
/**************************************************************************
 *
 * Filename:  mas5_scale.c
 *
 * Purpose: 
 *
 * Creation: 
 *
 * Author:    Steven Eschrich
 *
 * Copyright: Copyright (C) 2007, Moffitt Cancer Center 
 *            All rights reserved. 
 *
 * Update History
 * --------------
 * 04/14/05: Imported/repaired from old libaffy (AMH)
 * 03/07/08: New error handling scheme (AMH)
 * 09/20/10: Pooled memory allocator (AMH)
 * 10/22/10: Use new AFFY_COMBINED_FLAGS instead of AFFY_MAS5_FLAGS
 * 01/06/12: Added global scaling for quantile pre-normalized probesets
 * 03/06/14: Only include points >= 0 in the averages
 *
 **************************************************************************/

#include <affy_mas5.h>

static double trimmed_mean(double *values, int n, double lo, double hi);

int affy_mas5_scale(AFFY_CHIPSET *c, AFFY_COMBINED_FLAGS *f, AFFY_ERROR *err)
{
  int               i, n, num_probesets;
  int               n_above_zero;
  double           *signal, sf = 0;
  double            avg_sf;
  double            value;
  LIBUTILS_PB_STATE pbs;

  if (c == 0 || c->num_chips == 0 || c->cdf == 0)
    return (-1);

  pb_init(&pbs);
  pb_begin(&pbs, c->num_chips, "Scaling probeset values to %5.lf", f->scale_target);

  if (f == NULL)
    f = affy_mas5_get_defaults(err);

  AFFY_CHECK_ERROR(err, -1);

  num_probesets = c->cdf->numprobesets;

  signal = h_calloc(num_probesets, sizeof(double));
  if (signal == NULL)
    AFFY_HANDLE_ERROR("calloc failed", AFFY_ERROR_OUTOFMEM, err, -1);

  /* scale all chips by the same factor, average of all individual factors */
  if (f->use_quantile_normalization)
  {
    avg_sf = 0;

    /* calculate average scaling factor across all chips */
    for (n = 0; n < c->num_chips; n++)
    {
      /* Load signal and get scaling factor */
      n_above_zero = 0;
      for (i = 0; i < num_probesets; i++)
      {
        value = c->chip[n]->probe_set[i];
        if (value > 0)
          signal[n_above_zero++] = value;
      }
      sf =
        f->scale_target / trimmed_mean(signal, n_above_zero,
                                       f->trimmed_mean_low,
                                       f->trimmed_mean_high);
      avg_sf += sf;
    }
    
    if (c->num_chips)
      avg_sf /= c->num_chips;

    /* scale each individual chip */
    for (n = 0; n < c->num_chips; n++)
    {
      pb_tick(&pbs,1,"Scaling chip %d",n+1);
      info(" Sample %s scaled, sf=%f ", c->chip[n]->filename, sf);

      /* Store adjusted value */
      for (i = 0; i < num_probesets; i++)
      {
        c->chip[n]->probe_set[i] *= avg_sf;
      }
    }
  }
  /* scale each separately as per standard MAS5 */
  else
  {
    for (n = 0; n < c->num_chips; n++)
    {
      pb_tick(&pbs,1,"Scaling chip %d",n+1);

      /* Load signal and get scaling factor */
      n_above_zero = 0;
      for (i = 0; i < num_probesets; i++)
      {
        value = c->chip[n]->probe_set[i];
        if (value > 0)
          signal[n_above_zero++] = value;
      }
      sf =
        f->scale_target / trimmed_mean(signal, n_above_zero,
                                       f->trimmed_mean_low,
                                       f->trimmed_mean_high);

      info(" Sample %s scaled, sf=%f ", c->chip[n]->filename, sf);

      /* Store adjusted value */
      for (i = 0; i < num_probesets; i++)
      {
        c->chip[n]->probe_set[i] *= sf;
      }
    }
  }

  pb_finish(&pbs, "Finished MAS5 linear probeset scaling");

  h_free(signal);

  return (0);
}


int affy_mas5_scale_iron(AFFY_CHIPSET *c, AFFY_CHIPSET *model_chipset,
                         AFFY_COMBINED_FLAGS *f, AFFY_ERROR *err)
{
  int               i, n, num_probesets;
  int               n_above_zero;
  double           *signal, sf = 0;
  double            value;
  LIBUTILS_PB_STATE pbs;

  if (c == 0 || c->num_chips == 0 || c->cdf == 0)
    return (-1);
  if (model_chipset == 0 || model_chipset->num_chips == 0 ||
      model_chipset->cdf == 0)
  {
      return (-1);
  }

  pb_init(&pbs);
  pb_begin(&pbs, c->num_chips, "Scaling probeset values to %5.lf", f->scale_target);

  if (f == NULL)
    f = affy_mas5_get_defaults(err);

  AFFY_CHECK_ERROR(err, -1);

  num_probesets = c->cdf->numprobesets;

  signal = h_calloc(num_probesets, sizeof(double));
  if (signal == NULL)
    AFFY_HANDLE_ERROR("calloc failed", AFFY_ERROR_OUTOFMEM, err, -1);

  /* scale all chips by the same factor, calculated from the model chipset */
  n_above_zero = 0;
  for (i = 0; i < num_probesets; i++)
  {
    value = model_chipset->chip[0]->probe_set[i];
    if (value > 0)
      signal[n_above_zero++] = value;
  }
  sf = f->scale_target / trimmed_mean(signal, n_above_zero,
                                      f->trimmed_mean_low,
                                      f->trimmed_mean_high);

  for (n = 0; n < c->num_chips; n++)
  {
    pb_tick(&pbs,1,"Scaling chip %d",n+1);
    info(" Sample %s scaled, sf=%f ", c->chip[n]->filename, sf);

    /* Store adjusted value */
    for (i = 0; i < num_probesets; i++)
    {
      c->chip[n]->probe_set[i] *= sf;
    }
  }

  pb_finish(&pbs, "Finished MAS5 linear probeset scaling");

  h_free(signal);

  return (0);
}


/* 
 * Return the trimmed mean of a list of numbers, trimmed by the
 * top lo and hi percentages.
 */
static double trimmed_mean(double *values, int n, double lo, double hi)
{
  double sum = 0;
  int    min_value, max_value, i;

  assert(values != NULL);

  /* Calculate number to start/stop at in sorted list */
  min_value = (int)(n * lo);
  max_value = (int)(n * hi);

  /* Sort values */
  qsort(values, n, sizeof(double), dcompare);

  /* Accumulate sums */
  for (i = min_value; i <= max_value; i++)
  {
    sum += values[i];
  }

  /* Divide by num signals */
  sum = sum / (max_value - min_value + 1);

  return (sum);
}
