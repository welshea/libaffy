
/**************************************************************************
 *
 * Filename:  mas5_signal.c
 *
 * Purpose:   Calculate signal values using the MAS5 algorithm.
 *
 * Creation: 
 *
 * Author:    Steven Eschrich
 *
 * Copyright: Copyright (C) 2007, Moffitt Cancer Center.
 *            All rights reserved.
 *
 * Update History
 * --------------
 * 04/14/05: Imported/repaired from old libaffy (AMH)
 * 04/19/05: AFFY_CELFILE uses AFFY_CELL's now (AMH)
 * 03/07:08: New error handling scheme (AMH)
 * 09/16/10: Add support for P/A calls (EAW)
 * 09/20/10: Pooled memory allocator (AMH)
 * 10/22/10: Use new AFFY_COMBINED_FLAGS instead of AFFY_MAS5_FLAGS
 * 09/19/12: Handle some rare special Tukey's Biweight cases (EAW)
 * 03/06/14: Skip MM subtraction when chip is missing MM probes (EAW)
 *
 **************************************************************************/

#include "affy_mas5.h"

/* LOG2 == ln(2) */
static const double LOG2    = 0.69314718055994530941;
static const int    c       = 5;
static const double epsilon = 0.0001;

/* Private routines */

static int    double_sort_compare(const void *p1, const void *p2);
static double median(double *x, int n, double *range_ptr, AFFY_ERROR *err);
static double tukey_biweight(double *x, int n, AFFY_ERROR *err);
static double calculate_specific_background(double *pm, 
                                            double *mm, 
                                            int n, 
                                            AFFY_ERROR *err);
static double calculate_probeset_signal(AFFY_CHIP *c, 
                                        int probeset_num,
					AFFY_COMBINED_FLAGS *f,
                                        AFFY_ERROR *err);


/* *** DO NOT USE *** -- experimental development code
 * also, haven't updated the pm/mm probe to deal with cdf->seen_xy[x][y] yet...
 */
static double calculate_probeset_signal_iron(AFFY_CHIP *c, 
                                             int probeset_num,
					     AFFY_COMBINED_FLAGS *f,
                                             AFFY_ERROR *err)
{
  AFFY_PROBESET *p;
  AFFY_CELL    **data;
  double        *pm, *mm;
  double         signal, signal_log_value_pm, signal_log_value_mm;
  int            n, i, j, x, y, *mempool;
  double         r;
  
  /* Some shortcuts */
  p    = &(c->cdf->probeset[probeset_num]);
  n    = p->numprobes;
  data = c->cel->data;

  mempool = h_malloc(sizeof(int));
  if (mempool == NULL)
    AFFY_HANDLE_ERROR("malloc failed", AFFY_ERROR_OUTOFMEM, err, 0.0);

  /* Assign pm/mm values for this chip (based on existing layout) */
  pm = h_subcalloc(mempool, n, sizeof(double));
  if (pm == NULL)
  {
    h_free(mempool);
    AFFY_HANDLE_ERROR("calloc failed", AFFY_ERROR_OUTOFMEM, err, 0.0);
  }

  mm = h_subcalloc(mempool, n, sizeof(double));
  if (mm == NULL)
  {
    h_free(mempool);
    AFFY_HANDLE_ERROR("calloc failed", AFFY_ERROR_OUTOFMEM, err, 0.0);
  }
  

  /*
     Assign pm/mm values based on probes. Caveat: masked probes (either 
     in PM or MM) are ignored in the computation.
   */
  for (i = 0, j = 0; i < n; i++)
  {
    x = p->probe[i].pm.x;
    y = p->probe[i].pm.y;
    if (affy_ismasked(c, x, y))
      continue;
    pm[j] = log(max_macro(data[x][y].value, f->delta)) / LOG2;

    x = p->probe[i].mm.x;
    y = p->probe[i].mm.y;
    if (affy_ismasked(c, x, y))
      continue;
    mm[j] = log(max_macro(data[x][y].value, f->delta)) / LOG2;
    j++;
  }

  /* Uh oh, all probes are masked!  We'll have to use them all... */
  if (j == 0)
  {
    for (i = 0, j = 0; i < n; i++)
    {
      x = p->probe[i].pm.x;
      y = p->probe[i].pm.y;
      pm[j++] = log(max_macro(data[x][y].value, f->delta)) / LOG2;;
    }
  }

  /* In case the total number of probes is less than n */
  if (n > j)
  {
    debug("%s:Adjusting for %d probes instead of %d", p->name, j, n);
    n = j;
  }
  
  /* FIXME -- cases where #mm == 0 are not handled properly */

  /* Tukey's Biweight signals (in log2 space) */
  signal_log_value_pm = tukey_biweight(pm, n, err);
  signal_log_value_mm = tukey_biweight(mm, n, err);
  
  /* correlate PM/MM vectors */
  r = calculate_pearson_r_double(pm, mm, n);
  /*  if (r < 0) r = 0; */

  /* subtract MM probeset signal in non-log space, weighted by r */
  signal = pow(2.0, signal_log_value_pm) -
               0.5 * (1 - r) * pow(2.0, signal_log_value_mm);

  signal = max_macro(signal, f->delta);

  h_free(mempool);

  AFFY_CHECK_ERROR(err, 0.0);

  return (signal);
}


static double calculate_probeset_signal(AFFY_CHIP *c, 
                                        int probeset_num,
					AFFY_COMBINED_FLAGS *f,
                                        AFFY_ERROR *err)
{
  AFFY_PROBESET *p;
  AFFY_CELL    **data;
  double        *pm, *pv;
  double         signal, signal_log_value;
  int            n, i, j, x, y, *mempool;
  
  /* Some shortcuts */
  p    = &(c->cdf->probeset[probeset_num]);
  n    = p->numprobes;
  data = c->cel->data;

  mempool = h_malloc(sizeof(int));
  if (mempool == NULL)
    AFFY_HANDLE_ERROR("malloc failed", AFFY_ERROR_OUTOFMEM, err, 0.0);

  /* Assign pm values for this chip (based on existing layout) */
  pm = h_subcalloc(mempool, n, sizeof(double));
  if (pm == NULL)
  {
    h_free(mempool);
    AFFY_HANDLE_ERROR("calloc failed", AFFY_ERROR_OUTOFMEM, err, 0.0);
  }

  pv = h_subcalloc(mempool, n, sizeof(double));
  if (pv == NULL)
  {
    h_free(mempool);
    AFFY_HANDLE_ERROR("calloc failed", AFFY_ERROR_OUTOFMEM, err, 0.0);
  }


  /*
     Assign pm values based on probes. Caveat: masked probes are ignored
     in the computation.
   */
  for (i = 0, j = 0; i < n; i++)
  {
    x = p->probe[i].pm.x;
    y = p->probe[i].pm.y;
    if (affy_ismasked(c, x, y))
      continue;
    pm[j++] = data[x][y].value;
  }
  
  /* Uh oh, all probes are masked.  We'll have to use all probes instead... */
  if (j == 0)
  {
    for (i = 0, j = 0; i < n; i++)
    {
      x = p->probe[i].pm.x;
      y = p->probe[i].pm.y;
      pm[j++] = data[x][y].value;
    }
  }

  /* In case the total number of probes is less than n */
  if (n > j)
  {
    debug("%s:Adjusting for %d probes instead of %d", p->name, j, n);
    n = j;
  }

  /* Calculate Probe Values (PV) */
  for (j = 0; j < n; j++)
    pv[j] = log(max_macro(pm[j], f->delta)) / LOG2;

  /* Signal log value is here */
  signal_log_value = tukey_biweight(pv, n, err);

  h_free(mempool);

  AFFY_CHECK_ERROR(err, 0.0);

  /* Then take the antilog and we're done */
  signal = pow(2, signal_log_value);

  return (signal);
}

/*----------------------------------------------------------------------
   The specific background is the Tukey Biweight of log differences
 -----------------------------------------------------------------------*/
static double calculate_specific_background(double *pm, 
                                            double *mm, 
                                            int n, 
                                            AFFY_ERROR *err)
{
  double *d;
  double  r;
  int     j;

  d = h_malloc(n * sizeof(double));
  if (d == NULL)
    AFFY_HANDLE_ERROR("malloc failed", AFFY_ERROR_OUTOFMEM, err, 0.0);

  /* Difference of log values */
  for (j = 0; j < n; j++)
    d[j] = (log(pm[j]) / LOG2 - log(mm[j]) / LOG2);

  /* Compute this biweight */
  r = tukey_biweight(d, n, err);

  h_free(d);

  AFFY_CHECK_ERROR(err, 0.0);

  return (r);
}

/*--------------------------------------------------------------------
  Calculate Tukey's Biweighted Average, per the Affy docs.
  --------------------------------------------------------------------*/
static double tukey_biweight(double *x, int n, AFFY_ERROR *err)
{
  double *u;
  int     i, *mempool;
  double  M, S;
  double  Tbi_num = 0, Tbi_denom = 0;
  double *diffs;
  double  range;
  
  /* special case for n == 1, to avoid any potential issues */
  if (n == 1)
  {
    return x[0];
  }
  
  /* special case for n == 2, to avoid potential epsilon-related errors */
  if (n == 2)
  {
    return (0.5 * (x[0] + x[1]));
  }

  mempool = h_malloc(sizeof(int));
  if (mempool == NULL)
    AFFY_HANDLE_ERROR("malloc failed", AFFY_ERROR_OUTOFMEM, err, 0.0);

  u = h_suballoc(mempool, n * sizeof(double));
  if (u == NULL)
  {
    h_free(mempool);
    AFFY_HANDLE_ERROR("malloc failed", AFFY_ERROR_OUTOFMEM, err, 0.0);
  }

  /* Calculate median */
  M = median(x, n, &range, err);
  if (err->type != AFFY_ERROR_NONE)
  {
    h_free(mempool);
    return (0.0);
  }
  
  /* zero variance, return first value */
  if (range <= DBL_EPSILON)
  {
    h_free(mempool);
    return x[0];
  }

  /* Calculate S, median of absolute differences from M */
  diffs = h_suballoc(mempool, n * sizeof(double));
  if (diffs == NULL)
  {
    h_free(mempool);
    AFFY_HANDLE_ERROR("malloc failed", AFFY_ERROR_OUTOFMEM, err, 0.0);
  }

  for (i = 0; i < n; i++)
    diffs[i] = fabs(x[i] - M);

  S = median(diffs, n, &range, err);

  if (err->type != AFFY_ERROR_NONE)
  {
    h_free(mempool);

    return (0.0);
  }

  /* Calculate distance measure */
  for (i = 0; i < n; i++)
    u[i] = (x[i] - M) / (c * S + epsilon);

  /* Finally, calculate the result */
  for (i = 0; i < n; i++)
  {
    double usquared, w;

    /* Function w(u) is 0 for all |u|>1 */
    if (fabs(u[i]) > 1)
      continue;

    usquared = u[i] * u[i];
    w = (1.0 - usquared) * (1.0 - usquared);
    Tbi_num += w * x[i];
    Tbi_denom += w;
  }
  
  /* If all points are distant, take a flat average instead */
  if (Tbi_denom <= DBL_EPSILON)
  {
    Tbi_denom = n;
    Tbi_num = 0.0;

    for (i = 0; i < n; i++)
      Tbi_num += x[i];
  }

  h_free(mempool);

  return (Tbi_num / Tbi_denom);
}

/*
 * Calculate the median of n numbers (x[0]..x[n-1]) without touching
 *  the x array.
 */
static double median(double *x, int n, double *range_ptr, AFFY_ERROR *err)
{
  double *d;
  int     i;
  double  M;

  /* Copy the array */
  d = h_malloc(n * sizeof(double));
  if (d == NULL)
    AFFY_HANDLE_ERROR("malloc failed", AFFY_ERROR_OUTOFMEM, err, 0.0);

  for (i = 0; i < n; i++)
    d[i] = x[i];

  /* Sort in increasing order */
  qsort(d, n, sizeof(double), double_sort_compare);

  /* Median is middle value, or mean of two middle values */
  if (n % 2 == 1)
    M = d[n / 2];
  else
    M = (d[(n / 2) - 1] + d[n / 2]) / 2.0;

  *range_ptr = d[n-1] - d[0];

  h_free(d);

  return (M);
}

static int double_sort_compare(const void *p1, const void *p2)
{
  double d1 = *((double *)p1);
  double d2 = *((double *)p2);

  if (d1 < d2)
    return (-1);
  if (d1 > d2)
    return (1);
  if (d1 == d2)
    return (0);

  return (-1);
}

int affy_mas5_signal(AFFY_CHIPSET *c, AFFY_COMBINED_FLAGS *f, AFFY_ERROR *err)
{
  int               i, n;
  int               num_probesets;
  LIBUTILS_PB_STATE pbs;

  assert(c      != NULL);
  assert(c->cdf != NULL);
  assert(f      != NULL);

  if (c->num_chips == 0)
    return (-1);

  pb_init(&pbs);

  num_probesets = c->cdf->numprobesets;

  pb_begin(&pbs, 
           c->num_chips*num_probesets, 
           "Calculating signal for probesets using Tukey's biweight method");
  for (n = 0; n < c->num_chips; n++)
  {
    c->chip[n]->probe_set = h_subcalloc(c->chip[n], 
                                        num_probesets, 
                                        sizeof(double));
    if (c->chip[n]->probe_set == NULL)
      AFFY_HANDLE_ERROR_GOTO("calloc failed", AFFY_ERROR_OUTOFMEM, err, err);

    c->chip[n]->numprobesets = num_probesets;

    for (i = 0; i < num_probesets; i++)
    {
      pb_tick(&pbs, 1, "Calculating probeset signal");
      c->chip[n]->probe_set[i] = calculate_probeset_signal(c->chip[n], 
                                                           i, 
                                                           f, 
                                                           err);
      AFFY_CHECK_ERROR_GOTO(err, err);
    }
  }

  pb_finish(&pbs, "Finished Tukey's Biweight probeset summarization");

  return (0);

err:
  return (-1);
}


int affy_iron_signal(AFFY_CHIPSET *c, AFFY_COMBINED_FLAGS *f, AFFY_ERROR *err)
{
  int               i, n;
  int               num_probesets;
  LIBUTILS_PB_STATE pbs;

  assert(c      != NULL);
  assert(c->cdf != NULL);
  assert(f      != NULL);

  if (c->num_chips == 0)
    return (-1);

  pb_init(&pbs);

  num_probesets = c->cdf->numprobesets;

  pb_begin(&pbs, 
           c->num_chips*num_probesets, 
           "Calculating signal for chip using IRON method");
  for (n = 0; n < c->num_chips; n++)
  {
    c->chip[n]->probe_set = h_subcalloc(c->chip[n], 
                                        num_probesets, 
                                        sizeof(double));
    if (c->chip[n]->probe_set == NULL)
      AFFY_HANDLE_ERROR_GOTO("calloc failed", AFFY_ERROR_OUTOFMEM, err, err);

    c->chip[n]->numprobesets = num_probesets;

    for (i = 0; i < num_probesets; i++)
    {
      pb_tick(&pbs, 1, "Calculating probeset signal");
      c->chip[n]->probe_set[i] = calculate_probeset_signal_iron(c->chip[n], 
                                                           i, 
                                                           f, 
                                                           err);
      AFFY_CHECK_ERROR_GOTO(err, err);
    }
  }

  pb_finish(&pbs, "Finished IRON probeset summarization");

  return (0);

err:
  return (-1);
}


int affy_mas5_subtract_mm_signal_probe(AFFY_CHIP *c,
                                       AFFY_COMBINED_FLAGS *f,
                                       AFFY_ERROR *err)
{
  AFFY_PROBESET *p;
  AFFY_CELL    **data;
  double        *pm = NULL, *mm = NULL;
  int            probeset_num, n, i, x, y;
  int            numprobesets;
  double         SB, im;
  LIBUTILS_PB_STATE pbs;
  
  /* chip is missing MM probes, abort */
  if (c->cdf->no_mm_flag == 1)
    return 1;
  
  pb_init(&pbs);
  data = c->cel->data;
  numprobesets = c->cdf->numprobesets;

  pb_begin(&pbs, 2, "MM Probe subtraction");

  for (probeset_num = 0; probeset_num < numprobesets; probeset_num++)
  {
    p    = &(c->cdf->probeset[probeset_num]);
    n    = p->numprobes;

    /* Assign pm/mm values for this chip (based on existing layout) */
    pm = realloc(pm, n * sizeof(double));
    if (pm == NULL)
      AFFY_HANDLE_ERROR("calloc failed", AFFY_ERROR_OUTOFMEM, err, 0);

    mm = realloc(mm, n * sizeof(double));
    if (mm == NULL)
    {
      free(mm);
      AFFY_HANDLE_ERROR("calloc failed", AFFY_ERROR_OUTOFMEM, err, 0);
    }

    /* Calculate the SB. NOTE: For some reason, it is calculated
     * on all PM/MM probe pairs, even if they are masked.
     */
    for (i = 0; i < n; i++)
    {
      x = p->probe[i].pm.x;
      y = p->probe[i].pm.y;
      pm[i] = data[x][y].value;

      /* hack for missing MM probes, where MM coords == PM coords
       */
      if (p->probe[i].pm.x == p->probe[i].mm.x &&
          p->probe[i].pm.y == p->probe[i].mm.y)
      {
        if (affy_ismasked(c, x, y))
          continue;
        mm[i] = 0;
      }
      else
      {
        x = p->probe[i].mm.x;
        y = p->probe[i].mm.y;
        mm[i] = data[x][y].value;
      }
    }

    SB = calculate_specific_background(pm, mm, n, err);

    if (err->type != AFFY_ERROR_NONE)
    {
      free(pm);
      free(mm);

      return (0);
    }

    /* Calculate Probe Values (PV) */
    for (i = 0; i < n; i++)
    {
      /* Easiest case, mm < pm */
      if (pm[i] > mm[i])
      {
        im = mm[i];
      }
      else
      {
        /* Based on SB, calculate IM */
        if (SB - f->contrast_tau > 0)
          im = pm[i] / pow(2, SB);
        else
          im =
            pm[i] / pow(2,
                        (f->contrast_tau /
                         (1 + ((f->contrast_tau - SB) / f->scale_tau))));
      }

      /* hack for missing MM probes, where MM coords == PM coords
       * leave background unchanged
       */
      if (p->probe[i].pm.x == p->probe[i].mm.x &&
          p->probe[i].pm.y == p->probe[i].mm.y)
      {
        continue;
      }
      /* handle MM as usual */
      else
      {
        x = p->probe[i].pm.x;
        y = p->probe[i].pm.y;
        data[x][y].value = pm[i] - im;
        
        x = p->probe[i].mm.x;
        y = p->probe[i].mm.y;
        data[x][y].value = 0;
      }
    }
  }

  pb_finish(&pbs, "Finished MM probe subtraction");
  
  free(pm);
  free(mm);

  AFFY_CHECK_ERROR(err, 0);
  
  return 1;
}
