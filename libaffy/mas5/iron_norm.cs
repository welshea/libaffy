/**************************************************************************
 *
 * Filename:  iron_norm.c
 *
 * Purpose:   libaffy wrappers around pairwise normalization code.
 *
 * Creation:  2012/03/16
 *
 * Author:    Eric Welsh
 *
 * Copyright: Copyright (C) 2012, Moffitt Cancer Center
 *            All rights reserved.
 *
 * Update History
 * --------------
 * 2012/03/16: Creation (EAW)
 * 2013/01/10: fixed typo in numprobes2 assignment that resulted in it using
 *             an uninitialized value during PM-only normaliztation (EAW)
 * 2014/03/14: fixed to work with exon arrays (EAW)
 * 2018/06/01: added support for probeset exclusions during IRON training (EAW)
 * 2018/09/11: extend exclusions to probeset-level normalization (EAW)
 * 2018/09/11: support leaving exclusions unscaled (EAW)
 * 2018/09/14: split unscaled support into new spikeins related code (EAW)
 * 2019/03/13: skip AFFX probesets during probeset-level normalization (EAW)
 * *** TODO -- fix 1:many probe:probeset stuff ***
 *
 **************************************************************************/

#include "affy.h"

#define MIN_SIGNAL       1E-5
#define DO_FLOOR         0	/* do not enable, is handled elsewhere now */

#define RANK_FRACTION    0.01           /* 0.01 */
#define RANK_FRACTION2   0.10           /* 0.10 */

void affy_pairwise_normalization(AFFY_CHIPSET *cs,
                                 AFFY_CHIP *model_chip,
                                 unsigned int opts,
                                 AFFY_COMBINED_FLAGS *f,
                                 AFFY_ERROR *err)
{
  affy_int32    x, y, numprobes, numprobes2;
  double       *model_signals, *input_signals, *scale_factors;
  double        rmsd, frac;
  char         *mask_model, *mask, mask_char;
  int          *mempool;
  unsigned int  i, j;
  affy_uint32   p;
  AFFY_CDFFILE *cdf;
  AFFY_CELL   **chip_data;
  AFFY_CELFILE *model_cel;
  char         *filestem;

  assert(cs              != NULL);
  assert(cs->cdf         != NULL);
  assert(cs->chip        != NULL);
  assert(model_chip      != NULL);
  assert(model_chip->cel != NULL);

  if (cs->num_chips == 0)
    return;

  cdf        = cs->cdf;
  numprobes  = cdf->numprobes;
  numprobes2 = (opts & AFFY_PAIRWISE_PM_ONLY) ? numprobes : numprobes << 1;
  model_cel  = model_chip->cel;

  mempool = h_malloc(sizeof(int));
  if (mempool == NULL)
    AFFY_HANDLE_ERROR_VOID("malloc failed", AFFY_ERROR_OUTOFMEM, err);

  model_signals = h_suballoc(mempool, numprobes2 * sizeof(double));
  if (model_signals == NULL)
    AFFY_HANDLE_ERROR_GOTO("malloc failed",
                           AFFY_ERROR_OUTOFMEM,
                           err,
                           cleanup);

  input_signals = h_suballoc(mempool, numprobes2 * sizeof(double));
  if (input_signals == NULL)
    AFFY_HANDLE_ERROR_GOTO("malloc failed",
                           AFFY_ERROR_OUTOFMEM,
                           err,
                           cleanup);
  mask = h_suballoc(mempool, numprobes2 * sizeof(double));
  if (mask == NULL)
    AFFY_HANDLE_ERROR_GOTO("malloc failed",
                           AFFY_ERROR_OUTOFMEM,
                           err,
                           cleanup);

  mask_model = h_suballoc(mempool, numprobes2 * sizeof(double));
  if (mask_model == NULL)
    AFFY_HANDLE_ERROR_GOTO("malloc failed",
                           AFFY_ERROR_OUTOFMEM,
                           err,
                           cleanup);

  scale_factors = h_suballoc(mempool, numprobes2 * sizeof(double));
  if (scale_factors == NULL)
    AFFY_HANDLE_ERROR_GOTO("malloc failed",
                           AFFY_ERROR_OUTOFMEM,
                           err,
                           cleanup);

  /* NOTE -- there will be some minor issues with exclusions of probes
   *  that belong to multiple probesets.  If even one copy of a probe
   *  is not excluded, one copy is left unmasked.
   *
   *  Also, if a probe has duplicates that are in both spikein and
   *   non-spikein probesets, it *will* be normalized, rather
   *   than left unnormalized.  This is not likely to occur,
   *   since spikein probesets should generally not contain any
   *   probes that are present in non-spikein probesets.
   *
   *  However, if the pm[] vector is used, to keep separate copies of
   *   all the probes, then only the non-spikein copy will be normalized.
   *   The pm[] vector optimization is currently not used, though, so this
   *   exception currently does not trigger.
   */

  if (opts & AFFY_PAIRWISE_PM_ONLY)
  {
    model_signals = model_chip->pm;

    /* initialize stuff to deal with duplicate probes */
    memset(cdf->seen_xy[0], 0, cdf->numrows*cdf->numcols*sizeof(affy_uint8));
   
    for (p = 0; p < numprobes; p++)
    {
      x = cdf->probe[p]->pm.x;
      y = cdf->probe[p]->pm.y;
      
      mask_char = (bit_test(model_cel->mask[x], y)) ||
        (cdf->cell_type[x][y] == AFFY_UNDEFINED_LOCATION) ||
        (cdf->cell_type[x][y] == AFFY_QC_LOCATION);
      
      /* skip AFFX/control probesets */
      if (affy_is_control_probe(&cdf->probe[p]))
        mask_char = 1;

      /* mask low intensity points */
      if (model_signals[p] < 1)
        mask_char = 1;

      /* mask probesets that we want to exclude from training */
      if (f->use_exclusions && cdf->exclusions)
      {
        if (bsearch(&cdf->probe[p]->ps->name, &cdf->exclusions[0],
            cdf->numexclusions, sizeof(char *), compare_string))
        {
          mask_char = 1;
        }
      }
      /* we also want to exclude spikins */
      if (f->use_spikeins && cdf->spikeins)
      {
        if (bsearch(&cdf->probe[p]->ps->name, &cdf->spikeins[0],
            cdf->numspikeins, sizeof(char *), compare_string))
        {
          mask_char = 1;
        }
      }

      /* mask all but the first unmasked copy of duplicate probes */
      if (mask_char == 0 && cdf->seen_xy[x][y] == 0)
        cdf->seen_xy[x][y] = 1;
      else if (cdf->seen_xy[x][y])
        mask_char = 1;

      mask_model[p] = mask_char;
    }
  }
  else
  {
    memset(cdf->seen_xy[0], 0, cdf->numrows*cdf->numcols*sizeof(affy_uint8));

    for (p = 0, j = 0; p < numprobes; p++)
    {
      x = cdf->probe[p]->pm.x;
      y = cdf->probe[p]->pm.y;
      
      model_signals[j] = model_cel->data[x][y].value;
      mask_char = (bit_test(model_cel->mask[x], y)) ||
        (cdf->cell_type[x][y] == AFFY_UNDEFINED_LOCATION) ||
        (cdf->cell_type[x][y] == AFFY_QC_LOCATION);

      /* skip AFFX/control probesets */
      if (affy_is_control_probe(&cdf->probe[p]))
        mask_char = 1;

      /* mask low intensity points */
      if (model_signals[j] < 1)
        mask_char = 1;

      /* mask probesets that we want to exclude from training */
      if (f->use_exclusions && cdf->exclusions)
      {
        if (bsearch(&cdf->probe[p]->ps->name, &cdf->exclusions[0],
            cdf->numexclusions, sizeof(char *), compare_string))
        {
          mask_char = 1;
        }
      }
      /* we also want to exclude spikeins */
      if (f->use_spikeins && cdf->spikeins)
      {
        if (bsearch(&cdf->probe[p]->ps->name, &cdf->spikeins[0],
            cdf->numspikeins, sizeof(char *), compare_string))
        {
          mask_char = 1;
        }
      }

      /* mask all but the first unmasked copy of duplicate probes */
      if (mask_char == 0 && cdf->seen_xy[x][y] == 0)
        cdf->seen_xy[x][y] = 1;
      else if (cdf->seen_xy[x][y])
        mask_char = 1;

      mask_model[j] = mask_char;

      j++;

     
      /* hack for missing MM probes, where MM coords == PM coords
       */
      if (cdf->probe[p]->pm.x == cdf->probe[p]->mm.x &&
          cdf->probe[p]->pm.y == cdf->probe[p]->mm.y)
      {
        continue;
      }
      else
      {
        x = cdf->probe[p]->mm.x;
        y = cdf->probe[p]->mm.y;

        /* do not train on MM probes */
        model_signals[j] = model_cel->data[x][y].value;
        mask_model[j] = 1;

        j++;
      }
    }
  }

  for (i = 0; i < cs->num_chips; i++)
  {
    chip_data = cs->chip[i]->cel->data;
    filestem  = stem_from_filename_safer(cs->chip[i]->filename);

    if (opts & AFFY_PAIRWISE_PM_ONLY)
    {
      input_signals = cs->chip[i]->pm;
      /*
       * as of right now, our RMA ignores masks, should we honor them here?
       * Yes (EAW)
       */
      for (p = 0; p < numprobes; p++)
      {
        mask[p] = mask_model[p] | bit_test(cs->chip[i]->cel->mask[x], y);

        /* mask low intensity points */
        if (input_signals[p] < 1)
          mask[p] = 1;
      }
    }
    else
    {
      for (p = 0, j = 0; p < numprobes; p++)
      {
        x = cdf->probe[p]->pm.x;
        y = cdf->probe[p]->pm.y;

        input_signals[j] = chip_data[x][y].value;
        mask[j] = mask_model[j] | bit_test(cs->chip[i]->cel->mask[x], y);

        /* mask low intensity points */
        if (input_signals[j] < 1)
          mask[j] = 1;

        j++;
       
        /* hack for missing MM probes, where MM coords == PM coords
         */
        if (cdf->probe[p]->pm.x == cdf->probe[p]->mm.x &&
            cdf->probe[p]->pm.y == cdf->probe[p]->mm.y)
        {
          continue;
        }
        else
        {
          x = cdf->probe[p]->mm.x;
          y = cdf->probe[p]->mm.y;

          /* do not train on MM probes */
          input_signals[j] = chip_data[x][y].value;
          mask[j] = 1;

          j++;
        }
      }

      numprobes2 = j;
    }
    
    /* Rank Fraction Cutoff = 0.01 */
    fill_normalization_scales(filestem,
                              model_signals,
                              input_signals,
                              scale_factors,
                              mask,
                              numprobes2,
                              RANK_FRACTION,
                              RANK_FRACTION2,
                              f,
                              &frac,
                              &rmsd,
                              err);
    AFFY_CHECK_ERROR_GOTO(err, cleanup);

    fprintf(stderr, "pairwise\tprobe-level\t%s\t%s\t%f\t%f\t%f\n",
            model_chip->filename,
            cs->chip[i]->filename,
            frac, rmsd, (rmsd + 1E-5) / (frac + 1E-5));

    if (opts & AFFY_PAIRWISE_PM_ONLY)
    {
      for (p = 0; p < numprobes; p++)
      {
        /* preserve missing data */
        if (cs->chip[i]->pm[p])
        {
          /* HACK -- set scaling factors for spikein probesets to 1 */
          if (f->use_spikeins && cdf->spikeins)
          {
            if (bsearch(&cdf->probe[p]->ps->name, &cdf->spikeins[0],
                cdf->numspikeins, sizeof(char *), compare_string))
            {
              scale_factors[p] = 1.0;
            }
          }

          if (scale_factors[p] > 0)
            cs->chip[i]->pm[p] *= scale_factors[p];

#if DO_FLOOR
          if (cs->chip[i]->pm[p] < MIN_SIGNAL)
            cs->chip[i]->pm[p] = MIN_SIGNAL;
#endif
        }
      }
    }
    else
    {
      memset(cdf->seen_xy[0], 0, cdf->numrows*cdf->numcols*sizeof(affy_uint8));

      for (p = 0, j = 0; p < numprobes; p++)
      {
        x = cdf->probe[p]->pm.x;
        y = cdf->probe[p]->pm.y;

        /* preserve missing data */
        if (chip_data[x][y].value)
        {
          /* HACK -- set scaling factors for spikein probesets to 1 */
          if (f->use_spikeins && cdf->spikeins)
          {
            if (bsearch(&cdf->probe[p]->ps->name, &cdf->spikeins[0],
                cdf->numspikeins, sizeof(char *), compare_string))
            {
              scale_factors[j] = 1.0;
            }
          }
          
          if (scale_factors[j] > 0 && scale_factors[j] != 1.0)
          {
            /* only scale a duplicate probe once */
            if (cdf->seen_xy[x][y] == 0)
            {
              chip_data[x][y].value *= scale_factors[j];
              cdf->seen_xy[x][y] = 1;
            }
          }

#if DO_FLOOR
          if (chip_data[x][y].value < MIN_SIGNAL)
            chip_data[x][y].value = MIN_SIGNAL;
#endif
        }

        j++;

        /* hack for missing MM probes, where MM coords == PM coords
         */
        if (cdf->probe[p]->pm.x == cdf->probe[p]->mm.x &&
            cdf->probe[p]->pm.y == cdf->probe[p]->mm.y)
        {
          continue;
        }
        else
        {
          x = cdf->probe[p]->mm.x;
          y = cdf->probe[p]->mm.y;

          /* preserve missing data */
          if (chip_data[x][y].value)
          {
            /* HACK -- set scaling factors for spikein probesets to 1 */
            if (f->use_spikeins && cdf->spikeins)
            {
              if (bsearch(&cdf->probe[p]->ps->name, &cdf->spikeins[0],
                  cdf->numspikeins, sizeof(char *), compare_string))
              {
                scale_factors[j] = 1.0;
              }
            }


            if (scale_factors[j] > 0 && scale_factors[j] != 1.0)
            {
              /* only scale a duplicate probe once */
              if (cdf->seen_xy[x][y] == 0)
              {
                chip_data[x][y].value *= scale_factors[j];
                cdf->seen_xy[x][y] = 1;
              }
            }

#if DO_FLOOR
            if (chip_data[x][y].value < MIN_SIGNAL)
              chip_data[x][y].value = MIN_SIGNAL;
#endif
          }

          j++;
        }
      }
    }
  }

cleanup:
  h_free(mempool);
}


void affy_pairwise_normalization_probeset(AFFY_CHIPSET *cs,
                                          AFFY_CHIP *model_chip,
                                          int unlog_flag,
                                          AFFY_COMBINED_FLAGS *f,
                                          AFFY_ERROR *err)
{
  affy_int32    numprobesets;
  double       *model_signals, *input_signals, *scale_factors;
  double        rmsd, frac;
  double        log2 = log(2.0);
  char         *mask;
  int          *mempool;
  unsigned int  i;
  affy_uint32   p;
  AFFY_CDFFILE *cdf;
  char         *filestem;

#if DEBUG_SKIP_PROBESET_NORM
  /* skip probeset normalization, since we're generating probe graphs */
  return;
#endif

  assert(cs              != NULL);
  assert(cs->cdf         != NULL);
  assert(cs->chip        != NULL);
  assert(model_chip      != NULL);
  assert(model_chip->cel != NULL);

  if (cs->num_chips == 0)
    return;

  cdf           = cs->cdf;
  numprobesets  = cdf->numprobesets;

  mempool = h_malloc(sizeof(int));
  if (mempool == NULL)
    AFFY_HANDLE_ERROR_VOID("malloc failed", AFFY_ERROR_OUTOFMEM, err);

  mask = h_suballoc(mempool, numprobesets * sizeof(double));
  if (mask == NULL)
    AFFY_HANDLE_ERROR_GOTO("malloc failed",
                            AFFY_ERROR_OUTOFMEM,
                            err,
                            cleanup);

  scale_factors = h_suballoc(mempool, numprobesets * sizeof(double));
  if (scale_factors == NULL)
    AFFY_HANDLE_ERROR_GOTO("malloc failed",
                           AFFY_ERROR_OUTOFMEM,
                           err,
                           cleanup);

  memset(mask, 0, numprobesets * sizeof(char));

  model_signals  = model_chip->probe_set;
  if (unlog_flag)
    for (p = 0; p < numprobesets; p++)
      model_signals[p] = pow(2.0, model_signals[p]);

  for (i = 0; i < cs->num_chips; i++)
  {
    input_signals = cs->chip[i]->probe_set;
    filestem      = stem_from_filename_safer(cs->chip[i]->filename);

    if (unlog_flag)
      for (p = 0; p < numprobesets; p++)
        input_signals[p] = pow(2, input_signals[p]);

    memset(mask, 0, numprobesets * sizeof(char));
    for (p = 0; p < numprobesets; p++)
    {
      /* mask obvious AFFX control probesets */
      if (affy_is_control_probeset(&cdf->probeset[p]))
        mask[p] = 1;
    
      /* mask low intensity points */
      if (model_signals[p] < 1 || input_signals[p] < 1)
        mask[p] = 1;

      /* mask probesets that we want to exclude from training */
      if (f->use_exclusions && cdf->exclusions)
      {
        if (bsearch(&cdf->probeset[p].name, &cdf->exclusions[0],
            cdf->numexclusions, sizeof(char *), compare_string))
        {
          mask[p] = 1;
        }
      }
      /* we also want to exclude spikeins */
      if (f->use_spikeins && cdf->spikeins)
      {
        if (bsearch(&cdf->probeset[p].name, &cdf->spikeins[0],
            cdf->numspikeins, sizeof(char *), compare_string))
        {
          mask[p] = 1;
        }
      }
    }

    /* Rank Fraction Cutoff = 0.01 */
    fill_normalization_scales(filestem,
                              model_signals,
                              input_signals,
                              scale_factors,
                              mask,
                              numprobesets,
                              RANK_FRACTION,
                              RANK_FRACTION2,
                              f,
                              &frac,
                              &rmsd,
                              err);
    AFFY_CHECK_ERROR_GOTO(err, cleanup);

    fprintf(stderr, "pairwise\tprobeset-level\t%s\t%s\t%f\t%f\t%f\n",
            model_chip->filename,
            cs->chip[i]->filename,
            frac, rmsd, (rmsd + 1E-5) / (frac + 1E-5));

    for (p = 0; p < numprobesets; p++)
    {
      /* HACK -- set scaling factors for spikein probesets to 1 */
      if (f->use_spikeins && cdf->spikeins)
      {
        if (bsearch(&cdf->probeset[p].name, &cdf->spikeins[0],
            cdf->numspikeins, sizeof(char *), compare_string))
        {
            scale_factors[p] = 1.0;
        }
      }
    
      if (scale_factors[p] > 0)
        input_signals[p] *= scale_factors[p];

#if 0
      if (input_signals[p] < MIN_SIGNAL)
        input_signals[p] = MIN_SIGNAL;
#else
      if (input_signals[p] < 1.0)
        input_signals[p] = 1.0;
#endif
    }

    if (unlog_flag)
      for (p = 0; p < numprobesets; p++)
        input_signals[p] = log(input_signals[p]) / log2;
  }

cleanup:
  h_free(mempool);
}


void affy_floor_probe(AFFY_CHIPSET *cs,
                      double floor_value,
                      AFFY_ERROR *err)
{
  affy_int32    x, y, numprobes;
  unsigned int  i;
  affy_uint32   p;
  AFFY_CDFFILE *cdf;
  AFFY_CELL   **chip_data;

  assert(cs              != NULL);
  assert(cs->cdf         != NULL);
  assert(cs->chip        != NULL);

  if (cs->num_chips == 0)
    return;

  cdf       = cs->cdf;
  numprobes = cdf->numprobes;

  for (i = 0; i < cs->num_chips; i++)
  {
    chip_data = cs->chip[i]->cel->data;

    if (cs->chip[i]->pm)
    {
      for (p = 0; p < numprobes; p++)
      {
        if (cs->chip[i]->pm[p] < floor_value)
          cs->chip[i]->pm[p] = floor_value;
      }
    }
    if (cdf->probe && chip_data)
    {
      memset(cdf->seen_xy[0], 0, cdf->numrows*cdf->numcols*sizeof(affy_uint8));

      for (p = 0; p < numprobes; p++)
      {
        x = cdf->probe[p]->pm.x;
        y = cdf->probe[p]->pm.y;

        if (cdf->seen_xy[x][y] == 0)
        {
          if (chip_data[x][y].value < floor_value)
            chip_data[x][y].value = floor_value;
        }
        cdf->seen_xy[x][y] = 1;

        /* hack for missing MM probes, where MM coords == PM coords
         */
        if (cdf->probe[p]->pm.x == cdf->probe[p]->mm.x &&
            cdf->probe[p]->pm.y == cdf->probe[p]->mm.y)
        {
          continue;
        }
        else
        {
          x = cdf->probe[p]->mm.x;
          y = cdf->probe[p]->mm.y;

          if (cdf->seen_xy[x][y] == 0)
          {
            if (chip_data[x][y].value < floor_value)
              chip_data[x][y].value = floor_value;
          }
          cdf->seen_xy[x][y] = 1;
        }
      }
    }
  }
}


void affy_floor_probeset(AFFY_CHIPSET *cs,
                         double floor_value,
                         AFFY_ERROR *err)
{
  affy_int32    numprobesets;
  double       *input_signals;
  unsigned int  i;
  affy_uint32   p;
  AFFY_CDFFILE *cdf;

  assert(cs              != NULL);
  assert(cs->cdf         != NULL);
  assert(cs->chip        != NULL);

  if (cs->num_chips == 0)
    return;

  cdf           = cs->cdf;
  numprobesets  = cdf->numprobesets;

  for (i = 0; i < cs->num_chips; i++)
  {
    input_signals = cs->chip[i]->probe_set;

    for (p = 0; p < numprobesets; p++)
      if (input_signals[p] < floor_value)
        input_signals[p] = floor_value;
  }
}


/* considers <= MIN_SIGNAL to be zero */
void affy_floor_probeset_to_min_non_zero(AFFY_CHIPSET *cs,
                                         AFFY_ERROR *err)
{
  affy_int32    numprobesets;
  double       *input_signals;
  double        min;
  unsigned int  i;
  affy_uint32   p;
  AFFY_CDFFILE *cdf;

  assert(cs              != NULL);
  assert(cs->cdf         != NULL);
  assert(cs->chip        != NULL);

  if (cs->num_chips == 0)
    return;

  cdf           = cs->cdf;
  numprobesets  = cdf->numprobesets;

  for (i = 0; i < cs->num_chips; i++)
  {
    input_signals = cs->chip[i]->probe_set;
    min = 9E99;

    for (p = 0; p < numprobesets; p++)
      if (input_signals[p] > MIN_SIGNAL && input_signals[p] < min)
        min = input_signals[p];

    for (p = 0; p < numprobesets; p++)
      if (input_signals[p] < min)
        input_signals[p] = min;
  }
}


void affy_unlog_probeset(AFFY_CHIPSET *cs, AFFY_ERROR *err)
{
  affy_int32    numprobesets;
  double       *input_signals;
  unsigned int  i;
  affy_uint32   p;
  AFFY_CDFFILE *cdf;

  assert(cs              != NULL);
  assert(cs->cdf         != NULL);
  assert(cs->chip        != NULL);

  if (cs->num_chips == 0)
    return;

  cdf           = cs->cdf;
  numprobesets  = cdf->numprobesets;

  for (i = 0; i < cs->num_chips; i++)
  {
    input_signals = cs->chip[i]->probe_set;

    for (p = 0; p < numprobesets; p++)
      input_signals[p] = pow(2.0, input_signals[p]);
  }
}


void affy_probeset_weighted_mean(AFFY_CHIP *chip, AFFY_ERROR *err)
{
  struct weight_stuff
  {
    double     weight;
    affy_int32 n_windows;
  };

  struct weight_stuff *weights = NULL;
  double              *log_signals = NULL;
  double              *input_signals;
  double               avg, weight, diff, min_weight, max_weight;
  double               weight_sum;
  affy_int32           numprobesets;
  affy_uint32          p, i, w;
  

  assert(chip              != NULL);
  assert(chip->cdf         != NULL);
/*
  assert(cs->chip        != NULL);
  if (cs->num_chips == 0)
    return;
*/


  numprobesets  = chip->cdf->numprobesets;
  w = (int) (numprobesets * 0.01 + 0.5);
  
  weights     = (struct weight_stuff *) calloc(numprobesets,
                                               sizeof(struct weight_stuff));
  log_signals = (double *) calloc(numprobesets, sizeof(double));

  input_signals = chip->probe_set;
  
  for (p = 0; p < numprobesets; p++)
    log_signals[p] = log(input_signals[p]);

  /* slide bin windows */
  for (p = 0; p <= numprobesets - w; p++)
  {
    avg = 0;
    for (i = p; i < p + w; i++)
      avg += log_signals[p];
    avg /= w;
    
    weight = 0;
    for (i = p; i < p + w; i++)
    {
      diff = log_signals[p] - avg;
      weight += diff * diff;
    }
    weight = sqrt(weight / w);

    for (i = p; i < p + w; i++)
    {
      weights[i].weight += weight;
      weights[i].n_windows++;
    }
  }
  
  min_weight =  9E99;
  max_weight = -9E99;
  for (p = 0; p < numprobesets; p++)
  {
    weights[p].weight /= weights[p].n_windows;
    
    if (weights[p].weight > 0 && weights[p].weight < min_weight)
      min_weight = weights[p].weight;

    if (weights[p].weight > max_weight)
      max_weight = weights[p].weight;
  }
  
  for (p = 0; p < numprobesets; p++)
  {
    if (weights[p].weight == 0.0)
      weights[p].weight = min_weight;
    
    weights[p].weight = pow(min_weight / weights[p].weight, 4.0);
  }
  
  avg = 0;
  weight_sum = 0.0;

  for (p = 0; p < numprobesets; p++)
  {
    weight      = weights[p].weight;
    avg        += weight * log_signals[p];
    weight_sum += weight;
  }
  
  avg /= weight_sum;
  
  printf("WeightedAvg:\t%f\n", exp(avg));
  
  if (weights) free(weights);
  if (log_signals) free(log_signals);
}
