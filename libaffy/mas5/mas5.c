
/**************************************************************************
 *
 * Filename:  mas5.c
 *
 * Purpose:   Top-level wrapper to tie together all the phases of MAS5.
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
 * 04/19/05: Use affy_free_cel_file() instead of manual code (AMH)
 * 10/17/06: Added option for quantile normalization (AMH)
 * 03/07/08: New error handling scheme (AMH)
 * 06/04/08: Revised CHIPSET handling (AMH)
 * 09/20/10: Pooled memory allocator (AMH)
 * 10/07/10: Pairwise normalization (AMH)
 * 10/20/10: Added 2nd pairwise normalization pass at probeset level (EAW)
 * 10/22/10: Use new AFFY_COMBINED_FLAGS instead of AFFY_MAS5_FLAGS
 * 03/04/11: Quantile normalization actually normalizes now (EAW)
 * 06/09/11: Moved MAS5 MM background subtraction to after normalization (EAW)
 * 06/10/11: Moved P/M/A calls to after normalization (EAW)
 * 08/05/11: Began adding support for separate PM/MM RMA bg (EAW)
 * 03/12/13: Added probeset-level quantile normalization pass (EAW)
 * 05/10/13: Replaced affy_free_cel_file() with
 *           affy_mostly_free_cel_file(), to allow for corrupt chip
 *           warnings at program end if --salvage was used (EAW)
 * 03/11/14: Added checks for no-MM chips (EAW)
 * 03/12/14: No longer use any pm[] array stuff, due to exon arrays (EAW)
 * 03/14/19: added support for exclusions and spikeins files (EAW)
 * 05/22/19: added --ignore-chip-mismatch support (EAW)
 * 08/12/20: pass flags to affy_create_chipset() (EAW)
 * 08/18/20: add flags to enable/disable iron or quantile probeset norm after
 *           probe norm (EAW)
 * 01/10/24: pass flags to affy_mean_normalization() (EAW)
 *
 **************************************************************************/

#include "affy_mas5.h"
#include "affy_rma.h"

#if TRAP_FLOAT_ERRORS
#define _GNU_SOURCE
#include <fenv.h>
#endif

extern void affy_iron_background_correction_probeset(AFFY_CHIPSET *cs,
                                              AFFY_COMBINED_FLAGS *f,
                                              AFFY_ERROR *err);

static void mas5_to_rma_pm(AFFY_CHIP *cp, AFFY_ERROR *err)
{
  affy_int32 numprobes, k;

  assert(cp            != NULL);
  assert(cp->cdf       != NULL);
  assert(cp->cel       != NULL);
  assert(cp->cel->data != NULL);

  numprobes = cp->cdf->numprobes;

  if (cp->pm == NULL)
    cp->pm = h_subcalloc(cp, numprobes, sizeof(double));
  if (cp->pm == NULL)
    AFFY_HANDLE_ERROR_VOID("malloc failed", AFFY_ERROR_OUTOFMEM, err);
  
  for (k = 0; k < numprobes; k++)
  {
    int x = cp->cdf->probe[k]->pm.x;
    int y = cp->cdf->probe[k]->pm.y;
    
    cp->pm[k] = cp->cel->data[x][y].value;
  }
}


/* free cel->data after converting to RMA PM-only array */
static void mas5_to_rma_pm_free(AFFY_CHIP *cp, AFFY_ERROR *err)
{
  affy_int32 numprobes, k;

  assert(cp            != NULL);
  assert(cp->cdf       != NULL);
  assert(cp->cel       != NULL);
  assert(cp->cel->data != NULL);

  numprobes = cp->cdf->numprobes;

  if (cp->pm == NULL)
    cp->pm = h_subcalloc(cp, numprobes, sizeof(double));
  if (cp->pm == NULL)
    AFFY_HANDLE_ERROR_VOID("malloc failed", AFFY_ERROR_OUTOFMEM, err);
  
  for (k = 0; k < numprobes; k++)
  {
    int x = cp->cdf->probe[k]->pm.x;
    int y = cp->cdf->probe[k]->pm.y;
    
    cp->pm[k] = cp->cel->data[x][y].value;
  }
  
  h_free(cp->cel->data);
  cp->cel->data = NULL;
}


static void rma_to_mas5_pm(AFFY_CHIP *cp, AFFY_ERROR *err)
{
  affy_int32 numprobes, k;

  assert(cp            != NULL);
  assert(cp->cdf       != NULL);
  assert(cp->cel       != NULL);
  assert(cp->cel->data != NULL);
  assert(cp->pm        != NULL);

  numprobes = cp->cdf->numprobes;

  for (k = 0; k < numprobes; k++)
  {
    int x, y;
    
    x = cp->cdf->probe[k]->pm.x;
    y = cp->cdf->probe[k]->pm.y;

    cp->cel->data[x][y].value = cp->pm[k];

    /* hack for missing MM probes, where MM coords == PM coords
     */
    if (cp->cdf->probe[k]->pm.x == cp->cdf->probe[k]->mm.x &&
        cp->cdf->probe[k]->pm.y == cp->cdf->probe[k]->mm.y)
    {
      continue;
    }
    else
    {
      x = cp->cdf->probe[k]->mm.x;
      y = cp->cdf->probe[k]->mm.y;

      cp->cel->data[x][y].value = 0;
    }
  }
}

AFFY_CHIPSET *affy_mas5(char **filelist, AFFY_COMBINED_FLAGS *f, AFFY_ERROR *err)
{
  AFFY_CHIPSET         *result, *temp, *model_chipset = NULL;
  AFFY_CHIP            *model_chip = NULL;
  AFFY_COMBINED_FLAGS  default_flags;
  int                  i, max_chips, chips_processed;
  char                 *chip_type = NULL, **p;

  assert(filelist != NULL);

#if TRAP_FLOAT_ERRORS
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  result = temp = NULL;
  chips_processed = 0;

  /* Get the array type. */
  chip_type = affy_get_cdf_name_from_cel(filelist[0], err);
  AFFY_CHECK_ERROR(err, NULL);

  /* In case flags not used, create default entry */
  if (f == NULL)
  {
    affy_rma_set_defaults(&default_flags);
    affy_mas5_set_defaults(&default_flags);
    f = &default_flags;
  }

  i = (f->use_quantile_normalization == true) 
      + (f->use_pairwise_normalization == true) 
      + (f->use_mean_normalization == true);
  if (i > 1)
  {
    AFFY_HANDLE_ERROR_GOTO("ERROR - Multiple normalization methods selected",
                           AFFY_ERROR_NOTSUPP, err, cleanup);
  }

  /* can't calculate P/M/A calls without mis-match probes */
  if (f->bg_rma && f->output_present_absent)
  {
    AFFY_HANDLE_ERROR_GOTO("ERROR - Can not calculate P/M/A calls with RMA background subtraction",
                           AFFY_ERROR_NOTSUPP, err, cleanup);
  }
  
  /* sanity checks */
  if (f->bg_iron)
  {
    f->use_mm_probe_subtraction = false;
  }

  result = affy_create_chipset(1, chip_type, f->cdf_directory, f, err);
  AFFY_CHECK_ERROR_GOTO(err, cleanup);


  /* load in probesets to exclude from IRON training */
  /* I haven't retrofitted this file to use much in the way of h_alloc(), etc.
   *  so I should eventually replace "result" with mempool as elsewhere
   */
  if (f->use_exclusions)
  {
    affy_load_exclusions_file(f->exclusions_filename, result->cdf,
                              result, err);
  }
  if (f->use_spikeins)
    affy_load_spikeins_file(f->spikeins_filename, result->cdf,
                            result, err);


  temp = affy_clone_chipset(result, err);
  AFFY_CHECK_ERROR_GOTO(err, cleanup);

  if (f->use_pairwise_normalization)
  {
    info("Loading pairwise normalization model from %s",
         f->pairwise_model_filename);

    model_chipset = affy_clone_chipset(result, err);
    AFFY_CHECK_ERROR_GOTO(err, cleanup);

    hattach(model_chipset, temp);

    affy_load_chipset_single(model_chipset, f->pairwise_model_filename,
                             f->ignore_chip_mismatch, err);
    AFFY_CHECK_ERROR_GOTO(err, cleanup);

    /* abort on corrupt CEL files, unless --salvage is used */
    if (model_chipset->chip[0]->cel->corrupt_flag &&
        f->salvage_corrupt == false)
          AFFY_HANDLE_ERROR_GOTO("corrupt CEL file",
                                 AFFY_ERROR_BADFORMAT,
                                 err,
                                 cleanup);

    model_chip = model_chipset->chip[0];

    if (f->use_background_correction)
    {
      if (f->bg_mas5)
      {
        affy_mas5_background_correction(model_chipset, f, err);
        AFFY_CHECK_ERROR_GOTO(err, cleanup);
      }
      else if (f->bg_rma || f->bg_rma_both)
      {
        if (f->use_mm_probe_subtraction &&
            model_chipset->chip[0]->cdf->no_mm_flag == 0)
        {
          affy_rma_background_correct_pm_mm_together(model_chipset,0,0,err);
        }
        else
        {
          if (f->bg_rma)
          {
            affy_rma_background_correct_pm_mm_together(model_chipset,0,1,err);
          }
          else
          {
            affy_rma_background_correct_pm_mm_together(model_chipset,0,0,err);
          }
        }
      }
      else if (f->bg_iron)
      {
        if (model_chipset->chip[0]->cdf->no_mm_flag == 0)
        {
          affy_rma_background_correct_pm_mm_together(model_chipset,0,0,err);
        }
        else
        {
          affy_rma_background_correct_pm_mm_together(model_chipset,0,0,err);
        }
      }
      else if (f->bg_global)
      {
        affy_global_background_correct(model_chipset,0,err);
      }
    }

    info("Pairwise reference sample loaded");
  }

  /* Count up chips */
  for (p = filelist, max_chips = 0; *p != NULL; p++)
    max_chips++;

  result = affy_resize_chipset(result, max_chips, err);
  AFFY_CHECK_ERROR_GOTO(err, cleanup);

  /* Load each chip */
  for (i = 0; i < max_chips; i++)
  {
    /*
     * The idea here is to load each chip into the result chipset,
     * copy its ptr to the temporary "singleton" chipset which is then
     * passed to the various processing routines (to operate on one
     * chip at a time).  We're left with the finished product in
     * `result'.
     */

    /* Load chip, skipping CEL files that can't be loaded */
    affy_load_chipset_single(result, filelist[i],
                             f->ignore_chip_mismatch, err);
    if (err->type != AFFY_ERROR_NONE)
      continue;

    /* Temp chipset now contains the most recently loaded chip */
    temp->chip[0] = result->chip[(result->num_chips) - 1];
    temp->num_chips = 1;

    /* abort on corrupt CEL files, unless --salvage is used */
    if (temp->chip[0]->cel->corrupt_flag &&
        f->salvage_corrupt == false)
          AFFY_HANDLE_ERROR_GOTO("corrupt CEL file",
                                 AFFY_ERROR_BADFORMAT,
                                 err,
                                 cleanup);

    /* Process chip according to various flags */

    if (f->use_background_correction)
    {
      if (f->bg_mas5)
      {
        affy_mas5_background_correction(temp, f, err);
        AFFY_CHECK_ERROR_GOTO(err, cleanup);
      }
      else if (f->bg_rma || f->bg_rma_both)
      {
        if (f->use_mm_probe_subtraction &&
            temp->chip[0]->cdf->no_mm_flag == 0)
        {
          affy_rma_background_correct_pm_mm_together(temp,0,0,err);
        }
        else
        {
          if (f->bg_rma)
          {
            affy_rma_background_correct_pm_mm_together(temp,0,1,err);
          }
          else
          {
            affy_rma_background_correct_pm_mm_together(temp,0,0,err);
          }
        }
      }
      else if (f->bg_iron)
      {
        if (temp->chip[0]->cdf->no_mm_flag == 0)
        {
          affy_rma_background_correct_pm_mm_together(temp,0,0,err);
        }
        else
        {
          affy_rma_background_correct_pm_mm_together(temp,0,0,err);
        }
      }
      else if (f->bg_global)
      {
        affy_global_background_correct(temp,0,err);
      }
    }

    if (f->use_pairwise_normalization)
    {
      info("Performing pairwise probe normalization...");
      affy_pairwise_normalization(temp, 
                                  model_chip, 
                                  AFFY_PAIRWISE_DEFAULT,
                                  f, err);

      AFFY_CHECK_ERROR_GOTO(err, cleanup);

      affy_floor_probe(temp, 1E-5, err);
      AFFY_CHECK_ERROR_GOTO(err, cleanup);

      info("done.\n");
    }

    /* we can only save memory and summarize one at a time if no quantile */
    if (f->use_quantile_normalization == 0)
    {
      /* Make present/absent calls after scaling/normalization takes place */
      if (f->output_present_absent &&
          (f->use_background_correction == 0 || f->bg_mas5 ||
           f->bg_rma_both))
      {
        if (temp->chip[0]->cdf->no_mm_flag == 0)
        {
          affy_mas5_call(temp, f, err);
          AFFY_CHECK_ERROR_GOTO(err, cleanup);
        }
      }

      /* subtract MAS5 MM signals after normalization */
      if (f->use_background_correction && f->use_mm_probe_subtraction)
      {
        if (temp->chip[0]->cdf->no_mm_flag == 0)
        {
          affy_mas5_subtract_mm_signal_probe(temp->chip[0], f, err);
          AFFY_CHECK_ERROR_GOTO(err, cleanup);
        }
      }

      /* renormalize if MM was subtracted... */
      /* ERROR -- model chipset has not been MM subtracted yet!!! */
      if (0 && f->use_background_correction && f->use_mm_probe_subtraction)
      {
        if (f->use_mean_normalization)
        {
          affy_mean_normalization(temp, f->mean_normalization_target_mean, f);
        }
        else if (f->use_pairwise_normalization)
        {
          info("Performing pairwise probe normalization again...");
          affy_pairwise_normalization(temp, 
                                      model_chip, 
                                      AFFY_PAIRWISE_DEFAULT,
                                      f, err);
          AFFY_CHECK_ERROR_GOTO(err, cleanup);

          affy_floor_probe(temp, 1E-5, err);
          AFFY_CHECK_ERROR_GOTO(err, cleanup);

          info("done.\n");
        }
      }

      if (f->use_background_correction && f->bg_iron)
      {
        affy_iron_signal(temp, f, err);
        AFFY_CHECK_ERROR_GOTO(err, cleanup);
      }
      else if (f->use_tukey_biweight)
      {
        affy_mas5_signal(temp, f, err);
        AFFY_CHECK_ERROR_GOTO(err, cleanup);
      }


      /* Free chip space, if we don't need to keep it for printing probes
       * or median polish probeset summarization
       */
      if (f->dump_probe_values == false &&
          f->use_median_polish == false)
      {
        affy_mostly_free_cel_file(temp->chip[0]->cel);
/*      temp->chip[0]->cel = NULL; */
      }
    }

    info("Finished one-at-a-time processing: %s\n", filelist[i]);

    chips_processed++;
  }

  /* Option to use mean normalization */
  /* Must go after all chips are loaded now, so that mean of means can be
   * calculated if target mean = 0.
   */
  if (f->use_normalization && f->use_mean_normalization)
  {
    affy_mean_normalization(result, f->mean_normalization_target_mean, f);
  }

  /* apply postponed quantile normalization */
  if (f->use_quantile_normalization)
  {
    if (f->bg_rma) /* use PM, skipping MM which are now all zeroed */
    {
      affy_quantile_normalization(result, 1, err);
      AFFY_CHECK_ERROR_GOTO(err, cleanup);
    }
    else    /* use both PM and MM */
    {
      affy_quantile_normalization(result, 0, err);
      AFFY_CHECK_ERROR_GOTO(err, cleanup);
    }

    /* Make present/absent calls after scaling/normalization takes place */
    if (f->output_present_absent &&
        (f->use_background_correction == 0 || f->bg_mas5 ||
         f->bg_rma_both))
    {
      if (result->chip[0]->cdf->no_mm_flag == 0)
      {
        affy_mas5_call(result, f, err);
        AFFY_CHECK_ERROR_GOTO(err, cleanup);
      }
    }

    /* subtract MAS5 MM signals after normalization */
    if (f->use_background_correction && f->use_mm_probe_subtraction)
    {
      temp->num_chips = 1;

      for (i = 0; i < max_chips; i++)
      {
        temp->chip[0] = result->chip[i];

        if (temp->chip[0]->cdf->no_mm_flag == 0)
        {
          affy_mas5_subtract_mm_signal_probe(temp->chip[0], f, err);
          AFFY_CHECK_ERROR_GOTO(err, cleanup);
        }
      }
    }
    
    /* renormalize if MM was subtracted... */
    if (0 && f->use_background_correction && f->use_mm_probe_subtraction)
    {
      if (f->bg_rma) /* use PM, skipping MM which are now all zeroed */
      {
        affy_quantile_normalization(result, 1, err);
        AFFY_CHECK_ERROR_GOTO(err, cleanup);
      }
      else    /* use both PM and MM */
      {
        affy_quantile_normalization(result, 0, err);
        AFFY_CHECK_ERROR_GOTO(err, cleanup);
      }
    }

    if (f->use_background_correction && f->bg_iron)
    {
      affy_iron_signal(result, f, err);
      AFFY_CHECK_ERROR_GOTO(err, cleanup);
    }
    else if (f->use_tukey_biweight)
    {
      affy_mas5_signal(result, f, err);
      AFFY_CHECK_ERROR_GOTO(err, cleanup);
    }
  }

  /* apply floors to model chip */
  if (f->use_pairwise_normalization)
    affy_floor_probe(model_chipset, 1E-5, err);

  /* summarize/normalize probesets if we haven't already */
  if (f->use_median_polish)
  {
    /* discard PMs entirely, since we don't need them anymore */
    if (f->dump_probe_values == false)
    {
        for (i = 0; i < max_chips; i++)
            mas5_to_rma_pm_free(result->chip[i], err);

        if (f->use_pairwise_normalization)
            mas5_to_rma_pm_free(model_chip, err);
    }
    else
    {
        for (i = 0; i < max_chips; i++)
            mas5_to_rma_pm(result->chip[i], err);

        if (f->use_pairwise_normalization)
            mas5_to_rma_pm(model_chip, err);
    }

    /* flag to store affinities for later use on model chip */
    if (f->use_pairwise_normalization)
      f->reuse_affinities = true;

    affy_rma_signal(result, f, 0, err);
    
    /* apply median polish trained on other chips to the model chip */
    if (f->use_pairwise_normalization)
    {
      model_chipset->affinities = result->affinities;
      model_chipset->t_values = result->t_values;
      model_chipset->mp_allocated_flag = result->mp_allocated_flag;
      model_chipset->mp_populated_flag = result->mp_populated_flag;
    
      info("Performing probeset summarization on reference sample...");
      affy_rma_signal(model_chipset, f, 0, err);
      AFFY_CHECK_ERROR_GOTO(err, cleanup);
    }
    
#if 1
    /* unlog the data, since the MAS5 routines assume unlogged data */
    affy_unlog_probeset(result, err);
    AFFY_CHECK_ERROR_GOTO(err, cleanup);

    if (f->use_pairwise_normalization)
    {
      affy_unlog_probeset(model_chipset, err);
      AFFY_CHECK_ERROR_GOTO(err, cleanup);
    }
#endif
  }


  if (f->use_tukey_biweight && f->use_pairwise_normalization)
  {
    /* apply postponed MM subtraction to model chipset, calculate probesets */
    if (f->use_background_correction && f->use_mm_probe_subtraction)
    {
        if (model_chip->cdf->no_mm_flag == 0)
        {
          affy_mas5_subtract_mm_signal_probe(model_chip, f, err);
          AFFY_CHECK_ERROR_GOTO(err, cleanup);
        }
    }

    if (f->use_background_correction && f->bg_iron)
    {
      affy_iron_signal(model_chipset, f, err);
      AFFY_CHECK_ERROR_GOTO(err, cleanup);
    }
    else if (f->use_tukey_biweight)
    {
      affy_mas5_signal(model_chipset, f, err);
      AFFY_CHECK_ERROR_GOTO(err, cleanup);
    }
  }

  /* apply median scaling to probesets */
  if (f->use_probeset_scaling &&
      !f->use_quantile_normalization &&
      !f->use_pairwise_normalization)
  {
    affy_mas5_scale(result, f, err);
    AFFY_CHECK_ERROR_GOTO(err, cleanup);
  }

  if (f->use_normalization && f->normalize_probesets)
  {
    /* quantile probeset normalization */
    if (f->use_quantile_normalization)
    {
      info("Performing quantile probeset normalization...");
      affy_quantile_normalization_probeset(result, err);
      AFFY_CHECK_ERROR_GOTO(err, cleanup);
      info("done.\n");
    }

    /* pairwise probeset normalization */
    if (f->use_pairwise_normalization)
    {
      info("Performing pairwise probeset normalization...");
      affy_pairwise_normalization_probeset(result, model_chip, 0,
                                         f, err);
      AFFY_CHECK_ERROR_GOTO(err, cleanup);
      info("done.\n");

      /* Apply scaling factor from model chipset to all chips */
/*
      affy_mas5_scale_iron(result, model_chipset, f, err);
      AFFY_CHECK_ERROR_GOTO(err, cleanup);
*/
    }
  }


  /* floor probeset signal intensity  */
  if (f->bg_iron)
  {
    if (f->use_pairwise_normalization ||
        f->bioconductor_compatability == 0)
    {
      if (f->floor_non_zero_to_one)
        affy_floor_probeset_non_zero_to_one(model_chipset, err);
      else if (f->floor_to_min_non_zero)
        affy_floor_probeset_to_min_non_zero(model_chipset, err);
      else
        affy_floor_probeset(model_chipset, 1.0, err);
    }

    if (f->floor_non_zero_to_one)
      affy_floor_probeset_non_zero_to_one(result, err);
    else if (f->floor_to_min_non_zero)
      affy_floor_probeset_to_min_non_zero(result, err);
    else
      affy_floor_probeset(result, 1.0, err);
  }

  
  if (f->dump_probe_values)
  {
    unsigned int ci;

    affy_write_probe_values(result, f->probe_filename, 0, err);

    for (ci = 0; ci < result->num_chips; ci++)
    {
      affy_mostly_free_cel_file(result->chip[ci]->cel);
/*    result->chip[ci]->cel = NULL; */
    }

    AFFY_CHECK_ERROR_GOTO(err, cleanup);
  }
  
  /* Free up the temp & model chipsets */
  h_free(temp);
  h_free(chip_type);

  info("MAS5/IRON finished on %d samples", chips_processed);

  return (result);

cleanup:
  h_free(chip_type);
  h_free(temp);
  h_free(result);

  return (NULL);
}
