
/**************************************************************************
 *
 * Filename: rma.c
 *
 * Purpose:  Top-level function to perform RMA; ties together all the
 *           subsidiary processing phases.
 *
 * Creation: 
 *
 * Author:   Steven Eschrich
 *
 *
 * Update History
 * --------------
 * 04/14/05: Imported/repaired from old libaffy (AMH)
 * 04/19/05: AFFY_CELFILE now uses AFFY_CELL's (AMH)
 * 12/11/06: Added various internal assertions (AMH)
 * 06/05/08: New error handling and chipset API (AMH)
 * 09/20/10: Pooled memory allocator (AMH)
 * 10/18/10: Pairwise normalization (AMH)
 * 10/18/10: Do not sort filelist alphabetically (EAW)
 * 10/20/10: Added 2nd pairwise normalization pass at probeset level (EAW)
 * 10/22/10: Use new AFFY_COMBINED_FLAGS instead of AFFY_RMA_FLAGS (EAW)
 * 03/10/11: Mean normalization works again, probably broken years ago? (EAW)
 * 05/10/13: replaced affy_free_cel_file() with
 *           affy_mostly_free_cel_file(), to allow for corrupt chip
 *           warnings at program end if --salvage was used (EAW)
 * 05/22/19: --ignore-chip-mismatch support (EAW)
 * 08/12/20: pass flags to affy_create_chipset() (EAW)
 * 09/05/23: change utils_getline() to fgets_strip_realloc() (EAW)
 *
 **************************************************************************/

#include "affy_rma.h"
#include "affy_mas5.h"

extern void affy_iron_background_correction_probeset(AFFY_CHIPSET *cs,
                                              AFFY_COMBINED_FLAGS *f,
                                              AFFY_ERROR *err);

static void load_pm(AFFY_CHIP *cp, AFFY_ERROR *err)
{
  affy_int32 numprobes, k;

  assert(cp            != NULL);
  assert(cp->cdf       != NULL);
  assert(cp->cel       != NULL);
  assert(cp->cel->data != NULL);

  numprobes = cp->cdf->numprobes;

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

AFFY_CHIPSET *affy_rma(char **filelist, AFFY_COMBINED_FLAGS *f,
                       AFFY_ERROR *err)
{
  AFFY_CHIPSET         *result = NULL, *model_chipset = NULL, *temp = NULL;
  AFFY_CHIP            *model_chip = NULL;
  AFFY_COMBINED_FLAGS  default_flags;
  int                  i;
  char                 *chip_type, **p;
  int                  max_chips, numprobes, *mempool = NULL;
  double               *mean = NULL;
  int                  safe_to_write_affinities_flag = 0;

  assert(filelist != NULL);

  /* In case flags not used, create default entry */
  if (f == NULL)
  {
    affy_mas5_set_defaults(&default_flags);
    affy_rma_set_defaults(&default_flags);
    f = &default_flags;
  }
  
  /* sanity check, can not mix-and-match quantile with post-norm non-rma bg */
  if (f->use_background_correction && f->normalize_before_bg &&
      !f->bg_rma)
  {
    AFFY_HANDLE_ERROR_GOTO("Normalization before BG unsupported with selected bg method",
                           AFFY_ERROR_NOTSUPP,
                           err, 
                           cleanup);
  }
  if (f->use_background_correction && f->normalize_before_bg &&
      !f->use_pairwise_normalization)
  {
    warn("WARNING - non-IRON normalization before BG may yield odd results\n");
  }


  mempool = h_malloc(sizeof(int));
  if (mempool == NULL)
    AFFY_HANDLE_ERROR("malloc failed", AFFY_ERROR_OUTOFMEM, err, NULL);

  /* Get the array type. */
  chip_type = affy_get_cdf_name_from_cel(filelist[0], err);
  if (chip_type == NULL)
    AFFY_CHECK_ERROR_GOTO(err, cleanup);

  hattach(chip_type, mempool);

  /* Count up chips */
  for (p = filelist, max_chips = 0; *p != NULL; p++)
    max_chips++;

  /* Create structure */
  result = affy_create_chipset(max_chips, chip_type, f->cdf_directory, f, err);
  AFFY_CHECK_ERROR_GOTO(err, cleanup);

  temp = affy_clone_chipset(result, err);
  AFFY_CHECK_ERROR_GOTO(err, cleanup);

  /* abort if any MM probes are missing */
  if (result->cdf->dupe_probes_flag)
    AFFY_HANDLE_ERROR("multiple probesets share same probe, use 'iron --norm-quantile --median-polish' instead", AFFY_ERROR_NOTSUPP, err, NULL);

  numprobes = result->cdf->numprobes;
  
  /* sanity checks for affinity reuse flag */
  if (f->use_rma_probeset_singletons)
    f->reuse_affinities = false;
  if (f->use_saved_affinities)
    f->reuse_affinities = false;
  
  if (f->use_pairwise_normalization)
  {
    info("Loading pairwise normalization model from %s",
         f->pairwise_model_filename);
    
    model_chipset = affy_clone_chipset(result, err);
    AFFY_CHECK_ERROR_GOTO(err, cleanup);

    hattach(model_chipset, mempool);

    affy_load_chipset_single(model_chipset, f->pairwise_model_filename,
                             f->ignore_chip_mismatch, err);
    AFFY_CHECK_ERROR_GOTO(err, cleanup);

    model_chip = model_chipset->chip[0];

    /* abort on corrupt CEL files, unless --salvage is used */
    if (model_chip->cel->corrupt_flag &&
        f->salvage_corrupt == false)
          AFFY_HANDLE_ERROR_GOTO("corrupt CEL file",
                                 AFFY_ERROR_BADFORMAT,
                                 err,
                                 cleanup);

    if (f->use_background_correction && !f->normalize_before_bg)
    {
      if (f->bg_mas5)
      {
        affy_mas5_background_correction(model_chipset, f, err);
        AFFY_CHECK_ERROR_GOTO(err, cleanup);

        if (f->use_mm_probe_subtraction)
        {
          affy_mas5_subtract_mm_signal_probe(model_chip, f, err);
          AFFY_CHECK_ERROR_GOTO(err, cleanup);
        }

        /* Load PM data, then free all data */
        load_pm(model_chipset->chip[0], err);
        AFFY_CHECK_ERROR_GOTO(err, cleanup);
      }
      else if (f->bg_rma)
      {
        /* Load PM data, then free all data */
        load_pm(model_chipset->chip[0], err);
        AFFY_CHECK_ERROR_GOTO(err, cleanup);

        affy_rma_background_correct(model_chipset, 0, err);
      }
    }
    else
    {
      /* Load PM data, then free all data */
      load_pm(model_chipset->chip[0], err);
      AFFY_CHECK_ERROR_GOTO(err, cleanup);
    }

    info("Pairwise reference sample loaded");
  }

  /* Load each chip */
  for (i = 0; i < max_chips; i++)
  {
    int cur_chip;

    /* Load chip */
    affy_load_chipset_single(result, filelist[i],
                             f->ignore_chip_mismatch, err);
    if (err->type != AFFY_ERROR_NONE)
      continue;

    cur_chip = result->num_chips - 1;

    /* abort on corrupt CEL files, unless --salvage is used */
    if (result->chip[cur_chip]->cel->corrupt_flag &&
        f->salvage_corrupt == false)
          AFFY_HANDLE_ERROR_GOTO("corrupt CEL file",
                                 AFFY_ERROR_BADFORMAT,
                                 err,
                                 cleanup);

    /* Background correct, before normalization */
    if (f->use_background_correction && !f->normalize_before_bg)
    {
      if (f->bg_mas5)
      {
        temp->chip[0] = result->chip[cur_chip];
        temp->num_chips = 1;

        affy_mas5_background_correction(temp, f, err);
        AFFY_CHECK_ERROR_GOTO(err, cleanup);

        if (f->use_mm_probe_subtraction)
        {
          affy_mas5_subtract_mm_signal_probe(temp->chip[0], f, err);
          AFFY_CHECK_ERROR_GOTO(err, cleanup);
        }

        /* Load PM data, then free all data */
        load_pm(result->chip[cur_chip], err);
        AFFY_CHECK_ERROR_GOTO(err, cleanup);
      }
      else if (f->bg_rma)
      {
        /* Load PM data, then free all data */
        load_pm(result->chip[cur_chip], err);
        AFFY_CHECK_ERROR_GOTO(err, cleanup);

        affy_rma_background_correct(result, cur_chip, err);
      }
    }
    else
    {
      /* Load PM data, then free all data */
      load_pm(result->chip[cur_chip], err);
      AFFY_CHECK_ERROR_GOTO(err, cleanup);
    }

    /* Normalize (partially) */
    if (f->use_normalization)
    {
      /* Option to use mean normalization */
      if (f->use_mean_normalization)
      {
        affy_mean_normalization(result, f->mean_normalization_target_mean);
      }
      else if (!f->use_pairwise_normalization)
      {
        /* Default is quantile normalization */
        if (mean == NULL)
        {
          mean = h_subcalloc(mempool, numprobes, sizeof(double));
          if (mean == NULL)
            AFFY_HANDLE_ERROR_GOTO("calloc failed",
                                   AFFY_ERROR_OUTOFMEM,
                                   err,
                                   cleanup);
        }

        affy_rma_quantile_normalization_chip(result, cur_chip, mean, f, err);
        AFFY_CHECK_ERROR_GOTO(err, cleanup);
      }
    }
  }

  /* Redistribute quantile means */
  if (f->use_normalization && !f->use_mean_normalization 
      && !f->use_pairwise_normalization)
  {
    if (f->use_saved_means)
    {
      FILE *fp;
      char *nl, *err_str;
      int   i = 0;
      int   max_string_len;
      
      nl             = NULL;
      max_string_len = 0;

      fp = fopen(f->means_filename, "rb");
      if (fp == NULL)
        AFFY_HANDLE_ERROR_GOTO("couldn't open saved means file",
                               AFFY_ERROR_NOTFOUND,
                               err,
                               cleanup);

      /* while ((nl = utils_getline(fp)) != NULL) */
      while (fgets_strip_realloc(&nl, &max_string_len, fp) != NULL)
      {
	mean[i] = strtod(nl, &err_str);
	
	if ((nl == err_str) && (mean[i] == 0))
        {
          if (nl) free(nl);
        
	  warn("error parsing mean value from %s, line %d\n", 
               f->means_filename, 
               i);
          AFFY_HANDLE_ERROR_GOTO("error parsing mean value", 
                                 AFFY_ERROR_BADFORMAT, 
                                 err, 
                                 cleanup);
	}

	i++;
      }
      fclose(fp);

      if (nl) free(nl);

      if (i != numprobes)
      {
	warn("expected %d means, found %d\n", numprobes, i);
        AFFY_HANDLE_ERROR_GOTO("incorrect number of saved means",
                               AFFY_ERROR_BADFORMAT,
                               err,
                               cleanup);
      }
    }
    else
    {
      for (i = 0; i < numprobes; i++)
	mean[i] /= result->num_chips;
    }
  }
  
  /* Save the means if such was requested. */
  if (f->dump_expression_means)
  {
    FILE *fp;
    int   i;
    char *mf_name = f->means_filename;
    
    assert(mean != NULL);
    
    fp = fopen(mf_name, "w");
    if (fp == NULL)
      AFFY_HANDLE_ERROR_GOTO("couldn't open means file for writing",
                             AFFY_ERROR_IO,
                             err,
                             cleanup);
    
    for (i = 0; i < numprobes; i++)
      fprintf(fp, "%.15e\n", mean[i]);
    
    fclose(fp);
  }

  /* XXX ask Steven about this. */
  if (f->use_normalization && !f->use_pairwise_normalization &&
      !f->use_mean_normalization)
  {
    affy_rma_quantile_normalization_chipset(result, mean, f);
  }
  else if (f->use_normalization && f->use_pairwise_normalization)
  {
    info("Performing pairwise probe normalization...");
    affy_pairwise_normalization(result, 
                                model_chip, 
                                AFFY_PAIRWISE_PM_ONLY,
                                f, err);
    AFFY_CHECK_ERROR_GOTO(err, cleanup);

    affy_floor_probe(result, 1E-5, err);
    AFFY_CHECK_ERROR_GOTO(err, cleanup);

    info("done.\n");
  }


  /* Background correct, after normalization */
  /* All non-RMA methods currently unsupported, due to PM-only stuff */
  /* This will only get fixed when I merge the MAS5 and RMA programs */
  if (f->use_background_correction && f->normalize_before_bg)
  {
    if (f->use_pairwise_normalization)
      affy_rma_background_correct(model_chipset, 0, err);

    for (i = 0; i < max_chips; i++)
    {
      int cur_chip = i;
    
      if (f->bg_rma)
      {
        affy_rma_background_correct(result, cur_chip, err);
      }
    }

    if (f->use_normalization && f->use_pairwise_normalization)
    {
      info("Performing 2nd pass post-BG pairwise probe normalization...");
      affy_pairwise_normalization(result, 
                                  model_chip, 
                                  AFFY_PAIRWISE_PM_ONLY,
                                  f, err);
      AFFY_CHECK_ERROR_GOTO(err, cleanup);

      affy_floor_probe(result, 1E-5, err);
      AFFY_CHECK_ERROR_GOTO(err, cleanup);

      info("done.\n");
    }
  }

  if (f->use_pairwise_normalization)
  {
    /* apply floors to model chip */
    affy_floor_probe(model_chipset, 1E-5, err);

    if (f->use_rma_probeset_singletons)
    {
      info("Performing probeset summarization on reference sample...");
      affy_rma_signal(model_chipset, f, safe_to_write_affinities_flag, err);
      AFFY_CHECK_ERROR_GOTO(err, cleanup);
    }
  }

  
  /* Dump out the raw pm values, if desired. */
  if (f->dump_probe_values)
  {
    affy_write_probe_values(result, f->probe_filename, AFFY_USE_PM, err);
    AFFY_CHECK_ERROR_GOTO(err, cleanup);
  }


  /* Actually calculate the expression */

  /* use each chip by itself when estimating affinities */
  if (f->use_rma_probeset_singletons)
  {
    /*
     * The idea here is to load each chip into the result chipset,
     * copy its ptr to the temporary "singleton" chipset which is then
     * passed to the various processing routines (to operate on one
     * chip at a time).  We're left with the finished product in
     * `result'.
     */
     for (i = 0; i < max_chips; i++)
     {
       temp->chip[0] = result->chip[i];
       temp->num_chips = 1;

       affy_rma_signal(temp, f, safe_to_write_affinities_flag, err);
       AFFY_CHECK_ERROR_GOTO(err, cleanup);
     }
  }
  /* use all the chips when estimating affinities */
  else
  {
    safe_to_write_affinities_flag = 1;

    affy_rma_signal(result, f, safe_to_write_affinities_flag, err);
    AFFY_CHECK_ERROR_GOTO(err, cleanup);
    
    safe_to_write_affinities_flag = 0;
  }

  if (f->use_normalization && f->use_pairwise_normalization &&
      !f->use_rma_probeset_singletons)
  {
    model_chipset->affinities = result->affinities;
    model_chipset->t_values = result->t_values;
    model_chipset->mp_allocated_flag = result->mp_allocated_flag;
    model_chipset->mp_populated_flag = result->mp_populated_flag;
    
    info("Performing probeset summarization on reference sample...");
    affy_rma_signal(model_chipset, f, safe_to_write_affinities_flag, err);
    AFFY_CHECK_ERROR_GOTO(err, cleanup);
  }

  if (f->use_background_correction && f->bg_iron)
  {
    affy_iron_background_correction_probeset(model_chipset, f, err);
    AFFY_CHECK_ERROR_GOTO(err, cleanup);

    affy_iron_background_correction_probeset(result, f, err);
    AFFY_CHECK_ERROR_GOTO(err, cleanup);
  }
  
#if 1
  if (f->use_normalization && f->use_pairwise_normalization)
  {
    info("Performing pairwise probeset normalization...");
    affy_pairwise_normalization_probeset(result, model_chip, 1, f, err);
    AFFY_CHECK_ERROR_GOTO(err, cleanup);
    info("done.\n");
  }
#endif

  /* floor to signal intensity of 1 (0 in log2 space) */
  if (!f->bioconductor_compatability)
    affy_floor_probeset(result, 0.0, err);

  info("RMA finished on %u samples", result->num_chips);
  
  /* Free up remaining space */
  for (i = 0; i < result->num_chips; i++)
  {
    affy_mostly_free_cel_file(result->chip[i]->cel);
/*  result->chip[i]->cel = NULL; */
  }
  
  h_free(temp);
  h_free(mempool);

  return (result);

cleanup:
  h_free(temp);
  h_free(mempool);
  affy_free_chipset(result);

  return (NULL);
}
