
/**************************************************************************
 *
 * Filename: illumina.c
 *
 * Purpose:  normalize Illumina data
 *
 *
 * Creation: 
 *
 * Author:   Eric A. Welsh
 *
 *
 * Update History
 * --------------
 * 04/06/11: Creation.  Modified from rma.c (EAW)
 * 01/10/13: renamed iron_linear_normalization flag to
 *           iron_global_scaling_normaliztion
 * 03/10/14: #ifdef out CEL qc fields to save memory (EAW)
 * 06/01/18: added affy_load_exclusions_file() (EAW)
 * 09/11/11: pass mempool to affy_load_exclusions_file()
 * 09/14/18: added affy_load_spikeins_file() (EAW)
 * 01/10/20: print column headers for various stderr output (EAW)
 * 09/05/23: change utils_getline() to fgets_strip_realloc() (EAW)
 * 10/05/23: change STDERR GlobalFitLine columns and headers (EAW)
 * 01/10/24: don't free cel data that doesn't exist (EAW)
 * 01/10/24: pass flags to affy_mean_normalization() (EAW)
 * 
 *
 **************************************************************************/

#include "affy_rma.h"
#include "affy_mas5.h"

extern void affy_iron_background_correction_probes(AFFY_CHIPSET *cs,
                                            AFFY_COMBINED_FLAGS *f,
                                            AFFY_ERROR *err);
extern void affy_iron_background_correction_probeset(AFFY_CHIPSET *cs,
                                              AFFY_COMBINED_FLAGS *f,
                                              AFFY_ERROR *err);

int strcmp_insensitive(const char *str1, const char *str2)
{
  char c1, c2;
  char *sptr1 = (char *) str1, *sptr2 = (char *) str2;
  
  assert(str1 != NULL);
  assert(str2 != NULL);
  
  while ((c1 = *sptr1) && (c2 = *sptr2))
  {
    c1 = tolower(c1);
    c2 = tolower(c2);
    
    if (c1 < c2) return -1;
    if (c2 < c1) return  1;
    
    sptr1++;
    sptr2++;
  }
  
  /* one string shorter than the other */
  if (*sptr2) return -1;
  if (*sptr1) return  1;
  
  return 0;
}

static void load_generic_spreadsheet(AFFY_CHIPSET *cs, char *filename,
                                AFFY_ERROR *err)
{
  FILE *data_file;
  AFFY_CHIP    *chip = NULL;
  AFFY_CELFILE *cf   = NULL;
  int max_string_len = 0;
  char *string       = NULL;
  char **fields      = NULL;
  char *sptr;
  int num_fields = 0;
  int max_field  = 0;
  int nbytes;

  int numprobes = 0;
  int ok_chip_flag, ok_probe_flag;
  int i, j;
  
  data_file = fopen(filename, "rb");
  if (!data_file)
    AFFY_HANDLE_ERROR_VOID("can not open data file", AFFY_ERROR_NOTFOUND, err);
  
  /* read header line */
  fgets_strip_realloc(&string, &max_string_len, data_file);
  num_fields = split_tabs(string, &fields, &max_field);
  
  cs->num_chips = 0;
  
  /* allocate chips */
  /* first field is probe, the rest are samples, skip empty columns */
  for (i = 1; i < num_fields; i++)
  {
    ok_chip_flag = 0;

    /* skip empty columns */
    for (sptr = fields[i]; sptr; sptr++)
    {
      if (!isspace(*sptr))
      {
        ok_chip_flag = 1;
        break;
      }
    }
    
    /* allocate chip */
    if (ok_chip_flag)
    {
      chip = h_calloc(1, sizeof(AFFY_CHIP));
      if (chip == NULL)
      {
        AFFY_HANDLE_ERROR_VOID("calloc failed", AFFY_ERROR_OUTOFMEM, err);
      }
      
      /* allocate CEL */
      cf = h_calloc(1, sizeof(AFFY_CELFILE));
      if (cf == NULL)
        AFFY_HANDLE_ERROR_VOID("calloc failed", AFFY_ERROR_OUTOFMEM, err);

      cf->filename = h_strdup(fields[i]);
      if (cf->filename == NULL)
      {
        AFFY_HANDLE_ERROR_VOID("strdup failed", AFFY_ERROR_OUTOFMEM, err);
      }
      hattach(cf->filename, cf);
      
      cf->numrows     = cs->numrows;
      cf->numcols     = cs->numcols;
      cf->nummasks    = 0;
      cf->numoutliers = 0;
      cf->mask        = NULL;
      cf->outlier     = NULL;

      nbytes = numbytes(cf->numcols);

      cf->data = h_subcalloc(cf, cf->numrows, sizeof(AFFY_CELL *));
      if (cf->data == NULL)
        AFFY_HANDLE_ERROR_VOID("calloc failed", AFFY_ERROR_OUTOFMEM, err);

      cf->data[0] = h_subcalloc(cf->data, 
                                cf->numcols * cf->numrows, 
                                sizeof(AFFY_CELL));
      if (cf->data[0] == NULL)
        AFFY_HANDLE_ERROR_VOID("calloc failed", AFFY_ERROR_OUTOFMEM, err);

      cf->mask = h_subcalloc(cf, cf->numrows, sizeof(bitstr_t *));
      if (cf->mask == NULL)
        AFFY_HANDLE_ERROR_VOID("calloc failed", AFFY_ERROR_OUTOFMEM, err);

      cf->mask[0] = h_subcalloc(cf->mask, nbytes * cf->numrows, sizeof(bitstr_t));
      if (cf->mask[0] == NULL)
        AFFY_HANDLE_ERROR_VOID("calloc failed", AFFY_ERROR_OUTOFMEM, err);

      cf->outlier = h_subcalloc(cf, cf->numrows, sizeof(bitstr_t *));
      if (cf->outlier == NULL)
        AFFY_HANDLE_ERROR_VOID("calloc failed", AFFY_ERROR_OUTOFMEM, err);

      cf->outlier[0] = h_subcalloc(cf->outlier, 
                                   nbytes * cf->numrows, 
                                   sizeof(bitstr_t));
      if (cf->outlier[0] == NULL)
        AFFY_HANDLE_ERROR_VOID("calloc failed", AFFY_ERROR_OUTOFMEM, err);
      
      for (j = 1; j < cf->numrows; j++)
      {
        cf->data[j]    = cf->data[0] + (j * cf->numcols);
        cf->mask[j]    = cf->mask[0] + (j * nbytes);
        cf->outlier[j] = cf->outlier[0] + (j * nbytes);
      }

      chip->cdf       = NULL;
      chip->cel       = cf;
      chip->dat       = NULL;
      chip->filename  = NULL;
      chip->probe_set = NULL;
      chip->pm        = NULL;

      hattach(chip->cel, chip);

      chip->filename = h_strdup(fields[i]);
      if (chip->filename == NULL)
      {
        affy_free_chip(chip);

        AFFY_HANDLE_ERROR_VOID("strdup failed", AFFY_ERROR_OUTOFMEM, err);
      }
      hattach(chip->filename, chip);

      cs->chip[cs->num_chips] = chip;
      hattach(cs->chip[cs->num_chips], cs->chip);
      cs->chip[cs->num_chips]->cdf = cs->cdf;
      cs->num_chips++;
    }
  }

  /* assume only a single line per probe */
  /* multiple identical probes will be treated as separate probes */
  /* skip blank probes */
  while(fgets_strip_realloc(&string, &max_string_len, data_file))
  {
    num_fields = split_tabs(string, &fields, &max_field);

    ok_probe_flag = 0;

    /* skip empty columns */
    for (sptr = fields[0]; sptr; sptr++)
    {
      if (!isspace(*sptr))
      {
        ok_probe_flag = 1;
        break;
      }
    }
    
    if (ok_probe_flag)
    {
      /* probeset name */
      cs->cdf->probeset[numprobes].name = h_strdup(fields[0]);
      if (cs->cdf->probeset[numprobes].name == NULL)
      {
        AFFY_HANDLE_ERROR_VOID("strdup failed", AFFY_ERROR_OUTOFMEM, err);
      }
      hattach(cs->cdf->probeset[numprobes].name, cs->cdf);


      cs->num_chips = 0;

      for (i = 1; i < num_fields; i++)
      {
        ok_chip_flag = 0;

        /* skip empty columns */
        for (sptr = fields[i]; sptr; sptr++)
        {
          if (!isspace(*sptr))
          {
            ok_chip_flag = 1;
            break;
          }
        }

        /* store intensity data */
        if (ok_chip_flag)
        {
          chip = cs->chip[cs->num_chips];
          cf   = chip->cel;

          cf->data[0][numprobes].value     = atof(fields[i]);
#ifdef STORE_CEL_QC
          cf->data[0][numprobes].stddev    = 0;
          cf->data[0][numprobes].numpixels = 0;
          cf->data[0][numprobes].pixels    = NULL;
#endif
          
          cs->num_chips++;
        }
      }

      numprobes++;
    }
  }
  
  fclose(data_file);

  if (string)
    free(string);
  if (fields)
    free(fields);
}

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

  if (cp->cel->data)
    h_free(cp->cel->data);
  cp->cel->data = NULL;
}


/* assume numprobes == numprobesets */
int fill_probesets_with_probes(AFFY_CHIPSET *cs,
                               AFFY_ERROR *err)
{
  affy_int32    x, y, numprobes, numprobesets;
  affy_uint32  i, p;
  AFFY_CDFFILE *cdf;
  AFFY_CELL   **chip_data;

  assert(cs              != NULL);
  assert(cs->cdf         != NULL);
  assert(cs->chip        != NULL);

  if (cs->num_chips == 0)
    return -1;

  cdf       = cs->cdf;
  numprobes = cdf->numprobes;
  numprobesets = cdf->numprobesets;

  for (i = 0; i < cs->num_chips; i++)
  {
    assert(numprobes == numprobesets);

    cs->chip[i]->probe_set = h_subcalloc(cs->chip[i], numprobesets,
                                         sizeof(double));
    if (cs->chip[i]->probe_set == NULL)
      AFFY_HANDLE_ERROR_GOTO("calloc failed", AFFY_ERROR_OUTOFMEM, err, err);

    cs->chip[i]->numprobesets = numprobesets;

    chip_data = cs->chip[i]->cel->data;

    if (cs->chip[i]->pm)
    {
      for (p = 0; p < numprobes; p++)
      {
        cs->chip[i]->probe_set[p] = cs->chip[i]->pm[p];
      }
    }
    if (cdf->probe && chip_data)
    {
      for (p = 0; p < numprobes; p++)
      {
        x = cdf->probe[p]->pm.x;
        y = cdf->probe[p]->pm.y;

        cs->chip[i]->probe_set[p] = chip_data[x][y].value;
      }
    }
  }
  
  return (0);

err:
  return (-1);
}


AFFY_CHIPSET *affy_illumina(char **filelist, AFFY_COMBINED_FLAGS *f,
                       AFFY_ERROR *err)
{
  AFFY_CHIPSET         *result = NULL, *model_chipset = NULL, *temp = NULL;
  AFFY_CHIP            *model_chip = NULL;
  AFFY_COMBINED_FLAGS  default_flags;
  int                  i;
  affy_uint32          max_chips, numprobes;
  int                  *mempool;
  double               *mean = NULL;
  char                 **sample_names = NULL;
  int                  model_chip_idx;

  assert(filelist != NULL);

  /* In case flags not used, create default entry */
  if (f == NULL)
  {
    affy_mas5_set_defaults(&default_flags);
    affy_rma_set_defaults(&default_flags);
    f = &default_flags;
  }
  

  mempool = h_malloc(sizeof(int));
  if (mempool == NULL)
    AFFY_HANDLE_ERROR("malloc failed", AFFY_ERROR_OUTOFMEM, err, NULL);

#if 0
  /* Get the array type. */
  chip_type = affy_get_cdf_name_from_cel(filelist[0], err);
  if (chip_type == NULL)
    AFFY_CHECK_ERROR_GOTO(err, cleanup);

  hattach(chip_type, mempool);
#endif

  /* scan for max probes and max chips */
  /* only read data in from the first file, ignore all others */
  get_generic_spreadsheet_bounds(filelist[0], &numprobes, &max_chips, err);
  AFFY_CHECK_ERROR_GOTO(err, cleanup);
  
  info("NumSamples:\t%d\tNumProbes:\t%d", max_chips, numprobes);


  /* Create structure */
  result = create_blank_generic_chipset(max_chips, numprobes, err);
  AFFY_CHECK_ERROR_GOTO(err, cleanup);

  /* Read in data */
  load_generic_spreadsheet(result, filelist[0], err);
  AFFY_CHECK_ERROR_GOTO(err, cleanup);

  /* allocate chip names */
  sample_names = h_subcalloc(mempool, max_chips, sizeof(char **));
  if (sample_names == NULL)
    AFFY_HANDLE_ERROR_GOTO("calloc failed",
                           AFFY_ERROR_OUTOFMEM,
                           err,
                           cleanup);
  
  /* store chip names */
  for (i = 0; i < max_chips; i++)
  {
    sample_names[i] = h_strdup(result->chip[i]->filename);
    if (sample_names[i] == NULL)
    {
      AFFY_HANDLE_ERROR("strdup failed", AFFY_ERROR_OUTOFMEM, err, NULL);
    }
    
    hattach(sample_names[i], mempool);
  }

  /* allocate temp chipset */
  temp = create_blank_generic_chipset(1, numprobes, err);
  AFFY_CHECK_ERROR_GOTO(err, cleanup);

  /* sanity check for affinity reuse flag */
  if (f->use_rma_probeset_singletons)
    f->reuse_affinities = false;
  if (f->use_saved_affinities)
    f->reuse_affinities = false;
  
  /* load in probesets to exclude from IRON training */
  if (f->use_exclusions)
    affy_load_exclusions_file(f->exclusions_filename, result->cdf,
                              mempool, err);
  if (f->use_spikeins)
    affy_load_spikeins_file(f->spikeins_filename, result->cdf,
                            mempool, err);
  
  if (f->use_pairwise_normalization)
  {
    info("Loading pairwise normalization model from %s",
         f->pairwise_model_filename);

    /* find model chip */
    model_chip_idx = -1;
    for (i = 0; i < max_chips; i++)
    {
      if (strcmp_insensitive(f->pairwise_model_filename,
                             result->chip[i]->filename) == 0)
      {
        model_chip_idx = i;
        break;
      }
    }

    if (model_chip_idx == -1)
    {
      AFFY_HANDLE_ERROR("can not find pairwise reference sample",
                        AFFY_ERROR_UNKNOWN, err, NULL);
    }

    /* clone model chipset */
    model_chipset = affy_clone_chipset_one_chip(result, model_chip_idx, err);
    AFFY_CHECK_ERROR_GOTO(err, cleanup);
    hattach(model_chipset, mempool);

    /* clone model chip */
    model_chipset->chip[0] = affy_clone_chip(result->chip[model_chip_idx],
                                             err);
    hattach(model_chipset->chip[0], model_chipset);

    model_chip             = model_chipset->chip[0];

#if 0
    /* perform a pass of normalization *before* background subtraction */
    if (f->use_normalization && f->use_pairwise_normalization)
    {
      affy_pairwise_normalization(result,
                                  model_chip,
                                  AFFY_PAIRWISE_DEFAULT,
                                  f, err);
    }
#endif

    if (f->use_background_correction)
    {
      if (f->bg_mas5)
      {
        affy_mas5_background_correction(model_chipset, f, err);
        AFFY_CHECK_ERROR_GOTO(err, cleanup);

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
      else if (f->bg_global)
      {
        /* Load PM data, then free all data */
        load_pm(model_chipset->chip[0], err);
        AFFY_CHECK_ERROR_GOTO(err, cleanup);

        affy_global_background_correct_pm_only(model_chipset, 0, err);
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
    result->num_chips = i+1;
    cur_chip = result->num_chips - 1;

    /* Background correct, before normalization */
    if (f->use_background_correction)
    {
      if (f->bg_mas5)
      {
        temp->chip[0] = result->chip[cur_chip];
        temp->num_chips = 1;

        affy_mas5_background_correction(temp, f, err);
        AFFY_CHECK_ERROR_GOTO(err, cleanup);

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
      else if (f->bg_global)
      {
        /* Load PM data, then free all data */
        load_pm(result->chip[cur_chip], err);
        AFFY_CHECK_ERROR_GOTO(err, cleanup);

        affy_global_background_correct_pm_only(result, cur_chip, err);
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
      if (!f->use_mean_normalization && !f->use_pairwise_normalization)
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

  /* Option to use mean normalization */
  /* Must go after all chips are loaded now, so that mean of means can be
   * calculated if target mean = 0.
   */
  if (f->use_normalization && f->use_mean_normalization)
  {
    affy_mean_normalization(result, f->mean_normalization_target_mean, f);
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
#if 0
  /* don't normalize the probes, normalize the probesets later */
  /* --iron-condense-training was changed to only work on probesets */
  else if (f->use_normalization && f->use_pairwise_normalization)
  {
    info("Performing pairwise probe normalization...");

    if (f->iron_global_scaling_normalization)
    {
        fprintf(stderr, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                        "GlobalScale:",
                        "SampleID",
                        "Scale",
                        "Log2Scale",
                        "TrainingSet",
                        "PresentBoth",
                        "PresentSample",
                        "PresentDataset",
                        "FractionTrain");
    }
    else if (f->iron_untilt_normalization)
    {
        fprintf(stderr, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                        "GlobalFitLine:",
                        "SampleID",
                        "Scale",
                        "Log2Scale",
                        "UnTiltDegrees",
                        "TrainingSet",
                        "PresentBoth",
                        "PresentSample",
                        "PresentDataset",
                        "FractionTrain");
    }
    
    affy_pairwise_normalization(result,
                                model_chip,
                                AFFY_PAIRWISE_PM_ONLY,
                                f, err);
    
    AFFY_CHECK_ERROR_GOTO(err, cleanup);
    info("done.\n");
  }
#endif

  /* Dump out the raw pm values, if desired. */
  if (f->dump_probe_values)
  {
    affy_write_probe_values(result, f->probe_filename, AFFY_USE_PM, err);
    AFFY_CHECK_ERROR_GOTO(err, cleanup);
  }

  fill_probesets_with_probes(result, err);

#if 1
  if (f->use_normalization && f->use_pairwise_normalization)
  {
    fill_probesets_with_probes(model_chipset, err);

    info("Performing pairwise probeset normalization...");

    if (f->iron_global_scaling_normalization)
    {
        fprintf(stderr, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                        "GlobalScale:",
                        "SampleID",
                        "Scale",
                        "Log2Scale",
                        "TrainingSet",
                        "PresentBoth",
                        "PresentSample",
                        "PresentDataset",
                        "FractionTrain");
    }
    else if (f->iron_untilt_normalization)
    {
        fprintf(stderr, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                        "GlobalFitLine:",
                        "SampleID",
                        "Scale",
                        "Log2Scale",
                        "UnTiltDegrees",
                        "TrainingSet",
                        "PresentBoth",
                        "PresentSample",
                        "PresentDataset",
                        "FractionTrain");
    }

    affy_pairwise_normalization_probeset(result, model_chip, 0, f, err);
    AFFY_CHECK_ERROR_GOTO(err, cleanup);
    info("done.\n");
  }
#endif

  info("IRON processing finished on %u samples", result->num_chips);
  
  /* Free up remaining space */
  for (i = 0; i < result->num_chips; i++)
  {
    affy_free_cel_file(result->chip[i]->cel);
    result->chip[i]->cel = NULL;
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
