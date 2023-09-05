
/**************************************************************************
 *
 * Filename: rma_signal.c
 *
 * Purpose:  Calculate the final expression values for a chipset.
 *
 * Creation: 
 *
 * Author:   Steven Eschrich
 *
 *
 * Update History
 * --------------
 * 04/14/05: Imported/repaired from old libaffy (AMH)
 * 04/12/07: Add support for incremental RMA (AMH)
 * 06/05/08: New error handling scheme (AMH)
 * 09/20/10: Pooled memory allocator (AMH)
 * 10/22/10: Use new AFFY_COMBINED_FLAGS instead of AFFY_RMA_FLAGS
 * 11/16/10: Added ability to reuse calculated median polish affinities
 * 11/19/10: Pass flags for bioconductor compatability (EAW)
 * 09/05/23: change utils_getline() to fgets_strip_realloc() (EAW)
 *
 **************************************************************************/

#include <affy_rma.h>

/* XXX remove strtok! */
static void read_affinities(FILE *infile, 
                            double *affinities, 
                            double *t,
			    unsigned int numprobes,
                            AFFY_ERROR *err)
{
  char *line = NULL, *val  = NULL;
  int   i;
  int   max_string_len = 0;

  assert(infile     != NULL);
  assert(affinities != NULL);
  assert(t          != NULL);

  /* Get the t-value. */
  /* if ((line = utils_getline(infile)) == NULL) */
  if (fgets_strip_realloc(&line, &max_string_len, infile) == NULL)
    AFFY_HANDLE_ERROR_GOTO("failed to parse affinity value from file",
                           AFFY_ERROR_BADFORMAT,
                           err,
                           cleanup);

  if (   ((val = strtok(line, " ")) == NULL)
      || ((val = strtok(NULL, " ")) == NULL))
    AFFY_HANDLE_ERROR_GOTO("failed to parse affinity value from file",
                           AFFY_ERROR_BADFORMAT,
                           err,
                           cleanup);

  if (parsefloat(val, t) == -1)
    AFFY_HANDLE_ERROR_GOTO("failed to parse affinity value from file",
                           AFFY_ERROR_BADFORMAT,
                           err,
                           cleanup);

  free(line);

  /* Now read an affinity for each probe. */
  for (i = 0; i < numprobes; i++)
  {
    int j;
    
    /* if ((line = utils_getline(infile)) == NULL) */
    if (fgets_strip_realloc(&line, &max_string_len, infile) == NULL)
      AFFY_HANDLE_ERROR_GOTO("failed to parse affinity value from file",
                             AFFY_ERROR_BADFORMAT,
                             err,
                             cleanup);
    
    /* Skip irrelevant stuff. */
    for (j = 0, val = strtok(line, " "); 
	 val && (j < 3); 
	 j++, val = strtok(NULL, " "))
      ;

    if ((j != 3) || (!val))
      AFFY_HANDLE_ERROR_GOTO("failed to parse affinity value from file",
                             AFFY_ERROR_BADFORMAT,
                             err,
                             cleanup);

    if (parsefloat(val, affinities+i) == -1)
      AFFY_HANDLE_ERROR_GOTO("failed to parse affinity value from file",
                             AFFY_ERROR_BADFORMAT,
                             err,
                             cleanup);
    
    free(line);
  }

  return;

cleanup:
  if (line) free(line);
}

/*
 * This is the PM-only array implementation of rma median polish:
 * useful for low-memory applications.
 */
void affy_rma_signal(AFFY_CHIPSET *c, AFFY_COMBINED_FLAGS *f,
                     int safe_to_write_affinities_flag, AFFY_ERROR *err)
{
  int               i, j, ps, *mempool, affinity_size = 0;
  unsigned int      numchips;
  double           *results = NULL, *affinities = NULL, t, **z = NULL;
  double            LOG2 = log(2.0);
  AFFY_CDFFILE     *cdf;
  FILE             *aff_file = NULL;
  LIBUTILS_PB_STATE pbs;

  /* Preconditions: The CHIPSET and the CDF exist */
  assert(c      != NULL);
  assert(c->cdf != NULL);
  assert(f      != NULL);

  /* Convenience variable */
  cdf      = c->cdf;
  numchips = c->num_chips;

  pb_init(&pbs);
  pb_begin(&pbs, cdf->numprobesets, "Calculating expressions");

  /* Allocate storage */
  mempool = h_malloc(sizeof(int));
  if (mempool == NULL)
    AFFY_HANDLE_ERROR_VOID("malloc failed", AFFY_ERROR_OUTOFMEM, err);

  for (i = 0; i < numchips; i++)
  {
    c->chip[i]->probe_set = h_subcalloc(c->chip[i], 
                                        cdf->numprobesets, 
                                        sizeof(double));
    if (c->chip[i]->probe_set == NULL)
      AFFY_HANDLE_ERROR_GOTO("calloc failed", 
                             AFFY_ERROR_OUTOFMEM, 
                             err, 
                             cleanup);

    c->chip[i]->numprobesets = cdf->numprobesets;
  }
  
  /* allocate arrays for reuse of median polish affinities */
  if (f->reuse_affinities &&
      !f->use_saved_affinities &&
      !c->mp_allocated_flag)
  {
    c->affinities = h_subcalloc(c,
                                cdf->numprobesets,
                                sizeof(double *));
    c->t_values   = h_subcalloc(c,
                                cdf->numprobesets,
                                sizeof(double *));
  }

  results = h_subcalloc(mempool, numchips, sizeof(double));
  if (results == NULL)
    AFFY_HANDLE_ERROR_GOTO("calloc failed", AFFY_ERROR_OUTOFMEM, err, cleanup);

  /* If we're using saved affinity values or dumping calculated ones,
     get the file opened first. */
  if (f->dump_probe_affinities && safe_to_write_affinities_flag)
  {
    aff_file = fopen(f->affinities_filename, "w");
    if (aff_file == NULL)
      AFFY_HANDLE_ERROR_GOTO("affinities file could not be written",
                             AFFY_ERROR_IO,
                             err,
                             cleanup);
  }
  else if (f->use_saved_affinities)
  {
    aff_file = fopen(f->affinities_filename, "rb");
    if (aff_file == NULL)
      AFFY_HANDLE_ERROR_GOTO("affinities file could not be read",
                             AFFY_ERROR_NOTFOUND,
                             err,
                             cleanup);
  }

  /* ----------- XXX */

  /* Median polish each probeset. The results are stored in probeset */
  for (ps = 0; ps < cdf->numprobesets; ps++)
  {
    AFFY_PROBESET  *p;
    int             numprobes;
    double          value;

    pb_tick(&pbs, 1, "Calculating signal for probe %d",ps+1);
    p = &(cdf->probeset[ps]);
    assert(p != NULL);

    numprobes = p->numprobes;

    z = NULL; /* for safety in case code gets reordered in here */

    /* Allocate storage for z matrix (for polishing) */
    z = create_matrix(numprobes, numchips);
    if (z == NULL)
      AFFY_HANDLE_ERROR_GOTO("create_matrix() out of memory",
                             AFFY_ERROR_OUTOFMEM,
                             err,
                             cleanup);

    /* Load logged data for this probe set */
    for (i = 0; i < numchips; i++)
    {
      /* Load probes from pm array or actual cel data */
      for (j = 0; j < numprobes; j++)
      {
        /* look out, CEL file can contain intensities <= 0 now !!
         *  HACK -- use the delta from the mas5 settings
         */
        value = c->chip[i]->pm[p->probe[j].index];
        if (value < f->delta)
        {
          value = f->delta;
        }
        
        z[j][i] = log(value) / LOG2;
      }
    }

    /* Even if affinities aren't stored, the median polish will record them */
    if (affinity_size < numprobes)
    {
      affinities = halloc(affinities, numprobes * sizeof(double));
      if (affinities == NULL)
        AFFY_HANDLE_ERROR_GOTO("realloc failed",
                               AFFY_ERROR_OUTOFMEM,
                               err,
                               cleanup);
      else
        /* XXX double check this is ok */
        hattach(affinities, mempool);

      affinity_size = numprobes;
    }

    if (f->use_saved_affinities ||
        (f->reuse_affinities && c->mp_populated_flag))
    {
      double t_p, t_g;
      double *affinities_ptr;

      if (f->use_saved_affinities)
      {
        read_affinities(aff_file, affinities, &t_g, numprobes, err);
        AFFY_CHECK_ERROR_GOTO(err, cleanup);
        
        affinities_ptr = affinities;
      }
      else
      {
        t_g = c->t_values[ps];
        affinities_ptr = c->affinities[ps];
      }

      /* At this point the affinity values have been successfully
	 loaded from disk, now we simply calculate the expression
	 value for this probeset on each chip. */
      for (i = 0; i < numchips; i++)
      {
	int probe_idx;

	/* First adjust the bg-corrected probe values for the affinity. */
	for (probe_idx = 0; probe_idx < numprobes; probe_idx++)
	  z[probe_idx][i] -= affinities_ptr[probe_idx];
	
	/* Next, median polish the corrected probe values. */
	affy_rma_median_polish(z, 0, i, numprobes, 1, NULL, NULL, &t_p, f,
	                       err);
        AFFY_CHECK_ERROR_GOTO(err, cleanup);

	/* Calculate and store probeset expression. */
	results[i] = t_p + t_g;
      }
    }
    else
    {
      affy_rma_median_polish(z, 0, 0, numprobes, numchips, results, 
			     affinities, &t, f, err);
      AFFY_CHECK_ERROR_GOTO(err, cleanup);
    }

    /* store affinities for median polish reuse */
    if (f->reuse_affinities &&
        !f->use_saved_affinities &&
        !c->mp_populated_flag)
    {
      double *dptr;
    
      dptr = c->affinities[ps] = h_subcalloc(c->affinities,
                                             numprobes,
                                             sizeof(double *));

      c->t_values[ps] = t;
      for (i = 0; i < numprobes; i++)
        dptr[i] = affinities[i];
    }

    /* Store results */
    for (i = 0; i < numchips; i++)
    {
      c->chip[i]->probe_set[ps] = results[i];
    }

    /* If desired, dump the probe affinities */
    if (f->dump_probe_affinities && safe_to_write_affinities_flag)
    {
      /* T-value and probeset name come first, on their own line. */
      fprintf(aff_file, "%s %.15e\n", p->name, t);

      /* Then each probe and its location, one per line. */
      for (i = 0; i < numprobes; i++)
      {
	AFFY_POINT pt = p->probe[i].pm;
	
	fprintf(aff_file, "%s %" AFFY_PRNd32 " %" AFFY_PRNd32 " %.15e\n",
		p->name, pt.x, pt.y, affinities[i]);
      }
    }

    free_matrix(z);
    z = NULL;
  }

  pb_finish(&pbs, "Finished median polish probeset summarization");

cleanup:
  /* mark reuse affinity arrays as allocated and populated */
  c->mp_allocated_flag = 1;
  c->mp_populated_flag = 1;

  /* Free up results storage */
  h_free(mempool);
  free_matrix(z);

  /* Close affinity file if it was opened. */
  if (aff_file)
    fclose(aff_file);
}
