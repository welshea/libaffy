
/**************************************************************************
 *
 * Filename:  quantile_normalizations.c
 *
 * Purpose:   Statistical routines.
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
 * 04/08/05: Imported/repaired from old libaffy (AMH)
 * 04/19/05: AFFY_CELFILE now uses AFFY_CELL's (AMH)
 * 03/07/08: New error handling scheme (AMH)
 * 09/20/10: Pooled memory allocator (AMH)
 * 03/12/13: added affy_quantile_normalize_probeset()
 * 03/14/14: fixed to work with exons arrays (EAW)
 *
 **************************************************************************/

#include <affy.h>
#include <utils.h>

/**********************************************************
 **
 ** file: qnorm.c
 **
 ** aim: A c implementation of the quantile normalization method 
 **
 ** Copyright (C) 2002-2003    Ben Bolstad
 **
 ** written by: B. M. Bolstad  <bolstad@stat.berkeley.edu>
 **
 ** written: Feb 2, 2002
 ** last modified: Apr 19, 2002
 ** 
 ** This c code implements the quantile normalization method
 ** for normalizing high density oligonucleotide data as discussed
 ** in
 **
 ** Bolstad, B. M., Irizarry R. A., Astrand, M, and Speed, T. P. (2003)(2003) 
 ** A Comparison of Normalization Methods for High 
 ** Density Oligonucleotide Array Data Based on Bias and Variance.
 ** Bioinformatics 19,2,pp 185-193
 **
 ** History
 ** Feb 2, 2002 - Intial c code version from original R code
 ** Apr 19, 2002 - Update to deal more correctly with ties (equal rank)
 ** Jan 2, 2003 - Documentation/Commenting updates reformating
 ** Feb 17, 2003 - add in a free(datvec) to qnorm(). clean up freeing of dimat
 ** Feb 25, 2003 - try to reduce or eliminate compiler warnings (with gcc -Wall)
 ** Feb 28, 2003 - update reference to normalization paper in comments
 ** Mar 25, 2003 - ability to use median, rather than mean in so called "robust" method
 ** sometime in 2012-2013? - added MAS5-style functions
 ** 2014-03-13   - modified to work with exon arrays
 **
 ***********************************************************/


/*
 * This normalizes both PM and MM values, which is different from the RMA
 * version, which only normalizes the PM values.
 */
void affy_quantile_normalization(AFFY_CHIPSET *d, affy_uint8 pm_only,
                                 AFFY_ERROR *err)
{
  int           p, i, j, x, y, number_of_probes;
  int           num_seen_probes = 0;
  double        sum;
  double       *rank, *mean;
  AFFY_CELFILE *cf;
  double     ***all_vals;
  int          *qnorm_pool;

  if (pm_only)
    info("Quantile normalization (PM-only)...");
  else
    info("Quantile normalization (PM, MM)...");

  number_of_probes = d->cdf->numprobes;

  /* Allocate a pool pointer to hang all our storage on */
  qnorm_pool = h_malloc(sizeof(int));
  if (qnorm_pool == NULL)
    AFFY_HANDLE_ERROR_VOID("malloc failed", AFFY_ERROR_OUTOFMEM, err);

  mean = (double *)h_subcalloc(qnorm_pool, number_of_probes, sizeof(double));
  if (mean == NULL)
  {
    h_free(qnorm_pool);
    AFFY_HANDLE_ERROR_VOID("calloc failed", AFFY_ERROR_OUTOFMEM, err);
  }

  all_vals = (double ***)h_subcalloc(qnorm_pool, 
                                     d->num_chips, 
                                     sizeof(double **));
  if (all_vals == NULL)
  {
    h_free(qnorm_pool);
    AFFY_HANDLE_ERROR_VOID("calloc failed", AFFY_ERROR_OUTOFMEM, err);
  }

  /* We start by creating a matrix of double pointers to data */
  for (i = 0; i < d->num_chips; i++)
  {
    cf = d->chip[i]->cel;
    j = 0;

    /* Allocate storage for this chip */
    all_vals[i] = (double **)h_subcalloc(qnorm_pool, 
                                         number_of_probes, 
                                         sizeof(double *));
    if (all_vals[i] == NULL)
    {
      h_free(qnorm_pool);
      AFFY_HANDLE_ERROR_VOID("calloc failed", AFFY_ERROR_OUTOFMEM, err);
    }

    /* Load in the pm for each probe */
    memset(d->cdf->seen_xy[0], 0,
           d->cdf->numrows * d->cdf->numcols * sizeof(affy_uint8));
    for (p = 0, j = num_seen_probes = 0; p < number_of_probes; p++)
    {
      x = d->cdf->probe[p]->pm.x;
      y = d->cdf->probe[p]->pm.y;
      
      if (d->cdf->seen_xy[x][y] == 0)
        all_vals[i][j++] = &(cf->data[x][y].value);
      d->cdf->seen_xy[x][y] = 1;
      
      if (pm_only)
        continue;

      /* hack for missing MM probes, where MM coords == PM coords
       */
      if (d->cdf->probe[p]->pm.x == d->cdf->probe[p]->mm.x &&
          d->cdf->probe[p]->pm.y == d->cdf->probe[p]->mm.y)
      {
        continue;
      }
      else
      {
        x = d->cdf->probe[p]->mm.x;
        y = d->cdf->probe[p]->mm.y;
      
        if (d->cdf->seen_xy[x][y] == 0)
          all_vals[i][j++] = &(cf->data[x][y].value);
        d->cdf->seen_xy[x][y] = 1;
      }
    }
    num_seen_probes = j;

    qsort(all_vals[i], num_seen_probes, sizeof(double *), affy_qnorm_compare);
  }

  /* Step two: calculate mean value at a given rank */
  for (j = 0; j < num_seen_probes; j++)
  {
    sum = 0;
    /* Calculate mean */
    for (i = 0; i < d->num_chips; i++)
      sum += *(all_vals[i][j]) / (double)d->num_chips;

    mean[j] = sum;
  }

  rank = (double *)h_subcalloc(qnorm_pool, num_seen_probes, sizeof(double));
  if (rank == NULL)
  {
    h_free(qnorm_pool);
    AFFY_HANDLE_ERROR_VOID("calloc failed", AFFY_ERROR_OUTOFMEM, err);
  }

  /* Step 3: Redistribute mean value to all chips */
  for (i = 0; i < d->num_chips; i++)
  {
    /* Rank order the intensities on this chip */
    affy_rank_order(rank, all_vals[i], num_seen_probes);
    for (j = 0; j < num_seen_probes; j++)
    {
      *(all_vals[i][j]) = mean[(int)floor(rank[j]) - 1];
    }
  }

  h_free(qnorm_pool);

  info("done.\n");
}


/*
 * This normalizes probesets
 */
void affy_quantile_normalization_probeset(AFFY_CHIPSET *d, AFFY_ERROR *err)
{
  int           i, j, x, y, number_of_probesets;
  double        sum;
  double       *rank, *mean;
  AFFY_CDFFILE *cdf;
  double     ***all_vals;
  int          *qnorm_pool;
  
  assert(d->cdf != NULL);
  cdf = d->cdf;

  info("Quantile normalization (probesets)...");

  number_of_probesets = d->cdf->numprobesets;

  /* Allocate a pool pointer to hang all our storage on */
  qnorm_pool = h_malloc(sizeof(int));
  if (qnorm_pool == NULL)
    AFFY_HANDLE_ERROR_VOID("malloc failed", AFFY_ERROR_OUTOFMEM, err);

  mean = (double *)h_subcalloc(qnorm_pool, number_of_probesets, sizeof(double));
  if (mean == NULL)
  {
    h_free(qnorm_pool);
    AFFY_HANDLE_ERROR_VOID("calloc failed", AFFY_ERROR_OUTOFMEM, err);
  }

  all_vals = (double ***)h_subcalloc(qnorm_pool, 
                                     d->num_chips, 
                                     sizeof(double **));
  if (all_vals == NULL)
  {
    h_free(qnorm_pool);
    AFFY_HANDLE_ERROR_VOID("calloc failed", AFFY_ERROR_OUTOFMEM, err);
  }

  /* We start by creating a matrix of double pointers to data */
  for (i = 0; i < d->num_chips; i++)
  {
    /* Allocate storage for this chip */
    all_vals[i] = (double **)h_subcalloc(qnorm_pool, 
                                         number_of_probesets,
                                         sizeof(double *));
    if (all_vals[i] == NULL)
    {
      h_free(qnorm_pool);
      AFFY_HANDLE_ERROR_VOID("calloc failed", AFFY_ERROR_OUTOFMEM, err);
    }

    /* store pointers to each probeset value */
    for (j = 0; j < number_of_probesets; j++)
    {
      all_vals[i][j] = &(d->chip[i]->probe_set[j]);
    }

    qsort(all_vals[i], number_of_probesets, sizeof(double *), affy_qnorm_compare);
  }

  /* Step two: calculate mean value at a given rank */
  for (j = 0; j < number_of_probesets; j++)
  {
    sum = 0;
    /* Calculate mean */
    for (i = 0; i < d->num_chips; i++)
      sum += *(all_vals[i][j]) / (double)d->num_chips;

    mean[j] = sum;
  }

  rank = (double *)h_subcalloc(qnorm_pool, number_of_probesets, sizeof(double));
  if (rank == NULL)
  {
    h_free(qnorm_pool);
    AFFY_HANDLE_ERROR_VOID("calloc failed", AFFY_ERROR_OUTOFMEM, err);
  }

  /* Step 3: Redistribute mean value to all chips */
  for (i = 0; i < d->num_chips; i++)
  {
    /* Rank order the intensities on this chip */
    affy_rank_order(rank, all_vals[i], number_of_probesets);
    for (j = 0; j < number_of_probesets; j++)
    {
      *(all_vals[i][j]) = mean[(int)floor(rank[j]) - 1];
    }
  }

  h_free(qnorm_pool);

  info("done.\n");
}


int affy_qnorm_compare(const void *p1, const void *p2)
{
  double x = **((double **)p1);
  double y = **((double **)p2);

  if (x < y)
    return (-1);
  if (x > y)
    return (1);
  return (0);
}

/************************************************************
 **
 ** double *affy_rank_order()
 **
 ** get ranks in the same manner as R does. Assume that *x is
 ** already sorted
 **
 *************************************************************/
void affy_rank_order(double *rank, double **x, int n)
{
  int i, j, k;

  i = 0;

  while (i < n)
  {
    j = i;

    while ((j < n - 1) && (*(x[j]) == *(x[j + 1])))
      j++;

    if (i != j)
    {
      for (k = i; k <= j; k++)
        rank[k] = (i + j + 2) / 2.0;
    }
    else
    {
      rank[i] = i + 1;
    }

    i = j + 1;
  }
}
