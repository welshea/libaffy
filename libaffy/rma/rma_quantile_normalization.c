
/**********************************************************
 **
 ** file: rma_quantile_normalization.c
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
 ** Feb 2,  2002 - Intial c code version from original R code
 ** Apr 19, 2002 - Update to deal more correctly with ties (equal rank)
 ** Jan 2,  2003 - Documentation/Commenting updates reformating
 ** Feb 17, 2003 - add in a free(datvec) to qnorm(). clean up freeing of dimat
 ** Feb 25, 2003 - try to reduce or eliminate compiler warnings (with gcc -Wall)
 ** Feb 28, 2003 - update reference to normalization paper in comments
 ** Mar 25, 2003 - ability to use median, rather than mean in so 
 **                called "robust" method
 **
 ** 10/22/10: Use new AFFY_COMBINED_FLAGS instead of AFFY_RMA_FLAGS
 **
 ***********************************************************/

#include <affy_rma.h>

typedef struct
{
  double data;
  int index;
} dataitem;

int qnorm_compare(const void *p1, const void *p2)
{
  double x;
  double y;

  assert(p1 != NULL);
  assert(p2 != NULL);

  x = ((dataitem *) p1)->data;
  y = ((dataitem *) p2)->data;

  if (x < y)
    return (-1);
  if (x > y)
    return (1);

  return (0);
}

/************************************************************
 **
 ** double *rank_order()
 **
 ** get ranks in the same manner as R does. Assume that *x is
 ** already sorted
 ** NB: Can this be in-place now?
 *************************************************************/
void rank_order(double *rank, dataitem *x, unsigned int n)
{
  int i, j, k;

  assert(rank != NULL);
  assert(x    != NULL);

  i = 0;

  while (i < n)
  {
    j = i;

    while ((j < n - 1) && (x[j].data == x[j + 1].data))
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
  /*return rank; */
}

void affy_rma_quantile_normalization_chip(AFFY_CHIPSET *c, 
					  int chipnum,
					  double *mean, 
					  AFFY_COMBINED_FLAGS *f,
                                          AFFY_ERROR *err)
{
  int               i, np, *mempool;
  double           *rank = NULL;
  dataitem         *vals = NULL;
  int               num_probes;
  AFFY_PROBE      **probe_arr;
  LIBUTILS_PB_STATE pbs;

  assert(c             != NULL);
  assert(f             != NULL);
  assert(mean          != NULL);
  assert(c->cdf        != NULL);
  assert(c->cdf->probe != NULL);

  num_probes = c->cdf->numprobes;
  probe_arr = c->cdf->probe;

  pb_init(&pbs);
  pb_begin(&pbs, 2, "Quantile Normalization");

  /* Allocate storage */
  mempool = h_malloc(sizeof(int));
  if (mempool == NULL)
    AFFY_HANDLE_ERROR_VOID("malloc failed", AFFY_ERROR_OUTOFMEM, err);

  vals = h_subcalloc(mempool, num_probes, sizeof(dataitem));
  if (vals == NULL)
    AFFY_HANDLE_ERROR_GOTO("calloc failed", AFFY_ERROR_OUTOFMEM, err, cleanup);

  /* Load in pm for each probe */
  np = 0;

  for (i = 0; i < num_probes; i++)
  {
    if ((f->normalize_affx_probes) 
	|| !(affy_is_control_probe(probe_arr[i])))
    {
      assert(c->chip[chipnum]->pm != NULL);
      vals[np].data  = c->chip[chipnum]->pm[i];
      vals[np].index = i;
      np++;
    }
  }

  pb_tick(&pbs,1,"Accumulating means");
  qsort(vals, np, sizeof(dataitem), qnorm_compare);

  /* Step two: accumulate mean value at a given rank */
  if (!(f->use_saved_means))
  {
    for (i = 0; i < np; i++)
      mean[i] += vals[i].data;
  }

  pb_tick(&pbs,1,"Rank ordering");
  rank = h_subcalloc(mempool, np, sizeof(double));
  if (rank == NULL)
    AFFY_HANDLE_ERROR_GOTO("calloc failed", AFFY_ERROR_OUTOFMEM, err, cleanup);

  /* Rank order the intensities on this chip */
  rank_order(rank, vals, np);

  for (i = 0; i < num_probes; i++)
  {
    int r = vals[i].index;

   if ((f->normalize_affx_probes) 
       || !(affy_is_control_probe(probe_arr[i])))
     c->chip[chipnum]->pm[r] = floor(rank[i]) - 1;
  }

  pb_finish(&pbs,"Finished quantile normalization");

cleanup:
  h_free(mempool);
}

/**
 * Given rankings already computed per chip and an overall mean
 * profile, we can normalize each chip.
 */
void affy_rma_quantile_normalization_chipset(AFFY_CHIPSET *c, 
					     double *mean,
					     AFFY_COMBINED_FLAGS *f)
{
  int          i, j;
  int          numprobes;
  AFFY_PROBE **probe_arr;

  assert(c             != NULL);
  assert(mean          != NULL);
  assert(c->chip       != NULL);
  assert(f             != NULL);
  assert(c->cdf        != NULL);
  assert(c->cdf->probe != NULL);

  numprobes = c->cdf->numprobes;
  probe_arr = c->cdf->probe;

  for (i = 0; i < c->num_chips; i++)
  {
    /* The rank is stored in the pm array */
    for (j = 0; j < numprobes; j++)
    {
      if ((f->normalize_affx_probes) 
	  || !(affy_is_control_probe(probe_arr[j])))
      {
        assert(c->chip[i]->pm != NULL);

	c->chip[i]->pm[j] = mean[(int)c->chip[i]->pm[j]];
      }
    }
  }
}

