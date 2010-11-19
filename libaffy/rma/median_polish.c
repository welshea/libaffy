
/**************************************************************************
 *
 * Filename: median_polish.c
 *
 * Purpose: 
 *
 * Creation: 
 *
 * Author:   Steven Eschrich
 *
 *
 * Update History
 * --------------
 * 04/14/05: Imported/repaired from old libaffy (AMH)
 * 06/05/08: New error handling scheme (AMH)
 * 09/20/10: Pooled memory allocator (AMH)
 * 11/19/10: Pass flags to affy_median_polish() (EAW)
 *
 **************************************************************************/

/* This function was adapted from the bioconductor package */

/*****************************************************************************
 
  void median_polish(double *data, int rows, int cols, int *cur_rows, 
		     double *results, int nprobes, AFFY_COMBINED_FLAGS *f)

  double *data - a data matrix of dimension rows by cols (the entire PM matrix)
  int rows, cols - rows and columns dimensions of matrix
  int cur_rows - vector of length nprobes containg row indicies of *data 
                 matrix which apply for a particular probeset
  double *results - a vector of length cols already allocated. on output 
		   contains expression values
  int nprobes - number of probes in current probeset.

  a function to do median polish.

 *****************************************************************************/

#include <affy_rma.h>

static const int max_iterations = 10;

/* Add xdelta to x, for n elements */
static INLINE void vector_add(double *x, double *xdelta, unsigned int n)
{
  int i;
  
  assert(x      != NULL);
  assert(xdelta != NULL);

  for (i = 0; i < n; i++)
    x[i] += xdelta[i];
}

static INLINE double sum_abs(double **z, 
                             int startingrow, 
                             int startingcol, 
			     int numrows, 
                             int numcols)
{
  int    i, j;
  int    rowsleft, colsleft;
  double sum = 0.0;

  assert(z != NULL);

  for (i = startingrow, rowsleft = numrows; rowsleft > 0; i++, rowsleft--)
    for (j = startingcol, colsleft = numcols; colsleft > 0; j++, colsleft--)
      sum += fabs(z[i][j]);

  return (sum);

/*   int i, j; */
/*   double sum = 0.0; */

/*   for (i = 0; i < numrows; i++) */
/*     for (j = 0; j < numcols; j++) */
/*       sum += fabs(z[i][j]); */

/*   return sum; */

}

static INLINE void subtract_by_row(double **z, 
                                   double *rdelta, 
				   int startingrow, 
                                   int startingcol, 
				   int numrows, 
                                   int numcols)
{
  int i, j, rowsleft, colsleft;

  assert(z != NULL);

  for (i = startingrow, rowsleft = numrows; rowsleft > 0; i++, rowsleft--)
    for (j = startingcol, colsleft = numcols; colsleft > 0; j++, colsleft--)
      z[i][j] -= rdelta[numrows-rowsleft];

/*   int i, j; */

/*   for (i = 0; i < numrows; i++) */
/*     for (j = 0; j < numcols; j++) */
/*       z[i][j] -= rdelta[i]; */
/*   return; */

}

static INLINE void subtract_by_col(double **z, 
                                   double *cdelta, 
				   int startingrow, 
                                   int startingcol, 
				   int numrows, 
                                   int numcols)
{
  int i, j, rowsleft, colsleft;

  assert(z      != NULL);
  assert(cdelta != NULL);

/*   fprintf(stderr, "startingrow = %d, startingcol = %d, numrows = %d, " */
/* 	  "numcols = %d\n", startingrow, startingcol, numrows, numcols); */
/*   fflush(stderr); */

  for (j = startingcol, colsleft = numcols; colsleft > 0; j++, colsleft--)
  {
/*     fprintf(stderr, "colsleft = %d, j = %d\n", colsleft, j); */
/*     fflush(stderr); */
    for (i = startingrow, rowsleft = numrows; rowsleft > 0; i++, rowsleft--)
    {
/*       fprintf(stderr, "rowsleft = %d, i = %d\n", rowsleft, i); */
/*       fprintf(stderr, "z[%d][%d]  cdelta[%d]\n", i, j, numcols-colsleft); */
/*       fflush(stderr); */
      z[i][j] -= cdelta[numcols-colsleft];
    }
  }

/*   int i, j; */

/*   for (j = 0; j < numcols; j++) */
/*   { */
/*     for (i = 0; i < numrows; i++) */
/*     { */
/*       fprintf(stderr, "z[%d][%d]  cdelta[%d]\n", i, j, j); */
     
/*  z[i][j] -= cdelta[j]; */
/*     } */
/*   } */
}

/*
  Do a median polish on the chipset c, using probe set # ps. The results
  are stored in the corresponding probe_set array for each chip.
*/
void affy_rma_median_polish(double **z, 
                            int startingprobe, 
                            int startingchip, 
			    int numprobes, 
                            int numchips, 
                            double *results, 
			    double *affinities,
                            double *t_val,
                            AFFY_COMBINED_FLAGS *f,
                            AFFY_ERROR *err)
{
  double *rdelta = NULL, *r = NULL, *cdelta = NULL, *c = NULL;
  double eps = 0.01;
  double oldsum = 0.0, newsum = 0.0;
  double t = 0.0;
  double delta;
  int    i, j, it, *mempool;

  assert(z    != NULL);
  assert(z[0] != NULL);

  mempool = h_malloc(sizeof(int));
  if (mempool == NULL)
    AFFY_HANDLE_ERROR_VOID("malloc failed", AFFY_ERROR_OUTOFMEM, err);

  rdelta = h_subcalloc(mempool, numprobes, sizeof(double));
  if (rdelta == NULL)
    AFFY_HANDLE_ERROR_GOTO("calloc failed", AFFY_ERROR_OUTOFMEM, err, cleanup);

  cdelta = h_subcalloc(mempool, numchips, sizeof(double));
  if (rdelta == NULL)
    AFFY_HANDLE_ERROR_GOTO("calloc failed", AFFY_ERROR_OUTOFMEM, err, cleanup);

  r = h_subcalloc(mempool, numprobes, sizeof(double));
  if (rdelta == NULL)
    AFFY_HANDLE_ERROR_GOTO("calloc failed", AFFY_ERROR_OUTOFMEM, err, cleanup);

  c = h_subcalloc(mempool, numchips, sizeof(double));
  if (rdelta == NULL)
    AFFY_HANDLE_ERROR_GOTO("calloc failed", AFFY_ERROR_OUTOFMEM, err, cleanup);

  /* Do the actual median polish */
  for (it = 0; it < max_iterations; it++)
  {
    /* Rows first */
    affy_get_row_median(z, 
                        rdelta, 
                        startingprobe, 
                        startingchip, 
                        numprobes, 
			numchips,
			f,
                        err);
    AFFY_CHECK_ERROR_GOTO(err, cleanup);

    subtract_by_row(z, 
                    rdelta, 
                    startingprobe, 
                    startingchip, 
                    numprobes, 
		    numchips);
    vector_add(r, rdelta, numprobes);

    delta = affy_median_save(c, numchips, f, err);
    AFFY_CHECK_ERROR_GOTO(err, cleanup);

    for (j = 0; j < numchips; j++)
      c[j] -= delta;
    t += delta;

    /* Then columns */
    affy_get_column_median(z, 
                           cdelta, 
                           startingprobe, 
                           startingchip, 
                           numprobes, 
			   numchips,
			   f,
                           err);
    AFFY_CHECK_ERROR_GOTO(err, cleanup);

    subtract_by_col(z, cdelta, startingprobe, startingchip, numprobes, 
		    numchips);
    vector_add(c, cdelta, numchips);

    delta = affy_median_save(r, numprobes, f, err);
    AFFY_CHECK_ERROR_GOTO(err, cleanup);

    for (i = 0; i < numprobes; i++)
      r[i] -= delta;

    t += delta;

    /* If changes are small, then quit */
    newsum = sum_abs(z, startingprobe, startingchip, numprobes, numchips);

    if ((newsum == 0.0) || (fabs(1.0 - oldsum / newsum) < eps))
      break;

    oldsum = newsum;
  }

  /* The result is the expression value */
  if (results != NULL)
  {
    for (j = 0; j < numchips; j++)
      results[j] = t + c[j];
  }

  if (affinities != NULL)
  {
    for (j = 0; j < numprobes; j++)
      affinities[j] = t + r[j];
  }

  if (t_val != NULL)
    *t_val = t;

cleanup:
  h_free(mempool);
}
