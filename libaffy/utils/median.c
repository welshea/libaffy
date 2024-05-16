
/**************************************************************************
 *
 * Filename:  median.c
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
 * 04/14/05: Imported/repaired from old libaffy (AMH)
 * 04/19/07: Added support for incremental RMA (AMH)
 * 03/07/08: New error handling scheme (AMH)
 * 09/20/10: Pooled memory allocator (AMH)
 * 11/19/10: Fixed median calculations, added bioconductor compatability (EAW)
 * 08/12/20: removed bioconductor compatiblity, both gave same results (EAW)
 * 05/16/24: optimized median math (EAW)
 *
 **************************************************************************/

#include "utils.h"
#include "affy.h"

/**************************************************************************
 **
 ** double affy_median(double *x, int length, AFFY_COMBINED_FLAGS *f)
 **
 ** double *x - vector
 ** int length - length of *x
 **
 ** returns the median of *x
 **
 *************************************************************************/

/*
 * Copy the results over to a temporary array, then get the median.
 * This preserves the initial ordering.
 */
double affy_median_save(double *x, int length, AFFY_COMBINED_FLAGS *f,
                        AFFY_ERROR *err)
{
  int     i;
  double *buffer;
  double  result;

  buffer = h_malloc(length * sizeof(double));
  if (buffer == NULL)
    AFFY_HANDLE_ERROR("malloc failed", AFFY_ERROR_OUTOFMEM, err, 0.0);

  for (i = 0; i < length; i++)
    buffer[i] = x[i];

  result = affy_median(buffer, length, f);

  h_free(buffer);

  return (result);
}

/* 
 *  This function will sort x and find the median. Note this is 
 *  destructive, therefore you should call affy_median_save() if you 
 *  care about this.
 */
double affy_median(double *x, int length, AFFY_COMBINED_FLAGS *f)
{
  int    half;
  double med;

  assert(x != NULL);

  qsort(x, length, sizeof(double), affy_median_sort);


  half = length >> 1;

  if (length % 2)
    med = x[half];
  else
    med = 0.5 * (x[half - 1] + x[half]);

  return (med);
}

/*****************************************************************************
 **
 ** void affy_get_row_median(double *z, double *rdelta, int startrow,
 **                          int startcol, int rows, int cols,
 **                          AFFY_COMBINED_FLAGS *f)
 **
 ** double *z - matrix of dimension  rows*cols
 ** double *rdelta - on output will contain row medians (vector of length rows)
 ** int rows, cols - dimension of matrix
 ** int startrow, startcol - beginning row/col for working with submatrices
 **
 ** get the row medians of a matrix 
 **
 *****************************************************************************/
void affy_get_row_median(double **z, 
                         double *rdelta, 
                         int startrow, 
			 int startcol, 
                         int numrows, 
                         int numcolumns,
                         AFFY_COMBINED_FLAGS *f,
                         AFFY_ERROR *err)
{
  int     i, j;
  int     rowsleft, colsleft;
  double *buffer;

  assert(z      != NULL);
  assert(rdelta != NULL);

  buffer = h_malloc(numcolumns * sizeof(double));
  if (buffer == NULL)
    AFFY_HANDLE_ERROR_VOID("malloc failed", AFFY_ERROR_OUTOFMEM, err);

  for (i = startrow, rowsleft = numrows; rowsleft > 0; i++, rowsleft--)
  {
    for (j = startcol, colsleft = numcolumns; colsleft > 0; j++, colsleft--)
      buffer[numcolumns-colsleft] = z[i][j];

    rdelta[numrows-rowsleft] = affy_median(buffer, numcolumns, f);
  }

  h_free(buffer);
}

/*****************************************************************************
 **
 ** void affy_get_col_median(double *z, double *cdelta, int startrow, 
 **                          int startcol, int rows, int cols,
 **                          AFFY_COMBINED_FLAGS *f)
 **
 ** double *z - matrix of dimension  rows*cols
 ** double *cdelta - on output will contain col medians (vector of length cols)
 ** int rows, cols - dimension of matrix
 ** int startrow, startcol - beginning row/col for working with submatrices
 **
 ** get the col medians of a matrix 
 **
 *****************************************************************************/
void affy_get_column_median(double **z, 
                            double *cdelta, 
                            int startrow, 
			    int startcol, 
                            int numrows, 
                            int numcolumns,
                            AFFY_COMBINED_FLAGS *f,
                            AFFY_ERROR *err)
{
  int     i, j;
  int     rowsleft, colsleft;
  double *buffer;

  assert(z      != NULL);
  assert(cdelta != NULL);

  buffer = h_malloc(numrows * sizeof(double));
  if (buffer == NULL)
    AFFY_HANDLE_ERROR_VOID("malloc failed", AFFY_ERROR_OUTOFMEM, err);

  for (j = startcol, colsleft = numcolumns; colsleft > 0; j++, colsleft--)
  {
    for (i = startrow, rowsleft = numrows; rowsleft > 0; i++, rowsleft--)
      buffer[numrows-rowsleft] = z[i][j];

    cdelta[numcolumns-colsleft] = affy_median(buffer, numrows, f);
  }

  h_free(buffer);
}

/* 
 * Sorting function for qsort(). Given two double pointers, return the
 * ordering.
 */
int affy_median_sort(const void *p1, const void *p2)
{
  double n1 = *((double *)p1);
  double n2 = *((double *)p2);

  if (n1 < n2)
    return (-1);
  if (n1 > n2)
    return (1);

  return (0);
}
