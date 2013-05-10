
/**************************************************************************
 *
 * Filename:  cel_file_util.c
 *
 * Purpose:   Utility functions for CEL files.
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
 * 04/08/05: Imported/fixed up from old libaffy (AMH)
 * 04/19/05: Now we use AFFY_CELL's, plus a utility function (AMH)
 * 03/11/08: New error handling scheme (AMH)
 * 09/20/10: Pooled memory allocator (AMH)
 * 05/10/13: Added affy_mostly_free_cel_file() (EAW)
 *
 **************************************************************************/

#include <affy.h>
#include <utils.h>

void affy_free_cel_file(AFFY_CELFILE *cf)
{
  h_free(cf);
}

/* free the matricies, keeping some useful information about the CEL file */
void affy_mostly_free_cel_file(AFFY_CELFILE *cf)
{
  if (cf)
  {
    if (cf->data)    h_free(cf->data);
    if (cf->mask)    h_free(cf->mask);
    if (cf->outlier) h_free(cf->outlier);
    
    cf->data    = NULL;
    cf->mask    = NULL;
    cf->outlier = NULL;
  }
}

/* 
 * affy_matrix_from_cel(): Extract cell value matrix from an AFFY_CELFILE
 *
 * Inputs: A pointer to a valid AFFY_CELFILE structure.
 *
 * Outputs: A pointer to the 2D array of values on success, NULL on failure
 *
 * Side effects: Memory is allocated which must be freed by the caller.
 *
 * Comments: For various processing it is desireable to have quick,
 *           direct access to cell values which are ordinarily hidden
 *           within AFFY_CELL's.  Thus this handy utility function to
 *           pull out those values and place them in a block of doubles.
 */

double **affy_matrix_from_cel(AFFY_CELFILE *cf, AFFY_ERROR *err)
{
  int      x, y;
  double **matrix;

  assert(cf != NULL);

  /* Allocate the matrix. */
  if ((matrix = create_matrix(cf->numrows, cf->numcols)) == NULL)
    AFFY_HANDLE_ERROR("out of memory creating matrix", 
		      AFFY_ERROR_OUTOFMEM, 
		      err, 
		      NULL);

  /* 
   * Copy the values over. XXX -- use pointers to data within AFFY_CELL
   * instead?
   */
  for (y = 0; y < cf->numrows; y++)
    for (x = 0; x < cf->numcols; x++)
      matrix[y][x] = cf->data[y][x].value;

  return (matrix);
}
