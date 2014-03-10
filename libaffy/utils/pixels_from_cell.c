
/**************************************************************************
 *
 * Filename:  pixels_from_cell.c
 *
 * Purpose:   Obtain an AFFY_PIXREGION for a given cell.
 *
 * Creation:  21 Apr 2005 
 *
 * Author:    Andrew M. Hoerter
 *
 * Copyright: Copyright (C) 2007, Moffitt Cancer Center.
 *            All rights reserved.
 *
 * Update History
 * --------------
 * 04/21/05: Creation (AMH)
 * 03/07/08: New error handling scheme (AMH)
 * 09/20/10: Misc updates (AMH)
 * 03/10/14: #ifdef out CEL qc fields to save memory (EAW)
 *
 **************************************************************************/

#include <affy.h>

#ifdef STORE_CEL_QC
/* 
 * affy_pixels_from_cell(): Return the AFFY_PIXREGION corresponding to a 
 *                          cell.
 *
 * Inputs: cp points to an initialized/loaded AFFY_CHIP, x/y specify the
 *         cell coordinates. 
 *        
 *
 * Outputs: Pointer to an AFFY_PIXREGION or NULL on error.
 *
 * Side effects: The AFFY_CELL structure is updated to cache the result.
 *
 * TODO: Boundary conditions.
 *
 *
 */
AFFY_PIXREGION *affy_pixels_from_cell(AFFY_CHIP *cp, 
                                      int x, 
                                      int y, 
                                      AFFY_ERROR *err)
{
  int             i;
  AFFY_POINT      px, px_r,px_l;
  AFFY_PIXREGION *px_result;

  assert(cp      != NULL);
  assert(cp->dat != NULL);

  /* Check if already cached */
  if (   (cp->cel != NULL) 
      && (cp->cel->data != NULL) 
      && (cp->cel->data[x][y].pixels != NULL))
    return (cp->cel->data[x][y].pixels);

  /* 
   * Otherwise, need to make a new one.  Note that freeing a DAT file 
   * structure will free all associated PIXREGION's.
   */
  px_result = h_subcalloc(cp->dat, 1, sizeof(AFFY_PIXREGION));
  if (px_result == NULL)
    AFFY_HANDLE_ERROR("calloc failed", AFFY_ERROR_OUTOFMEM, err, NULL);

  /* 
   * First get the starting position in image, and the right/lower
   * boundaries
   */
  px   = affy_cell_to_pixel(cp, x, y);
  px_r = affy_cell_to_pixel(cp, x+1, y);
  px_l = affy_cell_to_pixel(cp, x, y+1);

  px_result->numcols = px_r.x-px.x;
  px_result->numrows = px_l.y-px.y;
  px_result->data    = h_subcalloc(px_result, 
                                   px_result->numrows, 
                                   sizeof(unsigned int *));
  if (px_result->data == NULL)
  {
    h_free(px_result);
    AFFY_HANDLE_ERROR("calloc failed", AFFY_ERROR_OUTOFMEM, err, NULL);
  }

  for (i = 0; i < px_result->numrows; i++) 
    px_result->data[i] = &(cp->dat->pixels.data[px.y+i][px.x]);
  
  if (cp->cel != NULL)
    cp->cel->data[x][y].pixels = px_result;

  return (px_result);
}
#endif    /* STORE_CEL_QC */

/*
 * Convert a cell location (a cell is a collection of pixels) to image
 * coordinates. This requires rotation of the image which is done in
 * real time via a bilinear algorithm
 * 
 * Detecting and Correcting Misalignment in Affymetrix Data
 * Keith A. Baggerly
 * Technical Report
 * MD Anderson Cancer Center
 */
AFFY_POINT affy_cell_to_pixel(AFFY_CHIP *cp, int x, int y)
{
  int        ax, bx, cx, dx, ay, by, cy, dy;
  double     rows, cols;
  double     xn, yn, xp, yp, newx, newy;
  AFFY_POINT newp;

  assert(cp != NULL);
	
  /* Account for corners being off by one */
  ax = cp->dat->grid_ul.x;
  ay = cp->dat->grid_ul.y;
  bx = cp->dat->grid_ur.x + 1;
  by = cp->dat->grid_ur.y;
  cx = cp->dat->grid_ll.x;
  cy = cp->dat->grid_ll.y + 1;
  dx = cp->dat->grid_lr.x + 1;
  dy = cp->dat->grid_lr.y + 1;

  /* Get chip-specific rows/cols */
  rows = cp->cdf->numrows;
  cols = cp->cdf->numcols;
	
  /* 
   * Calculate x/y coordinates for current position, and one below and
   * to right 
   */
  xn = ax*((cols-x)/cols)*((rows-y)/rows) 
    + bx*(x/cols)*((rows-y)/rows) 
    + cx*((cols-x)/cols)*(y/rows) 
    + dx*(x/cols)*(y/rows);
			
  yn = ay*((cols-x)/cols)*((rows-y)/rows) 
    + by*(x/cols)*((rows-y)/rows) 
    + cy*((cols-x)/cols)*(y/rows) 
    + dy*(x/cols)*(y/rows);
  y++;
  xp = ax*((cols-x)/cols)*((rows-y)/rows) 
    + bx*(x/cols)*((rows-y)/rows) 
    + cx*((cols-x)/cols)*(y/rows) 
    + dx*(x/cols)*(y/rows);
  y--;
  x++;
  yp = ay*((cols-x)/cols)*((rows-y)/rows) 
    + by*(x/cols)*((rows-y)/rows) 
    + cy*((cols-x)/cols)*(y/rows) 
    + dy*(x/cols)*(y/rows);
  x--;

  /* Take the average of the current position and next one */
  newx = (xn+xp)/2.0;
  newy = (yn+yp)/2.0;

  /* Assign new values to a point, to return */
  newp.x = (int)(newx+0.5);
  newp.y = (int)(newy+0.5);

  return (newp);
}
