/**************************************************************************
 *
 * Filename: mean.c
 *
 * Purpose:  Statistical routines.
 *
 * Creation: 18 November, 2010
 *
 * Author:    Eric A. Welsh
 *
 * Copyright: Copyright (C) 2010, Moffitt Cancer Center.
 *            All rights reserved.
 *
 * Update History
 * --------------
 * 11/18/10: File creation (EAW)
 * 03/11/11: Added geometric mean function (EAW)
 *
 **************************************************************************/

#include <utils.h>
#include <affy.h>

/**************************************************************************
 **
 ** double affy_mean(double *x, int length)
 **
 ** double *x - vector
 ** int length - length of *x
 **
 ** returns the mean of *x
 **
 *************************************************************************/

/* 
 *  This function will find the mean.
 */
double affy_mean(double *x, int length)
{
  double mean = 0.0;
  int i;

  assert(x != NULL);
  
  for (i = 0; i < length; i++)
    mean += x[i];
  if (length)
    mean /= length;

  return (mean);
}


/**************************************************************************
 **
 ** double affy_mean_geometric_floor_1(double *x, int length)
 **
 ** double *x - vector
 ** int length - length of *x
 **
 ** returns the geometric mean of *x
 **
 *************************************************************************/

/* 
 *  This function will find the geometric mean, flooring values to 1
 */
double affy_mean_geometric_floor_1(double *x, int length)
{
  double mean = 0.0;
  double value;
  int i;

  assert(x != NULL);
  
  for (i = 0; i < length; i++)
  {
    value = x[i];
    if (value < 1.0)
      value = 1.0;

    mean += log(value);
  }
  if (length)
    mean /= length;

  return (exp(mean));
}
