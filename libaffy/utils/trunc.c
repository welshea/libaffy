
/**************************************************************************
 *
 * Filename:  trunc.c
 *
 * Purpose:   Portable implementation of trunc()
 *
 * Creation: 
 *
 * Author:    Andrew M. Hoerter
 *
 * Copyright: Copyright (C) 2008, Moffitt Cancer Center.
 *            All rights reserved.
 *
 * Update History
 * --------------
 * 10/28/08: Factored out from pnorm.c (AMH)
 *
 **************************************************************************/

#include <math.h>

#include <utils.h>

/*
 * trunc() is not available in C89.  Provide an implementation of it
 * here.
 */
double affy_trunc(double x)
{
  double result = 0;

  if (x < 0) 
    result = ceil(x);
  else         
    result = floor(x);
  
  return (result);
}
