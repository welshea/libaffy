
/**************************************************************************
 *
 * Filename:  isundefined.c
 *
 * Purpose:   Cell-type testing predicate.
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
 * 4/14/05: Imported/repaired from old libaffy (AMH)
 *
 **************************************************************************/

#include <utils.h>
#include <affy.h>

bool affy_isundefined(AFFY_CHIP *chip, int x, int y)
{
  assert(chip      != NULL);
  assert(chip->cdf != NULL);
  assert(x < chip->cdf->numcols);
  assert(y < chip->cdf->numrows);
  assert(x >= 0);
  assert(y >= 0);

  if (chip->cdf->cell_type[x][y] == AFFY_UNDEFINED_LOCATION)
    return (true);
  else
    return (false);
}
