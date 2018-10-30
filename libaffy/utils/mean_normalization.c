
/**************************************************************************
 *
 * Filename:  mean_normalization.c
 *
 * Purpose:   Statistical normalization routines.
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
 * 1/04/07: Add various assertions for debugging (AMH)
 * 4/14/05: Imported/repaired from old libaffy (AMH)
 * 4/19/05: AFFY_CELFILE now uses AFFY_CELL's (AMH)
 * 3/10/11: Works on both PM/MM and PM-only data, rather than segfault (EAW)
 * 3/06/14: Fixed x/y numcols/numrows assertions (they were swapped) (EAW)
 * 3/06/14: Further fixed missing MM functionality, so that it works correctly now (EAW)
 * 3/06/14: affy_mean_normalization only uses points >0 now
 * 3/14/14: fixed to work with exon arrays (EAW)
 * 10/29/18: handle 1:many probe:probeset
 *
 **************************************************************************/

#include <affy.h>
#include <utils.h>

/* 
 *  affy_mean_normalization()
 *
 *  Normalize the chips to the same constant mean intensity
 */
void affy_mean_normalization(AFFY_CHIPSET *d, double target_mean)
{
  int           i, j, x, y, n;
  double        mean, value;
  AFFY_CELFILE *cf;
  AFFY_CDFFILE *cdf = d->cdf;
  int           number_of_probes;

  assert(d   != NULL); 
  assert(cdf != NULL);

  info("Performing mean normalization...");

  number_of_probes = cdf->numprobes;

  for (i = 0; i < d->num_chips; i++)
  {
    assert(d->chip[i]      != NULL);
    assert(d->chip[i]->cel != NULL);

    cf = d->chip[i]->cel;
    mean = 0;

    /* Calculate mean */
    
    /* both PM and MM */
    n = 0;
    if (cf->data)
    {
      memset(cdf->seen_xy[0], 0, cdf->numrows*cdf->numcols*sizeof(affy_uint8));
      for (j = 0; j < number_of_probes; j++)
      {
        assert(cdf->probe[j] != NULL);

        x = cdf->probe[j]->pm.x;
        y = cdf->probe[j]->pm.y;

        assert(x >= 0);
        assert(x < cf->numcols);
        assert(y >= 0);
        assert(y < cf->numrows); 

#if 0
        debug("row ptr = %p\n", cf->data + x); 
        debug("AFFY_CELL %p\n", cf->data[x] + y);
#endif

        if (cdf->seen_xy[x][y] == 0)
        {
          value = cf->data[x][y].value;
          if (value > 0)
          {
            mean += value;
            n++;
          }
        }
        cdf->seen_xy[x][y] = 1;

        /* hack for missing MM probes, where MM coords == PM coords
         */
        if (cdf->probe[j]->pm.x == cdf->probe[j]->mm.x &&
            cdf->probe[j]->pm.y == cdf->probe[j]->mm.y)
        {
          continue;
        }
        else
        {
          x = cdf->probe[j]->mm.x;
          y = cdf->probe[j]->mm.y;

          assert(x >= 0);
          assert(x < cf->numcols);
          assert(y >= 0);
          assert(y < cf->numrows);

          if (cdf->seen_xy[x][y] == 0)
          {
            value = cf->data[x][y].value;
            if (value > 0)
            {
              mean += value;
              n++;
            }
          }
          cdf->seen_xy[x][y] = 1;
        }
      }
      
      if (n)
        mean /= n;

      /* Now shift all values by the same factor */
      memset(cdf->seen_xy[0], 0, cdf->numrows*cdf->numcols*sizeof(affy_uint8));
      for (j = 0; j < number_of_probes; j++)
      {
        x = cdf->probe[j]->pm.x;
        y = cdf->probe[j]->pm.y;

        if (cdf->seen_xy[x][y] == 0)
          cf->data[x][y].value *= target_mean / mean;
        cdf->seen_xy[x][y] = 1;

        /* hack for missing MM probes, where MM coords == PM coords
         */
        if (cdf->probe[j]->pm.x == cdf->probe[j]->mm.x &&
            cdf->probe[j]->pm.y == cdf->probe[j]->mm.y)
        {
          continue;
        }

        x = cdf->probe[j]->mm.x;
        y = cdf->probe[j]->mm.y;

        if (cdf->seen_xy[x][y] == 0)
          cf->data[x][y].value *= target_mean / mean;
        cdf->seen_xy[x][y] = 1;
      }
    }
    /* only PM */
    else
    {
      assert(d->chip[i]->pm != NULL);

      /* deal with duplicate probes */
      if (cdf->dupe_probes_flag)
      {
        memset(cdf->seen_xy[0], 0, cdf->numrows*cdf->numcols*sizeof(affy_uint8));
        for (j = 0; j < number_of_probes; j++)
        {
          x = cdf->probe[j]->pm.x;
          y = cdf->probe[j]->pm.y;

          if (cdf->seen_xy[x][y] == 0)
          {
            value = d->chip[i]->pm[j];
            if (value > 0)
            {
              mean += value;
              n++;
            }
          }
          cdf->seen_xy[x][y] = 1;
        }
      }
      else
      {
        for (j = 0; j < number_of_probes; j++)
        {
          value = d->chip[i]->pm[j];
          if (value > 0)
          {
            mean += value;
            n++;
          }
        }
      }

      if (n)
        mean /= n;

      /* Now shift all values by the same factor */
      for (j = 0; j < number_of_probes; j++)
        d->chip[i]->pm[j] *= target_mean / mean;
    }
  }

  info("done.\n");
}
