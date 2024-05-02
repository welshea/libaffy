
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
 *  1/04/07: Add various assertions for debugging (AMH)
 *  4/14/05: Imported/repaired from old libaffy (AMH)
 *  4/19/05: AFFY_CELFILE now uses AFFY_CELL's (AMH)
 *  3/10/11: Works on both PM/MM and PM-only data, rather than segfault (EAW)
 *  3/06/14: Fixed x/y numcols/numrows assertions (they were swapped) (EAW)
 *  3/06/14: Further fixed missing MM functionality, so that it works correctly now (EAW)
 *  3/06/14: affy_mean_normalization only uses points >0 now (EAW)
 *  3/14/14: fixed to work with exon arrays (EAW)
 * 10/29/18: handle 1:many probe:probeset (EAW)
 *  1/10/24: use geometric mean instead of arithmetic mean (EAW)
 *  4/25/24: add affy_median_normalization() function (EAW)
 *
 **************************************************************************/

#include <affy.h>
#include <utils.h>

/* does not take intensity into account, only "technical" reasons */
char is_masked_probe(AFFY_CDFFILE *cdf, AFFY_CELFILE *cf,
                     int x, int y, int p, AFFY_COMBINED_FLAGS *f)
{
  char mask_char;

  mask_char = (bit_test(cf->mask[x], y)) ||
    (cdf->cell_type[x][y] == AFFY_UNDEFINED_LOCATION) ||
    (cdf->cell_type[x][y] == AFFY_QC_LOCATION);

  /* skip AFFX/control probesets */
  if (affy_is_control_string(cdf->probe[p]->ps->name))
    mask_char = 1;

  /* mask probesets that we want to exclude from training */
  if (f->use_exclusions && cdf->exclusions)
  {
    if (bsearch(&cdf->probe[p]->ps->name, &cdf->exclusions[0],
        cdf->numexclusions, sizeof(char *), compare_string))
    {
      mask_char = 1;
    }
  }
  /* we also want to exclude spikins */
  if (f->use_spikeins && cdf->spikeins)
  {
    if (bsearch(&cdf->probe[p]->ps->name, &cdf->spikeins[0],
        cdf->numspikeins, sizeof(char *), compare_string))
    {
      mask_char = 1;
    }
  }
  
  return mask_char;
}


/* 
 *  affy_mean_normalization()
 *
 *  Normalize the chips to the same constant mean intensity
 */
void affy_mean_normalization(AFFY_CHIPSET *d, double target_mean,
                             AFFY_COMBINED_FLAGS *f)
{
  AFFY_CDFFILE *cdf = d->cdf;
  AFFY_CELFILE *cf;
  double       *mean_array = NULL;
  double        mean, value, min;
  int           number_of_probes;
  int           i, j, x, y, n;
  char          mask_char;

  assert(d   != NULL); 
  assert(cdf != NULL);
  
  mean_array = (double *) calloc(d->num_chips, sizeof(double));

  info("Performing mean normalization...");

  number_of_probes = cdf->numprobes;

  for (i = 0; i < d->num_chips; i++)
  {
    assert(d->chip[i]      != NULL);
    assert(d->chip[i]->cel != NULL);

    cf = d->chip[i]->cel;

    /* find minimum value per chip */
    min = 9.99E99;

    /* both PM and MM */
    if (cf->data)
    {
      for (j = 0; j < number_of_probes; j++)
      {
        assert(cdf->probe[j] != NULL);

        x = cdf->probe[j]->pm.x;
        y = cdf->probe[j]->pm.y;

        assert(x >= 0);
        assert(x < cf->numcols);
        assert(y >= 0);
        assert(y < cf->numrows); 

        mask_char = is_masked_probe(cdf, cf, x, y, j, f);

        /* don't skip < 0 in the min value calculation */
        value = cf->data[x][y].value;
        if (mask_char == 0 && value < min)
          min = value;

        /* hack for missing MM probes, where MM coords == PM coords
         */
        if (cdf->probe[j]->pm.x == cdf->probe[j]->mm.x &&
            cdf->probe[j]->pm.y == cdf->probe[j]->mm.y)
        {
          continue;
        }

        x = cdf->probe[j]->mm.x;
        y = cdf->probe[j]->mm.y;

        assert(x >= 0);
        assert(x < cf->numcols);
        assert(y >= 0);
        assert(y < cf->numrows);

        mask_char = is_masked_probe(cdf, cf, x, y, j, f);

        /* don't skip < 0 in the min value calculation */
        value = cf->data[x][y].value;
        if (mask_char == 0 && value < min)
          min = value;
      }
    }
    /* only PM */
    else
    {
      assert(d->chip[i]->pm != NULL);

      for (j = 0; j < number_of_probes; j++)
      {
        x = cdf->probe[j]->pm.x;
        y = cdf->probe[j]->pm.y;

        assert(x >= 0);
        assert(x < cf->numcols);
        assert(y >= 0);
        assert(y < cf->numrows); 

        mask_char = is_masked_probe(cdf, cf, x, y, j, f);

        /* don't skip < 0 in the min value calculation */
        value = d->chip[i]->pm[j];
        if (mask_char == 0 && value < min)
          min = value;
      }
    }

    /* Calculate mean */
    mean = 0;
    
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
          mask_char = is_masked_probe(cdf, cf, x, y, j, f);

          value = cf->data[x][y].value;
          if (mask_char == 0 && value > 0 &&
              (value > min || f->m_include_min))
          {
            mean += log(value);
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

        x = cdf->probe[j]->mm.x;
        y = cdf->probe[j]->mm.y;

        assert(x >= 0);
        assert(x < cf->numcols);
        assert(y >= 0);
        assert(y < cf->numrows);

        if (cdf->seen_xy[x][y] == 0)
        {
          mask_char = is_masked_probe(cdf, cf, x, y, j, f);

          value = cf->data[x][y].value;
          if (mask_char == 0 && value > 0 &&
              (value > min || f->m_include_min))
          {
            mean += log(value);
            n++;
          }
        }
        cdf->seen_xy[x][y] = 1;
      }
      
      if (n)
        mean /= n;
      
      mean_array[i] = mean;
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
            mask_char = is_masked_probe(cdf, cf, x, y, j, f);

            value = d->chip[i]->pm[j];
            if (mask_char == 0 && value > 0 &&
                (value > min || f->m_include_min))
            {
              mean += log(value);
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
          x = cdf->probe[j]->pm.x;
          y = cdf->probe[j]->pm.y;

          mask_char = is_masked_probe(cdf, cf, x, y, j, f);

          value = d->chip[i]->pm[j];
          if (mask_char == 0 && value > 0 &&
              (value > min || f->m_include_min))
          {
            mean += log(value);
            n++;
          }
        }
      }

      if (n)
        mean /= n;

      mean_array[i] = mean;
    }
  }


  /* if target mean is zero, set target mean to mean of means */
  if (target_mean == 0)
  {
    for (i = 0; i < d->num_chips; i++)
      target_mean += mean_array[i];

    if (d->num_chips)
      target_mean /= d->num_chips;
    target_mean = exp(target_mean);
  }


  for (i = 0; i < d->num_chips; i++)
  {
    assert(d->chip[i]      != NULL);
    assert(d->chip[i]->cel != NULL);

    cf = d->chip[i]->cel;
    
    mean = exp(mean_array[i]);
    
    /* both PM and MM */
    if (cf->data)
    {
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
      /* DO NOT DO THIS, since the pm[] array already has separate copies */
      if (0 && cdf->dupe_probes_flag)
      {
        memset(cdf->seen_xy[0], 0, cdf->numrows*cdf->numcols*sizeof(affy_uint8));
        for (j = 0; j < number_of_probes; j++)
        {
          x = cdf->probe[j]->pm.x;
          y = cdf->probe[j]->pm.y;

          if (cdf->seen_xy[x][y] == 0)
            d->chip[i]->pm[j] *= target_mean / mean;
          cdf->seen_xy[x][y] = 1;
        }
      }
      else
      {
        /* Now shift all values by the same factor */
        /* safe to scale them all, since dupe probes have their own copies */
        for (j = 0; j < number_of_probes; j++)
          d->chip[i]->pm[j] *= target_mean / mean;
      }
    }
  }

  info("done.\n");
  
  if (mean_array)
      free(mean_array);
}


/* 
 *  affy_median_normalization()
 *
 *  Normalize the chips to the same constant median intensity
 */
void affy_median_normalization(AFFY_CHIPSET *d, double target_median,
                             AFFY_COMBINED_FLAGS *f)
{
  AFFY_CDFFILE *cdf = d->cdf;
  AFFY_CELFILE *cf;
  double       *median_array = NULL;
  double       *value_array  = NULL;
  double        median, value, min, min_higher;
  int           number_of_probes;
  int           i, j, x, y, n;
  char          mask_char;

  assert(d   != NULL); 
  assert(cdf != NULL);
  
  info("Performing median normalization...");

  number_of_probes = cdf->numprobes;

  median_array = (double *) calloc(d->num_chips, sizeof(double));
  value_array  = (double *) calloc(2 * number_of_probes, sizeof(double));

  for (i = 0; i < d->num_chips; i++)
  {
    assert(d->chip[i]      != NULL);
    assert(d->chip[i]->cel != NULL);

    cf = d->chip[i]->cel;

    /* find minimum value per chip */
    min        = 9.99E99;
    min_higher = 9.99E99;

    /* both PM and MM */
    if (cf->data)
    {
      for (j = 0; j < number_of_probes; j++)
      {
        assert(cdf->probe[j] != NULL);

        x = cdf->probe[j]->pm.x;
        y = cdf->probe[j]->pm.y;

        assert(x >= 0);
        assert(x < cf->numcols);
        assert(y >= 0);
        assert(y < cf->numrows); 

        mask_char = is_masked_probe(cdf, cf, x, y, j, f);

        /* don't skip < 0 in the min value calculation */
        value = cf->data[x][y].value;
        if (mask_char == 0 && value < min)
          min = value;

        /* hack for missing MM probes, where MM coords == PM coords
         */
        if (cdf->probe[j]->pm.x == cdf->probe[j]->mm.x &&
            cdf->probe[j]->pm.y == cdf->probe[j]->mm.y)
        {
          continue;
        }

        x = cdf->probe[j]->mm.x;
        y = cdf->probe[j]->mm.y;

        assert(x >= 0);
        assert(x < cf->numcols);
        assert(y >= 0);
        assert(y < cf->numrows);

        mask_char = is_masked_probe(cdf, cf, x, y, j, f);

        /* don't skip < 0 in the min value calculation */
        value = cf->data[x][y].value;
        if (mask_char == 0 && value < min)
          min = value;
      }
    }
    /* only PM */
    else
    {
      assert(d->chip[i]->pm != NULL);

      for (j = 0; j < number_of_probes; j++)
      {
        x = cdf->probe[j]->pm.x;
        y = cdf->probe[j]->pm.y;

        assert(x >= 0);
        assert(x < cf->numcols);
        assert(y >= 0);
        assert(y < cf->numrows); 

        mask_char = is_masked_probe(cdf, cf, x, y, j, f);

        /* don't skip < 0 in the min value calculation */
        value = d->chip[i]->pm[j];
        if (mask_char == 0 && value < min)
          min = value;
      }
    }

    /* Calculate median */
    median = min;
    
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
          mask_char = is_masked_probe(cdf, cf, x, y, j, f);

          value = cf->data[x][y].value;
          if (mask_char == 0 && value > 0 &&
              (value > min || f->m_include_min))
          {
            value_array[n++] = value;
            
            if (value > min && value < min_higher)
              min_higher = value;
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

        x = cdf->probe[j]->mm.x;
        y = cdf->probe[j]->mm.y;

        assert(x >= 0);
        assert(x < cf->numcols);
        assert(y >= 0);
        assert(y < cf->numrows);

        if (cdf->seen_xy[x][y] == 0)
        {
          mask_char = is_masked_probe(cdf, cf, x, y, j, f);

          value = cf->data[x][y].value;
          if (mask_char == 0 && value > 0 &&
              (value > min || f->m_include_min))
          {
            value_array[n++] = value;

            if (value > min && value < min_higher)
              min_higher = value;
          }
        }
        cdf->seen_xy[x][y] = 1;
      }
      
      median = affy_median(value_array, n, f);

      /* HACK -- min value isn't trustworthy, use next-higher value */
      if (median == min && min_higher != 9.99E99)
      {
          median = min_higher;
      }

      median_array[i] = median;
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
            mask_char = is_masked_probe(cdf, cf, x, y, j, f);

            value = d->chip[i]->pm[j];
            if (mask_char == 0 && value > 0 &&
                (value > min || f->m_include_min))
            {
              value_array[n++] = value;

              if (value > min && value < min_higher)
                min_higher = value;
            }
          }
          cdf->seen_xy[x][y] = 1;
        }
      }
      else
      {
        for (j = 0; j < number_of_probes; j++)
        {
          x = cdf->probe[j]->pm.x;
          y = cdf->probe[j]->pm.y;

          mask_char = is_masked_probe(cdf, cf, x, y, j, f);

          value = d->chip[i]->pm[j];
          if (mask_char == 0 && value > 0 &&
              (value > min || f->m_include_min))
          {
            value_array[n++] = value;

            if (value > min && value < min_higher)
              min_higher = value;
          }
        }
      }

      median = affy_median(value_array, n, f);

      /* HACK -- min value isn't trustworthy, use next-higher value */
      if (median == min && min_higher != 9.99E99)
      {
          median = min_higher;
      }

      median_array[i] = median;
    }
  }


  /* if target median is zero, set target median to geometric mean of medians */
  if (target_median == 0)
  {
    for (i = 0; i < d->num_chips; i++)
    {
      target_median += log(median_array[i]);
    }

    if (d->num_chips)
      target_median /= d->num_chips;
    target_median = exp(target_median);
  }


  for (i = 0; i < d->num_chips; i++)
  {
    assert(d->chip[i]      != NULL);
    assert(d->chip[i]->cel != NULL);

    cf = d->chip[i]->cel;
    
    median = median_array[i];
    
    /* both PM and MM */
    if (cf->data)
    {
      /* Now shift all values by the same factor */
      memset(cdf->seen_xy[0], 0, cdf->numrows*cdf->numcols*sizeof(affy_uint8));
      for (j = 0; j < number_of_probes; j++)
      {
        x = cdf->probe[j]->pm.x;
        y = cdf->probe[j]->pm.y;

        if (cdf->seen_xy[x][y] == 0)
          cf->data[x][y].value *= target_median / median;
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
          cf->data[x][y].value *= target_median / median;
        cdf->seen_xy[x][y] = 1;
      }
    }
    /* only PM */
    else
    {
      assert(d->chip[i]->pm != NULL);

      /* deal with duplicate probes */
      /* DO NOT DO THIS, since the pm[] array already has separate copies */
      if (0 && cdf->dupe_probes_flag)
      {
        memset(cdf->seen_xy[0], 0, cdf->numrows*cdf->numcols*sizeof(affy_uint8));
        for (j = 0; j < number_of_probes; j++)
        {
          x = cdf->probe[j]->pm.x;
          y = cdf->probe[j]->pm.y;

          if (cdf->seen_xy[x][y] == 0)
            d->chip[i]->pm[j] *= target_median / median;
          cdf->seen_xy[x][y] = 1;
        }
      }
      else
      {
        /* Now shift all values by the same factor */
        /* safe to scale them all, since dupe probes have their own copies */
        for (j = 0; j < number_of_probes; j++)
          d->chip[i]->pm[j] *= target_median / median;
      }
    }
  }

  info("done.\n");
  
  if (median_array)
      free(median_array);
  if (value_array)
      free(value_array);
}
