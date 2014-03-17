
/**************************************************************************
 *
 * Filename:  write_probe_values.c
 *
 * Purpose:   Write individual probe values from a set of chips to a
 *            given file
 *
 * Creation: 
 *
 * Author:    Andrew M. Hoerter
 *
 * Copyright: Copyright (C) 2009, Moffitt Cancer Center.
 *            All rights reserved.
 *
 * Update History
 * --------------
 * 07/01/09: Created (AMH)
 * 11/20/10: Support printing values < 0.000001 (EAW)
 * 03/17/14: changed a few variable types to be more correct,
 *            should not have any effect on results (EAW)
 *
 **************************************************************************/

#include <affy.h>
#include <utils.h>

void affy_write_probe_values(AFFY_CHIPSET *cs, 
                             char *filename, 
                             int opts,
                             AFFY_ERROR *err)
{
  FILE        *values_file;
  affy_int32   ps_idx, p_idx;
  unsigned int c_idx;

  if ((values_file = fopen(filename, "w")) == NULL)
    AFFY_HANDLE_ERROR_VOID("couldn't open probe values file for writing",
                           AFFY_ERROR_IO,
                           err);

  /* print header */
  fprintf(values_file, "%s\t", filename);
  for (c_idx = 0; c_idx < cs->num_chips; c_idx++)
  {
    fprintf(values_file, "%s", cs->chip[c_idx]->filename);
      
    if (c_idx < (cs->num_chips - 1))
      fprintf(values_file, "\t");
  }

  fprintf(values_file, "\n");

  /* print probe values */
  for (ps_idx = 0; ps_idx < cs->cdf->numprobesets; ps_idx++)
  {
    AFFY_PROBESET *ps = cs->cdf->probeset + ps_idx;

    for (p_idx = 0; p_idx < ps->numprobes; p_idx++)
    {
      affy_int32 x_loc = ps->probe[p_idx].pm.x;
      affy_int32 y_loc = ps->probe[p_idx].pm.y;

      fprintf(values_file, "%s.%" AFFY_PRNd32, ps->name, p_idx);

      for (c_idx = 0; c_idx < cs->num_chips; c_idx++)
      {
        AFFY_CHIP *c = cs->chip[c_idx];
        affy_float32 val;

        if (opts & AFFY_USE_PM)
        {
          val = c->pm[ps->probe[p_idx].index];
        }
        else
        {
          val = c->cel->data[x_loc][y_loc].value;
        }

        if (val < 0.000001)
          fprintf(values_file, "\t%e", val);
        else
          fprintf(values_file, "\t%f", val);
      }

      fprintf(values_file, "\n");
    }
  }
    
  fclose(values_file);
}
