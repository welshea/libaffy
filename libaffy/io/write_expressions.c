
/**************************************************************************
 *
 * Filename:  write_expressions.c
 *
 * Purpose:   Write expression data from a set of chips to a given file
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
 * 04/08/05: Imported/repaired from old libaffy (AMH)
 * 09/04/07: Factor out filestem code into libutils (AMH)
 * 03/14/08: New error handling scheme (AMH)
 * 09/16/10: Add support for P/A calls (EAW)
 * 10/19/10: Add unlog option (AMH)
 *
 **************************************************************************/

#include "affy.h"
#include "utils.h"
#include "affy_mas5.h"

void affy_write_expressions(AFFY_CHIPSET *c, 
                            char         *filename, 
                            unsigned int  opts,
                            AFFY_ERROR   *err)
{
  unsigned int n, i, print_pa, unlog_flag, log_flag;
  int          ret;
  FILE        *fp;
  double       log2 = log(2.0);

  assert(filename != NULL);
  assert(c        != NULL);

  print_pa   = opts & AFFY_WRITE_EXPR_PA;
  unlog_flag = opts & AFFY_WRITE_EXPR_UNLOG;
  log_flag   = opts & AFFY_WRITE_EXPR_LOG;

  fp = fopen(filename, "w");
  if (fp == NULL)
    AFFY_HANDLE_ERROR_VOID("couldn't open output file", AFFY_ERROR_IO, err);

  /* Print the header first */
  if (fprintf(fp, "%s\t", filename) < 0)
  {
    fclose(fp);
    AFFY_HANDLE_ERROR_VOID("I/O error writing expressions", 
                           AFFY_ERROR_IO, 
                           err);
  }

  for (i = 0; i < c->num_chips; i++) 
  {
    char *filestem = stem_from_filename_safer(c->chip[i]->filename);
    char  sep;

    if (filestem == NULL)
    {
      fclose(fp);
      AFFY_HANDLE_ERROR_VOID("stem_from_filename_safer failed", 
                             AFFY_ERROR_OUTOFMEM,
                             err);
    }
  
    sep = (i < c->num_chips - 1) ? '\t' : '\n';

    if (print_pa)
    {
      ret = fprintf(fp, "%s_EXTR\t%s_CALL\t%s_PVAL%c", 
                    filestem,
                    filestem,
                    filestem,
                    sep);
    }
    else
    {
      ret = fprintf(fp, "%s%c", filestem, sep);
    }

    if (ret < 0)
    {
      fclose(fp);
      free(filestem);
      AFFY_HANDLE_ERROR_VOID("I/O error writing expressions", 
                             AFFY_ERROR_IO, 
                             err);
    }

    free(filestem);
  }
  
  /* Now the data */
  for (i = 0; i < c->cdf->numprobesets; i++)
  {
    /* checking for write errors once per row should be quite sufficient */
    if (fprintf(fp, "%s", c->cdf->probeset[i].name) < 0)
      goto err;

    for (n = 0; n < c->num_chips; n++) 
    {
      double data = c->chip[n]->probe_set[i];
      
      /* preserve missing data */
      if (data)
      {
        if (unlog_flag && !log_flag)
          data = pow(2.0, data);
        else if (log_flag && !unlog_flag)
          data = log(data) / log2;
      }
      
      if (fprintf(fp, "\t%f", data) < 0)
        goto err;

      if (print_pa)
      {
        double pvalue = c->chip[n]->probe_set_call_pvalue[i];

        if (fprintf(fp, "\t%c", affy_mas5_pvalue_call(pvalue)) < 0)
          goto err;
        if (fprintf(fp, "\t%e", pvalue) < 0)
          goto err; 
      }
    }

    if (fprintf(fp, "\n") < 0)
      goto err;
  }
  
  fclose(fp);

  return;

err:
  fclose(fp);
  AFFY_HANDLE_ERROR_VOID("I/O error writing expressions", AFFY_ERROR_IO, err);
}

