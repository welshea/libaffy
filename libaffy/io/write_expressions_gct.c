
/**************************************************************************
 *
 * Filename:  write_expressions_gct.c
 *
 * Purpose:   Write expression data in GCT format from a set of chips to a 
 *            given file.  GCT is used by GenePattern pipelines.  For
 *            documentation on the GCT format, see the following URL:
 *   
 *         http://www.broad.mit.edu/cancer/software/genepattern/tutorial/
 *         gp_fileformats_text.html#gct
 *
 * Creation:  09/04/2007
 *
 * Author:    Andrew Hoerter
 *
 * Copyright: Copyright (C) 2007, Moffitt Cancer Center.
 *            All rights reserved.
 *
 * Update History
 * --------------
 * 09/04/07: File creation, copied from write_expressions.c (AMH)
 * 03/14/08: New error handling scheme (AMH)
 *
 **************************************************************************/

#include <affy.h>
#include <utils.h>

void affy_write_expressions_gct(AFFY_CHIPSET *c, 
                                char *filename, 
                                AFFY_ERROR *err)
{
  unsigned int n, i;
  FILE        *fp;

  assert(c        != NULL);
  assert(filename != NULL);

  fp = fopen(filename, "w");
  if (fp == NULL)
    AFFY_HANDLE_ERROR_VOID("couldn't open output file", AFFY_ERROR_IO, err);

  /* First, the version string, always the same. */
  if (fprintf(fp, "#1.2\t%s\n", c->cdf->array_type) < 0)
  {
    fclose(fp);
    AFFY_HANDLE_ERROR_VOID("I/O error writing expressions", 
                           AFFY_ERROR_IO, 
                           err);
  }

  /* 
   * Next, the number of rows (probesets/genes) and columns 
   * (samples/CEL files) 
   */
  if (fprintf(fp, "%d\t%d\n", c->cdf->numprobesets, c->num_chips) < 0)
  {
    fclose(fp);
    AFFY_HANDLE_ERROR_VOID("I/O error writing expressions", 
                           AFFY_ERROR_IO, 
                           err);
  }

  /* 
   * Next a list of the sample identifiers.  In this case we will simply
   * reuse the code from write_expressions() and print the CEL filename
   * stems, which seems reasonable.
   */
  if (fprintf(fp, "Name\tDescription\t") < 0)
  {
    fclose(fp);
    AFFY_HANDLE_ERROR_VOID("I/O error writing expressions", 
                           AFFY_ERROR_IO, 
                           err);
  }

  for (i = 0; i < c->num_chips; i++) 
  {
    char *filestem = stem_from_filename_safer(c->chip[i]->filename);
    if (filestem == NULL)
    {
      fclose(fp);
      AFFY_HANDLE_ERROR_VOID("stem_from_filename_safer failed", 
                             AFFY_ERROR_OUTOFMEM,
                             err);
    }
	       
    if (fprintf(fp, "%s%c", filestem, (i < c->num_chips - 1) ? '\t' : '\n') < 0)
    {
      fclose(fp);
      free(filestem);
      AFFY_HANDLE_ERROR_VOID("I/O error writing expressions", 
                             AFFY_ERROR_IO, 
                             err);
    }
    
    free(filestem);
  }
  
  /* Now the data. */
  for (i = 0; i < c->cdf->numprobesets; i++)
  {
    char *ps_name = c->cdf->probeset[i].name;

    /* 
     * Note the difference from write_expressions() here; in GCT there
     * is a second column with a description, it isn't used by any
     * programs though so we just repeat the probeset name.
     */ 
    if (fprintf(fp, "%s\t%s", ps_name, ps_name) < 0)
    {
      fclose(fp);
      AFFY_HANDLE_ERROR_VOID("I/O error writing expressions", 
                             AFFY_ERROR_IO, 
                             err);
    }

    for (n = 0; n < c->num_chips; n++)
    {
      if (fprintf(fp, "\t%f", c->chip[n]->probe_set[i]) < 0)
      {
        fclose(fp);
        AFFY_HANDLE_ERROR_VOID("I/O error writing expressions", 
                               AFFY_ERROR_IO, 
                               err);
        
      }
    }

    if (fprintf(fp, "\n") < 0)
    {
      fclose(fp);
      AFFY_HANDLE_ERROR_VOID("I/O error writing expressions", 
                             AFFY_ERROR_IO, 
                             err);
    }
  }
  
  fclose(fp);
}
