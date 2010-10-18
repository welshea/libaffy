/**************************************************************************
 *
 * Filename:  write_pixel_region.c
 *
 * Purpose:   Write pixels from a DAT to a given file. See function at
 *            bottom of file for attempt at generic interface.
 *
 * Creation:  07/20/05
 *
 * Author:    Steven Eschrich
 *
 * Copyright: Copyright (C) 2007, Moffitt Cancer Center.
 *            All rights reserved.
 *
 * Update History
 * --------------
 * 07/20/05: Creation date
 * 04/04/08: New error handling scheme, other cleanups (AMH)
 * 09/20/10: Pooled memory allocator (AMH)
 *
 **************************************************************************/

#include "affy.h"

#ifdef AFFY_TIFF_SUPPORT
# include <tiffio.h>
#endif

void affy_pixregion2tiff(AFFY_PIXREGION *p, char *filename, AFFY_ERROR *err)
{
#ifdef AFFY_TIFF_SUPPORT
  int     i, j;
  TIFF   *fp;
  uint16 *x = NULL;

  /* Suppress noise to stderr */
  (void)TIFFSetErrorHandler(NULL);

  fp = TIFFOpen(filename, "w");
	
  if (fp) 
  {
    x = h_calloc(p->numcols, sizeof(uint16));
    if (x == NULL)
      AFFY_HANDLE_ERROR_GOTO("calloc failed", 
                             AFFY_ERROR_OUTOFMEM, 
                             err, 
                             cleanup);
    
    TIFFSetField(fp, TIFFTAG_BITSPERSAMPLE, 16);
    TIFFSetField(fp, TIFFTAG_IMAGELENGTH, p->numrows);
    TIFFSetField(fp, TIFFTAG_IMAGEWIDTH, p->numcols);
    TIFFSetField(fp, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);

    for (i = 0; i < p->numrows; i++) 
    {
      for (j = 0; j < p->numcols; j++) 
        x[j] = p->data[i][j];

      if (TIFFWriteScanline(fp, x, i, 1) == -1)
        AFFY_HANDLE_ERROR_GOTO("I/O error writing TIFF file", 
                               AFFY_ERROR_IO,
                               err,
                               cleanup);
    }
  }
  else
  {
    AFFY_HANDLE_ERROR_VOID("couldn't open file for TIFF output", 
                           AFFY_ERROR_SYSPERM, 
                           err);
  }

cleanup:
  TIFFClose(fp);
  h_free(x);	
#else
  AFFY_HANDLE_ERROR_VOID("no TIFF support available", AFFY_ERROR_NOTSUPP, err);
#endif
}

void affy_pixregion2text(AFFY_PIXREGION *p, char *filename, AFFY_ERROR *err)
{
  FILE *fp;
  int   i, j;
	
  fp = fopen(filename,"w");
	
  if (fp) 
  {
    for (i = 0; i < p->numrows; i++) 
    {
      for (j = 0; j < p->numcols; j++) 
      {
        fprintf(fp, "%d", p->data[i][j]);

        if (j < p->numcols - 1) 
          fprintf(fp, "\t");
      }
      
      fprintf(fp, "\n");
      
      if (ferror(fp) != 0)
        AFFY_HANDLE_ERROR_GOTO("I/O error writing pixel region",
                               AFFY_ERROR_IO,
                               err,
                               cleanup);
    }    
  }
  else
  {
    AFFY_HANDLE_ERROR_VOID("couldn't open file for output", 
                           AFFY_ERROR_SYSPERM, 
                           err);
  }

cleanup:
  fclose(fp);
}

/*
 * This function could be used as a default preference of what
 * image writer to use at compile time. Could be more sophisticated.
 */
void affy_write_pixel_region(AFFY_PIXREGION *pr, 
                             char *filename, 
                             AFFY_ERROR *err)
{
#ifdef AFFY_TIFF_SUPPORT
  affy_pixregion2tiff(pr, filename, err);
#else
  affy_pixregion2text(pr, filename, err);
#endif
}
