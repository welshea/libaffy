
/**************************************************************************
 *
 * Filename:  write_binary_cel_file.c
 *
 * Purpose:   Write an AFFY_CELFILE structure to a binary on-disk file.
 *
 * Creation:  8 Oct 2010
 *
 * Author:    Andrew M. Hoerter
 *
 * Copyright: Copyright (C) 2010, Moffitt Cancer Center.
 *            All rights reserved.
 *
 * Update History
 * --------------
 * 10/08/10: Initial creation (AMH)
 * 03/10/14: #ifdef out CEL qc fields to save memory (EAW)
 *
 **************************************************************************/

#include <affy.h>

static void write_header_section(FILE *fp, 
                                 AFFY_CHIP *cp,  
                                 AFFY_ERROR *err);
static void write_intensity_section(FILE *fp, 
                                    AFFY_CHIP *cp, 
                                    LIBUTILS_PB_STATE *pbs, 
                                    AFFY_ERROR *err);
static void write_mask_section(FILE *fp, 
                               AFFY_CHIP *cp,  
                               AFFY_ERROR *err);
static void write_outlier_section(FILE *fp, 
                                  AFFY_CHIP *cp,  
                                  AFFY_ERROR *err);

void affy_write_binary_cel_file(FILE *fp, 
                                AFFY_CHIP *cp,
                                AFFY_ERROR *err)
{
  affy_int32        magic = AFFY_CEL_BINARYFILE_MAGIC;
  LIBUTILS_PB_STATE pbs;

  assert(cp->cel             != NULL);
  assert(cp->cdf             != NULL);
  assert(cp->cdf->array_type != NULL);
  assert(fp                  != NULL);

  pb_init(&pbs);

  /* Process by section */
  /* Magic Number */
  if (affy_write32_le(fp, &magic) != 0)
    AFFY_HANDLE_ERROR_GOTO("I/O error writing binary CEL file", 
                           AFFY_ERROR_IO, 
                           err,
                           cleanup);

  write_header_section(fp, cp, err);
  AFFY_CHECK_ERROR_GOTO(err, cleanup);

  write_intensity_section(fp, cp, &pbs, err);
  AFFY_CHECK_ERROR_GOTO(err, cleanup);

  write_mask_section(fp, cp, err);
  AFFY_CHECK_ERROR_GOTO(err, cleanup);

  write_outlier_section(fp, cp, err);
  AFFY_CHECK_ERROR_GOTO(err, cleanup);

  /* No support for subgrids, we don't keep track of them */
cleanup:
  pb_cleanup(&pbs);
}

static void write_header_section(FILE *fp, 
                                 AFFY_CHIP *cp, 
                                 AFFY_ERROR *err)
{
  int         i;
  affy_int32  tmp_ai;
  char        dh[MAXBUF];

  /* Version */
  tmp_ai = 4;
  if (affy_write32_le(fp, &tmp_ai) != 0)
    AFFY_HANDLE_ERROR_VOID("I/O error writing binary CEL file",
                           AFFY_ERROR_IO,
                           err);

  if (affy_write32_le(fp, &cp->cel->numcols) != 0)
    AFFY_HANDLE_ERROR_VOID("I/O error writing binary CEL file",
                           AFFY_ERROR_IO,
                           err);

  if (affy_write32_le(fp, &cp->cel->numrows) != 0)
    AFFY_HANDLE_ERROR_VOID("I/O error writing binary CEL file",
                           AFFY_ERROR_IO,
                           err);

  /* total number of cells */
  tmp_ai = cp->cel->numrows * cp->cel->numcols;
  if (affy_write32_le(fp, &tmp_ai) != 0)
    AFFY_HANDLE_ERROR_VOID("I/O error writing binary CEL file",
                           AFFY_ERROR_IO,
                           err);


  /* XXX: complete hack, do this in a better way some day */
  sprintf(dh, "DatHeader= %s.1sq", cp->cdf->array_type);
  tmp_ai = strlen(dh);

  if (affy_write32_le(fp, &tmp_ai) != 0)
    AFFY_HANDLE_ERROR_VOID("I/O error writing binary CEL file", 
                           AFFY_ERROR_IO, 
                           err);
  if (affy_writechars(fp, dh) != 0)
    AFFY_HANDLE_ERROR_VOID("I/O error writing binary CEL file", 
                           AFFY_ERROR_IO, 
                           err);
  
  /* Fake header strings */
  for (i = 0; i < 2; i++)
  {
    affy_int32 len   = 1;
    char      *empty = "0";
    
    if (affy_write32_le(fp, &len) != 0)
      AFFY_HANDLE_ERROR_VOID("I/O error writing binary CEL file", 
                             AFFY_ERROR_IO, 
                             err);
    if (affy_writechars(fp, empty) != 0)
      AFFY_HANDLE_ERROR_VOID("I/O error writing binary CEL file", 
                             AFFY_ERROR_IO, 
                             err);
    
  }
  
  /* Cell margin */
  tmp_ai = 0;
  if (affy_write32_le(fp, &tmp_ai) != 0)
    AFFY_HANDLE_ERROR_VOID("I/O error writing binary CEL file",
                           AFFY_ERROR_IO,
                           err);
  
  if (affy_write32_le(fp, &cp->cel->numoutliers) != 0)
    AFFY_HANDLE_ERROR_VOID("I/O error writing binary CEL file",
                           AFFY_ERROR_IO,
                           err);
  
  if (affy_write32_le(fp, &cp->cel->nummasks) != 0)
    AFFY_HANDLE_ERROR_VOID("I/O error writing binary CEL file",
                           AFFY_ERROR_IO,
                           err);
    
  /* Number of subgrids */
  tmp_ai = 0;
  if (affy_write32_le(fp, &tmp_ai) != 0)
    AFFY_HANDLE_ERROR_VOID("I/O error writing binary CEL file",
                           AFFY_ERROR_IO,
                           err);
    
}

static void write_intensity_section(FILE *fp, 
                                    AFFY_CHIP *cp, 
                                    LIBUTILS_PB_STATE *pbs, 
                                    AFFY_ERROR *err)
{
  affy_int32  x, y;
  affy_uint32 num_cells;

  num_cells = cp->cel->numrows * cp->cel->numcols;

  pb_begin(pbs, num_cells, "Writing CEL file");

  for (y = 0; y < cp->cel->numrows; y++)
    for (x = 0; x < cp->cel->numcols; x++)
    {
      affy_float32 tmp_f;

#ifndef STORE_CEL_QC
      affy_int16   tmp_i16;
#endif

      tmp_f = cp->cel->data[x][y].value;
      if (affy_write32_le(fp, &tmp_f) != 0)
        AFFY_HANDLE_ERROR_VOID("I/O error writing binary CEL file",
                               AFFY_ERROR_IO,
                               err);

#ifdef STORE_CEL_QC
      tmp_f = cp->cel->data[x][y].stddev;
#else
      tmp_f = 0.0;
#endif
      if (affy_write32_le(fp, &tmp_f) != 0)
        AFFY_HANDLE_ERROR_VOID("I/O error writing binary CEL file",
                               AFFY_ERROR_IO,
                               err);

#ifdef STORE_CEL_QC
      if (affy_write16_le(fp, &(cp->cel->data[x][y].numpixels)) != 0)
        AFFY_HANDLE_ERROR_VOID("I/O error writing binary CEL file",
                               AFFY_ERROR_IO,
                               err);
#else
      tmp_i16 = 1;
      if (affy_write16_le(fp, &tmp_i16) != 0)
        AFFY_HANDLE_ERROR_VOID("I/O error writing binary CEL file",
                               AFFY_ERROR_IO,
                               err);
#endif
      
      pb_tick(pbs, 1, "");
    }

  pb_finish(pbs, "%" AFFY_PRNu32 " cells", num_cells);
}

static void write_mask_section(FILE *fp, 
                               AFFY_CHIP *cp,
                               AFFY_ERROR *err)
{
  affy_int16 x, y;

  if (cp->cel->nummasks == 0)
    return;

  for (x = 0; x < cp->cel->numcols; x++)
    for (y = 0; y < cp->cel->numrows; y++)
      if (affy_ismasked(cp, x, y))
      {
        if (affy_write16_le(fp, &x) != 0)
          AFFY_HANDLE_ERROR_VOID("I/O error writing binary CEL file",
                                 AFFY_ERROR_IO,
                                 err);

        if (affy_write16_le(fp, &y) != 0)
          AFFY_HANDLE_ERROR_VOID("I/O error writing binary CEL file",
                                 AFFY_ERROR_IO,
                                 err);
      }
}

static void write_outlier_section(FILE *fp, 
                                  AFFY_CHIP *cp, 
                                  AFFY_ERROR *err)
{
  affy_int16 x, y;

  if (cp->cel->numoutliers == 0)
    return;

  for (x = 0; x < cp->cel->numcols; x++)
    for (y = 0; y < cp->cel->numrows; y++)
      if (affy_isoutlier(cp, x, y))
      {
        if (affy_write16_le(fp, &x) != 0)
          AFFY_HANDLE_ERROR_VOID("I/O error writing binary CEL file",
                                 AFFY_ERROR_IO,
                                 err);

        if (affy_write16_le(fp, &y) != 0)
          AFFY_HANDLE_ERROR_VOID("I/O error writing binary CEL file",
                                 AFFY_ERROR_IO,
                                 err);
      }
}
