
/**************************************************************************
 *
 * Filename:  load_binary_cel_file.c
 *
 * Purpose:   Parse a CEL file and initialize an accompanying structure.
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
 * 04/19/05: Now we use AFFY_CELL and collect all cell data (AMH)
 * 10/03/07: Minor beautification/cleanups (AMH)
 * 01/25/08: Convert to use readmulti() (AMH)
 * 03/12/08: New error handling scheme (AMH)
 * 09/20/10: Pooled memory allocator (AMH)
 * 05/13/13: Check for negative coords in mask/outliers (EAW)
 * 05/13/13: Partial support for salvaging corrupt CEL files (EAW)
 * 03/10/14: #ifdef out CEL qc fields to save memory (EAW)
 * 03/17/14: fixed row/col memory allocation errors, the dimensions were swapped (EAW)
 *
 **************************************************************************/

#include <affy.h>

static void process_header_section(FILE *fp, 
                                   AFFY_CELFILE *cf, 
                                   AFFY_ERROR *err);
static void process_intensity_section(FILE *fp, 
                                      AFFY_CELFILE *cf,
                                      LIBUTILS_PB_STATE *pbs,
                                      AFFY_ERROR *err);
static void process_mask_section(FILE *fp, 
                                 AFFY_CELFILE *cf, 
                                 LIBUTILS_PB_STATE *pbs,
                                 AFFY_ERROR *err);
static void process_outlier_section(FILE *fp, 
                                    AFFY_CELFILE *cf, 
                                    LIBUTILS_PB_STATE *pbs,
                                    AFFY_ERROR *err);

void affy_load_binary_cel_file(FILE *fp, 
                               AFFY_CELFILE *cf,
                               LIBUTILS_PB_STATE *pbs,
                               AFFY_ERROR *err)
{
  affy_int32        magic;

  /* Process by section */
  /* Magic Number */
  if (affy_read32_le(fp, &magic) != 0)
    AFFY_HANDLE_ERROR_VOID("I/O error in binary CEL file", AFFY_ERROR_IO, err);

  if (magic != AFFY_CEL_BINARYFILE_MAGIC)
    AFFY_HANDLE_ERROR_VOID("Bad magic in binary CEL file", 
                           AFFY_ERROR_BADFORMAT, err);

  process_header_section(fp, cf, err);
  AFFY_CHECK_ERROR_VOID(err);
  if (cf->corrupt_flag)
    return;

  process_intensity_section(fp, cf, pbs, err);
  AFFY_CHECK_ERROR_VOID(err);
  if (cf->corrupt_flag)
    return;

  process_mask_section(fp, cf, pbs, err);
  AFFY_CHECK_ERROR_VOID(err);
  if (cf->corrupt_flag)
    return;

  process_outlier_section(fp, cf, pbs, err);
  AFFY_CHECK_ERROR_VOID(err);
  if (cf->corrupt_flag)
    return;

  /* Skip the sub-grid section */
}

/* 
 * Read an integer from the filestream, and seek ahead using it as an
 * offset.
 */
static void read_offset_and_skip(FILE *fp, AFFY_ERROR *err)
{
  affy_int32 ofs;

  assert(fp != NULL);

  if (affy_read32_le(fp, (void *)&ofs) != 0)
    AFFY_HANDLE_ERROR_VOID("I/O error", AFFY_ERROR_IO, err);

  if (fseek(fp, ofs, SEEK_CUR) != 0)
    AFFY_HANDLE_ERROR_VOID("fseek failed", AFFY_ERROR_IO, err);
}

static void process_header_section(FILE *fp, AFFY_CELFILE *cf, AFFY_ERROR *err)
{
  int             i;
  affy_int32      version;
      
  if (affy_readmulti(fp, "%3dl%x",
                     (void *)&version,
                     (void *)&cf->numcols,
                     (void *)&cf->numrows,
                     4) <= 0)
    AFFY_HANDLE_ERROR_VOID("I/O error in CEL header section", 
                           AFFY_ERROR_IO, 
                           err);
  
  info("Found XDA (binary) CEL version: %" AFFY_PRNd32, version);
       
  cf->data = h_suballoc(cf, cf->numcols * sizeof(AFFY_CELL *));
  if (cf->data == NULL)
    AFFY_HANDLE_ERROR_VOID("malloc failed", AFFY_ERROR_OUTOFMEM, err);

  cf->data[0] = h_subcalloc(cf->data, 
                            cf->numcols * cf->numrows, 
                            sizeof(AFFY_CELL));
  if (cf->data[0] == NULL)
    AFFY_HANDLE_ERROR_VOID("calloc failed", AFFY_ERROR_OUTOFMEM, err);

  cf->mask = h_suballoc(cf, cf->numcols * sizeof(affy_uint8 *));
  if (cf->mask == NULL)
    AFFY_HANDLE_ERROR_VOID("malloc failed", AFFY_ERROR_OUTOFMEM, err);

  cf->mask[0] = h_subcalloc(cf->mask, cf->numrows * cf->numcols, sizeof(affy_uint8));
  if (cf->mask[0] == NULL)
    AFFY_HANDLE_ERROR_VOID("calloc failed", AFFY_ERROR_OUTOFMEM, err);

  cf->outlier = h_suballoc(cf, cf->numcols * sizeof(affy_uint8 *));
  if (cf->outlier == NULL)
    AFFY_HANDLE_ERROR_VOID("malloc failed", AFFY_ERROR_OUTOFMEM, err);

  cf->outlier[0] = h_subcalloc(cf->outlier, 
                               cf->numrows * cf->numcols, 
                               sizeof(affy_uint8));
  if (cf->outlier[0] == NULL)
    AFFY_HANDLE_ERROR_VOID("calloc failed", AFFY_ERROR_OUTOFMEM, err);

  for (i = 1; i < cf->numcols; i++)
  {
    cf->data[i]    = cf->data[0] + (i * cf->numrows);
    cf->mask[i]    = cf->mask[0] + (i * cf->numrows);
    cf->outlier[i] = cf->outlier[0] + (i * cf->numrows);
  }

  for (i = 0; i < 3; i++)
  {
    read_offset_and_skip(fp, err);
    AFFY_CHECK_ERROR_VOID(err);
  }

  if (affy_readmulti(fp, "%x%2dl%x",
                     4,
                     (void *)&cf->numoutliers,
                     (void *)&cf->nummasks,
                     4) <= 0)
    AFFY_HANDLE_ERROR_VOID("I/O error in CEL header section", 
                           AFFY_ERROR_IO, 
                           err);
  
  info("CEL Dimensions: %" AFFY_PRNd32 "x%" AFFY_PRNd32, 
       cf->numcols, 
       cf->numrows);
}

static void process_intensity_section(FILE *fp, 
                                      AFFY_CELFILE *cf,
                                      LIBUTILS_PB_STATE *pbs,
                                      AFFY_ERROR *err)
{
  affy_int32 x, y, i, num_cells;

#ifndef STORE_CEL_QC
  double     stddev;
  affy_int16 numpixels;
#endif

  assert(fp != NULL);
  assert(cf != NULL);

  num_cells = cf->numrows * cf->numcols;

  x = 0;
  y = 0;

  pb_begin(pbs, num_cells, "Loading intensities");
        
  for (i = 0; i < num_cells; i++) 
  {
    /* Check for next row */
    if (x == cf->numcols) 
    {
      x = 0;
      y++;
    }
    
#ifdef STORE_CEL_QC
    if (affy_readmulti(fp, "%2Dl%hl", 
                       (void *)&(cf->data[x][y].value),
                       (void *)&(cf->data[x][y].stddev),
                       (void *)&(cf->data[x][y].numpixels)) <= 0)
#else
    if (affy_readmulti(fp, "%2Dl%hl", 
                       (void *)&(cf->data[x][y].value),
                       (void *)&(stddev),
                       (void *)&(numpixels)) <= 0)
#endif
      AFFY_HANDLE_ERROR_VOID("I/O error in CEL intensity section", 
                             AFFY_ERROR_IO, 
                             err);

    x++;

    pb_tick(pbs, 1, "");
  }

  pb_finish(pbs, "%" AFFY_PRNd32 " cells", num_cells);
}

static void process_mask_section(FILE *fp, 
                                 AFFY_CELFILE *cf,
                                 LIBUTILS_PB_STATE *pbs,
                                 AFFY_ERROR *err)
{
  affy_uint32 i, j;
  affy_int16  x, y;
  int corrupt_flag = 0;

  pb_begin(pbs, cf->nummasks, "Loading masks");

  j = 0;
  for (i = 0; i < cf->nummasks; i++) 
  {
    if (affy_readmulti(fp, "%2hl", (void *)&x, (void *)&y) <= 0)
    {
      if (corrupt_flag == 0)
      {
        cf->corrupt_flag = 1;
        fprintf(stderr, "\nCORRUPT_CEL_FILE: I/O error in CEL mask section:");
        fprintf(stderr, " %s\n", cf->filename);
      }
    
      corrupt_flag = 1;
    }

    if (corrupt_flag)
      continue;

    if ((x >= cf->numcols) || (y >= cf->numrows) || (x < 0) || (y < 0))
    {
      if (corrupt_flag == 0)
      {
        cf->corrupt_flag = 1;
        fprintf(stderr, "\nCORRUPT_CEL_FILE: Invalid mask location:");
        fprintf(stderr, " %s", cf->filename);
        fprintf(stderr, " %d", x);
        fprintf(stderr, " %d\n", y);
      }

      corrupt_flag = 1;
    }

    if (corrupt_flag)
      continue;

    bit_set(cf->mask[x], y);

    pb_tick(pbs, 1, "");
    
    j++;
  }
  
  cf->nummasks = j;

  pb_finish(pbs, "%" AFFY_PRNu32 " masks", cf->nummasks);
}

static void process_outlier_section(FILE *fp, 
                                    AFFY_CELFILE *cf,
                                    LIBUTILS_PB_STATE *pbs,
                                    AFFY_ERROR *err)
{
  affy_uint32 i, j;
  affy_int16  x, y;
  int corrupt_flag = 0;
        
  pb_begin(pbs, cf->numoutliers, "Loading outliers");

  j = 0;
  for (i = 0; i < cf->numoutliers; i++) 
  {
    if (affy_readmulti(fp, "%2hl", (void *)&x, (void *)&y) <= 0)
    {
      if (corrupt_flag == 0)
      {
        cf->corrupt_flag = 1;
        fprintf(stderr, "\nCORRUPT_CEL_FILE: I/O error in CEL outlier section:");
        fprintf(stderr, " %s\n", cf->filename);
      }
    
      corrupt_flag = 1;
    }

    if (corrupt_flag)
      continue;

    if ((x >= cf->numcols) || (y >= cf->numrows) || (x < 0) || (y < 0))
    {
      if (corrupt_flag == 0)
      {
        cf->corrupt_flag = 1;
        fprintf(stderr, "\nCORRUPT_CEL_FILE: Invalid outlier location:");
        fprintf(stderr, " %s", cf->filename);
        fprintf(stderr, " %d", x);
        fprintf(stderr, " %d\n", y);
      }

      corrupt_flag = 1;
    }

    if (corrupt_flag)
      continue;

    bit_set(cf->outlier[x], y);

    pb_tick(pbs, 1,"");
    
    j++;
  }
  
  cf->numoutliers = j;

  pb_finish(pbs, "%" AFFY_PRNu32 " outliers", cf->numoutliers);
}
