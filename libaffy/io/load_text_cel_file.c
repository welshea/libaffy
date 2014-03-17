
/**************************************************************************
 *
 * Filename:  load_text_cel_file.c
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
 * 01/23/08: Add checking for short reads (AMH)
 * 04/04/08: New error handling scheme, misc cleanups (AMH)
 * 09/20/10: Pooled memory allocator (AMH)
 * 05/13/13: Check for negative coords in mask/outliers (EAW)
 * 03/10/14: #ifdef out CEL qc fields to save memory (EAW)
 * 03/17/14: fixed row/col memory allocation errors, the dimensions were swapped (EAW)
 *
 **************************************************************************/

#include <affy.h>
#include <utils.h>

static void process_cel_section(AFFY_TEXTIO *tf, 
                                AFFY_CELFILE *cf,
                                LIBUTILS_PB_STATE *pbs,
                                AFFY_ERROR *err);
static void process_header_section(AFFY_TEXTIO *tf, 
                                   AFFY_CELFILE *cf,
                                   AFFY_ERROR *err);
static void process_intensity_section(AFFY_TEXTIO *tf, 
                                      AFFY_CELFILE *cf,
                                      LIBUTILS_PB_STATE *pbs,
                                      AFFY_ERROR *err);
static void process_mask_section(AFFY_TEXTIO *tf, 
                                 AFFY_CELFILE *cf,
                                 LIBUTILS_PB_STATE *pbs,
                                 AFFY_ERROR *err);
static void process_outlier_section(AFFY_TEXTIO *tf, 
                                    AFFY_CELFILE *cf,
                                    LIBUTILS_PB_STATE *pbs,
                                    AFFY_ERROR *err);
static void process_modified_section(AFFY_TEXTIO *tf, 
                                     AFFY_CELFILE *cf, 
                                     AFFY_ERROR *err);

void affy_load_text_cel_file(FILE *fp, 
                             AFFY_CELFILE *cf,
                             LIBUTILS_PB_STATE *pbs,
                             AFFY_ERROR *err)
{
  char             *str;
  AFFY_TEXTIO      *tf;

  assert(fp != NULL);
  assert(cf != NULL);

  tf = affy_textio_init(fp, err);
  AFFY_CHECK_ERROR_VOID(err);

  while ((str = affy_textio_get_next_line(tf)) != NULL)
  {
    /* Process by section */
    if (STREQ(str, "[CEL]"))
      process_cel_section(tf, cf, pbs, err);
    else if (STREQ(str, "[HEADER]"))
      process_header_section(tf, cf, err);
    else if (STREQ(str, "[INTENSITY]"))
      process_intensity_section(tf, cf, pbs, err);
    else if (STREQ(str, "[MASKS]"))
      process_mask_section(tf, cf, pbs, err);
    else if (STREQ(str, "[OUTLIERS]"))
      process_outlier_section(tf, cf, pbs, err);
    else if (STREQ(str, "[MODIFIED]"))
      process_modified_section(tf, cf, err);
    /* This case should probably not occur */
    else
    {
      info("(Skipping unknown section '%s'.)", str);
      affy_textio_skip_to_next_header(tf);
    }

    AFFY_CHECK_ERROR_GOTO(err, out);
  }

out:
  affy_textio_free(tf);
}

/*
  We have just read the [CEL] line. Process the remainder of the
  section until we get to eof or another []. 
*/
static void process_cel_section(AFFY_TEXTIO *tf, 
                                AFFY_CELFILE *cf,
                                LIBUTILS_PB_STATE *pbs,
                                AFFY_ERROR *err)
{
  char  *s;
  char  *kv[2];

  assert(tf != NULL);
  assert(cf != NULL);

  /* Read lines until eof or a new section */
  while ((s = affy_textio_get_next_line(tf)) != NULL)
  {
    if (*s == '[')
    {
      affy_textio_unget_next_line(tf);
      return;
    }

    if (split(s, kv, '=', 2) != 2)
      AFFY_HANDLE_ERROR_VOID("error parsing CEL section",
                             AFFY_ERROR_BADFORMAT,
                             err);

    if (STREQ(kv[0], "Version"))
      info("Found ASCII CEL version:  %s", kv[1]);

  }
}

static void process_header_section(AFFY_TEXTIO *tf, 
                                   AFFY_CELFILE *cf,
                                   AFFY_ERROR *err)
                                  
{
  char *s, *kv[2];
  int   i;

  assert(tf != NULL);
  assert(cf != NULL);

  /* Read lines until eof or a new section */
  while ((s = affy_textio_get_next_line(tf)) != NULL)
  {
    char *err_str = NULL;

    if (*s == '[')
    {
      affy_textio_unget_next_line(tf);
      break;
    }

    if (split(s, kv, '=', 2) != 2)
      AFFY_HANDLE_ERROR_VOID("error parsing CEL header section",
                             AFFY_ERROR_BADFORMAT,
                             err);

 
   if (STREQ(kv[0], "Cols"))
     cf->numcols = strtol(kv[1], &err_str, 10);
   else if (STREQ(kv[0], "Rows"))
     cf->numrows = strtol(kv[1], &err_str, 10);

   if (err_str == kv[1])
     AFFY_HANDLE_ERROR_VOID("error parsing CEL header section",
                            AFFY_ERROR_BADFORMAT,
                            err);
  }

  if ((cf->numcols <= 0) || (cf->numrows <= 0))
    AFFY_HANDLE_ERROR_VOID("invalid CEL file dimensions",
                           AFFY_ERROR_BADFORMAT,
                           err);

  cf->data = h_suballoc(cf, cf->numcols * sizeof(AFFY_CELL *));
  if (cf->data == NULL)
    AFFY_HANDLE_ERROR_VOID("malloc failed", AFFY_ERROR_OUTOFMEM, err);

  cf->data[0] = h_subcalloc(cf->data, 
                            cf->numcols * cf->numrows, 
                            sizeof(AFFY_CELL));
  if (cf->data[0] == NULL)
    AFFY_HANDLE_ERROR_VOID("calloc failed", AFFY_ERROR_OUTOFMEM, err);

  cf->mask = h_suballoc(cf, cf->numcols * sizeof(bitstr_t *));
  if (cf->mask == NULL)
    AFFY_HANDLE_ERROR_VOID("malloc failed", AFFY_ERROR_OUTOFMEM, err);

  cf->mask[0] = h_subcalloc(cf->mask, cf->numcols * cf->numrows, sizeof(bitstr_t));
  if (cf->mask[0] == NULL)
    AFFY_HANDLE_ERROR_VOID("calloc failed", AFFY_ERROR_OUTOFMEM, err);

  cf->outlier = h_suballoc(cf, cf->numcols * sizeof(bitstr_t *));
  if (cf->outlier == NULL)
    AFFY_HANDLE_ERROR_VOID("malloc failed", AFFY_ERROR_OUTOFMEM, err);

  cf->outlier[0] = h_subcalloc(cf->outlier, 
                               cf->numcols * cf->numrows, 
                               sizeof(bitstr_t));
  if (cf->outlier[0] == NULL)
    AFFY_HANDLE_ERROR_VOID("calloc failed", AFFY_ERROR_OUTOFMEM, err);

  for (i = 1; i < cf->numcols; i++)
  {
    cf->data[i]    = cf->data[0] + (i * cf->numrows);
    cf->mask[i]    = cf->mask[0] + (i * cf->numrows);
    cf->outlier[i] = cf->outlier[0] + (i * cf->numrows);
  }

  info("CEL Dimensions: %" AFFY_PRNd32 "x%" AFFY_PRNd32, 
       cf->numcols, 
       cf->numrows);
}

static void process_intensity_section(AFFY_TEXTIO *tf, 
                                      AFFY_CELFILE *cf,
                                      LIBUTILS_PB_STATE *pbs,
                                      AFFY_ERROR *err)
{
  char      *s, *kv[2];
  bool       read_cellheader = false;
  affy_int32 x, y, npixels, num_read = 0;
  double     val, stdv;

  assert(tf != NULL);
  assert(cf != NULL);

  pb_begin(pbs, cf->numrows * cf->numcols, "Loading intensities");

  /* Read lines until eof or a new section */
  while ((s = affy_textio_get_next_line(tf)) != NULL)
  {
    if (*s == '[')
    {
      affy_textio_unget_next_line(tf);
      break;
    }

    /* Only do streq's in case we haven't read the header */
    if (!read_cellheader)
    {
      if (split(s, kv, '=', 2) != 2)
        AFFY_HANDLE_ERROR_VOID("error parsing CEL intensity section",
                               AFFY_ERROR_BADFORMAT,
                               err);

      if (STREQ(kv[0], "NumberCells"))
        ;
      else if (STREQ(kv[0], "CellHeader"))
        read_cellheader = true;
    }
    else
    {
      /* Otherwise, this is a line of x,y coordinates and mean intensity */
      num_read++;
      pb_tick(pbs, 1,"");

      if (sscanf(s, "%" AFFY_PRNd32 " %" AFFY_PRNd32 " %lf %lf %" AFFY_PRNd32, 
                 &x, 
                 &y, 
                 &val, 
                 &stdv, 
                 &npixels) != 5)
        AFFY_HANDLE_ERROR_VOID("error parsing CEL intensity section",
                               AFFY_ERROR_BADFORMAT,
                               err);

      if ((x >= cf->numcols) || (y >= cf->numrows))
        AFFY_HANDLE_ERROR_VOID("Invalid intensity location",
                               AFFY_ERROR_BADFORMAT,
                               err);

      cf->data[x][y].value     = val;

#ifdef STORE_CEL_QC
      cf->data[x][y].stddev    = stdv;
      cf->data[x][y].numpixels = npixels;
#endif
    }
  }

  if (num_read < (cf->numrows * cf->numcols))
    AFFY_HANDLE_ERROR_VOID("truncated intensity section in CEL file",
                           AFFY_ERROR_BADFORMAT,
                           err);

  pb_finish(pbs, "%" AFFY_PRNd32 " cells", num_read);
}

static void process_mask_section(AFFY_TEXTIO *tf, 
                                 AFFY_CELFILE *cf,
                                 LIBUTILS_PB_STATE *pbs,
                                 AFFY_ERROR *err)
{
  char       *s, *kv[2];
  bool        read_maskheader = false;
  affy_int32  x, y;
  affy_uint32 num_masks = 0;

  assert(tf != NULL);
  assert(cf != NULL);

  /* Read lines until eof or a new section */
  while ((s = affy_textio_get_next_line(tf)) != NULL)
  {
    if (*s == '[')
    {
      affy_textio_unget_next_line(tf);
      break;
    }

    /* Only do streq's in case we haven't read the header */
    if (!read_maskheader)
    {
      if (split(s, kv, '=', 2) != 2)
        AFFY_HANDLE_ERROR_VOID("error parsing CEL mask section",
                               AFFY_ERROR_BADFORMAT,
                               err);

      if (STREQ(kv[0], "NumberCells"))
      {
        char *err_str;

        cf->nummasks = strtol(kv[1], &err_str, 10);
        if (err_str == kv[1])
          AFFY_HANDLE_ERROR_VOID("error parsing CEL mask section",
                                 AFFY_ERROR_BADFORMAT,
                                 err);
          
        pb_begin(pbs, cf->nummasks, "Loading masks");
      }
      else if (STREQ(kv[0], "CellHeader"))
      {
        read_maskheader = true;
      }
    }
    else
    {
      /* Otherwise, this is a line of x,y coordinates and mean intensity */
      if (sscanf(s, "%" AFFY_PRNd32 " %" AFFY_PRNd32, &x, &y) != 2)
        AFFY_HANDLE_ERROR_VOID("error parsing CEL mask section",
                               AFFY_ERROR_BADFORMAT,
                               err);

      if ((x >= cf->numcols) || (y >= cf->numrows) || (x < 0) || (y < 0))
        AFFY_HANDLE_ERROR_VOID("Invalid mask location",
                               AFFY_ERROR_BADFORMAT,
                               err);

      bit_set(cf->mask[x], y);
      num_masks++;

      pb_tick(pbs, 1,"");
    }
  }

  if (num_masks != cf->nummasks)
    warn("Mismatch on number of masks: %d actual, %d expected",
          num_masks, 
          cf->nummasks);

  pb_finish(pbs, "%" AFFY_PRNu32 " masks", num_masks);
}

static void process_outlier_section(AFFY_TEXTIO *tf, 
                                    AFFY_CELFILE *cf,
                                    LIBUTILS_PB_STATE *pbs,
                                    AFFY_ERROR *err)
{
  char       *s, *kv[2];
  bool        read_outlierheader = false;
  affy_int32  x, y;
  affy_uint32 num_outliers = 0;

  assert(tf != NULL);
  assert(cf != NULL);

  /* Read lines until eof or a new section */
  while ((s = affy_textio_get_next_line(tf)) != NULL)
  {
    if (*s == '[')
    {
      affy_textio_unget_next_line(tf);
      break;
    }

    /* Only do streq's in case we haven't read the header */
    if (!read_outlierheader)
    {
      if (split(s, kv, '=', 2) != 2)
        AFFY_HANDLE_ERROR_VOID("error parsing CEL outlier section",
                               AFFY_ERROR_BADFORMAT,
                               err);

      if (STREQ(kv[0], "NumberCells"))
      {
        char *err_str;

        cf->numoutliers = strtol(kv[1], &err_str, 10);
        if (err_str == kv[1])
          AFFY_HANDLE_ERROR_VOID("error parsing CEL outlier section",
                                 AFFY_ERROR_BADFORMAT,
                                 err);

        pb_begin(pbs, cf->numoutliers, "Loading outliers");
      }
      else if (STREQ(kv[0], "CellHeader"))
      {
        read_outlierheader = true;
      }
    }
    else
    {
      /* Otherwise, this is a line of x,y coordinates and mean intensity */
      if (sscanf(s, "%" AFFY_PRNd32 " %" AFFY_PRNd32, &x, &y) != 2)
        AFFY_HANDLE_ERROR_VOID("error parsing CEL outlier section",
                               AFFY_ERROR_BADFORMAT,
                               err);

      if ((x >= cf->numcols) || (y >= cf->numrows) || (x < 0) || (y < 0))
        AFFY_HANDLE_ERROR_VOID("Invalid outlier location",
                               AFFY_ERROR_BADFORMAT,
                               err);

      bit_set(cf->outlier[x], y);
      num_outliers++;

      pb_tick(pbs, 1,"");
    }
  }

  if (num_outliers != cf->numoutliers)
    warn("Mismatch on number of outliers: %" AFFY_PRNu32 
         " actual, %" AFFY_PRNu32 " expected",
         num_outliers, 
         cf->numoutliers);

  pb_finish(pbs, "%" AFFY_PRNu32 " outliers", num_outliers);
}

/* 
 * This section is supposed to be outdated post MAS 4, so it should be
 *  empty and ignored.
 */
static void process_modified_section(AFFY_TEXTIO *tf, 
                                     AFFY_CELFILE *cf, 
                                     AFFY_ERROR *err)
{
  assert(tf != NULL);
  assert(cf != NULL);

  affy_textio_skip_to_next_header(tf);
}
