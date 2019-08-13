
/**************************************************************************
 *
 * Filename:  load_calvin_cel_file.c
 *
 * Purpose:   Parse a Calvin CEL file and initialize an accompanying
 *            structure.
 *
 * Creation:  12/02/2008
 *
 * Author:    Andrew M Hoerter
 *
 * Copyright: Copyright (C) 2008, Moffitt Cancer Center.
 *            All rights reserved.
 *
 * Update History
 * --------------
 * 12/02/08: File creation (AMH)
 * 09/20/10: Pooled memory allocator (AMH)
 * 05/09/13: Fixed mask and outlier (harmless) memory allocation typos (EAW)
 * 05/13/13: Check for negative coords in mask/outliers (EAW)
 * 05/13/13: Fix regression due to AFFY_POINT int32 (EAW)
 * 05/13/13: Partial support for salvaging corrupt CEL files (EAW)
 * 03/06/14: fixed row/col memory allocation errors, the dimensions were swapped (EAW)
 * 03/10/14: #ifdef out CEL qc fields to save memory (EAW)
 * 03/17/14: fixed swapped row/col when reading data,
 *            only affected asymmetric chips (EAW)
 *
 **************************************************************************/

#include <affy.h>

static void process_intensity_dataset(AFFY_CALVINIO *cio,
                                      AFFY_CELFILE *cf,
                                      LIBUTILS_PB_STATE *pbs,
                                      AFFY_ERROR *err);
static void process_stddev_dataset(AFFY_CALVINIO *cio,
                                   AFFY_CELFILE *cf,
                                   LIBUTILS_PB_STATE *pbs,
                                   AFFY_ERROR *err);
static void process_cellpixel_dataset(AFFY_CALVINIO *cio,
                                      AFFY_CELFILE *cf,
                                      LIBUTILS_PB_STATE *pbs,
                                      AFFY_ERROR *err);
static void process_mask_dataset(AFFY_CALVINIO *cio,
                                 AFFY_CELFILE *cf,
                                 LIBUTILS_PB_STATE *pbs,
                                 AFFY_ERROR *err);
static void process_outlier_dataset(AFFY_CALVINIO *cio,
                                    AFFY_CELFILE *cf,
                                    LIBUTILS_PB_STATE *pbs,
                                    AFFY_ERROR *err);

static const AFFY_CALVIN_COLUMN_MAPPING point_map[] =
  {
    { "X", offsetof(AFFY_POINT16, x) },
    { "Y", offsetof(AFFY_POINT16, y) },
    { NULL, 0 }
  };

void affy_load_calvin_cel_file(FILE *fp,
                               AFFY_CELFILE *cf,
                               LIBUTILS_PB_STATE *pbs,
                               AFFY_ERROR *err)
{
  AFFY_CALVINIO          *cio = NULL;
  AFFY_CALVIN_FILEHEADER *fh = NULL;
  AFFY_CALVIN_DATAHEADER *dh = NULL;
  AFFY_CALVIN_PARAM      *param;
  affy_int32              i;

  assert(fp != NULL);
  assert(cf != NULL);

  cio = affy_calvinio_init(fp, err);
  AFFY_CHECK_ERROR_GOTO(err, cleanup);
  
  fh = affy_calvin_get_file_metadata(cio, err);
  AFFY_CHECK_ERROR_GOTO(err, cleanup);
  
  dh = affy_calvin_get_dataheader(cio, err);
  AFFY_CHECK_ERROR_GOTO(err, cleanup);

  info("Found Calvin (generic) CEL version: %" AFFY_PRNu8, fh->file_version);

  if ((param = affy_calvin_find_param(dh->params, 
                                      dh->num_params,
                                      "affymetrix-cel-cols")) == NULL)
    AFFY_HANDLE_ERROR_GOTO("CEL column parameter not found", 
                           AFFY_ERROR_BADFORMAT, 
                           err, 
                           cleanup);

  cf->numcols = param->value.int_val;

  if ((param = affy_calvin_find_param(dh->params, 
                                      dh->num_params,
                                      "affymetrix-cel-rows")) == NULL)
    AFFY_HANDLE_ERROR_GOTO("CEL row parameter not found", 
                           AFFY_ERROR_BADFORMAT, 
                           err, 
                           cleanup);

  cf->numrows = param->value.int_val;

  info("CEL Dimensions: %" AFFY_PRNd32 "x%" AFFY_PRNd32,
       cf->numcols,
       cf->numrows);

  cf->data = h_suballoc(cf, cf->numcols * sizeof(AFFY_CELL *));
  if (cf->data == NULL)
    AFFY_HANDLE_ERROR_GOTO("malloc failed", 
                           AFFY_ERROR_OUTOFMEM, 
                           err, 
                           cleanup);

  cf->data[0] = h_subcalloc(cf->data, 
                            cf->numcols * cf->numrows, 
                            sizeof(AFFY_CELL));
  if (cf->data[0] == NULL)
    AFFY_HANDLE_ERROR_GOTO("calloc failed", 
                           AFFY_ERROR_OUTOFMEM, 
                           err, 
                           cleanup);

  cf->mask = h_suballoc(cf, cf->numcols * sizeof(affy_uint8 *));
  if (cf->mask == NULL)
    AFFY_HANDLE_ERROR_GOTO("malloc failed", 
                           AFFY_ERROR_OUTOFMEM, 
                           err, 
                           cleanup);

  cf->mask[0] = h_subcalloc(cf->mask,
                            cf->numcols * cf->numrows,
                            sizeof(affy_uint8));
  if (cf->mask[0] == NULL)
    AFFY_HANDLE_ERROR_GOTO("calloc failed", 
                           AFFY_ERROR_OUTOFMEM, 
                           err, 
                           cleanup);

  cf->outlier = h_suballoc(cf, cf->numcols * sizeof(affy_uint8 *));
  if (cf->outlier == NULL)
    AFFY_HANDLE_ERROR_GOTO("malloc failed", 
                           AFFY_ERROR_OUTOFMEM, 
                           err, 
                           cleanup);

  cf->outlier[0] = h_subcalloc(cf->outlier, 
                               cf->numcols * cf->numrows, 
                               sizeof(affy_uint8));
  if (cf->outlier[0] == NULL)
    AFFY_HANDLE_ERROR_GOTO("calloc failed", 
                           AFFY_ERROR_OUTOFMEM, 
                           err, 
                           cleanup);

  for (i = 1; i < cf->numcols; i++)
  {
    cf->data[i]    = cf->data[0] + (i * cf->numrows);
    cf->mask[i]    = cf->mask[0] + (i * cf->numrows);
    cf->outlier[i] = cf->outlier[0] + (i * cf->numrows);
  }
  
  process_intensity_dataset(cio, cf, pbs, err);
  AFFY_CHECK_ERROR_VOID(err);
  if (cf->corrupt_flag)
    goto cleanup;

#if STORE_CEL_QC
  process_stddev_dataset(cio, cf, pbs, err);
  AFFY_CHECK_ERROR_VOID(err);
  if (cf->corrupt_flag)
    goto cleanup;
#endif

  process_mask_dataset(cio, cf, pbs, err);
  AFFY_CHECK_ERROR_VOID(err);
  if (cf->corrupt_flag)
    goto cleanup;

#if STORE_CEL_QC
  process_cellpixel_dataset(cio, cf, pbs, err);
  AFFY_CHECK_ERROR_VOID(err);
  if (cf->corrupt_flag)
    goto cleanup;
#endif

  process_outlier_dataset(cio, cf, pbs, err);
  AFFY_CHECK_ERROR_VOID(err);
  if (cf->corrupt_flag)
    goto cleanup;

 cleanup:
  affy_calvinio_free(cio);
  affy_free_calvin_fileheader(fh);
  affy_free_calvin_dataheader(dh);
}

static void process_intensity_dataset(AFFY_CALVINIO *cio,
                                      AFFY_CELFILE *cf,
                                      LIBUTILS_PB_STATE *pbs,
                                      AFFY_ERROR *err)
{
  AFFY_CALVIN_DATASET_IO    *dio;
  affy_uint32                ds_index, num_cells;
  AFFY_CALVIN_COLUMN_MAPPING ofs[] = { { "Intensity", 0 }, { NULL, 0 } };

  affy_int32 row, col;

  num_cells = cf->numrows * cf->numcols;

  pb_begin(pbs, num_cells, "Loading intensities");

  ds_index = affy_calvin_find_dataset_index(cio, 0, "Intensity", err);
  if (ds_index == -1)
    AFFY_HANDLE_ERROR_VOID("Intensity dataset not found",
                           AFFY_ERROR_BADFORMAT,
                           err);

  dio = affy_calvin_prepare_dataset(cio, 0, ds_index, err);
  AFFY_CHECK_ERROR_VOID(err);

  for (row = 0; row < cf->numrows; row++)
  {
    for (col = 0; col < cf->numcols; col++)
    {
      affy_float32 val;
      
      affy_calvin_read_dataset_rows(dio,
                                    pbs,
                                    col + (row * cf->numcols),
                                    1,
                                    (void *)(&val),
                                    sizeof(affy_float32),
                                    ofs,
                                    err);
      AFFY_CHECK_ERROR_GOTO(err, cleanup);
      
      cf->data[col][row].value = val;
      
      pb_tick(pbs, 1,"");
    }
  }
  
  pb_finish(pbs, "%" AFFY_PRNu32 " cells", num_cells);

cleanup:
  affy_calvin_close_dataset(dio);
}

#ifdef STORE_CEL_QC
static void process_stddev_dataset(AFFY_CALVINIO *cio,
                                   AFFY_CELFILE *cf,
                                   LIBUTILS_PB_STATE *pbs,
                                   AFFY_ERROR *err)
{
  AFFY_CALVIN_DATASET_IO    *dio;
  affy_uint32                ds_index, num_cells, i;
  AFFY_CALVIN_COLUMN_MAPPING ofs[] = { { "StdDev", 0 }, { NULL, 0 } };

  num_cells = cf->numrows * cf->numcols;

  ds_index = affy_calvin_find_dataset_index(cio, 0, "StdDev", err);
  if (ds_index == -1)
  {
    cf->corrupt_flag = 1;
    fprintf(stderr, "CORRUPT_CEL_FILE: Standard deviation dataset not found: %s\n",
            cf->filename);
    return;
  }

  pb_begin(pbs, num_cells, "Loading standard deviations");

  dio = affy_calvin_prepare_dataset(cio, 0, ds_index, err);
  AFFY_CHECK_ERROR_VOID(err);

  for (i = 0; i < num_cells; i++)
  {
    affy_float32 val;

    affy_calvin_read_dataset_rows(dio,
                                  pbs, 
                                  i, 
                                  1, 
                                  (void *)(&val),
                                  sizeof(affy_float32),
                                  ofs,
                                  err);
    AFFY_CHECK_ERROR_GOTO(err, cleanup);

    cf->data[0][i].stddev = val;
  }
  
  pb_finish(pbs, "%" AFFY_PRNu32 " cells", num_cells);

cleanup:
  affy_calvin_close_dataset(dio);
}

static void process_cellpixel_dataset(AFFY_CALVINIO *cio,
                                      AFFY_CELFILE *cf,
                                      LIBUTILS_PB_STATE *pbs,
                                      AFFY_ERROR *err)
{
  AFFY_CALVIN_DATASET_IO    *dio;
  affy_uint32                ds_index, num_cells;
  AFFY_CALVIN_COLUMN_MAPPING ofs[] = { { "Pixel", offsetof(AFFY_CELL, numpixels) },
                                       { NULL, 0 } };


  num_cells = cf->numrows * cf->numcols;

  ds_index = affy_calvin_find_dataset_index(cio, 0, "Pixel", err);
  if (ds_index == -1)
  {
    cf->corrupt_flag = 1;
    fprintf(stderr, "CORRUPT_CEL_FILE: Pixel dataset not found: %s\n",
            cf->filename);
    return;
  }

  pb_begin(pbs, num_cells, "Loading cell pixel counts");

  dio = affy_calvin_prepare_dataset(cio, 0, ds_index, err);
  AFFY_CHECK_ERROR_VOID(err);

  affy_calvin_read_dataset_rows(dio, 
                                pbs, 
                                0, 
                                num_cells, 
                                (void *)(cf->data[0]),
                                sizeof(AFFY_CELL),
                                ofs,
                                err);
  AFFY_CHECK_ERROR_GOTO(err, cleanup);
  
  pb_finish(pbs, "%" AFFY_PRNu32 " cells", num_cells);

cleanup:
  affy_calvin_close_dataset(dio);
}
#endif    /* STORE_CEL_QC */

static void process_mask_dataset(AFFY_CALVINIO *cio,
                                 AFFY_CELFILE *cf,
                                 LIBUTILS_PB_STATE *pbs,
                                 AFFY_ERROR *err)
{
  AFFY_CALVIN_DATASET_IO    *dio;
  AFFY_POINT16                 ap;
  affy_uint32                ds_index, i, j;
  int                        corrupt_flag = 0;

  ds_index = affy_calvin_find_dataset_index(cio, 0, "Mask", err);
  if (ds_index == -1)
  {
    cf->corrupt_flag = 1;
    fprintf(stderr, "CORRUPT_CEL_FILE: Mask dataset not found: %s\n",
            cf->filename);
    cf->nummasks = 0;
    return;
  }

  dio = affy_calvin_prepare_dataset(cio, 0, ds_index, err);
  AFFY_CHECK_ERROR_VOID(err);

  cf->nummasks = dio->metadata->num_rows;

  pb_begin(pbs, cf->nummasks, "Loading masks");

  j = 0;
  for (i = 0; i < cf->nummasks; i++)
  {
    ap.x = -9999;
    ap.y = -9999;

    affy_calvin_read_dataset_rows(dio, 
                                  pbs, 
                                  i,
                                  1, 
                                  (void *)&ap,
                                  sizeof(AFFY_POINT16),
                                  point_map,
                                  err);
    AFFY_CHECK_ERROR_GOTO(err, cleanup);

    if ((ap.x >= cf->numcols) || (ap.y >= cf->numrows) ||
        (ap.x < 0) || (ap.y < 0))
    {
      if (corrupt_flag == 0)
      {
        cf->corrupt_flag = 1;
        fprintf(stderr, "\nCORRUPT_CEL_FILE: Invalid mask location:");
        fprintf(stderr, " %s", cf->filename);
        fprintf(stderr, " %d", ap.x);
        fprintf(stderr, " %d\n", ap.y);
      }
      
      corrupt_flag = 1;
    }
    
    if (corrupt_flag)
      continue;

    bit_set(cf->mask[ap.x], ap.y);
    
    j++;
  }
  
  cf->nummasks = j;
  
  pb_finish(pbs, "%" AFFY_PRNu32 " masks", cf->nummasks);

cleanup:
  affy_calvin_close_dataset(dio);
}

static void process_outlier_dataset(AFFY_CALVINIO *cio,
                                    AFFY_CELFILE *cf,
                                    LIBUTILS_PB_STATE *pbs,
                                    AFFY_ERROR *err)
{
  AFFY_CALVIN_DATASET_IO    *dio;
  AFFY_POINT16                 ap;
  affy_uint32                ds_index, i, j;
  int                        corrupt_flag = 0;

  ds_index = affy_calvin_find_dataset_index(cio, 0, "Outlier", err);
  if (ds_index == -1)
  {
    cf->corrupt_flag = 1;
    fprintf(stderr, "CORRUPT_CEL_FILE: Outlier dataset not found: %s\n",
            cf->filename);
    cf->numoutliers = 0;
    return;
  }

  dio = affy_calvin_prepare_dataset(cio, 0, ds_index, err);
  AFFY_CHECK_ERROR_VOID(err);

  cf->numoutliers = dio->metadata->num_rows;

  pb_begin(pbs, cf->numoutliers, "Loading outliers");

  j = 0;
  for (i = 0; i < cf->numoutliers; i++)
  {
    ap.x = -9999;
    ap.y = -9999;
  
    affy_calvin_read_dataset_rows(dio, 
                                  pbs, 
                                  i,
                                  1, 
                                  (void *)&ap,
                                  sizeof(AFFY_POINT16),
                                  point_map,
                                  err);
    AFFY_CHECK_ERROR_GOTO(err, cleanup);

    if ((ap.x >= cf->numcols) || (ap.y >= cf->numrows) ||
        (ap.x < 0) || (ap.y < 0))
    {

      if (corrupt_flag == 0)
      {
        cf->corrupt_flag = 1;
        fprintf(stderr, "\nCORRUPT_CEL_FILE: Invalid outlier location:");
        fprintf(stderr, " %s", cf->filename);
        fprintf(stderr, " %d", ap.x);
        fprintf(stderr, " %d\n", ap.y);
      }
      
      corrupt_flag = 1;
    }
    
    if (corrupt_flag)
      continue;

    bit_set(cf->outlier[ap.x], ap.y);
    
    j++;
  }
  
  cf->numoutliers = j;

  pb_finish(pbs, "%" AFFY_PRNu32 " outliers", cf->numoutliers);

cleanup:
  affy_calvin_close_dataset(dio);
}
