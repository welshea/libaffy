
/**************************************************************************
 *
 * Filename:  dump_cel.c
 *
 * Purpose:   Dump out CEL file data in various formats.
 *
 * Creation:  4/20/2008 
 *
 * Author:    Andrew M Hoerter
 *
 * Copyright: Copyright (C) 2009, Moffitt Cancer Center.
 *            All rights reserved.
 *
 * Update History
 * --------------
 * 04/20/09: File creation (AMH)
 * 03/10/14: #ifdef out CEL qc fields to save memory (EAW)
 *
 **************************************************************************/

#include <affy.h>

#include "affydump.h"

#define JSON_BOOLEAN(x) ((x) ? "true" : "false")

static void print_cell_json(AFFY_CELL *cell, 
                            int ismasked, 
                            int isoutlier,
                            FILE *fp)
{
#ifdef STORE_CEL_QC
  const char *tmpl = "{\"intensity\":%f,\"stddev\":%f,"
                     "\"numpixels\":%" AFFY_PRNd16 ","
                     "\"masked\":%s,\"outlier\":%s}";
#else
  const char *tmpl = "{\"intensity\":%f,"
                     "\"masked\":%s,\"outlier\":%s}";
#endif

  assert(cell != NULL);
  assert(fp   != NULL);

#ifdef STORE_CEL_QC
  fprintf(fp, tmpl, cell->value, cell->stddev, cell->numpixels,
          JSON_BOOLEAN(ismasked), JSON_BOOLEAN(isoutlier));
#else
  fprintf(fp, tmpl, cell->value,
          JSON_BOOLEAN(ismasked), JSON_BOOLEAN(isoutlier));
#endif
}

void cel_to_json(void *vp, const char *output_name)
{
  affy_int32    i, j;
  AFFY_CELFILE *cf = (AFFY_CELFILE *)vp;
  FILE         *fp;

  assert(cf          != NULL);
  assert(output_name != NULL);

  fp = FOPEN(output_name, "w+");

  fprintf(fp, "{\"orig_filename\":\"%s\",", cf->filename);
  fprintf(fp, "\"numrows\":%" AFFY_PRNd32 ",", cf->numrows);
  fprintf(fp, "\"numcols\":%" AFFY_PRNd32 ",", cf->numcols);
  fprintf(fp, "\"nummasks\":%" AFFY_PRNu32 ",", cf->nummasks);
  fprintf(fp, "\"numoutliers\":%" AFFY_PRNu32 ",", cf->numoutliers);
  fprintf(fp, "\"cells\":[\n");
  
  for (i = 0; i < cf->numrows; i++)
    for (j = 0; j < cf->numcols; j++)
    {
      if ((i > 0) || (j > 0))
        fprintf(fp, ",\n    ");
      else
        fprintf(fp, "    ");
        
      print_cell_json(&(cf->data[j][i]), 
                      bit_test(cf->mask[j], i),
                      bit_test(cf->outlier[j], i),
                      fp);
    }
  
  fprintf(fp, "]}");
}

void cel_to_sexpr(void *vp, const char *output_name)
{
}

#if defined(AFFY_HAVE_NETCDF)
static void write_netcdf_bitmap(bitstr_t **bitmap, 
                                int ncid, 
                                int varid,
                                int rows,
                                int cols)
{
  unsigned int rec_count = 0, i, j;
  int          status;

  assert(bitmap != NULL);

  for (i = 0; i < rows; i++)
    for (j = 0; j < cols; j++)
    {
      if (bit_test(bitmap[j], i))
      {
        int count[] = { 1, 2 };
        int val[2], start[2];
        
        val[0] = j;
        val[1] = i;
        
        start[0] = rec_count++;
        start[1] = 0;
        
        if ((status = nc_put_vara_int(ncid, 
                                      varid, 
                                      start, 
                                      count, 
                                      val)) != NC_NOERR)
          die("Couldn't write NetCDF data (%s)", nc_strerror(status));
      }
    }
}

void cel_to_netcdf(void *vp, const char *output_name)
{
  AFFY_CELFILE *cf = (AFFY_CELFILE *)vp;
  int           ncid, status, i, j;
  int           cell_row_dim, cell_col_dim, xy_dim, record_dim;
  int           intensity_var, mask_var, outlier_var, stddev_var, 
                numpixels_var;
  int           cell_coord_dims[2], bitmap_dims[2];
  
  assert(output_name != NULL);
  assert(cf          != NULL);
  
  if (nc_create(output_name, 0, &ncid) != NC_NOERR)
    die("Couldn't write file %s", output_name);

  /* Dimensions */
  if ((status = nc_def_dim(ncid, 
                           "cell_row", 
                           cf->numrows, 
                           &cell_row_dim)) != NC_NOERR)
    die("Couldn't define NetCDF dimension (%s)", nc_strerror(status));

  if ((status = nc_def_dim(ncid, 
                           "cell_col", 
                           cf->numcols, 
                           &cell_col_dim)) != NC_NOERR)
    die("Couldn't define NetCDF dimension (%s)", nc_strerror(status));

  if ((status = nc_def_dim(ncid, 
                           "cell_rowcol", 
                           2, 
                           &xy_dim)) != NC_NOERR)
    die("Couldn't define NetCDF dimension (%s)", nc_strerror(status));

  if ((status = nc_def_dim(ncid, 
                           "record", 
                           NC_UNLIMITED, 
                           &record_dim)) != NC_NOERR)
    die("Couldn't define NetCDF dimension (%s)", nc_strerror(status));

  cell_coord_dims[0] = cell_row_dim;
  cell_coord_dims[1] = cell_col_dim;

  bitmap_dims[0] = record_dim;
  bitmap_dims[1] = xy_dim;

  /* Variables */
  if ((status = nc_def_var(ncid, 
                           "intensity", 
                           NC_DOUBLE, 
                           2, 
                           cell_coord_dims, 
                           &intensity_var)) != NC_NOERR)
    die("Couldn't define NetCDF variable (%s)", nc_strerror(status));

  if ((status = nc_def_var(ncid, 
                           "standard_deviation", 
                           NC_DOUBLE, 
                           2, 
                           cell_coord_dims, 
                           &stddev_var)) != NC_NOERR)
    die("Couldn't define NetCDF variable (%s)", nc_strerror(status));

  if ((status = nc_def_var(ncid, 
                           "number_of_pixels", 
                           NC_SHORT, 
                           2, 
                           cell_coord_dims, 
                           &numpixels_var)) != NC_NOERR)
    die("Couldn't define NetCDF variable (%s)", nc_strerror(status));

  if ((status = nc_def_var(ncid, 
                           "mask_coords", 
                           NC_INT, 
                           2, 
                           bitmap_dims, 
                           &mask_var)) != NC_NOERR)
    die("Couldn't define NetCDF variable (%s)", nc_strerror(status));

  if ((status = nc_def_var(ncid, 
                           "outlier_coords", 
                           NC_INT, 
                           2, 
                           bitmap_dims, 
                           &outlier_var)) != NC_NOERR)
    die("Couldn't define NetCDF variable (%s)", nc_strerror(status));

  /* Global attributes */
  if ((status = nc_put_att_text(ncid, 
                                NC_GLOBAL, 
                                "original_filename",
                                strlen(cf->filename),
                                cf->filename)) != NC_NOERR)
    die("Couldn't set NetCDF attribute (%s)", nc_strerror(status));

  /* Begin data phase */
  if ((status = nc_enddef(ncid)) != NC_NOERR)
    die("Couldn't leave NetCDF define mode (%s)", nc_strerror(status));

  /* Masks */
  write_netcdf_bitmap(cf->mask, 
                      ncid, 
                      mask_var, 
                      cf->numrows, 
                      cf->numcols);

  /* Outliers */
  write_netcdf_bitmap(cf->outlier, 
                      ncid, 
                      outlier_var, 
                      cf->numrows, 
                      cf->numcols);

  /* Intensity, standard deviation, number of pixels */
  for (i = 0; i < cf->numrows; i++)
    for (j = 0; j < cf->numcols; j++)
    {
      int        idx[2];
      double     intensity, stddev;
      affy_int16 numpixels;

      intensity = cf->data[j][i].value;
      stddev    = cf->data[j][i].stddev;
      numpixels = cf->data[j][i].numpixels;

      idx[0] = i;
      idx[1] = j;

      if ((status = nc_put_var1_double(ncid, 
                                       intensity_var, 
                                       idx, 
                                       &intensity)) != NC_NOERR)
        die("Couldn't write NetCDF data (%s)", nc_strerror(status));

      if ((status = nc_put_var1_short(ncid, 
                                      numpixels_var, 
                                      idx, 
                                      &numpixels)) != NC_NOERR)
        die("Couldn't write NetCDF data (%s)", nc_strerror(status));

      if ((status = nc_put_var1_double(ncid, 
                                       stddev_var, 
                                       idx, 
                                       &stddev)) != NC_NOERR)
        die("Couldn't write NetCDF data (%s)", nc_strerror(status));
    }

  /* All done, close the file */
  if (nc_close(ncid) != NC_NOERR)
    die("Couldn't close file %s", output_name);
}
#endif

