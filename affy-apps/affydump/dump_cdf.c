
/**************************************************************************
 *
 * Filename:  dump_cdf.c
 *
 * Purpose:   Dump out CDF file data in various formats.
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
 *
 **************************************************************************/

#include <affy.h>

#include "affydump.h"

/*
 * JSON routines.
 */

#define JSON_BOOLEAN(x) ((x) ? "true" : "false")

static void print_point_json(AFFY_POINT *p, FILE *fp)
{
  const char *tmpl = "{\"x\":%" AFFY_PRNd16 
                     ",\"y\":%" AFFY_PRNd16 "}";

  assert(p  != NULL);
  assert(fp != NULL);
  
  fprintf(fp, tmpl, p->x, p->y);
}

static void print_probe_json(AFFY_PROBE *p, FILE *fp)
{
  assert(p  != NULL);
  assert(fp != NULL);
  
  fprintf(fp, "{\"pm_loc\":");
  print_point_json(&(p->pm), fp);
  fprintf(fp, ",\"mm_loc\":");
  print_point_json(&(p->mm), fp);
  fprintf(fp, "}");
}

static void print_probeset_json(AFFY_PROBESET *ps, FILE *fp)
{
  int i;

  assert(ps != NULL);
  assert(fp != NULL);
  
  fprintf(fp, "{\"name\":\"%s\",\"numprobes\":%d,\"probes\":[",
          ps->name, ps->numprobes);
  
  for (i = 0; i < ps->numprobes; i++)
  {
    if (i > 0)
      fprintf(fp, ",");

    print_probe_json(&(ps->probe[i]), fp);
  }

  fprintf(fp, "]}");
}

void cdf_to_json(void *vp, const char *output_name)
{
  affy_int32    i;
  affy_uint32   n;
  AFFY_CDFFILE *cdf  = (AFFY_CDFFILE *)vp;
  FILE         *fp;
  const char   *tmpl = "{\"array_type\":\"%s\""
                       ",\"numrows\":%" AFFY_PRNu32
                       ",\"numcols\":%" AFFY_PRNu32
                       ",\"numprobes\":%" AFFY_PRNd32
                       ",\"numprobesets\":%" AFFY_PRNd32
                       ",\"numqcunits\":%" AFFY_PRNd32
                       ",\"probesets\":[\n";

  assert(cdf         != NULL);
  assert(output_name != NULL);

  fp = FOPEN(output_name, "w+");

  fprintf(fp, tmpl, cdf->array_type, cdf->numrows, cdf->numcols, 
          cdf->numprobes, cdf->numprobesets, cdf->numqcunits);

  for (i = 0; i < cdf->numprobesets; i++)
  {
    if (i > 0)
      fprintf(fp, ",\n    ");
    
    print_probeset_json(&(cdf->probeset[i]), fp);
  }
  
  fprintf(fp, "],\"cell_type\":[");
  
  for (n = 0; n < (cdf->numrows * cdf->numcols); n++)
  {
    if (n > 0)
      fprintf(fp, ",");

    fprintf(fp, "%" AFFY_PRNu8, cdf->cell_type[0][n]);
  }

  fprintf(fp, "]\n}\n");
}

/*
 * s-expression routines.
 */

void cdf_to_sexpr(void *vp, const char *output_name)
{
}

#if defined(AFFY_HAVE_NETCDF)
/*
 * NetCDF routines.
 */
void cdf_to_netcdf(void *vp, const char *output_name)
{
  AFFY_CDFFILE *cdf = (AFFY_CDFFILE *)vp;
  int           ncid, status;
  affy_int32    n, probecount;
  int           i;
  int           mm_pm_dim, probe_dim, psname_dim;
  int           psname_var, mm_pm_var;
  int           psname_dims[2], mm_pm_dims[2];
  
  assert(output_name != NULL);
  assert(cdf         != NULL);
  
  if (nc_create(output_name, 0, &ncid) != NC_NOERR)
    die("Couldn't write file %s", output_name);

  /* Dimensions */
  if ((status = nc_def_dim(ncid, 
                           "probeset_name", 
                           40,
                           &psname_dim)) != NC_NOERR)
    die("Couldn't define NetCDF dimension (%s)", nc_strerror(status));

  if ((status = nc_def_dim(ncid,
                           "mm_pm_dim", 
                           4, 
                           &mm_pm_dim)) != NC_NOERR)
    die("Couldn't define NetCDF dimension (%s)", nc_strerror(status));

  if ((status = nc_def_dim(ncid, 
                           "probe_id", 
                           NC_UNLIMITED, 
                           &probe_dim)) != NC_NOERR)
    die("Couldn't define NetCDF dimension (%s)", nc_strerror(status));

  mm_pm_dims[0] = probe_dim;
  mm_pm_dims[1] = mm_pm_dim;

  psname_dims[0] = probe_dim;
  psname_dims[1] = psname_dim;

  /* Variables */
  if ((status = nc_def_var(ncid, 
                           "probeset_name", 
                           NC_CHAR, 
                           2, 
                           psname_dims, 
                           &psname_var)) != NC_NOERR)
    die("Couldn't define NetCDF variable (%s)", nc_strerror(status));

  if ((status = nc_def_var(ncid, 
                           "mm_pm_location", 
                           NC_SHORT, 
                           2, 
                           mm_pm_dims, 
                           &mm_pm_var)) != NC_NOERR)
    die("Couldn't define NetCDF variable (%s)", nc_strerror(status));

  /* Global attributes */
/*   if ((status = nc_put_att_text(ncid,  */
/*                                 NC_GLOBAL,  */
/*                                 "array_type", */
/*                                 strlen(cdf->array_type), */
/*                                 cf->array_type)) != NC_NOERR) */
/*     die("Couldn't set NetCDF attribute (%s)", nc_strerror(status)); */

  /* Begin data phase */
  if ((status = nc_enddef(ncid)) != NC_NOERR)
    die("Couldn't leave NetCDF define mode (%s)", nc_strerror(status));

  /* Process each probe for each probeset */
  for (n = 0, probecount = 0; n < cdf->numprobesets; n++)
  {
    AFFY_PROBESET *ps;
    AFFY_PROBE    *ps_probes;

    ps        = &(cdf->probeset[n]);
    ps_probes = ps->probe;

    for (i = 0; i < ps->numprobes; i++)
    {
      affy_int16 loc_data[4];
      int        loc_start[2];
      int        loc_count[] = { 1, 4 };

      char psname_data[40];
      int  psname_start[2];
      int  psname_count[2];

      loc_data[0] = ps_probes[i].mm.x;
      loc_data[1] = ps_probes[i].mm.y;
      loc_data[2] = ps_probes[i].pm.x;
      loc_data[3] = ps_probes[i].pm.y;

      loc_start[0] = psname_start[0] = probecount++;
      loc_start[1] = 0;

      if ((status = nc_put_vara_short(ncid, 
                                      mm_pm_var, 
                                      loc_start,
                                      loc_count,
                                      loc_data)) != NC_NOERR)
        die("Couldn't write NetCDF data (%s)", nc_strerror(status));

      psname_start[1] = 0;
      psname_count[0] = 1;
      psname_count[1] = strlen(ps->name) + 1;
      strncpy(psname_data, ps->name, 40);

      if ((status = nc_put_vara_text(ncid, 
                                     psname_var, 
                                     psname_start,
                                     psname_count,
                                     psname_data)) != NC_NOERR)
        die("Couldn't write NetCDF data (%s)", nc_strerror(status));
    }

    
  }

  /* All done, close the file */
  if (nc_close(ncid) != NC_NOERR)
    die("Couldn't close file %s", output_name);
}

#endif

