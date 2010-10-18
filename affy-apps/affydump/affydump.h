/**************************************************************************
 *
 * Filename:  affydump.h
 *
 * Purpose:   Various definitions used by affydump.
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

#ifndef _AFFYDUMP_H
# define _AFFYDUMP_H

#if defined(AFFY_HAVE_NETCDF)
# include <netcdf.h>
#endif

struct write_format
{
  const char *input_type;
  const char *output_type;
  void      (*writer) (void *, const char *);
};

void cel_to_json(void *vp, const char *output_name);
void cdf_to_json(void *vp, const char *output_name);
void dat_to_json(void *vp, const char *output_name);

void cel_to_sexpr(void *vp, const char *output_name);
void cdf_to_sexpr(void *vp, const char *output_name);
void dat_to_sexpr(void *vp, const char *output_name);

#if defined(AFFY_HAVE_NETCDF)
void cel_to_netcdf(void *vp, const char *output_name);
void cdf_to_netcdf(void *vp, const char *output_name);
void dat_to_netcdf(void *vp, const char *output_name);
#endif

#endif
