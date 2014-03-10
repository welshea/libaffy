
/**************************************************************************
 *
 * Filename:  affydump.c
 *
 * Purpose:   Command-line driver for Affymetrix data extraction.
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
 * 03/10/14: abort cleanly when netcdf is requested but not available (EAW)
 *
 **************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include <affy.h>
#include <argp.h>

#include "affydump.h"

static struct write_format format_map[] =
  {
    { "cdf", "json", cdf_to_json },
    { "cel", "json", cel_to_json },
    { "cdf", "sexpr", cdf_to_sexpr },
    { "cel", "sexpr", cel_to_sexpr },
#if defined(AFFY_HAVE_NETCDF)
    { "cdf", "netcdf", cdf_to_netcdf },
    { "cel", "netcdf", cel_to_netcdf },
#endif
    { 0 }
  };

static char *output_format   = NULL;
static char *input_format    = NULL;
static char *output_filename = NULL;
static char *input_filename  = NULL;

const char *argp_program_version=affy_version;
const char *argp_program_bug_address="<Eric.Welsh@moffitt.org>";
bool opt_salvage_corrupt = false;

static struct argp_option options[] = {
  { "input-format", 't', "type", 0, 
    "Specify input type, one of `cel', `cdf'" },
  { "output-format", 'f', "type", 0,
    "Specify output format, one of `json', `netcdf' (where available)" },
  { "output-file", 'o', "filename", 0,
    "Specify output filename" },
  { "salvage", 24, 0, 0,
    "Attempt to salvage corrupt CEL files (may still result in corrupt data!)" },
  { 0 }
};

static error_t parse_opt(int key, char *arg, struct argp_state *state);
static struct argp argp = { options, parse_opt, 0, 
			    "affydump - Affymetrix data extraction"};

int main(int argc, char **argv)
{
  FILE                *infile;
  AFFY_CDFFILE        *cdf;
  AFFY_CELFILE        *cel;
  AFFY_ERROR          *err = NULL;
  void                *wrdata;
  struct write_format *wf;

  err = affy_get_default_error();

  argp_parse(&argp, argc, argv, 0, 0, 0);

  infile = FOPEN(input_filename, "r");
  
  /* Make a reasonable guess at input format if not specified */
  if (input_format == NULL)
  {
    if (endsWith(input_filename, ".cdf"))
    {
      printf("Assuming CDF file as input\n");
      input_format = "cdf";
    }
    else if (endsWith(input_filename, ".cel"))
    {
      printf("Assuming CEL file as input\n");
      input_format = "cel";
    }
    else
    {
      fprintf(stderr, "Please specify input type, couldn't guess it\n");

      if (err)
        free(err);

      exit(EXIT_FAILURE);
    }
  }

#if !defined(AFFY_HAVE_NETCDF)
fprintf(stderr, "Error: netcdf output unavailable (AFFY_HAVE_NETCDF is undefined)\n");
exit(AFFY_ERROR_NOTSUPP);
#endif

  /* Get the Affy data */
  if (strcmp(input_format, "cel") == 0)
  {
    cel = affy_load_cel_file(input_filename, err);
    wrdata = (void *)cel;

    if (cel->corrupt_flag && opt_salvage_corrupt == false)
      AFFY_HANDLE_ERROR_GOTO("corrupt CEL file",
                             AFFY_ERROR_BADFORMAT,
                             err,
                             cleanup);
  }
  else if (strcmp(input_format, "cdf") == 0)
  {
    cdf = affy_load_cdf_file_byname(input_filename, NULL, err);
    wrdata = (void *)cdf;
  }

  /* Look for an appropriate conversion function */
  for (wf = format_map; wf != NULL; wf++)
  {
    if ((strcmp(wf->input_type, input_format) == 0)
        && (strcmp(wf->output_type, output_format) == 0))
    {
      (wf->writer)(wrdata, output_filename);

      if (err)
        free(err);

      exit(EXIT_SUCCESS);
    }
  }
  
  /* If we got here, didn't find a usable function */
  fprintf(stderr, "Couldn't convert '%s' to '%s'\n", 
          input_format, output_format);

cleanup:

  if (err)
    free(err);

  exit(EXIT_FAILURE);
}

static error_t parse_opt(int key, char *arg, struct argp_state *state)
{
  switch (key)
  {
    case 24:
      opt_salvage_corrupt = true;
      break;
    case 'f':
      output_format = strdup(arg);
      break;
    case 't':
      input_format = strdup(arg);
      break;
    case 'o':
      output_filename = strdup(arg);
      break;
    case ARGP_KEY_ARG:
      /* consume at most one input filename */
      if (state->arg_num == 0)
        input_filename = strdup(arg);
      else
        argp_usage(state);
      break;
    case ARGP_KEY_END:
      if (state->arg_num != 1)
        argp_error(state, "exactly one input file required\n");
      if (output_format == NULL)
        argp_error(state, "an output format must be specified");
      if (output_filename == NULL)
        argp_error(state, "an output filename must be specified");
      break;
    default:
      return (ARGP_ERR_UNKNOWN);
  }
      
  return (0);
}
