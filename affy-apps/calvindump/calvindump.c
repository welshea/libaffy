
/**************************************************************************
 *
 * Filename: calvindump.c
 * 
 * Purpose:  Dump all metadata from an Affy Calvin file.
 *
 * Creation: 1 December, 2008
 *
 * Author:   Andrew M Hoerter
 *
 *
 * Update History
 * --------------
 * 12/01/08: Creation (AMH)
 *
 **************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <signal.h>

#include "affy.h"
#include "argp.h"

char **filelist = NULL;

/* Administrative options */
const char *argp_program_version="calvindump $Revision$";
const char *argp_program_bug_address="<Eric.Welsh@moffitt.org>";

static struct argp_option options[] = {
  {0}
};

static error_t parse_opt(int key, char *arg, struct argp_state *state);
static struct argp argp = { options, parse_opt, 0, 
			    "calvindump - Affymetrix Calvin debugging tool"};

void crash(AFFY_ERROR *err)
{
  raise(11);
}

int main(int argc, char **argv)
{
  int         i;
  AFFY_ERROR *err = NULL;

  err = affy_get_default_error();

  /*err->handler = crash;*/
        
  argp_parse(&argp, argc, argv, 0, 0, 0);

  /* Give up if we have no files to operate on */
  if (filelist == NULL)
  {
    fprintf(stderr, "no Calvin files specified, exiting\n");

    if (err)
      free(err);

    exit(EXIT_FAILURE);
  }
  
  for (i = 0; filelist[i] != NULL; i++)
  {
    AFFY_CALVINIO         *cio;
    AFFY_CALVIN_CONTAINER *cc;
    FILE                  *fp;

    if ((fp = fopen(filelist[i], "rb")) == NULL)
    {
      fprintf(stderr, "couldn't open %s for reading\n", filelist[i]);
      continue;
    }

    cio = affy_calvinio_init(fp, err);
    cc  = affy_calvin_load_container(cio, err);

    printf("\n----------\n%s\n----------\n", filelist[i]);
    affy_dump_calvin_container(cc);

    affy_free_calvin_container(cc);
    free(cc);

    affy_calvinio_free(cio);
  }

  if (err)
    free(err);
  
  exit(EXIT_SUCCESS);
}

static error_t parse_opt(int key, char *arg, struct argp_state *state)
{
  switch (key)
  {
    case ARGP_KEY_ARG: /* Do not try to get filenames one at a time, use all */
      return ARGP_ERR_UNKNOWN;
    case ARGP_KEY_ARGS:
      filelist = state->argv + state->next;
      break;
    default:
      return ARGP_ERR_UNKNOWN;
  }

  return (0);
}
