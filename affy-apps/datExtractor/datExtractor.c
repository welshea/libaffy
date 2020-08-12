
/* 2020/08/12:  modified to (hopefully) compile with f->cdf_filename stuff
 */

#include "affy.h"
#include "argp.h"

int write_probeset(AFFY_CHIP *cp, AFFY_PROBESET *ps, AFFY_ERROR *err);

/* Administrative options */
const char *argp_program_version=affy_version;
const char *argp_program_bug_address="<Eric.Welsh@moffitt.org>";
static struct argp_option options[]= {
  {"file",'f',"FILE",0,"Use FILE as the dat/cel file."},
  {"type",'t',"TYPE",0,"What file type as output (tiff|text)."},
  {0}
};

static error_t parse_opt(int key, char *arg, struct argp_state *state);
static struct argp argp={options, parse_opt,"DATFIILE ...","datExtractor - Affymetrix DAT file processing"};


#define OT_TEXT 1
#define OT_TIFF 2

char  *ext = "txt";
int    OUTPUT_TYPE   = OT_TEXT;
char **probeset_list = NULL;
char  *file          = NULL;

int main(int argc, char **argv)
{
  int         i, j;
  AFFY_CHIP  *cp;
  AFFY_ERROR *err = NULL;
  AFFY_COMBINED_FLAGS flags;

#ifndef STORE_CEL_QC
  fprintf(stderr, "Error: You must recompile with -DSTORE_CEL_QC for datExtractor to work.\n");
  fprintf(stderr, "       -DSTORE_CEL_QC is not defined by default, due to the increased memory\n");
  fprintf(stderr, "       overhead that it would incur for most of the other software.\n");
  fprintf(stderr, "       \"scons cflags_bundle=datextractor\" adds -DSTORE_CEL_QC to cflags_optimize\n");
  exit(AFFY_ERROR_NOTSUPP);
#endif
  
  err = affy_get_default_error();

  affy_rma_set_defaults(&flags);
  affy_mas5_set_defaults(&flags);

  argp_parse(&argp, argc, argv, 0,0,0);

  /* Try to load data. Errors should result in dying. */
  cp = calloc(1, sizeof(AFFY_CHIP));
  cp->dat = affy_load_dat_file(file, err);
  cp->cdf = affy_load_cdf_file(cp->dat->probe_array_type,
                               NULL, flags, err);

  /* Identify a probeset */
  for (i = 0; probeset_list[i] != NULL; i++)
  {
    /* Try and find the probeset in the chip */
    for (j = 0; j < cp->cdf->numprobesets; j++)
    {
      if (!strcmp(probeset_list[i], cp->cdf->probeset[j].name))
        write_probeset(cp, &(cp->cdf->probeset[j]), err);
    }
  }

  if (err)
    free(err);

  return (0);
}

static error_t parse_opt(int key, char *arg, struct argp_state *state)
{
   switch (key) {
    case 't':
      if (!strcmp(arg, "tif") || !strcmp(arg, "tiff"))
      {
        OUTPUT_TYPE = OT_TIFF;
        ext = "tiff";
      }
      else
      {
        OUTPUT_TYPE = OT_TEXT;
        ext = "txt";
      }
      break;
    case 'f':
      file = strdup(arg);
      break;
    case ARGP_KEY_ARG: /* Do not try to get probesets one at a time */
      return ARGP_ERR_UNKNOWN;
    case ARGP_KEY_ARGS:
      probeset_list=state->argv+state->next;
      break;
    case ARGP_KEY_NO_ARGS:
      fprintf(stderr,"You must supply probeset information\n");
      argp_usage(state);
      return ARGP_ERR_UNKNOWN; 
    default:
      return ARGP_ERR_UNKNOWN;
   }
   return 0;
}

/*
 * Given a probeset, write out data for each probe,
 * for mm and pm.
 */
int write_probeset(AFFY_CHIP *cp, AFFY_PROBESET *ps, AFFY_ERROR *err)
{
  int i;
  AFFY_PIXREGION *pr_pm, *pr_mm;
  char buf_pm[1024], buf_mm[1024];

#ifdef STORE_CEL_QC
  for (i = 0; i < ps->numprobes; i++)
  {
    AFFY_PROBE *p = &(ps->probe[i]);

    /* PM */
    pr_pm = affy_pixels_from_cell(cp, p->pm.x, p->pm.y, err);
    sprintf(buf_pm, "%s-pm-%d.%s", ps->name, i, ext);
    pr_mm = affy_pixels_from_cell(cp, p->mm.x, p->mm.y, err);
    sprintf(buf_mm, "%s-mm-%d.%s", ps->name, i, ext);

    switch (OUTPUT_TYPE)
    {
    case OT_TEXT:
      affy_pixregion2text(pr_pm, buf_pm, err);
      affy_pixregion2text(pr_mm, buf_mm, err);
      break;
    case OT_TIFF:
      affy_pixregion2tiff(pr_pm, buf_pm, err);
      affy_pixregion2tiff(pr_mm, buf_mm, err);
      break;
    }
  }
#endif

  return (0);
}

