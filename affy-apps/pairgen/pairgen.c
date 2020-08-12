
/**************************************************************************
 *
 * Filename: pairgen.c
 * 
 * Purpose:  Utility to create a 'model' chip for use with pairwise
 *           normalization..
 *
 * Creation: 7 October, 2010
 *
 * Author:   Andrew M. Hoerter
 *
 *
 * Update History
 * --------------
 * 10/07/10: File creation (AMH)
 * 11/19/10: Added average and median flags (EAW)
 * 03/11/11: Changed mean to geometric mean (EAW)
 * 03/11/11: Abort on chip load failure (EAW)
 * 05/22/19: added --ignore-chip-mismatch support (EAW)
 * 08/12/20: disable searching current working directory for CEL files (EAW)
 * 08/12/20: pass flags to affy_create_chipset() (EAW)
 *
 **************************************************************************/

#include <stdio.h>

#include "affy.h"
#include "argp.h"

#define SEARCH_WORKING_DIR 0

AFFY_COMBINED_FLAGS   flags;
char                 *output_file   = "median.CEL";
char                 *directory     = ".";
char                 *cdf_directory = ".";
char                **filelist;
int                  *mempool = NULL;
double              **all_probes;
int                   debug_level = 2;
int                   opt_average_flag = 0;

/* Administrative options */
const char *argp_program_version     = affy_version;
const char *argp_program_bug_address = "<Eric.Welsh@moffitt.org>";
static struct argp_option options[]  = 
{
  { "dir",      'd', "DIRECTORY",  0, "Use DIRECTORY as working directory"  },
  { "cdf",      'c', "CDFDIR",     0, "Use CDFDIR as location for CDF file" },
  { "output",   'o', "OUTPUTFILE", 0, "Output model chip to OUTPUTFILE"     },
  { "average",  'a',           0,  0, "Use geometric mean of probes"        },
  { "median",   'm',           0,  0, "Use median of probes"                },
  { "salvage",  24,            0,  0, "Attempt to salvage corrupt CEL files (may still result in corrupt data!)" },
  { NULL }
};

static error_t parse_opt(int key, char *arg, struct argp_state *state);
static struct argp argp = { options, 
                            parse_opt, 
                            0, 
			    "pairgen - Pairwise normalization model chip generator"};


int main(int argc, char **argv)
{
  AFFY_CHIPSET     *cs = NULL;
  AFFY_CHIPSET     *temp = NULL;
  AFFY_CHIP         model_chip;
  AFFY_CELFILE      model_cel;
  FILE             *fp;
  int               status = EXIT_FAILURE, i;
  size_t            sz;
  unsigned int      max_chips;
  affy_uint16       row, col;
  affy_uint32       num_all_probes;
  affy_int32        probe_idx;
  char             *chip_type, **p;
  LIBUTILS_PB_STATE pbs;
  AFFY_ERROR       *err = NULL;

  pb_init(&pbs);

  mempool = h_malloc(sizeof(int));
  if (mempool == NULL)
  {
    fprintf(stderr, "malloc failed: out of memory\n");
    exit(EXIT_FAILURE);
  }

  err = affy_get_default_error();
        
  flags.bioconductor_compatability = false;
  argp_parse(&argp, argc, argv, 0, 0, 0);

  /* If files is NULL, open all CEL files in the current working directory */
  if (filelist == NULL && SEARCH_WORKING_DIR)
    filelist = affy_list_files(directory, ".cel", err);

  /* Give up if we have no files to operate on */
  if ((filelist == NULL) || (filelist[0] == NULL))
  {
    if (SEARCH_WORKING_DIR)
    {
      fprintf(stderr, 
              "no CEL files specified or found in current working directory, exiting\n");
    }
    else
    {
      fprintf(stderr, 
              "no CEL files specified, exiting\n");
    }

    goto cleanup;
  }

  /* Count up chips */
  for (p = filelist, max_chips = 0; *p != NULL; p++)
    max_chips++;
  
  /* Get the array type. */
  chip_type = affy_get_cdf_name_from_cel(filelist[0], err);
  AFFY_CHECK_ERROR_GOTO(err, cleanup);

  hattach(chip_type, mempool);

  /* Create temp chipset */
  cs = affy_create_chipset(1, chip_type, cdf_directory, &flags, err);
  AFFY_CHECK_ERROR_GOTO(err, cleanup);

  temp = affy_clone_chipset(cs, err);
  AFFY_CHECK_ERROR_GOTO(err, cleanup);

  cs = affy_resize_chipset(cs, max_chips, err);
  AFFY_CHECK_ERROR_GOTO(err, cleanup);

  hattach(cs, mempool);

  model_cel.filename    = output_file;
  
  model_cel.numrows     = cs->cdf->numrows;
  model_cel.numcols     = cs->cdf->numcols;
  model_cel.nummasks    = 0;
  model_cel.numoutliers = 0;
  model_cel.outlier     = NULL;
  model_cel.mask        = NULL;

  num_all_probes = cs->cdf->numrows * cs->cdf->numcols;

  model_cel.data = h_suballoc(mempool, 
                              cs->cdf->numrows * sizeof(AFFY_CELL *));
  if (model_cel.data == NULL)
  {
    fprintf(stderr, "malloc failed: out of memory");
    goto cleanup;
  }
  
  model_cel.data[0] = h_suballoc(model_cel.data,
                                 num_all_probes * sizeof(AFFY_CELL));
  if (model_cel.data[0] == NULL)
  {
    fprintf(stderr, "malloc failed: out of memory");
    goto cleanup;
  }

  for (row = 1; row < cs->cdf->numrows; row++)
    model_cel.data[row] = model_cel.data[row - 1] + cs->cdf->numcols;
  
  /* Set up the 2D array of chips/probes */
  all_probes = h_suballoc(mempool, num_all_probes * sizeof(double *));
  if (all_probes == NULL)
  {
    fprintf(stderr, "malloc failed: out of memory\n");
    goto cleanup;
  }

  sz = max_chips * num_all_probes * sizeof(double);
  all_probes[0] = h_suballoc(all_probes, sz);
  if (all_probes[0] == NULL)
  {
    fprintf(stderr, "malloc failed: out of memory\n");
    goto cleanup;
  }

  for (probe_idx = 1; probe_idx < num_all_probes; probe_idx++)
    all_probes[probe_idx] = all_probes[probe_idx - 1] + max_chips;

  for (i = 0; i < max_chips; i++)
  {
    /* Load each chip, abort on failure */
    affy_load_chipset_single(cs, filelist[i],
                             flags.ignore_chip_mismatch, err);
    if (err->type != AFFY_ERROR_NONE)
      AFFY_CHECK_ERROR_GOTO(err, cleanup);

    /* Temp chipset now contains the most recently loaded chip */
    temp->chip[0] = cs->chip[(cs->num_chips) - 1];
    temp->num_chips = 1;

    for (row = 0; row < cs->cdf->numrows; row++)
    {
      for (col = 0; col < cs->cdf->numcols; col++)
      {
        probe_idx = (cs->cdf->numcols * row) + col;
        
        all_probes[probe_idx][i] = temp->chip[0]->cel->data[row][col].value;
      }
    }
    
    affy_mostly_free_chip(temp->chip[0]);
  }

  info("Finished reading %u CEL files", max_chips);
  
  if (opt_average_flag)
    pb_begin(&pbs, num_all_probes, "Calculating averages for all probes");
  else
    pb_begin(&pbs, num_all_probes, "Calculating medians for all probes");
  
  for (row = 0; row < cs->cdf->numrows; row++)
  {
    for (col = 0; col < cs->cdf->numcols; col++)
    {
      probe_idx = (cs->cdf->numcols * row) + col;
    
      if (opt_average_flag)
        model_cel.data[row][col].value =
               affy_mean_geometric_floor_1(all_probes[probe_idx], max_chips);
      else
        model_cel.data[row][col].value = affy_median(all_probes[probe_idx],
                                                   max_chips, &flags);
#ifdef STORE_CEL_QC
      model_cel.data[row][col].numpixels = 0;
      model_cel.data[row][col].stddev    = 0.0;
#endif

      pb_tick(&pbs, 1, "");
    }
  }

  pb_finish(&pbs, "done");

  fp = fopen(output_file, "wb");
  if (fp == NULL)
    AFFY_HANDLE_ERROR_GOTO("couldn't fopen output file",
                           AFFY_ERROR_IO,
                           err,
                           cleanup);
  model_chip.cdf = cs->cdf;
  model_chip.cel = &model_cel;
  affy_write_binary_cel_file(fp, &model_chip, err);
  fclose(fp);
  AFFY_CHECK_ERROR_GOTO(err, cleanup);

  status = EXIT_SUCCESS;

  print_corrupt_chips_to_stderr(cs);

cleanup:  
  affy_free_chipset(cs);
  h_free(mempool);
  pb_cleanup(&pbs);

  if (err)
    free(err);

  exit(status);
}

static error_t parse_opt(int key, char *arg, struct argp_state *state)
{
  switch (key)
  {
    case 24:
      flags.salvage_corrupt = true;
      break;
    case 'd':
      directory = h_strdup(arg);
      hattach(directory, mempool);
      break;
    case 'c':
      cdf_directory = h_strdup(arg);
      hattach(cdf_directory, mempool);
      break;
    case 'o':
      output_file = h_strdup(arg);
      hattach(output_file, mempool);
      break;
    case 'a':
      opt_average_flag = 1;
      break;
    case 'm':
      opt_average_flag = 0;
      break;
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
