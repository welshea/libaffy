
/**************************************************************************
 *
 * Filename: rma.c
 * 
 * Purpose:  RMA command-line client.
 *
 * Creation: 4 April, 2005
 *
 * Author:   Steven Eschrich
 *
 *
 * Update History
 * --------------
 * 04/04/05: Import from old libaffy (AMH)
 * 04/08/05: DAT pixeldata now a 2D array of unsigned ints (AMH) 
 * 04/19/05: Add AFFY_CELL and AFFY_PIXREGION types (AMH)
 * 11/09/06: Add flag to ignore AFFX probes, freshen style a bit (AMH)
 * 04/06/07: Add support for incremental RMA (AMH)
 * 09/04/07: Add support for GCT output (AMH)
 * 09/11/07: Add check for empty filelist (AMH)
 * 03/03/08: Normalize filelist order (AMH)
 * 06/05/08: New error handling scheme (AMH)
 * 10/18/10: Pairwise normalization (AMH)
 * 10/19/10: Add unlog option (AMH)
 * 10/22/10: Use new AFFY_COMBINED_FLAGS instead of AFFY_RMA_FLAGS (EAW)
 * 08/05/11: Add MM subtraction flag (EAW)
 * 06/01/18: replaced compare_filename() with central compare_string() (EAW)
 * 05/22/19: added --ignore-chip-mismatch support (EAW)
 *
 **************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "affy.h"
#include "affy_mas5.h"
#include "affy_rma.h"
#include "argp.h"

char     *output_file   = "exprs-rma.txt";
char     *affinity_file = "affinities-rma.txt";
char     *means_file    = "means-rma.txt";
char     *probe_file    = "probes-rma.txt";
char     *directory     = ".";
bool      gct_format    = false;
char    **filelist;
int      *mempool;

AFFY_COMBINED_FLAGS flags;

int debug_level = 2;

/* Administrative options */
const char *argp_program_version=affy_version;
const char *argp_program_bug_address="<Eric.Welsh@moffitt.org>";
static struct argp_option options[] = {
  { "norm-quantile",1,0,0,"Quantile normalize probe data" },
  { "norm-none",2,0,0,"Disable normalization" },
#if UNSUPPORTED
  { "ignore-affx-probes",3,0,0, 
    "Ignore AFFX*/control probes when normalizing" },
#endif
  { "dump-affinities",'W',"dump-file",OPTION_ARG_OPTIONAL, 
    "Write affinity values to a dump file" },
  { "read-affinities",'A',"affinity_file",0,
    "Use saved affinities (incremental RMA)" },
  { "read-means",'M',"mean_file",0,
    "Use saved means (incremental RMA)" },
  { "dump-means",'w',"mean_file",OPTION_ARG_OPTIONAL,
    "Write mean values to a savefile" },
  { "dump-probes", 'p', "probe_file", OPTION_ARG_OPTIONAL,
    "Write raw probe values to a file" },
  { "bg-none",5,0,0,"Disable background correction" },
  { "gct-output-format",'g',0,0,"Output expressions in gct format" },
  { "norm-mean",'m',"TARGET",OPTION_ARG_OPTIONAL,
    "Normalize expression on chip to TARGET" },
  { "dir",'d',"DIRECTORY",0,"Use DIRECTORY as working directory" },
  { "cdf",'c',"CDFDIR",0,"Use CDFDIR as location for CDF file" },
  { "output",'o',"OUTPUTFILE",0,"Output expressions to OUTPUTFILE" },
#if UNSUPPORTED
  { "norm-iron",6,"MODEL-FILE",OPTION_ARG_OPTIONAL,"IRON normalize data" },
#endif
  { "unlog",7,0,0,"Output expressions in normal rather than log scale" },
#if UNSUPPORTED
  { "use-rma-probeset-singletons",8,0,0,"Summarize RMA probsets individually" },
  { "bg-mas5",9,0,0,"MAS5 Background-correct expression (without MM subtraction)" },
#endif
  { "bg-rma",10,0,0,"RMA Background-correct expression" },
#if UNSUPPORTED
  { "bg-iron",11,0,0,"IRON Background-correct expression" },
#endif
  { "probe-tab",12,"file.probe_tab",OPTION_ARG_OPTIONAL,"Probe seqs for sequence-specific background" },
#if UNSUPPORTED
  { "use-mm-subtraction",13,0,0,"Subtract MM from PM signal"},
  { "no-mm-subtraction",14,0,0,"Do NOT subtract MM from PM signal"},
#endif
  { "bioconductor-compatability",15,0,0,"Calculate exprs identical to bioconductor"},
#if UNSUPPORTED
  { "normalize-before-bg",16,0,0,"Normalize before background subtraction (only with bg-rma)" },
#endif
  { "log2",17,0,0,"Output log2 probesets"},
  { "salvage",24,0,0,
    "Attempt to salvage corrupt CEL files (may still result in corrupt data!)" },
  { "ignore-chip-mismatch", 137,   0, 0, "Do not abort when multiple chips types are detected" },
  {0}
};

static error_t parse_opt(int key, char *arg, struct argp_state *state);
static struct argp argp = { options, parse_opt, 0, 
                            "rma - RMA GeneChip Processing"};


int main(int argc, char **argv)
{
  AFFY_CHIPSET *c;
  int           numfiles;
  unsigned int  write_opts = AFFY_WRITE_EXPR_DEFAULT;
  AFFY_ERROR   *err = NULL;

  mempool = h_malloc(sizeof(int));
  if (mempool == NULL)
  {
    fprintf(stderr, "malloc failed: out of memory\n");
    exit(EXIT_FAILURE);
  }

  err = affy_get_default_error();
        
  affy_mas5_set_defaults(&flags);
  affy_rma_set_defaults(&flags);
  argp_parse(&argp, argc, argv, 0, 0, 0);

  /* Check for mutually exclusive options. */
  if (   (flags.use_saved_means && flags.dump_expression_means) 
      || (flags.use_saved_affinities && flags.dump_probe_affinities))
  {
    fprintf(stderr, "Error: probe affinities or mean values cannot be "
            "simultaneously saved and dumped\n");

    if (err)
      free(err);

    exit(EXIT_FAILURE);
  }

  /* If files is NULL, open all CEL files in the current directory */
  if (filelist == NULL) 
    filelist = affy_list_files(directory, ".cel", err);

  /* Give up if we have no files to operate on */
  if ((filelist == NULL) || (filelist[0] == NULL))
  {
    fprintf(stderr, 
            "no CEL files specified or found in current dir, exiting\n");

    if (err)
      free(err);

    exit(EXIT_FAILURE);
  }

  /* Normalize argument order. */
  for (numfiles = 0; filelist[numfiles] != NULL; numfiles++)
    ;

  /* Do not sort.  APT processes in input order. Also more user desirable. */
/*  qsort(filelist, numfiles, sizeof(char *), compare_string); */

  print_flags(&flags, output_file);
  
  c = affy_rma(filelist, &flags, err);

  if (flags.output_log2 == false)
    write_opts |= AFFY_WRITE_EXPR_UNLOG;
  
  /* Finally write out the entire expression array */
  if (gct_format)
    affy_write_expressions_gct(c, output_file, err);
  else
    affy_write_expressions(c, output_file, write_opts, err);

  print_corrupt_chips_to_stderr(c);

  affy_free_chipset(c);
  h_free(mempool);

  if (err)
    free(err);

  exit(EXIT_SUCCESS);
}

static error_t parse_opt(int key, char *arg, struct argp_state *state)
{
  switch (key)
  {
    case 1:
      flags.use_normalization = true;
      flags.use_quantile_normalization = true;
      flags.use_mean_normalization     = false;
      flags.use_probeset_scaling       = false;
      flags.use_pairwise_normalization = false;
      break;
    case 2:
      flags.use_normalization = false;
      flags.use_quantile_normalization = false;
      flags.use_mean_normalization     = false;
      flags.use_probeset_scaling       = false;
      flags.use_pairwise_normalization = false;
      break;
#if UNSUPPORTED
    case 3:
      flags.normalize_affx_probes = false;
      break;
#endif
    case 4:
      flags.use_background_correction = true;
      break;
    case 5:
      flags.use_background_correction = false;
      flags.bg_mas5 = false;
      flags.bg_rma  = false;
      flags.bg_iron = false;
      break;
    case 6:
      flags.use_pairwise_normalization = true;
      flags.reuse_affinities = true;
      if (arg != NULL)
      {
        flags.pairwise_model_filename = h_strdup(arg);
        hattach(flags.pairwise_model_filename, mempool);
      }
      break;
    case 7:
      flags.output_log2 = false;
      break;
    case 8:
      flags.use_saved_affinities = false;
      flags.use_rma_probeset_singletons = true;
      break;
    case 9:
      flags.use_background_correction = true;
      flags.bg_mas5 = true;
      flags.bg_rma = false;
      flags.bg_iron = false;
      break;
    case 10:
      flags.use_background_correction = true;
      flags.bg_mas5 = false;
      flags.bg_rma = true;
      flags.bg_iron = false;
      break;
    case 11:
      flags.use_background_correction = true;
      flags.bg_mas5 = false;
      flags.bg_rma = false;
      flags.bg_iron = true;
      break;
    case 12:
      if (arg != NULL) 
      {
        flags.probe_tab_filename = h_strdup(arg);
        hattach(flags.probe_tab_filename, mempool);
      }
      break;
    case 13:
      flags.use_mm_probe_subtraction = true;
      flags.use_mm_probeset_subtraction = false;
      break;
    case 14:
      flags.use_mm_probe_subtraction = false;
      break;
    case 15:
      flags.bioconductor_compatability = true;
      break;
    case 16:
      flags.normalize_before_bg = true;
      break;
    case 17:
      flags.output_log2 = true;
      break;
    case 24:
      flags.salvage_corrupt = true;
      break;
    case 137:
      flags.ignore_chip_mismatch = true;
      break;

    case 'd':
      directory = h_strdup(arg);
      hattach(directory, mempool);
      break;
    case 'm':
      flags.use_mean_normalization = true;
      flags.use_probeset_scaling   = true;
      if (arg != NULL) 
        flags.mean_normalization_target_mean = atof(arg);
      break;
    case 'W':
      flags.dump_probe_affinities=true;
      if (arg != NULL) 
      {
        flags.affinities_filename = h_strdup(arg);
        hattach(flags.affinities_filename, mempool);
      }
      break;
    case 'w':
      flags.dump_expression_means = true;
      if (arg != NULL) 
      {
        flags.means_filename = h_strdup(arg);
        hattach(flags.means_filename, mempool);
      }
      break;
    case 'p':
      flags.dump_probe_values = true;
      if (arg != NULL) 
      {
        flags.probe_filename = h_strdup(arg);
        hattach(flags.probe_filename, mempool);
      }
      break;
    case 'A':
      flags.use_saved_affinities = true;
      flags.affinities_filename = h_strdup(arg);
      hattach(flags.affinities_filename, mempool);
      break;
    case 'M':
      flags.use_saved_means = true;
      flags.means_filename = h_strdup(arg);
      hattach(flags.means_filename, mempool);
      break;
    case 'g':
      gct_format = true;
      break;
    case 'c':
      flags.cdf_directory = h_strdup(arg);
      hattach(flags.cdf_directory, mempool);
      break;
    case 'o':
      output_file = h_strdup(arg);
      hattach(output_file, mempool);
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
