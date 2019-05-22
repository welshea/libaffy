/**************************************************************************
 *
 * Filename: iron.c
 * 
 * Purpose: IRON command-line client, modified from MAS5 client.
 *
 * Creation date: 2011-03-16
 *
 * Author: Eric A. Welsh
 *
 *
 * Update History
 * --------------   
 * 04/04/05: Import from old libaffy (AMH)
 * 09/04/07: Add support for GCT output (AMH)
 * 09/11/07: Add check for empty filelist (AMH)
 * 03/03/08: Normalize filelist order (AMH)
 * 09/16/10: Add support for P/A calls (EAW)
 * 09/20/10: Pooled memory allocator (AMH)
 * 10/14/10: Pairwise normalization (EAW)
 * 08/05/11: Add MM subtraction flag (EAW)
 * 03/12/13: fixed --norm-none so as not to mean normalize probesets
 * 06/01/18: replaced compare_filename() with central compare_string() (EAW)
 * 05/22/19: added --ignore-chip-mismatch support (EAW)
 *
 **************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>

#include "affy.h"
#include "affy_mas5.h"
#include "affy_rma.h"
#include "argp.h"


char                 *output_file = "exprs-mas.txt";
AFFY_COMBINED_FLAGS   flags;
int                   gct_format = 0;
char                 *directory = ".";
char                **filelist;
int                  *mempool;

/* Administrative options */
const char *argp_program_version = affy_version;
const char *argp_program_bug_address = "<Eric.Welsh@moffitt.org>";
static struct argp_option options[] = {
  { "norm-mean",1,0,0,"Mean normalize probeset data" },
  { "norm-quantile",2,0,0,"Quantile normalize probe data" },
  { "bg-none",4,0,0,"Disable background correction" },
  { "directory",'d',"DIR",0,"Use directory as working directory" },
  { "bioconductor-compatability",5,0,0,"Calculate exprs identical to bioconductor"},
  { "dump-probes", 'p', "probe_file", OPTION_ARG_OPTIONAL,
    "Write raw probe values to a file" },
  { "gct-output-format",'g',0,0,"Write expressions in GCT format" },
  { "output-present-absent",'r',0,0,"Include present/absent calls in output" },
  { "output",'o',"FILE",0,"Write expressions to FILE" },
#if UNSUPPORTED
  { "scale",6,0,0,"Scale probesets to a constant factor" },
  { "no-scale",7,0,0,"Disable scaling probesets to a constant factor" },
#endif
  { "cdf",'c',"CDFDIR",0,"Use CDFDIR as location for CDF file" },
  { "norm-iron",8,"MODEL-FILE",OPTION_ARG_OPTIONAL,"IRON normalize data" },
  { "bg-mas5",9,0,0,"MAS5 background subtraction" },
  { "bg-rma",10,0,0,"RMA PM-only background subtraction" },
#if UNSUPPORTED
  { "bg-iron",11,0,0,"IRON background subtraction" },
  { "probe-tab",12,"file.probe_tab",OPTION_ARG_OPTIONAL,"Probe seqs for sequence-specific background" },
#endif
  { "use-mm-subtraction",13,0,0,"Subtract MM from PM signal"},
  { "no-mm-subtraction",14,0,0,"Do NOT subtract MM from PM signal"},
  { "tukey",15,0,0,"Tukey's Biweight probeset summarization"},
  { "median-polish",16,0,0,"Median Polish probeset summarization"},
  { "log2",17,0,0,"Output log2 probesets"},
  { "unlog",18,0,0,"Output non-logged probesets"},
  { "bg-rma-both",19,0,0,"RMA-like PM & MM background subtraction" },
  { "norm-none",20,0,0,"Disable normalization" },
  { "salvage",24,0,0,
    "Attempt to salvage corrupt CEL files (may still result in corrupt data!)" },
  { "ignore-chip-mismatch", 137,   0, 0, "Do not abort when multiple chips types are detected" },
  {0}
};

static error_t parse_opt(int key, char *arg, struct argp_state *state);
static struct argp argp =  { options, parse_opt, 0,
                             "MAS5 - MAS5 GeneChip Processing" };


int main(int argc, char **argv) 
{
  AFFY_CHIPSET *c;
  int           numfiles;
  AFFY_ERROR *err = NULL;

  mempool = h_malloc(sizeof(int));
  if (mempool == NULL)
  {
    fprintf(stderr, "malloc failed: out of memory\n");
    exit(EXIT_FAILURE);
  }

  err = affy_get_default_error();

  affy_rma_set_defaults(&flags);
  affy_mas5_set_defaults(&flags);

  argp_parse(&argp, argc, argv, 0, 0, 0);
                
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

  for (numfiles = 0; filelist[numfiles] != NULL; numfiles++)
    ;
  
 /* Do not sort.  APT processes in input order. Also more user desirable. */
/*  qsort(filelist, numfiles, sizeof(char *), compare_string); */

  print_flags(&flags, output_file);
  
  /* More efficient approach: do them one at a time */
  c = affy_mas5(filelist, &flags, err);

  if (gct_format)
  {
    affy_write_expressions_gct(c, output_file, err);
  }
  else
  {
    unsigned int writeopts = AFFY_WRITE_EXPR_DEFAULT;

    if (flags.output_present_absent)
      writeopts |= AFFY_WRITE_EXPR_PA;

    if (flags.output_log2)
      writeopts |= AFFY_WRITE_EXPR_LOG;

    affy_write_expressions(c, output_file, writeopts, err);  
  }

  print_corrupt_chips_to_stderr(c);

  h_free(mempool);
  affy_free_chipset(c);

  if (err)
    free(err);

  exit(EXIT_SUCCESS);
}

static error_t parse_opt(int key, char *arg, struct argp_state *state)
{
  switch (key) 
  {
    case 1:
      flags.use_normalization          = true;
      flags.use_mean_normalization     = true;
      flags.use_probeset_scaling       = true;
      flags.use_quantile_normalization = false;
      flags.use_pairwise_normalization = false;
      break;
    case 2:
      flags.use_normalization          = true;
      flags.use_quantile_normalization = true;
      flags.use_mean_normalization     = false;
      flags.use_pairwise_normalization = false;
      flags.use_probeset_scaling       = false;
      break;
    case 3:
      flags.use_background_correction = true;
      break;
    case 4:
      flags.use_background_correction = false;
      flags.bg_mas5 = false;
      flags.bg_rma  = false;
      flags.bg_iron = false;
      break;
    case 5:
      flags.bioconductor_compatability = true;
      break;
    case 6:
      flags.use_probeset_scaling = true;
      break;
    case 7:
      flags.use_probeset_scaling = false;
      break;
    case 8:
      flags.use_normalization          = true;
      flags.use_pairwise_normalization = true;
      if (arg != NULL)
      {
        flags.pairwise_model_filename = h_strdup(arg);
        hattach(flags.pairwise_model_filename, mempool);
      }
      flags.use_quantile_normalization = false;
      flags.use_mean_normalization     = false;
      flags.use_probeset_scaling       = false;
      break;
    case 9:
      flags.use_background_correction = true;
      flags.bg_mas5 = true;
      flags.bg_rma = false;
      flags.bg_rma_both = false;
      flags.bg_iron = false;
      flags.use_mm_probe_subtraction = true;
      break;
    case 10:
      flags.use_background_correction = true;
      flags.bg_mas5 = false;
      flags.bg_rma = true;
      flags.bg_rma_both = false;
      flags.bg_iron = false;
      flags.use_mm_probe_subtraction = false;
      break;
    case 11:
      flags.use_background_correction = true;
      flags.bg_mas5 = false;
      flags.bg_rma = false;
      flags.bg_rma_both = false;
      flags.bg_iron = true;
      flags.use_mm_probe_subtraction = false;
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
      flags.use_tukey_biweight = true;
      flags.use_median_polish = false;
      break;
    case 16:
      flags.use_tukey_biweight = false;
      flags.use_median_polish = true;
      break;
    case 17:
      flags.output_log2 = true;
      break;
    case 18:
      flags.output_log2 = false;
      break;
    case 19:
      flags.use_background_correction = true;
      flags.bg_mas5 = false;
      flags.bg_rma = false;
      flags.bg_rma_both = true;
      flags.bg_iron = false;
      flags.use_mm_probe_subtraction = false;
      break;
    case 20:
      flags.use_normalization          = false;
      flags.use_quantile_normalization = false;
      flags.use_pairwise_normalization = false;
      flags.use_mean_normalization     = false;
      flags.use_probeset_scaling       = false;
      break;
    case 24:
      flags.salvage_corrupt = true;
      break;
    case 137:
      flags.ignore_chip_mismatch = true;
      break;

    case 'g':
      gct_format = true;
      break;
    case 'r':
      flags.output_present_absent = true;
      break;
    case 'd':
      directory = h_strdup(arg);
      hattach(directory, mempool);
      break;
    case 'o':
      output_file = h_strdup(arg);
      hattach(output_file, mempool);
      break;
    case 'p':
      flags.dump_probe_values = true;
      if (arg != NULL) 
      {
        flags.probe_filename = h_strdup(arg);
        hattach(flags.probe_filename, mempool);
      }
      break;
    case 'c':
      flags.cdf_directory = h_strdup(arg);
      hattach(flags.cdf_directory, mempool);
      break;
    case ARGP_KEY_ARG: /* Do not try to get filenames one at a time */
      return (ARGP_ERR_UNKNOWN);
    case ARGP_KEY_ARGS:
      filelist = state->argv+state->next;
      break;
    default:
      return (ARGP_ERR_UNKNOWN);
  }

  return (0);
}
