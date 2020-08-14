
/**************************************************************************
 *
 * Filename: iron_generic.c
 * 
 * Purpose:  Generic IRON command-line client.
 *
 * Creation: 6 April, 2011
 *
 * Author:   Eric A. Welsh
 *
 *
 * Update History
 * --------------
 * 04/06/11: Creation.  Modified from rma.c (EAW)
 * 01/10/13: renamed --iron-linear flag to --iron-global-scaling
 * 04/08/13: removed attempt to load CEL files if no files given
 * 09/18/13: added iron_fit_both_x_y flag: better norm, may distrub rank order
 * 09/18/13: added iron_fit_window_frac (EAW)
 * 10/04/17: added iron_condense_training (EAW)
 * 06/01/18: replaced compare_filename() with compare_string() (EAW)
 * 06/01/18: added support for probeset exclusions during IRON training (EAW)
 * 09/11/18: extended IRON exclusion flags to support *not* scaling them (EAW)
 * 09/14/18: split off unscaled probesets into separate spikeins handling (EAW)
 * 10/11/19: --rnaseq now enables --iron-condense-training (EAW)
 *           this generally results in better fits, due to minimizing the
 *           contribution of large number of low read counts
 * 10/15/19: added --floor-non-zero-to-one flag
 * 05/16/20: change median polish to Tukey's Bi-weight probeset summarization
 *           change does not actually affect results, since 1:1 probeset:probe
 * 08/12/20: change description of --bioconductor-compatability (EAW)
 * 08/12/20: changed error message for no input files specified (EAW)
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
bool      unlog_expr    = false;
char    **filelist;
int      *mempool;

AFFY_COMBINED_FLAGS flags;

int debug_level = 2;

/* Administrative options */
const char *argp_program_version=affy_version;
const char *argp_program_bug_address="<Eric.Welsh@moffitt.org>";
static struct argp_option options[] = {
  { "norm-quantile",1,0,0,"Quantile normalize data" },
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
  { "norm-iron",6,"MODEL-FILE",OPTION_ARG_OPTIONAL,"Pairwise normalize data" },
  { "unlog",7,0,0,"Output expressions in normal rather than log scale" },
#if UNSUPPORTED
  { "bg-mas5",9,0,0,"MAS5 Background-correct expression (without MM subtraction)" },
#endif
  { "bg-rma",10,0,0,"RMA Background-correct expression" },
#if UNSUPPORTED
  { "bg-iron",11,0,0,"IRON Background-correct expression" },
#endif
  { "probe-tab",12,"file.probe_tab",OPTION_ARG_OPTIONAL,"Probe seqs for sequence-specific background" },
  { "bioconductor-compatability",13,0,0,"Calculate exprs more similar to bioconductor"},
#if UNSUPPORTED
  { "normalize-before-bg",14,0,0,"Normalize before background subtraction (only with bg-rma)" },
#endif
  { "log2",17,0,0,"Output log2 probesets (default)"},
  { "iron-global-scaling",21,0,0,"IRON: use per-chip global scaling"},
  { "iron-non-linear",22,0,0,"IRON: use per-chip non-linear scaling (default)"},
  { "iron-weight-exponent",23,"EXPONENT",0,
     "IRON: Pseudo-density weight exponent (microarray default: 4, unweighted: 0)" },
#if 1
  /* will include this later after testing it more */
  { "iron-untilt",25,0,0,"IRON: use per-chip linear-fit untilt scaling"},
#endif
  { "iron-fit-both-x-y",26,0,0,"IRON: Fit to both X and Y (better normalization, but may alter rank orders)"},
  { "iron-fit-only-y",27,0,0,"IRON: Fit only to Y (default)" },
  { "iron-fit-window-frac",28,"FRACTION",0,
     "IRON: Fit window width fraction (default: 0.10)" },
  { "proteomics",29,0,0,"Use defaults suitable for proteomics" },
  { "rnaseq",30,0,0,"Use defaults suitable for RNASeq" },
  { "floor-to-min",131,0,0,"Set final zero/near-zero values to min value per sample" },
  { "floor-none",132,0,0,"Do not apply any floors to final signals (default)" },
  { "iron-condense-training",133,0,0,"Condense identical X,Y prior to training" },
  { "iron-no-condense-training",134,0,0,"Do not condense identical X,Y prior to training (default)" },
  { "iron-exclusions",'x',"EXCLUSIONSFILE",0,"Ignore probesets from EXCLUSIONSFILE during curve fitting" },
  { "iron-spikeins",'S',"SPIKEINSFILE",0,"File with probesets, usually spikeins, to be left unnormalized" },
  { "bg-global",135,0,0,"Global background subtraction" },
  { "floor-non-zero-to-one",136,0,0,"Floor final non-zero values to 1.0" },
  { "microarray",137,0,0,"Use defaults suitable for microarrays (default)" },
  {0}
};

static error_t parse_opt(int key, char *arg, struct argp_state *state);
static struct argp argp = { options, parse_opt, 0, 
                            "iron_generic - IRON generic spreadsheet (un-logged intensities) processing"};


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

  /* IRON specific defaults */
  flags.use_pairwise_normalization = true;
  flags.use_mean_normalization     = false;
  flags.use_probeset_scaling       = false;
  flags.use_quantile_normalization = false;
  flags.bg_rma                     = true;
  flags.bg_mas5                    = false;
  flags.bg_iron                    = false;
  flags.bg_global                  = false;
  flags.use_mm_probe_subtraction   = false;
  flags.use_tukey_biweight         = true;
  flags.use_median_polish          = false;
  flags.output_log2                = true;

  argp_parse(&argp, argc, argv, 0, 0, 0);

  /* Check for mutually exclusive options. */
  if (   (flags.use_saved_means && flags.dump_expression_means) 
      || (flags.use_saved_affinities && flags.dump_probe_affinities))
  {
    fprintf(stderr, "Error: probe affinities or mean values cannot be "
            "simultaneously saved and dumped\n");

    h_free(mempool);

    if (err)
      free(err);

    exit(EXIT_FAILURE);
  }

  /* Give up if we have no files to operate on */
  if ((filelist == NULL) || (filelist[0] == NULL))
  {
    fprintf(stderr, "no input files specified, exiting\n");

    h_free(mempool);

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
  
  c = affy_illumina(filelist, &flags, err);

  if (flags.floor_non_zero_to_one)
    affy_floor_probeset_non_zero_to_one(c, err);
  
  if (flags.floor_to_min_non_zero)
    affy_floor_probeset_to_min_non_zero(c, err);


  if (flags.output_log2)
    write_opts |= AFFY_WRITE_EXPR_LOG;
  
  /* Finally write out the entire expression array */
  if (gct_format)
    affy_write_expressions_gct(c, output_file, err);
  else
    affy_write_expressions(c, output_file, write_opts, err);

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
      flags.iron_global_scaling_normalization = false;
      flags.iron_untilt_normalization = false;
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
      flags.bg_global = false;
      break;
    case 6:
      flags.use_normalization = true;
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
      flags.bg_global = false;
      break;
    case 10:
      flags.use_background_correction = true;
      flags.bg_mas5 = false;
      flags.bg_rma = true;
      flags.bg_iron = false;
      flags.bg_global = false;
      break;
    case 11:
      flags.use_background_correction = true;
      flags.bg_mas5 = false;
      flags.bg_rma = false;
      flags.bg_iron = true;
      flags.bg_global = false;
      break;
    case 12:
      if (arg != NULL) 
      {
        flags.probe_tab_filename = h_strdup(arg);
        hattach(flags.probe_tab_filename, mempool);
      }
      break;
    case 13:
      flags.bioconductor_compatability = true;
      break;
    case 14:
      flags.normalize_before_bg = true;
      break;
    case 17:
      flags.output_log2 = true;
      break;
    case 21:
      flags.iron_global_scaling_normalization = true;
      flags.iron_untilt_normalization = false;
      break;
    case 22:
      flags.iron_global_scaling_normalization = false;
      flags.iron_untilt_normalization = false;
      break;
    case 23:
      flags.iron_weight_exponent = atof(arg);
      break;
    case 25:
      flags.iron_untilt_normalization = true;
      flags.iron_global_scaling_normalization = false;
      break;
    case 26:
      flags.iron_fit_both_x_y = true;
      break;
    case 27:
      flags.iron_fit_both_x_y = false;
      break;
    case 28:
      flags.iron_fit_window_frac = atof(arg);
      break;
    /* --proteomics */
    case 29:
      flags.use_background_correction = false;
      flags.output_log2 = false;
      flags.iron_global_scaling_normalization = true;
      flags.iron_untilt_normalization = false;
      flags.iron_weight_exponent = 0;
      flags.iron_fit_both_x_y = true;
      flags.iron_condense_training = true;
      flags.floor_to_min_non_zero = false;
      flags.floor_non_zero_to_one = false;
      break;
    /* --rnaseq */
    case 30:
      flags.use_background_correction = false;
      flags.output_log2 = false;
      flags.iron_global_scaling_normalization = false;
      flags.iron_untilt_normalization = true;
      flags.iron_weight_exponent = 0;
      flags.iron_fit_both_x_y = false;
      flags.iron_condense_training = true;
      flags.floor_to_min_non_zero = false;
      flags.floor_non_zero_to_one = false;
      break;
    case 131:
      flags.floor_to_min_non_zero = true;
      flags.floor_non_zero_to_one = false;
      break;
    case 132:
      flags.floor_to_min_non_zero = false;
      flags.floor_non_zero_to_one = false;
      break;
    case 133:
      flags.iron_condense_training = true;
      break;
    case 134:
      flags.iron_condense_training = false;
      break;
    case 135:
      flags.use_background_correction = true;
      flags.bg_mas5     = false;
      flags.bg_rma      = false;
      flags.bg_rma_both = false;
      flags.bg_iron     = false;
      flags.bg_global   = true;
      flags.use_mm_probe_subtraction  = false;
      break;
    case 136:
      flags.floor_to_min_non_zero = false;
      flags.floor_non_zero_to_one = true;
      break;
    /* --microarray */
    case 137:
      flags.use_background_correction = true;
      flags.bg_mas5                   = false;
      flags.bg_rma                    = true;
      flags.bg_rma_both               = false;
      flags.bg_iron                   = false;
      flags.bg_global                 = false;
      flags.use_mm_probe_subtraction  = false;
      flags.output_log2 = true;
      flags.iron_global_scaling_normalization = false;
      flags.iron_untilt_normalization = false;
      flags.iron_weight_exponent = 4;
      flags.iron_fit_both_x_y = false;
      flags.iron_condense_training = false;
      flags.floor_to_min_non_zero = false;
      flags.floor_non_zero_to_one = false;
      break;
    
    case 'd':
      directory = h_strdup(arg);
      hattach(directory, mempool);
      break;
    case 'm':
      flags.use_normalization = true;
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
    case 'S':
      flags.use_spikeins = true;
      if (arg != NULL)
      {
        flags.spikeins_filename = h_strdup(arg);
        hattach(flags.spikeins_filename, mempool);
      }
      break;
    case 'x':
      flags.use_exclusions = true;
      if (arg != NULL)
      {
        flags.exclusions_filename = h_strdup(arg);
        hattach(flags.exclusions_filename, mempool);
      }
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
