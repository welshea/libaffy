
/**************************************************************************
 *
 * Filename: print_flags.c
 * 
 * Purpose:  Print various run-time parameters
 *
 * Creation: September 5th, 2012
 *
 * Author:   Eric A. Welsh, modified from work of Andrew M. Hoerter
 *
 *
 * Update History
 * --------------
 * 09/05/12: Creation.  (EAW)
 * 01/10/13: renamed iron_linear_normalization flag to
 *           iron_global_scaling_normalization
 * 01/10/18: added iron_weight_exponent
 * 09/18/13: added iron_fit_both_x_y flag
 * 09/18/13: added iron_fit_window_frac (EAW)
 * 10/04/17: added iron_condense_training (EAW)
 * 06/01/18: added support for probeset exclusions during IRON training (EAW)
 * 09/12/18: added support for not scaling excluded probesets (EAW)
 * 03/13/19: added support for iron_ignore_noise (EAW)
 * 09/13/23: added support for iron_check_saturated (EAW)
 * 09/13/23: added support for iron_ignore_low (EAW)
 * 04/26/25: added support for median normalization (EAW)
 * 04/12/17: added support for normalization before bg-sub (EAW)
 *
 **************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "affy.h"

#define boolstr(b) (b) ? "Yes" : "No"

void print_flags(AFFY_COMBINED_FLAGS *f, char *output_file_name)
{
  assert(f != NULL);

  printf("General flags for this run:\n");
  printf("======================================\n");
  printf("CDF Directory:                       %s\n", f->cdf_directory);
  printf("Output filename:                     %s\n", output_file_name);

  printf("BG Correction (global override):     %s\n", 
         boolstr(f->use_background_correction));
  printf("Normalize before BG Correction:      %s\n", 
         boolstr(f->normalize_before_bg));
  printf("MAS5 BG Correction:                  %s\n", 
         boolstr(f->bg_mas5));
  printf("RMA BG Correction:                   %s\n", 
         boolstr(f->bg_rma));
  printf("RMA-like PM+MM BG Correction:        %s\n", 
         boolstr(f->bg_rma_both));
  printf("Use MM Probe BG Subtraction:         %s\n", 
         boolstr(f->use_mm_probe_subtraction));

  printf("Quantile normalization:              %s\n", 
         boolstr(f->use_quantile_normalization));

  printf("IRON normalization:                  %s ", 
         boolstr(f->use_pairwise_normalization));
  if (f->use_pairwise_normalization)
  {
    printf("(model file: %s)\n",
           f->pairwise_model_filename);
   }
  else
    printf("\n");

  printf("Mean normalization:                  %s ",
         boolstr(f->use_mean_normalization));
  if (f->use_mean_normalization)
    printf("(target: %f)\n",
           f->mean_normalization_target_mean);
  else
    printf("\n");

  printf("Median normalization:                %s ",
         boolstr(f->use_median_normalization));
  if (f->use_median_normalization)
    printf("(target: %f)\n",
           f->median_normalization_target_median);
  else
    printf("\n");
  
  printf("Tukey's Biweight probesets:          %s\n",
         boolstr(f->use_tukey_biweight));
  printf("Median polish probesets:             %s\n",
         boolstr(f->use_median_polish));
  printf("Output log2 probesets:               %s\n",
         boolstr(f->output_log2));

  printf("Floor to min non-zero per sample:    %s\n", 
         boolstr(f->floor_to_min_non_zero));
  printf("Floor non-zero to one:               %s\n", 
         boolstr(f->floor_non_zero_to_one));

  printf("Dump probe values:                   %s ",
         boolstr(f->dump_probe_values));
  if (f->dump_probe_values)
    printf("(filename: %s)\n",
           f->probe_filename);
  else
    printf("\n");

  printf("Bioconductor compat:                 %s\n", 
         boolstr(f->bioconductor_compatability));
  printf("Output present/absent:               %s\n",
         boolstr(f->output_present_absent));
  printf("Salvage corrupt CEL files:           %s\n",
         boolstr(f->output_present_absent));


  printf("\n");
  printf("MAS5 specific flags for this run:\n");
  printf("======================================\n");

  printf("Probeset scaling:                    %s ",
         boolstr(f->use_probeset_scaling));
  if (f->use_probeset_scaling)
    printf("(target: %f)\n",
         f->scale_target);
  else
    printf("\n");

  printf("Trimmed mean low:                    %f\n", f->trimmed_mean_low);
  printf("Trimmed mean high:                   %f\n", f->trimmed_mean_high);
  printf("Number of zones (K):                 %d\n", f->K); 
  printf("Smoothing parameter:                 %d\n", f->smooth); 
  printf("Noise frac parameter:                %f\n", f->NoiseFrac); 
  printf("Delta parameter:                     %f\n", f->delta); 
  printf("Contrast tau parameter:              %f\n", f->contrast_tau); 
  printf("Scale tau parameter:                 %f\n", f->scale_tau); 

  printf("\n");
  printf("RMA specific flags for this run:\n");
  printf("======================================\n");

#if UNSUPPORTED
  printf("Normalize AFFX:                      %s\n",
         boolstr(f->normalize_affx_probes));
#endif

  printf("Dump probe affinities:               %s ",
         boolstr(f->dump_probe_affinities));
  if (f->dump_probe_affinities)
    printf("(filename: %s)\n",
           f->affinities_filename);
  else
    printf("\n");

  printf("Dump expression means:               %s ",
         boolstr(f->dump_expression_means));
  if (f->dump_expression_means)
    printf("(filename: %s)\n", f->means_filename);
  else
    printf("\n");

  printf("Use saved affinities:                %s ",
         boolstr(f->use_saved_affinities));
  if (f->use_saved_affinities)
    printf("(filename: %s)\n",
           f->affinities_filename);
  else
    printf("\n");

  printf("Use saved means:                     %s ",
         boolstr(f->use_saved_means));
  if (f->use_saved_means)
    printf("(filename: %s)\n", f->means_filename);
  else
    printf("\n");

  printf("\n");
  printf("IRON specific flags for this run:\n");
  printf("======================================\n");
  printf("Use global scaling factors instead:  %s\n",
         boolstr(f->iron_global_scaling_normalization));
  printf("Use single line untilting instead:   %s\n",
         boolstr(f->iron_untilt_normalization));
  printf("Exclude potentially 16-bit saturated during training: %s\n",
         boolstr(f->iron_check_saturated));
  printf("Exclude reference values <= 1 during training: %s\n",
         boolstr(f->iron_ignore_low));
  printf("Exclude noise-level during training: %s\n",
         boolstr(f->iron_ignore_noise));
  printf("Exclude probesets during training:   %s\n",
         boolstr(f->use_exclusions));
  if (f->use_exclusions)
    printf("Exclusion probeset filename:         %s\n",
           f->exclusions_filename);
  if (f->use_spikeins)
    printf("Spikeins probeset filename:          %s\n",
           f->spikeins_filename);
  printf("Fit to both X and Y:                 %s\n",
         boolstr(f->iron_fit_both_x_y));
  printf("Condense identical X,Y:              %s\n",
         boolstr(f->iron_condense_training));
  printf("Pseudo-density exponent:             %f\n",
          f->iron_weight_exponent);
  printf("Window width fraction:               %f\n",
          f->iron_fit_window_frac);
  printf("\n");
}
