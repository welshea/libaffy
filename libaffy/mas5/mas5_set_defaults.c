
/**************************************************************************
 *
 * Filename:  mas5_set_defaults.c
 *
 * Purpose:   Given an allocated MAS5 flags structure, lay down some
 *            default values.
 *
 * Creation: 
 *
 * Author:    Steven Eschrich
 *
 * Copyright: Copyright (C) 2007, Moffitt Cancer Center.
 *            All rights reserved.
 *
 * Update History
 * --------------
 * 04/14/05: Imported/repaired from old libaffy (AMH)
 * 10/17/06: Added option for quantile normalization (AMH)
 * 09/16/10: Add support for P/A calls (EAW)
 * 10/14/10: Pairwise normalization (EAW)
 * 10/22/10: Use new AFFY_COMBINED_FLAGS instead of AFFY_MAS5_FLAGS (EAW)
 * 11/18/10: Initialize bioconductor_compatablity flag (EAW)
 * 01/18/13: Initialize IRON-specific flags (EAW)
 * 09/18/13: added iron_fit_both_x_y flag: better norm, may distrub rank order
 * 09/18/13: added iron_fit_window_frac (EAW)
 * 10/04/17: added iron_condense_training (EAW)
 * 06/01/18: added support for probeset exclusions during IRON training (EAW)
 * 09/14/18: added spikeins support, similar to exclusions (EAW)
 * 08/12/20: add cdf_filename (EAW)
 * 08/18/20: add flags to enable/disable iron or quantile probeset norm after
 *           probe norm (EAW)
 *
 **************************************************************************/

#include <affy_mas5.h>

void affy_mas5_set_defaults(AFFY_COMBINED_FLAGS *f)
{
  f->scale_target = 500;
  f->trimmed_mean_low = 0.02;
  f->trimmed_mean_high = 0.98;
  f->K = 16;
  f->smooth = 100;
  f->NoiseFrac = 0.5;
  f->delta = pow(2.0,-20);
  f->contrast_tau = 0.03;
  f->scale_tau = 10;
  f->cdf_directory = ".";
  f->cdf_filename = "";
  f->probe_filename = "probe-values.txt";
  f->dump_probe_values = false;
  f->output_present_absent = false;
  f->use_background_correction = true;
  f->use_mean_normalization = true;
  f->use_probeset_scaling = true;
  f->use_quantile_normalization = false;
  f->use_pairwise_normalization = false;
  f->floor_to_min_non_zero      = false;
  f->floor_non_zero_to_one      = true;
  f->pairwise_model_filename = "median.CEL";
  f->mean_normalization_target_mean = 500;
  f->bioconductor_compatability = false;
  f->bg_mas5 = true;
  f->bg_rma = false;
  f->bg_rma_both = false;
  f->bg_iron = false;
  f->use_mm_probe_subtraction = true;
  f->use_mm_probeset_subtraction = false;
  f->probe_tab_filename = "probe_tab.txt";
  f->use_tukey_biweight = true;
  f->use_median_polish = false;
  f->normalize_before_bg = false;
  f->bioconductor_compatability = false;
  f->output_log2 = false;
  f->iron_global_scaling_normalization = false;
  f->iron_fit_both_x_y = false;
  f->iron_weight_exponent = 4.0;
  f->iron_fit_window_frac     = 0.10;
  f->iron_condense_training = false;
  f->iron_ignore_noise = false;
  f->salvage_corrupt = false;
  f->use_exclusions = false;
  f->exclusions_filename = NULL;
  f->use_spikeins = false;
  f->spikeins_filename = NULL;
  f->ignore_chip_mismatch = false;
  f->normalize_probesets = false;
}
