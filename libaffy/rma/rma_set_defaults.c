
/**************************************************************************
 *
 * Filename: rma_set_defaults.c
 *
 * Purpose:  Fill out an RMA flags structure with default values.
 *
 * Creation: 
 *
 * Author:   Steven Eschrich
 *
 *
 * Update History
 * --------------
 * 04/14/05: Imported/repaired from old libaffy (AMH)
 * 11/07/06: Added option to ignore AFFX/control probes for normalization
 *           (AMH)
 * 04/10/07: Add support for incremental RMA (AMH)
 * 10/18/10: Pairwise normalization (AMH)
 * 10/22/10: Use new AFFY_COMBINED_FLAGS instead of AFFY_RMA_FLAGS
 * 11/19/10: Initialize bioconductor_compatablity flag (EAW)
 * 01/18/13: Initialize IRON-specific flags (EAW)
 * 09/18/13: added iron_fit_both_x_y flag: better norm, may distrub rank order
 * 09/18/13: added iron_fit_window_frac (EAW)
 * 10/04/17: added iron_condense_training (EAW)
 * 06/01/18: added support for probeset exclusions during IRON training (EAW)
 * 09/14/18: added spikeins support, similar to exclusions (EAW)
 * 08/12/20: add cdf_filename (EAW)
 * 08/18/20: add flags to enable/disable iron or quantile probeset norm after
 *           probe norm (EAW)
 * 09/13/23: added iron_check_saturated flag (EAW)
 * 09/13/23: added iron_ignore_low flag (EAW)
 *
 **************************************************************************/

#include <affy_rma.h>

void affy_rma_set_defaults(AFFY_COMBINED_FLAGS *f)
{
  f->use_background_correction         = true;
  f->use_normalization                 = true;
  f->normalize_affx_probes             = true;
  f->use_mean_normalization            = false;
  f->mean_normalization_target_mean    = 500;
  f->use_pairwise_normalization        = false;
  f->pairwise_model_filename           = "median.CEL";
  f->dump_probe_affinities             = false;
  f->dump_expression_means             = false;
  f->use_saved_affinities              = false;
  f->use_saved_means                   = false;
  f->affinities_filename               = "affinities.txt";
  f->means_filename                    = "mean-values.txt";
  f->probe_filename                    = "probe-values.txt";
  f->cdf_directory                     = ".";
  f->cdf_filename                      = "";
  f->bg_mas5                           = false;
  f->bg_rma                            = true;
  f->bg_rma_both                       = false;
  f->bg_iron                           = false;
  f->use_mm_probe_subtraction          = false;
  f->use_mm_probeset_subtraction       = false;
  f->use_rma_probeset_singletons       = false;
  f->probe_tab_filename                = "probe_tab.txt";
  f->use_tukey_biweight                = false;
  f->use_median_polish                 = true;
  f->normalize_before_bg               = false;
  f->reuse_affinities                  = false;
  f->bioconductor_compatability        = false;
  f->output_log2                       = true;
  f->iron_global_scaling_normalization = false;
  f->iron_fit_both_x_y                 = false;
  f->iron_weight_exponent              = 4.0;
  f->iron_fit_window_frac              = 0.10;
  f->iron_condense_training            = false;
  f->iron_check_saturated              = true;
  f->iron_ignore_low                   = true;
  f->iron_ignore_noise                 = false;
  f->salvage_corrupt                   = false;
  f->floor_to_min_non_zero             = false;
  f->floor_non_zero_to_one             = false;

  f->use_probeset_scaling              = false;
  f->use_quantile_normalization        = true;
  f->use_pairwise_normalization        = false;
  f->bg_mas5                           = false;
  f->bg_rma                            = true;
  f->bg_rma_both                       = false;
  f->bg_iron                           = false;
  f->use_mm_probe_subtraction          = false;
  f->use_mm_probeset_subtraction       = false;
  f->use_tukey_biweight                = false;
  f->use_median_polish                 = true;
  f->output_log2                       = true;
  f->salvage_corrupt                   = false;
  f->use_exclusions                    = false;
  f->exclusions_filename               = NULL;
  f->use_spikeins                      = false;
  f->spikeins_filename                 = NULL;
  f->ignore_chip_mismatch              = false;
  f->normalize_probesets               = false;
}
