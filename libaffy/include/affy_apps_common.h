/**************************************************************************
 *
 * Filename:  affy_apps_common.h
 *
 * Purpose:   MAS5 and RMA option structures, so they can cross-communicate
 *
 * Creation:
 *
 * Author:    Eric A. Welsh
 *
 * Copyright: Copyright (C) 2010, Moffitt Cancer Center.
 *            All rights reserved.
 *
 * Update History
 * --------------
 * 10/22/10: Initial version (EAW)
 * 01/10/13: renamed iron_linear_normalization flag to
 *           iron_global_scaling_normalization
 * 04/18/13: added iron_untilt_normalization, for fitting a straight line
 * 09/18/13: added iron_fit_both_x_y flag: better norm, may distrub rank order
 * 09/18/13: added iron_fit_window_frac (EAW)
 * 10/04/17: added iron_condense_training (EAW)
 * 06/01/18: added support for probeset exclusions during IRON training (EAW)
 * 03/13/19: added iron_ignore_noise (EAW)
 * 05/22/19: added ignore_chip_mismatch flag (EAW)
 *
 **************************************************************************/

/*
 * Structure to hold all relevant information
 * about running an instance of MAS5.0 and/or RMA. The
 * defaults for these fields are in parenthesis, but
 * the ultimate authority for these values is defined 
 * within the files mas5/mas5_set_defaults.c and
 * rma/rma_set_defaults.c
 */
typedef struct affy_combined_flag_struct
{
  /* *** Options common to both MAS5 and RMA *** */

  /** (.) Default location to look for CDF file */
  char *cdf_directory;

  /** (true) Run MAS5.0 background correction */
  bool use_background_correction;

  /** (false) Use a mean normalization of all probes prior to processing.
   *          This option is mutually exclusive with 
   *          use_quantile_normalization and use_pairwise_normalization.
   */
  bool use_mean_normalization;

  /** (false) Use a pairwise normalization of all probes prior to processing.
   *          This option is mutually exclusive with use_mean_normalization
   *          and use_quantile_normalization.
   */
  bool use_pairwise_normalization;

  /* Model chip to use for pairwise normalization. */
  char *pairwise_model_filename;

  /** (500) Mean normalize all probes to target mean value */
  double mean_normalization_target_mean;

  /* floor final zero values to the minimum non-zero value per chip */
  double floor_to_min_non_zero;

  /* floor non-zero values to 1 */
  double floor_non_zero_to_one;

  /** (false) Dump raw probe values to a file */
  bool dump_probe_values;
  /** ("probe-values.txt") Optional probe values filename */
  char *probe_filename;

  /* output log2 probesets */
  bool output_log2;
  
  /** flags for MAS5 or RMA background subtraction */
  bool bg_mas5;
  bool bg_rma;
  bool bg_rma_both;
  bool bg_iron;
  bool bg_global;
  
  bool use_mm_probe_subtraction;
  bool use_mm_probeset_subtraction;
  
  /** probe_tab file */
  char *probe_tab_filename;
  
  bool use_median_polish;
  bool use_tukey_biweight;
  bool normalize_before_bg;
  
  bool salvage_corrupt;



  /* *** MAS5 specific options */

  /** (false) Use a quantile normalization of all probes prior to processing.
   *          This option is mutually exclusive with use_mean_normalization
   *          and use_pairwise_normalization.
   */
  bool use_quantile_normalization;

  /** (true) Scale probesets to a constant value, determined as scale_target */
  bool use_probeset_scaling;

  /** (500) If probeset scaling is used, this should be the target trimmed mean */
  double scale_target;
  /** (0.02) Percentage below which to remove from trimmed mean */
  double trimmed_mean_low;
  /** (0.98) Percentage above which to remove from trimmed mean */
  double trimmed_mean_high;
  
  /** (false) Use Bioconductor compatability mode */
  bool bioconductor_compatability;

  /** (false) Output present/absent calls in expression output file */
  bool output_present_absent;

  /** (16) The number of rectangular zones on the chip */
  int K; 
  /** (100) MAS5.0 smooth parameter */
  int smooth;
  /** (0.5) MAS5.0 noiseFrac parameter */
  double NoiseFrac;
  /** (2.0e-20) MAS5.0 delta parameter */
  double delta;
  /** (0.03) MAS5.0 Contrast tau parameter */
  double contrast_tau;
  /** (10) MAS5.0 Scale tau parameter */
  double scale_tau;



  /* *** RMA specific options *** */

  /** (true) Run normalization step */
  bool use_normalization;

  /** (true) Use AFFX (control) probes when normalizing */
  bool normalize_affx_probes;

  /** (false) Dump probe affinities to a file */
  bool dump_probe_affinities;

  /** ("affinities.txt") Optional affinity filename */
  char *affinities_filename;

  /** (false) Dump mean values to a file */
  bool dump_expression_means;

  /** ("mean-values.txt") Optional mean values filename */
  char *means_filename;

  /** (false) Use previously stored probe affinities rather than
      calculating them as normal. */
  bool use_saved_affinities;

  /** (false) Use previously stored means rather than calculating them
      as normal. */
  bool use_saved_means;
  
  /** (false) Peform probeset summarization at the single-chip level */
  bool use_rma_probeset_singletons;
  
  /** (false) Store and reuse calculated median polish parameters */
  bool reuse_affinities;

  /** (false) do not abort when multiple chip types are found
   ** This is needed for Decipher commercial FFPE prostate chips that
   ** contain incorret chip type information.
   ** Some claim to be HuEx-1_0-st-v2, while others claim to be Deciper.
   ** All chips are actually HuEx-1_0-st-v2.
   ** This option is a workaround for faulty CEL files.
   **
   ** Don't be surprised if it crashes whenever the chips truly are different.
   **/
   
  bool ignore_chip_mismatch;  
  
  
  /* IRON specific options */
  bool   iron_global_scaling_normalization;
  bool   iron_fit_both_x_y;
  bool   iron_untilt_normalization;
  bool   iron_condense_training;
  bool   iron_ignore_noise;
  double iron_weight_exponent;
  double iron_fit_window_frac;
  
  /* currently IRON-specific, but may expand to other methods eventually */
  bool   use_exclusions;
  bool   use_spikeins;
  char   *exclusions_filename;
  char   *spikeins_filename;

} AFFY_COMBINED_FLAGS;
