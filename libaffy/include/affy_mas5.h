
/**************************************************************************
 *
 * Filename:  affy_mas5.h
 *
 * Purpose:   MAS5 declarations.
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
 * 04/08/05: Imported from old libaffy (AMH)
 * 10/17/06: Added option for quantile normalization (AMH)
 * 03/07/08: New error handling scheme (AMH)
 * 09/16/10: Add support for P/A calls (EAW)
 * 10/14/10: Pairwise normalization (EAW)
 * 10/22/10: Use new AFFY_COMBINED_FLAGS instead of AFFY_MAS5_FLAGS (EAW)
 * 08/12/20: add cdf_filename (EAW)
 *
 **************************************************************************/

#ifndef _AFFY_MAS5_H_
#define _AFFY_MAS5_H_

#include <math.h>
#include <utils.h>

#include "affy.h"

/*
 * Structure to hold all relevant information
 * about running an instance of MAS5.0. The
 * defaults for these fields are in parenthesis, but
 * the ultimate authority for these values is defined 
 * within the file mas5/mas5_set_defaults.c.
 */
typedef struct affy_mas5_flag_struct
{
  /** (.) Default location to look for CDF file */
  char *cdf_directory;

  /* full path to CDF file */
  char *cdf_filename;

  /** (true) Run MAS5.0 background correction */
  bool use_background_correction;

  /** (false) Use a mean normalization of all probes prior to processing.
   *          This option is mutually exclusive with 
   *          use_quantile_normalization and use_pairwise_normalization.
   */
  bool use_mean_normalization;

  /** (false) Use a quantile normalization of all probes prior to processing.
   *          This option is mutually exclusive with use_mean_normalization
   *          and use_pairwise_normalization.
   */
  bool use_quantile_normalization;

  /** (false) Use a pairwise normalization of all probes prior to processing.
   *          This option is mutually exclusive with use_mean_normalization
   *          and use_quantile_normalization.
   */
  bool use_pairwise_normalization;

  /* Model chip to use for pairwise normalization. */
  char *pairwise_model_filename;

  /** (500) Mean normalize all probes to target mean value */
  double mean_normalization_target_mean;

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

  /** (false) Dump raw probe values to a file */
  bool dump_probe_values;
  /** ("probe-values.txt") Optional probe values filename */
  char *probe_filename;

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

} AFFY_MAS5_FLAGS;

#ifdef __cplusplus
extern "C"
{
#endif

  AFFY_COMBINED_FLAGS *affy_mas5_get_defaults(AFFY_ERROR *err);
  void                 affy_mas5_set_defaults(AFFY_COMBINED_FLAGS *f);

  AFFY_CHIPSET *affy_mas5(char **filelist, AFFY_COMBINED_FLAGS *f, AFFY_ERROR *err);

  int affy_mas5_background_correction(AFFY_CHIPSET *c, 
                                      AFFY_COMBINED_FLAGS *f, 
                                      AFFY_ERROR *err);
  int  affy_mas5_signal(AFFY_CHIPSET *c, AFFY_COMBINED_FLAGS *f, AFFY_ERROR *err);
  int  affy_mas5_scale(AFFY_CHIPSET *c, AFFY_COMBINED_FLAGS *f, AFFY_ERROR *err);
  int  affy_mas5_scale_iron(AFFY_CHIPSET *c, AFFY_CHIPSET *model_chipset,
                            AFFY_COMBINED_FLAGS *f, AFFY_ERROR *err);
  int  affy_mas5_call(AFFY_CHIPSET *c, AFFY_COMBINED_FLAGS *f, AFFY_ERROR *err);
  char affy_mas5_pvalue_call(double pvalue);
  int  affy_mas5_subtract_mm_signal_probe(AFFY_CHIP *c,
                                          AFFY_COMBINED_FLAGS *f,
                                          AFFY_ERROR *err);
  int affy_iron_signal(AFFY_CHIPSET *c, AFFY_COMBINED_FLAGS *f,
                       AFFY_ERROR *err);


#ifdef __cplusplus
};
#endif

#endif
