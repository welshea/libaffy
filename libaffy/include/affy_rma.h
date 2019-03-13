
/**************************************************************************
 *
 * Filename: affy_rma.h
 *
 * Purpose:  RMA library declarations.
 *
 * Creation:
 *
 * Author:   Steven Eschrich
 *
 *
 * Update History
 * --------------
 * 04/08/05: Imported from old libaffy (AMH)
 * 11/07/06: Added option to ignore AFFX/control probes for normalization
 *           (AMH)
 * 04/06/07: Added incremental RMA support (AMH)
 * 02/26/08: Namespace cleanup, new error handling (AMH)
 * 06/05/08: Finish error handling cleanups (AMH)
 * 10/18/10: Pairwise normalization (AMH)
 * 10/22/10: Use new AFFY_COMBINED_FLAGS instead of AFFY_RMA_FLAGS (EAW)
 * 11/19/10: Pass flags to affy_median_polish() (EAW)
 * 04/06/11: HACK -- Added affy_illumina() entry point (EAW)
 * 03/13/19: add estimate_global_bg_sub() (EAW)
 *
 **************************************************************************/


#ifndef _AFFY_RMA_H_
#define _AFFY_RMA_H_

#include <math.h>
#include "affy.h"
#include "utils.h"

/** 
 * Structure to hold all relevant information
 * about running an instance of RMA. The defaults 
 * for these fields are in parenthesis, but the
 * ultimate authority for these values is defined
 * within the file rma/rma_set_defaults.c.
 */
typedef struct affy_rma_flag_struct
{
  /** (.) Default location to look for CDF file */
  char *cdf_directory;

	/** (true) Run background correction */
  bool use_background_correction;

  /** (true) Run normalization step */
  bool use_normalization;

  /** (true) Use AFFX (control) probes when normalizing */
  bool normalize_affx_probes;

  /** (false) Use mean normalization instead of quantile normalization */
  bool use_mean_normalization;

  /** (500) If using mean normalization, set mean equal to target */
  double mean_normalization_target_mean;

  /** (false) Use a pairwise normalization of all probes prior to processing.
   *          This option is mutually exclusive with other normalization
   *          options.
   */
  bool use_pairwise_normalization;

  /* Model chip to use for pairwise normalization. */
  char *pairwise_model_filename;

  /** (false) Dump probe affinities to a file */
  bool dump_probe_affinities;

  /** ("affinities.txt") Optional affinity filename */
  char *affinities_filename;

  /** (false) Dump mean values to a file */
  bool dump_expression_means;

  /** (false) Dump raw probe values to a file */
  bool dump_probe_values;

  /** ("probe-values.txt") Optional probe values filename */
  char *probe_filename;

  /** ("mean-values.txt") Optional mean values filename */
  char *means_filename;

  /** (false) Use previously stored probe affinities rather than
      calculating them as normal. */
  bool use_saved_affinities;

  /** (false) Use previously stored means rather than calculating them
      as normal. */
  bool use_saved_means;
} AFFY_RMA_FLAGS;

#ifdef __cplusplus
extern "C"
{
#endif

  void                 affy_rma_set_defaults(AFFY_COMBINED_FLAGS *f);
  AFFY_COMBINED_FLAGS *affy_rma_get_defaults(AFFY_ERROR *err);

/* 
 * General entrypoint for a typical RMA. A null pointer can be
 * passed for the flags argument (defaults are selected).
 */
  AFFY_CHIPSET *affy_rma(char **filelist, AFFY_COMBINED_FLAGS *f, AFFY_ERROR *err);
  AFFY_CHIPSET *affy_illumina(char **filelist, AFFY_COMBINED_FLAGS *f, AFFY_ERROR *err);

  /* Individual functions used in RMA */
  void affy_rma_background_correct(AFFY_CHIPSET *c, 
                                   unsigned int chipnum, 
                                   AFFY_ERROR *err);
  void affy_rma_background_correct_pm_mm_separately(AFFY_CHIPSET *c, 
                                                    unsigned int chipnum, 
                                                    AFFY_ERROR *err);
  void affy_rma_background_correct_pm_mm_together(AFFY_CHIPSET *c, 
                                                  unsigned int chipnum,
                                                  unsigned char pm_only,
                                                  AFFY_ERROR *err);
  void affy_rma_quantile_normalization_chip(AFFY_CHIPSET *c, 
					    int chipnum,
					    double *mean, 
					    AFFY_COMBINED_FLAGS *f,
                                            AFFY_ERROR *err);
  void affy_global_background_correct(AFFY_CHIPSET *c, 
                                      unsigned int chipnum, 
                                      AFFY_ERROR *err);
  void affy_global_background_correct_pm_only(AFFY_CHIPSET *c, 
                                      unsigned int chipnum, 
                                      AFFY_ERROR *err);
  void affy_rma_quantile_normalization_chipset(AFFY_CHIPSET *c, 
					       double *mean,
					       AFFY_COMBINED_FLAGS *f);
  void affy_rma_signal(AFFY_CHIPSET *c, AFFY_COMBINED_FLAGS *f,
                       int safe_to_write_affinities_flag, AFFY_ERROR *err);
  void affy_rma_median_polish(double **z, 
			      int startingprobe, 
			      int startingchip, 
			      int numprobes, 
			      int numchips, 
			      double *results, 
			      double *affinities, 
			      double *t_val,
			      AFFY_COMBINED_FLAGS *f,
                              AFFY_ERROR *err);
  double estimate_global_bg_sub(double *pm, 
                                int n,
                                int already_logged_flag,
                                double *ret_log_peak, double *ret_log_sd,
                                AFFY_ERROR *err);


#ifdef __cplusplus
};
#endif

#endif
