
/******************************************************************************
 **
 ** file: rma_background.c
 **
 ** Copyright (C) 2002 - 2003 B. M. Bolstad
 ** 
 ** Written by: B. M. Bolstad  <bolstad@stat.berkeley.edu>
 ** Implementation dates: 2002-2003
 **
 ** 2012-03-20 added check for negative/zero sigma prior to taking sqrt (EAW)
 ** 2013-05-24 handle cases where max density point is below data range (EAW)
 ** 2013-05-28 do not adjust intensities that are zero to begin with (EAW)
 ** 2013-08-15 require points used in bg estimation to be > 0 intensity (EAW)
 ** 2014-03-11 sigma was being calculated off the unfiltered pm[],
 **            which was only a problem when there are intensities <= 0,
 **            so we never saw it until using it on such data (EAW)
 ** 2014-03-11 added additional functions for comparison to RNAExpress (EAW)
 ** 2014-03-12 updated to work with exon arrays (EAW)
 ** 2015-02-25 added functions for simple global bg subtraction (EAW)
 ** 2018-09-14 added pm-only global background correction (EAW)
 ** 2018-10-29 fix affy_rma_background_correct to handle 1:many probe:probeset
 ** 2019-03-13 changed calling convention for estimate_global_bg_sub() (EAW)
 ** 2019-03-13 modified estimate_global_bg_sub() behavior a bit (EAW)
 */

#define TINY_VALUE 1E-16

#include <affy_rma.h>
#include <float.h>

/* Local functions */

/* phi- compute the standard normal density. */
static INLINE double phi(double x)
{
  return (1 / sqrt(2 * AFFY_PI) * exp(-0.5 * x * x));
}

/* Phi - Compute the standard normal distribution function */
static INLINE double Phi(double x)
{
  return affy_pnorm5(x, 0.0, 1.0, 1, 0);
}

/***************************************************************
 **
 ** double estimate_alpha(double *PM,double PMmax, int rows,int
 ** cols,int column)
 **
 ** estimate the alpha parameter given vector x value of maximum of density
 ** of x, dimensions of x matrix and column of interest using method proposed
 ** in affy2
 **
 **
 ***************************************************************/
static double estimate_alpha(double *x, 
                             unsigned int length, 
                             double max, 
                             AFFY_ERROR *err)
{
  double alpha;
  int    i;

  assert(x != NULL);

  /* Adjust x by the max value */
  for (i = 0; i < length; i++)
    x[i] -= max;

  alpha = affy_max_density(x, length, err);
  AFFY_CHECK_ERROR(err, 0.0);

  alpha = 1.0 / alpha;

  return (alpha);
}

/* old algorithm from RMAExpress affy_1.1.1 */
static double estimate_alpha_2002(double *x, 
                                  unsigned int length, 
                                  double max, 
                                  AFFY_ERROR *err)
{
  double alpha, tmpsum = 0.0;
  int    i, numtop = 0;

  assert(x != NULL);

  for (i = 0; i < length; i++)
  {
    if (x[i] > max)
    {
      tmpsum += x[i] - max;
      numtop++;
    }
  }

  alpha = numtop / tmpsum;

  return (alpha);
}

/***************************************************************
 **
 ** double get_sd()
 ** 
 ** estimate the sigma parameter given vector x, value of maximum 
 ** of density of x, and length of x
 ***************************************************************/
static double get_sd(double *x, unsigned int n, double max)
{
  double sigma, tmpsum = 0.0;
  int    i, numtop = 0;

  assert(x != NULL);

  for (i = 0; i < n; i++)
  {
    if (x[i] < max)
    {
      tmpsum += (x[i] - max) * (x[i] - max);
      numtop++;
    }
  }

  sigma = tmpsum / (numtop - 1);
  
  /* will result in NaNs later if sigma <= 0 */
  if (sigma <   0)     sigma = 0;
  if (sigma == -0)     sigma = 0;
  if (sigma < DBL_MIN) sigma = DBL_MIN;

  sigma = sqrt(sigma) * sqrt(2.0) / 0.85;

  return (sigma);
}

/* do not scale by sqrt(2)/0.85 */
static double get_sd_no_scale(double *x, unsigned int n, double max)
{
  double sigma, tmpsum = 0.0;
  int    i, numtop = 0;

  assert(x != NULL);

  for (i = 0; i < n; i++)
  {
    if (x[i] < max)
    {
      tmpsum += (x[i] - max) * (x[i] - max);
      numtop++;
    }
  }

  sigma = tmpsum / (numtop - 1);
  
  /* will result in NaNs later if sigma <= 0 */
  if (sigma <   0)     sigma = 0;
  if (sigma == -0)     sigma = 0;
  if (sigma < DBL_MIN) sigma = DBL_MIN;

  sigma = sqrt(sigma);

  return (sigma);
}

/* sd of log intensities */
static double get_sd_log(double *x, unsigned int n, double max)
{
  double sigma, tmpsum = 0.0;
  double log_max, diff;
  int    i, numtop = 0;

  assert(x != NULL);
  
  log_max = log(max);

  for (i = 0; i < n; i++)
  {
    if (x[i] < max)
    {
      /* log(x[i]) - log(max) */
      diff = log(x[i] / max);

      tmpsum += diff * diff;
      numtop++;
    }
  }

  sigma = tmpsum / (numtop - 1);
  
  /* will result in NaNs later if sigma <= 0 */
  if (sigma <   0)     sigma = 0;
  if (sigma == -0)     sigma = 0;
  if (sigma < DBL_MIN) sigma = DBL_MIN;

  sigma = sqrt(sigma);

  return (sigma);
}


double estimate_global_bg_sub(double *pm, 
                              int n,
                              int already_logged_flag,
                              double *ret_log_peak, double *ret_log_sd,
                              AFFY_ERROR *err)
{
  int    i = 0, numx = 0;
  double max, *x, min_bigger = 0;
  double max_unlog;
  double sd, diff, bg;

  x = h_calloc(n, sizeof(double));
  if (x == NULL)
    AFFY_HANDLE_ERROR_VOID_ZERO("calloc failed", AFFY_ERROR_OUTOFMEM, err);

  /* use unloged data for peak estimation */

  /* Get all values > 0 */
  numx = 0;
  for (i = 0; i < n; i++)
  {
    if (pm[i] > TINY_VALUE)
    {
      if (already_logged_flag)
        x[numx++] = exp(pm[i]);
      else
        x[numx++] = pm[i];
    }
  }
  /* if none are > 0, use the entire vector */
  if (numx == 0)
  {
    for (i = 0; i < n; i++)
    {
      if (already_logged_flag)
        x[numx++] = exp(pm[i]);
      else
        x[numx++] = pm[i];
    }
  }

  max_unlog = affy_max_density(x, numx, err);
  AFFY_CHECK_ERROR_VOID_ZERO(err);

  max = 0;
  if (max_unlog > 0)
    max = log(max_unlog);


  /* use log values for SD estimation */
  /* Get all values less than max */
  numx = 0;
  for (i = 0; i < n; i++)
  {
    if (pm[i] < max_unlog && pm[i] > TINY_VALUE)
    {
      if (already_logged_flag)
        x[numx++] = pm[i];
      else
        x[numx++] = log(pm[i]);
    }
  }
  /* If none were found, then get all values less than or equal to max */
  if (numx == 0)
  {
    for (i = 0; i < n; i++)
    {
      if (pm[i] <= max_unlog && pm[i] > TINY_VALUE)
      {
        if (already_logged_flag)
          x[numx++] = pm[i];
        else
          x[numx++] = log(pm[i]);
      }

      /* keep track of min_bigger, in case we still can't find anything */
      else if (pm[i] > 0 && (pm[i] < min_bigger || min_bigger == 0))
        min_bigger = pm[i];
    }
  }
  /* Uh oh, we still didn't find anything, so take drastic measures */
  if (numx == 0)
  {
    for (i = 0; i < n; i++)
    {
      if (pm[i] <= min_bigger && pm[i] > TINY_VALUE)
      {
        if (already_logged_flag)
          x[numx++] = pm[i];
        else
          x[numx++] = log(pm[i]);
      }
    }
  }


  /* calculate sd */
  sd = 0.0;
  for (i = 0; i < numx; i++)
  {
    diff = max - x[i];
    sd += diff * diff;
  }
  if (sd)
  {
    sd = sqrt(sd / numx);
  }


  bg = max - 2.0 * sd;
  bg = exp(bg) - 1.0;

  if (bg < 0) bg = 0;


  if (ret_log_peak) *ret_log_peak = max;
  if (ret_log_sd)   *ret_log_sd   = sd;
  
  return bg;

cleanup:
  if (x)
    h_free(x);
  
  return 0;
}

static void estimate_bg_parameters(double *pm, 
                                   int n, 
                                   double *alpha, 
				   double *mu, 
                                   double *sigma,
                                   AFFY_ERROR *err)
{
  int    i = 0, numx = 0;
  double max, *x, min_bigger = 0;

  x = h_calloc(n, sizeof(double));
  if (x == NULL)
    AFFY_HANDLE_ERROR_VOID("calloc failed", AFFY_ERROR_OUTOFMEM, err);

  /* Get all values > 0 */
  numx = 0;
  for (i = 0; i < n; i++)
    if (pm[i] > TINY_VALUE)
      x[numx++] = pm[i];

  /* if none are > 0, use the entire vector */
  if (numx == 0)
    max = affy_max_density(pm, n, err);
  /* otherwise, use only the points > 0 */
  else
    max = affy_max_density(x, numx, err);
  AFFY_CHECK_ERROR_VOID(err);

  /* Get all values less than max */
  numx = 0;
  for (i = 0; i < n; i++)
  {
    if (pm[i] < max && pm[i] > TINY_VALUE)
      x[numx++] = pm[i];
  }
  
  /* If none were found, then get all values less than or equal to max */
  if (numx == 0)
  {
    for (i = 0; i < n; i++)
    {
      if (pm[i] <= max && pm[i] > TINY_VALUE)
        x[numx++] = pm[i];

      /* keep track of min_bigger, in case we still can't find anything */
      else if (pm[i] > 0 && (pm[i] < min_bigger || min_bigger == 0))
        min_bigger = pm[i];
    }
  }
  
  /* Uh oh, we still didn't find anything, so take drastic measures */
  if (numx == 0)
  {
    for (i = 0; i < n; i++)
    {
      if (pm[i] <= min_bigger && pm[i] > TINY_VALUE)
        x[numx++] = pm[i];
    }
  }

  /* Get max density and sd */
  max = affy_max_density(x, numx, err);
  AFFY_CHECK_ERROR_GOTO(err, cleanup);

  *sigma = get_sd(x, numx, max) * 0.85;

  /* Get all values greater than max */
  numx = 0;
  for (i = 0; i < n; i++)
  {
    if (pm[i] > max)
      x[numx++] = pm[i];
  }

  /* the 0.85 is to fix up constant in above */
  *alpha = estimate_alpha(x, numx, max, err);
  *mu = max;

cleanup:
  if (x)
    h_free(x);
}

/* old functions, from RMAExpress affy_1.1.1 */
static void estimate_bg_parameters_2002(double *pm, 
                                        int n, 
                                        double *alpha, 
				        double *mu, 
                                        double *sigma,
                                        AFFY_ERROR *err)
{
  int    i = 0, numx = 0;
  double max, *x, min_bigger = 0;

  x = h_calloc(n, sizeof(double));
  if (x == NULL)
    AFFY_HANDLE_ERROR_VOID("calloc failed", AFFY_ERROR_OUTOFMEM, err);

  /* Get all values > 0 */
  numx = 0;
  for (i = 0; i < n; i++)
    if (pm[i] > TINY_VALUE)
      x[numx++] = pm[i];

  /* if none are > 0, use the entire vector */
  if (numx == 0)
  {
    max = affy_max_density(pm, n, err);
    *sigma = get_sd(pm, n, max);
    *alpha = estimate_alpha_2002(pm, n, max, err);
  }
  /* otherwise, use only the points > 0 */
  else
  {
    max = affy_max_density(x, numx, err);
    *sigma = get_sd(x, numx, max);
    *alpha = estimate_alpha_2002(x, numx, max, err);
  }
  AFFY_CHECK_ERROR_VOID(err);

  *mu = max;

cleanup:
  if (x)
    h_free(x);
}

/*
 * Background correct the values in pm 
 */
void affy_rma_background_correct(AFFY_CHIPSET *c, 
                                 unsigned int chipnum, 
                                 AFFY_ERROR *err)
{
  double             b, alpha, mu, sigma, *pm;
  double            *pm_nodupes = NULL;
  int                j, n;
  int                p, n2 = 0, x, y;
  LIBUTILS_PB_STATE  pbs;

  AFFY_CDFFILE      *cdf;

  assert(c                    != NULL);
  assert(c->cdf               != NULL);
  assert(c->chip              != NULL);
  assert(c->chip[chipnum]->pm != NULL);
  assert(chipnum < c->num_chips);

  cdf        = c->cdf;
  n          = cdf->numprobes;
  pm         = c->chip[chipnum]->pm;

  /* deal with duplicate probes */
  if (cdf->dupe_probes_flag)
  {
    pm_nodupes = h_calloc(n, sizeof(double));
    if (pm_nodupes == NULL)
      AFFY_HANDLE_ERROR_GOTO("calloc failed", AFFY_ERROR_OUTOFMEM, err, cleanup);

    memset(cdf->seen_xy[0], 0, cdf->numrows*cdf->numcols*sizeof(affy_uint8));
    for (p = 0, n2 = j = 0; p < n; p++)
    {
      x = cdf->probe[p]->pm.x;
      y = cdf->probe[p]->pm.y;

      if (cdf->seen_xy[x][y] == 0)
        pm_nodupes[j++] = pm[p];
      cdf->seen_xy[x][y] = 1;
    }
    n2 = j;
  }

  pb_init(&pbs);
  pb_begin(&pbs, 2, "RMA Background correction");
  pb_tick(&pbs,1, "Estimating background parameters");

  if (cdf->dupe_probes_flag)
  {
    estimate_bg_parameters(pm_nodupes, n2, &alpha, &mu, &sigma, err);
  }
  else
  {
    estimate_bg_parameters(pm, n, &alpha, &mu, &sigma, err);
    /* estimate_bg_parameters_2002(pm, n, &alpha, &mu, &sigma, err); */
  }

  AFFY_CHECK_ERROR_VOID(err);

  /* Adjust values. */
  b = mu + (alpha * sigma * sigma);

  pb_tick(&pbs,1,"Calculating PM values");
  for (j = 0; j < n; j++)
  {
    double a = pm[j] - b;

    if (pm[j] >= TINY_VALUE)
    {
      pm[j] = a + sigma * phi(a / sigma) / Phi(a / sigma);
/*
      pm[j] = a + sigma *
              (phi(a / sigma) - phi((pm[j] - a) / sigma)) /
              (Phi(a / sigma) + Phi((pm[j] - a) / sigma) - 1.0);
*/
    }
    
    if (pm[j] < TINY_VALUE)
      pm[j] = 0;
  }

  pb_finish(&pbs, "Finished background correction");

cleanup:
  if (pm_nodupes)
    h_free(pm_nodupes);
}

/*
 * Background correct the values, pm and mm in their original structures
 *
 * Treat them as one single distribution
 */
void affy_rma_background_correct_pm_mm_together(AFFY_CHIPSET *c, 
                                                unsigned int chipnum,
                                                unsigned char pm_only,
                                                AFFY_ERROR *err)
{
  double            b, alpha, mu, sigma;
  double            *mmpm = NULL;
  int               p, j, n, n2, x, y;
  LIBUTILS_PB_STATE pbs;
  AFFY_CDFFILE      *cdf;
  AFFY_CELFILE      *cel;

  assert(c                           != NULL);
  assert(c->cdf                      != NULL);
  assert(c->cdf->probe               != NULL);
  assert(chipnum < c->num_chips);
  assert(c->chip                     != NULL);
  assert(c->chip[chipnum]            != NULL);
  assert(c->chip[chipnum]->cel       != NULL);
  assert(c->chip[chipnum]->cel->data != NULL);

  cdf = c->cdf;
  cel = c->chip[chipnum]->cel;
  n   = cdf->numprobes;
  n2  = 2*n;

  pb_init(&pbs);
  if (pm_only)
    pb_begin(&pbs, 2, "RMA Background correction");
  else
    pb_begin(&pbs, 2, "RMA Background correction (PM/MM, together)");
  pb_tick(&pbs, 1, "Estimating background parameters");

  mmpm = h_calloc(n2, sizeof(double));
  if (mmpm == NULL)
    AFFY_HANDLE_ERROR_GOTO("calloc failed", AFFY_ERROR_OUTOFMEM, err, cleanup);

  /* store PM and MM into single working array */
  memset(cdf->seen_xy[0], 0, cdf->numrows*cdf->numcols*sizeof(affy_uint8));
  for (p = 0, n2 = j = 0; p < n; p++)
  {
    x = cdf->probe[p]->pm.x;
    y = cdf->probe[p]->pm.y;

    if (cdf->seen_xy[x][y] == 0)
      mmpm[j++] = cel->data[x][y].value;
    cdf->seen_xy[x][y] = 1;
    
    if (pm_only)
      continue;

    /* hack for missing MM probes, where MM coords == PM coords
     */
    if (cdf->probe[p]->pm.x == cdf->probe[p]->mm.x &&
        cdf->probe[p]->pm.y == cdf->probe[p]->mm.y)
    {
      continue;
    }
    else
    {
      x = cdf->probe[p]->mm.x;
      y = cdf->probe[p]->mm.y;

      if (cdf->seen_xy[x][y] == 0)
        mmpm[j++] = cel->data[x][y].value;
      cdf->seen_xy[x][y] = 1;
    }
  }
  n2 = j;


  /* *** background correction *** */

  estimate_bg_parameters(mmpm, n2, &alpha, &mu, &sigma, err);
/*  estimate_bg_parameters_2002(mmpm, n2, &alpha, &mu, &sigma, err); */
  AFFY_CHECK_ERROR_VOID(err);

  /* Adjust values. */
  b = mu + (alpha * sigma * sigma);

  if (pm_only)
    pb_tick(&pbs,1,"Calculating PM values");
  else
    pb_tick(&pbs,1,"Calculating PM+MM values");

  for (p = 0; p < n2; p++)
  {
    double a = mmpm[p] - b;

    if (mmpm[p] >= TINY_VALUE)
      mmpm[p] = a + sigma * phi(a / sigma) / Phi(a / sigma);
    
    if (mmpm[p] < TINY_VALUE)
      mmpm[p] = 0;
  }


  /* store corrected values back into original data structures */
  memset(cdf->seen_xy[0], 0, cdf->numrows*cdf->numcols*sizeof(affy_uint8));
  for (p = 0, j = 0; p < n; p++)
  {
    x = cdf->probe[p]->pm.x;
    y = cdf->probe[p]->pm.y;

    if (cdf->seen_xy[x][y] == 0)
      cel->data[x][y].value = mmpm[j++];
    cdf->seen_xy[x][y] = 1;

    if (pm_only)
      continue;

    /* hack for missing MM probes, where MM coords == PM coords
     */
    if (cdf->probe[p]->pm.x == cdf->probe[p]->mm.x &&
        cdf->probe[p]->pm.y == cdf->probe[p]->mm.y)
    {
      continue;
    }
    else
    {
      x = cdf->probe[p]->mm.x;
      y = cdf->probe[p]->mm.y;
      
      if (cdf->seen_xy[x][y] == 0)
        cel->data[x][y].value = mmpm[j++];
      cdf->seen_xy[x][y] = 1;
    }
  }


  pb_finish(&pbs, "Finished background correction");

cleanup:
  if (mmpm)
    h_free(mmpm);
}


/*
 * Background correct the values in pm by subtracting constant
 */
void affy_global_background_correct(AFFY_CHIPSET *c, 
                                    unsigned int chipnum,
                                    AFFY_ERROR *err)
{
  double            b;
  double            *mmpm = NULL;
  int               p, j, n, n2, x, y;
  LIBUTILS_PB_STATE pbs;
  AFFY_CDFFILE      *cdf;
  AFFY_CELFILE      *cel;
  unsigned char     pm_only;

  assert(c                           != NULL);
  assert(c->cdf                      != NULL);
  assert(c->cdf->probe               != NULL);
  assert(chipnum < c->num_chips);
  assert(c->chip                     != NULL);
  assert(c->chip[chipnum]            != NULL);
  assert(c->chip[chipnum]->cel       != NULL);
  assert(c->chip[chipnum]->cel->data != NULL);

  cdf     = c->cdf;
  cel     = c->chip[chipnum]->cel;
  n       = cdf->numprobes;
  n2      = 2*n;
  pm_only = cdf->no_mm_flag;

  pb_init(&pbs);
  pb_begin(&pbs, 2, "Global Background correction");
/*  pb_tick(&pbs, 1, "Estimating background parameters"); */

  mmpm = h_calloc(n2, sizeof(double));
  if (mmpm == NULL)
    AFFY_HANDLE_ERROR_GOTO("calloc failed", AFFY_ERROR_OUTOFMEM, err, cleanup);

  /* use only MM */
  memset(cdf->seen_xy[0], 0, cdf->numrows*cdf->numcols*sizeof(affy_uint8));
  for (p = 0, n2 = j = 0; p < n; p++)
  {
    x = cdf->probe[p]->mm.x;
    y = cdf->probe[p]->mm.y;

    if (cdf->seen_xy[x][y] == 0)
      mmpm[j++] = cel->data[x][y].value;
    cdf->seen_xy[x][y] = 1;
  }
  n2 = j;


  /* *** background correction *** */

  b = estimate_global_bg_sub(mmpm, n2, 0, NULL, NULL, err);
  AFFY_CHECK_ERROR_VOID(err);

  if (pm_only)
    pb_tick(&pbs,1,"Calculating PM values");
  else
    pb_tick(&pbs,1,"Calculating PM+MM values");


  /* store corrected values back into original data structures */
  memset(cdf->seen_xy[0], 0, cdf->numrows*cdf->numcols*sizeof(affy_uint8));
  for (p = 0, j = 0; p < n; p++)
  {
    x = cdf->probe[p]->pm.x;
    y = cdf->probe[p]->pm.y;

    if (cdf->seen_xy[x][y] == 0)
    {
      if (cel->data[x][y].value >= TINY_VALUE)
        cel->data[x][y].value -= b;
    
      if (cel->data[x][y].value < TINY_VALUE)
        cel->data[x][y].value = 0;
    }
    cdf->seen_xy[x][y] = 1;

    if (pm_only)
      continue;

    /* hack for missing MM probes, where MM coords == PM coords
     */
    if (cdf->probe[p]->pm.x == cdf->probe[p]->mm.x &&
        cdf->probe[p]->pm.y == cdf->probe[p]->mm.y)
    {
      continue;
    }
    else
    {
      x = cdf->probe[p]->mm.x;
      y = cdf->probe[p]->mm.y;
      
      if (cdf->seen_xy[x][y] == 0)
      {
        if (cel->data[x][y].value >= TINY_VALUE)
          cel->data[x][y].value -= b;
    
        if (cel->data[x][y].value < TINY_VALUE)
          cel->data[x][y].value = 0;
      }
      cdf->seen_xy[x][y] = 1;
    }
  }


  pb_finish(&pbs, "Finished background correction");

cleanup:
  if (mmpm)
    h_free(mmpm);
}


/*
 * Background correct the values in pm by subtracting constant
 */
void affy_global_background_correct_pm_only(AFFY_CHIPSET *c, 
                                            unsigned int chipnum,
                                            AFFY_ERROR *err)
{
  double            b, alpha, mu, sigma, *pm;
  int               j, n;
  LIBUTILS_PB_STATE pbs;

  assert(c                    != NULL);
  assert(c->cdf               != NULL);
  assert(c->chip              != NULL);
  assert(c->chip[chipnum]->pm != NULL);
  assert(chipnum < c->num_chips);

  n = c->cdf->numprobes;
  pb_init(&pbs);
  pm = c->chip[chipnum]->pm;
  pb_begin(&pbs, 2, "Global Background correction");
/*  pb_tick(&pbs,1, "Estimating background parameters"); */
/*  estimate_bg_parameters(pm, n, &alpha, &mu, &sigma, err); */
/*  estimate_bg_parameters_2002(pm, n, &alpha, &mu, &sigma, err); */
  AFFY_CHECK_ERROR_VOID(err);

  /* Adjust values. */
  b = estimate_global_bg_sub(pm, n, 0, NULL, NULL, err);

  pb_tick(&pbs,1,"Calculating PM values");

  for (j = 0; j < n; j++)
  {
    if (pm[j] >= TINY_VALUE)
    {
      pm[j] -= b;
    }
    
    if (pm[j] < TINY_VALUE)
      pm[j] = 0;
  }

  pb_finish(&pbs, "Finished background correction");
}
