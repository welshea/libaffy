
/**************************************************************************
 *
 * Filename:  mas5_wilcox.c
 *
 * Purpose:   assign P/M/A calls
 *
 * Creation:  01/10/08
 *
 * Author:    Eric A. Welsh
 *
 *
 * Update History
 * --------------
 * 09/16/10: initial version (EAW)
 * 09/20/10: Pooled memory allocator (AMH)
 * 12/22/11: replaced Abramowitz (1964) pnorm approximation with Hart (1968)
 *
 **************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "affy_wilcox.h"

#define TRUE 1

static int rcmp(double x, double y, signed char nalast)
{
  int nax = isnan(x), nay = isnan(y);

  if (nax && nay)
    return (0);
  if (nax)
    return (nalast ? 1 : -1);
  if (nay)
    return (nalast ? -1 : 1);
  if (x < y)
    return (-1);
  if (x > y)
    return (1);

  return (0);
}

static void rsort_with_index(double *x, int *indx, int n)
{
  double v;
  int    i, j, h, iv;

  for (h = 1; h <= n / 9; h = 3 * h + 1);
  for (; h > 0; h /= 3)
  {
    for (i = h; i < n; i++)
    {
      v = x[i];
      iv = indx[i];
      j = i;
      while (j >= h && rcmp(x[j - h], v, TRUE) > 0)
      {
        x[j] = x[j - h];
        indx[j] = indx[j-h];
        j -= h;
      }
      x[j] = v;
      indx[j] = iv;
    }
  }
}

/* Hart et al. 1968 approximation, taken from West 2009
 * "Better approximations to cumulative normal functions".
 *  Wilmott Magazine: 70.76.
 *   http://www.wilmott.com/pdfs/090721_west.pdf
 *
 * Should be roughly machine (double) precision.
 */
static double pnorm(double x)
{
  double cumnorm = 0;
  double xabs;
  double exponential;
  double build;

  xabs = fabs(x);
  
  if (xabs > 37)
  {
    cumnorm = 0;
  }
  else
  {
    exponential = exp(0.5 * -xabs*xabs);
    
    if (xabs < 7.07106781186547)
    {
      build = 3.52624965998911E-02 * xabs + 0.700383064443688;
      build = build * xabs + 6.37396220353165;
      build = build * xabs + 33.912866078383;
      build = build * xabs + 112.079291497871;
      build = build * xabs + 221.213596169931;
      build = build * xabs + 220.206867912376;
      
      cumnorm = exponential * build;
      
      build = 8.83883476483184E-02 * xabs + 1.75566716318264;
      build = build * xabs + 16.064177579207;
      build = build * xabs + 86.7807322029461;
      build = build * xabs + 296.564248779674;
      build = build * xabs + 637.333633378831;
      build = build * xabs + 793.826512519948;
      build = build * xabs + 440.413735824752;
      
      cumnorm = cumnorm / build;
    }
    else
    {
      build = xabs + 0.65;
      build = xabs + 4 / build;
      build = xabs + 3 / build;
      build = xabs + 2 / build;
      build = xabs + 1 / build;
      
      cumnorm = exponential / build / 2.506628274631;
    }
  }
  
  if (x > 0)
  {
    cumnorm = 1 - cumnorm;
  }
  
  return cumnorm;
}

/* Given a double array length nx, rank it, and put the results in 'r' */
static void rank(double *x, int nx, double *r)
{
  int i       = 0;
  int rank    = 1;
  int ranksum = 1;
  int ntie    = 1;
  int prev    = 0;

  r[0] = 1.0;

  for (i = 1; i < nx; i++)
  {
    if (x[i] == x[prev])
    {
      ntie++;
      rank++;
      ranksum += rank;
    }
    else
    {
      if (ntie > 1)
      {
        while (prev < i)
        {
          r[prev] = (double)ranksum / (double)ntie;
          prev++;
        }
      }

      rank++;
      ranksum = rank;
      r[i] = rank;
      prev = i;
      ntie = 1;
    }
  }

  if (ntie > 1)
  {
    while (prev < i)
    {
      r[prev] = (double)ranksum / (double)ntie;
      prev++;
    }
  }
}

/* a straight translation of relevant bits of the wilcox.test method
   in the R base library */
static double wilcox_approx(double *x, int n, double mu, AFFY_ERROR *err)
{
  int     i = 0, j = 0;
  double *r        = 0;
  double *absx     = 0;
  int    *xidx     = 0;
  double STATISTIC = 0;
  double NTIES_SUM = 0;
  int    prev      = 0;
  int    ntie      = 0;
  double z         = 0;
  double SIGMA     = 0;
  double PVAL      = 0;
  double nx        = n;
  int   *mempool;

  mempool = h_malloc(sizeof(int));
  if (mempool == NULL)
    AFFY_HANDLE_ERROR("malloc failed", AFFY_ERROR_OUTOFMEM, err, -DBL_MIN);

  for (i = 0; i < nx; i++)
  {
    x[j] = x[i] - mu;
    if (x[j] != 0)
      j++; /* eliminate zeros */
  }

  nx = j;
  r  = (double *)h_subcalloc(mempool, nx, sizeof(double));
  if (r == NULL)
  {
    h_free(mempool);
    AFFY_HANDLE_ERROR("calloc failed", AFFY_ERROR_OUTOFMEM, err, -DBL_MIN);
  }

  absx = (double *)h_subcalloc(mempool, nx, sizeof(double));
  if (absx == NULL)
  {
    h_free(mempool);
    AFFY_HANDLE_ERROR("calloc failed", AFFY_ERROR_OUTOFMEM, err, -DBL_MIN);
  }

  xidx = (int *)h_subcalloc(mempool, nx, sizeof(int));
  if (xidx == NULL)
  {
    h_free(mempool);
    AFFY_HANDLE_ERROR("calloc failed", AFFY_ERROR_OUTOFMEM, err, -DBL_MIN);
  }

  for(i = 0 ; i < nx; i++)
  {
    absx[i] = fabs(x[i]);
    xidx[i] = i;
  }

  rsort_with_index(absx, xidx, nx);
  rank(absx, nx, r);

  for(i = 0; i < nx; i++)
    r[i] = (x[xidx[i]] > 0) ? r[i] : -r[i];

  for(i = 0; i < nx; i++)
  {
    if (r[i] > 0)
      STATISTIC += r[i];
  }

  for (i = 1; i < nx; i++)
  {
    if (r[prev] == r[i])
    {
      ntie++;
    }
    else
    {
      if (ntie > 1)
        NTIES_SUM += ntie * ntie * ntie - ntie;

      ntie = 0;
      prev = i;
    }
  }

  NTIES_SUM += ntie * ntie * ntie - ntie; /* added by Crispin Noc 2005 */

  z       = STATISTIC - (nx * (nx + 1))/4;
  SIGMA   = sqrt((nx * (nx + 1) * (2 * nx + 1)) / 24 - (NTIES_SUM / 48));
  PVAL    = pnorm(z / SIGMA);
  PVAL    = 1 - PVAL;

  h_free(mempool);

  return (PVAL);
}

static int compare_abs_r(const void *keyval, const void *datum)
{
  struct affy_wilcox *ptr1, *ptr2;

  ptr1 = * (struct affy_wilcox **) keyval;
  ptr2 = * (struct affy_wilcox **) datum;

  if (ptr1->abs_r < ptr2->abs_r)
    return (-1);
  if (ptr1->abs_r > ptr2->abs_r)
    return (1);

  return (0);
}

static void assign_ranks(struct affy_wilcox **rset_sort, int n)
{
  int i, j;
  int rank      = 1;
  int tie_sum   = 0;
  int tie_count = 0;

  for (i = 0; i < n; i++, rank++)
  {
    /* check for tie in the future */
    if (i < n - 1 && rset_sort[i]->abs_r == rset_sort[i+1]->abs_r)
    {
      tie_sum += rank;
      tie_count++;
    }
    /* end of a tie */
    else if (tie_sum)
    {
      double tied_rank;

      tie_sum += rank;
      tie_count++;

      tied_rank = (double)tie_sum / (double)tie_count;

      for (j = i - tie_count + 1; j <= i; j++)
        rset_sort[j]->rank = tied_rank;

      tie_sum = 0;
      tie_count = 0;
    }
    else
    {
      rset_sort[i]->rank = rank;
    }
  }
}

static void recurse_sum(struct affy_wilcox *rset, int start, int stop,
                        double sum, double *pvalue, double S)
{
  /* negative rank */
  if (start == stop)
  {
    if (sum > S)
      *pvalue += 1.0;
    else if (sum == S)
      *pvalue += 0.5;
  }
  else
  {
    recurse_sum(rset, start+1, stop, sum, pvalue, S);
  }

  /* positive rank */
  sum += rset[start].rank;
  if (start == stop)
  {
    if (sum > S)
      *pvalue += 1.0;
    else if (sum == S)
      *pvalue += 0.5;
  }
  else
  {
    recurse_sum(rset, start+1, stop, sum, pvalue, S);
  }
}

double affy_mas5_calculate_wilcox_pvalue(struct affy_wilcox *rset, int n)
{
  int    combinations, i;
  double S = 0, pvalue = 0;

  for (i = 0; i < n; i++)
  {
    if (rset[i].r > 0)
      S += rset[i].rank;
  }

  combinations = pow(2, n);

  recurse_sum(rset, 0, n-1, 0.0, &pvalue, S);

  return (pvalue / combinations);
}

/* assume there are no zero points */
double affy_mas5_calculate_call_pvalue(double *values,
                                       int n,
                                       double tau,
                                       AFFY_ERROR *err)
{
  struct affy_wilcox  *rset      = NULL;
  struct affy_wilcox **rset_sort = NULL;
  int                  i, *mempool;
  double               pvalue;

  if (n == 0)
    return (1.0);

  /* Should be at least >= 20, set it as high as is feasible */
  /* HG-U133plus2 chip has one AFFX probeset with 69, 2nd biggest are 20 */
  if (n >= 21)
    return (wilcox_approx(values, n, tau, err));

  mempool = h_malloc(sizeof(int));
  if (mempool == NULL)
    AFFY_HANDLE_ERROR("malloc failed", AFFY_ERROR_OUTOFMEM, err, -DBL_MIN);

  rset = (struct affy_wilcox *)h_subcalloc(mempool,
                                           n,
                                           sizeof(struct affy_wilcox));
  if (rset == NULL)
  {
    h_free(mempool);
    AFFY_HANDLE_ERROR("calloc failed", AFFY_ERROR_OUTOFMEM, err, -DBL_MIN);
  }

  rset_sort = (struct affy_wilcox **)h_subcalloc(mempool,
                                                 n,
                                                 sizeof(struct affy_wilcox *));
  if (rset_sort == NULL)
  {
    h_free(mempool);
    AFFY_HANDLE_ERROR("calloc failed", AFFY_ERROR_OUTOFMEM, err, -DBL_MIN);
  }

  for (i = 0; i < n; i++)
    rset_sort[i] = &rset[i];

  for (i = 0; i < n; i++)
    rset[i].r = values[i] - tau;

  for (i = 0; i < n; i++)
    rset[i].abs_r = fabs(rset[i].r);

  qsort(rset_sort, n, sizeof(struct affy_wilcox *), compare_abs_r);

  assign_ranks(rset_sort, n);

  pvalue = affy_mas5_calculate_wilcox_pvalue(rset, n);

  h_free(mempool);

  return (pvalue);
}
