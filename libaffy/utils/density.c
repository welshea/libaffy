
/**************************************************************************
 *
 * Filename:  density.c
 *
 * Purpose:
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
 * 02/26/08: Moved/merged from rma/X_density.c (AMH)
 * 03/07/08: New error handling scheme (AMH)
 * 09/20/10: Pooled memory allocator (AMH)
 * 03/11/14: Added unweighted_massdist() function (EAW)
 *
 **************************************************************************/

#include <affy.h>
#ifdef SunOS
# include <ieeefp.h>
#endif

#define INVERSE_FFT 1
#define FFT         2

static int density_estimate_points = 16384;

/*** Old header from original source ***/
/*****************************************************************************
 **
 ** file: kernel_density.c
 **
 ** aim : compute weighted kernel density estimates
 **
 ** Copyright (C) 2003 Ben Bolstad
 **
 ** created on: Mar 9, 2003
 **
 ** Description
 ** 
 ** the aim here is to implement kernel density estimators, with the
 ** option to weight each observation. For speed we will use the FFT
 ** to convolve a weighted histogram with a kernel.
 ** 
 **
 ** History
 **
 ** Mar 9,  2003 - Initial version
 ** Mar 10, 2003 - Add in FFT framework, Push into AffyExtensions
 ** Mar 11, 2003 - Add ability to do kernel density with arbitrary weights
 ** Apr 21, 2003 - changes to work with RMAExpress
 ** Apr 22, 2003 - fix computation of bandwidth, add in linear
 **                interpolation step
 **
 ****************************************************************************/


/*********************************************************************
 ** void twiddle(int N, int i, double tf_real, double tf_imag)
 **
 ** N length of data series
 ** tf_real/tf_imag - on output contains real/imaginary part of twiddle factor 
 **
 ** twiddle factor in FFT 
 **
 ********************************************************************/
#define TWIDDLE(N, i, tf_real, tf_imag, what) \
   do \
   { \
      if ( i == 0 ) {\
	tf_real=1;\
	tf_imag=0;\
      } else {\
	tf_real=cos(2*AFFY_PI*(double)i/(double)N);  \
	tf_imag = sin(2*AFFY_PI*(double)i/(double)N); \
	if ( what == FFT ) tf_imag=-tf_imag;\
      }\
   } \
   while (0)

/* Private helper functions */
static double linear_interpolation(double v, double *x, double *y, int n);
static double bandwidth(double *x, int length, double iqr);
static double compute_sd(double *x, int length);
static void   fft_dif(double *f_real, double *f_imag, int p);
static void   fft_ditI(double *f_real, double *f_imag, int p);
static void   kernelize(double *data, int n, double bw, int kernel);
static void   fft_density_convolve(double *y, 
                                   double *kords, 
                                   int n, 
                                   AFFY_ERROR *err);
static void   weighted_massdist(double *x, 
                                int nx, 
                                double *w, 
                                double xlow,
				double xhigh, 
                                double *y, 
                                int ny);
static void   unweighted_massdist(double *x, 
                                int nx, 
                                double xlow,
				double xhigh, 
                                double *y, 
                                int ny);

/**********************************************************************
 **
 ** void kernel_density()
 **
 ** double *x - data vector
 ** int nx - length of x
 ** double *weights - a weight for each item of *x should be of length *nxxx
 ** double *output - place to output density values (x and Y)
 ** int N - length of output should be a power of two, preferably 512 or above
 **********************************************************************/
void affy_kernel_density(double *x, int nx, double *weights, double *dy,
			 double *dx, int N, AFFY_ERROR *err)
{

  int    i;
  double low, high, iqr, bw, from, to;
  double *kords;
  double *buffer;
  double *y;
  double *xords;
  int    *mempool;

  /* Allocate a memory pool handle */
  mempool = h_malloc(sizeof(int));
  if (mempool == NULL)
    AFFY_HANDLE_ERROR_VOID("malloc failed", AFFY_ERROR_OUTOFMEM, err);

  y = (double *)h_subcalloc(mempool, 2 * N, sizeof(double));
  if (y == NULL)
  {
    h_free(mempool);
    AFFY_HANDLE_ERROR_VOID("calloc failed", AFFY_ERROR_OUTOFMEM, err);
  }

  /* Calculate low/high points of data */
  buffer = h_suballoc(mempool, nx * sizeof(double));
  if (buffer == NULL)
  {
    h_free(mempool);
    AFFY_HANDLE_ERROR_VOID("malloc failed", AFFY_ERROR_OUTOFMEM, err);
  }

  for (i = 0; i < nx; i++)
  {
    buffer[i] = x[i];
  }

  qsort(buffer, nx, sizeof(double), dcompare);
  low  = buffer[0];
  high = buffer[nx - 1];
  iqr  = buffer[(int)(0.75 * nx + 0.5)] - buffer[(int)(0.25 * nx + 0.5)];
  bw   = bandwidth(x, nx, iqr);
  low  = low - 7 * bw;
  high = high + 7 * bw;

  kords = h_suballoc(mempool, 2 * N * sizeof(double));
  if (kords == NULL)
  {
    h_free(mempool);
    AFFY_HANDLE_ERROR_VOID("malloc failed", AFFY_ERROR_OUTOFMEM, err);
  }

  for (i = 0; i <= N; i++)
  {
    kords[i] = (double)i / (double)(2 * N - 1) * 2 * (high - low);
  }

  for (i = N + 1; i < 2 * N; i++)
  {
    kords[i] = -kords[2 * N - i];
  }

  kernelize(kords, 2 * N, bw, 2);

  weighted_massdist(x, nx, weights, low, high, y, N);
/*  unweighted_massdist(x, nx, low, high, y, N); */

  fft_density_convolve(y, kords, 2 * N, err);
  if (err->type != AFFY_ERROR_NONE)
  {
    h_free(mempool);

    return;
  }

  to  = high - 4 * bw;      /* corrections to get on correct output range */
  from = low + 4 * bw;

  xords = h_suballoc(mempool, N * sizeof(double));
  if (xords == NULL)
  {
    h_free(mempool);
    AFFY_HANDLE_ERROR_VOID("malloc failed", AFFY_ERROR_OUTOFMEM, err);
  }

  for (i = 0; i < N; i++)
  {
    xords[i] = (double)i / (double)(N - 1) * (high - low) + low;
    dx[i] = (double)i / (double)(N - 1) * (to - from) + from;
  }

  for (i = 0; i < N; i++)
  {
    kords[i] = kords[i] / (2 * N);
  }

  /* to get results that agree with R really need to do linear interpolation */
  for (i = 0; i < N; i++)
    dy[i] = linear_interpolation(dx[i], xords, kords, N);

  h_free(mempool);
}

/* This function was adapted from the Bioconductor RMA implementation */
/**************************************************************************
 **
 ** double max_density(double *z,int rows,int cols,int column, SEXP fn,SEXP rho)
 **
 ** double *z - matrix of dimension rows*cols
 ** int cols - matrix dimension
 ** int rows - matrix dimension
 ** int column - column of interest
 ** SEXP fn - R function for estimation of density
 ** SEXP rho - an R environment to work within
 **************************************************************************/
double affy_max_density(double *x, int n, AFFY_ERROR *err)
{
  int     i, imax;
  double *dx, *dy;
  double  final_result;
  double *weights;
  int    *mempool;

  /* Allocate a memory pool handle */
  mempool = h_malloc(sizeof(int));
  if (mempool == NULL)
    AFFY_HANDLE_ERROR("malloc failed", AFFY_ERROR_OUTOFMEM, err, 0);

  /* Assign standard weights */
  weights = h_subcalloc(mempool, n, sizeof(double));
  if (weights == NULL)
  {
    h_free(mempool);
    AFFY_HANDLE_ERROR("calloc failed", AFFY_ERROR_OUTOFMEM, err, 0);
  }

  for (i = 0; i < n; i++)
    weights[i] = 1.0;

  /* Space for output */
  dx = h_subcalloc(mempool, density_estimate_points, sizeof(double));
  if (dx == NULL)
  {
    h_free(mempool);
    AFFY_HANDLE_ERROR("calloc failed", AFFY_ERROR_OUTOFMEM, err, 0);
  }

  dy = h_subcalloc(mempool, density_estimate_points, sizeof(double));
  if (dy == NULL)
  {
    h_free(mempool);
    AFFY_HANDLE_ERROR("calloc failed", AFFY_ERROR_OUTOFMEM, err, 0);
  }

  affy_kernel_density(x, n, weights, dy, dx, density_estimate_points, err);
  if (err->type != AFFY_ERROR_NONE)
  {
    h_free(mempool);

    return(0.0);
  }

  /* Get the dx value corresponding to max dy */
  for (imax = 0, i = 1; i < density_estimate_points; i++)
    if (dy[imax] < dy[i])
      imax = i;
  final_result = dx[imax];

  h_free(mempool);

  return (final_result);
}

/*****************************************************************************
 **
 ** void weighted_massdist(double *x, int nx, double *w, double *xlow, double *xhigh, double *y, int *ny)
 **
 ** see AS R50 and AS 176  (AS = Applied Statistics)
 **
 ** idea is to discretize the data,  but have modified algorithm to put weights on each observation 
 **
 ** double *x - the data
 ** int nx - length of x
 ** double *w - weight for each one of x, vector should also be of length nx
 ** double *xlow - minimum value in x dimension
 ** double *xhigh - maximum value in x dimension
 ** double *y - on output will contain discretation scheme of data
 ** int ny - length of y
 **
 ****************************************************************************/
static void weighted_massdist(double *x, int nx, double *w, double xlow,
			      double xhigh, double *y, int ny)
{

  double fx, xdelta, xmass, xpos;
  int i, ix, ixmax, ixmin;

  ixmin = 0;
  ixmax = ny - 2;
  xmass = 0.0;
  xdelta = (xhigh - xlow) / (ny - 1);

  for (i = 0; i < ny; i++)
  {
    y[i] = 0.0;
  }

  for (i = 0; i < nx; i++)
  {
    xmass += w[i];
  }

  xmass = 1.0 / xmass;

  for (i = 0; i < nx; i++)
  {
    xpos = (x[i] - xlow) / xdelta;
    ix = (int)floor(xpos);
    fx = xpos - ix;
    if (ixmin <= ix && ix <= ixmax)
    {
      y[ix] += w[i] * (1 - fx);
      y[ix + 1] += w[i] * fx;
    }
    else if (ix == -1)
    {
      y[0] += w[i] * fx;
    }
    else if (ix == ixmax + 1)
    {
      y[ix] += w[i] * (1 - fx);
    }
  }

  for (i = 0; i < ny; i++)
    y[i] *= xmass;

}

/*****************************************************************************
 **
 ** void unweighted_massdist(double *x, int nx, double *w, double *xlow, double *xhigh, double *y, int *ny)
 **
 ** see AS R50 and AS 176  (AS = Applied Statistics)
 **
 ** idea is to discretize the data,  but have modified algorithm to put weights on each observation 
 **
 ** double *x - the data
 ** int nx - length of x
 ** double *w - weight for each one of x, vector should also be of length nx
 ** double *xlow - minimum value in x dimension
 ** double *xhigh - maximum value in x dimension
 ** double *y - on output will contain discretation scheme of data
 ** int ny - length of y
 **
 ****************************************************************************/
static void unweighted_massdist(double *x, int nx, double xlow,
			        double xhigh, double *y, int ny)
{

  double fx, xdelta, xpos;
  int i, ix, ixmax, ixmin;

  ixmin = 0;
  ixmax = ny - 2;
  xdelta = (xhigh - xlow) / (ny - 1);

  for (i = 0; i < ny; i++)
  {
    y[i] = 0.0;
  }

  for (i = 0; i < nx; i++)
  {
    xpos = (x[i] - xlow) / xdelta;
    ix = (int)floor(xpos);
    fx = xpos - ix;
    if (ixmin <= ix && ix <= ixmax)
    {
      y[ix] += (1 - fx);
      y[ix + 1] += fx;
    }
    else if (ix == -1)
    {
      y[0] += fx;
    }
    else if (ix == ixmax + 1)
    {
      y[ix] += (1 - fx);
    }
  }

  for (i = 0; i < ny; i++)
    y[i] /= (double) nx;
}

/*********************************************************************
 **
 ** void fft_dif(double *f_real, double *f_imag, int p){
 **
 ** compute the FFT using Decimation In Frequency of a data sequence of length 2^p
 **
 ** double *f_real - real component of data series
 ** double *f_imag - imaginary component of data series
 ** int p -  where 2^p is length of data series
 ** 
 ** computes the FFT in place, result is in reverse bit order.
 **
 ********************************************************************/
static void fft_dif(double *f_real, double *f_imag, int p)
{

  int BaseE, BaseO, i, j, k, Blocks, Points, Points2;
  double even_real, even_imag, odd_real, odd_imag;
  double tf_real, tf_imag;

  Blocks = 1;
  Points = 1 << p;

  for (i = 0; i < p; i++)
  {
    Points2 = Points >> 1;
    BaseE = 0;
    for (j = 0; j < Blocks; j++)
    {
      BaseO = BaseE + Points2;
      for (k = 0; k < Points2; k++)
      {
        even_real = f_real[BaseE + k] + f_real[BaseO + k];
        even_imag = f_imag[BaseE + k] + f_imag[BaseO + k];
        TWIDDLE(Points, k, tf_real, tf_imag, FFT);
        odd_real =
          (f_real[BaseE + k] - f_real[BaseO + k]) * tf_real -
          (f_imag[BaseE + k] - f_imag[BaseO + k]) * tf_imag;
        odd_imag =
          (f_real[BaseE + k] - f_real[BaseO + k]) * tf_imag +
          (f_imag[BaseE + k] - f_imag[BaseO + k]) * tf_real;
        f_real[BaseE + k] = even_real;
        f_imag[BaseE + k] = even_imag;
        f_real[BaseO + k] = odd_real;
        f_imag[BaseO + k] = odd_imag;

      }
      BaseE = BaseE + Points;
    }
    Blocks = Blocks << 1;
    Points = Points >> 1;
  }
}

/*********************************************************************
 **
 ** void fft_ditI(double *f_real, double *f_imag, int p){
 **
 ** compute the IFFT using Decimation In time of a data sequence of length 2^p
 **
 ** double *f_real - real component of data series
 ** double *f_imag - imaginary component of data series
 ** int p -  where 2^p is length of data series
 ** 
 ** computes the IFFT in place, where input is in reverse bit order.
 ** output is in normal order.
 **
 ********************************************************************/
static void fft_ditI(double *f_real, double *f_imag, int p)
{
  int i, j, k, Blocks, Points, Points2, BaseB, BaseT;
  double top_real, top_imag, bot_real, bot_imag, tf_real, tf_imag;

  Blocks = 1 << (p - 1);
  Points = 2;
  for (i = 0; i < p; i++)
  {
    Points2 = Points >> 1;
    BaseT = 0;
    for (j = 0; j < Blocks; j++)
    {
      BaseB = BaseT + Points2;
      for (k = 0; k < Points2; k++)
      {
        top_real = f_real[BaseT + k];
        top_imag = f_imag[BaseT + k];
        TWIDDLE(Points, k, tf_real, tf_imag, INVERSE_FFT);
        bot_real =
          f_real[BaseB + k] * tf_real - f_imag[BaseB + k] * tf_imag;
        bot_imag =
          f_real[BaseB + k] * tf_imag + f_imag[BaseB + k] * tf_real;
        f_real[BaseT + k] = top_real + bot_real;
        f_imag[BaseT + k] = top_imag + bot_imag;
        f_real[BaseB + k] = top_real - bot_real;
        f_imag[BaseB + k] = top_imag - bot_imag;
      }
      BaseT = BaseT + Points;
    }
    Blocks = Blocks >> 1;
    Points = Points << 1;
  }

}

static void fft_density_convolve(double *y, 
                                 double *kords, 
                                 int n, 
                                 AFFY_ERROR *err)
{
  int     i;
  /* ugly hack to stop rounding problems */
  int     nlog2 = (int)(log((double)n) / log(2.0) + 0.5);
  double *y_imag, *kords_imag, *conv_real, *conv_imag;
  int    *mempool;

  mempool = h_malloc(sizeof(int));
  if (mempool == NULL)
    AFFY_HANDLE_ERROR_VOID("malloc failed", AFFY_ERROR_OUTOFMEM, err);

  y_imag = (double *)h_subcalloc(mempool, n, sizeof(double));
  if (y_imag == NULL)
  {
    h_free(mempool);
    AFFY_HANDLE_ERROR_VOID("calloc failed", AFFY_ERROR_OUTOFMEM, err);
  }

  kords_imag = (double *)h_subcalloc(mempool, n, sizeof(double));
  if (kords_imag == NULL)
  {
    h_free(mempool);
    AFFY_HANDLE_ERROR_VOID("calloc failed", AFFY_ERROR_OUTOFMEM, err);
  }

  conv_real = y_imag;  
  conv_imag = kords_imag;

  fft_dif(y, y_imag, nlog2);
  fft_dif(kords, kords_imag, nlog2);

  /*
     We don't need y_imag and kords_imag after this, so
     reuse as conv_real/imag.
   */
  for (i = 0; i < n; i++)
  {
    double yi = y_imag[i];
    double ki = kords_imag[i];

    conv_real[i] = y[i] * kords[i] + yi * ki;
    conv_imag[i] = y[i] * (-1 * ki) + yi * kords[i];
  }

  fft_ditI(conv_real, conv_imag, nlog2);

  for (i = 0; i < n; i++)
  {
    kords[i] = conv_real[i];
  }

  h_free(mempool);
}

/**************************************************************
 **
 ** static void kernelize(double *data, int n, double bw, int kernel)
 **
 ** double *data - data to kernelize
 ** int n - length of data.
 ** double bw - bandwidth for Kernel
 ** int kernel - an integer specifying which kernel to use 
 **              1 is gaussian, 2 is Epanechnikov,
 **              3 is ...........
 **
 **
 ***************************************************************/
static void kernelize(double *data, int n, double bw, int kernel)
{
  double a = 0.0;
  int    i;

  if (kernel == 1)
  {
    /* Gaussian Kernel */
  }
  else if (kernel == 2)
  {
    /* Epanechnikov Kernel */
    a = bw * sqrt(5.0);

    for (i = 0; i < n; i++)
    {
      if (fabs(data[i]) < a)
      {
        data[i] =
          3.0 / (4.0 * a) * (1.0 -
                             (fabs(data[i]) / a) * (fabs(data[i]) / a));
      }
      else
      {
        data[i] = 0.0;
      }
    }
  }
}

/*****************************************************************
 **
 ** static double compute_sd(double *x, int length)
 **
 ** double *x - data vector
 ** int length - length of x
 **
 ** compute the standard deviation of a data vector
 **
 *****************************************************************/
static double compute_sd(double *x, int length)
{
  int    i;
  double sum = 0.0, sum2 = 0.0;

  for (i = 0; i < length; i++)
    sum += x[i];

  sum = sum / (double)length;

  for (i = 0; i < length; i++)
    sum2 += (x[i] - sum) * (x[i] - sum);

  return (sqrt(sum2 / (double)(length - 1)));
}

/*****************************************************************
 **
 ** static double bandwidth(double *x, int length, double iqr)
 **
 ** double *x - data vector
 ** int length - length of x
 ** double iqr - IQR of *x
 **
 ** compute the kernel bandwidth
 **
 *****************************************************************/
static double bandwidth(double *x, int length, double iqr)
{
  double hi;
  double lo;

  hi = compute_sd(x, length);

  if (hi > iqr)
    lo = iqr / 1.34;
  else
    lo = hi;

  if (lo == 0)
  {
    if (hi != 0)
    {
      lo = hi;
    }
    else if (fabs(x[1]) != 0)
    {
      lo = fabs(x[1]);
    }
    else
    {
      lo = 1.0;
    }
  }

  return (0.9 * lo * pow((double)length, -0.2));
}

/******************************************************************
 **
 ** linearly interpolate v given x and y.
 **
 **********************************************************************/
static double linear_interpolation(double v, double *x, double *y, int n)
{
  int i, j, ij;

  i = 0;
  j = n - 1;

  if (v < x[i])
    return (y[0]);
  if (v > x[j])
    return (y[n - 1]);

  /* find the correct interval by bisection */
  while (i < j - 1)
  {                             /* x[i] <= v <= x[j] */
    ij = (i + j) / 2;           /* i+1 <= ij <= j-1 */

    if (v < x[ij])
      j = ij;
    else
      i = ij;
    /* still i < j */
  }
  /* provably have i == j-1 */

  /* interpolation */
  if (v == x[j])
    return (y[j]);
  if (v == x[i])
    return (y[i]);
  /* impossible: if(x[j] == x[i]) return y[i]; */

  return (y[i] + (y[j] - y[i]) * ((v - x[i]) / (x[j] - x[i])));
}
