
/**************************************************************************
 *
 * Filename:  pnorm.c
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
 * 04/08/05: Imported/repaired from old libaffy (AMH)
 * 10/28/08: Refactor trunc() into its own source file (AMH)
 *
 **************************************************************************/

#include <affy.h>
#include <utils.h>

/*
 * Since log1p() is not available in C89, an implementation of it
 * has been copied from GSL, the GNU Scientific Library.  Full
 * license/attribution follows.
 */

/* sys/log1p.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA.
 */
static double _gsl_log1p (const double x)
{
  volatile double y, z;
  y = 1 + x;
  z = y - 1;
  return log(y) - (z-x)/y ;  /* cancels errors with IEEE arithmetic */
}

/* End inclusion from the GNU Scientific Library */

/* 
   This file is part of RMAExpress.

    RMAExpress is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    RMAExpress is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with RMAExpress; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA 
*/

/*********************************************************************
 **
 ** file: pnorm.c
 **
 ** aim: distribution function of normal distribution
 **
 ** note that method based upon noted reference
 **
 **  REFERENCE
 **
 **    Cody, W. D. (1993).
 **    ALGORITHM 715: SPECFUN - A Portable FORTRAN Package of
 **    Special Function Routines and Test Drivers".
 **    ACM Transactions on Mathematical Software. 19, 22-32.
 **
 **
 *********************************************************************/

/* S    REAL FUNCTION ANORM(ARG) */

/* D    DOUBLE PRECISION FUNCTION ANORM(ARG) */

/* ------------------------------------------------------------------ */

/* This function evaluates the normal distribution function: */

/*                              / x */

/*                     1       |       -t*t/2 */

/*          P(x) = ----------- |      e       dt */

/*                 sqrt(2 pi)  | */

/*                             /-oo */

/*   The main computation evaluates near-minimax approximations */

/*   derived from those in "Rational Chebyshev approximations for */

/*   the error function" by W. J. Cody, Math. Comp., 1969, 631-637. */

/*   This transportable program uses rational functions that */

/*   theoretically approximate the normal distribution function to */

/*   at least 18 significant decimal digits.  The accuracy achieved */

/*   depends on the arithmetic system, the compiler, the intrinsic */

/*   functions, and proper selection of the machine-dependent */

/*   constants. */

/* ******************************************************************* */

/* ******************************************************************* */

/* Explanation of machine-dependent constants.  Let */

/*   XMIN  = the smallest positive floating-point number. */

/* Then the following machine-dependent constants must be declared */

/*   in DATA statements.  IEEE values are provided as a default. */

/*   EPS   = argument below which anorm(x) may be represented by */

/*           0.5  and above which  x*x  will not underflow. */

/*           A conservative value is the largest machine number X */

/*           such that   1.0 + X = 1.0   to machine precision. */

/*   XLOW  = the most negative argument for which ANORM does not */

/*           vanish.  This is the negative of the solution to */

/*                    W(x) * (1-1/x**2) = XMIN, */

/*           where W(x) = exp(-x*x/2)/[x*sqrt(2*pi)]. */

/*   XUPPR = positive argument beyond which anorm = 1.0.  A */

/*           conservative value is the solution to the equation */

/*                    exp(-x*x/2) = EPS, */

/*           i.e., XUPPR = sqrt[-2 ln(eps)]. */

/*   Approximate values for some important machines are: */

/*                          XMIN        EPS        XLOW    XUPPR */

/*  CDC 7600      (S.P.)  3.13E-294   7.11E-15   -36.641   8.072 */

/*  CRAY-1        (S.P.)  4.58E-246   7.11E-157 -106.521  26.816 */

/*  IEEE (IBM/XT, */

/*    SUN, etc.)  (S.P.)  1.18E-38    5.96E-8    -12.949   5.768 */

/*  IEEE (IBM/XT, */

/*    SUN, etc.)  (D.P.)  2.23D-308   1.11D-16   -37.519   8.572 */

/*  IBM 195       (D.P.)  5.40D-79    1.39D-17   -18.781   8.811 */

/*  VAX D-Format  (D.P.)  2.94D-39    1.39D-17   -13.055   8.811 */

/*  VAX G-Format  (D.P.)  5.56D-309   1.11D-16   -37.556   8.572 */

/* ******************************************************************* */

/* ******************************************************************* */

/* Error returns */

/*  The program returns  ANORM = 0     for  ARG .LE. XLOW. */

/* Intrinsic functions required are: */

/*     ABS, AINT, EXP */

/*  Author: W. J. Cody */

/*          Mathematics and Computer Science Division */

/*          Argonne National Laboratory */

/*          Argonne, IL 60439 */

/*  Latest modification: March 15, 1992 */

/* ------------------------------------------------------------------ */

void affy_pnorm_both(double x, double *cum, double *ccum, int i_tail,
                     int log_p)
{
  /* Initialized data */

  static double a[5] = { 2.2352520354606839287, 161.02823106855587881,
    1067.6894854603709582, 18154.981253343561249,
    .065682337918207449113
  };
  static double b[4] = { 47.20258190468824187, 976.09855173777669322,
    10260.932208618978205, 45507.789335026729956
  };
  static double c[9] = { .39894151208813466764, 8.8831497943883759412,
    93.506656132177855979, 597.27027639480026226, 2494.5375852903726711,
    6848.1904505362823326, 11602.651437647350124, 9842.7148383839780218,
    1.0765576773720192317e-8
  };
  static double d[8] = { 22.266688044328115691, 235.38790178262499861,
    1519.377599407554805, 6485.558298266760755, 18615.571640885098091,
    34900.952721145977266, 38912.003286093271411, 19685.429676859990727
  };
  static double p[6] = { .21589853405795699, .1274011611602473639,
    .022235277870649807, .001421619193227893466, 2.9112874951168792e-5,
    .02307344176494017303
  };
  static double q[5] = { 1.28426009614491121, .468238212480865118,
    .0659881378689285515, .00378239633202758244, 7.29751555083966205e-5
  };
  static double SIXTEN = 16.;
  static double M_1_SQRT_2PI = .39894228040143267794;
  static double M_LN_SQRT_2PI = 0.918938533204672741780329736406;
  static double thrsh = 0.67448975;
  static double root32 = 5.656854248;
  static double eps = 1.11e-16;

  /* Local variables */
  static int i, lower, upper;
  static double y, del, xsq, xden, xnum, temp;

  lower = i_tail != 1;
  upper = i_tail != 0;

  y = fabs(x);

  if (y <= thrsh)
  {
    if (y > eps)
    {
      xsq = x * x;
      xnum = a[4] * xsq;
      xden = xsq;
      for (i = 0; i < 3; ++i)
      {
        xnum = (xnum + a[i]) * xsq;
        xden = (xden + b[i]) * xsq;
      }
    }
    else
      xnum = xden = 0.0;

    temp = x * (xnum + a[3]) / (xden + b[3]);
    if (lower)
      *cum = 0.5 + temp;
    if (upper)
      *ccum = 0.5 - temp;
    if (log_p)
    {
      if (lower)
        *cum = log(*cum);
      if (upper)
        *ccum = log(*ccum);
    }
    /* ------------------------------------------------------------------ */
    /*  Evaluate  anorm  for threshhold <= |X| <= sqrt(32) */
    /* ------------------------------------------------------------------ */
  }
  else if (y <= root32)
  {
    xnum = c[8] * y;
    xden = y;
    for (i = 0; i < 7; ++i)
    {
      xnum = (xnum + c[i]) * y;
      xden = (xden + d[i]) * y;
    }
    temp = (xnum + c[7]) / (xden + d[7]);

#define do_del(X)							\
	xsq = affy_trunc(X * SIXTEN) / SIXTEN;				\
	del = (X - xsq) * (X + xsq);					\
	if(log_p) {							\
	    *cum = (-xsq * xsq * 0.5) + (-del * 0.5) + log(temp);	\
	    if((lower && x > 0.) || (upper && x <= 0.))			\
		  *ccum = _gsl_log1p(-exp(-xsq * xsq * 0.5) * 		\
				exp(-del * 0.5) * temp);		\
	}								\
	else {								\
	    *cum = exp(-xsq * xsq * 0.5) * exp(-del * 0.5) * temp;	\
	    *ccum = 1.0 - *cum;						\
	}

#define swap_tail						\
	if (x > 0.) {/* swap  ccum <--> cum */			\
	    temp = *cum; if(lower) *cum = *ccum; *ccum = temp;	\
	}

    do_del(y);
    swap_tail;
    /* ------------------------------------------------------------------ */
    /*  Evaluate  anorm  for |X| > sqrt(32) */
    /* ------------------------------------------------------------------ */
  }
  else
   if ((-37.5193 < x) || (x < 8.2924))
  {                             /* originally had y < 50 */

    /* Evaluate pnorm for x in (-37.5, -5.657) union (5.657, 8.29) */

    xsq = 1.0 / (x * x);
    xnum = p[5] * xsq;
    xden = xsq;
    for (i = 0; i < 4; ++i)
    {
      xnum = (xnum + p[i]) * xsq;
      xden = (xden + q[i]) * xsq;
    }
    temp = xsq * (xnum + p[4]) / (xden + q[4]);
    temp = (M_1_SQRT_2PI - temp) / y;

    do_del(x);
    swap_tail;
  }
  else
  {                             /* x < -37.5193  OR  8.2924 < x */
    if (log_p)
    {                           /* be better than to just return log(0) or log(1) */
      xsq = x * x;
      if (xsq * eps < 1.)
        del = (1. - (1. - 5. / (xsq + 6.)) / (xsq + 4.)) / (xsq + 2.);
      else
        del = 0.;
      *cum = -.5 * xsq - M_LN_SQRT_2PI - log(y) + _gsl_log1p(del);
      *ccum = -0.;              /*log(1) */
      swap_tail;

    }
    else
    {
      if (x > 0)
      {
        *cum = 1.;
        *ccum = 0.;
      }
      else
      {
        *cum = 0.;
        *ccum = 1.;
      }
    }
  }

}

double affy_pnorm5(double x, double mu, double sigma, int lower_tail,
                   int log_p)
{
  double p, cp;

  x = (x - mu) / sigma;

  affy_pnorm_both(x, &p, &cp, (lower_tail ? 0 : 1), log_p);

  return (lower_tail ? p : cp);
}
