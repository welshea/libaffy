
/**************************************************************************
 *
 * Filename:  mas5_wilcox.h
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
 *
 **************************************************************************/

#include "affy.h"

struct affy_wilcox
{
  double r;
  double abs_r;
  double rank;
};

double affy_mas5_calculate_wilcox_pvalue(struct affy_wilcox *rset, int n);
double affy_mas5_calculate_call_pvalue(double *values, 
                                       int n, 
                                       double tau,
                                       AFFY_ERROR *err);
