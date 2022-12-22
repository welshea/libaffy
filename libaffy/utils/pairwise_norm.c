
/**************************************************************************
 *
 * Filename:  pairwise_norm.c
 *
 * Purpose:   Pairwise normalization method.
 *
 * Creation:  10/07/10
 *
 * Author:    Eric Welsh
 *
 * Copyright: Copyright (C) 2010, Moffitt Cancer Center
 *            All rights reserved.
 *
 * Update History
 * --------------
 * 10/07/10: Creation (EAW)
 * 01/10/13: Prevent iterative pruning from removing all points when N
 *           is small (EAW)
 * 01/10/13: Fixed rare case of reading beyond equation windows bounds,
 *           particularly when N points is small
 * 01/10/13: renamed linear normalization to global scaling, since that's
 *           actually what it is (true intensity-dependent linear
 *           linear normalization is not implemented yet)
 * 01/18/13: Pseudo-density weight exponent is now user definable.
 *           4 works well for microarray, 0 (unweighted) or at most 1 works
 *           best for proteomics.
 * 09/18/13: added iron_fit_both_x_y flag: better norm, may distrub rank order
 * 09/18/13: added iron_fit_window_frac (EAW)
 * 10/04/17: added iron_condense_training flag (EAW)
 * 06/01/18: added support for probeset exclusions during IRON training (EAW)
 * 10/11/19: print GlobalFitLine stats to stderr, similar to GlobalScale (EAW)
 *           stats are for the fit line, rather than the scaling adjustments
 *           (opposite direction/interpretation as GlobalScale)
 * 05/25/21: better handling of empty and near-empty samples (EAW)
 * 12/22/22: define M_PI if not already defined; for Ubuntu/Debian (EAW)
 *
 **************************************************************************/

#include "affy.h"

/* #include "pairwise_norm.h" */

#define MIN_SIGNAL       1E-5
#define DO_FLOOR         1

#if !defined(M_PI)
#define M_PI 3.14159265358979323846
#endif

/* (Unweighted): cisplatin breast subset works best with NO second pass,
 * rank fraction = 0.01, and window fraction = 0.05
 */

/* Window width fraction:
 *   0.05 too bumpy in UTSouthwestern Illumina 5058818012_E vs. 5067386018_F
 *   0.10 still a little bumpy, but much much smoother
 */

/* discard outliers from 1st pass and retrain with more lax rank order fit */
/* NOT RECOMMENDED, for experimental purposes only! */
#define SECOND_PASS_TRAIN 0


/*
 * flags for development and debug printing
 */
#define DEBUG_PRINT              1
#define DEBUG_FILE               0
#define DEBUG_COLOR_IRANK        0
#define DEBUG_DIE_EARLY          0
#define DEBUG_FIXED_RANK         0

/* should I replace with AFFY_POINT? */
struct xy_pair
{
  double x, y;
};


struct signal_pair
{
  int    index;
  double sig1;
  double sig2;
  int    rank1;
  int    rank2;
  int    rank_diff;
  char   initial_set_flag;
  char   irank_flag;
 
#if DEBUG_COLOR_IRANK
  double irank_frac_0;
  double irank_frac;
  double norm_err_scaled;
#endif
 
  double log_xy;
  double log_adjust;
  double fit_log_adjust;
  double norm_err;
  
  double weight;
  int n_windows;
};

struct eqn_window
{
  double slope;
  double offset;
  double start;
  double end;
};


static int compare_sort_xy_pair_by_x(const void *vptr1, const void *vptr2)
{
  struct xy_pair *pptr1 = (struct xy_pair *) vptr1;
  struct xy_pair *pptr2 = (struct xy_pair *) vptr2;

  assert(pptr1 != NULL);
  assert(pptr2 != NULL);

  if (pptr1->x < pptr2->x)
    return (-1);
  if (pptr1->x > pptr2->x)
    return (1);

  if (pptr1->y < pptr2->y)
    return (-1);
  if (pptr1->y > pptr2->y)
    return (1);

  return (0);
}


static int compare_sort_xy_pair_by_y(const void *vptr1, const void *vptr2)
{
  struct xy_pair *pptr1 = (struct xy_pair *) vptr1;
  struct xy_pair *pptr2 = (struct xy_pair *) vptr2;

  assert(pptr1 != NULL);
  assert(pptr2 != NULL);

  if (pptr1->y < pptr2->y)
    return (-1);
  if (pptr1->y > pptr2->y)
    return (1);
 
  if (pptr1->x < pptr2->x)
    return (-1);
  if (pptr1->x > pptr2->x)
    return (1);

  return (0);
}


static int compare_sort_xy_pair_by_x_plus_y(const void *vptr1,
                                            const void *vptr2)
{
  struct xy_pair *pptr1 = (struct xy_pair *) vptr1;
  struct xy_pair *pptr2 = (struct xy_pair *) vptr2;
  double val1, val2;

  assert(pptr1 != NULL);
  assert(pptr2 != NULL);

  val1 = pptr1->x + pptr1->y;
  val2 = pptr2->x + pptr2->y;

  if (val1 < val2)
    return (-1);
  if (val1 > val2)
    return (1);

  if (pptr1->y < pptr2->y)
    return (-1);
  if (pptr1->y > pptr2->y)
    return (1);
 
  if (pptr1->x < pptr2->x)
    return (-1);
  if (pptr1->x > pptr2->x)
    return (1);

  return (0);
}


static int compare_sort_sig1(const void *vptr1, const void *vptr2)
{
  struct signal_pair *pptr1 = *((struct signal_pair **) vptr1);
  struct signal_pair *pptr2 = *((struct signal_pair **) vptr2);

  assert(pptr1 != NULL);
  assert(pptr2 != NULL);
   
  if (pptr1->sig1 < pptr2->sig1)
    return (-1);
  if (pptr1->sig1 > pptr2->sig1)
    return (1);

  if (pptr1->sig2 < pptr2->sig2)
    return (-1);
  if (pptr1->sig2 > pptr2->sig2)
    return (1);

  if (pptr1->index < pptr2->index)
    return (-1);
  if (pptr1->index > pptr2->index)
    return (1);
   
  return (0);
}


static int compare_sort_sig2(const void *vptr1, const void *vptr2)
{
  struct signal_pair *pptr1 = *((struct signal_pair **) vptr1);
  struct signal_pair *pptr2 = *((struct signal_pair **) vptr2);

  assert(pptr1 != NULL);
  assert(pptr2 != NULL);
   
  if (pptr1->sig2 < pptr2->sig2)
    return (-1);
  if (pptr1->sig2 > pptr2->sig2)
    return (1);

  if (pptr1->sig1 < pptr2->sig1)
    return (-1);
  if (pptr1->sig1 > pptr2->sig1)
    return (1);

  if (pptr1->index < pptr2->index)
    return (-1);
  if (pptr1->index > pptr2->index)
    return (1);

  return (0);
}


static int compare_sort_log_xy(const void *vptr1, const void *vptr2)
{
  struct signal_pair *pptr1 = *((struct signal_pair **) vptr1);
  struct signal_pair *pptr2 = *((struct signal_pair **) vptr2);

  assert(pptr1 != NULL);
  assert(pptr2 != NULL);
   
  if (pptr1->log_xy < pptr2->log_xy)
    return (-1);
  if (pptr1->log_xy > pptr2->log_xy)
    return (1);

  if (pptr1->sig1 < pptr2->sig1)
    return (-1);
  if (pptr1->sig1 > pptr2->sig1)
    return (1);

  if (pptr1->sig2 < pptr2->sig2)
    return (-1);
  if (pptr1->sig2 > pptr2->sig2)
    return (1);
   
  if (pptr1->index < pptr2->index)
    return (-1);
  if (pptr1->index > pptr2->index)
    return (1);

  return (0);
}


static int fill_geometric_eqn_windows(struct eqn_window *eqn_windows,
                                      struct signal_pair **filt_ptrs,
                                      int num_pairs, double window_frac,
                                      int pass, double weight_exponent)
{
  struct signal_pair *pair_ptr;
  int                 n, i;
  int                 w = (int)(window_frac * num_pairs + 0.5);
  int                 wsmall = (int)(0.01 * num_pairs + 0.5);
  double              x, y;
  double              ss_xx = 0.0, ss_xy = 0.0;
  double              x_sum = 0.0, y_sum = 0.0;
  double              x_avg, y_avg;
  double              weight, weight_sum = 0.0;
  double              min_weight = 9E99, max_weight = -9E99;
   
  assert(filt_ptrs   != NULL);
  assert(eqn_windows != NULL);

  if (w < 100)
    w = 100;
  if (wsmall < 10)
    wsmall = 10;
  
  if (w > num_pairs)
    w = num_pairs;
  if (wsmall > num_pairs)
    wsmall = num_pairs;

  qsort(filt_ptrs,
        num_pairs,
        sizeof(struct signal_pair *),
        compare_sort_log_xy);


  /* weights */

  /* initialize weights */
  for (i = 0; i < num_pairs; i++)
  {
    filt_ptrs[i]->weight = 0;
    filt_ptrs[i]->n_windows = 0;
  }

  /* slide bin windows */
  for (n = 0; n <= num_pairs - wsmall; n++)
  {
    x_avg = 0;
    for (i = n; i < n + wsmall; i++)
      x_avg += filt_ptrs[i]->log_xy;
    x_avg /= wsmall;
    
    weight = 0;
    for (i = n; i < n + wsmall; i++)
    {
      x = filt_ptrs[i]->log_xy - x_avg;
      weight += x*x;
    }
    weight = sqrt(weight / wsmall);

    for (i = n; i < n + wsmall; i++)
    {
      filt_ptrs[i]->weight += weight;
      filt_ptrs[i]->n_windows++;
    }
  }
  for (i = 0; i < num_pairs; i++)
  {
    filt_ptrs[i]->weight /= filt_ptrs[i]->n_windows;
    
    if (filt_ptrs[i]->weight >= 1E-5 && filt_ptrs[i]->weight < min_weight)
      min_weight = filt_ptrs[i]->weight;
    
    if (filt_ptrs[i]->weight > max_weight)
      max_weight = filt_ptrs[i]->weight;
  }
  
  if (min_weight == 9E99)
    min_weight = max_weight;

#if DEBUG_PRINT
  fprintf(stderr, "Weights:\t%f\t%f\t%f\n",
          min_weight, max_weight, max_weight / min_weight);
#endif

  for (i = 0; i < num_pairs; i++)
  {
    if (filt_ptrs[i]->weight < 1E-5)
      filt_ptrs[i]->weight = min_weight;

    /*
     * w^4 = sigma^4 = variance^2; w^8 = (variance^2)^2
     */
     filt_ptrs[i]->weight = pow(filt_ptrs[i]->weight / max_weight, weight_exponent);
  }


  /* local windowed fits */
  n = 0;
  weight_sum = 0;

  for (i = 0; i < w; i++)
  {
    pair_ptr = filt_ptrs[i];
    weight = pair_ptr->weight;
    weight_sum += weight;
       
    x      = pair_ptr->log_xy;
    x_sum += weight * x;
    ss_xx += weight * x*x;

    y      = pair_ptr->log_adjust;
    y_sum += weight * y;
    ss_xy += weight * x*y;
  }
   
  while (1)
  {
    double temp;
 
    x_avg = x_sum / weight_sum;
    y_avg = y_sum / weight_sum;
   
    temp = ss_xx - weight_sum * x_avg * x_avg;
    if (temp)
      eqn_windows[n].slope = (ss_xy - weight_sum * x_avg * y_avg) / temp;
    else
      eqn_windows[n].slope = 0.0;
    eqn_windows[n].offset = y_avg - eqn_windows[n].slope * x_avg;
    eqn_windows[n].start  = filt_ptrs[n]->log_xy;
    eqn_windows[n].end    = filt_ptrs[n+w-1]->log_xy;

    if (n >= num_pairs - w)
    {
      n++;
      break;
    }

    /* add next point to sums */
    pair_ptr = filt_ptrs[n+w];
    weight = pair_ptr->weight;
    weight_sum += weight;
       
    x      = pair_ptr->log_xy;
    x_sum += weight * x;
    ss_xx += weight * x*x;

    y      = pair_ptr->log_adjust;
    y_sum += weight * y;
    ss_xy += weight * x*y;


    /* remove first point from sums */
    pair_ptr = filt_ptrs[n];
    weight = pair_ptr->weight;
    weight_sum -= weight;
    
    x      = pair_ptr->log_xy;
    x_sum -= weight * x;
    ss_xx -= weight * x*x;

    y      = pair_ptr->log_adjust;
    y_sum -= weight * y;
    ss_xy -= weight * x*y;
       
    n++;
  }

  return (n);
}


static void smooth_geometric_fits(struct eqn_window *eqn_windows,
                                  int num_eqn_windows,
                                  struct signal_pair **filt_ptrs,
                                  int num_pairs)
{
  struct signal_pair *pair_ptr;
  double              x, avg_adjust, sum_slope = 0, sum_offset = 0;
  int                 i, j, min_eqn_idx = 0, end_eqn_idx = 0;
  int                 old_min_eqn_idx = 0, old_end_eqn_idx = 0;

  assert(filt_ptrs   != NULL);
  assert(eqn_windows != NULL);

  for (i = 0; i < num_pairs; i++)
  {
    pair_ptr = filt_ptrs[i];
    x        = pair_ptr->log_xy;

#if 0
    if (isinf(x))
      printf("NAN %d %f %f\n", i, pair_ptr->sig1, pair_ptr->sig2);
#endif
   
    old_min_eqn_idx = min_eqn_idx;
    old_end_eqn_idx = end_eqn_idx;
       
    /* find first possibly overlapping eqn window */
    while (min_eqn_idx < num_eqn_windows && eqn_windows[min_eqn_idx].end < x)
      min_eqn_idx++;

    if (end_eqn_idx < min_eqn_idx)
      end_eqn_idx = min_eqn_idx;

    while (end_eqn_idx < num_eqn_windows &&
           x >= eqn_windows[end_eqn_idx].start &&
           x <= eqn_windows[end_eqn_idx].end)
      end_eqn_idx++;

    /* sum slope and offset */
    if (i == 0)
    {
      for (j = min_eqn_idx; j < end_eqn_idx; j++)
      {
        sum_slope += eqn_windows[j].slope;
        sum_offset += eqn_windows[j].offset;
      }
    }
    /* subtract out old window, add in new window */
    else
    {
      for (j = old_min_eqn_idx; j < min_eqn_idx; j++)
      {
        sum_slope -= eqn_windows[j].slope;
        sum_offset -= eqn_windows[j].offset;
      }

      for (j = old_end_eqn_idx; j < end_eqn_idx; j++)
      {
        sum_slope += eqn_windows[j].slope;
        sum_offset += eqn_windows[j].offset;
      }
    }

    /* average the fit adjust value over all overlapping eqns */
    avg_adjust = (sum_slope * x + sum_offset) / (end_eqn_idx - min_eqn_idx);

    pair_ptr->fit_log_adjust = avg_adjust;
    pair_ptr->norm_err       = avg_adjust - pair_ptr->log_adjust;

#if 0
    if (isnan(pair_ptr->fit_log_adjust))
      printf("NAN_LOGADJUST %f\n", x);
    if (isnan(avg_adjust))
      printf("NAN_AVGADJUST %f\n", x);
#endif
  }
}


/* WARNING -- can go haywire on rapid changes, DO NOT USE */
static double lagrange_interp(double x,
                              struct xy_pair *xy_pairs,
                              int last_idx,
                              int idx)
{
  struct xy_pair *xy;
  double          a, sum = 0;

  assert(xy_pairs != NULL);
   
  /* cubic lagrange */
  if (idx >= 2 && idx < last_idx)
  {
    xy = xy_pairs + idx - 2;
   
    a    = (x - xy[1].x) / (xy[0].x - xy[1].x);
    a   *= (x - xy[2].x) / (xy[0].x - xy[2].x);
    a   *= (x - xy[3].x) / (xy[0].x - xy[3].x);
    sum += a * xy[0].y;

    a    = (x - xy[0].x) / (xy[1].x - xy[0].x);
    a   *= (x - xy[2].x) / (xy[1].x - xy[2].x);
    a   *= (x - xy[3].x) / (xy[1].x - xy[3].x);
    sum += a * xy[1].y;

    a    = (x - xy[0].x) / (xy[2].x - xy[0].x);
    a   *= (x - xy[1].x) / (xy[2].x - xy[1].x);
    a   *= (x - xy[3].x) / (xy[2].x - xy[3].x);
    sum += a * xy[2].y;

    a    = (x - xy[0].x) / (xy[3].x - xy[0].x);
    a   *= (x - xy[1].x) / (xy[3].x - xy[1].x);
    a   *= (x - xy[2].x) / (xy[3].x - xy[2].x);
    sum += a * xy[3].y;

    return (sum);
  }
  /* linear */
  else if (idx >= 1 && xy_pairs[idx].x != xy_pairs[idx-1].x)
  {
    a = (xy_pairs[idx].x - x) / (xy_pairs[idx].x - xy_pairs[idx-1].x);

    return (a * xy_pairs[idx-1].y + (1.0 - a) * xy_pairs[idx].y);
  }

  /* use the value at idx */
  return (xy_pairs[idx].y);
}


static double linear_interp(double x,
                            struct xy_pair *xy_pairs,
                            int last_idx,
                            int idx)
{
  double          a;

  assert(xy_pairs != NULL);

  if (idx >= 1 && xy_pairs[idx].x != xy_pairs[idx-1].x)
  {
    a = (xy_pairs[idx].x - x) / (xy_pairs[idx].x - xy_pairs[idx-1].x);

    return (a * xy_pairs[idx-1].y + (1.0 - a) * xy_pairs[idx].y);
  }

  /* use the value at idx */
  return (xy_pairs[idx].y);
}


static void interpolate_final_scales(struct signal_pair **pair_ptrs,
                                     int num_pairs,
                                     struct signal_pair **filt_ptrs_train,
                                     int num_pairs_train,
                                     int fit_both_x_y_flag,
                                     AFFY_ERROR *err)
{
  struct xy_pair    *xy_pairs = NULL;
  struct xy_pair    *xy_pairs2 = NULL;
  struct signal_pair *pair_ptr;
  double              x, y, old_x, old_y, sum;
  int                 num_pairs_train2 = 0;
  int                 num_pairs_train3 = 0;
  int                 i, num;
  int                 min_idx, last_idx;

  assert(pair_ptrs       != NULL);
  assert(filt_ptrs_train != NULL);
 

  xy_pairs = (struct xy_pair *) calloc(num_pairs_train,
                                       sizeof(struct xy_pair));
  if (xy_pairs == NULL)
    AFFY_HANDLE_ERROR_VOID("malloc failed", AFFY_ERROR_OUTOFMEM, err);

  qsort(pair_ptrs,
        num_pairs,
        sizeof(struct signal_pair *),
        compare_sort_sig2);
  qsort(filt_ptrs_train,
        num_pairs_train,
        sizeof(struct signal_pair *),
        compare_sort_sig2);

  /* condense overlapping points */
  /* project onto average best fit line */
  old_x = -9E99;
  old_y = -9E99;
  for (i = 0, num_pairs_train2 = 0; i < num_pairs_train; i++)
  {
    pair_ptr = filt_ptrs_train[i];
    
    /* project x/y onto best fit avg line using norm_err */
    x = log(pair_ptr->sig1) + 0.5 * pair_ptr->norm_err;
    y = log(pair_ptr->sig2) - 0.5 * pair_ptr->norm_err;
    
    /* is new point different from last point? */
    /* take round off error into account */
    if (fabs(x - old_x) > 1E-14 || fabs(y - old_y) > 1E-14)
    {
      xy_pairs[num_pairs_train2].x   = x;
      xy_pairs[num_pairs_train2++].y = y;
    }
    
    old_x = x;
    old_y = y;
  }

  qsort(xy_pairs,
        num_pairs_train2,
        sizeof(struct xy_pair),
        compare_sort_xy_pair_by_y);

  xy_pairs2 = (struct xy_pair *) calloc(num_pairs_train,
                                        sizeof(struct xy_pair));
  if (xy_pairs2 == NULL)
    AFFY_HANDLE_ERROR_VOID("malloc failed", AFFY_ERROR_OUTOFMEM, err);


  /* store as [log(y), fit_log_adjust] pairs */
  old_x = -9E99;
  sum = 0;
  num = 0;
  for (i = 0, num_pairs_train3 = -1; i < num_pairs_train2; i++)
  {
    x = xy_pairs[i].y;
    
#if 0
    if (old_x == x)
    {
      fprintf(stderr, "SAMEXY     %f %f %f     %f %f %f\n",
              xy_pairs[i-1].x,
              xy_pairs[i-1].y,
              xy_pairs[i-1].x - xy_pairs[i-1].y,
              xy_pairs[i].x,
              xy_pairs[i].y,
              xy_pairs[i].x - xy_pairs[i].y);
      printf("SAMEXY     %f %f %f     %f %f %f\n",
              xy_pairs[i-1].x,
              xy_pairs[i-1].y,
              xy_pairs[i-1].x - xy_pairs[i-1].y,
              xy_pairs[i].x,
              xy_pairs[i].y,
              xy_pairs[i].x - xy_pairs[i].y);
    }
#endif
    
    if (x != old_x)
    {
      num_pairs_train3++;
      sum = 0;
      num = 0;
    }
    old_x = x;

    sum += xy_pairs[i].x - xy_pairs[i].y;
    num++;

    xy_pairs2[num_pairs_train3].x = x;
    
    /* average scales for identical points */
    if (num > 1)
      xy_pairs2[num_pairs_train3].y = sum / num;
    else
      xy_pairs2[num_pairs_train3].y = sum;
  }

  num_pairs_train3++;
  last_idx = num_pairs_train3 - 1;
  min_idx = 0;
  old_x = -9E99;


#if DEBUG_PRINT
  fprintf(stderr, "TrainingY\t%d\t%d\t%d\n",
          num_pairs_train, num_pairs_train2, num_pairs_train3);
#endif

  for (i = 0; i < num_pairs; i++)
  {
    pair_ptr = pair_ptrs[i];
    x        = log(pair_ptr->sig2);

    /* same X as last time, use the previous fit */
    if (x == old_x && i)
    {
      pair_ptr->fit_log_adjust = sum;
      continue;
    }

    old_x = x;

    /* before first training point */
    /* use average of end 10 adjustments */
    if (x < xy_pairs2[0].x)
    {
      int j;
    
      sum = 0.0;
      for (j = 0; j < 10 && j < last_idx; j++)
        sum += xy_pairs2[j].y;
      if (j)
        sum /= j;
      pair_ptr->fit_log_adjust = sum;

      continue;
    }
        
    /* after last training point */
    /* use average of end 10 adjustments */
    if (x > xy_pairs2[last_idx].x)
    {
      int j;

      sum = 0.0;
      for (j = 0; j < 10 && j < last_idx; j++)
        sum += xy_pairs2[last_idx - j].y;
      if (j)
        sum /= j;
      pair_ptr->fit_log_adjust = sum;

      continue;
    }
        
    /* skip to at or just after current point */
    while(xy_pairs2[min_idx].x < x)
      min_idx++;

    /* Use linear interpolation.
     * WARNING -- curve can have occasional too-rapid changes which
     * greatly confuse Lagrange interpolation and lead to big overshoots.
     * This tends to happen when the adjustment is very close to zero
     * and it oscillates a bit around zero.  Linear interpolation has
     * no issues, since it doesn't require smoothness.
     */
    sum = linear_interp(x, xy_pairs2, last_idx, min_idx);
    pair_ptr->fit_log_adjust = sum;

#if 0
    printf("DEBUG\t%f\t%d\t%f\t%f\t%f\n",
           x, min_idx, xy_pairs2[min_idx].x, xy_pairs2[min_idx].y, pair_ptr->fit_log_adjust);
#endif
  }


  /* Fit scale vs. X, fit scale vs. Y, average the scaling from both fits.
   *
   * Not usually recommended, since it will result in altered intensity rank
   * orders.
   *
   * However, for especially ill-behaved data, it can result in overall better
   * normalizations.
   */
  if (fit_both_x_y_flag)
  {
    qsort(pair_ptrs,
          num_pairs,
          sizeof(struct signal_pair *),
          compare_sort_sig1);
    qsort(filt_ptrs_train,
          num_pairs_train,
          sizeof(struct signal_pair *),
          compare_sort_sig1);

    /* condense overlapping points */
    /* project onto average best fit line */
    old_x = -9E99;
    old_y = -9E99;
    for (i = 0, num_pairs_train2 = 0; i < num_pairs_train; i++)
    {
      pair_ptr = filt_ptrs_train[i];
    
      /* project x/y onto best fit avg line using norm_err */
      x = log(pair_ptr->sig1) + 0.5 * pair_ptr->norm_err;
      y = log(pair_ptr->sig2) - 0.5 * pair_ptr->norm_err;
    
      /* is new point different from last point? */
      /* take round off error into account */
      if (fabs(x - old_x) > 1E-14 || fabs(y - old_y) > 1E-14)
      {
        xy_pairs[num_pairs_train2].x   = x;
        xy_pairs[num_pairs_train2++].y = y;
      }
    
      old_x = x;
      old_y = y;
    }

    qsort(xy_pairs,
          num_pairs_train2,
          sizeof(struct xy_pair),
          compare_sort_xy_pair_by_x);

    /* store as [log(x), fit_log_adjust] pairs */
    old_x = -9E99;
    sum = 0;
    num = 0;
    for (i = 0, num_pairs_train3 = -1; i < num_pairs_train2; i++)
    {
      x = xy_pairs[i].x;
    
#if 0
      if (old_x == x)
      {
        fprintf(stderr, "SAMEXY     %f %f %f     %f %f %f\n",
                xy_pairs[i-1].x,
                xy_pairs[i-1].y,
                xy_pairs[i-1].x - xy_pairs[i-1].y,
                xy_pairs[i].x,
                xy_pairs[i].y,
                xy_pairs[i].x - xy_pairs[i].y);
        printf("SAMEXY     %f %f %f     %f %f %f\n",
                xy_pairs[i-1].x,
                xy_pairs[i-1].y,
                xy_pairs[i-1].x - xy_pairs[i-1].y,
                xy_pairs[i].x,
                xy_pairs[i].y,
                xy_pairs[i].x - xy_pairs[i].y);
      }
#endif
    
      if (x != old_x)
      {
        num_pairs_train3++;
        sum = 0;
        num = 0;
      }
      old_x = x;

      sum += xy_pairs[i].x - xy_pairs[i].y;
      num++;

      xy_pairs2[num_pairs_train3].x = x;
    
      /* average scales for identical points */
      if (num > 1)
        xy_pairs2[num_pairs_train3].y = sum / num;
      else
        xy_pairs2[num_pairs_train3].y = sum;
    }

    num_pairs_train3++;
    last_idx = num_pairs_train3 - 1;
    min_idx = 0;
    old_x = -9E99;


#if DEBUG_PRINT
    fprintf(stderr, "TrainingX\t%d\t%d\t%d\n",
            num_pairs_train, num_pairs_train2, num_pairs_train3);
#endif

    for (i = 0; i < num_pairs; i++)
    {
      pair_ptr = pair_ptrs[i];
      x        = log(pair_ptr->sig1);

      /* same X as last time, use the previous fit */
      if (x == old_x && i)
      {
        pair_ptr->fit_log_adjust = 0.5 * (pair_ptr->fit_log_adjust + sum);
        continue;
      }

      old_x = x;

      /* before first training point */
      /* use average of end 10 adjustments */
      if (x < xy_pairs2[0].x)
      {
        int j;
    
        sum = 0.0;
        for (j = 0; j < 10 && j < last_idx; j++)
          sum += xy_pairs2[j].y;
        if (j)
          sum /= j;
        pair_ptr->fit_log_adjust = 0.5 * (pair_ptr->fit_log_adjust + sum);

        continue;
      }
        
      /* after last training point */
      /* use average of end 10 adjustments */
      if (x > xy_pairs2[last_idx].x)
      {
        int j;

        sum = 0.0;
        for (j = 0; j < 10 && j < last_idx; j++)
          sum += xy_pairs2[last_idx - j].y;
        if (j)
          sum /= j;
        pair_ptr->fit_log_adjust = 0.5 * (pair_ptr->fit_log_adjust + sum);

        continue;
      }

      /* skip to at or just after current point */
      while(xy_pairs2[min_idx].x < x)
        min_idx++;

      /* Use linear interpolation.
       * WARNING -- curve can have occasional too-rapid changes which
       * greatly confuse Lagrange interpolation and lead to big overshoots.
       * This tends to happen when the adjustment is very close to zero
       * and it oscillates a bit around zero.  Linear interpolation has
       * no issues, since it doesn't require smoothness.
       */
      sum = linear_interp(x, xy_pairs2, last_idx, min_idx);
      pair_ptr->fit_log_adjust = 0.5 * (pair_ptr->fit_log_adjust + sum);

#if 0
      printf("DEBUG\t%f\t%d\t%f\t%f\t%f\n",
             x, min_idx, xy_pairs2[min_idx].x, xy_pairs2[min_idx].y, pair_ptr->fit_log_adjust);
#endif
    }
  }


  if (xy_pairs)
    free(xy_pairs);

  if (xy_pairs2)
    free(xy_pairs2);
}


/* out of date, needs to be updated with latest stuff from regular function */
int refine_training_set(struct signal_pair **pair_ptrs,
                        int num_pairs,
                        struct signal_pair **filt_ptrs_train,
                        int num_pairs_train,
                        AFFY_ERROR *err)
{
  struct xy_pair    *xy_pairs = NULL;
  struct xy_pair    *xy_pairs2 = NULL;
  struct signal_pair *pair_ptr;
  double              x, y, old_x, old_y, sum;
  int                 num_pairs_train2 = 0;
  int                 num_pairs_train3 = 0;
  int                 i, num;
  int                 min_idx = 0, last_idx;
  double              rmsd;

  assert(pair_ptrs       != NULL);
  assert(filt_ptrs_train != NULL);

  xy_pairs = (struct xy_pair *) calloc(num_pairs_train,
                                       sizeof(struct xy_pair));
  if (xy_pairs == NULL)
    AFFY_HANDLE_ERROR("malloc failed", AFFY_ERROR_OUTOFMEM, err, 0);

  qsort(pair_ptrs,
        num_pairs,
        sizeof(struct signal_pair *),
        compare_sort_log_xy);
  qsort(filt_ptrs_train,
        num_pairs_train,
        sizeof(struct signal_pair *),
        compare_sort_log_xy);

  /* condense overlapping points */
  /* project onto average best fit line */
  old_x = -9E99;
  old_y = -9E99;
  for (i = 0, num_pairs_train2 = 0; i < num_pairs_train; i++)
  {
    pair_ptr = filt_ptrs_train[i];
    
    /* project x/y onto best fit avg line using norm_err */
    x = log(pair_ptr->sig1) + 0.5 * pair_ptr->norm_err;
    y = log(pair_ptr->sig2) - 0.5 * pair_ptr->norm_err;
    
    /* is new point different from last point? */
    /* take round off error into account */
    if (fabs(x - old_x) > 1E-14 || fabs(y - old_y) > 1E-14)
    {
      xy_pairs[num_pairs_train2].x   = x;
      xy_pairs[num_pairs_train2++].y = y;
    }
    
    old_x = x;
    old_y = y;
  }

  qsort(xy_pairs,
        num_pairs_train2,
        sizeof(struct xy_pair),
        compare_sort_xy_pair_by_x_plus_y);

  xy_pairs2 = (struct xy_pair *) calloc(num_pairs_train2,
                                        sizeof(struct xy_pair));
  if (xy_pairs2 == NULL)
    AFFY_HANDLE_ERROR("malloc failed", AFFY_ERROR_OUTOFMEM, err, 0);


  /* store as [log(x*y), log(x/y)] pairs */
  old_x = -9E99;
  sum = 0;
  num = 0;
  for (i = 0, num_pairs_train3 = -1; i < num_pairs_train2; i++)
  {
    x = xy_pairs[i].x + xy_pairs[i].y;
    
    if (x != old_x)
    {
      num_pairs_train3++;
      sum = 0;
      num = 0;
    }
    old_x = x;

    sum += xy_pairs[i].x - xy_pairs[i].y;	/* log(x*y), log(x/y) */
    num++;

    xy_pairs2[num_pairs_train3].x = x;
    
    /* average scales for identical points */
    if (num > 1)
      xy_pairs2[num_pairs_train3].y = sum / num;
    else
      xy_pairs2[num_pairs_train3].y = sum;
  }

  num_pairs_train3++;
  last_idx = num_pairs_train3 - 1;
  old_x = -9E99;


#if DEBUG_PRINT
  fprintf(stderr, "Training\t%d\t%d\t%d\n",
          num_pairs_train, num_pairs_train2, num_pairs_train3);
#endif

  for (i = 0; i < num_pairs; i++)
  {
    pair_ptr = pair_ptrs[i];
    x        = pair_ptr->log_xy;

    /* same X as last time, use the previous fit */
    if (x == old_x && i)
    {
      pair_ptr->fit_log_adjust = pair_ptrs[i-1]->fit_log_adjust;
      pair_ptr->norm_err = pair_ptr->fit_log_adjust - pair_ptr->log_adjust;

      continue;
    }

    old_x = x;

    /* before first training point */
    /* use average of end 10 adjustments */
    if (x < xy_pairs2[0].x)
    {
      int j;
    
      sum = 0.0;
      for (j = 0; j < 10 && j < last_idx; j++)
        sum += xy_pairs2[j].y;
      if (j)
        sum /= j;
      pair_ptr->fit_log_adjust = sum;
      pair_ptr->norm_err = pair_ptr->fit_log_adjust - pair_ptr->log_adjust;

      continue;
    }
        
    /* after last training point */
    /* use average of end 10 adjustments */
    if (x > xy_pairs2[last_idx].x)
    {
      int j;

      sum = 0.0;
      for (j = 0; j < 10 && j < last_idx; j++)
        sum += xy_pairs2[last_idx - j].y;
      if (j)
        sum /= j;
      pair_ptr->fit_log_adjust = sum;
      pair_ptr->norm_err = pair_ptr->fit_log_adjust - pair_ptr->log_adjust;

      continue;
    }
        
    /* skip to at or just after current point */
    while(xy_pairs2[min_idx].x < x)
      min_idx++;

    /* Use linear interpolation.
     * WARNING -- curve can have occasional too-rapid changes which
     * greatly confuse Lagrange interpolation and lead to big overshoots.
     * This tends to happen when the adjustment is very close to zero
     * and it oscillates a bit around zero.  Linear interpolation has
     * no issues, since it doesn't require smoothness.
     */
    pair_ptr->fit_log_adjust = linear_interp(x, 
                                             xy_pairs2, 
                                             last_idx, 
                                             min_idx);

    pair_ptr->norm_err = pair_ptr->fit_log_adjust - pair_ptr->log_adjust;

#if 0
    printf("DEBUG\t%f\t%d\t%f\t%f\t%f\n",
           x, min_idx, xy_pairs2[min_idx].x, xy_pairs2[min_idx].y, pair_ptr->fit_log_adjust);
#endif
  }


#if 1
  /* now that we've calculated norm_err for all points, calculate the RMSD */
  rmsd = 0.0;
  for (i = 0; i < num_pairs; i++)
    rmsd += pair_ptrs[i]->norm_err * pair_ptrs[i]->norm_err;

  /* GSM467598, GSM467552, GSM467594 */
  /* 5.0 fits diagonal, but doesn't prune much and hoses GSM467526 */
  /* sqrt(2) sd ensures >= 50% coverage (Chebyshev inequality) */
  rmsd = sqrt(2.0) * sqrt(rmsd / num_pairs);
#else
  rmsd = 0.0;
  for (i = 0; i < num_pairs_train; i++)
    rmsd += filt_ptrs_train[i]->norm_err * filt_ptrs_train[i]->norm_err;
  rmsd = 5*sqrt(rmsd / num_pairs_train);
#endif
  
#if 1
  /* fill the new training set with points < RMSD cutoff */
  for (i = 0, num_pairs_train = 0; i < num_pairs; i++)
  {
#if DEBUG_COLOR_IRANK
    pair_ptrs[i]->norm_err_scaled = pair_ptrs[i]->norm_err / rmsd;
#endif
  
    if (rmsd < 1E-5 || fabs(pair_ptrs[i]->norm_err) < rmsd)
      filt_ptrs_train[num_pairs_train++] = pair_ptrs[i];
  }
#else
  /* filter original training set of points > RMSD cutoff */
  for (i = 0, j = 0; i < num_pairs_train; i++)
    if (rmsd < 1E-5 || fabs(filt_ptrs_train[i]->norm_err) < rmsd)
      filt_ptrs_train[j++] = filt_ptrs_train[i];
  num_pairs_train = j;
#endif
  
  if (xy_pairs)
    free(xy_pairs);

  if (xy_pairs2)
    free(xy_pairs2);

  return num_pairs_train;
}


void fill_normalization_scales(char *filestem,
                               double *signals1,
                               double *signals2,
                               double *signals2_scales,
                               char *mask_array,
                               int num_spots,
                               double rank_frac_cutoff,
                               double rank_frac_cutoff2,
                               int condense_training_flag,
                               AFFY_COMBINED_FLAGS *f,
                               double *return_training_frac,
                               double *return_rmsd,
                               AFFY_ERROR *err)
{
  double min_sig1 = 9.0E8, min_sig2 = 9.0E8;
  int    bit16_flag1 = 1, bit16_flag2 = 1, i, *mempool;
   
  struct signal_pair  *signal_pairs = NULL;
  struct signal_pair  *pair_ptr;
  struct signal_pair **tmp_ptrs1 = NULL, **tmp_ptrs2 = NULL;
  struct signal_pair **tmp_ptrs3 = NULL, **tmp_ptrs4 = NULL;
  struct signal_pair **old_filt1, **old_filt2, **filt1, **filt2;

  int old_num_filtered  = -42;
  int num_filtered      = 0;
  int num_unpruned      = 0;
  int orig_num_unpruned = 0;
  int num_not_weak      = 0;
  int num_both_not_weak = 0;
   
  int    max_rank_diff;
  int    old_rank_diff_cutoff = num_spots + 1;
  int    rank_diff_cutoff = num_spots + 1;
  int    start, end;
  double rank_diff_cutoff_frac = 999;
  double old_rank_diff_cutoff_frac = 999;
   
  struct eqn_window *eqn_windows = NULL;
  int    num_eqns = 0;
  double rmsd;
  double global_scale = 0.0;
  int    global_scaling_flag      = f->iron_global_scaling_normalization;
  int    fit_both_x_y_flag        = f->iron_fit_both_x_y;
  double weight_exponent          = f->iron_weight_exponent;
   
#if DEBUG_FILE
  FILE *debug_file;
  double temp;
#endif

#if DEBUG_FIXED_RANK
  rank_frac_cutoff = 0.005;
#endif

  for (i = 0; i < num_spots; i++)
  {
    if (signals2[i] > MIN_SIGNAL)
    {
      num_not_weak++;
      
      if (signals1[i] > MIN_SIGNAL)
        num_both_not_weak++;
    }
  }

  /* check to see if we are normalizing vs. self
   * or if there are no good points
   */
  for (i = 0; i < num_spots; i++)
  {
    if (fabs(signals1[i] - signals2[i]) > 1E-5)
      break;
  }
  if (i == num_spots || num_both_not_weak == 0)
  {
    for (i = 0; i < num_spots; i++)
      signals2_scales[i] = 1.0;
   
    *return_training_frac = 1.0;
    *return_rmsd = 0.0;

    if (global_scaling_flag)
      fprintf(stderr, "GlobalScale:\t%s\t%f\t%f\t%d\t%d\t%d\t%d\t%f\n",
                      filestem,
                      1.0, 0.0,
                      num_not_weak,
                      num_both_not_weak,
                      num_not_weak, num_spots,
                      1.0);
    else if (f->iron_untilt_normalization)
      fprintf(stderr, "GlobalFitLine:\t%s\t%f\t%f\t%f\t%d\t%d\t%d\n",
                      filestem,
                      1.0, 0.0, 0.0,
                      num_both_not_weak, num_not_weak, num_spots);

    return;
  }

  mempool = h_malloc(sizeof(int));
  if (mempool == NULL)
    AFFY_HANDLE_ERROR_VOID("malloc failed", AFFY_ERROR_OUTOFMEM, err);

  signal_pairs = h_subcalloc(mempool, num_spots, sizeof(struct signal_pair));
  if (signal_pairs == NULL)
    AFFY_HANDLE_ERROR_GOTO("calloc failed",
                           AFFY_ERROR_OUTOFMEM,
                           err,
                           cleanup);

  tmp_ptrs1 = h_subcalloc(mempool, num_spots, sizeof(struct signal_pair *));
  if (tmp_ptrs1 == NULL)
    AFFY_HANDLE_ERROR_GOTO("calloc failed",
                           AFFY_ERROR_OUTOFMEM,
                           err,
                           cleanup);
 
  tmp_ptrs2 = h_subcalloc(mempool, num_spots, sizeof(struct signal_pair *));
  if (tmp_ptrs2 == NULL)
    AFFY_HANDLE_ERROR_GOTO("calloc failed",
                           AFFY_ERROR_OUTOFMEM,
                           err,
                           cleanup);

  tmp_ptrs3 = h_subcalloc(mempool, num_spots, sizeof(struct signal_pair *));
  if (tmp_ptrs3 == NULL)
    AFFY_HANDLE_ERROR_GOTO("calloc failed",
                           AFFY_ERROR_OUTOFMEM,
                           err,
                           cleanup);

  tmp_ptrs4 = h_subcalloc(mempool, num_spots, sizeof(struct signal_pair *));
  if (tmp_ptrs4 == NULL)
    AFFY_HANDLE_ERROR_GOTO("calloc failed",
                           AFFY_ERROR_OUTOFMEM,
                           err,
                           cleanup);
   
  old_filt1 = tmp_ptrs1;
  old_filt2 = tmp_ptrs2;
  filt1     = tmp_ptrs3;
  filt2     = tmp_ptrs4;

  /* assume 16-bit scanner if channel <= 65536 */
  for (i = 0; i < num_spots; i++)
  {
    if (signals1[i] > 65536)
      bit16_flag1 = 0;

    if (mask_array[i])
      continue;

    if (signals1[i] >= MIN_SIGNAL && signals1[i] < min_sig1)
      min_sig1 = signals1[i];
  }

  for (i = 0; i < num_spots; i++)
  {
    if (signals2[i] > 65536)
      bit16_flag2 = 0;

    if (mask_array[i])
      continue;

    if (signals2[i] >= MIN_SIGNAL && signals2[i] < min_sig2)
      min_sig2 = signals2[i];
  }

  for (i = 0; i < num_spots; i++)
  {
    pair_ptr            = signal_pairs + i;
    pair_ptr->index     = i;
    pair_ptr->sig2      = signals2[i];
    pair_ptr->sig1      = signals1[i];
    pair_ptr->weight    = 0;
    pair_ptr->n_windows = 0;
    
    /* avoid log(0) */
    if (pair_ptr->sig1 < MIN_SIGNAL)
    {
      pair_ptr->sig1 = MIN_SIGNAL;
    }

    if (pair_ptr->sig2 < MIN_SIGNAL)
    {
      pair_ptr->sig2 = MIN_SIGNAL;
    }

    pair_ptr->log_xy    = log(pair_ptr->sig1 * pair_ptr->sig2);

    /* filter control spots from training set */
    if (mask_array[i])
      continue;

    /* skip spots that are extremely dark in either channel */
#if DO_FLOOR
    if (pair_ptr->sig1 <= MIN_SIGNAL || pair_ptr->sig2 <= MIN_SIGNAL)
      continue;
#endif
    if (pair_ptr->sig1 <= MIN_SIGNAL && pair_ptr->sig2 > MIN_SIGNAL)
      continue;
    if (pair_ptr->sig2 <= MIN_SIGNAL && pair_ptr->sig1 > MIN_SIGNAL)
      continue;
    if (pair_ptr->sig1 <= min_sig1 || pair_ptr->sig2 <= min_sig2)
      continue;

    /* skip spots that are likely saturated in at least one channel */
    if ((bit16_flag1 && pair_ptr->sig1 >= 64000) ||
        (bit16_flag2 && pair_ptr->sig2 >= 64000))
      continue;

    /* store initial set of points */
    filt1[num_filtered] = filt2[num_filtered] = pair_ptr;

    pair_ptr->initial_set_flag = 1;

    num_filtered++;
  }

  /* uh oh, no good points left after first pass of filtering;
   * remove the minimum observed non-weak value and saturation filters
   */
  if (num_filtered == 0)
  {
    for (i = 0; i < num_spots; i++)
    {
      pair_ptr            = signal_pairs + i;
      pair_ptr->index     = i;
      pair_ptr->sig2      = signals2[i];
      pair_ptr->sig1      = signals1[i];
      pair_ptr->weight    = 0;
      pair_ptr->n_windows = 0;
      
      /* avoid log(0) */
      if (pair_ptr->sig1 < MIN_SIGNAL)
      {
        pair_ptr->sig1 = MIN_SIGNAL;
      }

      if (pair_ptr->sig2 < MIN_SIGNAL)
      {
        pair_ptr->sig2 = MIN_SIGNAL;
      }

      pair_ptr->log_xy    = log(pair_ptr->sig1 * pair_ptr->sig2);

      /* filter control spots from training set */
      if (mask_array[i])
        continue;

      /* skip spots that are extremely dark in either channel */
#if DO_FLOOR
      if (pair_ptr->sig1 <= MIN_SIGNAL || pair_ptr->sig2 <= MIN_SIGNAL)
        continue;
#endif
      if (pair_ptr->sig1 <= MIN_SIGNAL && pair_ptr->sig2 > MIN_SIGNAL)
        continue;
      if (pair_ptr->sig2 <= MIN_SIGNAL && pair_ptr->sig1 > MIN_SIGNAL)
        continue;
#if 0
      if (pair_ptr->sig1 <= min_sig1 || pair_ptr->sig2 <= min_sig2)
        continue;

      /* skip spots that are likely saturated in at least one channel */
      if ((bit16_flag1 && pair_ptr->sig1 >= 64000) ||
          (bit16_flag2 && pair_ptr->sig2 >= 64000))
        continue;
#endif

      /* store initial set of points */
      filt1[num_filtered] = filt2[num_filtered] = pair_ptr;

      pair_ptr->initial_set_flag = 1;

      num_filtered++;
    }
  }

  /* still no training points, exit without normalizing */
  if (num_filtered == 0)
  {
    for (i = 0; i < num_spots; i++)
      signals2_scales[i] = 1.0;
   
    *return_training_frac = 1.0;
    *return_rmsd = 0.0;

    if (global_scaling_flag)
      fprintf(stderr, "GlobalScale:\t%s\t%f\t%f\t%d\t%d\t%d\t%d\t%f\n",
                      filestem,
                      1.0, 0.0,
                      num_not_weak,
                      num_both_not_weak,
                      num_not_weak, num_spots,
                      1.0);
    else if (f->iron_untilt_normalization)
      fprintf(stderr, "GlobalFitLine:\t%s\t%f\t%f\t%f\t%d\t%d\t%d\n",
                      filestem,
                      1.0, 0.0, 0.0,
                      num_both_not_weak, num_not_weak, num_spots);

    /* free the memory we've allocated thus far */
    h_free(mempool);

    return;
  }

  /* condense identical points; they cause too many problems */
  if (condense_training_flag)
  {
    qsort(filt1, num_filtered, sizeof(struct signal_pair *), compare_sort_sig2);
    old_num_filtered = num_filtered;
    filt2[0] = filt1[0];
    for (i = 1, num_filtered = 1; i < old_num_filtered; i++)
    {
        /* skip if identical to old point */
        if (filt1[i]->sig1 == filt1[i-1]->sig1 &&
            filt1[i]->sig2 == filt1[i-1]->sig2)
        {
            continue;
        }

        filt2[num_filtered++] = filt1[i];
    }
    /* copy new filtered pointers into old filtered pointers */
    memcpy(filt1, filt2, num_filtered * sizeof(struct signal_pair *));
  }

  qsort(filt1, num_filtered, sizeof(struct signal_pair *), compare_sort_sig1);
  qsort(filt2, num_filtered, sizeof(struct signal_pair *), compare_sort_sig2);

  num_unpruned = num_filtered;
  orig_num_unpruned = num_unpruned;

  /* iteratively prune training spots */
  while(num_filtered * rank_diff_cutoff_frac > 1.0 + 1E-5 &&
        (num_filtered != old_num_filtered ||
         rank_diff_cutoff_frac >= rank_frac_cutoff + 1E-5))
  {
    old_rank_diff_cutoff = rank_diff_cutoff;
    old_num_filtered     = num_filtered;

    if (old_filt1 == tmp_ptrs1)
    {
      old_filt1 = tmp_ptrs3;
      old_filt2 = tmp_ptrs4;
      filt1     = tmp_ptrs1;
      filt2     = tmp_ptrs2;
    }
    else
    {
      old_filt1 = tmp_ptrs1;
      old_filt2 = tmp_ptrs2;
      filt1     = tmp_ptrs3;
      filt2     = tmp_ptrs4;
    }

    /* store ranks */
    for (i = 0; i < old_num_filtered; i++)
    {
      old_filt1[i]->rank1 = i;
      old_filt2[i]->rank2 = i;
    }
   
    /* find max rank diff */
    for (max_rank_diff = 0, i = 0; i < old_num_filtered; i++)
    {
      old_filt1[i]->rank_diff =
                        labs(old_filt1[i]->rank1 - old_filt1[i]->rank2);
           
#if DEBUG_COLOR_IRANK
      old_filt1[i]->irank_frac = old_filt1[i]->rank_diff / (double) old_num_filtered;
      if (old_num_filtered == num_unpruned)
        old_filt1[i]->irank_frac_0 = old_filt1[i]->irank_frac;
#endif
       
      if (old_filt1[i]->rank_diff > max_rank_diff)
        max_rank_diff = old_filt1[i]->rank_diff;
    }
       
    /* set cutoff to max observed minus 0.5% */
    old_rank_diff_cutoff_frac = rank_diff_cutoff_frac;
    rank_diff_cutoff_frac = (double) max_rank_diff /
                            old_num_filtered - 0.005;
   
    /* floor the cutoff at the defined floor (0.01 recommended) */
    if (rank_diff_cutoff_frac < rank_frac_cutoff)
      rank_diff_cutoff_frac = rank_frac_cutoff;

#if DEBUG_FIXED_RANK
    /* just one pass of rank order filtering, for demonstration */
    rank_diff_cutoff_frac = rank_frac_cutoff;
#endif   

    rank_diff_cutoff = (int)(old_num_filtered * rank_diff_cutoff_frac + 0.5);

    /* prune by rank diff */
    for (num_filtered = 0, end = i = 0, start = -1; i < old_num_filtered; i++)
    {
      if (old_filt1[i]->rank_diff >= rank_diff_cutoff)
      {
        if (start >= 0)
        {
          int n = end - start + 1;

          memcpy(filt1 + num_filtered, old_filt1 + start,
            n * sizeof(struct signal_pair *));
          num_filtered += n;
        }
       
        start = -1;

        continue;
      }
     
      if (start < 0)
        start = i;

      end = i;
    }
    if (start >= 0)
    {
      memcpy(filt1 + num_filtered, old_filt1 + start,
        (end - start + 1) * sizeof(struct signal_pair *));
    }

    for (num_filtered = 0, end = i = 0, start = -1; i < old_num_filtered; i++)
    {
      if (old_filt2[i]->rank_diff >= rank_diff_cutoff)
      {
        if (start >= 0)
        {
          int n = end - start + 1;

          memcpy(filt2 + num_filtered, old_filt2 + start,
            n * sizeof(struct signal_pair *));
          num_filtered += n;
        }
       
        start = -1;

        continue;
      }
     
      if (start < 0)
        start = i;

      end = i;
    }
    if (start >= 0)
    {
      int n = end - start + 1;

      memcpy(filt2 + num_filtered, old_filt2 + start,
        n * sizeof(struct signal_pair *));
      num_filtered += n;
    }
  }

  /* we've pruned too much, back up an iteration */
  if (num_filtered * rank_diff_cutoff_frac < 1.0 + 1E-5)
  {
    rank_diff_cutoff      = old_rank_diff_cutoff;
    num_filtered          = old_num_filtered;
    rank_diff_cutoff_frac = old_rank_diff_cutoff_frac;

    if (old_filt1 == tmp_ptrs1)
    {
      old_filt1 = tmp_ptrs3;
      old_filt2 = tmp_ptrs4;
      filt1     = tmp_ptrs1;
      filt2     = tmp_ptrs2;
    }
    else
    {
      old_filt1 = tmp_ptrs1;
      old_filt2 = tmp_ptrs2;
      filt1     = tmp_ptrs3;
      filt2     = tmp_ptrs4;
    }
  }
  
#if DEBUG_PRINT
  fprintf(stderr,
          "IRank:\t%d\t%d\t%d\t%f\t%f\n",
          num_spots, num_unpruned, num_filtered,
          rank_diff_cutoff_frac,
          (double)num_filtered / (double)num_unpruned);
#endif

  eqn_windows = (struct eqn_window *)h_subcalloc(mempool,
                                                 num_filtered,
                                                 sizeof(struct eqn_window));
  if (eqn_windows == NULL)
    AFFY_HANDLE_ERROR_GOTO("calloc failed",
                           AFFY_ERROR_OUTOFMEM,
                           err,
                           cleanup);
   
  /* initialize values for smoothed piecewise linear fit (geometric) */
  /* calculate for ALL points, since we will eventually use them later */
  for (i = 0; i < num_spots; i++)
  {
    pair_ptr = &signal_pairs[i];
    pair_ptr->log_adjust = log(pair_ptr->sig1 / pair_ptr->sig2);
  }

  /* windowed linear fits, log(x/y) vs. log(x*y) */
  /* use filt2, since filt1 is already sorted on X which we'll use later */

  num_eqns = fill_geometric_eqn_windows(eqn_windows, filt2, num_filtered,
                                        f->iron_fit_window_frac, 1,
                                        weight_exponent);
  smooth_geometric_fits(eqn_windows, num_eqns, filt2, num_filtered);

#if 0
  /* calculate post-normalized rmsd from center line */
  rmsd = 0.0;
  for (i = 0; i < num_filtered; i++)
    rmsd += filt1[i]->norm_err * filt1[i]->norm_err;

  if (num_filtered)
    rmsd = sqrt(rmsd / num_filtered);
 
  /* filter training points beyond 1 rmsd, store in filt2 */
  /* filt2 will contain final training set, filt1 contains previous set */
  old_num_filtered = num_filtered;

  num_filtered = 0;
  for (i = 0; i < old_num_filtered; i++)
  {
    if (fabs(filt1[i]->norm_err) < 5*rmsd + 1E-5)
      filt2[num_filtered++] = filt1[i];
  }

  /* refit on final reduced training set */
  num_eqns = fill_geometric_eqn_windows(eqn_windows, filt2,
                                        num_filtered, f->iron_fit_window_frac,
                                        2, weight_exponent);
  smooth_geometric_fits(eqn_windows, num_eqns, filt2, num_filtered);
#endif

#if SECOND_PASS_TRAIN && DEBUG_FIXED_RANK == 0
  /* refit with larger swath around best fit line */
  /* store pointers to all non-bad points into filt1 */
  num_unpruned = 0;
  for (i = 0; i < num_spots; i++)
    if (signal_pairs[i].initial_set_flag)
      filt1[num_unpruned++] = signal_pairs + i;

  num_filtered = refine_training_set(filt1, num_unpruned, filt2, num_filtered,
                                     err);

#if DEBUG_COLOR_IRANK
  for (i = 0; i < num_spots; i++)
  {
    signal_pairs[i].initial_set_flag = 0;
    signal_pairs[i].irank_frac = 0;
    signal_pairs[i].irank_frac_0 = 0;
  }
  for (i = 0; i < num_filtered; i++)
    filt2[i]->initial_set_flag = 1;
#endif

  for (i = 0; i < num_filtered; i++)
    filt1[i] = filt2[i];


  qsort(filt1, num_filtered, sizeof(struct signal_pair *), compare_sort_sig1);
  qsort(filt2, num_filtered, sizeof(struct signal_pair *), compare_sort_sig2);

  old_num_filtered = -42;
  num_unpruned = num_filtered;

  /* iteratively prune training spots */
  while(num_filtered * rank_diff_cutoff_frac > 1.0 + 1E-5 &&
        (num_filtered != old_num_filtered ||
         rank_diff_cutoff_frac >= rank_frac_cutoff + 1E-5))
  {
    old_rank_diff_cutoff = rank_diff_cutoff;
    old_num_filtered     = num_filtered;
       
    if (old_filt1 == tmp_ptrs1)
    {
      old_filt1 = tmp_ptrs3;
      old_filt2 = tmp_ptrs4;
      filt1     = tmp_ptrs1;
      filt2     = tmp_ptrs2;
    }
    else
    {
      old_filt1 = tmp_ptrs1;
      old_filt2 = tmp_ptrs2;
      filt1     = tmp_ptrs3;
      filt2     = tmp_ptrs4;
    }

    /* store ranks */
    for (i = 0; i < old_num_filtered; i++)
    {
      old_filt1[i]->rank1 = i;
      old_filt2[i]->rank2 = i;
    }
   
    /* find max rank diff */
    for (max_rank_diff = 0, i = 0; i < old_num_filtered; i++)
    {
      old_filt1[i]->rank_diff =
                        labs(old_filt1[i]->rank1 - old_filt1[i]->rank2);
           
#if DEBUG_COLOR_IRANK
      old_filt1[i]->irank_frac = rank_diff / (double) old_num_filtered;
      if (old_num_filtered == num_unpruned)
        old_filt1[i]->irank_frac_0 = old_filt1[i]->irank_frac;
#endif
       
      if (old_filt1[i]->rank_diff > max_rank_diff)
        max_rank_diff = old_filt1[i]->rank_diff;
    }
       
    /* set cutoff to max observed minus 0.5% */
    old_rank_diff_cutoff_frac = rank_diff_cutoff_frac;
    rank_diff_cutoff_frac = (double) max_rank_diff /
                            old_num_filtered - 0.005;
   
    /* floor the cutoff at the defined floor (0.1 recommended) */
    if (rank_diff_cutoff_frac < rank_frac_cutoff2)
      rank_diff_cutoff_frac = rank_frac_cutoff2;

#if DEBUG_FIXED_RANK
    /* just one pass of rank order filtering, for demonstration */
    rank_diff_cutoff_frac = rank_frac_cutoff2;
#endif   

    rank_diff_cutoff = (int)(old_num_filtered * rank_diff_cutoff_frac + 0.5);

    /* prune by rank diff */
    for (num_filtered = 0, end = i = 0, start = -1; i < old_num_filtered; i++)
    {
      if (old_filt1[i]->rank_diff >= rank_diff_cutoff)
      {
        if (start >= 0)
        {
          int n = end - start + 1;

          memcpy(filt1 + num_filtered, old_filt1 + start,
            n * sizeof(struct signal_pair *));
          num_filtered += n;
        }
       
        start = -1;

        continue;
      }
     
      if (start < 0)
        start = i;

      end = i;
    }
    if (start >= 0)
    {
      memcpy(filt1 + num_filtered, old_filt1 + start,
        (end - start + 1) * sizeof(struct signal_pair *));
    }

    for (num_filtered = 0, end = i = 0, start = -1; i < old_num_filtered; i++)
    {
      if (old_filt2[i]->rank_diff >= rank_diff_cutoff)
      {
        if (start >= 0)
        {
          int n = end - start + 1;

          memcpy(filt2 + num_filtered, old_filt2 + start,
            n * sizeof(struct signal_pair *));
          num_filtered += n;
        }
       
        start = -1;

        continue;
      }
     
      if (start < 0)
        start = i;

      end = i;
    }
    if (start >= 0)
    {
      int n = end - start + 1;

      memcpy(filt2 + num_filtered, old_filt2 + start,
        n * sizeof(struct signal_pair *));
      num_filtered += n;
    }
  }

  /* we've pruned too much, back up an iteration */
  if (num_filtered * rank_diff_cutoff_frac < 1.0 + 1E-5)
  {
    rank_diff_cutoff      = old_rank_diff_cutoff;
    num_filtered          = old_num_filtered;
    rank_diff_cutoff_frac = old_rank_diff_cutoff_frac;

    if (old_filt1 == tmp_ptrs1)
    {
      old_filt1 = tmp_ptrs3;
      old_filt2 = tmp_ptrs4;
      filt1     = tmp_ptrs1;
      filt2     = tmp_ptrs2;
    }
    else
    {
      old_filt1 = tmp_ptrs1;
      old_filt2 = tmp_ptrs2;
      filt1     = tmp_ptrs3;
      filt2     = tmp_ptrs4;
    }
  }

#if DEBUG_PRINT
  fprintf(stderr,
          "IRank:\t%d\t%d\t%d\t%f\t%f\n",
          num_spots, num_unpruned, num_filtered,
          rank_diff_cutoff_frac,
          (double)num_filtered / (double)num_unpruned);
#endif

  /* reallocate equation windows to hold new larger training set */
  fprintf(stderr, "NumFiltered\t%d\n", num_filtered);
  eqn_windows = (struct eqn_window *)h_realloc(eqn_windows,
                                               num_filtered *
                                               sizeof(struct eqn_window));
  if (eqn_windows == NULL)
    AFFY_HANDLE_ERROR_GOTO("calloc failed",
                           AFFY_ERROR_OUTOFMEM,
                           err,
                           cleanup);

  /* refit on final reduced training set */
  num_eqns = fill_geometric_eqn_windows(eqn_windows, filt2,
                                        num_filtered, f->iron_fit_window_frac,
                                        2, weight_exponent);
  smooth_geometric_fits(eqn_windows, num_eqns, filt2, num_filtered);
#endif

  /* remember which points were in the irank training set */
  for (i = 0; i < num_filtered; i++)
    filt2[i]->irank_flag = 1;

  /* calculate adjustments for ALL points */
  /* store pointers to ALL points in filt1 */
  for (i = 0; i < num_spots; i++)
    filt1[i] = signal_pairs + i;

  interpolate_final_scales(filt1, num_spots, filt2, num_filtered,
                           fit_both_x_y_flag, err);
  
  /* use a single global scaling factor, rather than non-linear scaling */
  if (global_scaling_flag)
  {
#if 0
    /* refit linear line on final reduced training set, single full window */
    /* appears to give a worse fit than using non-linear fit (?) */
    num_eqns = fill_geometric_eqn_windows(eqn_windows, filt2,
                                          num_filtered, 1.0, 1,
                                          weight_exponent);
    smooth_geometric_fits(eqn_windows, num_eqns, filt2, num_filtered);
    interpolate_final_scales(filt1, num_spots, filt2, num_filtered,
                             fit_both_x_y_flag, err);
#endif

    affy_int32 count = 0;
    
    for (i = 0; i < num_spots; i++)
    {
      if (signal_pairs[i].irank_flag)
      {
        global_scale += signal_pairs[i].fit_log_adjust;
        count++;
      }
    }

    global_scale = exp(global_scale / count);

    fprintf(stderr, "GlobalScale:\t%s\t%f\t%f\t%d\t%d\t%d\t%d\t%f\n",
                    filestem,
                    global_scale, log(global_scale) / log(2.0),
                    count, num_both_not_weak, num_not_weak, num_spots,
                    (double) count / (double) num_both_not_weak);
  }
  /* use a single line fit to the entire training set */
  else if (f->iron_untilt_normalization)
  {
    affy_int32 count = 0;

    /* refit linear line on final reduced training set, single full window */
    num_eqns = fill_geometric_eqn_windows(eqn_windows, filt2,
                                          num_filtered, 1.0, 1,
                                          weight_exponent);
    smooth_geometric_fits(eqn_windows, num_eqns, filt2, num_filtered);
    interpolate_final_scales(filt1, num_spots, filt2, num_filtered,
                             fit_both_x_y_flag, err);

    /* average adjustements together for printing QC info */
    global_scale = 0.0;
    for (i = 0; i < num_spots; i++)
    {
      if (signal_pairs[i].irank_flag)
      {
        global_scale += signal_pairs[i].fit_log_adjust;
        count++;
      }
    }

    global_scale = exp(global_scale / count);

    fprintf(stderr, "GlobalFitLine:\t%s\t%f\t%f\t%f\t%d\t%d\t%d\n",
                    filestem,
                    1.0 / global_scale,
                    -log(global_scale) / log(2.0),
                    -(180.0 * atan(eqn_windows[0].slope) / M_PI),
                    num_both_not_weak, num_not_weak, num_spots);
  }

  /* store scaling multipliers in signals2_scales array */
  for (i = 0; i < num_spots; i++)
  {
#if DO_FLOOR
    /* zero out extremely weak signals */
    if (signal_pairs[i].sig2 <= MIN_SIGNAL)
      signals2_scales[i] = 0;
    else
#endif

    if (global_scaling_flag)
      signals2_scales[i] = global_scale;
    else
      signals2_scales[i] = exp(signal_pairs[i].fit_log_adjust);
  }


#if DEBUG_FILE
  debug_file = fopen("irank_set.txt", "wb");
  if (!debug_file)
    fprintf(stderr, "ERROR -- Can't open irank_set.txt\n");
 
  fprintf(debug_file, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
          "ProbeID", "log10_X", "log10_Y", "log10_Y_norm",
          "Weight",
          "IRankSet", "InitialRankSet", "IRankFrac", "IRankFrac_0",
          "X0_proj", "Y0_proj", "Y0_proj_norm",
          "log10(X*Y)", "log10(X/Y)", "log10(X/Y_norm)",
          "log10(X0_proj*Y0_proj)", "log10(X0_proj/Y0_proj)",
          "log10(X0_proj_norm/Y0_proj_norm)");

  /* find lowest non-zero rank_diff */
  temp = 1.0;
  for (i = 0; i < num_spots; i++)
  {
      if (signal_pairs[i].irank_frac > 0 &&
          signal_pairs[i].irank_frac < temp)
      {
          temp = signal_pairs[i].irank_frac;
      }
      if (signal_pairs[i].irank_frac_0 > 0 &&
          signal_pairs[i].irank_frac_0 < temp)
      {
          temp = signal_pairs[i].irank_frac_0;
      }
  }

  for (i = 0; i < num_spots; i++)
  {
    double temp_signal, irank_frac, irank_frac_0;
    double x_train, y_train, y_train_norm;
 
    if (mask_array[i])
      continue;
    
    temp_signal = signal_pairs[i].sig2 * signals2_scales[i];
    if (temp_signal < MIN_SIGNAL)
      temp_signal = MIN_SIGNAL;

    x_train = log(signal_pairs[i].sig1) / log(10.0);
    y_train = log(signal_pairs[i].sig2) / log(10.0);
    y_train_norm = log(temp_signal) / log(10.0);

    if (signal_pairs[i].irank_flag)
    {
      x_train = (log(signal_pairs[i].sig1) + 0.5 * signal_pairs[i].norm_err) /
              log(10.0);
      y_train = (log(signal_pairs[i].sig2) - 0.5 * signal_pairs[i].norm_err) /
              log(10.0);
      y_train_norm = x_train;
    }
    
    irank_frac = signal_pairs[i].irank_frac;
    irank_frac_0 = signal_pairs[i].irank_frac_0;
    
    if (irank_frac < temp) irank_frac = temp;
    if (irank_frac_0 < temp) irank_frac_0 = temp;
 
    fprintf(debug_file, "%d\t%f\t%f\t%f\t%e\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
            i,
            log(signal_pairs[i].sig1) / log(10.0),
            log(signal_pairs[i].sig2) / log(10.0),
            log(temp_signal) / log(10.0),
#if 0
            fabs(signal_pairs[i].norm_err_scaled),
#else
            signal_pairs[i].weight,
#endif

#if DEBUG_COLOR_IRANK
            signal_pairs[i].irank_flag,
            signal_pairs[i].initial_set_flag,
            log10(irank_frac),
            log10(irank_frac_0),
#else
            0.0, 0.0, 0.0, 0.0,
#endif
            x_train,
            y_train,
            y_train_norm,
            log(signal_pairs[i].sig1 * signal_pairs[i].sig2) / log(10.0),
            log(signal_pairs[i].sig1 / signal_pairs[i].sig2) / log(10.0),
            log(signal_pairs[i].sig1 / temp_signal) / log(10.0),
            x_train + y_train,
            x_train - y_train,
            x_train - y_train_norm);
    /*
     * some editors' auto-indentation freaks out unless you remove the
     * ); from the #ifdef
     */
  }
  fclose(debug_file);
#endif
   
  /* return array similarity metrics */
  rmsd = 0;
  for (i = 0; i < num_spots; i++)
    if (signal_pairs[i].initial_set_flag)
      rmsd += signal_pairs[i].fit_log_adjust * signal_pairs[i].fit_log_adjust;
 
  if (orig_num_unpruned)
    rmsd = sqrt(rmsd / orig_num_unpruned);

  *return_rmsd          = rmsd / log(10.0);
  *return_training_frac = (double) num_filtered / (double) orig_num_unpruned;

#if DEBUG_PRINT
  fprintf(stderr, "SimilarityMetrics:\tTrain\t%f\tRMSD\t%f\n",
     *return_training_frac, *return_rmsd);
#endif

cleanup:
  h_free(mempool);

#if DEBUG_DIE_EARLY
  exit(0);
#endif
}
