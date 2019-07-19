#include <math.h>

double calculate_pearson_r_float(float *array1, float *array2, int n)
{
  double x_avg = 0, y_avg = 0;
  double sum_xy_diff = 0, sum_x2_diff = 0, sum_y2_diff = 0;
  double x_diff, y_diff;
  double temp;
  int i;

  for (i = 0; i < n; i++)
  {
    x_avg += array1[i];
    y_avg += array2[i];
  }
  x_avg /= n;
  y_avg /= n;

  for (i = 0; i < n; i++)
  {
    x_diff = array1[i] - x_avg;
    sum_x2_diff += x_diff * x_diff;

    y_diff = array2[i] - y_avg;
    sum_y2_diff += y_diff * y_diff;

    sum_xy_diff += x_diff * y_diff;
  }

  if (sum_x2_diff && sum_y2_diff)
  {
    sum_x2_diff = sqrt(sum_x2_diff);
    sum_y2_diff = sqrt(sum_y2_diff);
    temp = sum_x2_diff * sum_y2_diff;
  
    if (temp)
    {
      temp = sum_xy_diff / temp;
      
      /* round off errors can lead to slightly greater than 1 */
      if (temp >  1.0)
        return    1.0;
      if (temp < -1.0)
        return   -1.0;
      
      return temp;
    }
  }

  return 0;
}


double calculate_pearson_r_float_skip_weak(float *array1, float *array2,
                                           int n)
{
  double x_avg = 0, y_avg = 0;
  double sum_xy_diff = 0, sum_x2_diff = 0, sum_y2_diff = 0;
  double x_diff, y_diff;
  double temp;
  int i;
  int count = 0;

  for (i = 0; i < n; i++)
  {
    if (array1[i] <= 0 || array2[i] <= 0)
      continue;
  
    x_avg += array1[i];
    y_avg += array2[i];
    
    count++;
  }
  x_avg /= count;
  y_avg /= count;

  for (i = 0; i < n; i++)
  {
    if (array1[i] <= 0 || array2[i] <= 0)
      continue;

    x_diff = array1[i] - x_avg;
    sum_x2_diff += x_diff * x_diff;

    y_diff = array2[i] - y_avg;
    sum_y2_diff += y_diff * y_diff;

    sum_xy_diff += x_diff * y_diff;
  }

  if (sum_x2_diff && sum_y2_diff)
  {
    sum_x2_diff = sqrt(sum_x2_diff);
    sum_y2_diff = sqrt(sum_y2_diff);
    temp = sum_x2_diff * sum_y2_diff;
  
    if (temp)
    {
      temp = sum_xy_diff / temp;
      
      /* round off errors can lead to slightly greater than 1 */
      if (temp >  1.0)
        return    1.0;
      if (temp < -1.0)
        return   -1.0;
      
      return temp;
    }
  }

  return 0;
}


double calculate_pearson_r_double(double *array1, double *array2, int n)
{
  double x_avg = 0, y_avg = 0;
  double sum_xy_diff = 0, sum_x2_diff = 0, sum_y2_diff = 0;
  double x_diff, y_diff;
  double temp;
  int i;

  for (i = 0; i < n; i++)
  {
    x_avg += array1[i];
    y_avg += array2[i];
  }
  x_avg /= n;
  y_avg /= n;

  for (i = 0; i < n; i++)
  {
    x_diff = array1[i] - x_avg;
    sum_x2_diff += x_diff * x_diff;

    y_diff = array2[i] - y_avg;
    sum_y2_diff += y_diff * y_diff;

    sum_xy_diff += x_diff * y_diff;
  }

  if (sum_x2_diff && sum_y2_diff)
  {
    sum_x2_diff = sqrt(sum_x2_diff);
    sum_y2_diff = sqrt(sum_y2_diff);
    temp = sum_x2_diff * sum_y2_diff;
  
    if (temp)
    {
      temp = sum_xy_diff / temp;
      
      /* round off errors can lead to slightly greater than 1 */
      if (temp >  1.0)
        return    1.0;
      if (temp < -1.0)
        return   -1.0;
      
      return temp;
    }
  }

  return 0;
}
