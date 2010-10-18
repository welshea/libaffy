#include "utils.h"

/**
 * Convenient wrapper for create a 2D matrix without excessive 
 * calloc.
 */

double **create_matrix(unsigned int rows, unsigned int columns)
{
  double     **x;
  unsigned int i;

  x = calloc(rows, sizeof(double *));
  if (x == NULL)
    return (NULL);

  x[0] = calloc(rows * columns, sizeof(double));
  if (x[0] == NULL)
  {
    free(x);
    return (NULL);
  }

  for (i = 1; i < rows; i++)
  {
    x[i] = x[0] + (i * columns);
  }

  return (x);
}
