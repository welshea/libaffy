#include "utils.h"

void free_matrix(double **x)
{
  if (x == NULL)
    return;

  free(x[0]);
  free(x);
}
