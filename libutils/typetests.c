/* 
 * Simple sanity checks for the runtime environment type sizes. 
 */

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
  if (sizeof(float) != 4)
  {    
    printf("This system does not use 32-bit floats! "
           "Some assumptions may be invalid.\n");
    exit(EXIT_FAILURE);
  }

  if (sizeof(double) != 8)
  {
    printf("This system does not use 64-bit doubles! "
           "Some assumptions may be invalid.\n");
    exit(EXIT_FAILURE);
  }

  exit(EXIT_SUCCESS);
}
