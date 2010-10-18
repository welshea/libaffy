
/**
 * A generic routine to compare two floats numerically. Useful
 * for qsort() routines.
 * 
 * NB: Should probably consider equivalence a bit more loosely (e.g. epsilon)
 */

int fcompare(const void *p1, const void *p2)
{
  float n1 = *((float *)p1);
  float n2 = *((float *)p2);

  if (n1 < n2)
    return -1;
  if (n1 > n2)
    return 1;
  return 0;
}
