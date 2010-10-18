
/**
 * A generic routine to compare two doubles numerically. Useful
 * for qsort() routines.
 */
int dcompare(const void *p1, const void *p2)
{
  double n1 = *((double *)p1);
  double n2 = *((double *)p2);

  if (n1 < n2)
    return -1;
  if (n1 > n2)
    return 1;
  return 0;
}
