
/**
 * A generic routine to compare two ints numerically. Useful
 * for qsort() routines.
 */
int icompare(const void *p1, const void *p2)
{
  int n1 = *((int *)p1);
  int n2 = *((int *)p2);

  if (n1 < n2)
    return -1;
  if (n1 > n2)
    return 1;
  return 0;
}
