
#include <utils.h>

/* Split the line into key, value based on the = */
int split(char *str, char **kv, int split_val, int maxsplit)
{
  char *p, *s;
  int   i = 1;

  assert(maxsplit > 0);
  assert(str   != NULL);
  assert(kv    != NULL);

  p = s = kv[0] = str;

  while (((p = strchr(s, split_val)) != NULL) && (i < maxsplit))
  {
    s  = p + 1;
    *p = '\0';
    kv[i++] = s;
  }
 
  return (i);
}

