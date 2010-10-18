
/*
  Trim a string of whitespace from both ends.
*/
#include <string.h>
#include <ctype.h>

#include "utils.h"

char *trim(char *p)
{
  register char *l;

  if (p == NULL || *p == '\0')
    return p;

  l = ltrim(p);
  rtrim(l);

  return (l);
}
