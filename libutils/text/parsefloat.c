
#include "utils.h"

/*
 * Simplistic routine to parse a floating point number while hiding the
 * annoying strtod() API.
 *
 * The converted value from p_str is placed in p_dest.  0 is returned on
 * success, -1 on failure.  The contents of p_dest are undefined on error.
 */
int parsefloat(char *p_str, double *p_dest)
{
  char *err;

  assert(p_dest != NULL);
  assert(p_str  != NULL);

  *p_dest = strtod(p_str, &err);

  if ((*p_dest == 0) && (err == p_str))
    return (-1);
  else
    return (0);
}
