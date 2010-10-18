#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>

#include "util_log.h"

void status(const char *msg, ...)
{
  va_list ap;

  va_start(ap, msg);
  vfprintf(stderr, msg, ap);
  va_end(ap);

  fprintf(stderr, "\n");
}
