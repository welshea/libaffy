#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "util_log.h"

#define PB_NUM_TICKS 20

void pb_init(LIBUTILS_PB_STATE *pbs)
{
  if (pbs == NULL)
    return;

  pbs->depth = 0;
}

void pb_cleanup(LIBUTILS_PB_STATE *pbs)
{
  return;
}

void pb_begin(LIBUTILS_PB_STATE *pbs, unsigned int max, char *title, ...)
{
  va_list ap;

  /* Return safely if there is no state available. */
  if (pbs == NULL)
    return;

  assert(pbs->depth < LIBUTILS_MAX_PB_DEPTH);

  pbs->cur_ticks[pbs->depth]     = 0;
  pbs->tick_interval[pbs->depth] = ceil((double)max / PB_NUM_TICKS);
  pbs->max[pbs->depth]=max;
  pbs->depth++;
  
  fprintf(stderr, "[");

  if (title != NULL)
  {
    va_start(ap, title);
  
    vfprintf(stderr, title, ap);
    
    va_end(ap);
  }
}

void pb_tick(LIBUTILS_PB_STATE *pbs, unsigned int tick_sz, char *msg, ...)
{
  int i;

  if (pbs == NULL)
    return;

  assert(pbs->depth > 0);

  i = pbs->depth - 1;

  pbs->cur_ticks[i] += tick_sz;

  while (pbs->cur_ticks[i] > pbs->tick_interval[i])
  {
    fprintf(stderr, ".");
    pbs->cur_ticks[i] -= pbs->tick_interval[i];
  }
}

void pb_msg(LIBUTILS_PB_STATE *pbs, char *msg, ...)
{
  va_list ap;

  assert(msg != NULL);
  
  /* Return safely if there is no state available. */
  if (pbs == NULL)
    return;

  va_start(ap, msg);

  fprintf(stderr, "(");
  vfprintf(stderr, msg, ap);
  fprintf(stderr, ")");

  va_end(ap);
}

void pb_finish(LIBUTILS_PB_STATE *pbs, char *msg, ...)
{
  va_list ap;

  assert(msg != NULL);

  /* Return safely if there is no state available. */
  if (pbs == NULL)
    return;

  assert(pbs->depth > 0);

  if (msg != NULL)
  {
    va_start(ap, msg);

    fprintf(stderr, "(");
    vfprintf(stderr, msg, ap);
    fprintf(stderr, ")");
    
    va_end(ap);
  }

  fprintf(stderr, "]");

  if (pbs->depth-- == 1)
    fprintf(stderr, "\n");
}
