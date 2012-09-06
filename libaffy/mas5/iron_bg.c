#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "affy.h"

#define MIN_SIGNAL       1E-5	/* no lower than 0.05 */

void affy_iron_background_correction_probeset(AFFY_CHIPSET *cs,
                                              AFFY_COMBINED_FLAGS *f,
                                              AFFY_ERROR *err)
{
  AFFY_CDFFILE       *cdf;
  AFFY_CELL         **chip_data;
  AFFY_CHIP          *chip;
  int                *mempool;
  int                 numprobes, numprobesets;
  int                 i;
  double              log2;

  assert(cs              != NULL);
  assert(cs->cdf         != NULL);
  assert(cs->chip        != NULL);

  if (cs->num_chips == 0)
    return;

  cdf          = cs->cdf;
  numprobes    = cdf->numprobes;
  numprobesets = cdf->numprobesets;

  mempool = h_malloc(sizeof(int));
  if (mempool == NULL)
    AFFY_HANDLE_ERROR_GOTO("malloc failed", AFFY_ERROR_OUTOFMEM, err, cleanup);

  log2 = log(2.0);
  for (i = 0; i < cs->num_chips; i++)
  {
    double *probe_set_ptr;
  
    chip = cs->chip[i];
    chip_data = chip->cel->data;
    probe_set_ptr = chip->probe_set;
  }

cleanup:
  h_free(mempool);
}
