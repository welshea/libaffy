
/**************************************************************************
 *
 * Filename: print_corrupt_chips.c
 * 
 * Purpose:  Output list of corrupt CEL files
 *
 * Creation: May 9th, 2013
 *
 * Author:   Eric A. Welsh
 *
 *
 * Update History
 * --------------
 * 05/09/13: Creation.  (EAW)
 *
 **************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "affy.h"

void print_corrupt_chips_to_stderr(AFFY_CHIPSET *cs)
{
  unsigned int i;

  assert(cs              != NULL);
  assert(cs->chip        != NULL);

  if (cs->num_chips == 0)
    return;

  for (i = 0; i < cs->num_chips; i++)
  {
    if (cs->chip[i] && cs->chip[i]->cel && cs->chip[i]->cel->corrupt_flag)
    {
      fprintf(stderr, "Corrupt CEL file: %s\n", cs->chip[i]->cel->filename);
    }
  }
}
