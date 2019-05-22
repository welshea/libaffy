
/**************************************************************************
 *
 * Filename:  load_chipset.c
 *
 * Purpose:   Routines for loading chipset data.
 *
 * Creation: 
 *
 * Author:    Steven Eschrich
 *
 * Copyright: Copyright (C) 2007, Moffitt Cancer Center.
 *            All rights reserved.
 *
 * Update History
 * --------------
 * 04/08/05: Imported/repaired from old libaffy (AMH)
 * 03/14/08: New error handling scheme (AMH)
 * 04/15/08: Refactorization (AMH)
 * 06/04/08: Rename file to be more descriptive, other changes (AMH)
 * 09/20/10: Pooled memory allocator (AMH)
 * 10/26/10: Fixed missing warn %s arg in affy_load_chipset_single() (EAW)
 * 05/22/19: Added --ignore-chip-mismatch support (EAW)
 *
 **************************************************************************/

#include <affy.h>
#include <utils.h>

void affy_load_chipset_single(AFFY_CHIPSET *cs,
                              char *pathname,
                              bool ignore_chip_mismatch,
                              AFFY_ERROR *err)
{
  char *chip_type;

  assert(cs       != NULL);
  assert(pathname != NULL);

  assert(cs->num_chips <= cs->max_chips);

  if (cs->num_chips == cs->max_chips)
    AFFY_HANDLE_ERROR_VOID("chipset is full", AFFY_ERROR_LIMITREACHED, err);
  
  /* 
   * Check the array type against the type recorded for this
   * chipset.
   */
  chip_type = affy_get_cdf_name_from_cel(pathname, err);
  AFFY_CHECK_ERROR_VOID(err);
  
  assert(cs->array_type != NULL);
  if (strcmp(chip_type, cs->array_type) != 0 &&
      ignore_chip_mismatch == 0)
  {
    warn("Array type mismatch for CEL file %s.  Expected %s, "
         "found %s", 
         pathname,
         cs->array_type, 
         chip_type);

    AFFY_HANDLE_ERROR_GOTO("CEL file array type does not match chipset", 
                           AFFY_ERROR_WRONGTYPE, 
                           err, 
                           done);
  }
  else
  {
    /* Everything is in order, attempt to load the chip. */
    assert(cs->chip != NULL);
    cs->chip[cs->num_chips] = affy_load_chip(pathname, err);
    if (err->type == AFFY_ERROR_NONE)
    {
      hattach(cs->chip[cs->num_chips], cs->chip);
      cs->chip[cs->num_chips]->cdf = cs->cdf;
      cs->num_chips++;
    }
  }

done:
  h_free(chip_type);
}

void affy_load_chipset(AFFY_CHIPSET *cs, char **filelist,
                       bool ignore_chip_mismatch)
{
  char     **p;
  AFFY_ERROR err;

  assert(cs       != NULL);
  assert(filelist != NULL);

  /* Load each chip. */
  for (p = filelist; (*p != NULL) && (cs->num_chips < cs->max_chips); p++)
    affy_load_chipset_single(cs, *p, ignore_chip_mismatch, &err);
}
