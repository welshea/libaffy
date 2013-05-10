
/**************************************************************************
 *
 * Filename:  load_chip.c
 *
 * Purpose:   Initialize an AFFY_CHIP structure from disk.
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
 * 09/20/10: Pooled memory allocator (AMH)
 *
 **************************************************************************/

#include <affy.h>
#include <utils.h>

AFFY_CHIP *affy_load_chip(char *filename, AFFY_ERROR *err)
{
  AFFY_CELFILE *c;
  AFFY_CHIP    *chip;

  assert(filename != NULL);

  /* First try and open the file, if not quit now */
  c = affy_load_cel_file(filename, err);
  AFFY_CHECK_ERROR(err, NULL);

  /* OK, allocate storage */
  chip = h_malloc(sizeof(AFFY_CHIP));
  if (chip == NULL)
  {
    affy_free_cel_file(c);

    AFFY_HANDLE_ERROR("malloc failed", AFFY_ERROR_OUTOFMEM, err, NULL);
  }

  chip->cdf                   = NULL;
  chip->cel                   = c;
  chip->dat                   = NULL;
  chip->filename              = NULL;
  chip->probe_set             = NULL;
  chip->probe_set_call_pvalue = NULL;
  chip->pm                    = NULL;

  hattach(chip->cel, chip);

  chip->filename = h_strdup(c->filename);
  if (chip->filename == NULL)
  {
    affy_free_chip(chip);

    AFFY_HANDLE_ERROR("strdup failed", AFFY_ERROR_OUTOFMEM, err, NULL);
  }

  hattach(chip->filename, chip);

  return (chip);
}
