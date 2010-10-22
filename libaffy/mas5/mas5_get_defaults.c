
/**************************************************************************
 *
 * Filename:  mas5_get_defaults.c
 *
 * Purpose:   Initialize a fresh MAS5 flags structure with appropriate
 *            defaults and return it.
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
 * 04/14/05: Imported/repaired from old libaffy (AMH)
 * 03/07:08: New error handling scheme (AMH)
 * 09/20/10: Pooled memory allocator (AMH)
 * 10/22/10: Use new AFFY_COMBINED_FLAGS instead of AFFY_MAS5_FLAGS
 *
 **************************************************************************/

#include <affy_mas5.h>

AFFY_COMBINED_FLAGS *affy_mas5_get_defaults(AFFY_ERROR *err)
{
  AFFY_COMBINED_FLAGS *f;

  f = h_malloc(sizeof(AFFY_COMBINED_FLAGS));

  if (f == NULL)
    AFFY_HANDLE_ERROR("malloc failed", AFFY_ERROR_OUTOFMEM, err, NULL);

  affy_mas5_set_defaults(f);

  return (f);
}
