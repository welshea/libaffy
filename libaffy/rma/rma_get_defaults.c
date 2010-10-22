
/**************************************************************************
 *
 * Filename: rma_get_defaults.c
 *
 * Purpose:  Allocate a new RMA flags structure and initialize it with
 *           appropriate defaults.
 *
 * Creation: 
 *
 * Author:   Steven Eschrich
 *
 *
 * Update History
 * -------------
 * 04/14/05: Imported/repaired from old libaffy (AMH)
 * 06/05/08: New error handling scheme (AMH)
 * 09/20/10: Pooled memory allocator (AMH)
 * 10/22/10: Use new AFFY_COMBINED_FLAGS instead of AFFY_RMA_FLAGS
 *
 **************************************************************************/

#include <affy_rma.h>

AFFY_COMBINED_FLAGS *affy_rma_get_defaults(AFFY_ERROR *err)
{
  AFFY_COMBINED_FLAGS *f = h_malloc(sizeof(AFFY_COMBINED_FLAGS));

  if (f == NULL)
    AFFY_HANDLE_ERROR("malloc failed", AFFY_ERROR_OUTOFMEM, err, NULL);

  affy_rma_set_defaults(f);

  return (f);
}
