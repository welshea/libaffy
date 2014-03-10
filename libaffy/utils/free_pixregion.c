
/**************************************************************************
 *
 * Filename:  free_pixregion.c
 *
 * Purpose:   Free an AFFY_PIXREGION.
 *
 * Creation:  20 Sep 2010
 *
 * Author:    Andrew M. Hoerter
 *
 * Copyright: Copyright (C) 2010, Moffitt Cancer Center.
 *            All rights reserved.
 *
 * Update History
 * --------------
 * 09/20/10: Creation (AMH)
 * 03/10/14: #ifdef out CEL qc fields to save memory (EAW)
 *
 **************************************************************************/

#include <affy.h>

#ifdef STORE_CEL_QC
/* 
 * affy_free_pixregion(): Free an AFFY_PIXREGION.
 *
 * Inputs: p_reg points to an AFFY_PIXREGION
 *
 * Outputs: None.
 *
 * Side effects: Associated memory is freed and the parent AFFY_CELL 
 *               pointer to this pixregion (if any) is zeroed.
 *
 */
void affy_free_pixregion(AFFY_PIXREGION *p_reg)
{
  if (p_reg && p_reg->cell)
    p_reg->cell->pixels = NULL;

  h_free(p_reg);
}
#endif
