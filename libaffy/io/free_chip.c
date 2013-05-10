
/**************************************************************************
 *
 * Filename:  free_chip.c
 * 
 * Purpose:   Free routines for AFFY_CHIPs and AFFY_CHIPSETs.
 *
 * Creation:  05/27/05
 *
 * Author:    Andrew M. Hoerter
 *
 * Copyright: Copyright (C) 2007, Moffitt Cancer Center.
 *            All rights reserved.
 *
 * Update History
 * --------------
 * 05/27/05: Creation (AMH)
 * 09/20/10: Pooled memory allocator (AMH)
 * 05/10/13: Added affy_mostly_free_chip() (EAW)
 *
 **************************************************************************/

#include <affy.h>

/* 
 * affy_free_chip():
 *
 * Inputs: *ch contains an initialized AFFY_CHIP structure.
 * Outputs: None.
 * Side effects: Associated memory is freed.
 */
void affy_free_chip(AFFY_CHIP *ch)
{
  h_free(ch);
}

/* 
 * affy_free_chip():
 *
 * Inputs: *ch contains an initialized AFFY_CHIP structure.
 * Outputs: None.
 * Side effects: Associated memory is freed.
 *
 * Only *mostly* frees the data structure.  Keeps the informative stuff from
 *  the cel data structure, as well as the filename and cdf pointer.
 */
void affy_mostly_free_chip(AFFY_CHIP *ch)
{
  if (ch->cel)                   affy_mostly_free_cel_file(ch->cel);
  if (ch->dat)                   h_free(ch->dat);
  if (ch->probe_set)             h_free(ch->probe_set);
  if (ch->probe_set_call_pvalue) h_free(ch->probe_set_call_pvalue);
  if (ch->pm)                    h_free(ch->pm);
  
  ch->dat                        = NULL;
  ch->probe_set                  = NULL;
  ch->probe_set_call_pvalue      = NULL;
  ch->pm                         = NULL;
}


/* 
 * affy_free_chipset():
 *
 * Inputs: *ch contains an initialized AFFY_CHIPSET structure.
 * Outputs: None.
 * Side effects: Associated memory is freed.
 */
void affy_free_chipset(AFFY_CHIPSET *cs)
{
  h_free(cs);
}
