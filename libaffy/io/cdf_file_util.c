
/**************************************************************************
 *
 * Filename:  cdf_file_util.c
 * 
 * Purpose:   CDF file utility routines.
 *
 * Creation:  03/27/05
 *
 * Author:    Andrew M. Hoerter
 *
 * Copyright: Copyright (C) 2007, Moffitt Cancer Center.
 *            All rights reserved.
 *
 * Update History
 * --------------
 * 05/27/05: Creation (AMH)
 * 03/11/08: New error handling scheme (AMH)
 * 09/20/10: Pooled memory allocator (AMH)
 *
 **************************************************************************/

#include <affy.h>

/* 
 * affy_free_cdf_file(): Free a CDF file structure and any associated
 * storage.
 *
 * Inputs: *cdf contains an initialized CDF file structure.
 * Outputs: None.
 * Side effects: Associated memory is freed.
 */
void affy_free_cdf_file(AFFY_CDFFILE *cdf)
{
  h_free(cdf);
}
