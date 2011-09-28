
/**************************************************************************
 *
 * Filename:  affy_iron.h
 *
 * Purpose:   IRON declarations.
 *
 * Creation:
 *
 * Author:    Eric A. Welsh
 *
 * Copyright: Copyright (C) 2011, Moffitt Cancer Center.
 *            All rights reserved.
 *
 * Update History
 * --------------
 * 9/28/11: file creation (EAW)
 *
 **************************************************************************/

#ifndef _AFFY_IRON_H_
#define _AFFY_IRON_H_

#include <math.h>
#include <utils.h>

#include "affy.h"

AFFY_CHIPSET *affy_iron(char **filelist, AFFY_COMBINED_FLAGS *f, AFFY_ERROR *err);
