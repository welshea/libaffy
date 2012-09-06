/**************************************************************************
 *
 * Filename:  chipset.c
 *
 * Purpose:   Operations on an AFFY_CHIPSET.
 *
 * Creation:  1/8/2009
 *
 * Author:    Andrew M. Hoerter
 *
 * Copyright: Copyright (C) 2008, Moffitt Cancer Center.
 *            All rights reserved.
 *
 * Update History
 * --------------
 * 01/08/09: Creation (SAE)
 * 09/20/10: Pooled memory allocator (AMH)
 * 04/11/11: Added affy_clone_chipset_one_chip() (EAW)
 *
 **************************************************************************/

#include <affy.h>

/**
 * Create a shallow copy of the chipset. The CDF is shared with the
 * clone, however the chip * is not. Both chipsets point to the
 * same chips initially.  It is important to note that it's entirely
 * ok to h_free() the clone chipset, it won't affect the original one
 * (since all the sub-structures are attached to the original CS pointer,
 * not the clone pointer).  However, free'ing the original chipset will
 * invalidate the clone.
 */
AFFY_CHIPSET *affy_clone_chipset(AFFY_CHIPSET *cur_chip, 
                                 AFFY_ERROR *err)
{
  AFFY_CHIPSET *cs;
  int           i;

  assert(cur_chip != NULL);

  cs = h_malloc(sizeof(AFFY_CHIPSET));
  if (cs == NULL)
    AFFY_HANDLE_ERROR("malloc failed", AFFY_ERROR_OUTOFMEM, err, NULL);

  /* Initialize the chipset */
  cs->cdf       = cur_chip->cdf;
  cs->max_chips = cur_chip->max_chips;
  cs->num_chips = cur_chip->num_chips;

  cs->chip = h_subcalloc(cs, cs->max_chips, sizeof(AFFY_CHIP *));
  if (cs->chip == NULL )
    AFFY_HANDLE_ERROR_GOTO("calloc failed", AFFY_ERROR_OUTOFMEM, err, cleanup);
  for (i = 0; i < cs->max_chips; i++) 
    cs->chip[i] = cur_chip->chip[i];

  cs->numrows    = cur_chip->numrows;
  cs->numcols    = cur_chip->numcols;
  cs->array_type = cur_chip->array_type;

  return (cs);

 cleanup:
  h_free(cs);

  return (NULL);
}

/**
 * Create a shallow copy of the chipset. The CDF is shared with the
 * clone, however the chip * is not. Both chipsets point to the
 * same chips initially.  It is important to note that it's entirely
 * ok to h_free() the clone chipset, it won't affect the original one
 * (since all the sub-structures are attached to the original CS pointer,
 * not the clone pointer).  However, free'ing the original chipset will
 * invalidate the clone.
 *
 * Only a single selected chip (chip_index) is copied into the new chipset.
 */
AFFY_CHIPSET *affy_clone_chipset_one_chip(AFFY_CHIPSET *cur_chip,
                                          int chip_idx,
                                          AFFY_ERROR *err)
{
  AFFY_CHIPSET *cs;

  assert(cur_chip != NULL);
  assert(chip_idx < cur_chip->max_chips);

  cs = h_malloc(sizeof(AFFY_CHIPSET));
  if (cs == NULL)
    AFFY_HANDLE_ERROR("malloc failed", AFFY_ERROR_OUTOFMEM, err, NULL);

  /* Initialize the chipset */
  cs->cdf       = cur_chip->cdf;
  cs->max_chips = 1;
  cs->num_chips = 1;

  cs->chip = h_subcalloc(cs, 1, sizeof(AFFY_CHIP *));
  if (cs->chip == NULL )
    AFFY_HANDLE_ERROR_GOTO("calloc failed", AFFY_ERROR_OUTOFMEM, err, cleanup);
  cs->chip[0] = cur_chip->chip[chip_idx];

  cs->numrows    = cur_chip->numrows;
  cs->numcols    = cur_chip->numcols;
  cs->array_type = cur_chip->array_type;

  return (cs);

 cleanup:
  h_free(cs);

  return (NULL);
}

/**
 * Resize a chipset to the desired number of chips, keeping as
 * many of the pointers to original chips as fit in the new max_chips.
 */
AFFY_CHIPSET *affy_resize_chipset(AFFY_CHIPSET *cs, 
                                  unsigned int max_chips, 
                                  AFFY_ERROR *err)
{
  AFFY_CHIP **ac;
  int         i;

  assert(cs != NULL);

  ac = h_subcalloc(cs, max_chips, sizeof(AFFY_CHIP *));
  if (ac == NULL)
    AFFY_HANDLE_ERROR("calloc failed", AFFY_ERROR_OUTOFMEM, err, NULL);
    
  for (i = 0; i < cs->max_chips && i < max_chips; i++) 
  {
    ac[i] = cs->chip[i];
    if (ac[i] != NULL)
      hattach(cs->chip[i], ac);
  }

  h_free(cs->chip);

  cs->chip = ac;
  cs->max_chips = max_chips;
  cs->num_chips = (max_chips < cs->num_chips) ? max_chips : cs->num_chips;

  return (cs);
}
