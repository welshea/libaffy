/**************************************************************************
 *
 * Filename:  chip.c
 *
 * Purpose:   Operations on an AFFY_CHIP.
 *
 * Creation:  4/11/2011
 *
 * Author:    Eric A. Welsh
 *
 * Copyright: Copyright (C) 2011, Moffitt Cancer Center.
 *            All rights reserved.
 *
 * Update History
 * --------------
 * 04/11/11: Creation (EAW)
 *
 **************************************************************************/

#include <affy.h>

/**
 * Create a shallow copy of the chipset. The CDF is shared with the
 * clone, however the CEL * is not. Both chips point to the same CELs
 * initially.  It is important to note that it's entirely ok to h_free() the
 * clone chip, it won't affect the original one (since all the sub-structures
 * are attached to the original chip pointer, not the clone pointer).  However,
 * free'ing the original chipset will invalidate the clone.
 */
AFFY_CHIP *affy_clone_chip(AFFY_CHIP *cur_chip, AFFY_ERROR *err)
{
  AFFY_CELFILE *cur_cel = cur_chip->cel;
  AFFY_CHIP    *chip = NULL;
  AFFY_CELFILE *cf   = NULL;
  int           j, nbytes;

  assert(cur_chip != NULL);

  /* allocate chip */
  chip = h_calloc(1, sizeof(AFFY_CHIP));
  if (chip == NULL)
  {
    AFFY_HANDLE_ERROR_GOTO("calloc failed", AFFY_ERROR_OUTOFMEM, err, cleanup);
  }
      
  /* allocate CEL */
  cf = h_calloc(1, sizeof(AFFY_CELFILE));
  if (cf == NULL)
    AFFY_HANDLE_ERROR_GOTO("calloc failed", AFFY_ERROR_OUTOFMEM, err, cleanup);

  cf->filename = h_strdup(cur_cel->filename);
  if (cf->filename == NULL)
  {
    AFFY_HANDLE_ERROR_GOTO("strdup failed", AFFY_ERROR_OUTOFMEM, err, cleanup);
  }
  hattach(cf->filename, cf);
      
  cf->numrows     = cur_cel->numrows;
  cf->numcols     = cur_cel->numcols;
  cf->nummasks    = cur_cel->nummasks;
  cf->numoutliers = cur_cel->numoutliers;
  cf->mask        = NULL;
  cf->outlier     = NULL;

  nbytes = numbytes(cf->numcols);

  cf->data = h_subcalloc(cf, cf->numrows, sizeof(AFFY_CELL *));
  if (cf->data == NULL)
    AFFY_HANDLE_ERROR_GOTO("calloc failed", AFFY_ERROR_OUTOFMEM, err, cleanup);

  cf->data[0] = h_subcalloc(cf->data, 
                            cf->numcols * cf->numrows, 
                            sizeof(AFFY_CELL));
  if (cf->data[0] == NULL)
    AFFY_HANDLE_ERROR_GOTO("calloc failed", AFFY_ERROR_OUTOFMEM, err, cleanup);

  cf->mask = h_subcalloc(cf, cf->numrows, sizeof(bitstr_t *));
  if (cf->mask == NULL)
    AFFY_HANDLE_ERROR_GOTO("calloc failed", AFFY_ERROR_OUTOFMEM, err, cleanup);

  cf->mask[0] = h_subcalloc(cf->mask, nbytes * cf->numrows, sizeof(bitstr_t));
  if (cf->mask[0] == NULL)
    AFFY_HANDLE_ERROR_GOTO("calloc failed", AFFY_ERROR_OUTOFMEM, err, cleanup);

  cf->outlier = h_subcalloc(cf, cf->numrows, sizeof(bitstr_t *));
  if (cf->outlier == NULL)
    AFFY_HANDLE_ERROR_GOTO("calloc failed", AFFY_ERROR_OUTOFMEM, err, cleanup);

  cf->outlier[0] = h_subcalloc(cf->outlier, 
                               nbytes * cf->numrows, 
                               sizeof(bitstr_t));
  if (cf->outlier[0] == NULL)
    AFFY_HANDLE_ERROR_GOTO("calloc failed", AFFY_ERROR_OUTOFMEM, err, cleanup);
  
  for (j = 1; j < cf->numrows; j++)
  {
    cf->data[j]    = cf->data[0] + (j * cf->numcols);
    cf->mask[j]    = cf->mask[0] + (j * nbytes);
    cf->outlier[j] = cf->outlier[0] + (j * nbytes);
  }
  
  /* copy cel arrays */
  memcpy(cf->data[0], cur_cel->data[0],
         cf->numcols * cf->numrows * sizeof(AFFY_CELL));
  memcpy(cf->mask[0], cur_cel->mask[0],
         nbytes * cf->numrows * sizeof(bitstr_t));
  memcpy(cf->outlier[0], cur_cel->outlier[0],
         nbytes * cf->numrows * sizeof(bitstr_t));

  /* initialize/point chip stuff */
  chip->cdf                   = cur_chip->cdf;
  chip->cel                   = cf;
  chip->filename              = NULL;
  chip->dat                   = NULL;
  chip->probe_set             = NULL;
  chip->probe_set_call_pvalue = NULL;
  chip->pm                    = NULL;

  hattach(chip->cel, chip);

  chip->filename = h_strdup(cur_chip->filename);
  if (chip->filename == NULL)
  {
    affy_free_chip(chip);

    AFFY_HANDLE_ERROR_GOTO("strdup failed", AFFY_ERROR_OUTOFMEM, err, cleanup);
  }
  hattach(chip->filename, chip);

  return (chip);

 cleanup:
  h_free(chip);

  return (NULL);
}
