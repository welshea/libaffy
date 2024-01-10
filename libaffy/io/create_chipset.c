
/**************************************************************************
 *
 * Filename:  create_chipset.c
 *
 * Purpose:   Initialize an AFFY_CHIPSET.
 *
 * Creation:  04/11/08
 *
 * Author:    Andrew M. Hoerter
 *
 * Copyright: Copyright (C) 2008, Moffitt Cancer Center.
 *            All rights reserved.
 *
 * Update History
 * --------------
 * 04/11/08: Creation (AMH)
 * 06/02/08: Refactor things a bit (AMH)
 * 09/20/10: Pooled memory allocator (AMH)
 * 11/16/10: Added reuse median polish parameters (EAW)
 * 04/07/11: Added generic chipset/cdf functions (EAW)
 * 03/10/14: #ifdef out cdf->xy_ref code to save memory (EAW)
 * 02/25/15: set cdf->no_mm_flag in create_blank_generic_cdf() (EAW)
 * 03/15/19: generic cdf: allocate and initialize cdf->seen_xy (EAW)
 * 08/13/19: fix AFFY_HANDLE_ERROR_VOID where AFFY_HANDLE_ERROR needed (EAW)
 * 08/12/20: pass flags to affy_create_chipset() (EAW)
 * 01/10/24: swap create_blank_generic_cdf() numrows/numcols (EAW)
 *
 **************************************************************************/

#include <affy.h>

/*
 * affy_create_chipset()
 *
 * Creates a chipset with the provided maximum size and chip type.
 * `cdf_hint', if non-NULL, is a hint where to find the CDF file.
 *
 * Creating a CHIPSET is required before calling any other
 * CHIPSET-related functions.
 *
 * Note that the CDF file must exist at the time of CHIPSET creation,
 * or allocation will fail.
 */
AFFY_CHIPSET *affy_create_chipset(unsigned int max_chips,
                                  char *chip_type,
                                  char *cdf_hint,
                                  AFFY_COMBINED_FLAGS *f,
                                  AFFY_ERROR *err)
{
  AFFY_CHIPSET *cs;

  assert(chip_type != NULL);

  cs = h_calloc(1, sizeof(AFFY_CHIPSET));
  if (cs == NULL)
    AFFY_HANDLE_ERROR("calloc failed", AFFY_ERROR_OUTOFMEM, err, NULL);

  /* Initialize the chipset. */
  cs->cdf       = NULL;
  cs->max_chips = max_chips;
  cs->num_chips = 0;

  cs->chip = h_subcalloc(cs, max_chips, sizeof(AFFY_CHIP *));
  if (cs->chip == NULL)
    AFFY_HANDLE_ERROR_GOTO("calloc failed", AFFY_ERROR_OUTOFMEM, err, cleanup);

  cs->cdf = affy_load_cdf_file(chip_type, cdf_hint, f, err);
  AFFY_CHECK_ERROR_GOTO(err, cleanup);
  hattach(cs->cdf, cs);

  cs->numrows    = cs->cdf->numrows;
  cs->numcols    = cs->cdf->numcols;
  cs->array_type = cs->cdf->array_type;
  
  cs->affinities = NULL;
  cs->t_values   = NULL;
  cs->mp_allocated_flag = 0;
  cs->mp_populated_flag = 0;

  return (cs);

cleanup:
  affy_free_chipset(cs);

  return (NULL);
}


/*
 * Create a blank generic CDF structure.
 * There are no MM probes.
 */
AFFY_CDFFILE *create_blank_generic_cdf(unsigned int max_chips,
                                       unsigned int numprobes,
                                       AFFY_ERROR *err)
{
  AFFY_CDFFILE *cdf = NULL;
  int i;

  cdf = h_calloc(1, sizeof(AFFY_CDFFILE));
  if (cdf == NULL)
    AFFY_HANDLE_ERROR_GOTO("calloc failed", AFFY_ERROR_OUTOFMEM, err, cleanup);

  cdf->array_type = h_strdup("generic");
  if (cdf->array_type == NULL)
  {
    AFFY_HANDLE_ERROR_GOTO("strdup failed", 
                           AFFY_ERROR_OUTOFMEM,
                           err,
                           cleanup);
  }
  else
    hattach(cdf->array_type, cdf);

  
  cdf->numrows      = numprobes;   /* y is rows, not cols */
  cdf->numcols      = 1;           /* x is cols, not rows */
  cdf->numprobes    = numprobes;
  cdf->numprobesets = numprobes;
  cdf->numqcunits   = 0;
  cdf->no_mm_flag   = 1;

  cdf->cell_type    = NULL;
#ifdef STORE_XY_REF
  cdf->xy_ref       = NULL;
#endif
  cdf->probeset     = NULL;
  cdf->probe        = NULL;

/*    affy_load_text_cdf_file(fp, cdf, &pbs, err); */

  cdf->cell_type = h_subcalloc(cdf, cdf->numrows, sizeof(affy_uint8 *));
  if (cdf->cell_type == NULL)
    AFFY_HANDLE_ERROR("calloc failed", AFFY_ERROR_OUTOFMEM, err, NULL);

  cdf->cell_type[0] = h_subcalloc(cdf->cell_type, 
                                  cdf->numrows * cdf->numcols,
                                  sizeof(affy_uint8));
  if (cdf->cell_type[0] == NULL)
    AFFY_HANDLE_ERROR("calloc failed", AFFY_ERROR_OUTOFMEM, err, NULL);

  for (i = 1; i < cdf->numrows; i++)
    cdf->cell_type[i] = cdf->cell_type[i-1] + cdf->numcols;

#ifdef STORE_XY_REF
  cdf->xy_ref = h_subcalloc(cdf, cdf->numrows, sizeof(AFFY_PROBE **));
  if (cdf->xy_ref == NULL)
    AFFY_HANDLE_ERROR("calloc failed", AFFY_ERROR_OUTOFMEM, err, NULL);
  
  cdf->xy_ref[0] = h_subcalloc(cdf->xy_ref,
                               cdf->numrows * cdf->numcols,
                               sizeof(AFFY_PROBE *));
  if (cdf->xy_ref[0] == NULL)
    AFFY_HANDLE_ERROR("calloc failed", AFFY_ERROR_OUTOFMEM, err, NULL);

  for (i = 1; i < cdf->numrows; i++)
    cdf->xy_ref[i] = cdf->xy_ref[i-1] + cdf->numcols;
#endif

  cdf->probeset = h_subcalloc(cdf, cdf->numprobesets, sizeof(AFFY_PROBESET));
  if (cdf->probeset == NULL)
    AFFY_HANDLE_ERROR("calloc failed", AFFY_ERROR_OUTOFMEM, err, NULL);

  cdf->seen_xy = h_subcalloc(cdf, cdf->numcols, sizeof(affy_uint8 **));
  if (cdf->seen_xy == NULL)
    AFFY_HANDLE_ERROR("calloc failed", AFFY_ERROR_OUTOFMEM, err, NULL);
  
  cdf->seen_xy[0] = h_subcalloc(cdf->seen_xy, cdf->numrows * cdf->numcols,
                               sizeof(affy_uint8 *));
  if (cdf->seen_xy[0] == NULL)
    AFFY_HANDLE_ERROR("calloc failed", AFFY_ERROR_OUTOFMEM, err, NULL);

  for (i = 1; i < cdf->numcols; i++)
    cdf->seen_xy[i] = cdf->seen_xy[i-1] + cdf->numrows;

  cdf->probe = h_subcalloc(cdf, numprobes, sizeof(AFFY_PROBE *));
  if (cdf->probe == NULL)
    AFFY_HANDLE_ERROR("calloc failed", AFFY_ERROR_OUTOFMEM, err, NULL);


  for (i = 0; i < numprobes; i++)
  {
    /* Allocate enough storage for all probes in this probe set */
    cdf->probeset[i].probe = h_subcalloc(cdf, 1, sizeof(AFFY_PROBE));
    if (cdf->probeset[i].probe == NULL)
      AFFY_HANDLE_ERROR("calloc failed", AFFY_ERROR_OUTOFMEM, err, NULL);

    /* Name will be filled in elsewhere from this function */
    cdf->probeset[i].name      = NULL;

    cdf->probeset[i].numprobes = 1;
    cdf->probeset[i].index     = i;
    
    /* There are no MM probes, so set MM same as PM */
    cdf->probeset[i].probe[0].pm.x = 0;
    cdf->probeset[i].probe[0].pm.y = i;
    cdf->probeset[i].probe[0].mm.x = 0;
    cdf->probeset[i].probe[0].mm.y = i;
    
    /* store cell type (normal) */
    cdf->cell_type[0][i] = AFFY_NORMAL_LOCATION;

    /* point to probe */
#ifdef STORE_XY_REF
    cdf->xy_ref[0][i]    = &(cdf->probeset[i].probe[0]);
#endif
    cdf->probe[i]        = &(cdf->probeset[i].probe[0]);
    cdf->probe[i]->index = i;
    cdf->probe[i]->ps    = &(cdf->probeset[i]);
  }

cleanup:

  if (err->type != AFFY_ERROR_NONE)
  {
    affy_free_cdf_file(cdf);

    return (NULL);
  }
  
  return (cdf);
}



/*
 * create_blank_generic_chipset()
 *
 * Creates a generic chipset with the provided maximum size.
 * There are no MM probes.
 *
 * Creating a CHIPSET is required before calling any other
 * CHIPSET-related functions.
 */
AFFY_CHIPSET *create_blank_generic_chipset(unsigned int max_chips,
                                           unsigned int numprobes,
                                           AFFY_ERROR *err)
{
  AFFY_CHIPSET *cs;

  cs = h_calloc(1, sizeof(AFFY_CHIPSET));
  if (cs == NULL)
    AFFY_HANDLE_ERROR("calloc failed", AFFY_ERROR_OUTOFMEM, err, NULL);

  /* Initialize the chipset. */
  cs->cdf       = NULL;
  cs->max_chips = max_chips;
  cs->num_chips = 0;

  cs->chip = h_subcalloc(cs, max_chips, sizeof(AFFY_CHIP *));
  if (cs->chip == NULL)
    AFFY_HANDLE_ERROR_GOTO("calloc failed", AFFY_ERROR_OUTOFMEM, err, cleanup);

  cs->cdf = create_blank_generic_cdf(max_chips, numprobes, err);
  AFFY_CHECK_ERROR_GOTO(err, cleanup);
  hattach(cs->cdf, cs);

  cs->numrows    = cs->cdf->numrows;
  cs->numcols    = cs->cdf->numcols;
  cs->array_type = cs->cdf->array_type;
  
  cs->affinities = NULL;
  cs->t_values   = NULL;
  cs->mp_allocated_flag = 0;
  cs->mp_populated_flag = 0;

  return (cs);

cleanup:
  affy_free_chipset(cs);

  return (NULL);
}
