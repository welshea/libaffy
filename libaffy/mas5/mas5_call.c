
/**************************************************************************
 *
 * Filename:  mas5_call.c
 *
 * Purpose:   assign P/M/A calls
 *
 * Creation:  01/10/08
 *
 * Author:    Eric A. Welsh
 *
 *
 * Update History
 * --------------
 * 09/16/10: initial version (EAW)
 * 09/20/10: Pooled memory allocator (AMH)
 * 10/22/10: Use new AFFY_COMBINED_FLAGS instead of AFFY_MAS5_FLAGS
 * 06/09/11: modified algorithm to be closer to Affymetrix whitepaper
 * 03/06/14: Do not make calls if chip is missing MM probes (EAW)
 *
 **************************************************************************/

#include "affy_mas5.h"
#include "affy_wilcox.h"

/* Important tuning constants */
static const double TAU = 0.015;
static const double ALPHA1 = 0.04, ALPHA2 = 0.06;

/* Private routines */
static double calculate_probeset_call(AFFY_CHIP *c, 
                                      int probeset_num, 
                                      AFFY_ERROR *err);

static double calculate_probeset_call(AFFY_CHIP *c, 
                                      int probeset_num, 
                                      AFFY_ERROR *err)
{
  AFFY_PROBESET *p;
  AFFY_CELL    **data;
  double        *pm, *mm, *r, pvalue;
  int            n, i, j, x, y, *mempool;
  int            non_masked_count = 0, saturated_count = 0;

  /* Some shortcuts */
  p    = &(c->cdf->probeset[probeset_num]);
  n    = p->numprobes;
  data = c->cel->data;

  mempool = h_malloc(sizeof(int));
  if (mempool == NULL)
    AFFY_HANDLE_ERROR("malloc failed", AFFY_ERROR_OUTOFMEM, err, -DBL_MIN);

  /* Assign pm/mm values for this chip (based on existing layout) */
  pm = h_subcalloc(mempool, n, sizeof(double));
  if (pm == NULL)
  {
    h_free(mempool);
    AFFY_HANDLE_ERROR("calloc failed", AFFY_ERROR_OUTOFMEM, err, -DBL_MIN);
  }

  mm = h_subcalloc(mempool, n, sizeof(double));
  if (mm == NULL)
  {
    h_free(mempool);
    AFFY_HANDLE_ERROR("calloc failed", AFFY_ERROR_OUTOFMEM, err, -DBL_MIN);
  }

  r = h_subcalloc(mempool, n, sizeof(double));
  if (r == NULL)
  {
    h_free(mempool);
    AFFY_HANDLE_ERROR("calloc failed", AFFY_ERROR_OUTOFMEM, err, -DBL_MIN);
  }

  /*
    Assign pm/mm values based on probes. Caveat: masked probes (either 
    in PM or MM) are ignored in the computation.
  */
  for (i = 0, j = 0; i < n; i++)
  {
    double pm_i, mm_i;

    x = p->probe[i].pm.x;
    y = p->probe[i].pm.y;
    if (affy_ismasked(c, x, y))
      continue;
    pm_i = data[x][y].value;

    x = p->probe[i].mm.x;
    y = p->probe[i].mm.y;
    if (affy_ismasked(c, x, y))
      continue;
    mm_i = data[x][y].value;
    
    non_masked_count++;
    
    /* skip saturated mismatches */
    if (mm_i >= 46000)
    {
      saturated_count++;
      continue;
    }
    
    /* skip probes that are within TAU of each other */
    if (fabs(pm_i - mm_i) <= TAU)
      continue;

    /* Calculate R discrimination scores, skip (r - TAO == 0) */
    r[j] = (pm_i - mm_i) / (pm_i + mm_i);
    if (r[j] - TAU == 0)
      continue;

    pm[j] = pm_i;
    mm[j] = mm_i;
    j++;
  }
  
  /* assign completely saturated probesets as present (0.0) */
  if (saturated_count == non_masked_count)
  {
    return 0.0;
  }

  /* if no probes left, call 'A' (0.5) */
  if (j == 0)
    return 0.5;

  
  /* In case the total number of probes is less than n */
  if (n > j)
  {
#ifdef DEBUG
    info("%s:Adjusting for %d probes instead of %d", p->probe_set_id, j, n);
#endif
    n = j;
  }

  pvalue = affy_mas5_calculate_call_pvalue(r, n, TAU, err);

  h_free(mempool);

  AFFY_CHECK_ERROR(err, -DBL_MIN);

  return (pvalue);
}

char affy_mas5_pvalue_call(double pvalue)
{
  assert(pvalue >= 0.0);

  if (pvalue < ALPHA1)
    return ('P');
  
  if (pvalue < ALPHA2)
    return ('M');

  return ('A');
}

int affy_mas5_call(AFFY_CHIPSET *c, AFFY_COMBINED_FLAGS *f, AFFY_ERROR *err)
{
  int               i, n, num_probesets;
  LIBUTILS_PB_STATE pbs;

  assert(c            != NULL);
  assert(c->cdf       != NULL);
  assert(c->num_chips > 0);
  assert(f            != NULL);

  /* chip is missing MM probes, abort */
  for (n = 0; n < c->num_chips; n++)
  {
    /* chip is missing MM probes, abort */
    if (c->chip[n]->cdf->no_mm_flag == 1)
      return 0;
  }

  num_probesets = c->cdf->numprobesets;
  pb_init(&pbs);

  for (n = 0; n < c->num_chips; n++)
  {
    pb_begin(&pbs, num_probesets, 
             "Calculating calls for chip using Affymetrix method");
    c->chip[n]->probe_set_call_pvalue = h_subcalloc(c->chip[n],
                                                    num_probesets, 
                                                    sizeof(double));
    if (c->chip[n]->probe_set_call_pvalue == NULL)
      AFFY_HANDLE_ERROR_GOTO("calloc failed", AFFY_ERROR_OUTOFMEM, err, err);

    c->chip[n]->numprobesets = num_probesets;

    for (i = 0; i < num_probesets; i++)
    {
      pb_tick(&pbs, 1, "");
      c->chip[n]->probe_set_call_pvalue[i] = 
        calculate_probeset_call(c->chip[n], i, err);
    }

    pb_finish(&pbs, "Finished present/absent calls");
  }

  pb_cleanup(&pbs);

  return (0);

err:
  return (-1);
}
