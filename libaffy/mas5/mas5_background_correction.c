
/**************************************************************************
 *
 * Filename:  mas5_background_correction.c
 *
 * Purpose:   Perform background correction for MAS5.
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
 * 04/19/05: AFFY_CELFILE now uses AFFY_CELL's (AMH)
 * 10/03/07: Perform minor beautification/cleanup (AMH)
 * 03/07/08: New error handling scheme (AMH)
 * 09/20/10: Pooled memory allocator (AMH)
 * 10/22/10: Use new AFFY_COMBINED_FLAGS instead of AFFY_MAS5_FLAGS
 *
 **************************************************************************/

#include <affy_mas5.h>

#define START_X  0
#define START_Y  1
#define LENGTH_X 2
#define LENGTH_Y 3

static int    K;       /* The number of rectangular zones on the chip */
static int    smooth;
static double NoiseFrac;
static bool   bioconductor_compatability;

static AFFY_POINT *center;
static double     *bZ;
static double     *nZ;
static int         dim, default_grid_y_length, default_grid_x_length;

static int    find_centers(int rows, int cols);
static void   output_statistics();
static int    estimate_zone_background(AFFY_CHIP *chip, AFFY_ERROR *err);
static int    calculate_background();
static int    background(int x, int y, double *b, double *n);
static double w_k(int x, int y, int k);
static int    zone_information(int k, int type);

/*
 * This operation is detailed in the Affy white papers. It takes
 * a CELFILE object, background corrects it and returns.
 */
int affy_mas5_background_correction(AFFY_CHIPSET *c, 
                                    AFFY_COMBINED_FLAGS *f, 
                                    AFFY_ERROR *err)
{
  int n, *mempool;
  LIBUTILS_PB_STATE pbs;

  pb_init(&pbs);
  if (f == NULL)
    f = affy_mas5_get_defaults(err);

  AFFY_CHECK_ERROR(err, -1);

  K         = f->K;
  smooth    = f->smooth;
  NoiseFrac = f->NoiseFrac;
  bioconductor_compatability = f->bioconductor_compatability;

  pb_begin(&pbs, c->num_chips+2,"Background correction using Affymetrix method.");

  /* Compute these dynamically since we don't know K a priori */
  dim = sqrt(K);
  default_grid_x_length = c->numcols / dim;
  default_grid_y_length = c->numrows / dim;

  mempool = h_malloc(sizeof(int));
  if (mempool == NULL)
    AFFY_HANDLE_ERROR("malloc failed", AFFY_ERROR_OUTOFMEM, err, -1);

  /* Allocate space for these variables now */
  center = h_subcalloc(mempool, K, sizeof(AFFY_POINT));
  if (center == NULL)
  {
    h_free(mempool);
    AFFY_HANDLE_ERROR("calloc failed", AFFY_ERROR_OUTOFMEM, err, -1);
  }

  bZ = h_subcalloc(mempool, K, sizeof(double));
  if (bZ == NULL)
  {
    h_free(mempool);
    AFFY_HANDLE_ERROR("calloc failed", AFFY_ERROR_OUTOFMEM, err, -1);
  }

  nZ = h_subcalloc(mempool, K, sizeof(double));
  if (nZ == NULL)
  {
    h_free(mempool);
    AFFY_HANDLE_ERROR("calloc failed", AFFY_ERROR_OUTOFMEM, err, -1);
  }

  /*
   * Three steps: find grid centers, estimate zone background, then
   * finally compute the background intensities. The background is not
   * only calculated but intensities are adjusted relative to this
   * background.
   */
  pb_tick(&pbs,1,"Finding centers...");
  find_centers(c->numrows, c->numcols);
  pb_tick(&pbs,1,"Estimating zone background and calculating background correction: ");

  for (n = 0; n < c->num_chips; n++)
  {
    pb_tick(&pbs,1,"");
    estimate_zone_background(c->chip[n], err);
    if (err->type != AFFY_ERROR_NONE)
    {
      h_free(mempool);

      return (-1);
    }

    calculate_background(c->chip[n]);
  }

  pb_finish(&pbs,"Finished initial MAS5 background correction");
  h_free(mempool);

  return (0);
}

/* A simple calculation of grid centers, based on a square grid of size K */
static int find_centers(int rows, int cols)
{
  int running_x_offset, running_y_offset;
  int k;
  int lengthx = zone_information(0, LENGTH_X); /* first parameter N/A */
  int lengthy = zone_information(0, LENGTH_Y);
  int midy    = lengthy / 2;
  int midx    = lengthx / 2;

  /* 
   * Fill in the center points based on K. Bioconductor starts at 0-based
   * indexing and rounds up, whereas results suggest that MAS5.0 starts
   * at 1.
   */
  if (bioconductor_compatability)
    running_x_offset = running_y_offset = 0;
  else 
    running_x_offset = running_y_offset = 1;
  
  for (k = 0; k < K; k++)
  {
    /* Increment center if we wrap around a row */
    if (running_x_offset >= cols)
    {
      running_x_offset = 0;
      running_y_offset += lengthy;
    }

    /* Calculation is halfway into next region */
    center[k].x = running_x_offset + midx;
    running_x_offset += lengthx;

    center[k].y = running_y_offset + midy;
  }

  return (0);
}

/*
 * A zone is bit weird but here is the deal, thanks to perls of
 * wisdom gleaned from GNU Affymetrix code:
 * A zone cannot have split PM/MM pairs. Probe pairs appear to be
 * setup as ODDV then EVEN, so that the ODD is always the
 * smaller number.
 *        
 * What this means to zones is that the fixed length grid is not
 * fixed, since it must be adjusted to include an extra line at 
 * the beginning (which is QC anyway) and must be an extra pair
 * short at the end (the extra one from the top plus the QC row 
 * at the bottom).
 *    
 *
 * Given a particular zone
 */
static int zone_information(int k, int type) 
{
  switch (type) 
  {
    case START_X:
      /* The X coordinates are constant, so this is a simple calculation */
      return ((k % dim) * default_grid_x_length);
    case START_Y:
      return ((k / dim) * default_grid_y_length); 
    case LENGTH_X:
      return (default_grid_x_length);
    case LENGTH_Y:
      return (default_grid_y_length);
  }

  return (-1); 
}

/* This should calculate the lower 2% cutoff and stddev */
static int estimate_zone_background(AFFY_CHIP *chip, AFFY_ERROR *err)
{
  int           k, i, x, y;
  int           num_in_bg, num_bgvals, total_vals;
  double       *bgvals;
  AFFY_CELFILE *cf = chip->cel;

  /* What is 2% of a zone. This is the upper bound. */
  num_in_bg = 0.02 * (default_grid_x_length * (default_grid_y_length+1));
  bgvals    = h_calloc(num_in_bg, sizeof(double));
  if (bgvals == NULL)
    AFFY_HANDLE_ERROR("calloc failed", AFFY_ERROR_OUTOFMEM, err, -1);

  /* For each zone, copy the values to a separate array to work on */
  for (k = 0; k < K; k++)
  {
    int starty  = zone_information(k, START_Y);
    int startx  = zone_information(k, START_X);
    int lengthy = zone_information(k, LENGTH_Y);
    int lengthx = zone_information(k, LENGTH_X);

    num_bgvals = 0;
    total_vals = 0;

    /* IMPORTANT: Zero out this area for each zone */ 
    for (i = 0; i < num_in_bg; i++) 
      bgvals[i] = -1;

    /* Accumulate lower 2% of region */
    for (y = 0; y < lengthy; y++)
    {
      for (x = 0; x < lengthx; x++)
      {
        int ry = starty + y;    /* Real y offset */
        int rx = startx + x;    /* Real x offset */

        /* Skip over masked cells and undefined cells and QC cells */
        if (affy_ismasked(chip, rx, ry) 
            || affy_isundefined(chip, rx, ry) 
            || affy_isqc(chip, rx, ry))
          continue;

        total_vals++;

        /* Insert into bgvals list, in sorted order */
        if (num_bgvals < num_in_bg)
        {
          /* Insert into list */
          for (i = num_bgvals;
               i > 0 && bgvals[i - 1] > cf->data[rx][ry].value; i--)
            bgvals[i] = bgvals[i - 1];
          bgvals[i] = cf->data[rx][ry].value;
          num_bgvals++;
        }
        else if (bgvals[num_in_bg-1] < cf->data[rx][ry].value)
        {
          /* Do nothing */
        }
        else
        {
          /* Insert into list */
          for (i = num_in_bg - 1;
               i > 0 && bgvals[i - 1] > cf->data[rx][ry].value; i--)
            bgvals[i] = bgvals[i - 1];
          
          bgvals[i] = cf->data[rx][ry].value;
        }
      }
    }

    /* mean of the lower 2%: ASSUMPTION - 2% of usable cells in zone */
    num_bgvals = (int)(0.02 * total_vals);
    bZ[k]      = 0;

    for (i = 0; i < num_bgvals; i++)
      bZ[k] += bgvals[i];

    bZ[k] /= num_bgvals;

    /* Standard deviation */
    nZ[k] = 0;
 
    for (i = 0; i < num_bgvals; i++)
      nZ[k] += (bgvals[i] - bZ[k]) * (bgvals[i] - bZ[k]);

    nZ[k] = sqrt(nZ[k] / (num_bgvals - 1));
  }
#ifdef DEBUG
  /* Output the various statistics that come from this background */
  output_statistics();
#endif
  
  h_free(bgvals);

  return (0);
}

static void output_statistics() 
{
  double totalBG = 0;
  double totalN  = 0;
  int    k;
  
  for (k = 0; k < K; k++) 
  {
    totalBG += bZ[k];
    totalN  += nZ[k];
  }

  totalBG /= K;
  totalN  /= K;

  for (k = 0; k < K; k++) 
    info("Background[%d]=%f\n", k, bZ[k]);

  info("Average BG is %f\n",totalBG);

  for (k = 0; k < K; k++) 
    info("Noise[%d]=%f\n", k, nZ[k]);

  info("Average N is %f\n", totalN);
}

/* A Byzantine calculation involving many substeps */
static int calculate_background(AFFY_CHIP *chip)
{
  double        b, n, I_prime;
  int           x, y, pinterval, progress;
  AFFY_CELFILE *cf = chip->cel;

  assert(cf != NULL);

  pinterval = progress = cf->numrows / 10;

  for (y = 0; y < cf->numrows; y++)
  {
    if (y == progress)
      progress += pinterval;

    for (x = 0; x < cf->numcols; x++)
    {
      /* Skip over masked cells and undefined cells and QC cells */
      if (affy_ismasked(chip, x, y) 
          || affy_isundefined(chip, x, y) 
          || affy_isqc(chip, x, y))
        continue;

      /* Calculate both the b and n values */
      background(x, y, &b, &n);

      I_prime = max_macro(cf->data[x][y].value, 0.5);
      cf->data[x][y].value = max_macro(I_prime - b, NoiseFrac * n);
    }
  }

  return (0);
}

/* This is the b(x,y) calculation */
static int background(int x, int y, double *b, double *n)
{
  double denom = 0, n_n = 0, b_n = 0;
  double cur_weight;
  int    k;

  /* 
   * Although C operates on 0-based indexing, both Bioconductor and
   *  MAS5.0 assume that the grid is 1-based indexing. We increment
   *  both x and y coordinates before calculating background. Note
   *  this does not address the intensity at (x,y) from the cel file,
   *  which is 0-based indexing.
   */
  x++;
  y++;

  for (k = 0; k < K; k++)
  {
    cur_weight = w_k(x, y, k);

    denom += cur_weight;
    b_n   += cur_weight * bZ[k];
    n_n   += cur_weight * nZ[k];
  }

  *b = b_n / denom;
  *n = n_n / denom;

  return (0);
}

/* This is the w_k(row,col) calculation */
static double w_k(int x, int y, int k)
{
  double d;

  /* 
   * From my understanding of the Bioconductor code, an additional 0.5 
   * is added to the computation of the centers, which since they use 
   * floating point won't be truncated. Here, we (optionally) add that
   *  additional value back in, since our centers are integers.
   */
  if (bioconductor_compatability) 
  {
    d = ((x - center[k].x - 0.5) * (x - center[k].x - 0.5)) 
        + ((y - center[k].y - 0.5) * (y - center[k].y - 0.5));
  }
  else 
  {
    d = ((x - center[k].x) * (x - center[k].x)) 
        + ((y - center[k].y) * (y - center[k].y));
  }

  return (1.0 / (d + smooth));
}
