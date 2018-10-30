
/**************************************************************************
 *
 * Filename:  load_binary_cdf_file.c
 *
 * Purpose:   Parse a binary CDF file and initialize an accompanying structure.
 *
 * Creation:  11/10/05
 *
 * Author:    Steven Eschrich
 *
 * Copyright: Copyright (C) 2007, Moffitt Cancer Center.
 *            All rights reserved.
 *
 * Update History
 * --------------
 * 04/08/05: Imported/repaired from old libaffy (AMH)
 * 10/03/07: Minor beautification and cleanups (AMH)
 * 01/25/08: Convert to use readmulti() (AMH)
 * 03/11/08: New error handling scheme (AMH)
 * 09/20/10: Pooled memory allocator (AMH)
 * 03/06/14: fixed row/col memory allocation errors, the dimensions were swapped (EAW)
 * 03/10/14: #ifdef out cdf->xy_ref code to save memory (EAW)
 * 03/xx/14: around the same time, added support for no-MM CDF files (EAW)
 * 07/19/17: rewrote most of the binary CDF parser to fix numerous bugs(EAW)
 * 08/27/18: fixed some AFFY_ERROR return values that OSX clang caught (EAW)
 *
 **************************************************************************/

#include <affy.h>

static void process_qc_section(FILE *fp,
                               AFFY_CDFFILE *cdf,
                               LIBUTILS_PB_STATE *pbs,
                               AFFY_ERROR *err);

/* A probe section will be indexed by probe, and by x,y coordinates
 * Return values: 0,1 -- all_mm_flag status
 *               -1   -- did not read in any probesets
 */
static int process_probe_section(FILE *fp,
				 AFFY_CDFFILE *cdf,
				 int version,
				 int ps,
                                 LIBUTILS_PB_STATE *pbs,
				 AFFY_ERROR *err);

/*
 * The CDF file contains the description of the microarray chip in terms
 * of probes and probe sets. This information needs to be loaded into
 * memory so that calculations at the probe level can be matched
 * to (X,Y) coordinates in the CEL file.
 * 
 * The parameter is a chip_type as defined by the CEL file. This
 * should correspond to a cdf file somewhere.
 */
void affy_load_binary_cdf_file(FILE *fp, 
			       AFFY_CDFFILE *cdf,
                               LIBUTILS_PB_STATE *pbs,
			       AFFY_ERROR *err)
{
  affy_int32        version, numps;
  affy_int32        custom_len;
  affy_uint16       tmp_numcols, tmp_numrows;
  int               i, temp_int;
  char              ps_name[65];
  int               all_mm_flag = 1;

  assert(cdf != NULL);
  assert(fp  != NULL);

  /* 
   * *** Historic NOTE, up until 2017-07-19 ***
   * NOTE: it is assumed the magic number has already been read and
   * verified, and the file pointer is positioned after it.
   *
   * *** Updated NOTE, as of 2017-07-19 ***
   * NOTE: it appears that, at some time, affy_load_cdf_file() was changed
   *  to reopen fp, thus unreading the magic number, and causing the rest
   *  of this function to read everything in starting from the wrong file
   *  pointer position, causing lots of *bad things* to happen.
   *
   * Furthermore, a while ago, cdf->numcols and cdf->numrows were expanded
   *  to 32-bits, rather than 16-bits, to allow for various new things.
   *  This broke the readmulti() calls, which wrote to the structures
   *  directly, thus writing the wrong size data to the newly bigger variables,
   *  resulting in huge random garbage values in the high-order bits.
   *
   * Fixed this by reading to temp 16-bit variables first.
   *
   */

  /* skip past the magic number again, since we reopened the fp earlier */
  if (affy_readmulti(fp, "%x", 4) <= 0)
    AFFY_HANDLE_ERROR_VOID("I/O error, can't read CDF magic number",
                           AFFY_ERROR_IO,
                           err);
   
  if (affy_readmulti(fp, "%dl", (void *)&version) <= 0)
    AFFY_HANDLE_ERROR_VOID("I/O error in CDF header section: version number",
                           AFFY_ERROR_IO,
                           err);

  /* Report version */
  info("Found XDA (binary) CDF version %" AFFY_PRNd32, version);
  
  if (version == 4)
  {
    AFFY_HANDLE_ERROR_VOID("I/O error, binary CDF version 4 not supported",
                           AFFY_ERROR_IO,
                           err);
  }

  if (affy_readmulti(fp, "%2hl%3dl",
                     (void *)&tmp_numcols,
                     (void *)&tmp_numrows,
                     (void *)&numps,
                     (void *)&cdf->numqcunits,
                     (void *)&custom_len) <= 0)
    AFFY_HANDLE_ERROR_VOID("I/O error in CDF header section",
                           AFFY_ERROR_IO,
                           err);

  cdf->numcols = tmp_numcols;
  cdf->numrows = tmp_numrows;

#ifdef DEBUG
  fprintf(stderr, "NumCols: %d\n", (unsigned int) cdf->numcols);
  fprintf(stderr, "NumRows: %d\n", (unsigned int) cdf->numrows);
  fprintf(stderr, "NumPS:   %d\n", (unsigned int) numps);
  fprintf(stderr, "NumRows: %d\n", (unsigned int) cdf->numqcunits);
  fprintf(stderr, "CLen:    %d\n", (unsigned int) custom_len);
#endif

  /* skip past the CustomSeq reference sequence */
  if (custom_len && affy_readmulti(fp, "%x", custom_len) <= 0)
    AFFY_HANDLE_ERROR_VOID("I/O error, can't read past CustomSeq reference",
                           AFFY_ERROR_IO,
                           err);


  /* 
   * At this point we know the size of the chip, and the number
   * of probes.
   */
  cdf->cell_type = h_suballoc(cdf, cdf->numcols * sizeof(affy_uint8 *));
  if (cdf->cell_type == NULL)
    AFFY_HANDLE_ERROR_VOID("malloc failed",
                           AFFY_ERROR_OUTOFMEM,
                           err);

  cdf->cell_type[0] = h_subcalloc(cdf->cell_type,
                                  cdf->numrows * cdf->numcols, 
                                  sizeof(affy_uint8));
  if (cdf->cell_type[0] == NULL)
    AFFY_HANDLE_ERROR_VOID("calloc failed",
                           AFFY_ERROR_OUTOFMEM,
                           err);

#ifdef STORE_XY_REF
  cdf->xy_ref = h_suballoc(cdf, cdf->numcols * sizeof(AFFY_PROBE **));
  if (cdf->xy_ref == NULL)
    AFFY_HANDLE_ERROR_VOID("malloc failed",
                           AFFY_ERROR_OUTOFMEM,
                           err);
  
  cdf->xy_ref[0] = h_subcalloc(cdf->xy_ref,
                               cdf->numcols * cdf->numrows, 
                               sizeof(AFFY_PROBE *));
  if (cdf->xy_ref[0] == NULL)
    AFFY_HANDLE_ERROR_VOID("calloc failed",
                           AFFY_ERROR_OUTOFMEM,
                           err);
#endif

  for (i = 1; i < cdf->numcols; i++)
  {
    cdf->cell_type[i] = cdf->cell_type[i-1] + cdf->numrows;
#ifdef STORE_XY_REF
    cdf->xy_ref[i]    = cdf->xy_ref[i-1] + cdf->numrows;
#endif
  }

  cdf->probeset = h_subcalloc(cdf, numps, sizeof(AFFY_PROBESET));
  if (cdf->probeset == NULL)
    AFFY_HANDLE_ERROR_VOID("calloc failed",
                           AFFY_ERROR_OUTOFMEM,
                           err);

  cdf->seen_xy = h_subcalloc(cdf, cdf->numcols, sizeof(affy_uint8 **));
  if (cdf->seen_xy == NULL)
    AFFY_HANDLE_ERROR_VOID("calloc failed", AFFY_ERROR_OUTOFMEM, err);
  
  cdf->seen_xy[0] = h_subcalloc(cdf->seen_xy, cdf->numrows * cdf->numcols,
                               sizeof(affy_uint8 *));
  if (cdf->seen_xy[0] == NULL)
    AFFY_HANDLE_ERROR_VOID("calloc failed", AFFY_ERROR_OUTOFMEM, err);

  for (i = 1; i < cdf->numcols; i++)
    cdf->seen_xy[i] = cdf->seen_xy[i-1] + cdf->numrows;


  /* 
   * Since the number of probes cannot be known without first reading
   * in the entire CDF, we make a SWAG here to size the probe list.
   * Since each probe will have a PM/MM pair of cells (and some have a
   * quartet?), the number of cells divided by 2 provides an upper
   * limit for how many probes will exist.  Some space will likely go
   * unused.
   *
   * EAW -- Update: We can no longer assume PM/MM, now that Affy is making
   *                PM-only arrays.  So, don't divide by 2 for now, but
   *                realloc later to save space after CDF is read in.
   */
  cdf->probe = h_subcalloc(cdf,
                           (cdf->numrows * cdf->numcols),
                           sizeof(AFFY_PROBE *));
  if (cdf->probe == NULL)
    AFFY_HANDLE_ERROR_VOID("calloc failed",
                           AFFY_ERROR_OUTOFMEM,
                           err);
        
  /* Probeset names */
  /* These may not necessarily be guaranteed to be in the same order that
   * we read them in later?  So, don't store them here.
   */
  for (i = 0; i < numps; i++) 
  {
    if (affy_readmulti(fp, "%x", 64) <= 0)
      AFFY_HANDLE_ERROR_VOID("I/O error in CDF probeset name header section",
                             AFFY_ERROR_IO, err);

  }
  cdf->numprobesets = numps;
  
  /* Skip details about file positions */
  if (fseek(fp, cdf->numqcunits * 4 + cdf->numprobesets * 4, SEEK_CUR) != 0)
    AFFY_HANDLE_ERROR_VOID("couldn't seek within CDF file",
                           AFFY_ERROR_IO,
                           err);

  /* Process QC section */
  process_qc_section(fp, cdf, pbs, err);
  AFFY_CHECK_ERROR_VOID(err);

  /* Process Probe section */
  pb_begin(pbs, (cdf->numrows * cdf->numcols) / 2, "Loading probes");

  for (i = 0; i < cdf->numprobesets; i++) 
  {
    temp_int = process_probe_section(fp, cdf, version, i, pbs, err);
    AFFY_CHECK_ERROR_VOID(err);
    
    /* we read in some probes without a MM, so we can't realloc to save
     * memory
     */
    if (temp_int == 0)
      all_mm_flag = 0;
  }

  /* free extra memory that was pre-allocated for MM if we are missing some */
  if (all_mm_flag == 0)
  {
    cdf->probe = h_realloc(cdf->probe, cdf->numprobes * sizeof(AFFY_PROBE *));
  }
  
  /* flag CDF as missing MM probes */
  if (all_mm_flag == 0)
    cdf->no_mm_flag = 1;

  /* flag CDF as having duplicate probes */
  cdf->dupe_probes_flag = 0;
  for (i = 0; i < cdf->numrows * cdf->numcols; i++)
  {
    if (cdf->seen_xy[0][i] == 2)
    {
      cdf->dupe_probes_flag = 1;
      break;
    }
  }

  pb_finish(pbs, "%" AFFY_PRNd32 " probes", cdf->numprobes);
}

/* A quality control section */
static void process_qc_section(FILE *fp, 
                               AFFY_CDFFILE *cdf, 
                               LIBUTILS_PB_STATE *pbs, 
                               AFFY_ERROR *err)
{
  affy_int32 i, qc_counter, num_qc_cells = 0; 
  affy_int16 x, y;

  assert(fp  != NULL);
  assert(cdf != NULL);
        
  pb_begin(pbs, cdf->numqcunits, "Loading QC Units");

  /* For each qc unit, read in details */
  for (i = 0; i < cdf->numqcunits; i++) 
  {
    if (affy_readmulti(fp, "%x%dl", 2, (void *)&num_qc_cells) <= 0)
      AFFY_HANDLE_ERROR_VOID("couldn't read QC section", AFFY_ERROR_IO, err);

    for (qc_counter = 0; qc_counter < num_qc_cells; qc_counter++) 
    {
      if (affy_readmulti(fp, "%2hl%x", (void *)&x, (void *)&y, 3) <= 0)
	AFFY_HANDLE_ERROR_VOID("couldn't read QC section", AFFY_ERROR_IO, err);

      cdf->cell_type[x][y] = AFFY_QC_LOCATION;
#ifdef STORE_XY_REF
/*      cdf->xy_ref[x][y]    = NULL; */
#endif
    }
    
    pb_tick(pbs, 1,"Reading QC Unit %d",i+1);
  }

  pb_finish(pbs, "%" AFFY_PRNd32 " units", cdf->numqcunits);
}

/* A probe section will be indexed by probe, and by x,y coordinates */
/* 2017-07-19:
 *  EAW -- almost complete re-write to handle more modern binary files that
 *         have more complex unit/block structure
 */
static int process_probe_section(FILE *fp,
				 AFFY_CDFFILE *cdf,
				 int version,
				 int ps,
                                 LIBUTILS_PB_STATE *pbs,
				 AFFY_ERROR *err)
{
  affy_uint8    status = AFFY_NORMAL_LOCATION;
  affy_int32    atom, numprobes = 0;
  affy_uint16   x, y;
  int           j;
  affy_uint8    pbase, tbase;
  
  int rm_err;
  char          ps_name[65];
  affy_int32    block, cell;
  affy_uint16   unit_type;
  affy_int32    numblocks;
  affy_int32    numatoms;
  affy_int32    numcells;
  unsigned char cells_per_atom;
  
  int           all_mm_flag = 1;

  assert(fp  != NULL);
  assert(cdf != NULL);

  /* read unit type, numblocks */
  rm_err = affy_readmulti(fp, "%hl%x%dl%x",
                          (void *) &unit_type,
                          5,
                          (void *) &numblocks,
                          9);
  if (rm_err <= 0)
    AFFY_HANDLE_ERROR("probeset unit header read error", 
                      AFFY_ERROR_IO, err, -1);

#ifdef DEBUG
  fprintf(stderr, "UnitType:  %d\n", unit_type);
  fprintf(stderr, "NumBlocks: %d\n", numblocks);
#endif

  for (block = 0; block < numblocks; block++)
  {
    /* block header, up until block (probeset) name */
    rm_err = affy_readmulti(fp, "%2dl%c%x",
                            (void *) &numprobes,
                            (void *) &numcells,
                            (void *) &cells_per_atom,
                            9);
    if (rm_err <= 0)
        AFFY_HANDLE_ERROR("probeset block header read error", 
                          AFFY_ERROR_IO, err, -1);

    if (cells_per_atom < 2)
      all_mm_flag = 0;

#ifdef DEBUG
    fprintf(stderr, "NumProbes:    %d\n", numprobes);
    fprintf(stderr, "NumCells:     %d\n", numcells);
    fprintf(stderr, "CellsPerAtom: %d\n", cells_per_atom);
#endif

    cdf->probeset[ps].probe = NULL;
    if (affy_readchars(fp, ps_name, 65) == -1)
      AFFY_HANDLE_ERROR("couldn't read probeset name within block",
                        AFFY_ERROR_IO, err, -1);
    
    cdf->probeset[ps].name = h_strdup(ps_name);
    if (cdf->probeset[ps].name == NULL)
      AFFY_HANDLE_ERROR("strdup failed", AFFY_ERROR_OUTOFMEM, err, -1);
    hattach(cdf->probeset[ps].name, cdf);


#ifdef DEBUG
    fprintf(stderr, "%s\n", ps_name);
#endif
    
    if (version >= 2 && version <= 5)
    {
      /* Wobble situation, Allele code */
      if (affy_readmulti(fp, "%x", 4) <= 0)
        AFFY_HANDLE_ERROR("error in unused block section",
                          AFFY_ERROR_IO, err, -1);

      if (version >= 3)
      {
        /* Channel, RepType */
        if (affy_readmulti(fp, "%x", 2) <= 0)
          AFFY_HANDLE_ERROR("error in unused block section",
                            AFFY_ERROR_IO, err, -1);
      }
    }

    /* Allocate enough storage for all probes in this probe set */
    cdf->probeset[ps].probe = h_subcalloc(cdf, numprobes, sizeof(AFFY_PROBE));
    if (cdf->probeset[ps].probe == NULL)
      AFFY_HANDLE_ERROR("calloc failed", AFFY_ERROR_OUTOFMEM, err, -1);

    cdf->probeset[ps].numprobes = numprobes;
    cdf->probeset[ps].index     = ps;
    
    for (cell = 0; cell < numcells; cell++)
    {
      rm_err = affy_readmulti(fp, "%dl%2hl%x%2c",
                              (void *)&atom,
                              (void *)&x,
                              (void *)&y,
                              4,
                              (void *)&pbase,
                              (void *)&tbase);

      if (rm_err <= 0)
          AFFY_HANDLE_ERROR("probeset probe read error", 
                            AFFY_ERROR_IO, err, -1);

      if (version >= 2 && version <= 5)
      {
        /* Length, Physical grouping */
        if (affy_readmulti(fp, "%x", 4) <= 0)
          AFFY_HANDLE_ERROR("error in unused cell section",
                            AFFY_ERROR_IO, err, -1);

        if (version == 5)
        {
          /* Probe sequence ID */
          if (affy_readmulti(fp, "%x", 4) <= 0)
            AFFY_HANDLE_ERROR("error in unused cell section",
                              AFFY_ERROR_IO, err, -1);
        }
      }

      if (cdf->seen_xy[x][y] == 0)
        cdf->seen_xy[x][y] = 1;   /* seen it only once */
      else
        cdf->seen_xy[x][y] = 2;   /* seen it two or more times */

      cdf->cell_type[x][y] = status;

      /* HuEx-1_0-st-v2 text CDF file has non-letter cbase/pbase/tbase,
       *  treat them as PM.
       */
      if (pbase != tbase || !isalpha(pbase) || !isalpha(tbase))
      {
        cdf->probeset[ps].probe[cell].pm.x = x;
        cdf->probeset[ps].probe[cell].pm.y = y;
        
        /* fake MM coordinates if there are no MM */
        if (cells_per_atom == 1)
        {
          cdf->probeset[ps].probe[cell].mm.x = x;
          cdf->probeset[ps].probe[cell].mm.y = y;
        }
      }
      else
      {
        cdf->probeset[ps].probe[cell].mm.x = x;
        cdf->probeset[ps].probe[cell].mm.y = y;

        /* CDF says it is a MM, but it must really be a PM ? */
        /* store is as PM, otherwise things will fail badly */
        if (cells_per_atom == 1)
        {
          cdf->probeset[ps].probe[cell].pm.x = x;
          cdf->probeset[ps].probe[cell].pm.y = y;
        }
      }

#ifdef STORE_XY_REF
      cdf->xy_ref[x][y] = &(cdf->probeset[ps].probe[cell]);
#endif

      cdf->probeset[ps].probe[cell].ps = &(cdf->probeset[ps]);

      /* We've gone over the limit, due to shared probes between probeset.
       * Realloc'ing one-at-a-time is slow, but we're almost done by now.
       */
      if (cdf->numprobes >= cdf->numrows * cdf->numcols)
      {
        cdf->probe = h_realloc(cdf->probe,
                               (cdf->numprobes + 1) * sizeof(AFFY_PROBE *));
        if (cdf->probe == NULL)
          AFFY_HANDLE_ERROR_VOID_ZERO("realloc failed",
                                      AFFY_ERROR_OUTOFMEM, err);
      }

      /* Add probe to list of all probes */
      cdf->probe[cdf->numprobes]        = &(cdf->probeset[ps].probe[cell]);
      cdf->probe[cdf->numprobes]->index = cdf->numprobes;
      cdf->numprobes++;

      pb_tick(pbs, 1, "Reading probe %d", cell+1);
    }
  }
  
  if (numprobes)
  {
    if (all_mm_flag) return 1;
    return 0;
  }
  
  return -1;
}
