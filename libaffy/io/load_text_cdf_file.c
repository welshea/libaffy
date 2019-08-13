
/**************************************************************************
 *
 * Filename:  load_text_cdf_file.c
 *
 * Purpose:   Parse a CDF file and initialize an accompanying structure.
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
 * 04/08/05: Imported/repaired from old libaffy (AMH)
 * 03/14/08: New error handling scheme (AMH)
 * 04/18/08: Text I/O changes (AMH)
 * 09/20/10: Pooled memory allocator (AMH)
 * 06/20/12: added support for Human Exon arrays (EAW)
 * 09/19/12: added proper support for multiple banks per probeset (EAW)
 * 09/20/12: Fixed: pbase, tbase, atom fields were all off by one field (EAW)
 *            behavior is probably unchanged due to luck regarding code and
 *            the surrounding fields all working out to the same end result
 * 09/20/12: handled more problems in official unofficial HuEx CDF file
 * 03/06/14: fixed row/col memory allocation, dimensions were swapped (EAW)
 * 03/10/14: cdf->seen_xy, cdf->no_mm_flag, cdf->dupe_probes_flag (EAW)
 * 03/10/14: #ifdef out cdf->xy_ref code to save memory (EAW)
 * 03/12/14: print Number of Probests after probes are loaded (EAW)
 * 10/22/18: support broken BrainArray CDFs with all MM probes missing (EAW)
 * 05/23/19: workaround more errors in HuEx-1_0-st-v2.text.cdf controls (EAW)
 *
 **************************************************************************/

#include <affy.h>
#include <utils.h>

static void process_chip_section(AFFY_TEXTIO *tf, 
                                 AFFY_CDFFILE *cdf,
                                 AFFY_ERROR *err);
static void process_qc_section(AFFY_TEXTIO *tf, 
                               AFFY_CDFFILE *cdf, 
                               int qcnum,
                               LIBUTILS_PB_STATE *pbs,
                               AFFY_ERROR *err);
static int  process_probe_section(AFFY_TEXTIO *tf, 
                                  AFFY_CDFFILE *cdf, 
                                  int *probe_set_num_ptr,
                                  LIBUTILS_PB_STATE *pbs,
                                  char *old_probeset_name,
                                  AFFY_ERROR *err);

void affy_load_text_cdf_file(FILE *fp, 
                             AFFY_CDFFILE *cdf,
                             LIBUTILS_PB_STATE *pbs,
                             AFFY_ERROR *err)
{
  char              *kv[2], *str;
  char              *old_probeset_name = NULL;
  int                unit_no, block_no, probe_set_num = 0;
  AFFY_TEXTIO       *tf;
  bool               pb_started = false;
  int                all_mm_flag = 1, temp_int;
  affy_uint32        p;

  assert(fp  != NULL);
  assert(cdf != NULL);
  
  /* Initialize text parser. */
  tf = affy_textio_init(fp, err);
  AFFY_CHECK_ERROR_VOID(err);
  
  /* Parse file, looking at [] section headings */
  affy_textio_reset_next_line(tf);

  while ((str = affy_textio_get_next_line(tf)) != NULL)
  {
    /* Process by section */
    if (STREQ(str, "[CDF]"))
    {
      str = affy_textio_get_next_line(tf);
      if (str == NULL)
        AFFY_HANDLE_ERROR_GOTO("error parsing CDF", 
                               AFFY_ERROR_BADFORMAT, 
                               err, 
                               out);
        
      split(str, kv, '=', 2);
      if (STREQ(kv[0], "Version"))
        info("Found ASCII CDF version %s", kv[1]);
    }
    else if (STREQ(str, "[Chip]"))
    {
      process_chip_section(tf, cdf, err);
      AFFY_CHECK_ERROR_GOTO(err, out);

    }
    else if (!strncmp(str, "[QC", 3))
    {
      char *err_str;
      int   qcnum;

      qcnum = strtol(str + 3, &err_str, 10);
      if (err_str == (str + 3))
        AFFY_HANDLE_ERROR_GOTO("couldn't parse QC unit number in CDF",
                               AFFY_ERROR_BADFORMAT,
                               err,
                               out);
      
      process_qc_section(tf, cdf, qcnum, pbs, err);
      AFFY_CHECK_ERROR_GOTO(err, out);
    }
    else if (sscanf(str, "[Unit%d_Block%d]", &unit_no, &block_no) == 2)
    {
      if ( !pb_started ) {
         pb_begin(pbs, cdf->numprobesets, "Loading CDF File");
         pb_msg(pbs, "Loading Probesets...");
         pb_started = true;
      }
      temp_int = process_probe_section(tf, cdf, &probe_set_num, pbs,
                                       old_probeset_name, err);
      AFFY_CHECK_ERROR_GOTO(err, out);
      pb_tick(pbs, 1,"Loaded probeset %d",probe_set_num);
      
      /* we didn't read any probes, so skip this probeset */
      if (temp_int < 0)
        continue;

      /* remember, probe_set_num was incremented in process_probe_section() */
      old_probeset_name = cdf->probeset[probe_set_num-1].name;

      /* we read in some probes without a MM, so we can't realloc to save
       * memory
       */
      if (temp_int == 0)
        all_mm_flag = 0;
    }
    else
    {
      affy_textio_skip_to_next_header(tf);
    }
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
  for (p = 0; p < cdf->numrows * cdf->numcols; p++)
  {
    if (cdf->seen_xy[0][p] == 2)
    {
      cdf->dupe_probes_flag = 1;
      break;
    }
  }

out:
  affy_textio_free(tf);

  /* So that the free function works properly... */
  cdf->numprobesets = probe_set_num;

  if (pb_started)
    pb_finish(pbs, "%" AFFY_PRNd32 " probes", cdf->numprobes);
  
  info("Number of Probesets: %d", cdf->numprobesets);
}

static void process_chip_section(AFFY_TEXTIO *tf, 
                                 AFFY_CDFFILE *cdf, 
                                 AFFY_ERROR *err)
{
  char *s, *err_str = NULL, *kv[2];
  int   i;

  assert(tf  != NULL);
  assert(cdf != NULL);

  /* Read lines until eof or a new section */
  while ((s = affy_textio_get_next_line(tf)) != NULL)
  {
    int numsplit;

    if (*s == '[')
    {
      affy_textio_unget_next_line(tf);
      break;
    }

    numsplit = split(s, kv, '=', 2);
    if (numsplit == 2)
    {
      if (STREQ(kv[0], "Rows"))
        cdf->numrows = strtol(kv[1], &err_str, 10);
      else if (STREQ(kv[0], "Cols"))
        cdf->numcols = strtol(kv[1], &err_str, 10);
      else if (STREQ(kv[0], "NumQCUnits"))
        cdf->numqcunits = strtol(kv[1], &err_str, 10);
      else if (STREQ(kv[0], "NumberOfUnits"))
        cdf->numprobesets = strtol(kv[1], &err_str, 10);
      
      if (kv[1] == err_str)
        AFFY_HANDLE_ERROR_VOID("bad chip section in CDF file", 
                               AFFY_ERROR_BADFORMAT, 
                               err);
    }
  }

  /* 
   * At this point we know the size of the chip, and the number of
   * probes.
   */
  cdf->cell_type = h_subcalloc(cdf, cdf->numcols, sizeof(affy_uint8 *));
  if (cdf->cell_type == NULL)
    AFFY_HANDLE_ERROR_VOID("calloc failed", AFFY_ERROR_OUTOFMEM, err);

  cdf->cell_type[0] = h_subcalloc(cdf->cell_type, 
                                  cdf->numrows * cdf->numcols,
                                  sizeof(affy_uint8));
  if (cdf->cell_type[0] == NULL)
    AFFY_HANDLE_ERROR_VOID("calloc failed", AFFY_ERROR_OUTOFMEM, err);

  for (i = 1; i < cdf->numcols; i++)
    cdf->cell_type[i] = cdf->cell_type[i-1] + cdf->numrows;

#ifdef STORE_XY_REF
  cdf->xy_ref = h_subcalloc(cdf, cdf->numcols, sizeof(AFFY_PROBE **));
  if (cdf->xy_ref == NULL)
    AFFY_HANDLE_ERROR_VOID("calloc failed", AFFY_ERROR_OUTOFMEM, err);
  
  cdf->xy_ref[0] = h_subcalloc(cdf->xy_ref, cdf->numrows * cdf->numcols,
                               sizeof(AFFY_PROBE *));
  if (cdf->xy_ref[0] == NULL)
    AFFY_HANDLE_ERROR_VOID("calloc failed", AFFY_ERROR_OUTOFMEM, err);

  for (i = 1; i < cdf->numcols; i++)
    cdf->xy_ref[i] = cdf->xy_ref[i-1] + cdf->numrows;
#endif

  cdf->probeset = h_subcalloc(cdf, cdf->numprobesets, sizeof(AFFY_PROBESET));
  if (cdf->probeset == NULL)
    AFFY_HANDLE_ERROR_VOID("calloc failed", AFFY_ERROR_OUTOFMEM, err);

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
    AFFY_HANDLE_ERROR_VOID("calloc failed", AFFY_ERROR_OUTOFMEM, err);
}

/* A quality control section */
static void process_qc_section(AFFY_TEXTIO *tf, 
                               AFFY_CDFFILE *cdf, 
                               int qcnum,
                               LIBUTILS_PB_STATE *pbs,
                               AFFY_ERROR *err)
{
  char *s, *err_str, *kv[2];
  int   num_qc_cells = 0, qc_counter;
  int   x, y;

  /* Find number of cells, then cell header */
  while ((s = affy_textio_get_next_line(tf)) != NULL)
  {
    split(s, kv, '=', 2);
 
    if (STREQ(kv[0], "CellHeader"))
      break;
    if (STREQ(kv[0], "NumberCells"))
    {
      num_qc_cells = strtol(kv[1], &err_str, 10);
      if (err_str == kv[1])
        AFFY_HANDLE_ERROR_VOID("bad QC section in CDF file", 
                               AFFY_ERROR_BADFORMAT, 
                               err);
    }
  }

  for (qc_counter = 0; qc_counter < num_qc_cells; qc_counter++)
  {
    int ret;

    s = affy_textio_get_next_line(tf);
    if (s == NULL)
      AFFY_HANDLE_ERROR_VOID("bad QC section in CDF file",
                             AFFY_ERROR_BADFORMAT,
                             err);
    split(s, kv, '=', 2);

    ret = sscanf(kv[1], "%d%d", &x, &y);
    
    if ((ret != 2) 
        || (x < 0) 
        || (y < 0) 
        || (x >= cdf->numcols)
        || (y >= cdf->numrows))
      AFFY_HANDLE_ERROR_VOID("bad QC section in CDF file", 
                             AFFY_ERROR_BADFORMAT, 
                             err);

    cdf->cell_type[x][y] = AFFY_QC_LOCATION;
#ifdef STORE_XY_REF
/*    cdf->xy_ref[x][y]    = NULL; */
#endif
  }
}

/* A probe section will be indexed by probe, and by x,y coordinates
 * Return values: 0,1 -- all_mm_flag status
                 -1   -- did not read in any probesets
 */
static int process_probe_section(AFFY_TEXTIO *tf, 
                                 AFFY_CDFFILE *cdf, 
                                 int *probe_set_num_ptr,
                                 LIBUTILS_PB_STATE *pbs,
                                 char *old_probeset_name,
                                 AFFY_ERROR *err)
{
  int    ps = *probe_set_num_ptr;
  int    status = AFFY_NORMAL_LOCATION, numprobes = 0;
  int    numatoms=0, numcells = 0, cells_per_atom = 2;
  int    read_in_a_probe_flag = 0;
  int    read_in_a_probeset_flag = 0;
  int    all_mm_flag = 1;
  char  *kv[2], *s;
  int    x, y, i, j;
  int    old_numprobes, new_numprobes;
  int    mm_count, pm_count;
  char   pbase, tbase;
  int    atom, old_atom;

  assert(tf  != NULL);
  assert(cdf != NULL);

  /* Either the CDF file header lied and there are extra probesets,
   * or we have a multiblock probeset at the very end of the chip
   * (which would cause an off-by-one error due to how things are coded).
   * 
   * If multiblock, then things are fine and will be dealt with appropriately,
   * but if not, then the CDF file header lied and we'll just have to hope
   * that things are otherwise kosher.
   */
  if (ps >= cdf->numprobesets)
  {
    cdf->probeset = h_realloc(cdf->probeset, (cdf->numprobesets + 1) *
                              sizeof(AFFY_PROBESET));
    memset(&cdf->probeset[cdf->numprobesets], 0, sizeof(AFFY_PROBESET));
    cdf->numprobesets++;
  }

  /* 
   * Ensure that dynamically allocated stuff has a default NULL
   * value. 
   */
  cdf->probeset[ps].name  = NULL;
  cdf->probeset[ps].probe = NULL;

  /* Find number of cells, then cell header */
  while ((s = affy_textio_get_next_line(tf)) != NULL)
  {
    split(s, kv, '=', 2);

    /* Name is important, also sets the location type */
    if (STREQ(kv[0], "Name"))
    {
      cdf->probeset[ps].name = h_strdup(kv[1]);
      
      if (cdf->probeset[ps].name == NULL)
        AFFY_HANDLE_ERROR("strdup failed", AFFY_ERROR_OUTOFMEM, err, -1);
      hattach(cdf->probeset[ps].name, cdf);
    }
    else if (STREQ(kv[0], "NumAtoms"))
    {
      char *err_str;

      numatoms = strtol(kv[1], &err_str, 10);
      if (err_str == kv[1])
        AFFY_HANDLE_ERROR("couldn't parse probeset probe count",
                          AFFY_ERROR_BADFORMAT, err, -1);
      numprobes = numatoms;
    }
    else if (STREQ(kv[0], "NumCells"))
    {
      char *err_str;

      numcells = strtol(kv[1], &err_str, 10);
      if (err_str == kv[1])
        AFFY_HANDLE_ERROR("couldn't parse probeset probe count",
                           AFFY_ERROR_BADFORMAT, err, -1);
    }
    else if (STREQ(kv[0], "CellHeader"))
      break;
  }

  /* By this point, numprobes should have a valid value. */
  if (numprobes <= 0 && numcells <= 0)
    AFFY_HANDLE_ERROR("bad number of probes in probeset section",
                      AFFY_ERROR_BADFORMAT, err, -1);

  /* Handle NumAtoms=0 control probes, which are stored oddly */
  if (numprobes <= 0 && numcells > 0)
    numprobes = numcells;

  /* Handle more control probe errors in officially unsupported
   *  HuEx-1_0-st-v2.text.cdf file
   *
   * Example: Unit4057134_Block1; NumAtoms=1 NumCells=162
   */
  if (numprobes == 1 && numcells / cells_per_atom > numprobes)
    numprobes = numcells;

  /* Usually == 2, since each PM has a MM.  Can == 1 for exon arrays... */
  cells_per_atom = numcells / numprobes;
  
  
  if (cells_per_atom < 2)
    all_mm_flag = 0;

  /* deal with multiple sequential blocks per probeset */
  old_numprobes = 0;
  new_numprobes = numprobes;
  if (old_probeset_name &&
      strcmp(cdf->probeset[ps].name, old_probeset_name) == 0)
  {
    /* de-increment the ps, since it was erroneously pre-incremented */
    ps--;
  
    old_numprobes = cdf->probeset[ps].numprobes;
    new_numprobes = old_numprobes + numprobes;
    
#if 0
    printf("MultiBlock\t%s\t%d\t%d\t%d\n",
           cdf->probeset[ps].name, old_numprobes, numprobes, new_numprobes);
#endif
  }

  /* Allocate enough storage for all probes in this probe set */
  cdf->probeset[ps].probe = h_subcalloc(cdf, new_numprobes,
                                        sizeof(AFFY_PROBE));
  if (cdf->probeset[ps].probe == NULL)
    AFFY_HANDLE_ERROR("calloc failed", AFFY_ERROR_OUTOFMEM, err, -1);

  cdf->probeset[ps].numprobes = new_numprobes;
  cdf->probeset[ps].index     = ps;

  /* Get each probe */
  read_in_a_probeset_flag = 0;
  mm_count = 0;
  pm_count = 0;
  for (i = 0; i < numprobes; i++)
  {
    read_in_a_probe_flag = 0;

    atom = old_atom = -1;
    for (j = 0; j < cells_per_atom; j++)
    {
      char *err_str;
      char *f[16];

      s = affy_textio_get_next_line(tf);

      /* check for corrupt BrainArray CDF file where they lied about
       * having MM probes but actually didn't provide any
       */
      if (s == NULL || *s == '[')
      {
          /* we've read everything in that matters, back up and move on */
          if (pm_count == numprobes && mm_count == 0)
          {
              all_mm_flag = 0;
              affy_textio_unget_next_line(tf);
              break;
          }
      }

      if (s == NULL)
        AFFY_HANDLE_ERROR("bad probeset section in CDF",
                          AFFY_ERROR_BADFORMAT, err, -1);

      split(s, kv, '=', 2);
      split(kv[1], f, '\t', 16);

      old_atom = atom;
      atom = strtol(f[10], &err_str, 10);
      if (err_str == f[10])
        AFFY_HANDLE_ERROR("bad probeset section in CDF",
                          AFFY_ERROR_BADFORMAT, err, -1);

      /* check for corrupt BrainArray CDF file where they lied about
       * having MM probes but actually didn't provide any
       */
      if (j && atom != old_atom)
      {
          /* all MM are missing so far, back up and read in next probe */
          if (pm_count && mm_count == 0)
          {
              all_mm_flag = 0;
              affy_textio_unget_next_line(tf);
              continue;
          }

          AFFY_HANDLE_ERROR("bad probeset section in CDF",
                            AFFY_ERROR_BADFORMAT, err, -1);
      }
      
      x = strtol(f[0], &err_str, 10);
      if (err_str == f[0])
        AFFY_HANDLE_ERROR("bad probeset section in CDF",
                          AFFY_ERROR_BADFORMAT, err, -1);

      y = strtol(f[1], &err_str, 10);
      if (err_str == f[1])
        AFFY_HANDLE_ERROR("bad probeset section in CDF",
                          AFFY_ERROR_BADFORMAT, err, -1);
      
      if (cdf->seen_xy[x][y] == 0)
        cdf->seen_xy[x][y] = 1;   /* seen it only once */
      else
        cdf->seen_xy[x][y] = 2;   /* seen it two or more times */

      pbase = *f[8];
      tbase = *f[9];

      cdf->cell_type[x][y] = status;

      /* HuEx-1_0-st-v2 CDF file has non-letter cbase/pbase/tbase,
       *  treat them as PM.
       */
      if (pbase != tbase || !isalpha(pbase) || !isalpha(tbase))
      {
        cdf->probeset[ps].probe[old_numprobes + i].pm.x = x;
        cdf->probeset[ps].probe[old_numprobes + i].pm.y = y;
        
        /* fake MM coordinates if there are no MM */
        if (cells_per_atom == 1)
        {
          cdf->probeset[ps].probe[old_numprobes + i].mm.x = x;
          cdf->probeset[ps].probe[old_numprobes + i].mm.y = y;
        }
        
        pm_count++;
      }
      else
      {
        cdf->probeset[ps].probe[old_numprobes + i].mm.x = x;
        cdf->probeset[ps].probe[old_numprobes + i].mm.y = y;

        /* CDF says it is a MM, but it must really be a PM ? */
        /* store is as PM, otherwise things will fail badly */
        if (cells_per_atom == 1)
        {
          cdf->probeset[ps].probe[old_numprobes + i].pm.x = x;
          cdf->probeset[ps].probe[old_numprobes + i].pm.y = y;
          
          pm_count++;
        }
        else
        {
          mm_count++;
        }
      }

#ifdef STORE_XY_REF
      cdf->xy_ref[x][y] = &(cdf->probeset[ps].probe[old_numprobes + i]);
#endif

      read_in_a_probe_flag = 1;
    }

    /* we didn't read in any probes */
    if (read_in_a_probe_flag == 0)
      continue;
    
    cdf->probeset[ps].probe[old_numprobes + i].ps = &(cdf->probeset[ps]);

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
    cdf->probe[cdf->numprobes] =
      &(cdf->probeset[ps].probe[old_numprobes + i]);
    cdf->probe[cdf->numprobes]->index = cdf->numprobes;
    cdf->numprobes++;

    read_in_a_probeset_flag = 1;
  }

  if (pm_count != numprobes)
  {
      fprintf(stderr, "Problematic probeset: %s %d %d\n", cdf->probeset[ps].name, numprobes, pm_count);
      AFFY_HANDLE_ERROR("bad probeset section in CDF, not enough probes to fill probeset",
                        AFFY_ERROR_BADFORMAT, err, -1);
  }
  
  if (read_in_a_probeset_flag)
  {
    /* increment probeset index, whether or not there are multiblocks next */
    *probe_set_num_ptr = ps + 1;

    return all_mm_flag;
  }

  return -1;
}
