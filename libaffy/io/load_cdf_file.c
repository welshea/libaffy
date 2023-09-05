
/**************************************************************************
 *
 * Filename:  load_cdf_file.c
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
 * 10/05/07: Beautification/cleanup, new binary I/O functions (AMH)
 * 01/25/08: Repair binary I/O functions again (AMH)
 * 03/14/08: New error handling scheme (AMH)
 * 04/23/09: Add a routine to read CDF files by name (AMH)
 * 09/20/10: Pooled memory allocator (AMH)
 * 03/10/14: #ifdef out cdf->xy_ref code to save memory (EAW)
 * 03/14/14: updated to deal with exon arrays (EAW)
 * 06/01/18: change cdf malloc to calloc, so it is initialized to zeroes (EAW)
 * 08/12/20: store full path to CDF file in flags, so we can print later (EAW)
 * 09/05/23: fopen() everything as "rb" (EAW)
 *
 **************************************************************************/

#include <affy.h>

static bool file_readable(const char *filename)
{
  FILE *fp;

  if ((fp = fopen(filename, "rb")) == NULL)
  {
    return false;
  }
  else
  {
    fclose(fp);
    return true;
  }
}

/*
 * Load a CDF file given only a filename.  chip_type, if not NULL, is
 * used to initialize the array type field of the CDF structure.
 */
AFFY_CDFFILE *affy_load_cdf_file_byname(char *cdf_filename, 
                                        char *chip_type,
                                        AFFY_ERROR *err)
{
  AFFY_CDFFILE     *cdf = NULL;
  FILE             *fp;
  affy_int32        magic;
  LIBUTILS_PB_STATE pbs;

  assert(cdf_filename != NULL);

  pb_init(&pbs);

  /* Open file */
  fp = fopen(cdf_filename, "rb");
  if (fp == NULL)
    AFFY_HANDLE_ERROR("error opening CDF file", AFFY_ERROR_IO, err, NULL);

  info("Loading %s CDF file...", cdf_filename);

  /* Do a detection of the type, then call the appropriate code. */
  if (affy_read32_le(fp, (void *)&magic) != 0)
    AFFY_HANDLE_ERROR_GOTO("error reading CDF magic", 
                           AFFY_ERROR_IO, 
                           err, 
                           cleanup);

  cdf = h_calloc(1, sizeof(AFFY_CDFFILE));
  if (cdf == NULL)
    AFFY_HANDLE_ERROR_GOTO("malloc failed", AFFY_ERROR_OUTOFMEM, err, cleanup);

  if (chip_type != NULL)
  {
    cdf->array_type = h_strdup(chip_type);
    if (cdf->array_type == NULL)
      AFFY_HANDLE_ERROR_GOTO("strdup failed", 
                             AFFY_ERROR_OUTOFMEM,
                             err,
                             cleanup);
    hattach(cdf->array_type, cdf);
  }
  else
  {
    cdf->array_type = NULL;
  }

  cdf->cell_type    = NULL;
#ifdef STORE_XY_REF
  cdf->xy_ref       = NULL;
#endif
  cdf->probeset     = NULL;
  cdf->probe        = NULL;
  cdf->numprobesets = 0;
  cdf->numprobes    = 0;
  cdf->no_mm_flag   = 0;

  if (magic == AFFY_CDF_BINARYFILE_MAGIC) 
  {
    /* Reopen cdf file in binary mode - for PC's */
    fclose(fp);

    if ((fp = fopen(cdf_filename, "rb")) == NULL)
      AFFY_HANDLE_ERROR_GOTO("couldn't reopen CDF in binary mode", 
                             AFFY_ERROR_IO,
                             err, 
                             cleanup);

    /* Load the file */
    affy_load_binary_cdf_file(fp, cdf, &pbs, err);
  } 
  else 
  {
    rewind(fp);
    affy_load_text_cdf_file(fp, cdf, &pbs, err);
  }

cleanup:
  fclose(fp);
  pb_cleanup(&pbs);

  if (err->type != AFFY_ERROR_NONE)
  {
    affy_free_cdf_file(cdf);

    return (NULL);
  }
  
  return (cdf);
}

/*
 * The CDF file contains the description of the microarray chip in terms
 * of probes and probe sets. This information needs to be loaded into
 * memory so that calculations at the probe level can be matched
 * to (X,Y) coordinates in the CEL file.
 * 
 * The parameter is a chip_type as defined by the CEL file. This should
 * correspond to a cdf file somewhere.
 */
AFFY_CDFFILE *affy_load_cdf_file(char *chip_type, char *dir,
                                 AFFY_COMBINED_FLAGS *f,
                                 AFFY_ERROR *err)
{
  char cdf_filename[MAXBUF];

  assert(chip_type != NULL);

  /* Try to find cdf file. */
  if (dir != NULL)
  {
    /* Check if dir is actually a cdf file */
    if (endsWith(dir, ".CDF") && file_readable(dir))
    {
      sprintf(cdf_filename, "%s", dir);
      goto FOUND;
    }
    /* Check if dir is actually a cdf file */
    if (endsWith(dir, ".cdf") && file_readable(dir))
    {
      sprintf(cdf_filename, "%s", dir);
      goto FOUND;
    }

    /* Otherwise, treat it as a directory */
    sprintf(cdf_filename, "%s/%s.CDF", dir, chip_type);
    if (file_readable(cdf_filename))
      goto FOUND;

    sprintf(cdf_filename, "%s/%s.cdf", dir, chip_type);
    if (file_readable(cdf_filename))
      goto FOUND;
  }

  /* Even if not asked, check the current directory */
  sprintf(cdf_filename, "%s.CDF", chip_type);
  if (file_readable(cdf_filename))
    goto FOUND;

  sprintf(cdf_filename, "%s.cdf", chip_type);
  if (file_readable(cdf_filename))
    goto FOUND;

  /* If we got here, can't find it and die */
  AFFY_HANDLE_ERROR("can't locate CDF file", AFFY_ERROR_NOTFOUND, err, NULL);


FOUND:

  f->cdf_filename = strdup(cdf_filename);
  
  /* continuation of print_flags(), since we cannot know it ahead of time */
  printf("Path to CDF file:                    %s\n\n", f->cdf_filename);

  return (affy_load_cdf_file_byname(cdf_filename, chip_type, err));
}
