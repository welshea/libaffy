
/**************************************************************************
 *
 * Filename:  load_cel_file.c
 *
 * Purpose:   Parse a CEL file and initialize an accompanying structure.
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
 * 04/19/05: Now we use AFFY_CELL and collect all cell data (AMH)
 * 10/01/07: Add support for a new format: Calvin (AMH)
 * 03/14/08: New error handling scheme (AMH)
 * 09/20/10: Pooled memory allocator (AMH)
 * 09/19/12: Added sanity checker for NaN and Inf (EAW)
 *
 **************************************************************************/

#include <affy.h>

#if PARANOID_CEL_LOADER
int affy_cel_sanity_fix(AFFY_CELFILE *cf)
{
    AFFY_CELL **data = cf->data;
    int numrows = cf->numrows;
    int numcols = cf->numcols;
    int row, col;
    int num_bogus = 0;
    double value;
    
    for (row = 0; row < numrows; row++)
    {
        for (col = 0; col < numcols; col++)
        {
            value = data[row][col].value;
            
            if (isinf(value) || isnan(value))
            {
                data[row][col].value = 0.0;
                num_bogus++;
            }
        }
    }

    return num_bogus;
}
#endif

AFFY_CELFILE *affy_load_cel_file(char *filename, AFFY_ERROR *err)
{
  FILE             *fp;
  AFFY_CELFILE     *cf = NULL;
  affy_int32        int_magic;
  affy_uint8        byte_magic;
  LIBUTILS_PB_STATE pbs;
#if PARANOID_CEL_LOADER
  int num_bogus = 0;
#endif
        
  assert(filename != NULL);

  pb_init(&pbs);

  /* Open file. */
  fp = fopen(filename, "rb");
  if (fp == NULL)
    AFFY_HANDLE_ERROR("couldn't open CEL file", AFFY_ERROR_NOTFOUND, err, NULL);

  info("Loading CEL file %s", filename);

  cf = h_malloc(sizeof(AFFY_CELFILE));
  if (cf == NULL)
    AFFY_HANDLE_ERROR_GOTO("malloc failed", AFFY_ERROR_OUTOFMEM, err, done);

  cf->filename = h_strdup(filename);
  if (cf->filename == NULL)
    AFFY_HANDLE_ERROR_GOTO("strdup failed", AFFY_ERROR_OUTOFMEM, err, done);
  hattach(cf->filename, cf);
  
  cf->corrupt_flag = 0;

  /* 
   * Check the file magic and call appropriate loading function.  The
   * reading of two magic values is needed since the old binary format
   * magic is stored in a full int, which may be backwards due to
   * endianness, bla bla...
   */
  if (affy_read32_le(fp, (void *)(&int_magic)) < 0)
    AFFY_HANDLE_ERROR_GOTO("I/O error reading CEL magic", 
                           AFFY_ERROR_IO, 
                           err,
                           done);

  rewind(fp);

  if (affy_read8(fp, (void *)(&byte_magic)) < 0)
    AFFY_HANDLE_ERROR_GOTO("I/O error reading CEL magic", 
                           AFFY_ERROR_IO, 
                           err, 
                           done);

  rewind(fp);

  /* Check for Calvin magic first. */
  if (byte_magic == AFFY_CALVIN_FILEMAGIC)
  {
    /* Calvin (Command Console) "generic" format. */
    affy_load_calvin_cel_file(fp, cf, &pbs, err);
  }
  else if (int_magic == AFFY_CEL_BINARYFILE_MAGIC)
  {
    /* "Old" binary (GC) format. */
    affy_load_binary_cel_file(fp, cf, &pbs, err);
  } 
  else 
  {
    /* Assumed to be text; reopen in text format (probably not necessary) */
    fclose(fp);
    fp = fopen(filename, "r");
    if (fp == NULL)
      AFFY_HANDLE_ERROR_GOTO("error reopening CEL file", 
                             AFFY_ERROR_IO, 
                             err, 
                             done);

    affy_load_text_cel_file(fp, cf, &pbs, err);
  }
  
#if PARANOID_CEL_LOADER
  num_bogus = affy_cel_sanity_fix(cf);
  if (num_bogus)
  {
    info("Zeroed %d nan/inf values in CEL file %s", num_bogus, filename);
  }
#endif

done:        
  fclose(fp);
  pb_cleanup(&pbs);

  if (err->type != AFFY_ERROR_NONE)
    affy_free_cel_file(cf);

  return (cf);
}
