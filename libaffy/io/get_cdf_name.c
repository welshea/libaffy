
/**************************************************************************
 *
 * Filename:  get_cdf_name.c
 *
 * Purpose:   Extract CDF filename from DatHeader in CEL file.
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
 * 12/11/06: Add assertions (AMH)
 * 04/08/05: Imported/repaired from old libaffy (AMH)
 * 10/05/07: Update to use new binary I/O functions (AMH)
 * 03/11/08: New error handling scheme (AMH)
 * 09/20/10: Pooled memory allocator (AMH)
 * 09/05/23: change fgets() calls to EOL-safe functions (EAW)
 *
 **************************************************************************/

#include <affy.h>

#define HDR_OFFSET 20

/*
 * Given a string which is the DatHeader line of a CEL file, figure
 * out the corresponding CDF file.
 */

char *affy_get_cdf_name(const char *buf, AFFY_ERROR *err)
{
  char *p;
  char *result;

  assert(buf != NULL);

  /* Get the cdf name */
  p = strstr(buf, ".1sq");

  if (p == NULL)
    AFFY_HANDLE_ERROR("bad DatHeader format", AFFY_ERROR_BADFORMAT, err, NULL);

  *p = 0;

  while (p > buf && *(p - 1) != ' ' && *(p - 1) != 024)
    p--;

  result = h_strdup(p);
  if (result == NULL)
    AFFY_HANDLE_ERROR("out of memory for strdup", 
		      AFFY_ERROR_OUTOFMEM, 
		      err, 
		      NULL);

  return (result);
}

/*
 * Given a cel file name, open it and find the corresponding
 * cdf name.
 */
char *affy_get_cdf_name_from_cel(const char *filename, AFFY_ERROR *err)
{
  FILE       *fp;
  affy_int32  int_magic;
  affy_uint8  byte_magic;
  int         max_string_len = 0;
  char        buf[MAXBUF];
  char       *fgets_buffer = NULL;
  char       *result       = NULL;

  assert(filename != NULL);

  /* Figure out if file is txt (V3) or binary (V4) */
  fp = fopen(filename, "rb");
  if (fp == NULL)
  {
    fprintf(stderr, "CEL file: %s\n", filename);
    AFFY_HANDLE_ERROR("CEL file fopen failed", AFFY_ERROR_NOTFOUND, err, NULL);
  }

  if (affy_read32_le(fp, (void *)(&int_magic)) < 0)
    AFFY_HANDLE_ERROR("I/O error reading CEL magic", 
                      AFFY_ERROR_IO, 
                      err,
                      NULL);
  
  rewind(fp);

  if (affy_read8(fp, (void *)(&byte_magic)) < 0)
    AFFY_HANDLE_ERROR("I/O error reading CEL magic", 
                      AFFY_ERROR_IO, 
                      err, 
                      NULL);

  rewind(fp);

  if (int_magic == AFFY_CEL_BINARYFILE_MAGIC)
  {
    affy_int32 hdrlen;
    int        readin = 0;

    if (fseek(fp, HDR_OFFSET, SEEK_SET) == -1)
    {
      fclose(fp);
      AFFY_HANDLE_ERROR("seek failed", AFFY_ERROR_IO, err, NULL);
    }
   
    if (affy_read32_le(fp, (void *)&hdrlen) != 0)
    {
      fclose(fp);
      AFFY_HANDLE_ERROR("couldn't read CEL header", AFFY_ERROR_IO, err, NULL);
    }

    while (affy_readchars(fp, buf, MAXBUF) && readin < hdrlen)
    {
      if (!strncmp(buf, "DatHeader=", 10))
        break;
      readin = ftell(fp) - HDR_OFFSET + 4;
    }

    result = affy_get_cdf_name(buf, err);

    fclose(fp);
  }
  else if (byte_magic == AFFY_CALVIN_FILEMAGIC)
  {
    AFFY_CALVINIO          *cio;
    AFFY_CALVIN_DATAHEADER *dh;
    AFFY_CALVIN_PARAM      *cp;
    
    cio = affy_calvinio_init(fp, err);
    AFFY_CHECK_ERROR(err, NULL);

    dh = affy_calvin_get_dataheader(cio, err);
    AFFY_CHECK_ERROR(err, NULL);

    cp = affy_calvin_find_param(dh->params, 
                                dh->num_params, 
                                "affymetrix-array-type");

    if (cp != NULL)
    {
      result = h_strdup(cp->value.string_val);

      affy_free_calvin_dataheader(dh);
      affy_calvinio_free(cio);
      fclose(fp);
      
      if (result == NULL)
        AFFY_HANDLE_ERROR("strdup failed", AFFY_ERROR_OUTOFMEM, err, NULL);
    }
    else
    {
      affy_free_calvin_dataheader(dh);
      affy_calvinio_free(cio);
      fclose(fp);
      
      AFFY_HANDLE_ERROR("couldn't determine Calvin array type",
                        AFFY_ERROR_BADFORMAT,
                        err,
                        NULL);
    }
  }
  else /* text CEL file */
  {
    /* Find the dat header line */
    while (fgets_strip_realloc(&fgets_buffer, &max_string_len, fp) != NULL)
      if (!strncmp(fgets_buffer, "DatHeader=", 10))
        break;
    
    result = affy_get_cdf_name(fgets_buffer, err);

    fclose(fp);
    
    if (fgets_buffer)
        free(fgets_buffer);

    AFFY_CHECK_ERROR(err, NULL);
  }

  return (result);
}
