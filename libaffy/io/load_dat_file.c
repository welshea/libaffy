
/**************************************************************************
 *
 * Filename:  load_dat_file.c
 * 
 * Purpose:   Binary DAT file reading routines.
 *
 * Creation:  04/01/05
 *
 * Author:    Andrew M. Hoerter
 *
 * Copyright: Copyright (C) 2007, Moffitt Cancer Center.
 *            All rights reserved.
 *
 * Update History
 * --------------
 * 04/01/05: Creation (AMH)
 * 04/19/05: Now we use an AFFY_PIXREGION, plus other loose ends (AMH)
 * 10/03/07: Minor beautification/cleanups (AMH)
 * 10/05/07: New binary I/O functions (AMH)
 * 01/25/08: Convert to  use readmulti() (AMH)
 * 03/14/08: New error handling scheme (AMH)
 * 09/20/10: Pooled memory allocator (AMH)
 *
 **************************************************************************/

#include <affy.h>

#define TMPBUF_SIZE 250
#define FIELD_START "\x14 "
#define FIELD_END   " \x14"

/* 
 * affy_load_dat_file(): read in a DAT file and fill out a new AFFY_DATFILE
 *                       structure.
 *
 * Inputs: *filename contains the desired filename.
 * Outputs: pointer to AFFY_DATFILE structure on success, NULL on failure
 */
AFFY_DATFILE *affy_load_dat_file(char *filename, AFFY_ERROR *err)
{
  AFFY_DATFILE     *newdat;
  FILE             *input;
  affy_uint8        magic;
  affy_uint32       i;
  char              tmpbuf[TMPBUF_SIZE];
  char             *walk, *save;
  LIBUTILS_PB_STATE pbs;

  assert(filename != NULL);

  pb_init(&pbs);

  if ((input = fopen(filename, "rb")) == NULL)
    AFFY_HANDLE_ERROR("couldn't open DAT file", AFFY_ERROR_NOTFOUND, err, NULL);

  newdat = (AFFY_DATFILE *)h_malloc(sizeof(AFFY_DATFILE));
  if (newdat == NULL)
    AFFY_HANDLE_ERROR_GOTO("malloc failed", AFFY_ERROR_OUTOFMEM, err, cleanup);

  newdat->experiment_name  = 
  newdat->probe_array_type = NULL;

  newdat->pixels.data = NULL;

  /* Begin actual parsing of the file. */

  /* Check magic. */
  if (affy_read8(input, &magic) != 0)
    AFFY_HANDLE_ERROR_GOTO("error reading DAT file magic", 
                           AFFY_ERROR_IO, 
                           err, 
                           cleanup);

  if (magic != AFFY_DAT_FILEMAGIC)
    AFFY_HANDLE_ERROR_GOTO("bad DAT file magic", 
                           AFFY_ERROR_BADFORMAT, 
                           err, 
                           cleanup);

  if (affy_readmulti(input, "%2hl%3dl%2fl%x",
                     (void *)&newdat->pixels.numcols,
                     (void *)&newdat->pixels.numrows,
                     (void *)&newdat->numpixels,
                     (void *)&newdat->minpixel,
                     (void *)&newdat->maxpixel,
                     (void *)&newdat->meanpixel,
                     (void *)&newdat->std_dev_pixel,
                     18) <= 0)
    AFFY_HANDLE_ERROR_GOTO("I/O error reading DAT file", 
                           AFFY_ERROR_IO, 
                           err, 
                           cleanup);

  if (affy_readchars(input, tmpbuf, 8) != 0)
    AFFY_HANDLE_ERROR_GOTO("I/O error reading DAT file", 
                           AFFY_ERROR_IO, 
                           err, 
                           cleanup);

  newdat->pixel_width = strtoul(tmpbuf + 4, NULL, 10);

  if (affy_readchars(input, tmpbuf, 8) != 0)
    AFFY_HANDLE_ERROR_GOTO("I/O error reading DAT file", 
                           AFFY_ERROR_IO, 
                           err, 
                           cleanup);

  newdat->pixel_height = strtoul(tmpbuf + 4, NULL, 10);

  if (affy_readchars(input, tmpbuf, 7) != 0)
    AFFY_HANDLE_ERROR_GOTO("I/O error reading DAT file", 
                           AFFY_ERROR_IO, 
                           err, 
                           cleanup);

  newdat->scanspeed = strtoul(tmpbuf + 3, NULL, 10);

  if (affy_readchars(input, tmpbuf, 8) != 0)
    AFFY_HANDLE_ERROR_GOTO("I/O error reading DAT file", 
                           AFFY_ERROR_IO, 
                           err, 
                           cleanup);

  newdat->temperature = strtod(tmpbuf, NULL);

  if (affy_readchars(input, tmpbuf, 5) != 0)
    AFFY_HANDLE_ERROR_GOTO("I/O error reading DAT file", 
                           AFFY_ERROR_IO, 
                           err, 
                           cleanup);

  newdat->laser_power = strtod(tmpbuf, NULL);

  if (affy_readchars(input, tmpbuf, sizeof(newdat->timestamp)) != 0)
    AFFY_HANDLE_ERROR_GOTO("I/O error reading DAT file", 
                           AFFY_ERROR_IO, 
                           err, 
                           cleanup);

  strcpy(newdat->timestamp, tmpbuf);

  /* Variable-size subfields, read them as one unit. */
  if (affy_readchars(input, tmpbuf, 221) != 0)
    AFFY_HANDLE_ERROR_GOTO("I/O error reading DAT file", 
                           AFFY_ERROR_IO, 
                           err, 
                           cleanup);

  if ((walk = strstr(tmpbuf, FIELD_START)) == NULL)
    AFFY_HANDLE_ERROR_GOTO("bad DAT file magic", 
                           AFFY_ERROR_BADFORMAT, 
                           err, 
                           cleanup);
  else
    *walk++ = '\0';

  newdat->scannerid = h_strdup(tmpbuf);
  if (newdat->scannerid == NULL)
    AFFY_HANDLE_ERROR_GOTO("strdup failed", AFFY_ERROR_OUTOFMEM, err, cleanup);
  hattach(newdat->scannerid, newdat);
  

  if ((walk = strstr(walk, FIELD_START)) == NULL)
    AFFY_HANDLE_ERROR_GOTO("bad DAT file magic", 
                           AFFY_ERROR_BADFORMAT, 
                           err, 
                           cleanup);
  else
    save = walk + strlen(FIELD_START);

  if ((walk = strstr(walk, ".1sq")) == NULL)
    AFFY_HANDLE_ERROR_GOTO("bad DAT file magic", 
                           AFFY_ERROR_BADFORMAT, 
                           err, 
                           cleanup);
  else
    *walk = '\0';

  newdat->probe_array_type = h_strdup(save);
  if (newdat->probe_array_type == NULL)
    AFFY_HANDLE_ERROR_GOTO("strdup failed", AFFY_ERROR_OUTOFMEM, err, cleanup);
  hattach(newdat->probe_array_type, newdat);

  if (affy_readmulti(input, "%2fl%dl%9hl",
                     (void *)&newdat->avg_dc_offset,
                     (void *)&newdat->std_dev_dc_offset,
                     (void *)&newdat->numsamples_dc_offset,
                     (void *)&newdat->grid_ul.x,
                     (void *)&newdat->grid_ul.y,
                     (void *)&newdat->grid_ur.x,
                     (void *)&newdat->grid_ur.y,
                     (void *)&newdat->grid_lr.x,
                     (void *)&newdat->grid_lr.y,
                     (void *)&newdat->grid_ll.x,
                     (void *)&newdat->grid_ll.y,
                     (void *)&newdat->cellmargin) <= 0)
    AFFY_HANDLE_ERROR_GOTO("I/O error reading DAT file", 
                           AFFY_ERROR_IO, 
                           err, 
                           cleanup);

  if (affy_readchars(input, tmpbuf, 155) != 0)
    AFFY_HANDLE_ERROR_GOTO("I/O error reading DAT file", 
                           AFFY_ERROR_IO, 
                           err, 
                           cleanup);

  newdat->experiment_name = h_strdup(tmpbuf);
  if (newdat->experiment_name == NULL)
    AFFY_HANDLE_ERROR_GOTO("strdup failed", AFFY_ERROR_OUTOFMEM, err, cleanup);
  hattach(newdat->experiment_name, newdat);

  /* Allocate pixelmap memory. */
  newdat->pixels.data =
    (unsigned int **)h_suballoc(newdat, 
                                sizeof(unsigned int *) 
                                * newdat->pixels.numrows);
  if (newdat->pixels.data == NULL)
    AFFY_HANDLE_ERROR_GOTO("malloc failed", AFFY_ERROR_OUTOFMEM, err, cleanup);

  newdat->pixels.data[0] =
    (unsigned int *)h_suballoc(newdat, 
                               sizeof(unsigned int) * newdat->numpixels);
  if (newdat->pixels.data[0] == NULL)
    AFFY_HANDLE_ERROR_GOTO("malloc failed", AFFY_ERROR_OUTOFMEM, err, cleanup);

  /* Fixup the pointers for the "lines" dimension of the pixelmap. */
  for (i = 1; i < newdat->pixels.numrows; i++)
    newdat->pixels.data[i] =
      newdat->pixels.data[i - 1] + newdat->pixels.numcols;

  if (affy_readmulti(input, "%*hl", 
                     newdat->numpixels, (void *)newdat->pixels.data[0]) <= 0)
    AFFY_HANDLE_ERROR_GOTO("I/O error reading DAT file", 
                           AFFY_ERROR_IO, 
                           err, 
                           cleanup);

  fclose(input);

  return (newdat);

cleanup:
  fclose(input);
  affy_free_dat_file(newdat);

  return (NULL);
}
