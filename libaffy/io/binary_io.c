
/**************************************************************************
 *
 * Filename:  binary_io.c
 * 
 * Purpose:   Routines for reading/writing discrete data types from a
 *            binary file.  Endianness is accounted for and corrected.
 *
 * Creation:  04/04/05
 *
 * Author:    Andrew M. Hoerter
 *
 * Copyright: Copyright (C) 2007, Moffitt Cancer Center.
 *            All rights reserved.
 *
 * Update History
 * --------------
 * 04/04/05: Creation (AMH)
 * 10/05/07: Slight revamp; add support for both LE/BE I/O (AMH)
 * 01/10/08: Convert the native-endian read macros to functions (AMH)
 * 01/23/08: Remove some obsolete native-type-size I/O functions (AMH)
 * 10/08/10: Add write routines (AMH)
 *
 **************************************************************************/

#include "affy.h"
#include "endian_config.h"

/* Note that the endianness #define comes from the libutils package. */

/* Endian-swapping macros (x = (void *)). */
#define AFFY_SWAP16(x) *(affy_uint16 *)(x) =		               \
                         (((*(affy_uint16 *)(x)) & 0x00ff) << 8) |     \
			 (((*(affy_uint16 *)(x)) >> 8))

#define AFFY_SWAP32(x) *(affy_uint32 *)(x) =			       \
                         ((*(affy_uint32 *)(x) << 24) |	      	       \
                         ((*(affy_uint32 *)(x) << 8) & 0x00ff0000UL) | \
                         ((*(affy_uint32 *)(x) >> 8) & 0x0000ff00UL) | \
			 (*(affy_uint32 *)(x) >> 24))

#define AFFY_SWAPFL(x) AFFY_SWAP32(x)
#define AFFY_SWAPDBL(x) AFFY_SWAP64(x)

static void AFFY_SWAP64(void *buf)
{
  affy_uint8 *byte = (affy_uint8 *)buf;
  affy_uint8 tmp;

  tmp = byte[0];
  byte[0] = byte[7];
  byte[7] = tmp;

  tmp = byte[1];
  byte[1] = byte[6];
  byte[6] = tmp;

  tmp = byte[2];
  byte[2] = byte[5];
  byte[5] = tmp;

  tmp = byte[3];
  byte[3] = byte[4];
  byte[4] = tmp;
}

/* 
 * affy_read8(): Read an unsigned byte off the file stream. 
 *
 * Inputs: *input is an open file input stream, *buf points to allocated
 *          storage.
 * Outputs: 0 on success, -1 on failure (error or EOF).
 * Side effects: None.
 */
int affy_read8(FILE *input, void *buf)
{
  assert(input != NULL);
  assert(buf   != NULL);

  if (fread(buf, 1, 1, input) < 1)
    return (-1);

  /* For bytes, that's it. */
  return (0);
}

/* 
 * affy_write8(): Write an unsigned byte to the file stream. 
 *
 * Inputs: *output is an open file output stream, *buf points to the desired 
 *         value.
 * Outputs: 0 on success, -1 on failure (error or EOF).
 * Side effects: None.
 */
int affy_write8(FILE *output, void *buf)
{
  assert(output != NULL);
  assert(buf    != NULL);

  if (fwrite(buf, 1, 1, output) < 1)
    return (-1);

  /* For bytes, that's it. */
  return (0);
}

/* 
 * affy_read16_le(): Read a 16-bit little-endian integer off the file stream. 
 *
 * Inputs: *input is an open file input stream, *buf points to allocated
 *          storage.
 * Outputs: 0 on success, -1 on failure (error or EOF).
 * Side effects: None.
 */
int affy_read16_le(FILE *input, void *buf)
{
  assert(input != NULL);
  assert(buf   != NULL);

  if (fread(buf, 2, 1, input) < 1)
    return (-1);

#ifdef LIBUTILS_BIG_ENDIAN
  AFFY_SWAP16(buf);
#endif

  return (0);
}

/* 
 * affy_write16_le(): Write a 16-bit little-endian integer to the file stream.
 *
 * Inputs: *output is an open file output stream, *buf points to the desired 
 *         value.
 * Outputs: 0 on success, -1 on failure (error or EOF).
 * Side effects: None.
 */
int affy_write16_le(FILE *output, void *buf)
{
  affy_uint16 tmp = *((affy_uint16 *)buf);

  assert(output != NULL);
  assert(buf    != NULL);

#ifdef LIBUTILS_BIG_ENDIAN
  AFFY_SWAP16(&tmp);
#endif

  if (fwrite((void *)&tmp, 2, 1, output) < 1)
    return (-1);

  return (0);
}

/* 
 * affy_read16_be(): Read a 16-bit big-endian integer off the file stream. 
 *
 * Inputs: *input is an open file input stream, *buf points to allocated
 *          storage.
 * Outputs: 0 on success, -1 on failure (error or EOF).
 * Side effects: None.
 */
int affy_read16_be(FILE *input, void *buf)
{
  assert(input != NULL);
  assert(buf   != NULL);

  if (fread(buf, 2, 1, input) < 1)
    return (-1);

#ifndef LIBUTILS_BIG_ENDIAN
  AFFY_SWAP16(buf);
#endif

  return (0);
}

/* 
 * affy_write16_be(): Write a 16-bit big-endian integer to the file stream.
 *
 * Inputs: *output is an open file output stream, *buf points to the desired 
 *         value.
 * Outputs: 0 on success, -1 on failure (error or EOF).
 * Side effects: None.
 */
int affy_write16_be(FILE *output, void *buf)
{
  affy_uint16 tmp = *((affy_uint16 *)buf);

  assert(output != NULL);
  assert(buf    != NULL);

#ifdef LIBUTILS_LITTLE_ENDIAN
  AFFY_SWAP16(&tmp);
#endif

  if (fwrite((void *)&tmp, 2, 1, output) < 1)
    return (-1);

  return (0);
}

/* 
 * affy_read32_le(): Read a 32-bit little-endian integer or float off the 
 *                   file stream. 
 *
 * Inputs: *input is an open file input stream, *buf points to allocated
 *          storage.
 * Outputs: 0 on success, -1 on failure (error or EOF).
 * Side effects: None.
 */
int affy_read32_le(FILE *input, void *buf)
{
  assert(input != NULL);
  assert(buf   != NULL);

  if (fread(buf, 4, 1, input) < 1)
    return (-1);

#ifdef LIBUTILS_BIG_ENDIAN
  AFFY_SWAP32(buf);
#endif

  return (0);
}

/* 
 * affy_write32_le(): Write a 32-bit little-endian integer or float to the 
 *                   file stream. 
 *
 * Inputs: *output is an open file output stream, *buf points to the desired 
 *         value.
 * Outputs: 0 on success, -1 on failure (error or EOF).
 * Side effects: None.
 */
int affy_write32_le(FILE *output, void *buf)
{
  assert(output != NULL);
  assert(buf    != NULL);

#ifdef LIBUTILS_BIG_ENDIAN
  AFFY_SWAP32(&tmp);
#endif

  if (fwrite(buf, 4, 1, output) < 1)
    return (-1);

  return (0);
}

/* 
 * affy_read32_be(): Read a 32-bit big-endian integer or float off the 
 *                   file stream. 
 *
 * Inputs: *input is an open file input stream, *buf points to allocated
 *          storage.
 * Outputs: 0 on success, -1 on failure (error or EOF).
 * Side effects: None.
 */
int affy_read32_be(FILE *input, void *buf)
{
  assert(input != NULL);
  assert(buf   != NULL);

  if (fread(buf, 4, 1, input) < 1)
    return (-1);

#ifndef LIBUTILS_BIG_ENDIAN
  AFFY_SWAP32(buf);
#endif

  return (0);
}

/* 
 * affy_write32_be(): Write a 32-bit big-endian integer or float to the 
 *                   file stream. 
 *
 * Inputs: *output is an open file output stream, *buf points to the desired 
 *         value.
 * Outputs: 0 on success, -1 on failure (error or EOF).
 * Side effects: None.
 */
int affy_write32_be(FILE *output, void *buf)
{
  affy_uint32 tmp = *((affy_uint32 *)buf);

  assert(output != NULL);
  assert(buf    != NULL);

#ifdef LIBUTILS_LITTLE_ENDIAN
  AFFY_SWAP32(&tmp);
#endif

  if (fwrite(buf, 4, 1, output) < 1)
    return (-1);

  return (0);
}

/* 
 * affy_read64_le(): Read a 64-bit little-endian integer or float off the 
 *                   file stream. 
 *
 * Inputs: *input is an open file input stream, *buf points to allocated
 *          storage.
 * Outputs: 0 on success, -1 on failure (error or EOF).
 * Side effects: None.
 */
int affy_read64_le(FILE *input, void *buf)
{
  assert(input != NULL);
  assert(buf   != NULL);

  if (fread(buf, 8, 1, input) < 1)
    return (-1);

#ifdef LIBUTILS_BIG_ENDIAN
  AFFY_SWAP64(buf);
#endif

  return (0);
}

/* 
 * affy_write64_le(): Write a 64-bit little-endian integer or float to the 
 *                    file stream. 
 *
 * Inputs: *output is an open file output stream, *buf points to the desired 
 *         value.
 * Outputs: 0 on success, -1 on failure (error or EOF).
 * Side effects: None.
 */
int affy_write64_le(FILE *output, void *buf)
{
  affy_uint8 tmp[8];

  assert(output != NULL);
  assert(buf    != NULL);

  memcpy(&tmp, buf, 8);

#ifdef LIBUTILS_BIG_ENDIAN
  AFFY_SWAP64(&tmp);
#endif

  if (fwrite((void *)&tmp, 8, 1, output) < 1)
    return (-1);

  return (0);
}

/* 
 * affy_read64_be(): Read a 64-bit big-endian integer or float off the 
 *                   file stream. 
 *
 * Inputs: *input is an open file input stream, *buf points to allocated
 *          storage.
 * Outputs: 0 on success, -1 on failure (error or EOF).
 * Side effects: None.
 */
int affy_read64_be(FILE *input, void *buf)
{
  assert(input != NULL);
  assert(buf   != NULL);

  if (fread(buf, 8, 1, input) < 1)
    return (-1);

#ifndef LIBUTILS_BIG_ENDIAN
  AFFY_SWAP64(buf);
#endif

  return (0);
}

/* 
 * affy_write64_be(): Write a 64-bit big-endian integer or float to the 
 *                    file stream. 
 *
 * Inputs: *output is an open file output stream, *buf points to the desired 
 *         value.
 * Outputs: 0 on success, -1 on failure (error or EOF).
 * Side effects: None.
 */
int affy_write64_be(FILE *output, void *buf)
{
  affy_uint8 tmp[8];

  assert(output != NULL);
  assert(buf    != NULL);

  memcpy(&tmp, buf, 8);

#ifdef LIBUTILS_LITTLE_ENDIAN
  AFFY_SWAP64(&tmp);
#endif

  if (fwrite((void *)&tmp, 8, 1, output) < 1)
    return (-1);

  return (0);
}

/* 
 * affy_read64(): Read a 64-bit native-endian integer or float off the 
 *                file stream. 
 *
 * Inputs: *input is an open file input stream, *buf points to allocated
 *          storage.
 * Outputs: 0 on success, -1 on failure (error or EOF).
 * Side effects: None.
 */
int affy_read64(FILE *input, void *buf)
{
  assert(input != NULL);
  assert(buf   != NULL);

  if (fread(buf, 8, 1, input) < 1)
    return (-1);

  return (0);
}

/* 
 * affy_write64(): Write a 64-bit native-endian integer or float to the 
 *                 file stream. 
 *
 * Inputs: *output is an open file output stream, *buf points to the 
 *         desired value.
 * Outputs: 0 on success, -1 on failure (error or EOF).
 * Side effects: None.
 */
int affy_write64(FILE *output, void *buf)
{
  assert(output != NULL);
  assert(buf    != NULL);

  if (fwrite(buf, 8, 1, output) < 1)
    return (-1);

  return (0);
}

/* 
 * affy_read32(): Read a 32-bit native-endian integer or float off the 
 *                file stream. 
 *
 * Inputs: *input is an open file input stream, *buf points to allocated
 *          storage.
 * Outputs: 0 on success, -1 on failure (error or EOF).
 * Side effects: None.
 */
int affy_read32(FILE *input, void *buf)
{
  assert(input != NULL);
  assert(buf   != NULL);

  if (fread(buf, 4, 1, input) < 1)
    return (-1);

  return (0);
}

/* 
 * affy_write32(): Write a 32-bit native-endian integer or float to the 
 *                 file stream. 
 *
 * Inputs: *output is an open file output stream, *buf points to the 
 *         desired value.
 * Outputs: 0 on success, -1 on failure (error or EOF).
 * Side effects: None.
 */
int affy_write32(FILE *output, void *buf)
{
  assert(output != NULL);
  assert(buf    != NULL);

  if (fwrite(buf, 4, 1, output) < 1)
    return (-1);

  return (0);
}
/* 
 * affy_read16(): Read a 16-bit native-endian integer or float off the 
 *                file stream. 
 *
 * Inputs: *input is an open file input stream, *buf points to allocated
 *          storage.
 * Outputs: 0 on success, -1 on failure (error or EOF).
 * Side effects: None.
 */
int affy_read16(FILE *input, void *buf)
{
  assert(input != NULL);
  assert(buf   != NULL);

  if (fread(buf, 2, 1, input) < 1)
    return (-1);

  return (0);
}

/* 
 * affy_write16(): Write a 16-bit native-endian integer or float to the 
 *                 file stream. 
 *
 * Inputs: *output is an open file output stream, *buf points to the 
 *         desired value.
 * Outputs: 0 on success, -1 on failure (error or EOF).
 * Side effects: None.
 */
int affy_write16(FILE *output, void *buf)
{
  assert(output != NULL);
  assert(buf    != NULL);

  if (fwrite(buf, 2, 1, output) < 1)
    return (-1);

  return (0);
}

/* 
 * affy_readchars(): Read a series of characters off the file stream. 
 *
 * Inputs: *input is an open file input stream, *buf points to allocated
 *          storage large enough to hold the desired number of bytes.
 *          At most numbytes-1 characters will be read, and a NUL byte
 *          will be appended.
 * Outputs: 0 on success, -1 on failure (error or EOF).
 * Side effects: None.
 *
 * Comments: Provided for ease of reading character strings.  Almost
 *           equivalent to fgets() but ignoring newlines.
 */
int affy_readchars(FILE *input, char *buf, size_t numbytes)
{
  size_t i;
  int    next;

  assert(input != NULL);
  assert(buf   != NULL);

  for (i = 0; i < (numbytes - 1); i++)
  {
    next = fgetc(input);

    if (next == EOF)
    {
      buf[i] = '\0';
      return (-1);
    }
    else
    {
      buf[i] = next;
    }
  }

  buf[i] = '\0';

  return (0);
}

/* 
 * affy_writechars(): Write a series of NUL terminated characters to the 
 *                    file stream. 
 *
 * Inputs: *output is an open file output stream, *buf points to the 
 *         desired string.
 * Outputs: 0 on success, -1 on failure (error or EOF).
 * Side effects: None.
 *
 * Comments:
 */
int affy_writechars(FILE *output, char *buf)
{
  char *c = buf;

  assert(output != NULL);
  assert(buf    != NULL);

  while (*c != '\0')
  {
    if (fputc(*c++, output) == EOF)
      return (-1);
  }

  return (0);
}
