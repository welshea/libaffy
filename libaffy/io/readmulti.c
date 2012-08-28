
/**************************************************************************
 *
 * Filename:  readmulti.c
 * 
 * Purpose:   Convenience routine for allowing scanf-style reading of
 *            many input values from a file stream.
 *
 * Creation:  01/09/08
 *
 * Author:    Andrew M. Hoerter
 *
 * Copyright: Copyright (C) 2008, Moffitt Cancer Center.
 *            All rights reserved.
 *
 * Update History
 * --------------
 * 01/09/08: Creation (AMH)
 *
 **************************************************************************/

#include <affy.h>
    
/* 
 * affy_readmulti(): "Vectorized" read function to conveniently read
 *                   multiple values in a concise fashion.
 *
 * This function provides a somewhat similar interface to fscanf()
 * with a couple of twists.  The format string contained in fmt
 * directs a series of binary values to be read from fp.  The
 * structure of the format string is as follows:
 *
 *    %[* | N]{c,h,d,f,D}[l,b]...
 *
 * or, %x
 *
 * Each piece of the format string consists of an optional "repeat"
 * specifier, followed by a size specifier (c/h/d/f/D), followed by an
 * optional endianness specifier.  c corresponds to 8 bits, h to 16, d
 * to 32, f to 64.  Endianness may be big or little, with the
 * corresponding swapping performed as needed.  If the endianness
 * specifier is omitted, no swapping is performed whatsoever.  'D' is
 * a special size specifier; internally, libaffy uses doubles (usually
 * 64-bit in size).  However, it is commonly necessary to read 32-bit
 * values from disk which are then widened to 64-bits for internal
 * use.  The 'D' specifier reads a 32-bit quantity from disk into a
 * float, widens it to a double, and stores it in the desired location.
 *
 * A repetition specifier may be a '*', indicating that an extra
 * integer argument will be consumed to specify how many times to read
 * to the supplied destination address; or it may be a number, in
 * which case that many destination address arguments will be consumed
 * for each read.  This is to allow for compile-time vs. run-time
 * flexibility.
 *
 * The special size specifier 'x' indicates that a seek should be
 * performed relative to the current file position.  A single long
 * argument is consumed from the arg list as an offset.
 *
 * Items in the format string are separated by the % symbol.  Remember
 * that unlike the scanf() family of functions, affy_readmulti() deals
 * with binary files so the lexical arrangement of its format string is
 * not relevant, the '%' convention was only chosen for familiarity.
 * Whitespace, if present in the format string, is ignored.
 *
 * The variable arguments to this function consist of a series of
 * unsigned (long) integers and void pointers.  For each piece of the
 * format string, if a repeat flag or repeat count is present, the
 * function will perform a number of reads of the given size.  The
 * number is determined either by a single additional variable
 * argument (in the case of a repeat flag '*'), or by the number
 * present in the format string directive (e.g. '%20dl').  In the case
 * of a numerical repeat count, each read consumes a destination
 * pointer argument.
 *
 * Then, a (void *) pointing to the destination for the data is
 * consumed.
 *
 * The integer return value communicates two pieces of information.
 * First, its absolute value is equal to the number of format string
 * directives which were successfully processed (note that a single
 * directive repeated a number of times still only contributes 1 to
 * this total).  Second, if the value is less than or equal to 0, an
 * error occurred while reading (or an empty format string was provided).
 */
int affy_readmulti(FILE *fp, const char *fmt, ...)
{
  va_list     ap;
  int         numread = 0;
  const char *p_cur;

  assert(fp  != NULL);
  assert(fmt != NULL);

  va_start(ap, fmt);

  p_cur = fmt;

  while (*p_cur != '\0')
  {
    unsigned int i, readsz, repeat_count = 1, repeat_flag = 0, widen_float = 0;
    char         read_type, *p_destarg = NULL;
    int        (*read_function)(FILE *, void *);

    /* Advance to the next % marker. */
    if (*p_cur++ != '%')
      continue;

    /* Handle seeks. */
    if (*p_cur == 'x')
    {
      long ofs = va_arg(ap, long);

      if (fseek(fp, ofs, SEEK_CUR) != 0)
        return (-numread);
      else
        p_cur++;
      
      numread++;
      continue;
    }

    /* Handle repeat count marker, if present. */
    if (*p_cur == '*')
    {
      repeat_count = va_arg(ap, int);
      p_cur++;
    }
    else if (isdigit(*p_cur))
    {
      /* Numerical repeat specifier. */
      repeat_flag = 1;

      repeat_count = strtoul(p_cur, (char **) &p_cur, 10);
    }

    /* Select the type of read to be performed. */
    read_type = *p_cur++;

    switch (read_type)
    {
      case 'c':
        readsz = 1;

        read_function = affy_read8;
        break;

      case 'h':
        readsz = 2;
        
        if (*p_cur == 'l')
          read_function = affy_read16_le;
        else if (*p_cur == 'b')
          read_function = affy_read16_be;
        else
          read_function = affy_read16;
        break;

      case 'D':
        widen_float = 1;
      case 'd':
        readsz = 4;

        if (*p_cur == 'l')
          read_function = affy_read32_le;
        else if (*p_cur == 'b')
          read_function = affy_read32_be;
        else
          read_function = affy_read32;
        break;

      case 'f':
        readsz = 8;

        if (*p_cur == 'l')
          read_function = affy_read64_le;
        else if (*p_cur == 'b')
          read_function = affy_read64_be;
        else
          read_function = affy_read64;
        break;
        
      default:
        return (-numread);
    }
    
    /* Perform the actual read. */
    for (i = 0; i < repeat_count; i++)
    {
      affy_float32  tmp;
      affy_float64 *p_dest64;
      void         *p_dest;

      if ((p_destarg == NULL) || (repeat_flag == 1))
        p_destarg = va_arg(ap, void *);

      p_dest = (widen_float) ? (void *)(&tmp) : (void *)p_destarg;

      if (read_function(fp, p_dest) != 0)
        return (-numread);

      if (widen_float)
      {
        p_dest64  = (affy_float64 *)p_destarg;
        *p_dest64 = tmp;
      }

      if (repeat_flag == 0)
        p_destarg += readsz;
    }

    numread++;
  }

  /* Finished with no errors. */
  return (numread);
}
