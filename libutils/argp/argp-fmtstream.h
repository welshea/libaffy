/* Word-wrapping and line-truncating streams.
   Copyright (C) 1997, 2006 Free Software Foundation, Inc.
   This file is part of the GNU C Library.
   Written by Miles Bader <miles@gnu.ai.mit.edu>.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along
   with this program; if not, write to the Free Software Foundation,
   Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA. */

/* This package emulates glibc `line_wrap_stream' semantics for systems that
   don't have that.  If the system does have it, it is just a wrapper for
   that.  This header file is only used internally while compiling argp, and
   shouldn't be installed.  */

#ifndef _ARGP_FMTSTREAM_H
#define _ARGP_FMTSTREAM_H

#include <stdio.h>
#include <string.h>

#include "utils.h"

#ifndef __attribute__
/* This feature is available in gcc versions 2.5 and later.  */
# if __GNUC__ < 2 || (__GNUC__ == 2 && __GNUC_MINOR__ < 5) || __STRICT_ANSI__
#  define __attribute__(Spec) /* empty */
# endif
/* The __-protected variants of `format' and `printf' attributes
   are accepted by gcc versions 2.6.4 (effectively 2.7) and later.  */
# if __GNUC__ < 2 || (__GNUC__ == 2 && __GNUC_MINOR__ < 7) || __STRICT_ANSI__
#  define __format__ format
#  define __printf__ printf
# endif
#endif


/* Guess we have to define our own version.  */

struct argp_fmtstream
{
  FILE *stream;			/* The stream we're outputting to.  */

  size_t lmargin, rmargin;	/* Left and right margins.  */
  int wmargin;		/* Margin to wrap to, or -1 to truncate.  */

  /* Point in buffer to which we've processed for wrapping, but not output.  */
  size_t point_offs;
  /* Output column at POINT_OFFS, or -1 meaning 0 but don't add lmargin.  */
  int point_col;

  char *buf;			/* Output buffer.  */
  char *p;			/* Current end of text in BUF. */
  char *end;			/* Absolute end of BUF.  */
};

typedef struct argp_fmtstream *argp_fmtstream_t;


int 
argp_fmtstream_printf (struct argp_fmtstream *fs, const char *fmt, ...);

void argp_fmtstream_free (argp_fmtstream_t fs);
argp_fmtstream_t
argp_make_fmtstream (FILE *stream,
                     size_t lmargin, size_t rmargin, int wmargin);

/* Access macros for various bits of state.  */
#define argp_fmtstream_lmargin(__fs) ((__fs)->lmargin)
#define argp_fmtstream_rmargin(__fs) ((__fs)->rmargin)
#define argp_fmtstream_wmargin(__fs) ((__fs)->wmargin)

/* Internal routines.  */
 void argp_fmtstream_update (argp_fmtstream_t fs);
 int argp_fmtstream_ensure (argp_fmtstream_t fs, size_t amount);


#ifndef ARGP_FS_EI
#define ARGP_FS_EI static INLINE
#endif

ARGP_FS_EI size_t
argp_fmtstream_write (argp_fmtstream_t fs,
			const char *str, size_t len)
{
  if (fs->p + len <= fs->end || argp_fmtstream_ensure (fs, len))
    {
      memcpy (fs->p, str, len);
      fs->p += len;
      return len;
    }
  else
    return 0;
}

ARGP_FS_EI int
argp_fmtstream_puts (argp_fmtstream_t fs, const char *str)
{
  size_t len = strlen (str);
  if (len)
    {
      size_t wrote = argp_fmtstream_write (fs, str, len);
      return wrote == len ? 0 : -1;
    }
  else
    return 0;
}

ARGP_FS_EI int
argp_fmtstream_putc (argp_fmtstream_t fs, int ch)
{
  if (fs->p < fs->end || argp_fmtstream_ensure (fs, 1))
    return *fs->p++ = ch;
  else
    return EOF;
}

/* Set FS's left margin to LMARGIN and return the old value.  */
ARGP_FS_EI size_t
argp_fmtstream_set_lmargin (argp_fmtstream_t fs, size_t lmargin)
{
  size_t old;
  if ((size_t) (fs->p - fs->buf) > fs->point_offs)
    argp_fmtstream_update (fs);
  old = fs->lmargin;
  fs->lmargin = lmargin;
  return old;
}

/* Set FS's right margin to RMARGIN and return the old value.  */
ARGP_FS_EI size_t
argp_fmtstream_set_rmargin (argp_fmtstream_t fs, size_t rmargin)
{
  size_t old;
  if ((size_t) (fs->p - fs->buf) > fs->point_offs)
    argp_fmtstream_update (fs);
  old = fs->rmargin;
  fs->rmargin = rmargin;
  return old;
}

/* Set FS's wrap margin to WMARGIN and return the old value.  */
ARGP_FS_EI size_t
argp_fmtstream_set_wmargin (argp_fmtstream_t fs, size_t wmargin)
{
  size_t old;
  if ((size_t) (fs->p - fs->buf) > fs->point_offs)
    argp_fmtstream_update (fs);
  old = fs->wmargin;
  fs->wmargin = wmargin;
  return old;
}

/* Return the column number of the current output point in FS.  */
ARGP_FS_EI size_t
argp_fmtstream_point (argp_fmtstream_t fs)
{
  if ((size_t) (fs->p - fs->buf) > fs->point_offs)
    argp_fmtstream_update (fs);
  return fs->point_col >= 0 ? fs->point_col : 0;
}

#endif /* argp-fmtstream.h */
