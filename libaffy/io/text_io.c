
/**************************************************************************
 *
 * Filename:  text_io.c
 *
 * Purpose:   Routines for dealing with plain text I/O.
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
 * 03/18/08: Thread safety improvements (AMH)
 * 09/20/10: Pooled memory allocator (AMH)
 * 09/05/23: modifications to replace fgets() with EOL-safe function
 *
 **************************************************************************/

#include <affy.h>
#include <utils.h>

/*==================================================================
  Simple line reader with undo
*===================================================================*/

/*
 * Free up an AFFY_TEXTIO structure.
 */
void affy_textio_free(AFFY_TEXTIO *tf)
{
  assert(tf != NULL);
  
  /* free buffer, which was allocated with malloc() instead of h_free() */
  if (tf->buf)
    free(tf->buf);
  tf->buf         = NULL;
  tf->max_buf_len = 0;

  h_free(tf);
}

/*
 * Initialize an AFFY_TEXTIO structure and return it.  This must be
 * done before using the text I/O routines.
 */
AFFY_TEXTIO *affy_textio_init(FILE *fp, AFFY_ERROR *err)
{
  AFFY_TEXTIO *tf;
  
  assert(fp != NULL);

  tf = h_malloc(sizeof(AFFY_TEXTIO));
  if (tf == NULL)
    AFFY_HANDLE_ERROR("malloc failed", AFFY_ERROR_OUTOFMEM, err, NULL);

  tf->skip_read   = false;
  tf->fp          = fp;
  tf->buf         = calloc(MAXBUF, sizeof(char));   /* don't use h_alloc() */
  tf->max_buf_len = MAXBUF;   /* needed for fgets_strip_realloc() */

  if (tf->buf == NULL)
  {
    h_free(tf);
    
    AFFY_HANDLE_ERROR("malloc failed", AFFY_ERROR_OUTOFMEM, err, NULL);
  }

  return (tf);
}

/* Read a single line from the data file, trim it, and return it */
char *affy_textio_get_next_line(AFFY_TEXTIO *tf)
{
  char *p = NULL;

  assert(tf != NULL);

  while (1)
  {
    /* A special flag allows us to unget a line */
    if (tf->skip_read)
    {
      tf->skip_read = false;
    }
    else if (fgets_strip_realloc(&(tf->buf),
                                 &(tf->max_buf_len), tf->fp) == NULL)
    {
      return (NULL);
    }

    /* Trim the whitespace */
    p = trim(tf->buf);

    /* If not empty, return it */
    if (strlen(p) > 0)
      return (p);
  }

  /* Not reached */
}

void affy_textio_unget_next_line(AFFY_TEXTIO *tf)
{
  assert(tf != NULL);

  tf->skip_read = true;
}

void affy_textio_reset_next_line(AFFY_TEXTIO *tf)
{
  assert(tf != NULL);

  tf->skip_read = false;
}

void affy_textio_skip_to_next_header(AFFY_TEXTIO *tf)
{
  char *s;

  assert(tf != NULL);

  while ((s = affy_textio_get_next_line(tf)) != NULL)
  {
    if (*s == '[')
    {
      affy_textio_unget_next_line(tf);

      return;
    }
  }
}
