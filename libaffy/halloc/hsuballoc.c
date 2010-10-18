/**************************************************************************
 *
 * Filename:  hsuballoc.c
 * 
 * Purpose:   Wrapper routines for halloc to implement a pooled malloc().
 *
 * Creation:  09/20/2010
 *
 * Author:    Andrew M. Hoerter
 *
 * Copyright: Copyright (C) 2010, Moffitt Cancer Center.
 *            All rights reserved.
 *
 * Update History
 * --------------
 * 09/20/10: Initial creation (AMH)
 *
 **************************************************************************/

#include <assert.h>

#include "halloc.h"

/*
 * Same as h_malloc(), except it automatically attaches the new storage to
 * parent and refuses to malloc 0 bytes.
 */
void *h_suballoc(void *parent, size_t len)
{
  void *new;

  assert(len > 0);

  new = h_malloc(len);
  if (new != NULL)
    hattach(new, parent);
  
  return (new);
}

/*
 * Same as h_calloc(), except it automatically attaches the new storage to
 * parent and refuses to malloc 0 bytes.
 */
void *h_subcalloc(void *parent, size_t nmemb, size_t sz)
{
  void *new;

  assert(sz    > 0);
  assert(nmemb > 0);

  new = h_calloc(nmemb, sz);
  if (new != NULL)
    hattach(new, parent);
  
  return (new);
}
