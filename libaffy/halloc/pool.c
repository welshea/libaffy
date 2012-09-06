/**************************************************************************
 *
 * Filename:  pool.c
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

#include "affy.h"
#include "halloc.h"

affy_mempool *affy_pool_create(void)
{
  return (h_malloc(sizeof(affy_mempool)));
}

void affy_pool_destroy(affy_mempool *pool)
{
  h_free(pool);
}

void *affy_pool_alloc(affy_mempool *pool, size_t len)
{
  void *new;

  assert(pool != NULL);
  assert(len > 0);
  
  new = h_malloc(len);
  if (new != NULL)
    hattach(new, pool);
  
  return (new);
}
 
void affy_pool_free(void *mem)
{
  h_free(mem);
}

void affy_pool_attach(affy_mempool *child, affy_mempool *parent)
{
  assert(child  != NULL);
  assert(parent != NULL);

  hattach(child, parent);
}
