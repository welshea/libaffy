
/**************************************************************************
 *
 * Filename:  load_calvin_file.c
 *
 * Purpose:   Parse a Calvin (Command Console "generic") file and initialize 
 *            the appropriate structures.  This routine makes no attempt to
 *            interpret the data or produce a DAT/CEL file object, that is
 *            handled by other functions.
 *
 * Creation:  10/01/07
 *
 * Author:    Andrew M. Hoerter
 *
 * Copyright: Copyright (C) 2007, Moffitt Cancer Center.
 *            All rights reserved.
 *
 * Update History
 * --------------
 * 10/01/07: Creation (AMH)
 * 04/18/08: New error handling scheme (AMH)
 * 04/23/08: Refactor in terms of new lower level I/O routines (AMH)
 * 12/01/08: More refactoring (AMH)
 * 09/20/10: Pooled memory allocator (AMH)
 *
 **************************************************************************/

#include <affy.h>

AFFY_CALVIN_CONTAINER *affy_calvin_load_container(AFFY_CALVINIO *cio, 
                                                  AFFY_ERROR *err)
{
  AFFY_CALVIN_CONTAINER  *cc;
  AFFY_CALVIN_FILEHEADER *fh;
  AFFY_CALVIN_DATAHEADER *dh;
  affy_uint32             i;

  assert(cio != NULL);

  cc = h_malloc(sizeof(AFFY_CALVIN_CONTAINER));
  if (cc == NULL)
    AFFY_HANDLE_ERROR("malloc failed", AFFY_ERROR_OUTOFMEM, err, NULL);

  /* Perform basic structure initialization */
  cc->file_header->num_datagroups    = 0;
  cc->data_groups                    = NULL;

  /* Get the file header. */
  cc->file_header = fh = affy_calvin_get_file_metadata(cio, err);
  AFFY_CHECK_ERROR_GOTO(err, cleanup);

  hattach(fh, cc);

  cc->data_groups = h_subcalloc(cc, 
                                fh->num_datagroups, 
                                sizeof(AFFY_CALVIN_DATAGROUP *));
  if (cc->data_groups == NULL)
    AFFY_HANDLE_ERROR_GOTO("calloc failed", AFFY_ERROR_OUTOFMEM, err, cleanup);

  hattach(cc->data_groups, cc);

  /* Dataheader. */
  cc->data_header = dh = affy_calvin_get_dataheader(cio, err);
  AFFY_CHECK_ERROR_GOTO(err, cleanup);

  hattach(dh, cc);

  /* Data groups. */
  for (i = 0; i < fh->num_datagroups; i++)
  {
    affy_uint32            j;
    AFFY_CALVIN_DATAGROUP *cur_dg;

    cc->data_groups[i] = affy_calvin_get_datagroup_metadata(cio, i, err);
    AFFY_CHECK_ERROR_GOTO(err, cleanup);
    cur_dg = cc->data_groups[i];

    hattach(cur_dg, cc->data_groups);

    /* Allocate the dataset array. */
    cur_dg->datasets = h_subcalloc(cur_dg,
                                   cur_dg->num_datasets,
                                   sizeof(AFFY_CALVIN_DATASET *));
    if (cur_dg->datasets == NULL)
      AFFY_HANDLE_ERROR_GOTO("calloc failed", 
                             AFFY_ERROR_OUTOFMEM, 
                             err, 
                             cleanup);

    /* Process each dataset in this group. */
    for (j = 0; j < cur_dg->num_datasets; j++)
    {
      cur_dg->datasets[j] = affy_calvin_get_dataset_metadata(cio, i, j, err);
      AFFY_CHECK_ERROR_GOTO(err, cleanup);

      hattach(cur_dg->datasets[j], cur_dg->datasets);

      /* data reading functions TBD */
    }
  }

  return (cc);

cleanup:
  affy_free_calvin_container(cc);

  return (NULL);
}
