
/**************************************************************************
 *
 * Filename:  load_exclusions_file.c
 *
 * Purpose:   Load list of probes/probesets to exclude from IRON training set
 *
 * Creation: 
 *
 * Author:    Eric A. Welsh
 *
 * Copyright: Copyright (C) 2018, Moffitt Cancer Center.
 *            All rights reserved.
 *
 * Update History
 * --------------
 * 2018-06-01: file creation (EAW)
 * 2018-08-29: fixed memory pool reallocation errors (EAW)
 * 2018-09-11: more memory pool leaks (EAW)
 * 2018-09-14: added related new spikein code here as well (EAW)
 * 2019-03-14: changed int mempool to void mempool (EAW)
 * 2020-03-20: handle empty files without crashing (EAW)
 *
 **************************************************************************/

#include <affy.h>
#include "halloc.h"

/* store the list in the CDF structure
 * seemed like a reasonable place to put it
 *
 * CDF structure should already be initialized prior to calling this
 */
void affy_load_exclusions_file(char *filename, AFFY_CDFFILE *cdf,
                               void *mempool,
                               AFFY_ERROR *err)
{
  FILE *data_file;
  int max_string_len = 0;
  char *string = NULL;
  char **fields = NULL;
  char *sptr;
  int num_fields = 0;
  int max_field = 0;
  
  int max_count = 0;
  int count = 0;

  if (cdf->exclusions)
    h_free(cdf->exclusions);
  
  cdf->exclusions = NULL;
  
  data_file = fopen(filename, "rb");
  if (!data_file)
    AFFY_HANDLE_ERROR_VOID("can not open data file", AFFY_ERROR_NOTFOUND, err);

  /* assume the file is just a single column list, no header line */
  while(fgets_strip_realloc(&string, &max_string_len, data_file))
  {
      num_fields = split_tabs(string, &fields, &max_field);
      
      if (num_fields && fields[0][0])
      {
          /* allocate first exclusion */
          if (cdf->exclusions == NULL)
          {
              cdf->exclusions = h_suballoc(cdf, sizeof(char *));
              max_count = 1;
          }
          
          /* increase the size of the pointer array if we need to */
          if (count + 1 > max_count)
          {
              /* allocate a little extra, to avoid some memcpy */
              max_count = 1.01 * (count + 1);
              cdf->exclusions = h_realloc(cdf->exclusions,
                                          max_count * sizeof(char *));
          }
          
          /* store the string */
          cdf->exclusions[count] = h_strdup(fields[0]);
          
          if (mempool)
              hattach(cdf->exclusions[count], mempool);

          count++;
      }
  }
  
  fclose(data_file);
  
  cdf->numexclusions = count;
  
  /* shrink any over-allocated memory */
  if (max_count > count)
      cdf->exclusions = h_realloc(cdf->exclusions, count * sizeof(char *));

  /* sort exclusions strings */
  qsort(cdf->exclusions, count, sizeof(char *), compare_string);

  if (cdf->exclusions)
    hattach(cdf->exclusions, cdf);
  
  if (string)
    free(string);

  if (fields)
    free(fields);
}


/* store the list in the CDF structure
 * seemed like a reasonable place to put it
 *
 * CDF structure should already be initialized prior to calling this
 */
void affy_load_spikeins_file(char *filename, AFFY_CDFFILE *cdf,
                             void *mempool,
                             AFFY_ERROR *err)
{
  FILE *data_file;
  int max_string_len = 0;
  char *string = NULL;
  char **fields = NULL;
  char *sptr;
  int num_fields = 0;
  int max_field = 0;
  
  int max_count = 0;
  int count = 0;

  if (cdf->spikeins)
    h_free(cdf->spikeins);
  
  cdf->spikeins = NULL;
  
  data_file = fopen(filename, "rb");
  if (!data_file)
    AFFY_HANDLE_ERROR_VOID("can not open data file", AFFY_ERROR_NOTFOUND, err);

  /* assume the file is just a single column list, no header line */
  while(fgets_strip_realloc(&string, &max_string_len, data_file))
  {
      num_fields = split_tabs(string, &fields, &max_field);
      
      if (num_fields && fields[0][0])
      {
          /* allocate first exclusion */
          if (cdf->spikeins == NULL)
          {
              cdf->spikeins = h_suballoc(cdf, sizeof(char *));
              max_count = 1;
          }
          
          /* increase the size of the pointer array if we need to */
          if (count + 1 > max_count)
          {
              /* allocate a little extra, to avoid some memcpy */
              max_count = 1.01 * (count + 1);
              cdf->spikeins = h_realloc(cdf->spikeins,
                                        max_count * sizeof(char *));
          }
          
          /* store the string */
          cdf->spikeins[count] = h_strdup(fields[0]);
          
          if (mempool)
              hattach(cdf->spikeins[count], mempool);

          count++;
      }
  }
  
  fclose(data_file);
  
  cdf->numspikeins = count;
  
  /* shrink any over-allocated memory */
  if (max_count > count)
      cdf->spikeins = h_realloc(cdf->spikeins, count * sizeof(char *));

  /* sort spikeins strings */
  qsort(cdf->spikeins, count, sizeof(char *), compare_string);

  if (cdf->spikeins)
    hattach(cdf->spikeins, cdf);
  
  if (string)
    free(string);

  if (fields)
    free(fields);
}
