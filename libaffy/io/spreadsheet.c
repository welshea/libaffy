/**************************************************************************
 *
 * Filename:  spreadsheet.c
 *
 * Purpose:   Operations on a signal/sample spreadsheet.
 *
 * Creation:  4/11/2011
 *
 * Author:    Eric A. Welsh
 *
 * Copyright: Copyright (C) 2011, Moffitt Cancer Center.
 *            All rights reserved.
 *
 * Update History
 * --------------
 * 04/11/11: Creation (EAW)
 *
 **************************************************************************/

#include <affy.h>

/* returns number of rows and columns in spreadsheet */
void get_generic_spreadsheet_bounds(char *filename,
                                    affy_uint32 *return_max_rows,
                                    affy_uint32 *return_max_cols,
                                    AFFY_ERROR *err)
{
  FILE *data_file;
  int max_string_len = 0;
  char *string = NULL;
  char **fields = NULL;
  char *sptr;
  int num_fields = 0;
  int max_field = 0;

  int ok_flag;
  int max_cols = 0;
  int max_rows = 0;
  int i;
  
  data_file = fopen(filename, "rb");
  if (!data_file)
    AFFY_HANDLE_ERROR_VOID("can not open data file", AFFY_ERROR_NOTFOUND, err);
  
  /* read header line */
  fgets_strip_realloc(&string, &max_string_len, data_file);
  num_fields = split_tabs(string, &fields, &max_field);
  
  /* first field is probe, the rest are samples, skip empty columns */
  for (i = 1; i < num_fields; i++)
  {
    ok_flag = 0;

    for (sptr = fields[i]; sptr; sptr++)
    {
      if (!isspace(*sptr))
      {
        ok_flag = 1;
        break;
      }
    }
    
    if (ok_flag)
      max_cols++;
  }
  

  /* assume only a single line per probe */
  /* multiple identical probes will be treated as separate probes */
  /* skip blank probes */
  while(fgets_strip_realloc(&string, &max_string_len, data_file))
  {
    num_fields = split_tabs(string, &fields, &max_field);

    ok_flag = 0;

    for (sptr = fields[0]; sptr; sptr++)
    {
      if (!isspace(*sptr))
      {
        ok_flag = 1;
        break;
      }
    }
    
    if (ok_flag)
      max_rows++;
  }
  
  fclose(data_file);

  *return_max_cols = max_cols;
  *return_max_rows = max_rows;
  
  if (string)
    free(string);
  if (fields)
    free(fields);
}
