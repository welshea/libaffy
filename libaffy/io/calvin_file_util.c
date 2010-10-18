
/**************************************************************************
 *
 * Filename:  calvin_file_util.c
 * 
 * Purpose:   Calvin (command console) file utility routines.
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
 * 04/18/08: Add free routine (AMH)
 * 09/20/10: Pooled memory allocator (AMH)
 *
 **************************************************************************/

#include <affy.h>

/* 
 * Be sure to keep the ordering of this table in sync with the
 * AFFY_CALVIN_DATA_TYPE enum.
 */
static const char *calvin_type_labels[] =
  {
    "BYTE",    /* AFFY_CALVIN_BYTE     */
    "UBYTE",   /* AFFY_CALVIN_UBYTE    */
    "SHORT",   /* AFFY_CALVIN_SHORT    */
    "USHORT",  /* AFFY_CALVIN_USHORT   */
    "INT",     /* AFFY_CALVIN_INT      */
    "UINT",    /* AFFY_CALVIN_UINT     */
    "FLOAT",   /* AFFY_CALVIN_FLOAT    */
    "DOUBLE",  /* AFFY_CALVIN_DOUBLE   */
    "STRING",  /* AFFY_CALVIN_STRING   */
    "WSTRING", /* AFFY_CALVIN_WSTRING  */
    "UNKNOWN"  /* AFFY_CALVIN_UNKNOWN  */
  };

static void dump_dataheader(AFFY_CALVIN_DATAHEADER *dh, int depth)
{
  affy_uint32 i;
  char       *padding;

  assert(dh != NULL);
  assert(depth >= 1);
  
  padding = (char *)h_malloc(depth + 1);
  if (padding == NULL)
    return;

  memset((void *)padding, '\t', depth);
  padding[depth] = '\0';

  printf("%sData Header\n", padding);
  printf("%s-----------\n", padding);

  printf("%sType identifier: %s\n", padding, dh->type_identifier);
  printf("%sFile identifier: %s\n", padding, dh->file_identifier);
  printf("%sTimestamp:       %s\n", padding, dh->timestamp);
  printf("%sLocale:          %s\n", padding, dh->locale);

  printf("\n%s\tParameters\n", padding);
  printf("%s\t----------\n", padding);
  
  for (i = 0; i < dh->num_params; i++)
  {  
    printf("%s\t", padding);
    affy_print_calvin_param(dh->params + i);
    printf("\n");
  }

  printf("\n%s\tParent Data Headers\n", padding);
  printf("%s\t-------------------\n\n", padding);

  for (i = 0; i < dh->num_parent_headers; i++)
    dump_dataheader(dh->parent_headers[i], depth + 2);

  h_free(padding);
}

static void dump_dataset(AFFY_CALVIN_DATASET *ds)
{
  affy_uint32  i;
  
  assert(ds != NULL);

  printf("\t\t-- Data Set \"%s\": %" AFFY_PRNu32 " rows, %" AFFY_PRNu32 " cols\n",
         ds->name,
         ds->num_rows,
         ds->num_cols);

  printf("\n\t\t\tParameters\n");
  printf("\t\t\t----------\n");
  
  for (i = 0; i < ds->num_params; i++)
  {
    printf("\t\t\t");
    affy_print_calvin_param(ds->params + i);
    printf("\n");
  }

  printf("\n\n");
}
         
static void dump_datagroup(AFFY_CALVIN_DATAGROUP *dg)
{
  affy_uint32 i;

  assert(dg != NULL);

  printf("\n\t++ Data Group \"%s\", %" AFFY_PRNu32 " dataset(s)\n\n",
         dg->name,
         dg->num_datasets);

  for (i = 0; i < dg->num_datasets; i++)
    dump_dataset(dg->datasets[i]);
}

void affy_free_calvin_dataheader(AFFY_CALVIN_DATAHEADER *dh)
{
  h_free(dh);
}

void affy_free_calvin_column(AFFY_CALVIN_COLUMN *col)
{
  h_free(col);
}

void affy_free_calvin_dataset(AFFY_CALVIN_DATASET *ds)
{
  h_free(ds);
}

void affy_free_calvin_fileheader(AFFY_CALVIN_FILEHEADER *fh)
{
  h_free(fh);
}

void affy_free_calvin_datagroup(AFFY_CALVIN_DATAGROUP *dg)
{
  h_free(dg);
}

void affy_free_calvin_container(AFFY_CALVIN_CONTAINER *cc)
{
  h_free(cc);
}

void affy_calvin_close_dataset(AFFY_CALVIN_DATASET_IO *dio)
{
  h_free(dio);
}

void affy_print_calvin_value(AFFY_CALVIN_DATA data, 
                             AFFY_CALVIN_DATA_TYPE type)
{
  switch (type)
  {
    case AFFY_CALVIN_BYTE:
      printf("%" AFFY_PRNd8, data.byte_val);
      break;

    case AFFY_CALVIN_UBYTE:
      printf("%" AFFY_PRNu8, data.ubyte_val);
      break;

    case AFFY_CALVIN_SHORT:
      printf("%" AFFY_PRNd16, data.short_val);
      break;

    case AFFY_CALVIN_USHORT:
      printf("%" AFFY_PRNu16, data.ushort_val);
      break;

    case AFFY_CALVIN_INT:
      printf("%" AFFY_PRNd32, data.int_val);
      break;

    case AFFY_CALVIN_UINT:
      printf("%" AFFY_PRNu32, data.uint_val);
      break;

    case AFFY_CALVIN_FLOAT:
      printf("%f", data.float_val);
      break;

    case AFFY_CALVIN_DOUBLE:
      printf("%e", data.double_val);
      break;

    case AFFY_CALVIN_WSTRING:   /* WSTRING's were converted to (char *) */
    case AFFY_CALVIN_STRING:
      printf("\"%s\"", data.string_val);
      break;

    default:
      printf("(unknown data type %d)", type);
  }
}

void affy_print_calvin_param(AFFY_CALVIN_PARAM *cp)
{
  assert(cp != NULL);

  printf("%s (%s) = ", cp->name, calvin_type_labels[cp->type]);
  affy_print_calvin_value(cp->value, cp->type);
}

void affy_dump_calvin_container(AFFY_CALVIN_CONTAINER *cc)
{
  affy_int32 i;

  assert(cc != NULL);

  printf("Calvin container version %" AFFY_PRNu8, 
	 cc->file_header->file_version);
  printf(", %" AFFY_PRNu32 " data group(s)\n", 
         cc->file_header->num_datagroups);
  printf("-------------------------------------------------\n\n");

  dump_dataheader(cc->data_header, 1);

  printf("Data Groups\n");
  printf("-----------\n");

  for (i = 0; i < cc->file_header->num_datagroups; i++)
    dump_datagroup(cc->data_groups[i]);

  printf("-------------\nEnd container\n-------------\n");
}

