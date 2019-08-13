
/**************************************************************************
 *
 * Filename:  calvin_io.c
 *
 * Purpose:   Provide low-level primitives for reading Affymetrix Calvin
 *            ("generic") files
 *
 * Creation:  04/18/08
 *
 * Author:    Andrew M. Hoerter
 *
 * Copyright: Copyright (C) 2008, Moffitt Cancer Center.
 *            All rights reserved.
 *
 * Update History
 * --------------
 * 04/18/08: Creation (AMH)
 * 11/26/08: Rework the I/O layer significantly (AMH)
 * 09/20/10: Pooled memory allocator (AMH)
 * 05/13/13: Changes to support salvaging of corrupt CEL files (EAW)
 *
 **************************************************************************/

#include <affy.h>

static char *read_string(FILE *fp, AFFY_ERROR *err);
static char *read_wstring(FILE *fp, AFFY_ERROR *err);
static void  process_dataheader(AFFY_CALVINIO *cio, 
                                AFFY_CALVIN_DATAHEADER *dh,
                                AFFY_ERROR *err);
static void  read_parameter(AFFY_CALVINIO *cio, 
                            AFFY_CALVIN_PARAM *param,
                            AFFY_ERROR *err);
static void  move_to_datagroup(AFFY_CALVINIO *cio, 
                               affy_uint32 dg_idx, 
                               AFFY_ERROR *err);
static int   move_to_dataset(AFFY_CALVINIO *cio, 
                             affy_uint32 dg_idx, 
                             affy_uint32 ds_idx,
                             AFFY_ERROR *err);


typedef int(*affy_binary_reader)(FILE *, void *);

struct calvin_type_info
{ 
  const char        *label; 
  affy_binary_reader read_func; 
  size_t             size;
};

/* 
 * Be sure to keep the ordering of this table in sync with the
 * AFFY_CALVIN_DATA_TYPE enum.
 */
static const struct calvin_type_info calvin_type_table[] =
  {
    /* AFFY_CALVIN_BYTE    */
    { "text/x-calvin-integer-8",           affy_read8,     1 },
    /* AFFY_CALVIN_UBYTE   */
    { "text/x-calvin-unsigned-integer-8",  affy_read8,     1 },
    /* AFFY_CALVIN_SHORT   */
    { "text/x-calvin-integer-16",          affy_read16_be, 2 },
    /* AFFY_CALVIN_USHORT  */
    { "text/x-calvin-unsigned-integer-16", affy_read16_be, 2 },
    /* AFFY_CALVIN_INT     */
    { "text/x-calvin-integer-32",          affy_read32_be, 4 },
    /* AFFY_CALVIN_UINT    */
    { "text/x-calvin-unsigned-integer-32", affy_read32_be, 4 },
    /* AFFY_CALVIN_FLOAT   */
    { "text/x-calvin-float",               affy_read32_be, 4 },
    /* AFFY_CALVIN_DOUBLE  */
    { NULL,                                affy_read64_be, 8 },
    /* AFFY_CALVIN_STRING  */
    { "text/ascii",                        NULL,           0 },
    /* AFFY_CALVIN_WSTRING */
    { "text/plain",                        NULL,           0 }
  };

/* Translate parameter type strings into an internal type enum. */
static AFFY_CALVIN_DATA_TYPE paramtype_from_string(char *s)
{
  size_t tab_size;
  int    i;

  tab_size = sizeof(calvin_type_table) / sizeof(struct calvin_type_info);

  for (i = 0; i < tab_size; i++)
  {
    if ((calvin_type_table[i].label != NULL) &&
        (strcmp(calvin_type_table[i].label, s) == 0))
      return (i);
  }

  /* Else, it's unknown. */
  return (AFFY_CALVIN_UNKNOWN);
}

static const AFFY_CALVIN_COLUMN_MAPPING *mapping_for_column(const AFFY_CALVIN_COLUMN_MAPPING *mappings,
                                                      const char *name)
{
  const AFFY_CALVIN_COLUMN_MAPPING *cmap;

  assert(name     != NULL);
  assert(mappings != NULL);
  
  for (cmap = mappings; cmap->name != NULL; cmap++)
  {
    if (strcasecmp(name, cmap->name) == 0)
      return (cmap);
  }

  return (NULL);
}

/*
 * A reader for simple, ASCII/C strings, with the slight added twist
 * of a length field at the beginning.
 */
static char *read_string(FILE *fp, AFFY_ERROR *err)
{
  affy_int32 sz;
  char      *result;
  int        i;

  assert(fp != NULL);

  if (affy_read32_be(fp, (void *)(&sz)) != 0)
    AFFY_HANDLE_ERROR("I/O error reading Calvin file", 
		      AFFY_ERROR_IO, 
		      err, 
		      NULL);

  if (sz < 0)
    AFFY_HANDLE_ERROR("corrupt string size in Calvin file", 
		      AFFY_ERROR_BADFORMAT,
		      err,
		      NULL);

  result = h_malloc(sz + 1);
  if (result == NULL)
    AFFY_HANDLE_ERROR("malloc failed", AFFY_ERROR_OUTOFMEM, err, NULL);

  for (i = 0; i < sz; i++)
  {
    if (affy_read8(fp, (void *)(result + i)) != 0)
    {
      h_free(result);

      AFFY_HANDLE_ERROR("I/O error reading Calvin file", 
			AFFY_ERROR_IO, 
			err, 
			NULL);
    }
  }

  result[sz] = '\0';

  return (result);
}

/*
 * A sort of pseudo-reader for Unicode type strings.  Basically we
 * fake it and discard every other byte.  While by far not a
 * desireable solution, this will handle probably 99% of the data
 * files likely to be encountered until a more correct mechanism is
 * implemented.
 */
static char *read_wstring(FILE *fp, AFFY_ERROR *err)
{
  affy_int32  len;
  char       *result_buf;
  int         i;

  assert(fp != NULL);

  /* Get the string length first. */
  if (affy_read32_be(fp, (void *)(&len)) != 0)
    AFFY_HANDLE_ERROR("I/O error reading Calvin file", 
		      AFFY_ERROR_IO, 
		      err, 
		      NULL);

  result_buf = h_malloc(len + 1);
  if (result_buf == NULL)
    AFFY_HANDLE_ERROR("malloc failed", AFFY_ERROR_OUTOFMEM, err, NULL);

  for (i = 0; i < len; i++)
  {
    affy_uint8 tmp[2];

    if (affy_read16(fp, (void *)(&tmp)) != 0)
    {
#if 0
      h_free(result_buf);
      AFFY_HANDLE_ERROR("I/O error reading Calvin file",
			AFFY_ERROR_IO,
			err,
			NULL);
#else
      /* rather than exit with an error, return NULL
       * in order to support salvaging of corrupt CEL files
       */
      result_buf[0] = '\0';
      return result_buf;
#endif
    }

    result_buf[i] = tmp[1];
  }

  result_buf[len] = '\0';

  return (result_buf);
}

/*
 * Function to do the actual work of reading the dataheader; assumes
 * the filepointer is positioned correctly to read the next header.
 */
static void process_dataheader(AFFY_CALVINIO *cio, 
                               AFFY_CALVIN_DATAHEADER *dh,
                               AFFY_ERROR *err)
{
  affy_int32 i;

  assert(cio != NULL);
  assert(dh  != NULL);

  dh->type_identifier    = NULL;
  dh->file_identifier    = NULL;
  dh->timestamp          = NULL;
  dh->locale             = NULL;
  dh->num_params         = 0;
  dh->num_parent_headers = 0;
  dh->parent_headers     = NULL;
  dh->params             = NULL;

  /* Type identifier (GUID) */
  dh->type_identifier = read_string(cio->fp, err);
  AFFY_CHECK_ERROR_VOID(err);
  hattach(dh->type_identifier, dh);

  /* File identifier (GUID) */
  dh->file_identifier = read_string(cio->fp, err);
  AFFY_CHECK_ERROR_VOID(err);
  hattach(dh->file_identifier, dh);

  /* Timestamp (WSTRING) */
  dh->timestamp = read_wstring(cio->fp, err);
  AFFY_CHECK_ERROR_VOID(err);
  hattach(dh->timestamp, dh);

  /* Locale (WSTRING?) */
  dh->locale = read_wstring(cio->fp, err);
  AFFY_CHECK_ERROR_VOID(err);
  hattach(dh->locale, dh);

  /* Number of params (INT) */
  if (affy_read32_be(cio->fp, (void *)(&dh->num_params)) != 0)
    AFFY_HANDLE_ERROR_VOID("I/O error reading Calvin file",
			   AFFY_ERROR_IO,
			   err);

  if (dh->num_params < 0)
    AFFY_HANDLE_ERROR_VOID("corrupt field in dataheader",
                           AFFY_ERROR_BADFORMAT,
                           err);

  if (dh->num_params > 0)
  {
    dh->params = h_subcalloc(dh, dh->num_params, sizeof(AFFY_CALVIN_PARAM));
    if (dh->params == NULL)
      AFFY_HANDLE_ERROR_VOID("calloc failed", 
			     AFFY_ERROR_OUTOFMEM,
			     err);
    
    /* Parameter array (WSTRING, STRING, WSTRING) */
    for (i = 0; i < dh->num_params; i++)
    {
      AFFY_CALVIN_PARAM *param = dh->params + i;

      read_parameter(cio, param, err);
      AFFY_CHECK_ERROR_VOID(err);

      hattach(param->name, dh->params);
      if (param->type == AFFY_CALVIN_STRING)
        hattach(param->value.string_val, dh->params);
    }
  }
  
  /* Process parent data headers. */
  if (affy_read32_be(cio->fp, (void *)(&dh->num_parent_headers)) != 0)
    AFFY_HANDLE_ERROR_VOID("I/O error reading Calvin file",
			   AFFY_ERROR_IO,
			   err);

  if (dh->num_parent_headers < 0)
    AFFY_HANDLE_ERROR_VOID("corrupt field in dataheader",
                           AFFY_ERROR_BADFORMAT,
                           err);

  if (dh->num_parent_headers > 0)
  {
    dh->parent_headers = h_subcalloc(dh, 
                                     dh->num_parent_headers, 
                                     sizeof(AFFY_CALVIN_DATAHEADER *));
    if (dh->parent_headers == NULL)
      AFFY_HANDLE_ERROR_VOID("calloc failed",
                             AFFY_ERROR_OUTOFMEM,
                             err);
    
    hattach(dh->parent_headers, dh);

    for (i = 0; i < dh->num_parent_headers; i++)
    {
      dh->parent_headers[i] = h_suballoc(dh->parent_headers,
                                         sizeof(AFFY_CALVIN_DATAHEADER));
      if (dh->parent_headers[i] == NULL)
        AFFY_HANDLE_ERROR_VOID("malloc failed", AFFY_ERROR_OUTOFMEM, err);
      
      process_dataheader(cio, 
                         dh->parent_headers[i], 
                         err);
      AFFY_CHECK_ERROR_VOID(err);
    }
  }
}

/*
 * Parse a Calvin parameter.
 */
static void read_parameter(AFFY_CALVINIO *cio, 
                           AFFY_CALVIN_PARAM *param, 
                           AFFY_ERROR *err)
{
  affy_int32            len;
  char                 *type;
  long                  value_pos, end_pos;
  AFFY_CALVIN_DATA_TYPE input_type;

  assert(cio   != NULL);
  assert(param != NULL);

  param->type = AFFY_CALVIN_UNKNOWN;

  /* Parameter name (WSTRING) */
  param->name = read_wstring(cio->fp, err);
  AFFY_CHECK_ERROR_VOID(err);

  /* First, fetch the length of the value string, saving our position. */
  if ((value_pos = ftell(cio->fp)) == -1)
    AFFY_HANDLE_ERROR_GOTO("I/O error reading Calvin file", 
			   AFFY_ERROR_IO,
			   err,
                           cleanup);
  
  if (affy_read32_be(cio->fp, (void *)(&len)) != 0)
    AFFY_HANDLE_ERROR_GOTO("I/O error reading Calvin file", 
			   AFFY_ERROR_IO,
			   err,
                           cleanup);

  if (len < 0)
    AFFY_HANDLE_ERROR_GOTO("corrupt parameter in Calvin file", 
                           AFFY_ERROR_BADFORMAT,
                           err,
                           cleanup);

  /* Seek ahead to the type string. */
  if (fseek(cio->fp, len, SEEK_CUR) != 0)
    AFFY_HANDLE_ERROR_GOTO("I/O error reading Calvin file", 
			   AFFY_ERROR_IO,
			   err,
                           cleanup);

  type = read_wstring(cio->fp, err);
  AFFY_CHECK_ERROR_GOTO(err, cleanup);

  input_type = paramtype_from_string(type);
  h_free(type);

  /* Save end-of-parameter position. */
  if ((end_pos = ftell(cio->fp)) == -1)
    AFFY_HANDLE_ERROR_GOTO("I/O error reading Calvin file", 
			   AFFY_ERROR_IO,
			   err,
                           cleanup);

  /* Now rewind and read the actual value. */
  if (fseek(cio->fp, value_pos, SEEK_SET) != 0)
    AFFY_HANDLE_ERROR_GOTO("I/O error reading Calvin file", 
			   AFFY_ERROR_IO,
			   err,
                           cleanup);

  /* XXX But... if it isn't a string, skip the length INT. */
  if ((input_type != AFFY_CALVIN_STRING) 
      && (input_type != AFFY_CALVIN_WSTRING))
    if (fseek(cio->fp, 4, SEEK_CUR) != 0)
      AFFY_HANDLE_ERROR_GOTO("I/O error reading Calvin file", 
			     AFFY_ERROR_IO,
			     err,
                             cleanup);

  affy_calvin_read_data(cio, &param->value, input_type, err);
  AFFY_CHECK_ERROR_GOTO(err, cleanup);

  /* If the input type was a WSTRING, we converted it to a STRING. */
  if (input_type == AFFY_CALVIN_WSTRING)
    param->type = AFFY_CALVIN_STRING;
  else
    param->type = input_type;
  
  if (fseek(cio->fp, end_pos, SEEK_SET) != 0)
    AFFY_HANDLE_ERROR_GOTO("I/O error reading Calvin file",
                           AFFY_ERROR_IO,
                           err,
                           cleanup);

  return;

cleanup:
  h_free(param->name);
  if (param->type == AFFY_CALVIN_STRING)
    h_free(param->value.string_val);
}

/* 
 * Seek to the Nth datagroup in the container, leaving the file
 * pointer resting at the first field.
 */  
static void move_to_datagroup(AFFY_CALVINIO *cio, 
                              affy_uint32 dg_idx, 
                              AFFY_ERROR *err) 
{
  affy_uint32 next_ofs;
  affy_uint32 cur_dg = 0;

  assert(cio != NULL);

  if (dg_idx >= cio->num_datagroups)
    AFFY_HANDLE_ERROR_VOID("index out of range",
                           AFFY_ERROR_BADPARAM,
                           err);

  next_ofs = cio->first_datagroup;

  do
  {
    affy_uint32 ofs_save = next_ofs;

    if (fseek(cio->fp, next_ofs, SEEK_SET) != 0)
      AFFY_HANDLE_ERROR_VOID("I/O error reading Calvin file", 
                             AFFY_ERROR_IO,
                             err);

    /* Get next offset and return to beginning. */
    if (affy_read32_be(cio->fp, (void *)(&next_ofs)) != 0)
      AFFY_HANDLE_ERROR_VOID("I/O error reading Calvin file", 
                             AFFY_ERROR_IO,
                             err);

    if (fseek(cio->fp, ofs_save, SEEK_SET) != 0)
      AFFY_HANDLE_ERROR_VOID("I/O error reading Calvin file", 
                             AFFY_ERROR_IO,
                             err);
  } while (cur_dg++ != dg_idx);
}

/* 
 * Seek to the Mth dataset within the Nth datagroup in the container,
 * leaving the file pointer resting at the first field.
 */  
static int move_to_dataset(AFFY_CALVINIO *cio, 
                           affy_uint32 dg_idx, 
                           affy_uint32 ds_idx,
                           AFFY_ERROR *err)
{
  affy_uint32 num_datasets, first_dataset, next_ofs, cur_ds = 0;
  affy_int32  signed_num_datasets;

  assert(cio != NULL);

  move_to_datagroup(cio, dg_idx, err);
  AFFY_CHECK_ERROR_VOID_ZERO(err);

  if (fseek(cio->fp, 4, SEEK_CUR) != 0)
    AFFY_HANDLE_ERROR_VOID_ZERO("I/O error reading Calvin file", 
                                AFFY_ERROR_IO,
                                err);

  if (affy_read32_be(cio->fp, (void *)(&first_dataset)) != 0)
    AFFY_HANDLE_ERROR_VOID_ZERO("I/O error reading Calvin file", 
                                AFFY_ERROR_IO,
                                err);

  if (affy_read32_be(cio->fp, (void *)(&signed_num_datasets)) != 0)
    AFFY_HANDLE_ERROR_VOID_ZERO("I/O error reading Calvin file", 
                                AFFY_ERROR_IO, 
                                err);

  if (signed_num_datasets <= 0)
    AFFY_HANDLE_ERROR_VOID_ZERO("empty data group", AFFY_ERROR_BADFORMAT, err);

  num_datasets = signed_num_datasets; /* guaranteed positive */

#if 0
  if (ds_idx >= num_datasets)
    AFFY_HANDLE_ERROR_VOID("index out of range", AFFY_ERROR_BADPARAM, err);
#else
  /* HACK -- file is likely truncated, try to salvage somewhat gracefully */
  if (ds_idx >= num_datasets)
    return -1;
#endif
  
  next_ofs = first_dataset;

  do
  {
    affy_uint32 ofs_save = next_ofs;

    if (fseek(cio->fp, next_ofs + 4, SEEK_SET) != 0)
      AFFY_HANDLE_ERROR_VOID_ZERO("I/O error reading Calvin file", 
                                  AFFY_ERROR_IO,
                                  err);

    /* Get next offset and return to beginning. */
    if (affy_read32_be(cio->fp, (void *)(&next_ofs)) != 0)
      AFFY_HANDLE_ERROR_VOID_ZERO("I/O error reading Calvin file", 
                                  AFFY_ERROR_IO,
                                  err);

    if (fseek(cio->fp, ofs_save, SEEK_SET) != 0)
      AFFY_HANDLE_ERROR_VOID_ZERO("I/O error reading Calvin file", 
                                  AFFY_ERROR_IO,
                                  err);
  } while (cur_ds++ != ds_idx);
  
  return 0;
}

/* 
 * Free a Calvin I/O context.
 */
void affy_calvinio_free(AFFY_CALVINIO *cio)
{
  h_free(cio);
}

/*
 * Reset a Calvin I/O context (currently just a stub).
 */
void affy_calvinio_reset(AFFY_CALVINIO *cio)
{
  assert(cio != NULL);
}

/* Allocate and initialize a Calvin I/O context. */
AFFY_CALVINIO *affy_calvinio_init(FILE *fp, AFFY_ERROR *err)
{
  AFFY_CALVINIO *result;
  affy_uint8     magic, version;
  affy_int32     num_dgs;
  affy_uint32    first_dg;

  assert(fp != NULL);

  /* First, verify the magic number. */
  if (affy_read8(fp, (void *)(&magic)) != 0)
    AFFY_HANDLE_ERROR("I/O error reading Calvin file", 
                      AFFY_ERROR_IO,
                      err,
                      NULL);

  if (magic != AFFY_CALVIN_FILEMAGIC)
    AFFY_HANDLE_ERROR("bad Calvin file magic", 
                      AFFY_ERROR_BADFORMAT,
                      err,
                      NULL);

  if (affy_read8(fp, (void *)(&version)) != 0)
    AFFY_HANDLE_ERROR("I/O error reading Calvin file", 
                      AFFY_ERROR_IO,
                      err,
                      NULL);

  if (affy_read32_be(fp, (void *)(&num_dgs)) != 0)
    AFFY_HANDLE_ERROR("I/O error reading Calvin file", 
                      AFFY_ERROR_IO,
                      err,
                      NULL);

  if (num_dgs < 0)
    AFFY_HANDLE_ERROR("Corrupt Calvin header",
                      AFFY_ERROR_BADFORMAT,
                      err,
                      NULL);
  
  if (affy_read32_be(fp, (void *)(&first_dg)) != 0)
    AFFY_HANDLE_ERROR("I/O error reading Calvin file", 
                      AFFY_ERROR_IO,
                      err,
                      NULL);

  /* Create and initialize the structure. */
  result = h_malloc(sizeof(AFFY_CALVINIO));
  if (result == NULL)
    AFFY_HANDLE_ERROR("malloc failed", AFFY_ERROR_OUTOFMEM, err, NULL);

  result->fp              = fp;
  result->first_datagroup = first_dg;
  result->num_datagroups  = num_dgs;
  result->file_version    = version;

  return (result);
}

/* 
 * A generic reader for Calvin data entities.
 */
void affy_calvin_read_data(AFFY_CALVINIO *cio,
                           AFFY_CALVIN_DATA *dest,
                           AFFY_CALVIN_DATA_TYPE type,
                           AFFY_ERROR *err)
{
  assert(dest != NULL);
  assert(cio  != NULL);

  if (type == AFFY_CALVIN_STRING)
  {
    dest->string_val = read_string(cio->fp, err);
    AFFY_CHECK_ERROR_VOID(err);
  }
  else if (type == AFFY_CALVIN_WSTRING)
  {
    dest->string_val = read_wstring(cio->fp, err);
    AFFY_CHECK_ERROR_VOID(err);
  }
  else
  {
    void              *dest_field;
    affy_binary_reader read_func;

    switch (type)
    {
      case (AFFY_CALVIN_BYTE):
        dest_field = (void *)(&dest->byte_val);
        break;

      case (AFFY_CALVIN_UBYTE):
        dest_field = (void *)(&dest->ubyte_val);
        break;
        
      case (AFFY_CALVIN_SHORT):
        dest_field = (void *)(&dest->short_val);
        break;
        
      case (AFFY_CALVIN_USHORT):
        dest_field = (void *)(&dest->ushort_val);
        break;
        
      case (AFFY_CALVIN_INT):
        dest_field = (void *)(&dest->int_val);
        break;

      case (AFFY_CALVIN_UINT):
        dest_field = (void *)(&dest->uint_val);
        break;

      case (AFFY_CALVIN_FLOAT):
        dest_field = (void *)(&dest->float_val);
        break;

      case (AFFY_CALVIN_DOUBLE):
        dest_field = (void *)(&dest->double_val);
        break;
    
      default:
        AFFY_HANDLE_ERROR_VOID("unknown calvin data type",
                               AFFY_ERROR_BADPARAM,
                               err);
    }

    /* 
     * Unknown types were already trapped in the case statement, so
     * it's safe to do the read. 
     */
    read_func = calvin_type_table[type].read_func;
    if (read_func(cio->fp, dest_field) == -1)
      AFFY_HANDLE_ERROR_VOID("I/O error reading Calvin file",
                             AFFY_ERROR_IO,
                             err);
  }
}

/* 
 * Fetch metadata for the file itself.
 */
AFFY_CALVIN_FILEHEADER *affy_calvin_get_file_metadata(AFFY_CALVINIO *cio,
                                                      AFFY_ERROR *err)
{
  AFFY_CALVIN_FILEHEADER *fh = NULL;

  assert(cio != NULL);

  fh = h_malloc(sizeof(AFFY_CALVIN_FILEHEADER));
  if (fh == NULL)
    AFFY_HANDLE_ERROR("malloc failed", AFFY_ERROR_OUTOFMEM, err, NULL);

  fh->file_version   = cio->file_version;
  fh->num_datagroups = cio->num_datagroups;

  return (fh);
}

/*
 * Fetch metadata for the given datagroup.
 */
AFFY_CALVIN_DATAGROUP *affy_calvin_get_datagroup_metadata(AFFY_CALVINIO *cio,
                                                          affy_uint32 dg_index,
                                                          AFFY_ERROR *err)
{
  affy_int32             signed_num_datasets;
  AFFY_CALVIN_DATAGROUP *dg = NULL;

  assert(cio != NULL);

  move_to_datagroup(cio, dg_index, err);
  AFFY_CHECK_ERROR(err, NULL);

  if (fseek(cio->fp, 8, SEEK_CUR) != 0)
    AFFY_HANDLE_ERROR("I/O error reading Calvin file", 
                      AFFY_ERROR_IO,
                      err,
                      NULL);

  if (affy_read32_be(cio->fp, (void *)(&signed_num_datasets)) != 0)
    AFFY_HANDLE_ERROR("I/O error reading Calvin file", 
                      AFFY_ERROR_IO, 
                      err,
                      NULL);
  
  if (signed_num_datasets <= 0)
    AFFY_HANDLE_ERROR("empty data group", AFFY_ERROR_BADFORMAT, err, NULL);

  dg = h_malloc(sizeof(AFFY_CALVIN_DATAGROUP));
  if (dg == NULL)
    AFFY_HANDLE_ERROR("malloc failed", AFFY_ERROR_OUTOFMEM, err, NULL);

  dg->datasets     = NULL;
  dg->num_datasets = signed_num_datasets;
  dg->name         = read_wstring(cio->fp, err);
  AFFY_CHECK_ERROR_GOTO(err, cleanup);

  hattach(dg->name, dg);

  return (dg);

cleanup:
  h_free(dg);

  return (NULL);
}

/*
 * Fetch metadata for the given dataset.
 */
AFFY_CALVIN_DATASET *affy_calvin_get_dataset_metadata(AFFY_CALVINIO *cio,
                                                      affy_uint32 dg_index,
                                                      affy_uint32 ds_index,
                                                      AFFY_ERROR *err)
{
  affy_int32           p;
  affy_uint32          col;
  AFFY_CALVIN_DATASET *ds;

  assert(cio != NULL);

  if (move_to_dataset(cio, dg_index, ds_index, err))
    AFFY_CHECK_ERROR(err, NULL);

  ds = h_malloc(sizeof(AFFY_CALVIN_DATASET));
  if (ds == NULL)
    AFFY_HANDLE_ERROR("malloc failed", AFFY_ERROR_OUTOFMEM, err, NULL);

  ds->num_params = 0;
  ds->num_cols   = 0;
  ds->columns    = NULL;
  ds->data       = NULL;
  ds->params     = NULL;
  ds->name       = NULL;

  if (fseek(cio->fp, 8, SEEK_CUR) != 0)
    AFFY_HANDLE_ERROR_GOTO("I/O error reading Calvin file", 
                           AFFY_ERROR_IO,
                           err,
                           cleanup);

  /* Name (WSTRING) */
  ds->name = read_wstring(cio->fp, err);
  AFFY_CHECK_ERROR_GOTO(err, cleanup);

  hattach(ds->name, ds);

  /* Number of params (INT) */
  if (affy_read32_be(cio->fp, (void *)(&ds->num_params)) != 0)
    AFFY_HANDLE_ERROR_GOTO("I/O error reading Calvin file",
                           AFFY_ERROR_IO,
                           err,
                           cleanup);
  
  if (ds->num_params < 0)
    AFFY_HANDLE_ERROR_GOTO("corrupt field in dataset header",
                           AFFY_ERROR_BADFORMAT,
                           err,
                           cleanup);

  if (ds->num_params > 0)
  {
    ds->params = h_subcalloc(ds, ds->num_params, sizeof(AFFY_CALVIN_PARAM));
    if (ds->params == NULL)
      AFFY_HANDLE_ERROR_GOTO("calloc failed",
                             AFFY_ERROR_OUTOFMEM,
                             err,
                             cleanup);
    
    /* Parameter array (WSTRING, STRING, WSTRING) */
    for (p = 0; p < ds->num_params; p++)
    {
      AFFY_CALVIN_PARAM *param = ds->params + p;

      read_parameter(cio, param, err);
      AFFY_CHECK_ERROR_GOTO(err, cleanup);

      hattach(param->name, ds->params);
      if (param->type == AFFY_CALVIN_STRING)
        hattach(param->value.string_val, ds->params);
    }
  }

  /* Number of columns (UINT) */
  if (affy_read32_be(cio->fp, (void *)(&ds->num_cols)) != 0)
    AFFY_HANDLE_ERROR_GOTO("I/O error reading Calvin file",
                           AFFY_ERROR_IO,
                           err,
                           cleanup);

  ds->columns = (AFFY_CALVIN_COLUMN *)h_subcalloc(ds,
                                                  ds->num_cols,
                                                  sizeof(AFFY_CALVIN_COLUMN));
  if (ds->columns == NULL)
    AFFY_HANDLE_ERROR_GOTO("calloc failed",
                           AFFY_ERROR_OUTOFMEM,
                           err,
                           cleanup);

  /* Column metadata */
  for (col = 0; col < ds->num_cols; col++)
  {
    affy_uint8 tmp;

    ds->columns[col].name = read_wstring(cio->fp, err);
    AFFY_CHECK_ERROR_GOTO(err, cleanup);

    hattach(ds->columns[col].name, ds->columns);

    if (affy_read8(cio->fp, (void *)(&tmp)) != 0)
      AFFY_HANDLE_ERROR_GOTO("I/O error reading Calvin file",
			     AFFY_ERROR_IO,
			     err,
			     cleanup);

    if (affy_read32_be(cio->fp, (void *)(&(ds->columns[col].size))) != 0)
      AFFY_HANDLE_ERROR_GOTO("I/O error reading Calvin file",
			     AFFY_ERROR_IO,
			     err,
			     cleanup);

    ds->columns[col].type = tmp;
  }

  /* Number of rows (UINT) */
  if (affy_read32_be(cio->fp, (void *)(&(ds->num_rows))) != 0)
    AFFY_HANDLE_ERROR_GOTO("I/O error reading Calvin file",
                           AFFY_ERROR_IO,
                           err,
                           cleanup);

  return (ds);

cleanup:
  affy_free_calvin_dataset(ds);
  
  return (NULL);
}

/* 
 * Given a datagroup name, locate the index, if it is present.
 */
affy_int32 affy_calvin_find_datagroup_index(AFFY_CALVINIO *cio,
                                            const char *datagroup_name,
                                            AFFY_ERROR *err)
{
  affy_uint32 dg_idx;

  assert(cio            != NULL);
  assert(datagroup_name != NULL);

  for (dg_idx = 0; dg_idx < cio->num_datagroups; dg_idx++)
  {
    char *name;
    int   compare;

    move_to_datagroup(cio, dg_idx, err);
    AFFY_CHECK_ERROR(err, -1);

    /* Seek to datagroup name */
    if (fseek(cio->fp, 12, SEEK_CUR) != 0)
      AFFY_HANDLE_ERROR("I/O error reading Calvin file", 
                        AFFY_ERROR_IO,
                        err,
                        -1);
        
    name = read_wstring(cio->fp, err);
    AFFY_CHECK_ERROR(err, -1);

    compare = strcasecmp(datagroup_name, name);
    h_free(name);

    if (compare == 0)
      return (dg_idx);
  }

  return (-1);
}

/* 
 * Given a dataset name and a datagroup index, locate the dataset
 * index, if it is present.
 */
affy_int32 affy_calvin_find_dataset_index(AFFY_CALVINIO *cio,
                                          affy_uint32 dg_index,
                                          const char *dataset_name,
                                          AFFY_ERROR *err)
{
  affy_uint32 ds_index;

  assert(cio          != NULL);
  assert(dataset_name != NULL);

  for (ds_index = 0;; ds_index++)
  {
    char *name;
    int   compare;

    if (move_to_dataset(cio, dg_index, ds_index, err))
    {
      return (-1);
/*
    AFFY_CHECK_ERROR(err, -1);
*/
    }

    /* Seek to dataset name */
    if (fseek(cio->fp, 8, SEEK_CUR) != 0)
      AFFY_HANDLE_ERROR("I/O error reading Calvin file", 
                        AFFY_ERROR_IO,
                        err,
                        -1);

    name = read_wstring(cio->fp, err);
    AFFY_CHECK_ERROR(err, -1);

    compare = strcasecmp(dataset_name, name);
    h_free(name);

    if (compare == 0)
      return (ds_index);
  }
}

affy_int32 affy_calvin_find_column_index(AFFY_CALVIN_DATASET_IO *dio,
                                         const char *column_name,
                                         AFFY_ERROR *err)
{
  AFFY_CALVIN_DATASET *metadata;
  affy_uint32          i;

  assert(dio         != NULL);
  assert(column_name != NULL);

  metadata = dio->metadata;

  for (i = 0; i < metadata->num_cols; i++)
  {
    if (strcasecmp(column_name, metadata->columns[i].name) == 0)
      return (i);
  }

  return (-1);
}

/* 
 * Find the specified parameter by name in an array of parameters of
 * the given size and return a pointer to it, or NULL if the name was
 * not found.  No storage is allocated here, the reference returned is
 * to an existing allocation.
 */
AFFY_CALVIN_PARAM *affy_calvin_find_param(AFFY_CALVIN_PARAM *params,
                                          affy_uint32 num_params,
                                          const char *name)
{
  affy_uint32 i;

  assert(params != NULL);

  for (i = 0; i < num_params; i++)
  {
    if (strcasecmp(name, params[i].name) == 0)
      return (params + i);
  }
  
  return (NULL);
}

AFFY_CALVIN_DATASET_IO *affy_calvin_prepare_dataset(AFFY_CALVINIO *cio,
                                                    affy_uint32 dg_index,
                                                    affy_uint32 ds_index,
                                                    AFFY_ERROR *err)
{
  affy_uint32             i;
  AFFY_CALVIN_DATASET_IO *dio = NULL;

  assert(cio != NULL);

  dio = h_malloc(sizeof(AFFY_CALVIN_DATASET_IO));
  if (dio == NULL)
    AFFY_HANDLE_ERROR("malloc failed", AFFY_ERROR_OUTOFMEM, err, NULL);

  dio->calvin_io = cio;

  if (move_to_dataset(cio, dg_index, ds_index, err))
    AFFY_CHECK_ERROR_GOTO(err, cleanup);

  /* Save offset to first data element in the array */
  if (affy_read32_be(cio->fp, (void *)(&(dio->initial_offset))) != 0)
    AFFY_HANDLE_ERROR_GOTO("I/O error reading Calvin file",
                           AFFY_ERROR_IO,
                           err,
                           cleanup);

  /* Cache metadata */
  dio->metadata = affy_calvin_get_dataset_metadata(cio,
                                                   dg_index, 
                                                   ds_index, 
                                                   err);
  AFFY_CHECK_ERROR_GOTO(err, cleanup);
  
  hattach(dio->metadata, dio);

  /* Cache the length of a row in bytes. */
  for (i = 0, dio->row_length = 0; i < dio->metadata->num_cols; i++)
    dio->row_length += dio->metadata->columns[i].size;

  return (dio);

cleanup:
  h_free(dio);

  return (NULL);
}


/* 
 * Read a single column from the given dataset I/O context.
 */
void affy_calvin_read_dataset_col(AFFY_CALVIN_DATASET_IO *dio,
                                  LIBUTILS_PB_STATE *pbs,
                                  affy_uint32 col_index,
                                  void *dest,
                                  AFFY_ERROR *err)
{
  AFFY_CALVIN_DATASET  *metadata;
  AFFY_CALVIN_DATA_TYPE type;
  AFFY_CALVINIO        *cio;
  affy_uint32           skip_distance, column_offset, i;
  char                **string_dest = dest, *byte_dest = dest;
  affy_binary_reader    read_func = NULL;
 
  assert(dio  != NULL);
  assert(dest != NULL);

  metadata = dio->metadata;
  cio      = dio->calvin_io;

  if (col_index >= metadata->num_cols)
    AFFY_HANDLE_ERROR_VOID("column index out of range",
                           AFFY_ERROR_BADPARAM,
                           err);

  type = metadata->columns[col_index].type;

  switch (type)
  {
    case AFFY_CALVIN_STRING:
    case AFFY_CALVIN_WSTRING:
      break;
    case AFFY_CALVIN_UNKNOWN:  
      AFFY_HANDLE_ERROR_VOID("unknown column type",
                             AFFY_ERROR_BADFORMAT,
                             err);
    default:
      read_func = calvin_type_table[type].read_func;
  }

  /* Calculate offsets. */
  for (i = column_offset = 0; i < metadata->num_cols; i++)
  {
    if (i < col_index)
      column_offset += metadata->columns[i].size;
  }

  skip_distance = dio->row_length - metadata->columns[col_index].size;

  if (fseek(cio->fp, dio->initial_offset + column_offset, SEEK_SET) != 0)
    AFFY_HANDLE_ERROR_VOID("I/O error reading Calvin file",
                           AFFY_ERROR_IO,
                           err);

  for (i = 0; i < metadata->num_rows; i++)
  {

    if (type == AFFY_CALVIN_WSTRING)
    {
      *(string_dest++) = read_wstring(cio->fp, err);
      AFFY_CHECK_ERROR_VOID(err);
    }
    else if (type == AFFY_CALVIN_STRING)
    {
      *(string_dest++) = read_string(cio->fp, err);
      AFFY_CHECK_ERROR_VOID(err);
    }
    else
    {
      if (read_func(cio->fp, (void *)byte_dest) != 0)
        AFFY_HANDLE_ERROR_VOID("I/O error reading Calvin file",
                               AFFY_ERROR_IO,
                               err);
      byte_dest += calvin_type_table[type].size;
    }
    
    /* Move to the same column on the next row. */
    if (fseek(cio->fp, skip_distance, SEEK_CUR) != 0)
      AFFY_HANDLE_ERROR_VOID("I/O error reading Calvin file",
                             AFFY_ERROR_IO,
                             err);
    pb_tick(pbs, 1,"");
  }
}

void affy_calvin_read_dataset_rows(AFFY_CALVIN_DATASET_IO *dio,
                                   LIBUTILS_PB_STATE *pbs,                  
                                   affy_uint32 start_row,
                                   affy_uint32 num_rows,
                                   void *base,
                                   size_t base_sz,
                                   const AFFY_CALVIN_COLUMN_MAPPING *offsets,
                                   AFFY_ERROR *err)
{
  affy_uint32                       ds_rows, i, j;
  char                             *byte_base = base;
  FILE                             *fp;
  AFFY_CALVIN_DATASET              *metadata;
  const AFFY_CALVIN_COLUMN_MAPPING *cmap;

  assert(dio     != NULL);
  assert(base    != NULL);
  assert(offsets != NULL);

  metadata = dio->metadata;
  fp       = dio->calvin_io->fp;
  ds_rows  = metadata->num_rows;

  if ((start_row + num_rows) > ds_rows)
    AFFY_HANDLE_ERROR_VOID("too many rows requested",
                           AFFY_ERROR_BADPARAM,
                           err);

  if (fseek(fp, 
            dio->initial_offset + (start_row * dio->row_length), 
            SEEK_SET) != 0)
    AFFY_HANDLE_ERROR_VOID("I/O error reading Calvin file",
                           AFFY_ERROR_IO,
                           err);

  for (i = start_row; i < start_row + num_rows; i++)
  {
    for (j = 0; j < metadata->num_cols; j++)
    {
      cmap = mapping_for_column(offsets, metadata->columns[j].name);

      if (cmap != NULL)
      {
        void              *dest = byte_base + cmap->offset;
        affy_binary_reader read_func;
        
        read_func = calvin_type_table[metadata->columns[j].type].read_func;

        if (read_func(fp, dest) != 0)
          AFFY_HANDLE_ERROR_VOID("I/O error reading Calvin file",
                                 AFFY_ERROR_IO,
                                 err);
      }
      else /* skip column */
      {
        if (fseek(fp, metadata->columns[j].size, SEEK_CUR) != 0)
          AFFY_HANDLE_ERROR_VOID("I/O error reading Calvin file",
                                 AFFY_ERROR_IO,
                                 err); 
      }
    }

    byte_base += base_sz;

    pb_tick(pbs, 1,"");
  }
}

/* 
 * Reads an entire dataset into memory at once.  Not likely to be very
 * useful if you need to pass the result directly to some other
 * function since the results are stored in a union type.  See
 * affy_calvin_read_dataset_col() if you want to read a given column
 * (homogeneous type) into your own storage.
 *
 * You must have previously read in the dataset metadata, in order
 * for the column type information to be available.
 */
/* void affy_calvin_read_dataset_rows(AFFY_CALVINIO *cio,  */
/*                                    AFFY_CALVIN_DATASET *ds, */
/*                                    AFFY_ERROR *err) */
/* { */
/*   affy_uint32 i, j; */
  
/*   assert(cio != NULL); */
/*   assert(ds  != NULL); */

/*   /\* The actual data. *\/ */
/*   ds->data = calloc(ds->num_rows, sizeof(AFFY_CALVIN_DATA *)); */
/*   if (ds->data == NULL) */
/*     AFFY_HANDLE_ERROR_VOID("calloc failed", AFFY_ERROR_OUTOFMEM, err); */

/*   ds->data[0] = calloc(ds->num_rows * ds->num_cols, sizeof(AFFY_CALVIN_DATA)); */
/*   if (ds->data[0] == NULL) */
/*     AFFY_HANDLE_ERROR_VOID("calloc failed", AFFY_ERROR_OUTOFMEM, err); */

/*   for (i = 1; i < ds->num_rows; i++) */
/*     ds->data[i] = ds->data[i-1] + ds->num_cols; */

/*   for (i = 0; i < ds->num_rows; i++, ds->rows_read++) */
/*     for (j = 0; j < ds->num_cols; j++, ds->cols_read++) */
/*     { */
/*       /\*  */
/*        * Detect the case where the column claims to be a STRING, but is */
/*        * actually a WSTRING.  Self-describing formats are fun when */
/*        * they lie to you. */
/*        *\/ */
/*       if (ds->columns[j].type == AFFY_CALVIN_STRING) */
/*       { */
/*         affy_int32 strlen; */

/*         if (affy_read32_be(cio->fp, (void *)&strlen) != 0) */
/*           AFFY_HANDLE_ERROR_VOID("I/O error reading Calvin file",  */
/*                                  AFFY_ERROR_IO, */
/*                                  err); */

/*         if (fseek(cio->fp, -4, SEEK_CUR) == -1) */
/*           AFFY_HANDLE_ERROR_VOID("I/O error reading Calvin file",  */
/*                                  AFFY_ERROR_IO, */
/*                                  err); */

/*         if (((ds->columns[j].size - 4) / 2) == strlen) */
/*           /\* Actually a Unicode string, correct the type field. *\/ */
/*           ds->columns[j].type = AFFY_CALVIN_WSTRING; */
/*       } */

/*       affy_calvin_read_data(cio, &(ds->data[i][j]), ds->columns[j].type, err); */
/*       AFFY_CHECK_ERROR_VOID(err); */
/*     } */
/* } */

/*
 * Read the Calvin dataheader and all its parents, recursively.
 */
AFFY_CALVIN_DATAHEADER *affy_calvin_get_dataheader(AFFY_CALVINIO *cio, 
                                                   AFFY_ERROR *err)
{
  AFFY_CALVIN_DATAHEADER *dh;

  assert(cio != NULL);

  if (fseek(cio->fp, 10, SEEK_SET) != 0)
    AFFY_HANDLE_ERROR("I/O error reading Calvin file", 
                      AFFY_ERROR_IO,
                      err,
                      NULL);

  dh = h_malloc(sizeof(AFFY_CALVIN_DATAHEADER));
  if (dh == NULL)
    AFFY_HANDLE_ERROR("malloc failed", AFFY_ERROR_OUTOFMEM, err, NULL);

  process_dataheader(cio, dh, err);
  if (err->type != AFFY_ERROR_NONE)
  {
    affy_free_calvin_dataheader(dh);
   
    return (NULL);
  }

  return (dh);
}
