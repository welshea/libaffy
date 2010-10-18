
/**************************************************************************
 *
 * Filename:  error.c
 *
 * Purpose:   Error-handling utility routines.  To some extent this
 *            style of error handling was inspired by the GNU
 *            Scientific Library.
 *
 * Creation:  03/05/08
 *
 * Author:    Andrew M Hoerter
 *
 * Copyright: Copyright (C) 2007, Moffitt Cancer Center 
 *            All rights reserved. 
 *
 * Update History
 * --------------
 * 03/05/08: Creation (AMH)
 *
 **************************************************************************/

#include <utils.h>
#include <affy.h>

#define ERRBUF_SIZE 100

static void affy_die(AFFY_ERROR *err);

static char errbuf[ERRBUF_SIZE];

static void affy_die(AFFY_ERROR *err)
{
  assert(err != NULL);
  
  fprintf(stderr, 
          "ERROR: %s (%s) [%s:%d]\n", 
          err->descr, 
          affy_strerror(err->type), 
          err->module, 
          err->location);
  exit(EXIT_FAILURE);
}

void affy_clone_error(AFFY_ERROR *e1, AFFY_ERROR *e2)
{
  assert(e1 != NULL);
  assert(e2 != NULL);

  e1->type      = e2->type;
  e1->timestamp = e2->timestamp;
  e1->descr     = e2->descr;
  e1->module    = e2->module;
  e1->location  = e2->location;
  /* Typically it isn't useful to clone the handler */
}

AFFY_ERROR *affy_get_default_error(void)
{
  AFFY_ERROR *result;

  result = malloc(sizeof(AFFY_ERROR));
  if (result == NULL)
    /* How ironic. */
    return (NULL);

  result->type    = AFFY_ERROR_NONE;
  /* Other fields are undefined until filled in for an actual error. */
  result->handler = affy_die;

  return (result);
}

const char *affy_strerror(AFFY_ERROR_TYPE err)
{
  /* 
   * The values of extended/user-defined error codes are simply
   * printed verbatim.  If something more sophisticated is needed,
   * add it here.
   */
  if (err >= 100)
  {
    portable_snprintf(errbuf, ERRBUF_SIZE, "User-defined error %d", err);

    return (errbuf);
  }

  switch (err)
  {
    case AFFY_ERROR_NONE:         return "No error";
    case AFFY_ERROR_NOTFOUND:     return "File not found";
    case AFFY_ERROR_SYSPERM:      return "Permission denied";
    case AFFY_ERROR_NOTREADY:     return "Resource not ready";
    case AFFY_ERROR_LIMITREACHED: return "Limit/quota reached";
    case AFFY_ERROR_IO:           return "I/O error";
    case AFFY_ERROR_WRONGTYPE:    return "Type error";
    case AFFY_ERROR_OUTOFMEM:     return "Out of memory";
    case AFFY_ERROR_BADPARAM:     return "Bad parameter";
    case AFFY_ERROR_BADFORMAT:    return "Bad format";
    case AFFY_ERROR_UNKNOWN:      return "Unknown error";
    default:                      return "Undefined error";
  }
}
