
/**************************************************************************
 *
 * Filename:  is_control_probe.c
 *
 * Purpose:   Probe testing predicate for control probes.
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
 * 2005-04-14: Imported/repaired from old libaffy (AMH)
 * 2019-03-15: Added probeset and string functions, more control strings (EAW)
 * 2019-08-13: removed "unsigned" from affy_is_control_string() chars (EAW)
 *
 **************************************************************************/

#include <utils.h>
#include <affy.h>


bool affy_is_control_string(char *string)
{
  char *str_lower = NULL;
  char *sptr, *sptr2;
  unsigned int length;

  assert(string != NULL);
  
  length = strlen(string);

  /* strcasestr() and strcasecmp() are non-standard BSD/POSIX extensions,
   * so the simplest solution is to make a lowercase copy ahead of time
   */

  /* allocate lowercase string */
  if (length)
  {
    str_lower = (char *) malloc((length+1) * sizeof(char));

    /* this should never happen, we must have run out of memory */
    if (!str_lower)
      return false;
    
    str_lower[length] = '\0';
  }
  else
    return(false);
  
  /* lowercase the string */
  sptr  = str_lower;
  sptr2 = string;

  while (*sptr2 != '\0')
    *sptr++ = tolower(*sptr2++);

  /* AFFX */
  if (strncmp(string, "AFFX", 4) == 0)
  {
    if (str_lower)
      free(str_lower);
    return (true);
  }

  /* spike in */
  if ((sptr = strstr(str_lower, "spike")))
  {
    sptr += 5;
    if (strncmp(sptr, "in",  2) == 0 ||
        strncmp(sptr, "-in", 3) == 0 ||
        strncmp(sptr, "_in", 3) == 0 ||
        strncmp(sptr, " in", 3) == 0)
    {
      if (str_lower)
        free(str_lower);
      return (true);
    }
  }
  /* control */
  if (strstr(str_lower, "control"))
  {
    if (str_lower)
      free(str_lower);
    return (true);
  }

  if (str_lower)
    free(str_lower);
  
  return (false);
}

bool affy_is_control_probe(AFFY_PROBE *p_probe)
{
  assert(p_probe     != NULL);
  assert(p_probe->ps != NULL);

  /* 
   * Any probeset name beginning with AFFX is designated a control
   * probeset, and probes belonging to it are control probes.
   */
  if (affy_is_control_string(p_probe->ps->name))
    return (true);
  else
    return (false);
}


bool affy_is_control_probeset(AFFY_PROBESET *ps)
{
  assert(ps != NULL);

  /* 
   * Any probeset name beginning with AFFX is designated a control
   * probeset, and probes belonging to it are control probes.
   */
  if (affy_is_control_string(ps->name))
  {
    return (true);
  }
  else
    return (false);
}
