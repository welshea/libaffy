
/**************************************************************************
 *
 * Filename:  list_files.c
 *
 * Purpose:   Find files with a given extension.
 *
 * Creation:  04/08/05
 *
 * Author:    Steven Eschrich
 *
 * Copyright: Copyright (C) 2007, Moffitt Cancer Center.
 *            All rights reserved.
 *
 * Update History
 * --------------
 * 04/08/05: Imported/repaired from old libaffy (AMH)
 * 09/11/07: Cosmetic updates/cleanups (AMH)
 * 03/11/08: New error handling scheme (AMH)
 * 09/20/10: Pooled memory allocator (AMH)
 * 09/06/12: Swapped POSIX/WIN32 #ifdef order to better handle cygwin (EAW)
 * 09/12/19: fixed *new assignment typo in WIN32 list_files() (EAW)
 *
 **************************************************************************/

#include <affy.h>
#include <utils.h>

struct name_n
{
  char          *name;
  struct name_n *next;
};

#if defined(AFFY_WIN32_ENV)
/* Code mostly courtesy of Dan Maas from Joel on Software. */
#include <windows.h>
char **affy_list_files(char *directory, char *extension, AFFY_ERROR *err)
{
  WIN32_FIND_DATA fdata;
  HANDLE          dhandle;
  struct name_n  *names = NULL, *np;
  char            buf[MAX_PATH];
  char          **filelist = NULL;
  int             num_chips = 0, i;

  assert(directory != NULL);
  assert(extension != NULL);
  
  portable_snprintf(buf, sizeof(buf), "%s\\*%s", directory, extension);
  
  if ((dhandle = FindFirstFile(buf, &fdata)) == INVALID_HANDLE_VALUE) 
  {
    AFFY_HANDLE_ERROR("FindFirstFile() failed",
                      AFFY_ERROR_UNKNOWN,
                      err,
                      NULL);
  }

  do
  {
    /* If this is reached, FindFirstFile() found at least one match */
    struct name_n *new  = h_malloc(sizeof(struct name_n));
    unsigned int   slen = strlen(fdata.cFileName) + strlen(directory) + 2;

    if (new == NULL)
      AFFY_HANDLE_ERROR_GOTO("malloc failed", AFFY_ERROR_OUTOFMEM, err, err);

    hattach(names, new);
    new->next = names;
    names     = new;

    new->name = h_suballoc(new, slen);
    if (new->name == NULL)
      AFFY_HANDLE_ERROR_GOTO("malloc failed", AFFY_ERROR_OUTOFMEM, err, err);
    
    sprintf(new->name, "%s\\%s", directory, fdata.cFileName);
    num_chips++;
  } while (FindNextFile(dhandle, &fdata));

  if (GetLastError() != ERROR_NO_MORE_FILES) 
  {
    AFFY_HANDLE_ERROR_GOTO("FindNextFile() failed",
                           AFFY_ERROR_UNKNOWN,
                           err,
                           err);
  }
   
  i = 0;
  for (np = names; np != NULL; np = np->next)
  {
    assert(i < num_chips);
    
    filelist[i++] = np->name;
    hattach(np->name, filelist);
  }
 
  filelist[num_chips] = NULL;

  FindClose(dhandle);
  h_free(names);

  return (filelist);

err:
  FindClose(dhandle);

  h_free(filelist);
  h_free(names);

  return (NULL);
}
#else
# include <sys/types.h>
# include <dirent.h>
# include <unistd.h>

char **affy_list_files(char *directory, char *extension, AFFY_ERROR *err)
{
  struct dirent  *de;
  struct name_n  *names = NULL, *np;
  DIR            *d;
  char          **filelist = NULL;
  unsigned int    num_chips = 0, i;

  assert(directory != NULL);
  assert(extension != NULL);

  d = opendir(directory);
  if (d == NULL)
    AFFY_HANDLE_ERROR("opendir() failed", AFFY_ERROR_NOTFOUND, err, NULL);

  while ((de = readdir(d)) != NULL)
  {
    if (endsWith(de->d_name, extension))
    {
      struct name_n *new  = h_malloc(sizeof(struct name_n));
      unsigned int   slen = strlen(de->d_name) + strlen(directory) + 2;

      if (new == NULL)
	AFFY_HANDLE_ERROR_GOTO("malloc failed", 
                               AFFY_ERROR_OUTOFMEM, 
                               err, 
                               cleanup);

      if (names != NULL)
        hattach(names, new);

      new->next = names;
      names     = new;

      new->name = h_suballoc(new, slen);
      if (new->name == NULL)
	AFFY_HANDLE_ERROR_GOTO("malloc failed", 
                               AFFY_ERROR_OUTOFMEM, 
                               err, 
                               cleanup);

      sprintf(new->name, "%s/%s", directory, de->d_name);
      num_chips++;
    }
  }

  if (num_chips == 0)
    goto cleanup;

  filelist = h_malloc((num_chips + 1) * sizeof(char *));
  if (filelist == NULL)
    AFFY_HANDLE_ERROR_GOTO("malloc failed", 
                           AFFY_ERROR_OUTOFMEM, 
                           err, 
                           cleanup);

  
  i = 0;
  for (np = names; np != NULL; np = np->next)
  {
    assert(i < num_chips);
    
    filelist[i++] = np->name;
    hattach(np->name, filelist);
  }

  /* Last spot in list is null */
  filelist[num_chips] = NULL;

  closedir(d);
  h_free(names);

  return (filelist);

cleanup:
  closedir(d);
  h_free(filelist);
  h_free(names);

  return (NULL);
}
#endif
