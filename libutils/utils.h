
#ifndef __UTILS_H_
#define __UTILS_H_

/* 
 *  utils.h
 *
 *  Adapted from the original libutils, written by:
 *
 *       Steven Eschrich
 *	 eschrich@csee.usf.edu
 *	 
 *  ANSI-ified by Andrew M. Hoerter.  Where possible I have kept to the
 *  original libutils interfaces, but in some cases changes have been
 *  made to make the code standards-conforming.  Any resulting bugs are
 *  mine alone.
 *
 *  04/08/13 -- added stem_from_filename_safer(), to use instead of the old
 *              stem_from_filename(), since sample names in spreadsheets
 *              can have .A01, etc. extensions as part of the true sample
 *              name, and thus should not be removed. (EAW)
 *
 */

#include "endian_config.h"
#include "sys_config.h"

/* System defines needed somewhere below */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <assert.h>
#include <stdarg.h>

/* Extra includes */
/***********************************************************
   bitstring.h may be included in some systems, so we will  
   check that no one else has included it, then include our
   own copy of it.
 **********************************************************/
#ifndef bitstr_t
#include "bitstring_.h"
#endif

/********************************************************
 * Compiler-specific knobs.
 * 
 ********************************************************/

#if defined(__GNUC__) && !defined(__STRICT_ANSI__)
/* GCC-isms. */
#define INLINE inline
/* Apparently G++ doesn't allow the unused attribute on functions. */
 #if !defined(__cplusplus)
  #define UNUSED __attribute ((unused))
  #include "stdbool_.h"
 #else
  #define UNUSED
 #endif

#elif (defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 199901L))
/* C99-isms. */
# define INLINE inline
# define UNUSED
 /* Include the stdbool directly */
# include <stdbool.h>

#else
/* The rest of the world is C89. */
# define INLINE
# define UNUSED
#   if !defined(__cplusplus)
#     include "stdbool_.h"
#   endif
#endif

/********************************************************
 * Constants defined within libutils.
 * 
 ********************************************************/

/* The buffer size generally to create statically, a large number */
#define MAXBUF        (10*1024)

/* If arrays grow, by what increment */
#define REALLOC_INCR  (1024)

/* 
 * The character used for directory separators in filenames.  This should
 * really be detected somehow, but at least this is more portable than
 * what it replaces.
 */
#define DIRECTORY_SEPARATOR '/'

/********************************************************
 * Variables defined within the library.
 * 
 ********************************************************/

/* An indicator of the endian-ness of the system */
extern const int  libutils_big_endian;
extern const char libutils_version[];

/********************************************************
 * Function prototypes.
 * 
 ********************************************************/
#ifdef __cplusplus
extern "C"
{
#endif

  char *trim(char *p);
  int   split(char *str, char **kv, int split_val, int maxsplit);

  /* Get the next line from an open file. */
  char *utils_getline(FILE *fp);

  /* Simple interface for parsing floats. */
  int parsefloat(char *p_str, double *p_dest);

  /* Compare two doubles (for sorting) */
  int dcompare(const void *p1, const void *p2);
  int icompare(const void *p1, const void *p2);
  int fcompare(const void *p1, const void *p2);

  /* Matrix allocation */
  double **create_matrix(unsigned int rows, unsigned int columns);
  void     free_matrix(double **m);

#ifdef HAVE_STRDUP
  char *strdup(const char *s);
#endif
#ifdef HAVE_STRCASECMP
  int   strcasecmp(const char *s1, const char *s2);
#endif
#ifdef HAVE_STRNDUP
  char *strndup(const char *s,size_t n);
#endif
#ifdef HAVE_STRNLEN
  char *strnlen(const char *string, size_t maxlen);
#endif
  char *strchrnul(const char *s, int c_in);
  void *mempcpy(void *dest, const void *src, size_t n);
  int   portable_snprintf(char *str, 
                          size_t str_m, 
                          const char *fmt, /*args*/ ...);
  int   portable_vsnprintf(char *str, 
                           size_t str_m, 
                           const char *fmt, 
                           va_list ap);
#ifdef __cplusplus
};
#endif

/*****************************************************************************
   Macros
 
 	 Various useful macros (die, warn, info, etc). These macros ultimately
 	 call implementations after doing some massaging.

*****************************************************************************/

/* Sensible string comparison */
#define STREQ(a,b)   ( !strcasecmp(a,b) )

/* Many logging utilities are defined separately */
#include "util_log.h"

#define MALLOC(size)         MALLOC_impl(size)
#define CALLOC(n, size)      CALLOC_impl(n, size)
#define REALLOC(q, n)        REALLOC_impl(q, n)
#define FOPEN(fname, mode)   FOPEN_impl(fname, mode)
#ifdef ZLIB_H
# define GZOPEN(fname, mode) GZOPEN_impl(fname, mode)
#endif

/* These simpler min/max macros are NOT side-effect safe! */
#define min_macro(X, Y)                     \
        ((X < Y) ? X : Y)

#define max_macro(X, Y)                     \
        ((X > Y) ? X : Y )

/* A perlism - chomp the last character off of the string */
#define chomp(q) \
      do { \
         register char *__p=q; \
         register int __plen=strlen(__p); \
         if ( __p != 0 && __plen > 0 ) \
	   __p[__plen-1]=0; \
      } while(0)

/* A perlism, push a new element on the given array */
#define push(array, element, array_count, array_type) \
      do { \
         register int __tmp_count=array_count; \
         set_locus(__FILE__, __LINE__);         \
         if ( __tmp_count % REALLOC_INCR == 0 ) \
             array=(array_type *)realloc(array, \
                                            (__tmp_count+REALLOC_INCR)* \
                                              sizeof(array_type)); \
         if ( ! array ) { \
           die("Memory allocation failure");\
         }\
         array[__tmp_count]=element; \
      } while (0)

/* A perlism, pop the last element off of the array */
#define pop(array, array_count, array_type) \
        ({ \
          register int __tmp_count=array_count; \
          array_type __retval=0; \
          if ( __tmp_count-1 >= 0 ) \
             __retval=array[__tmp_count-1]; \
          __retval; \
        })

/*****************************************************************************
  SECTION2: Inline functions

  If you have GCC or a C99 compiler, these ought to be real inlines;
  otherwise they will be small static functions which a smart compiler
  will hopefully inline anyway.  If all else fails, oh well.
*****************************************************************************/

/* MALLOC space, exiting the program on failure */
static INLINE void *MALLOC_impl(size_t size)
{
  register void *result;

  if ((result = malloc(size)) == NULL)
  {
    die("error allocating memory");
    return (NULL);              /* To silence compiler diagnostics. */
  }
  else
    return (result);
} UNUSED

/* CALLOC space, exiting the program on failure */
static INLINE void *CALLOC_impl(size_t n, size_t size)
{
  register void *result;

  if ((result = calloc(n, size)) == NULL)
  {
    die("error allocating memory");
    return (NULL);              /* To silence compiler diagnostics. */
  }
  else
    return (result);
} UNUSED

/* REALLOC space, exiting the program on failure */
static INLINE void *REALLOC_impl(void *q, size_t n)
{
  register void *result;

  if ((result = realloc(q, n)) == NULL)
  {
    die("error re-allocating memory");
    return (NULL);              /* To silence compiler diagnostics. */
  }
  else
    return (result);
} UNUSED

/* Open a file, exiting program on failure. */
static INLINE void *FOPEN_impl(const char *filename, const char *mode)
{
  FILE *result;

  if ((result = fopen(filename, mode)) == NULL)
  {
    die("Error opening %s (mode %s): %s", filename, mode,
             strerror(errno));
    return (NULL);              /* To silence compiler diagnostics. */
  }
  else
    return (result);
} UNUSED

#ifdef ZLIB_H
/* Open a file, exiting program on failure. */
static INLINE void *GZOPEN_impl(const char *filename, const char *mode)
{
  gzFile *result;

  if ((result = gzopen(filename, mode)) == NULL)
  {
    die("Error opening %s (mode %s): %s", filename, mode,
             strerror(errno));
    return (NULL);              /* To silence compiler diagnostics. */
  }
  else
    return (result);
} UNUSED
#endif

/* The bit-string based operations have been taken over by
   bitstring.h, therefore we only provide backward compatability
   for old functions.
*/
#define numbytes		bitstr_size
#define bget(v,index)		(bit_test(v,index)?1:0)
#define bset(v,index,value)	do { \
					if (value!=0){bit_set(v,index);} \
					else {bit_clear(v,index);} \
				} while (0)

/*
 * Given a pathname, return the filestem, i.e. the result of stripping
 * the directory prefix and file extension (if any)
 *
 * Example: /a/b/c/foo.txt  -->  foo
 *
 * Storage allocated for the resulting string must be freed by the caller.
 */
static INLINE char *stem_from_filename(const char *p)
{
  const char *q;
  char       *result, *r;
                
  /* Strip out leading path information */
  q = strrchr(p, DIRECTORY_SEPARATOR);

  /* Either we found it (increment) or not */
  if (q != NULL)
    q++;
  else
    q = p;

  /* 
   * q now points to the filestem plus its extension, if any.  Next we
   * make a copy in preparation for stripping the extension and returning
   * the finished result.
   */
  result = (char *)MALLOC(strlen(q) + 1);
  strcpy(result, q);
                
  /* Strip out extension as well */
  r = strrchr(result, '.');
  if (r != NULL) 
    *r = '\0';

  return (result);
} UNUSED


/*
 * Given a pathname, return the filestem, i.e. the result of stripping
 * the directory prefix and file extension (if any)
 *
 * *ONLY* removes .CEL and .TXT/.TEXT file extensions at the end, since
 *  some samplenames have .A01, etc. at the end when reading from
 *  spreadsheets
 *
 * Example: /a/b/c/foo.txt  -->  foo
 *
 * Storage allocated for the resulting string must be freed by the caller.
 */
static INLINE char *stem_from_filename_safer(const char *p)
{
  const char *q;
  char       *result, *r, *s, *lc_ext;
  
  /* Strip out leading path information */
  q = strrchr(p, DIRECTORY_SEPARATOR);

  /* Either we found it (increment) or not */
  if (q != NULL)
    q++;
  else
    q = p;

  /* 
   * q now points to the filestem plus its extension, if any.  Next we
   * make a copy in preparation for stripping the extension and returning
   * the finished result.
   */
  result = (char *)MALLOC(strlen(q) + 1);
  strcpy(result, q);

  /* Find file extension */
  r = strrchr(result, '.');

  /* Only strip .CEL and .TXT */
  if (r != NULL)
  {
    /* copy file extension and lowercase it for later comparison */
    lc_ext = strdup(r);
    s = lc_ext;
    while (*s)
    {
      *s = tolower(*s);
      s++;    /* shut the mac compiler up, instead of *s++ = tolower(*s) */
    }
    
    /* remove if .cel or .txt */
    if (strcmp(lc_ext, ".cel") == 0 ||
        strcmp(lc_ext, ".txt") == 0 ||
        strcmp(lc_ext, ".text") == 0)
    {
      *r = '\0';
    }
    
    free(lc_ext);
  }

  return (result);
} UNUSED


/* Remove comments from line, given comment char c */
static INLINE int strip_comments(char *p, char c)
{
  if (p==NULL)
    return -1;
  for (; *p != '\0'; p++)
  {
    if (*p == c)
    {
      *p = '\0';
      break;
    }
  }
  return 0;
} UNUSED

/* Trim leading spaces (whitespace) */
static INLINE char *ltrim(char *p)
{
  if (p==NULL)
    return p;
  for (; *p  && isspace((int)(*p)); p++);
  return p;
} UNUSED

/* Seek from the end of p backwards to the first non-blank character.
   This is a bit difficult since array bounds must be kept in mind. */
static INLINE char *rtrim(char *p)
{
  char *end;
  if (p==NULL)
    return 0;
  /* Seek back to beginning or first non-blank */
  for (end=p+strlen(p)-1; end > p && isspace((int)(*end)); end--)
     *end='\0';
  /* Are we at the beginning? If so, we might either have a blank
     or a non-blank. */
  if (end == p)
  {
    if (isspace((int)(*end)))
      *end = '\0';                 /* \n or the like */
    else if (*end != '\0')
      *(end + 1) = '\0';           /* non-empty string, zero next char */
  }
  return end;
} UNUSED

/*
  endsWith(). Like the java sort of version, does it end with the
  given string. See also startsWith(). Does case insensitive comparison.
*/
static INLINE int endsWith(char *str, char *sub)
{
  int N = strlen(str);
  int last = strlen(sub);

  if (last > N)
    return 0;

  return (!strcasecmp(str + N - last, sub));
} UNUSED


#endif
