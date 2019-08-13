
/**************************************************************************
 *
 * Filename: affy_basetypes.h
 * 
 * Purpose: Define some basic types for use within libaffy.
 *
 * Creation date: 4 April, 2005
 *
 * Author: Andrew M. Hoerter
 *
 *
 * Update History
 * --------------
 * 4/4/05: Creation (AMH)
 *
 **************************************************************************/

#ifndef AFFY_BASETYPES_H
#define AFFY_BASETYPES_H

#include <affy_sys.h>

/* Basic types. */

/**************************************************************************/

#if (defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 199901L)) || defined(AFFY_HAVE_STDINT)

/* Use C99 types. */
#include <stdint.h>
#include <inttypes.h>

typedef int8_t   affy_int8;
typedef uint8_t  affy_uint8;
typedef int16_t  affy_int16;
typedef uint16_t affy_uint16;
typedef int32_t  affy_int32;
typedef uint32_t affy_uint32;

#define AFFY_PRNd8  PRId8
#define AFFY_PRNu8  PRIu8
#define AFFY_PRNu16 PRIu16
#define AFFY_PRNu32 PRIu32
#define AFFY_PRNd16 PRId16
#define AFFY_PRNd32 PRId32

#else

/* Calculate base types as best we can without explicit help. */
#if (UCHAR_MAX == 255U)
typedef unsigned char affy_uint8;
typedef signed char affy_int8;
#define AFFY_PRNd8 "d"
#define AFFY_PRNu8 "u"
#else
#  error No integral type is 8 bits on this platform.
#endif

#if (UINT_MAX == 65535U)
typedef unsigned int affy_uint16;
typedef int affy_int16;

#define AFFY_PRNu16 "u"
#define AFFY_PRNd16 "d"
#elif (USHRT_MAX == 65535U)
typedef unsigned short affy_uint16;
typedef short affy_int16;

#define AFFY_PRNu16 "u"
#define AFFY_PRNd16 "d"
#elif (UCHAR_MAX == 65535U)
typedef unsigned char affy_uint16;

#define AFFY_PRNu16 "u"
#else
#  error No integral type is 16 bits on this platform.
#endif

#if (UINT_MAX == 4294967295UL)
typedef unsigned int affy_uint32;
typedef int affy_int32;

#define AFFY_PRNu32 "u"
#define AFFY_PRNd32 "d"
#elif (ULONG_MAX == 4294967295UL)
typedef unsigned long affy_uint32;
typedef long affy_int32;

#define AFFY_PRNu32 "lu"
#define AFFY_PRNd32 "l"
#elif (USHRT_MAX == 4294967295UL)
typedef unsigned short affy_uint32;
typedef short affy_int32;

#define AFFY_PRNu32 "u"
#define AFFY_PRNd32 "d"
#elif (UCHAR_MAX == 4294967295UL)
typedef unsigned char affy_uint32;
typedef signed char affy_int32;

#define AFFY_PRNu32 "u"
#define AFFY_PRNd32 "d"
#else
#  error No integral type is 32 bits on this platform.
#endif

#endif

typedef time_t affy_timestamp;

/* XXX -- detect these somehow */
typedef float affy_float32;
typedef double affy_float64;


#if 0             /* get rid of it, so -Wall won't complain */
#if !defined(NDEBUG)
/* 
 * A "canary" declaration to test the type sizes.  If the compiler
 * barfs on this, there's big trouble in little China.
 */
static union
{
  char int8_t_incorrect[sizeof(affy_uint8) == 1];
  char int16_t_incorrect[sizeof(affy_int16) == 2];
  char uint16_t_incorrect[sizeof(affy_uint16) == 2];
  char int32_t_incorrect[sizeof(affy_int32) == 4];
  char uint32_t_incorrect[sizeof(affy_uint32) == 4];
  char float32_t_incorrect[sizeof(affy_float32) == 4];
  char float64_t_incorrect[sizeof(affy_float64) == 8];
} affy_typecanary;
#endif
#endif

#endif
