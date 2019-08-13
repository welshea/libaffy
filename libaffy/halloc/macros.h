/*
 *	Copyright (c) 2004-2010 Alex Pankratov. All rights reserved.
 *
 *	Hierarchical memory allocator, 1.2.1
 *	http://swapped.cc/halloc
 */

/*
 *	The program is distributed under terms of BSD license. 
 *	You can obtain the copy of the license by visiting:
 *	
 *	http://www.opensource.org/licenses/bsd-license.php
 *
 *      2019-08-13: (EAW) pointer arithmetic on void pointer is a GNU
 *       extension and is otherwise illegal.
 *       Changed offending (void *) to (char *) in structof()
 */

#ifndef _LIBP_MACROS_H_
#define _LIBP_MACROS_H_

#include <stddef.h>  /* offsetof */

/*
 	restore pointer to the structure by a pointer to its field
 */
#define structof(p,t,f) ((t*)(- offsetof(t,f) + (char*)(p)))

/*
 *	redefine for the target compiler
 */

/* Added by AMH for portability */
#if defined(__GNUC__) && !defined(__STRICT_ANSI__)
/* GCC-isms. */
# define INLINE inline
#elif (defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 199901L))
/* C99-isms. */
# define INLINE inline
#else
/* The rest of the world is C89. */
# define INLINE
#endif

#define static_inline static INLINE

#endif

