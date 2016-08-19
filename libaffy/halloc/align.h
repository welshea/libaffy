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
 *      2016-08-09 EAW -- renamed max_align_t to halloc_max_align_t
 *       so that it builds on OSX with clang and Xcode
 *       (clang/XCode defined their own conflicting max_align_t).
 */

#ifndef _LIBP_ALIGN_H_
#define _LIBP_ALIGN_H_

/*
 *	a type with the most strict alignment requirements
 */
union max_align
{
	char   c;
	short  s;
	long   l;
	int    i;
	float  f;
	double d;
	void * v;
	void (*q)(void);
};

typedef union max_align halloc_max_align_t;

#endif

