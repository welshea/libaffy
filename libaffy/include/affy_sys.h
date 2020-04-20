
/**************************************************************************
 *
 * Filename: affy_sys.h
 * 
 * Purpose: System-dependent glue.
 *
 * Creation date: 8 April, 2005
 *
 * Author: Andrew M. Hoerter
 *
 * Comments: 
 *
 *
 * Update History
 * --------------
 * 4/08/05: Creation (AMH)
 * 4/20/20: AFFY_POSIX commented as now obsolete (EAW)
 *
 **************************************************************************/

#ifndef AFFY_SYS_H
#define AFFY_SYS_H

/* 
 * Define this if you know your system has C99 integer types.  But
 * don't do it here, use the build system.
 */
/* #define AFFY_HAVE_STDINT */

/*
 * Defined if the library is being compiled on a POSIX-like system.
 * Note that this currently only matters for list_files.c, which is
 * used only by the command-line driver programs, and even they don't
 * truly require it.
 *
 * Also note that the build system is typically responsible for
 * setting this.
 *
 * Update: 2020-04-20
 *  removed AFFY_POSIX from list_files.c entirely, since POSIX is now
 *  assumed to be the default, with all others as special cases
 */
/* #define AFFY_POSIX */

#endif
