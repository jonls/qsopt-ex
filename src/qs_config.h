/* Taken from EGlib 10-09-26 */
/* EGlib "Efficient General Library" provides some basic structures and
 * algorithms commons in many optimization algorithms.
 *
 * Copyright (C) 2005 Daniel Espinoza and Marcos Goycoolea.
 * 
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation; either version 2.1 of the License, or (at your
 * option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public 
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA 
 * */
/* ========================================================================= */
/** Main Configuration for the library, as debug levels and so on
 * 
 * @par History:
 * - 2010-09-26
 * 						- Addapted for QSopt_ex
 * - 2010-08-13
 * 						- Add suport for autoconf and configure for header and feature
 * 						selection
 * - 2006-01-27
 *					- Handle some problems with stdint.h in SUN
 * - 2005-08-17
 * 					- Set memory aligment to 8 bits
 * - 2003-06-02
 * 					- First Implementation
 * @version 1.1.1 
 * */
/* ========================================================================= */
#ifndef __QS_CONFIG_H__
#define __QS_CONFIG_H__
#include "config.h"
#ifdef HAVE_STDIO_H
# include <stdio.h>
#endif
#ifdef HAVE_SYS_SOCKET_H
# include <sys/socket.h>
#endif
#ifdef HAVE_SYS_TYPES_H
# include <sys/types.h>
#endif
#ifdef HAVE_SYS_STAT_H
# include <sys/stat.h>
#endif
#ifdef STDC_HEADERS
# include <stdlib.h>
# include <stddef.h>
#else
# ifdef HAVE_STDLIB_H
# include <stdlib.h>
# endif
#endif
#ifdef HAVE_STRING_H
# if !defined STDC_HEADERS && defined HAVE_MEMORY_H
# include <memory.h>
# endif
# include <string.h>
#endif
#ifdef HAVE_STRINGS_H
# include <strings.h>
#endif
#ifdef HAVE_INTTYPES_H
# include <inttypes.h>
#endif
#ifdef HAVE_STDINT_H
# include <stdint.h>
#endif
#ifdef HAVE_UNISTD_H
# include <unistd.h>
#endif
#ifdef HAVE_ERRNO_H
# include <errno.h>
#endif
#ifdef HAVE_LIMITS_H
# include <limits.h>
#endif
#ifdef HAVE_MATH_H
# include <math.h>
#endif
#ifdef HAVE_FLOAT_H
# include <float.h>
#endif
#ifdef HAVE_GETOPT_H
# include <getopt.h>
#endif
#ifdef HAVE_NETINET_IN_H
# include <netinet/in.h>
#endif
#ifdef HAVE_NETINET_TCP_H
# include <netinet/tcp.h>
#endif
#ifdef HAVE_NETDB_H
# include <netdb.h>
#endif
#ifdef HAVE_FCNTL_H
# include <fcntl.h>
#endif
#ifdef HAVE_SYS_RESOURCE_H
# include <sys/resource.h>
#endif
#ifdef HAVE_SYS_PARAM_H
# include <sys/param.h>
#endif
#ifdef TIME_WITH_SYS_TIME
# include <sys/time.h>
# include <time.h>
#else
# ifdef HAVE_SYS_TIME_H
# include <sys/time.h>
# else
# include <time.h>
# endif
#endif
#ifdef HAVE_SYS_TIMES_H
# include <sys/times.h>
#endif
#ifdef HAVE_STDARG_H
# include <stdarg.h>
#endif
#ifdef HAVE_SYS_UTSNAME_H
# include <sys/utsname.h>
#endif
#ifdef HAVE_SIGNAL_H
# include <signal.h>
#endif
#ifdef HAVE_SETJMP_H
# include <setjmp.h>
#endif
/* ========================================================================= */
/** @brief if no gmp support, we do not include gmp.h, if on the otherhand, we
 * have libgmp, we MUST have gmp.h */
#ifdef HAVE_LIBGMP
# if HAVE_LIBGMP
#  ifdef HAVE_GMP_H
#   include <gmp.h>
#  else
#   error Must have gmp.h for compiling QSopt_ex
#  endif
# else
#  error Must have gmp.h for compiling QSopt_ex
# endif
# else
#  error Must have gmp.h for compiling QSopt_ex
#endif
/* ========================================================================= */
/** @brief assert Debug options definitions, by defoult set on */
#ifndef DEBUG
#warning you should define DEBUG, assuming it to be 1
#define DEBUG 1
#endif

/* ========================================================================= */
/** @brief assert Verbose options definition, by default set on */
#ifndef VERBOSE_LEVEL
#warning you should define VERBOSE_LEVEL, assuming it to be 1
#define VERBOSE_LEVEL 1
#endif
/* ========================================================================= */
#ifdef HAVE_EGLIB_H
# if HAVE_EGLIB_H
#  include "EGlib.h"
# else
#  error You must have EGlib.h for compilation
# endif
#else
#  error You must have EGlib.h for compilation
#endif
/* ========================================================================= */
/** @brief define version function name */
void QSopt_ex_version(void);
#endif
