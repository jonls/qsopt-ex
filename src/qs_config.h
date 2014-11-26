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

#if ! defined(_XOPEN_SOURCE) || _XOPEN_SOURCE < 600
# define _XOPEN_SOURCE 600
#endif

#define EG_NEWLINE "\n"

/* ========================================================================= */
/** @brief see if we have both pthread library and header file available, if
 * so, define HAVE_EG_THREAD */
#if defined HAVE_PTHREAD_H && HAVE_PTHREAD_H && defined HAVE_LIBPTHREAD && HAVE_LIBPTHREAD
#define HAVE_EG_THREAD 1
#else
#define HAVE_EG_THREAD 0
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
#endif
