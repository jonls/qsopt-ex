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
/** @defgroup EGtimer EGtimer 
 *
 * Here we implement types and functions for __timer functions
 *
 * @version 1.0.1 
 * @par History:
 * - 2013-04-18
 * 						- Add EGtimerZero static initializer
 * - 2006-01-25
 * 						- Fix compilation errors on sun, change includes accordingly and
 * 						code
 * - 2005-05-31
 * 						- Eliminate the definition of #EGwallClockTimer_t and replace it by
 * 						a macro definition that replace it by #EGtimer_t.
 * - 2004-01-20
 * 						- Add a 'wall clock' __timer (Renan-Marcos) type and functions
 * - 2003-05-08
 * 						- First Implementation
 * @note Up to now, this code will only work on linux machines, and maybe on
 * unix/posix systems.
 * */
/** @file 
 * @ingroup EGtimer */
/** @addtogroup EGtimer */
/** @{ */
/* ========================================================================= */

#ifndef __EG_TIMER_H__
#define __EG_TIMER_H__

#include <stdlib.h>
#include <errno.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "eg_macros.h"

/* ========================================================================= */
/** @brief Get system time.
 *
 * This function is for internal purposes only and should not be called from the
 * user space, it ask the (user) time from the system. 
 * @return the time (in seconds) stored as a double. */
#ifdef HAVE_GETRUSAGE
#define __EGzeit() ({\
	struct rusage __EGzeit_ru;\
	int __EGzeit_st = getrusage(RUSAGE_SELF,&__EGzeit_ru);\
	EXIT(__EGzeit_st,"getrusage failed with code error %d (%s)", errno, \
	strerror(errno));\
	(((double)__EGzeit_ru.ru_utime.tv_sec) + \
	((double)__EGzeit_ru.ru_utime.tv_usec)/1000000);})
#else
#ifdef HAVE_TIMES
#ifdef CLK_TCK
#define MACHINE_FREQ CLK_TCK
#else
#define MACHINE_FREQ HZ
#endif
#define __EGzeit() ({\
	struct tms __EGzeit_now;\
	times(&__EGzeit_now);\
	((double) __EGzeit_now.tms_utime)/((double) MACHINE_FREQ);})
#else
#error Your system does not have (or the configure script could not find)\
 getrusage nor times functions, and thus we are unable to provide \
 timing functions. Without them this library will not compile in this system
#endif
#endif

/* ========================================================================= */
/** @brief this structure holds a __timer structure */
typedef struct
{
	double time;	/**< hold the accumulated time */
	double stime;	/**< hols the last time when we start counting, this is only 
										 for internal purposes, the user should only use the 
										 field 'time' */
}
EGtimer_t;
/* ========================================================================= */
/** @brief static initializer */
#define EGtimerZero {.time=0,.stime=0}
/* ========================================================================= */
/** @brief This is done for backward compability, we used to define
 * EGwallClockTimer_t just as the normal __timer, so in reality we don't need
 * another type, but keep the name so that older code depending on this still
 * compiles. */
#define EGwallClockTimer_t EGtimer_t

/* ========================================================================= */
/** @brief Set a new starting time for the __timer.
 * @param __timer pointer to a EGtimer_t structure.
 * @return starting time (in seconds), and the type is a double. */
#define EGtimerStart(__timer) ({(__timer)->stime = __EGzeit();})

/* ========================================================================= */
/** @brief Stop a 'running' __timer and accumulate the run time.
 * @return the time elapsed since the last 'start' call (in seconds). */
#define EGtimerStop(__timer) ({(__timer)->time += __EGzeit() - (__timer)->stime;})

/* ========================================================================= */
/** @brief this function reset the accumulated time to zero */
#define EGtimerReset(__timer) ({(__timer)->time = 0;})

/* ========================================================================= */
/** @brief Set the starting time the current (wall) time.
 * @return the current wall time. */
#define EGwallClockTimerStart(__timer) ({(__timer)->stime = time(0);})

/* ========================================================================= */
/** @brief Stop a 'running' __timer and accumulate the (wall) runing time.
 * @return the wall time elapsed since the last initialization. */
#define EGwallClockTimerStop(__timer) ({\
	(__timer)->time += difftime(time(0),(__timer)->stime);})

/* ========================================================================= */
/** @brief Reset the accumulated time to zero */
#define EGwallClockTimerReset(__timer) EGtimerReset(__timer)

/* ========================================================================= */
/** @}
 * end of eg_timer.h */
#endif
