/****************************************************************************/
/*                                                                          */
/*  This file is part of QSopt_ex.                                          */
/*                                                                          */
/*  (c) Copyright 2006 by David Applegate, William Cook, Sanjeeb Dash,      */
/*  and Daniel Espinoza                                                     */
/*                                                                          */
/*  Sanjeeb Dash ownership of copyright in QSopt_ex is derived from his     */
/*  copyright in QSopt.                                                     */
/*                                                                          */
/*  This code may be used under the terms of the GNU General Public License */
/*  (Version 2.1 or later) as published by the Free Software Foundation.    */
/*                                                                          */
/*  Alternatively, use is granted for research purposes only.               */
/*                                                                          */
/*  It is your choice of which of these two licenses you are operating      */
/*  under.                                                                  */
/*                                                                          */
/*  We make no guarantees about the correctness or usefulness of this code. */
/*                                                                          */
/****************************************************************************/

/* RCSINFO $Id: zeit.c,v 1.2 2003/11/05 16:47:22 meven Exp $ */
/****************************************************************************/
/*                                                                          */
/*  This file is part of CONCORDE                                           */
/*                                                                          */
/*  (c) Copyright 1995--1999 by David Applegate, Robert Bixby,              */
/*  Vasek Chvatal, and William Cook                                         */
/*                                                                          */
/*  Permission is granted for academic research use.  For other uses,       */
/*  contact the authors for licensing options.                              */
/*                                                                          */
/*  Use at your own risk.  We make no guarantees about the                  */
/*  correctness or usefulness of this code.                                 */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/*                        TIMING FUNCTIONS                                  */
/*                                                                          */
/*                            TSP CODE                                      */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: Summer 1994  (cofeb16)                                            */
/*        December 1997 (dla)                                               */
/*                                                                          */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  double ILLutil_zeit (void)                                              */
/*        - To measure cpu time.                                            */
/*    To use this, set double t = ILLutil_zeit (), run the function you     */
/*    want to time, then compute ILLutil_zeit () - t.                       */
/*                                                                          */
/*  double ILLutil_real_zeit (void)                                         */
/*    - To measure wall clock time.                                         */
/*                                                                          */
/*    To use this, set double t = ILLutil_real_zeit (), run the function    */
/*    you want to time, then compute ILLutil_real_zeit () - t.              */
/*                                                                          */
/*  void ILLutil_init_timer (ILLutil_timer *t, const char *name)            */
/*    - Initializes a ILLutil_timer, and gives it a name.                   */
/*    - The name is silently truncated if it is too long.                   */
/*                                                                          */
/*  void ILLutil_start_timer (ILLutil_timer *t)                             */
/*    - Starts the timer.                                                   */
/*                                                                          */
/*  void ILLutil_suspend_timer (ILLutil_timer *t)                           */
/*    - Suspends the timer.  Similar to ILLutil_stop_timer, but doesn't     */
/*      count a call, and doesn't output.                                   */
/*                                                                          */
/*  void ILLutil_resume_timer (QSutil_timer *t)                             */
/*    - Resumes the timer after a suspend.                                  */
/*                                                                          */
/*  double ILLutil_stop_timer (ILLutil_timer *t, int printit)               */
/*    - Stops the timer, and returns the time since the last start.         */
/*    - if printit == 1, outputs the time spent.                            */
/*    - if printit == 2, outputs the time spent only if nonzero             */
/*    - if printit == 3,4, like 1,2, except brief, table-form output        */
/*                                                                          */
/*  double ILLutil_total_timer (ILLutil_timer *t, int printit)              */
/*    - Returns the cumulative time for this timer.                         */
/*    - if printit == 1, outputs the cumulative time.                       */
/*    - if printit == 2, outputs the cumulative time only if nonzero        */
/*    - if printit == 3,4, like 1,2, except brief, table-form output        */
/*                                                                          */
/****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "logging-private.h"

#include "util.h"

#ifdef HAVE_GETRUSAGE

#ifdef HAVE_SYS_RESOURCE_H
# include <sys/resource.h>
#endif

#define ZEIT_FCT "getrusage"

double ILLutil_zeit (
	void)
{
	struct rusage ru;
	double t;

	getrusage (RUSAGE_SELF, &ru);

	t = ((double) ru.ru_utime.tv_sec) +
		((double) ru.ru_utime.tv_usec) / 1000000.0;
	return t;
}
#else	/* HAVE_GETRUSAGE */

#ifdef HAVE_TIMES

#ifdef HAVE_SYS_PARAM_H
# include <sys/param.h>
#endif
#ifdef HAVE_SYS_TIMES_H
# include <sys/times.h>
#endif

#ifdef CLK_TCK
#define MACHINE_FREQ CLK_TCK
#else
#define MACHINE_FREQ HZ
#endif

#define ZEIT_FCT "times"

double ILLutil_zeit (
	void)
{
	struct tms now;

	times (&now);
	return ((double) now.tms_utime) / ((double) MACHINE_FREQ);
}
#else	/* HAVE_TIMES */

#ifdef HAVE_CLOCK

#ifndef CLOCKS_PER_SEC
#ifdef CLK_TCK
#define CLOCKS_PER_SEC CLK_TCK
#else
#define CLOCKS_PER_SEC 60
#endif
#endif

#define ZEIT_FCT "clock"

double ILLutil_zeit (
	void)
{
	return ((double) clock ()) / ((double) CLOCKS_PER_SEC);
}

#else	/* HAVE_CLOCK */

#define ZEIT_FCT "???"

double ILLutil_zeit (
	void)
{
	return 0.0;
}
#endif /* HAVE_CLOCK */
#endif /* HAVE_TIMES */
#endif /* HAVE_GETRUSAGE */

double ILLutil_real_zeit (
	void)
{
	return time (0);
}

void ILLutil_init_timer (
	ILLutil_timer * t,
	const char *name)
{
	t->szeit = -1.0;
	t->cum_zeit = 0.0;
	t->count = 0;
	if (name == (char *) NULL || name[0] == '\0')
	{
		strncpy (t->name, "ANONYMOUS", sizeof (t->name) - 1);
	}
	else
	{
		strncpy (t->name, name, sizeof (t->name) - 1);
	}
	t->name[sizeof (t->name) - 1] = '\0';
}

void ILLutil_start_timer (
	ILLutil_timer * t)
{
	if (t->szeit != -1.0)
	{
		QSlog("Warning: restarting running timer %s", t->name);
	}
	t->szeit = ILLutil_zeit ();
}

void ILLutil_suspend_timer (
	ILLutil_timer * t)
{
	if (t->szeit == -1.0)
	{
		QSlog("Warning: suspended non-running timer %s", t->name);
		return;
	}

	t->cum_zeit += ILLutil_zeit () - t->szeit;
	t->szeit = -1.0;
}

void ILLutil_resume_timer (
	ILLutil_timer * t)
{
	if (t->szeit != -1.0)
	{
		QSlog("Warning: resuming running timer %s", t->name);
		return;
	}
	t->szeit = ILLutil_zeit ();
}

static void ILL_print (
	ILLutil_timer * t,
	double z,
	int printit)
{
	if (printit == 1 || (printit == 2 && z > 0.0))
	{
		if (t->count > 1)
		{
			QSlog("Time for %s: %.2f seconds (%.2f total in %d calls).",
									t->name, z, t->cum_zeit, t->count);
		}
		else
		{
			QSlog("Time for %s: %.2f seconds.", t->name, z);
		}
	}
	else if (printit == 3 || (printit == 4 && z > 0.0))
	{
		QSlog("T %-34.34s %9.2f %9.2f %d (%s)",
								t->name, z, t->cum_zeit, t->count, ZEIT_FCT);
	}
}

double ILLutil_stop_timer (
	ILLutil_timer * t,
	int printit)
{
	double z;

	if (t->szeit == -1.0)
	{
		QSlog("Warning: stopping non-running timer %s", t->name);
		return 0.0;
	}
	z = ILLutil_zeit () - t->szeit;
	t->szeit = -1.0;
	t->cum_zeit += z;
	t->count++;
	ILL_print (t, z, printit);
	return z;
}

double ILLutil_total_timer (
	ILLutil_timer * t,
	int printit)
{
	double z = t->cum_zeit;

	if (t->szeit != -1.0)
		z += ILLutil_zeit () - t->szeit;
	ILL_print (t, z, printit);
	return z;
}
