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

/* RCS_INFO = "$Id: reporter.c,v 1.1 2003/11/05 16:49:52 meven Exp $"; */

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <stdio.h>

#include "reporter.h"
#include "except.h"


int ILL_fprintf (
	void *dest,
	const char *s)
{
	if (s != NULL)
		return fputs(s, (FILE *)dest);
	return 0;
}

void ILLstring_reporter_init (
	qsstring_reporter * reporter,
	qsreport_string_fct fct,
	void *dest)
{
	int rval = 0;

	ILL_FAILfalse (reporter != NULL, "Must get non NULL reporter");
	if (reporter != NULL)
	{
		reporter->report_fct = fct;
		reporter->dest = dest;
	}
CLEANUP:
	return;
}

void ILLstring_reporter_copy (
	qsstring_reporter * dest,
	qsstring_reporter * src)
{
	*dest = *src;
}


#ifdef REPORTER_MAIN
static int string_reporter_main (
	int ac,
	char **av)
{
	int i = 0;
	qsstring_reporter reporter;

	ILLstring_reporter_init (&reporter, ILL_fprintf, stdout);
	for (i = 0; i < ac; i++)
	{
		(void) ILLstring_report (av[i], &reporter);
		(void) ILLstring_report ("\n", &reporter);
	}
	(void) ILLstring_report (NULL, &reporter);

	return 0;
}

int main (
	int ac,
	char **av)
{
	return string_reporter_main (ac, av);
}
#endif
