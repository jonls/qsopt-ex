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

/* $RCSfile: reporter.h,v $ $Revision: 1.2 $ $Date: 2003/11/05 16:57:39 $ */
#ifndef REPORTER_FILE
#define REPORTER_FILE

#include "eg_io.h"

typedef int (
	*qsreport_string_fct) (
	void *dest,
	const char *s);

typedef struct qsstring_reporter
{
	qsreport_string_fct report_fct;
	void *dest;
}
qsstring_reporter;

extern int ILL_fprintf (
	void *dest,
	const char *s);

void ILLstring_reporter_init (
	qsstring_reporter * reporter,
	qsreport_string_fct fct,
	void *dest);

void ILLstring_reporter_copy (
	qsstring_reporter * dest,
	qsstring_reporter * src);

#define ILLstring_report(s, reporter)  \
	((reporter)->report_fct((reporter)->dest, s) < 0)
									 /* used from with ILL fct to report progress */

/* REPORTER_FILE */
#endif
