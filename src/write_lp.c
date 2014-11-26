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

/* RCS_INFO = "$RCSfile: wr_lp.c,v $ $Revision: 1.2 $ $Date: 2003/11/05 16:49:52 $"; */

/****************************************************************************/
/*                                                                          */
/*               Routines to support writing of LP files                    */
/*                                                                          */
/****************************************************************************/

/* 
 * -) anything after '\' is comment 
 * -) variables consist of a-z A-Z 0-9!"#$%(),;.?@_`'{}|~ 
 *    don't start with a digit or '.'
 */

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <stdlib.h>
#include <string.h>

#include "eg_lpnum.h"
#include "eg_io.h"

#include "iqsutil.h"
#include "lpdefs.h"
#include "write_lp.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

void ILLwrite_lp_state_init (
	ILLwrite_lp_state * line,
	const char *str)
{
	line->total = 0;
	line->p = line->buf;
	*line->p = '\0';
	if (str != NULL)
	{
		ILLwrite_lp_state_append (line, str);
	}
}

void ILLwrite_lp_state_append (
	ILLwrite_lp_state * line,
	const char *str)
{
	int len, rval = 0;

	ILL_FAILfalse (str, "Must have non NULL string");
	sprintf (line->p, "%s", str);
	len = strlen (line->p);
	line->total += len;
	line->p += len;
CLEANUP:
	return;
}

void ILLwrite_lp_state_append_coef (
	ILLwrite_lp_state * line,
	EGlpNum_t v,
	int cnt)
{
	EGlpNum_t ntmp;
	int len = 0;

	EGlpNumInitVar (ntmp);
	EGlpNumCopy (ntmp, v);
	if (EGlpNumIsLessZero (ntmp))
	{
		sprintf (line->p, " - ");
		len = 3;
		EGlpNumSign (ntmp);
	}
	else
	{
		if (cnt > 0)
		{
			sprintf (line->p, " + ");
			len = 3;
		}
		else
		{
			sprintf (line->p, " ");
			len = 1;
		}
	}
	line->p += len;
	line->total += len;
	if (EGlpNumIsNeqq (ntmp, oneLpNum))
	{
		ILLwrite_lp_state_append_number (line, ntmp);
	}
	EGlpNumClearVar (ntmp);
}

/* so that diff will not stumble over too many number format differences 
 * between c and java generated lp files we make here sure that doubles 
 * which are printed without a ".xxx" part get ".0" from us
 */
static void append_number (
	ILLwrite_lp_state * line,
	EGlpNum_t v);
void ILLwrite_lp_state_append_number (
	ILLwrite_lp_state * line,
	EGlpNum_t v)
{
	/* write a blank after 'inf' in case it is used as a coefficient and 
	 * a variable follows */
	if (EGlpNumIsEqqual (v, ILL_MAXDOUBLE))
	{
		ILLwrite_lp_state_append (line, "inf ");
	}
	else
	{
		if (EGlpNumIsEqqual (v, ILL_MINDOUBLE))
		{
			ILLwrite_lp_state_append (line, "-inf ");
		}
		else
			append_number (line, v);
	}
}

static void append_number (
	ILLwrite_lp_state * line,
	EGlpNum_t v)
{
	int len = 0;
	char *numstr = EGlpNumGetStr (v);

	sprintf (line->p, "%s%n", numstr, &len);
	EGfree (numstr);
	line->p += len;
	line->total += len;
}

#if 0
#define D_SCALE (1e9)
static void append_number (
	ILLwrite_lp_state * line,
	double x)
{
	/* Better code for writing rational problems */
	int i, k;
	int got = 0;
	double w, t;
	int nnum = 0;
	int nden = 0;

	if (x != 0.0)
	{
		if (x < 0.0)
			w = -x;
		else
			w = x;

		for (i = -9, t = 0.000000001; i <= 7; i++, t *= 10.0)
		{
			if (w >= t && w <= t * 10)
			{
				got = 1;
				break;
			}
		}
		if (got == 0)
		{
			fprintf (stderr, "Out-of-range number: %f\n", x);
			exit (1);
		}
	}

	if (x < 0.0)
	{
		sprintf (line->p, "-");
		line->p++;
		line->total++;
		x = -x;
	}

	while (x >= 10.0 * D_SCALE)
	{
		x /= 10.0;
		nnum++;
	}

	/* (x != (double) (int) x) is a hack to let small integers print nicely */

	while (x < D_SCALE && x != (double) (int) x)
	{
		x *= 10.0;
		nden++;
	}

	sprintf (line->p, "%.0f%n", x, &k);
	line->p += k;
	line->total += k;

	for (i = 0; i < nnum; i++)
	{
		sprintf (line->p, "0");
		line->p++;
		line->total++;
	}

	if (nden)
	{
		sprintf (line->p, "/1%n", &k);
		line->p += k;
		line->total += k;
		for (i = 0; i < nden; i++)
		{
			sprintf (line->p, "0");
			line->p++;
			line->total++;
		}
	}
}
#endif

void ILLwrite_lp_state_save_start (
	ILLwrite_lp_state * line)
{
	line->startlen = line->total;
}

void ILLwrite_lp_state_start (
	ILLwrite_lp_state * line)
{
	int j;

	for (j = 0; j < line->startlen; j++)
	{
		line->buf[j] = ' ';
	}
	line->buf[j] = '\0';
	line->p = line->buf + j;
	line->total = j;
}
