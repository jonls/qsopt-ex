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

/* RCS_INFO = "$RCSfile: solver.c,v $ $Revision: 1.2 $ $Date: 2003/11/05 16:49:52 $"; */

#include "qs_config.h"
#include "exact.h"
#include "solver.h"
#include "iqsutil.h"
#include "util.h"
#include "lpdefs.h"							/* for PRIMAL_SIMPLEX */
#include "qstruct.h"
#include "qsopt.h"
#include "binary.h"
#include "editor.h"
#include "price.h"
#include "lib.h"								/* for ILLmip_binary_dfs */

static char *fname = 0;
static int lpfile = 0;
static int solvemip = 0;
static int interactive = 0;
static int usescaling = 1;
static int showversion = 0;
static int simplexalgo = PRIMAL_SIMPLEX;
static int pstrategy = QS_PRICE_PSTEEP;
static int dstrategy = QS_PRICE_DSTEEP;
static unsigned precision = 128;
static int printsol = 0;
static char *readbasis = 0;
static char *writebasis = 0;

static void usage (
	char *s),
  get_ftype (
	char *name,
	int *ftype);

#ifdef TEST_FEEDBACK
static int feedback (
	FILE * dest,
	const char *str);
#endif
static int parseargs (
	int ac,
	char **av);

#ifndef WIN32
QSLIB_INTERFACE int main ( int ac, char **av);
int main ( int ac, char **av)
{
	int rval;
	QSopt_ex_version();
	QSexactStart();
	rval = solver_main (ac, av);
	QSexactClear();
	return rval;
}
#endif

QSLIB_INTERFACE int solver_main (
	int ac,
	char **av)
{
	int rval = 0;
	QSdata *p = 0;
	ILLutil_timer timer_solve;
	ILLutil_timer timer_read;
	int ftype = 0;								/* 0 mps, 1 lp */
	EGlpNum_t *x = 0;
	EGlpNum_t val;
	double szeit;

	rval = parseargs (ac, av);
	if (rval)
		ILL_CLEANUP;

	QSset_precision (precision);
	EGlpNumInitVar (val);
	if (showversion)
	{
		char *buf = 0;

		buf = QSversion ();
		if (buf == 0)
		{
			ILL_CLEANUP;
		}
		else
		{
			printf ("%s\n", buf);
			QSfree ((void *) buf);
		}
	}

	if (lpfile)
	{
		ftype = 1;
	}
	else
	{
		get_ftype (fname, &ftype);
	}

	ILLutil_init_timer (&timer_read, "SOLVER_READ");
	ILLutil_start_timer (&timer_read);
	if (ftype == 1)
	{
		p = QSread_prob ((const char *) fname, "LP");
		if (p == 0)
		{
			fprintf (stderr, "Could not read lp file.\n");
			rval = 1;
			ILL_CLEANUP_IF (rval);
		}
	}
	else
	{
		p = QSread_prob ((const char *) fname, "MPS");
		if (p == 0)
		{
			fprintf (stderr, "Could not read mps file.\n");
			rval = 1;
			ILL_CLEANUP_IF (rval);
		}
	}

	if (readbasis)
	{
		rval = QSread_and_load_basis (p, (const char *) readbasis);
		ILL_CLEANUP_IF (rval);
	}
	ILLutil_stop_timer (&timer_read, 1);
	rval = QSset_param (p, QS_PARAM_SIMPLEX_DISPLAY, 1)
		|| QSset_param (p, QS_PARAM_PRIMAL_PRICING, pstrategy)
		|| QSset_param (p, QS_PARAM_DUAL_PRICING, dstrategy)
		|| QSset_param (p, QS_PARAM_SIMPLEX_SCALING, usescaling);
	ILL_CLEANUP_IF (rval);

	if (interactive)
	{
		ILLeditor_init ();
		ILLeditor (p);
		goto CLEANUP;
	}

	szeit = ILLutil_zeit();
	ILLutil_init_timer (&timer_solve, "SOLVER");
	ILLutil_start_timer (&timer_solve);

#ifdef TEST_FEEDBACK
	QSset_reporter (p, (void *) feedback, stdout);
#endif
	rval = ILLeditor_solve (p, simplexalgo);
	ILLutil_stop_timer (&timer_solve, 1);
	ILL_CLEANUP_IF (rval);

	printf ("ZZ %s %.2f\n", p->lp->O->probname, ILLutil_zeit () - szeit);
 	fflush (stdout);

	if (printsol)
		x = EGlpNumAllocArray (p->lp->O->nstruct);

	if (solvemip)
	{
		rval = ILLmip_bfs (p->lp, &val, x, &(p->itcnt));
		ILL_CLEANUP_IF (rval);
		printf ("MIP Objective Value: %.6f\n", EGlpNumToLf (val));
		fflush (stdout);
		if (printsol && EGlpNumIsNeqq (val, ILL_MAXDOUBLE) &&
				EGlpNumIsNeqq (val, ILL_MINDOUBLE))
		{
			EGioFile_t*out = EGioOpenFILE(stdout);
			rval = ILLlib_print_x (out, p->lp, 0, x, 1);
			EGioClose(out);
			ILL_CLEANUP_IF (rval);
		}
	}
	else
	{
		if (writebasis)
		{
			rval = QSwrite_basis (p, 0, writebasis);
			ILL_CLEANUP_IF (rval);
		}
		if (printsol)
		{
			EGioFile_t*out = EGioOpenFILE(stdout);
			rval = ILLlib_print_x (out, p->lp, 0, 0, 1);
			EGioClose(out);
			ILL_CLEANUP_IF (rval);
		}
	}

CLEANUP:
	EGlpNumFreeArray (x);
	QSfree_prob (p);
	EGlpNumClearVar (val);
	return rval;									/* main return */
}

static void usage (
	char *s)
{
	char *buf = 0;

	buf = QSversion ();
	if (buf)
	{
		fprintf (stderr, "%s\n", buf);
		QSfree ((void *) buf);
	}

	fprintf (stderr, "Usage: %s [- below -] prob_file\n", s);
	fprintf (stderr, "   -b f  write basis to file f\n");
	fprintf (stderr, "   -B f  read initial basis from file f\n");
#if 0
	fprintf (stderr, "   -I    solve the MIP using BestBound\n");
	fprintf (stderr, "   -E    edit problem after solving initial version\n");
#endif
	fprintf (stderr, "   -L    input file is in lp format (default: mps)\n");
	fprintf (stderr, "   -O    print the final solution\n");
	fprintf (stderr, "   -p #  run primal simplex with pricing rule #\n");
	fprintf (stderr,
					 "         (%d-Dantzig, %d-Devex, %d-Steep (default), %d-Partial\n",
					 QS_PRICE_PDANTZIG, QS_PRICE_PDEVEX, QS_PRICE_PSTEEP,
					 QS_PRICE_PMULTPARTIAL);
	fprintf (stderr, "   -d #  run dual simplex with pricing rule #\n");
	fprintf (stderr, "         (%d-Dantzig, %d-Steep, %d-Partial, %d-Devex)\n",
					 QS_PRICE_DDANTZIG, QS_PRICE_DSTEEP, QS_PRICE_DMULTPARTIAL,
					 QS_PRICE_DDEVEX);
	fprintf (stderr, "   -S    do NOT scale the initial LP\n");
	fprintf (stderr, "   -v    print QSopt version number\n");
}

static int parseargs (
	int ac,
	char **av)
{
	int c;
	int boptind = 1;
	char *boptarg = 0;

	while ((c =
					ILLutil_bix_getopt (ac, av, "b:B:d:p:P:IELOSv", &boptind,
															&boptarg)) != EOF)
		switch (c)
		{
		case 'b':
			writebasis = boptarg;
			break;
		case 'B':
			readbasis = boptarg;
			break;
		case 'P':
			precision = atoi (boptarg);
			break;
		case 'd':
			simplexalgo = DUAL_SIMPLEX;
			dstrategy = atoi (boptarg);
			break;
#if 0
		case 'E':
			interactive = 1;
			break;
		case 'I':
			solvemip = 1;
			break;
#endif
		case 'L':
			lpfile = 1;
			break;
		case 'O':
			printsol = 1;
			break;
		case 'p':
			simplexalgo = PRIMAL_SIMPLEX;
			pstrategy = atoi (boptarg);
			break;
		case 'S':
			usescaling = 0;
			break;
		case 'v':
			showversion = 1;
			break;
		case '?':
		default:
			usage (av[0]);
			return 1;
		}
	if (boptind != (ac - 1))
	{
		usage (av[0]);
		return 1;
	}

	fname = av[boptind++];
	fprintf (stderr, "Reading problem from %s\n", fname);
	return 0;
}

static void get_ftype (
	char *name,
	int *ftype)
{
	char *q;

	q = strrchr (name, '.');
	if (q)
	{
		q++;
		if (!strcmp (q, "lp") || !strcmp (q, "LP"))
		{
			*ftype = 1;
		}
		else
		{
			*ftype = 0;
		}
	}
}

#ifdef TEST_FEEDBACK
static int feedback (
	FILE * dest,
	const char *str)
{
	if (str != NULL)
	{
		int rc = fprintf ((FILE *) dest, "FEEDBACK: %s", str);

		fflush (dest);
		return rc;
	}
	return 0;
}
#endif
