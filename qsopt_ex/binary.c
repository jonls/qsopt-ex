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

/* RCS_INFO = "$RCSfile: binary.c,v $ $Revision: 1.2 $ $Date: 2003/11/05 16:49:52 $"; */
static int TRACE = 0;

/****************************************************************************/
/*                                                                          */
/*                     Simple MIP Code to test LP Solver                    */
/*                                                                          */
/*  EXPORTED FUNCTIONS                                                      */
/*                                                                          */
/*    int ILLmip_bfs (lpinfo *lp, double *val, double *x)                   */
/*                                                                          */
/*  NOTES                                                                   */
/*                                                                          */
/*                                                                          */
/****************************************************************************/

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "qs_config.h"

#include "eg_lpnum.h"
#include "eg_io.h"

#include "priority.h"
#include "sortrus.h"
#include "iqsutil.h"
#include "lpdata.h"
#include "lpdefs.h"
#include "simplex.h"
#include "binary.h"
#include "price.h"
#include "lib.h"
#include "qstruct.h"
#include "qsopt.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

/*#define  ILL_INTTOL (0.000001)*/
#define ILL_INTTOL PFEAS_TOLER

#define  STRONG_PIVOTS     (50)
#define  STRONG_CANDIDATES (10)

#define ILL_BRANCH_STRONG_WEIGHT (10)
#define ILL_BRANCH_STRONG_VAL(v0,v1)                                 \
    (((v0) < (v1) ? (ILL_BRANCH_STRONG_WEIGHT * (v0) + (v1))         \
                  : (ILL_BRANCH_STRONG_WEIGHT * (v1) + (v0)))        \
                    / (ILL_BRANCH_STRONG_WEIGHT + 1.0))

#define ILL_BRANCH_PENALTY_WEIGHT (2)
#define ILL_BRANCH_PENALTY_VAL(v0,v1,f)                              \
    (((v0)*(f) < (v1)*(1.0-(f)) ?                                    \
        (ILL_BRANCH_PENALTY_WEIGHT * (v0)*(f) + (v1)*(1.0-(f)))    \
      : (ILL_BRANCH_PENALTY_WEIGHT * (v1)*(1.0-(f)) + (v0)*(f)))    \
                    / (ILL_BRANCH_PENALTY_WEIGHT + 1.0))



#define FIRSTBRANCH  1
#define MIDDLEBRANCH 2
#define STRONGBRANCH 3
#define PENALTYBRANCH 4


typedef struct bbnode
{
	struct bbnode *next;
	struct bbnode *prev;
	int id;
	int depth;
	int handle;
	EGlpNum_t bound;
	char *cstat;
	char *rstat;
	EGlpNum_t *rownorms;
	int rownorms_size;
	int bound_cnt;
	int *bound_indx;
	char *lu;
	EGlpNum_t *bounds;
	int bounds_size;
}
bbnode;

typedef struct mipinfo
{
	int branching_rule;
	int watch;
	int depth;
	int totalnodes;
	int activenodes;
	int totalpivots;
	int lastpivots;
	int objsense;
	EGlpNum_t objectivebound;
	EGlpNum_t value;
	EGlpNum_t *downpen;
	EGlpNum_t *uppen;
	EGlpNum_t *x;
	EGlpNum_t *bestx;
	EGlpNum_t *orig_lower;
	EGlpNum_t *orig_upper;
	EGlpNum_t *lower;
	EGlpNum_t *upper;
	int nstruct;									/* size of all EGlpNum_t arrays */
	lpinfo *lp;
	price_info *pinf;
	bbnode head_bbnode;
	ILLpriority *que;
	ILLptrworld ptrworld;
}
mipinfo;


ILL_PTRWORLD_ROUTINES (bbnode, bbnodealloc, bbnode_bulkalloc, bbnodefree)
ILL_PTRWORLD_LISTFREE_ROUTINE (bbnode, bbnode_listfree, bbnodefree)
ILL_PTRWORLD_LEAKS_ROUTINE (bbnode, bbnode_check_leaks, depth, int)
static void cleanup_mip ( mipinfo * minf), 
		choose_initial_price ( price_info * pinf), 
		best_bbnode ( mipinfo * minf, bbnode ** best),
		put_bbnode ( mipinfo * minf, bbnode * b),
		remove_bbnode ( bbnode * b),
		find_first_branch ( lpinfo * lp, EGlpNum_t * x, int *bvar),
		find_middle_branch ( lpinfo * lp, EGlpNum_t * x, int *bvar),
		check_integral ( lpinfo * lp, EGlpNum_t * x, int *yesno),
		copy_x ( int nstruct, EGlpNum_t * from_x, EGlpNum_t * to_x),
		init_mipinfo ( mipinfo * minf),
		free_mipinfo ( mipinfo * minf),
		init_bbnode ( bbnode * b),
		free_bbnode ( bbnode * b);

static int startup_mip ( mipinfo * minf, lpinfo * lp, price_info * pinf,
			EGlpNum_t * lpval, itcnt_t*itcnt),
		run_bfs ( mipinfo * minf, itcnt_t*itcnt),
		process_bfs_bbnode ( mipinfo * minf, bbnode * b, itcnt_t*itcnt),
		child_work ( mipinfo * minf, bbnode * active, int bvar, int bdir,
			EGlpNum_t * cval, int *cp, itcnt_t*itcnt),
		fix_variables ( lpinfo * lp, EGlpNum_t * bestval, bbnode * b,
			EGlpNum_t * wupper, EGlpNum_t * wlower, int *hit),
		find_branch ( mipinfo * minf, EGlpNum_t * x, EGlpNum_t * lpval,
			int *bvar, itcnt_t*itcnt),
		find_penalty_branch ( lpinfo * lp, price_info * pinf, EGlpNum_t * x,
			EGlpNum_t * downpen, EGlpNum_t * uppen, EGlpNum_t * lpval, int *bvar,
			itcnt_t*itcnt),
		find_strong_branch ( lpinfo * lp, price_info * pinf, EGlpNum_t * x,
			int *bvar, itcnt_t*itcnt),
		plunge ( mipinfo * minf, itcnt_t*itcnt),
		plunge_work ( mipinfo * minf, int depth, itcnt_t*itcnt),
		round_variables ( mipinfo * minf, int *count, EGlpNum_t * tol);

static void choose_initial_price ( price_info * pinf)
{
	pinf->pI_price = QS_PRICE_PSTEEP;
	pinf->pII_price = QS_PRICE_PSTEEP;
	pinf->dI_price = QS_PRICE_DSTEEP;
	pinf->dII_price = QS_PRICE_DSTEEP;
}

int ILLmip_bfs (
	lpinfo * lp,
	EGlpNum_t * val,
	EGlpNum_t * x,
	itcnt_t*itcnt)
{
	int tval, rval = 0;
	price_info pinf;
	mipinfo minf;
	bbnode *b;
	EGlpNum_t lpval;
	double szeit = ILLutil_zeit ();

	EGlpNumInitVar (lpval);
	EGlpNumInitVar (pinf.htrigger);

	ILLprice_init_pricing_info (&pinf);
	init_mipinfo (&minf);

	if (!lp)
	{
		fprintf (stderr, "ILLmip_bfs called without an LP\n");
		rval = 1;
		goto CLEANUP;
	}

	rval = startup_mip (&minf, lp, &pinf, &lpval, itcnt);
	ILL_CLEANUP_IF (rval);

	ILL_SAFE_MALLOC (minf.que, 1, ILLpriority);
	rval = ILLutil_priority_init (minf.que, lp->O->nstruct + 1);
	ILL_CLEANUP_IF (rval);

	b = bbnodealloc (&minf.ptrworld);
	init_bbnode (b);
	b->depth = 0;
	b->id = minf.totalnodes++;
	EGlpNumCopy (b->bound, lpval);
	ILL_SAFE_MALLOC (b->cstat, lp->O->nstruct, char);
	ILL_SAFE_MALLOC (b->rstat, lp->nrows, char);

	rval = ILLlib_getbasis (lp, b->cstat, b->rstat);
	ILL_CLEANUP_IF (rval);

	if (pinf.dII_price == QS_PRICE_DSTEEP)
	{
		b->rownorms = EGlpNumAllocArray (lp->nrows);
		tval = ILLlib_getrownorms (lp, &pinf, b->rownorms);
		if (tval)
		{
			printf ("Row norms not available\n");
			fflush (stdout);
			EGlpNumFreeArray (b->rownorms);
		}
	}

	rval = ILLutil_priority_insert (minf.que, (void *) b, &lpval, &(b->handle));
	ILL_CLEANUP_IF (rval);

	b->prev = &(minf.head_bbnode);
	b->next = 0;
	minf.head_bbnode.next = b;
	minf.activenodes++;

	minf.branching_rule = PENALTYBRANCH;

	rval = run_bfs (&minf, itcnt);
	ILL_CLEANUP_IF (rval);

	printf ("Total Number of Nodes: %d\n", minf.totalnodes);
	printf ("Total Number of Pivots: %d\n", minf.totalpivots);
	printf ("BFS MIP Runing Time: %.2f seconds\n", ILLutil_zeit () - szeit);
	fflush (stdout);

	EGlpNumCopy (*val, minf.value);
	if (minf.objsense == ILL_MAX)
		EGlpNumSign (*val);

	if (x && EGlpNumIsNeqq (minf.value, ILL_MAXDOUBLE))
	{
		copy_x (lp->O->nstruct, minf.bestx, x);
	}

CLEANUP:

	if (minf.que)
	{
		ILLutil_priority_free (minf.que);
		ILL_IFFREE (minf.que, ILLpriority);
	}
	cleanup_mip (&minf);
	free_mipinfo (&minf);
	ILLprice_free_pricing_info (&pinf);
	EGlpNumClearVar (lpval);
	EGlpNumClearVar (pinf.htrigger);
	ILL_RETURN (rval, "ILLmip_bfs");
}

static int startup_mip (
	mipinfo * minf,
	lpinfo * lp,
	price_info * pinf,
	EGlpNum_t * lpval,
	itcnt_t*itcnt)
{
	int rval = 0;
	int i, col, status, intcount = 0;
	EGlpNum_t val;
	ILLlpdata *qlp;

	EGlpNumInitVar (val);

	choose_initial_price (pinf);

	qlp = lp->O;

	rval = ILLlib_optimize (lp, 0, pinf, DUAL_SIMPLEX, &status, 0, itcnt);
	ILL_CLEANUP_IF (rval);

	minf->totalpivots += ILLlib_iter (lp);

	rval = ILLlib_objval (lp, 0, &val);
	ILL_CLEANUP_IF (rval);

	printf ("LP Value: %.6f\n", EGlpNumToLf (val));
	fflush (stdout);
	if (lpval)
		EGlpNumCopy (*lpval, val);

	if (qlp->intmarker)
	{
		for (i = 0; i < qlp->nstruct; i++)
		{
			if (qlp->intmarker[i])
			{
				col = qlp->structmap[i];
				intcount++;
				if (EGlpNumIsEqqual (qlp->lower[col], ILL_MINDOUBLE)
						|| EGlpNumIsEqqual (qlp->upper[col], ILL_MAXDOUBLE))
				{
					printf ("Instance has unbounded integer variable\n");
					fflush (stdout);
					rval = 1;
					goto CLEANUP;
				}
			}
		}
	}

	if (intcount == 0)
	{
		printf ("No integer variables\n");
		fflush (stdout);
		rval = 1;
		goto CLEANUP;
	}
	else
	{
		printf ("%d integer variables\n", intcount);
		fflush (stdout);
	}

	if (qlp->sinfo)
	{															/* Free the presolve LP and work with orginal */
		ILLlp_sinfo_free (qlp->sinfo);
		ILL_IFFREE (qlp->sinfo, ILLlp_sinfo);
	}


	minf->lp = lp;
	minf->pinf = pinf;
	minf->objsense = qlp->objsense;
	if (qlp->objsense == ILL_MAX)
	{															/* MIP codes work with min */
		for (i = 0; i < lp->ncols; i++)
		{
			EGlpNumCopyNeg (qlp->obj[i], qlp->obj[i]);
		}
		qlp->objsense = ILL_MIN;
	}

	minf->x = EGlpNumAllocArray (qlp->nstruct);
	minf->bestx = EGlpNumAllocArray (qlp->nstruct);
	minf->lower = EGlpNumAllocArray (qlp->nstruct);
	minf->upper = EGlpNumAllocArray (qlp->nstruct);
	minf->orig_lower = EGlpNumAllocArray (qlp->nstruct);
	minf->orig_upper = EGlpNumAllocArray (qlp->nstruct);
	minf->downpen = EGlpNumAllocArray (qlp->nstruct);
	minf->uppen = EGlpNumAllocArray (qlp->nstruct);
	minf->nstruct = qlp->nstruct;

	rval = ILLlib_get_x (lp, 0, minf->x);
	ILL_CLEANUP_IF (rval);

	for (i = 0; i < qlp->nstruct; i++)
	{
		EGlpNumCopy (minf->lower[i], qlp->lower[i]);
		EGlpNumCopy (minf->upper[i], qlp->upper[i]);
		EGlpNumCopy (minf->orig_lower[i], qlp->lower[i]);
		EGlpNumCopy (minf->orig_upper[i], qlp->upper[i]);
		EGlpNumOne (minf->downpen[i]);
		EGlpNumOne (minf->uppen[i]);
		EGlpNumSign (minf->downpen[i]);
		EGlpNumSign (minf->uppen[i]);
	}


CLEANUP:

	EGlpNumClearVar (val);
	ILL_RETURN (rval, "startup_mip");
}

static void cleanup_mip (
	mipinfo * minf)
{
	int i;
	ILLlpdata *qslp;

	if (minf && minf->lp)
	{
		qslp = minf->lp->O;
		if (minf->objsense == ILL_MAX)
		{
			for (i = 0; i < minf->lp->ncols; i++)
			{
				EGlpNumSign (qslp->obj[i]);
			}
			qslp->objsense = ILL_MIN;
		}
	}
}

static int run_bfs (
	mipinfo * minf,
	itcnt_t*itcnt)
{
	int rval = 0;
	bbnode *b;

	while (minf->head_bbnode.next)
	{
		best_bbnode (minf, &b);
		rval = process_bfs_bbnode (minf, b, itcnt);
		ILL_CLEANUP_IF (rval);
		remove_bbnode (b);
		free_bbnode (b);
		bbnodefree (&minf->ptrworld, b);
		minf->activenodes--;
	}

CLEANUP:

	ILL_RETURN (rval, "run_bfs");
}

static int process_bfs_bbnode (
	mipinfo * minf,
	bbnode * active,
	itcnt_t*itcnt)
{
	lpinfo *lp = minf->lp;
	ILLlp_basis B;
	int status, bvar = 0;
	int i, j, hit, dnp = 0, upp = 0;
	int nstruct = lp->O->nstruct;
	EGlpNum_t t, lpval, dnval, upval;
	EGlpNum_t *wupper = 0;
	EGlpNum_t *wlower = 0;
	int rval = 0;

	EGlpNumInitVar (t);
	EGlpNumInitVar (lpval);
	EGlpNumInitVar (dnval);
	EGlpNumInitVar (upval);

	ILLlp_basis_init (&B);

	if (minf->watch > 1)
	{
		printf ("Node %4d: %.3f", active->id, EGlpNumToLf (active->bound));
		if (EGlpNumIsNeqq (minf->value, ILL_MAXDOUBLE))
			printf (" %.3f", EGlpNumToLf (minf->value));
		else
			printf ("  None");
		printf (", Active %d ", minf->activenodes);
		fflush (stdout);
	}
	else if (minf->watch == 1)
	{
		if (minf->lastpivots > 1000)
		{
			minf->lastpivots = 0;
			printf ("Pivots %d, Active Nodes %d, Bound %.3f, Soln ",
							minf->totalpivots, minf->activenodes,
							EGlpNumToLf (active->bound));
			if (!EGlpNumIsLess (minf->value, ILL_MAXDOUBLE))
				printf ("%.3f", EGlpNumToLf (minf->value));
			else
				printf ("None\n");
		}
	}

	if (EGlpNumIsLeq (minf->objectivebound, active->bound))
	{
		if (minf->watch > 1)
		{
			printf ("  Node can be purged\n");
			fflush (stdout);
		}
		goto CLEANUP;
	}

	/*  Set the LP bounds for the node. */

	wlower = EGlpNumAllocArray (nstruct);
	wupper = EGlpNumAllocArray (nstruct);

	for (i = 0; i < nstruct; i++)
	{
		EGlpNumCopy (wlower[i], minf->orig_lower[i]);
		EGlpNumCopy (wupper[i], minf->orig_upper[i]);
	}
	for (i = 0; i < active->bound_cnt; i++)
	{
		j = active->bound_indx[i];
		if (active->lu[i] == 'L')
			EGlpNumCopy (wlower[j], active->bounds[i]);
		else
			EGlpNumCopy (wupper[j], active->bounds[i]);
	}

	if (active->bound_cnt > 0)
	{
		rval = ILLlib_chgbnds (lp, active->bound_cnt, active->bound_indx,
													 active->lu, active->bounds);
		ILL_CLEANUP_IF (rval);
	}

	/*  Solve the LP. */

	rval = ILLlib_loadbasis (&B, nstruct, lp->nrows, active->cstat,
													 active->rstat);
	ILL_CLEANUP_IF (rval);
	if (active->rownorms)
	{
		B.rownorms = EGlpNumAllocArray (lp->nrows);
		for (i = 0; i < lp->nrows; i++)
		{
			EGlpNumCopy (B.rownorms[i], active->rownorms[i]);
		}
	}

	rval = ILLlib_optimize (lp, &B, minf->pinf, DUAL_SIMPLEX, &status, 0, itcnt);
	ILL_CLEANUP_IF (rval);

	minf->totalpivots += ILLlib_iter (lp);
	minf->lastpivots += ILLlib_iter (lp);

	if (status == QS_LP_UNSOLVED)
	{
		printf ("Simplex did not solve the LP\n");
		fflush (stdout);
		rval = 1;
		ILL_CLEANUP;
	}

	if (status == QS_LP_INFEASIBLE)
	{
		printf ("  Infeasible LP, should have been purged earlier\n");
		fflush (stdout);
		rval = 1;
		ILL_CLEANUP;
	}

	if (active->depth < 0)
	{
		for (i = 0; i < nstruct; i++)
		{
			EGlpNumCopy (minf->lower[i], wlower[i]);
			EGlpNumCopy (minf->upper[i], wupper[i]);
		}
		rval = plunge (minf, itcnt);
		ILL_CLEANUP_IF (rval);
	}

	/*  Fix variables. */

	if (EGlpNumIsLess (minf->value, ILL_MAXDOUBLE))
	{
		rval = fix_variables (lp, &(minf->value), active, wupper, wlower, &hit);
		ILL_CLEANUP_IF (rval);

		if (hit)
		{
			rval = ILLlib_optimize (lp, &B, minf->pinf, DUAL_SIMPLEX, &status, 0, itcnt);
			ILL_CLEANUP_IF (rval);

			minf->totalpivots += ILLlib_iter (lp);
			minf->lastpivots += ILLlib_iter (lp);

			if (status == QS_LP_UNSOLVED)
			{
				printf ("Simplex did not solve the LP\n");
				fflush (stdout);
				rval = 1;
				ILL_CLEANUP;
			}

			if (status == QS_LP_INFEASIBLE)
			{
				printf ("  Infeasible LP after fixing\n");
				fflush (stdout);
				rval = 1;
				ILL_CLEANUP;
			}
		}
	}


	/*  Branch. */

	rval = ILLlib_get_x (lp, 0, minf->x);
	ILL_CLEANUP_IF (rval);

	rval = ILLlib_objval (lp, 0, &lpval);
	ILL_CLEANUP_IF (rval);

	rval = find_branch (minf, minf->x, &lpval, &bvar, itcnt);
	ILL_CLEANUP_IF (rval);

	if (bvar == -1)
	{
		printf ("Found integral solution: %f\n", EGlpNumToLf (lpval));
		if (EGlpNumIsLess (lpval, minf->value))
		{
			EGlpNumCopy (minf->value, lpval);
			EGlpNumCopy (minf->objectivebound, lpval);
			EGlpNumSubTo (minf->objectivebound, ILL_INTTOL);
			copy_x (nstruct, minf->x, minf->bestx);
		}
	}
	else
	{
		/* Create down child */

		rval = child_work (minf, active, bvar, 'D', &dnval, &dnp, itcnt);
		ILL_CLEANUP_IF (rval);

		/* Restore parent basis */

		rval = ILLlib_loadbasis (&B, nstruct, lp->nrows, active->cstat,
														 active->rstat);
		ILL_CLEANUP_IF (rval);
		if (active->rownorms)
		{
			B.rownorms = EGlpNumAllocArray (lp->nrows);
			for (i = 0; i < lp->nrows; i++)
			{
				EGlpNumCopy (B.rownorms[i], active->rownorms[i]);
			}
		}

		/* Create up child */

		rval = child_work (minf, active, bvar, 'U', &upval, &upp, itcnt);
		ILL_CLEANUP_IF (rval);

		if (minf->watch > 1)
		{
			if (EGlpNumIsEqqual (dnval, ILL_MAXDOUBLE))
			{
				printf ("DN->XXX");
			}
			else
			{
				printf ("DN->%.3f%c", EGlpNumToLf (dnval), dnp ? 'X' : ' ');
			}
			if (EGlpNumIsEqqual (upval, ILL_MAXDOUBLE))
			{
				printf ("UP->XXX\n");
			}
			else
			{
				printf ("UP->%.3f%c\n", EGlpNumToLf (upval), upp ? 'X' : ' ');
			}
			fflush (stdout);
		}
	}

	/* Set the LP bounds back to original values */

	for (i = 0; i < active->bound_cnt; i++)
	{
		if (active->lu[i] == 'L')
			EGlpNumCopy (t, minf->orig_lower[active->bound_indx[i]]);
		else
			EGlpNumCopy (t, minf->orig_upper[active->bound_indx[i]]);

		rval = ILLlib_chgbnd (lp, active->bound_indx[i], active->lu[i], t);
		ILL_CLEANUP_IF (rval);
	}

CLEANUP:

	EGlpNumFreeArray (wlower);
	EGlpNumFreeArray (wupper);
	ILLlp_basis_free (&B);
	EGlpNumClearVar (t);
	EGlpNumClearVar (lpval);
	EGlpNumClearVar (dnval);
	EGlpNumClearVar (upval);
	ILL_RETURN (rval, "process_bfs_bbnode");
}

static int child_work (
	mipinfo * minf,
	bbnode * active,
	int bvar,
	int bdir,
	EGlpNum_t * cval,
	int *cp,
	itcnt_t*itcnt)
{
	int tval, rval = 0;
	int i, status, intsol;
	EGlpNum_t t, oldt, lpval;
	EGlpNum_t *xi = &(minf->x[bvar]);
	lpinfo *lp = minf->lp;
	bbnode *b;

	EGlpNumInitVar (t);
	EGlpNumInitVar (lpval);
	EGlpNumInitVar (oldt);

	*cp = 0;

	if (bdir == 'D')
	{
		rval = ILLlib_getbnd (lp, bvar, 'U', &oldt);
		ILL_CLEANUP_IF (rval);
		EGlpNumFloor (t, *xi);
		rval = ILLlib_chgbnd (lp, bvar, 'U', t);
		ILL_CLEANUP_IF (rval);
	}
	else
	{
		rval = ILLlib_getbnd (lp, bvar, 'L', &oldt);
		ILL_CLEANUP_IF (rval);
		EGlpNumCeil (t, *xi);
		rval = ILLlib_chgbnd (lp, bvar, 'L', t);
		ILL_CLEANUP_IF (rval);
	}

	rval = ILLlib_optimize (lp, 0, minf->pinf, DUAL_SIMPLEX, &status, 0, itcnt);
	ILL_CLEANUP_IF (rval);

	minf->totalpivots += ILLlib_iter (lp);
	minf->lastpivots += ILLlib_iter (lp);

	if (status == QS_LP_UNSOLVED)
	{
		printf ("Simplex did not solve Child LP\n");
		fflush (stdout);
		rval = 1;
		ILL_CLEANUP;
	}

	if (status == QS_LP_INFEASIBLE)
	{
		EGlpNumCopy (*cval, ILL_MAXDOUBLE);
		*cp = 1;
	}
	else
	{
		rval = ILLlib_objval (lp, 0, &lpval);
		ILL_CLEANUP_IF (rval);
		EGlpNumCopy (*cval, lpval);

		/* What about the x vector?  Bico - 020531 */

		check_integral (lp, minf->x, &intsol);
		if (intsol)
		{
			if (EGlpNumIsLess (lpval, minf->value))
			{
				printf ("Found integral solution: %f\n", EGlpNumToLf (lpval));
				EGlpNumCopy (minf->value, lpval);
				EGlpNumCopy (minf->objectivebound, lpval);
				EGlpNumSubTo (minf->objectivebound, ILL_INTTOL);
				copy_x (lp->O->nstruct, minf->x, minf->bestx);
			}
		}

		if (EGlpNumIsLeq (minf->objectivebound, lpval))
		{
			*cp = 1;
		}
		else
		{
			b = bbnodealloc (&minf->ptrworld);
			init_bbnode (b);
			b->depth = active->depth + 1;
			b->id = minf->totalnodes;
			EGlpNumCopy (b->bound, lpval);
			ILL_SAFE_MALLOC (b->cstat, lp->O->nstruct, char);
			ILL_SAFE_MALLOC (b->rstat, lp->nrows, char);

			rval = ILLlib_getbasis (lp, b->cstat, b->rstat);
			ILL_CLEANUP_IF (rval);
			if (minf->pinf->dII_price == QS_PRICE_DSTEEP)
			{
				b->rownorms = EGlpNumAllocArray (lp->nrows);
				tval = ILLlib_getrownorms (lp, minf->pinf, b->rownorms);
				if (tval)
				{
					printf ("Row norms not available\n");
					fflush (stdout);
					printf ("A\n");
					exit (1);
					EGlpNumFreeArray (b->rownorms);
				}
			}
			ILL_SAFE_MALLOC (b->bound_indx, active->bound_cnt + 1, int);
			ILL_SAFE_MALLOC (b->lu, active->bound_cnt + 1, char);

			b->bounds = EGlpNumAllocArray (active->bound_cnt + 1);
			for (i = 0; i < active->bound_cnt; i++)
			{
				b->bound_indx[i] = active->bound_indx[i];
				b->lu[i] = active->lu[i];
				EGlpNumCopy (b->bounds[i], active->bounds[i]);
			}
			b->bound_indx[active->bound_cnt] = bvar;
			if (bdir == 'D')
				b->lu[active->bound_cnt] = 'U';
			else
				b->lu[active->bound_cnt] = 'L';
			EGlpNumCopy (b->bounds[active->bound_cnt], t);
			b->bound_cnt = active->bound_cnt + 1;

			rval = ILLutil_priority_insert (minf->que, (void *) b, &lpval,
																			&(b->handle));
			ILL_CLEANUP_IF (rval);

			put_bbnode (minf, b);
			minf->activenodes++;
		}
	}
	minf->totalnodes++;

	if (bdir == 'D')
	{
		rval = ILLlib_chgbnd (lp, bvar, 'U', oldt);
		ILL_CLEANUP_IF (rval);
	}
	else
	{
		rval = ILLlib_chgbnd (lp, bvar, 'L', oldt);
		ILL_CLEANUP_IF (rval);
	}

CLEANUP:

	EGlpNumClearVar (t);
	EGlpNumClearVar (lpval);
	EGlpNumClearVar (oldt);
	return rval;
}

static int fix_variables (
	lpinfo * lp,
	EGlpNum_t * bestval,
	bbnode * b,
	EGlpNum_t * wupper,
	EGlpNum_t * wlower,
	int *hit)
{
	int rval = 0;
	int i, nnew = 0;
	int nstruct = lp->O->nstruct;
	EGlpNum_t delta, lpval;
	int *new_indx = 0;
	char *new_lu = 0;
	EGlpNum_t *new_bounds = 0;
	EGlpNum_t *dj = 0;

	EGlpNumInitVar (delta);
	EGlpNumInitVar (lpval);

	*hit = 0;

	if (EGlpNumIsLess (*bestval, ILL_MAXDOUBLE))
	{
		rval = ILLlib_objval (lp, 0, &lpval);
		ILL_CLEANUP_IF (rval);
		//delta = bestval - lpval + ILL_INTTOL;
		EGlpNumCopy (delta, *bestval);
		EGlpNumSubTo (delta, lpval);
		EGlpNumAddTo (delta, ILL_INTTOL);

		ILL_SAFE_MALLOC (new_indx, nstruct, int);
		ILL_SAFE_MALLOC (new_lu, nstruct, char);

		dj = EGlpNumAllocArray (nstruct);
		new_bounds = EGlpNumAllocArray (nstruct);

		rval = ILLlib_solution (lp, 0, 0, 0, 0, 0, dj);
		ILL_CLEANUP_IF (rval);

		for (i = 0; i < nstruct; i++)
		{
			if (lp->O->intmarker[i])
			{
				if (EGlpNumIsNeqq (wlower[i], wupper[i]))
				{
					if (EGlpNumIsLess (delta, dj[i]))
					{
						EGlpNumSubTo (wupper[i], oneLpNum);
						rval = ILLlib_chgbnd (lp, i, 'U', wupper[i]);
						ILL_CLEANUP_IF (rval);
						new_indx[nnew] = i;
						new_lu[nnew] = 'U';
						EGlpNumCopy (new_bounds[nnew], wupper[i]);
						nnew++;
					}
					/*if (-dj[i] > delta) */
					EGlpNumSign (delta);
					if (EGlpNumIsLess (delta, dj[i]))
					{
						EGlpNumAddTo (wlower[i], oneLpNum);
						rval = ILLlib_chgbnd (lp, i, 'L', wlower[i]);
						ILL_CLEANUP_IF (rval);
						new_indx[nnew] = i;
						new_lu[nnew] = 'L';
						EGlpNumCopy (new_bounds[nnew], wlower[i]);
						nnew++;
					}
					EGlpNumSign (delta);
				}
			}
		}

		if (nnew)
		{
			b->bound_indx =
				EGrealloc (b->bound_indx, sizeof (int) * (b->bound_cnt + nnew));
			//rval = ILLutil_reallocrus_count ((void **) &(b->bound_indx),
			//                                 b->bound_cnt + nnew, sizeof (int));
			//ILL_CLEANUP_IF (rval);
			b->lu = EGrealloc (b->lu, sizeof (char) * (b->bound_cnt + nnew));
			//rval = ILLutil_reallocrus_count ((void **) &(b->lu),
			//                                 b->bound_cnt + nnew, sizeof (char));
			//ILL_CLEANUP_IF (rval);
			EGlpNumReallocArray (&(b->bounds), b->bound_cnt + nnew);
			for (i = 0; i < nnew; i++)
			{
				b->bound_indx[b->bound_cnt + i] = new_indx[i];
				b->lu[b->bound_cnt + i] = new_lu[i];
				EGlpNumCopy (b->bounds[b->bound_cnt + i], new_bounds[i]);
			}
			b->bound_cnt += nnew;
		}
	}

	*hit = nnew;

CLEANUP:

	ILL_IFFREE (new_indx, int);
	ILL_IFFREE (new_lu, char);

	EGlpNumFreeArray (dj);
	EGlpNumFreeArray (new_bounds);
	EGlpNumClearVar (delta);
	EGlpNumClearVar (lpval);
	return rval;
}

static void best_bbnode (
	mipinfo * minf,
	bbnode ** best)
{
#if 0
	bbnode *b;
	double bestval = ILL_MAXDOUBLE;

	for (b = minf->head_bbnode.next; b; b = b->next)
	{
		if (b->bound < bestval)
		{
			*best = b;
			bestval = b->bound;
		}
	}
#endif

	EGlpNum_t val;

	EGlpNumInitVar (val);
	ILLutil_priority_deletemin (minf->que, &val, (void **) best);
	EGlpNumClearVar (val);
}

static void put_bbnode (
	mipinfo * minf,
	bbnode * b)
{
	b->next = minf->head_bbnode.next;
	b->prev = &(minf->head_bbnode);
	if (b->next)
		b->next->prev = b;
	minf->head_bbnode.next = b;
}

static void remove_bbnode (
	bbnode * b)
{
	b->prev->next = b->next;
	if (b->next)
		b->next->prev = b->prev;
}

static int find_branch (
	mipinfo * minf,
	EGlpNum_t * x,
	EGlpNum_t * lpval,
	int *bvar,
	itcnt_t*itcnt)
{
	lpinfo *lp = minf->lp;
	int rval = 0;

	switch (minf->branching_rule)
	{
	case PENALTYBRANCH:
		rval = find_penalty_branch (lp, minf->pinf, x, minf->downpen,
																minf->uppen, lpval, bvar, itcnt);
		ILL_CLEANUP_IF (rval);
		break;
	case FIRSTBRANCH:
		find_first_branch (lp, x, bvar);
		break;
	case MIDDLEBRANCH:
		find_middle_branch (lp, x, bvar);
		break;
	case STRONGBRANCH:
		rval = find_strong_branch (lp, minf->pinf, x, bvar, itcnt);
		ILL_CLEANUP_IF (rval);
		break;
	default:
		fprintf (stderr, "Unknown branching rule.\n");
		rval = 1;
		goto CLEANUP;
	}

CLEANUP:

	ILL_RETURN (rval, "find_branch");
}

static void find_first_branch (
	lpinfo * lp,
	EGlpNum_t * x,
	int *bvar)
{
	int i, ibest = -1;
	ILLlpdata *qslp = lp->O;
	EGlpNum_t t;

	EGlpNumInitVar (t);

	for (i = 0; i < qslp->nstruct; i++)
	{
		if (qslp->intmarker[i])
		{
			/*t = ILLutil_our_frac (x[i]); */
			EGlpNumFloor (t, x[i]);
			EGlpNumSubTo (t, x[i]);
			EGlpNumSign (t);
			if ((EGlpNumIsNeqZero (t, ILL_INTTOL)) &&
					(EGlpNumIsNeq (t, oneLpNum, ILL_INTTOL)))
			{
				ibest = i;
				break;
			}
		}
	}
	*bvar = ibest;
	EGlpNumClearVar (t);
}

static void find_middle_branch (
	lpinfo * lp,
	EGlpNum_t * x,
	int *bvar)
{
	int i, ibest = -1;
	EGlpNum_t t, tbest;
	ILLlpdata *qlp = lp->O;

	EGlpNumInitVar (t);
	EGlpNumInitVar (tbest);
	EGlpNumSet (tbest, 0.5);

	for (i = 0; i < qlp->nstruct; i++)
	{
		if (qlp->intmarker[i])
		{
			/*t = ILLutil_our_frac (x[i]) - 0.5;
			 * if (t < 0.0)
			 * t = -t; */
			EGlpNumFloor (t, x[i]);
			EGlpNumMultUiTo (t, 2);
			EGlpNumSubTo (t, oneLpNum);
			EGlpNumDivUiTo (t, 2);
			if (EGlpNumIsLessZero (t))
				EGlpNumSign (t);
			/*if (t < tbest) */
			if (EGlpNumIsLess (t, tbest))
			{
				EGlpNumCopy (tbest, t);
				ibest = i;
			}
		}
	}

	/*if (tbest < (0.5 - ILL_INTTOL)) */
	EGlpNumAddTo (tbest, ILL_INTTOL);
	if (EGlpNumIsLessDbl (tbest, 0.5))
	{
		*bvar = ibest;
	}
	else
	{
		*bvar = -1;
	}
	EGlpNumClearVar (t);
	EGlpNumClearVar (tbest);
}

static int find_penalty_branch (
	lpinfo * lp,
	price_info * pinf,
	EGlpNum_t * x,
	EGlpNum_t * downpen,
	EGlpNum_t * uppen,
	EGlpNum_t * lpval,
	int *bvar,
	itcnt_t*itcnt)
{
	int rval = 0;
	int i, k, ibest = -1, ncand = 0, nneed = 0;
	ILLlpdata *qslp = lp->O;
	int *candidatelist = 0;
	int *needlist = 0;
	EGlpNum_t *fval = 0;
	EGlpNum_t *xlist = 0;
	EGlpNum_t *newdown = 0;
	EGlpNum_t *newup = 0;
	EGlpNum_t a, t, tbest;

	EGlpNumInitVar (a);
	EGlpNumInitVar (t);
	EGlpNumInitVar (tbest);
	EGlpNumCopy (tbest, ILL_MINDOUBLE);

	ILL_SAFE_MALLOC (candidatelist, qslp->nstruct, int);
	ILL_SAFE_MALLOC (needlist, qslp->nstruct, int);

	fval = EGlpNumAllocArray (qslp->nstruct);
	xlist = EGlpNumAllocArray (qslp->nstruct);
	for (i = 0; i < qslp->nstruct; i++)
	{
		if (qslp->intmarker[i])
		{
			/*fval[i] = x[i] - floor(x[i]); */
			EGlpNumFloor (fval[i], x[i]);
			EGlpNumSubTo (fval[i], x[i]);
			EGlpNumSign (fval[i]);
			if ((EGlpNumIsNeqZero (fval[i], ILL_INTTOL)) &&
					(EGlpNumIsNeq (fval[i], oneLpNum, ILL_INTTOL)))
			{
				candidatelist[ncand++] = i;
				/*if (downpen[i] == -1.0) */
				EGlpNumSign (downpen[i]);
				if (EGlpNumIsEqqual (downpen[i], oneLpNum))
				{
					EGlpNumCopy (xlist[nneed], x[i]);
					needlist[nneed++] = i;
				}
				EGlpNumSign (downpen[i]);
			}
		}
	}

	if (nneed > 0)
	{
		newdown = EGlpNumAllocArray (nneed);
		newup = EGlpNumAllocArray (nneed);
		rval = ILLlib_strongbranch (lp, pinf, needlist, nneed,
																0, newdown, newup,
																5 * STRONG_PIVOTS, ILL_MAXDOUBLE, itcnt);
		ILL_CLEANUP_IF (rval);

		for (i = 0; i < nneed; i++)
		{
			k = needlist[i];
			/*uppen[k] = (newup[i] - lpval) / (1.0 - fval[k]); */
			EGlpNumCopyDiff (uppen[k], newup[i], *lpval);
			EGlpNumCopyDiff (downpen[k], oneLpNum, fval[k]);
			EGlpNumDivTo (uppen[k], downpen[k]);
			/*downpen[k] = (newdown[i] - lpval) / fval[k]; */
			EGlpNumCopyDiffRatio (downpen[k], newdown[i], *lpval, fval[k]);

		}
	}

	for (i = 0; i < ncand; i++)
	{
		k = candidatelist[i];
		/*t = ILL_BRANCH_PENALTY_VAL (downpen[k], uppen[k], fval[k]); */
		EGlpNumCopy (t, downpen[k]);
		EGlpNumMultTo (t, fval[k]);
		EGlpNumCopyDiff (a, oneLpNum, fval[k]);
		EGlpNumMultTo (a, uppen[k]);
		if (EGlpNumIsLess (t, a))
		{
			EGlpNumMultUiTo (t, ILL_BRANCH_PENALTY_WEIGHT);
			EGlpNumAddTo (t, a);
		}
		else
		{
			EGlpNumMultUiTo (a, ILL_BRANCH_PENALTY_WEIGHT);
			EGlpNumAddTo (t, a);
		}
		EGlpNumDivUiTo (t, ILL_BRANCH_PENALTY_WEIGHT + 1);

		if (EGlpNumIsLess (tbest, t))
		{
			EGlpNumCopy (tbest, t);
			ibest = k;
		}
	}

	*bvar = ibest;

CLEANUP:

	EGlpNumClearVar (a);
	EGlpNumClearVar (t);
	EGlpNumClearVar (tbest);
	EGlpNumFreeArray (newdown);
	EGlpNumFreeArray (newup);
	EGlpNumFreeArray (fval);
	EGlpNumFreeArray (xlist);
	ILL_IFFREE (candidatelist, int);
	ILL_IFFREE (needlist, int);

	return rval;
}

static int find_strong_branch (
	lpinfo * lp,
	price_info * pinf,
	EGlpNum_t * x,
	int *bvar,
	itcnt_t*itcnt)
{
	int rval = 0;
	int i, ibest = -1, ncand = 0;
	int maxtrys = STRONG_CANDIDATES;
	EGlpNum_t t, tbest;
	ILLlpdata *qlp = lp->O;
	int *candidatelist = 0;
	int *newlist = 0;
	int *perm = 0;
	EGlpNum_t *tval = 0;
	EGlpNum_t *xlist = 0;
	EGlpNum_t *downpen = 0;
	EGlpNum_t *uppen = 0;
	ILLrandstate rstate;

	EGlpNumInitVar (t);
	EGlpNumInitVar (tbest);
	EGlpNumCopy (tbest, ILL_MINDOUBLE);

	ILLutil_sprand (999, &rstate);
	ILL_SAFE_MALLOC (candidatelist, qlp->nstruct, int);

	tval = EGlpNumAllocArray (qlp->nstruct);

	for (i = 0; i < qlp->nstruct; i++)
	{
		if (qlp->intmarker[i])
		{
			/*t = ILLutil_our_frac (x[i]) - 0.5;
			 * if (t < 0.0)
			 * t = -t; */
			EGlpNumFloor (t, x[i]);
			EGlpNumSubTo (t, x[i]);
			EGlpNumSign (t);
			EGlpNumMultUiTo (t, 2);
			EGlpNumSubTo (t, oneLpNum);
			if (EGlpNumIsLessZero (t))
				EGlpNumSign (t);
			/*if (t < (0.5 - ILL_INTTOL)) */
			if (EGlpNumIsNeq (t, oneLpNum, ILL_INTTOL))
			{
				candidatelist[ncand] = i;
				EGlpNumDivUiTo (t, 2);
				EGlpNumCopy (tval[ncand++], t);
			}
		}
	}

	if (ncand > 0)
	{
		if (ncand > maxtrys)
		{
			ILL_SAFE_MALLOC (perm, ncand, int);

			for (i = 0; i < ncand; i++)
			{
				perm[i] = i;
			}
			ILLutil_EGlpNum_rselect (perm, 0, ncand - 1, maxtrys, tval, &rstate);

			ILL_SAFE_MALLOC (newlist, maxtrys, int);

			for (i = 0; i < maxtrys; i++)
			{
				newlist[i] = candidatelist[perm[i]];
			}
			ILL_IFFREE (candidatelist, int);

			candidatelist = newlist;
			newlist = 0;
			ncand = maxtrys;
		}

		downpen = EGlpNumAllocArray (ncand);
		uppen = EGlpNumAllocArray (ncand);
		xlist = EGlpNumAllocArray (ncand);

		for (i = 0; i < ncand; i++)
		{
			EGlpNumCopy (xlist[i], x[candidatelist[i]]);
		}

		rval = ILLlib_strongbranch (lp, pinf, candidatelist, ncand,
																0, downpen, uppen, STRONG_PIVOTS,
																ILL_MAXDOUBLE, itcnt);
		ILL_CLEANUP_IF (rval);

		for (i = 0; i < ncand; i++)
		{
			/*t = ILL_BRANCH_STRONG_VAL (downpen[i], uppen[i]); */
			if (EGlpNumIsLess (downpen[i], uppen[i]))
			{
				EGlpNumCopy (t, downpen[i]);
				EGlpNumMultUiTo (t, ILL_BRANCH_STRONG_WEIGHT);
				EGlpNumAddTo (t, uppen[i]);
			}
			else
			{
				EGlpNumCopy (t, uppen[i]);
				EGlpNumMultUiTo (t, ILL_BRANCH_STRONG_WEIGHT);
				EGlpNumAddTo (t, downpen[i]);
			}
			EGlpNumDivUiTo (t, ILL_BRANCH_STRONG_WEIGHT + 1);
			if (EGlpNumIsLess (tbest, t))
			{
				EGlpNumCopy (tbest, t);
				ibest = candidatelist[i];
			}
		}
	}

	*bvar = ibest;


CLEANUP:

	EGlpNumClearVar (t);
	EGlpNumClearVar (tbest);
	EGlpNumFreeArray (tval);
	EGlpNumFreeArray (xlist);
	EGlpNumFreeArray (uppen);
	EGlpNumFreeArray (downpen);
	ILL_IFFREE (candidatelist, int);
	ILL_IFFREE (newlist, int);
	ILL_IFFREE (perm, int);

	ILL_RETURN (rval, "find_strong_branch");
}

static void check_integral (
	lpinfo * lp,
	EGlpNum_t * x,
	int *yesno)
{
	int i;
	EGlpNum_t t;
	ILLlpdata *qlp = lp->O;

	EGlpNumInitVar (t);

	for (i = 0; i < qlp->nstruct; i++)
	{
		if (qlp->intmarker[i])
		{
			/*t = ILLutil_our_frac (x[i]); */
			EGlpNumFloor (t, x[i]);
			EGlpNumSubTo (t, x[i]);
			EGlpNumSign (t);
			/*if (t > ILL_INTTOL && t < 1.0 - ILL_INTTOL) */
			if ((EGlpNumIsNeqZero (t, ILL_INTTOL)) &&
					(EGlpNumIsNeq (t, oneLpNum, ILL_INTTOL)))
			{
				*yesno = 0;
				EGlpNumClearVar (t);
				return;
			}
		}
	}

	*yesno = 1;
	EGlpNumClearVar (t);
}

static int plunge (
	mipinfo * minf,
	itcnt_t*itcnt)
{
	int rval = 0;
	int i, status;
	lpinfo *lp = minf->lp;
	ILLlpdata *qlp = minf->lp->O;
	EGlpNum_t *oldlower = 0;
	EGlpNum_t *oldupper = 0;

	if (minf->watch)
	{
		printf ("Plunging ...\n");
		fflush (stdout);
	}

	oldlower = EGlpNumAllocArray (qlp->nstruct);
	oldupper = EGlpNumAllocArray (qlp->nstruct);

	for (i = 0; i < qlp->nstruct; i++)
	{
		EGlpNumCopy (oldlower[i], minf->lower[i]);
		EGlpNumCopy (oldupper[i], minf->upper[i]);
	}

	rval = plunge_work (minf, 0, itcnt);
	ILL_CLEANUP_IF (rval);

	for (i = 0; i < qlp->nstruct; i++)
	{
		rval = ILLlib_chgbnd (lp, i, 'L', oldlower[i]);
		ILL_CLEANUP_IF (rval);
		rval = ILLlib_chgbnd (lp, i, 'U', oldupper[i]);
		ILL_CLEANUP_IF (rval);
		EGlpNumCopy (minf->lower[i], oldlower[i]);
		EGlpNumCopy (minf->upper[i], oldupper[i]);
	}

	rval = ILLlib_optimize (lp, 0, minf->pinf, DUAL_SIMPLEX, &status, 0, itcnt);
	ILL_CLEANUP_IF (rval);


CLEANUP:

	EGlpNumFreeArray (oldlower);
	EGlpNumFreeArray (oldupper);

	ILL_RETURN (rval, "plunge");
}

static int plunge_work (
	mipinfo * minf,
	int depth,
	itcnt_t*itcnt)
{
	int rval = 0;
	int bvar, status, count;
	EGlpNum_t lpval, val0, val1, int_tol;
	lpinfo *lp = minf->lp;

	EGlpNumInitVar (lpval);
	EGlpNumInitVar (val0);
	EGlpNumInitVar (val1);
	EGlpNumInitVar (int_tol);
	EGlpNumSet (int_tol, 0.001);

	rval = ILLlib_get_x (lp, 0, minf->x);
	ILL_CLEANUP_IF (rval);

	rval = round_variables (minf, &count, &int_tol /* 0.001 */ );
	ILL_CLEANUP_IF (rval);
	if (count)
	{
		rval = ILLlib_optimize (lp, 0, minf->pinf, DUAL_SIMPLEX, &status, 0, itcnt);
		ILL_CLEANUP_IF (rval);
		if (status != QS_LP_OPTIMAL)
		{
			goto CLEANUP;
		}
		rval = ILLlib_get_x (lp, 0, minf->x);
		ILL_CLEANUP_IF (rval);
	}

	find_middle_branch (lp, minf->x, &bvar);
	if (bvar == -1)
	{
		rval = ILLlib_objval (lp, 0, &lpval);
		ILL_CLEANUP_IF (rval);

		if (EGlpNumIsLess (lpval, minf->value))
		{
			printf ("Plunge Integral Solution: %.6f (Depth: %d)\n",
							EGlpNumToLf (lpval), depth);
			fflush (stdout);

			EGlpNumCopy (minf->value, lpval);
			EGlpNumCopyDiff (minf->objectivebound, lpval, ILL_INTTOL);
			copy_x (lp->O->nstruct, minf->x, minf->bestx);
		}
		goto CLEANUP;
	}

	EGlpNumOne (minf->lower[bvar]);
	rval = ILLlib_chgbnd (lp, bvar, 'L', oneLpNum);
	ILL_CLEANUP_IF (rval);
	rval = ILLlib_optimize (lp, 0, minf->pinf, DUAL_SIMPLEX, &status, 0, itcnt);
	ILL_CLEANUP_IF (rval);

	if (status == QS_LP_UNSOLVED)
	{
		printf ("Simplex did not solve the plunge LP\n");
		fflush (stdout);
		rval = 1;
		ILL_CLEANUP;
	}
	else if (status == QS_LP_INFEASIBLE)
	{
		EGlpNumCopy (val1, ILL_MAXDOUBLE);
	}
	else if (status == QS_LP_OPTIMAL)
	{
		rval = ILLlib_objval (lp, 0, &val1);
		ILL_CLEANUP_IF (rval);
	}
	else
	{
		ILL_CLEANUP;
	}

	rval = ILLlib_chgbnd (lp, bvar, 'L', zeroLpNum);
	ILL_CLEANUP_IF (rval);
	EGlpNumZero (minf->lower[bvar]);

	EGlpNumZero (minf->upper[bvar]);
	rval = ILLlib_chgbnd (lp, bvar, 'U', zeroLpNum);
	ILL_CLEANUP_IF (rval);
	rval = ILLlib_optimize (lp, 0, minf->pinf, DUAL_SIMPLEX, &status, 0, itcnt);
	ILL_CLEANUP_IF (rval);

	if (status == QS_LP_UNSOLVED)
	{
		printf ("Simplex did not solve the plunge LP\n");
		fflush (stdout);
		rval = 1;
		ILL_CLEANUP;
	}
	else if (status == QS_LP_INFEASIBLE)
	{
		EGlpNumCopy (val0, ILL_MAXDOUBLE);
	}
	else if (status == QS_LP_OPTIMAL)
	{
		rval = ILLlib_objval (lp, 0, &val0);
		ILL_CLEANUP_IF (rval);
	}
	else
	{
		ILL_CLEANUP;
	}

	rval = ILLlib_chgbnd (lp, bvar, 'U', oneLpNum);
	ILL_CLEANUP_IF (rval);
	EGlpNumCopy (minf->upper[bvar], oneLpNum);

	if (EGlpNumIsEqqual (val0, ILL_MAXDOUBLE) &&
			EGlpNumIsEqqual (val1, ILL_MAXDOUBLE))
	{
		ILL_CLEANUP;
	}

	if (EGlpNumIsLess (val0, val1))
	{
		EGlpNumZero (minf->upper[bvar]);
		rval = ILLlib_chgbnd (lp, bvar, 'U', zeroLpNum);
		ILL_CLEANUP_IF (rval);
		rval = ILLlib_optimize (lp, 0, minf->pinf, DUAL_SIMPLEX, &status, 0, itcnt);
		ILL_CLEANUP_IF (rval);
		rval = plunge_work (minf, depth + 1, itcnt);
		ILL_CLEANUP_IF (rval);
		rval = ILLlib_chgbnd (lp, bvar, 'U', oneLpNum);
		ILL_CLEANUP_IF (rval);
		EGlpNumOne (minf->upper[bvar]);
	}
	else
	{
		EGlpNumOne (minf->lower[bvar]);
		rval = ILLlib_chgbnd (lp, bvar, 'L', oneLpNum);
		ILL_CLEANUP_IF (rval);
		rval = ILLlib_optimize (lp, 0, minf->pinf, DUAL_SIMPLEX, &status, 0, itcnt);
		ILL_CLEANUP_IF (rval);
		rval = plunge_work (minf, depth + 1, itcnt);
		ILL_CLEANUP_IF (rval);
		rval = ILLlib_chgbnd (lp, bvar, 'L', zeroLpNum);
		ILL_CLEANUP_IF (rval);
		EGlpNumZero (minf->lower[bvar]);
	}

CLEANUP:

	EGlpNumClearVar (lpval);
	EGlpNumClearVar (val0);
	EGlpNumClearVar (val1);
	EGlpNumClearVar (int_tol);
	ILL_RETURN (rval, "plunge_work");
}

static int round_variables (
	mipinfo * minf,
	int *count,
	EGlpNum_t * tol)
{
	int rval = 0;
	int i, hit = 0;
	lpinfo *lp = minf->lp;
	ILLlpdata *qlp = lp->O;

	*count = 0;

	for (i = 0; i < qlp->nstruct; i++)
	{
		if (qlp->intmarker[i])
		{
			if (EGlpNumIsNeqq (minf->lower[i], minf->upper[i]))
			{
				if (EGlpNumIsLess (minf->x[i], *tol))
				{
					EGlpNumZero (minf->upper[i]);
					rval = ILLlib_chgbnd (lp, i, 'U', zeroLpNum);
					ILL_CLEANUP_IF (rval);
					hit++;
				}
				else if (EGlpNumIsEqual (minf->x[i], oneLpNum, *tol))
				{
					EGlpNumOne (minf->lower[i]);
					rval = ILLlib_chgbnd (lp, i, 'L', oneLpNum);
					ILL_CLEANUP_IF (rval);
					hit++;
				}
			}
		}
	}
	*count = hit;

CLEANUP:

	ILL_RETURN (rval, "round_variables");
}

static void copy_x (
	int nstruct,
	EGlpNum_t * from_x,
	EGlpNum_t * to_x)
{
	int j;

	for (j = 0; j < nstruct; j++)
	{
		EGlpNumCopy (to_x[j], from_x[j]);
	}
}

static void init_mipinfo (
	mipinfo * minf)
{
	if (minf)
	{
		minf->depth = 0;
		minf->totalnodes = 0;
		minf->activenodes = 0;
		minf->totalpivots = 0;
		minf->lastpivots = 0;
		minf->downpen = 0;
		minf->uppen = 0;
		minf->x = 0;
		minf->bestx = 0;
		minf->lower = 0;
		minf->upper = 0;
		minf->lp = 0;
		minf->pinf = 0;
		minf->head_bbnode.prev = 0;
		minf->head_bbnode.next = 0;
		minf->que = 0;
		minf->branching_rule = /* MIDDLEBRANCH */ STRONGBRANCH;
		minf->watch = 1;
		EGlpNumInitVar (minf->objectivebound);
		EGlpNumInitVar (minf->value);
		EGlpNumCopy (minf->objectivebound, ILL_MAXDOUBLE);
		EGlpNumCopy (minf->value, ILL_MAXDOUBLE);
		ILLptrworld_init (&minf->ptrworld);
	}
}

static void free_mipinfo (
	mipinfo * minf)
{
	int total, onlist;

	if (minf)
	{
		EGlpNumFreeArray (minf->downpen);
		EGlpNumFreeArray (minf->uppen);
		EGlpNumFreeArray (minf->x);
		EGlpNumFreeArray (minf->bestx);
		EGlpNumFreeArray (minf->lower);
		EGlpNumFreeArray (minf->upper);
		bbnode_listfree (&minf->ptrworld, minf->head_bbnode.next);
		if (bbnode_check_leaks (&minf->ptrworld, &total, &onlist))
		{
			fprintf (stderr, "WARNING: %d outstanding bbnodes\n", total - onlist);
		}
		ILLptrworld_delete (&minf->ptrworld);
		EGlpNumClearVar ((minf->objectivebound));
		EGlpNumClearVar ((minf->value));
		memset (minf, 0, sizeof (mipinfo));
		//init_mipinfo (minf);
	}
}

static void init_bbnode (
	bbnode * b)
{
	if (b)
	{
		b->next = 0;
		b->prev = 0;
		b->id = 0;
		b->depth = 0;
		b->handle = 0;
		b->cstat = 0;
		b->rstat = 0;
		b->rownorms = 0;
		b->bound_cnt = 0;
		b->bound_indx = 0;
		b->lu = 0;
		b->bounds = 0;
		EGlpNumInitVar ((b->bound));
		EGlpNumCopy (b->bound, ILL_MINDOUBLE);
	}
}

static void free_bbnode (
	bbnode * b)
{
	if (b)
	{
		EGlpNumFreeArray (b->rownorms);
		EGlpNumFreeArray (b->bounds);
		ILL_IFFREE (b->cstat, char);
		ILL_IFFREE (b->rstat, char);
		ILL_IFFREE (b->bound_indx, int);
		ILL_IFFREE (b->lu, char);

		EGlpNumClearVar ((b->bound));
		memset (b, 0, sizeof (bbnode));
	}
}
