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
/*    int EGLPNUM_TYPENAME_ILLmip_bfs (EGLPNUM_TYPENAME_lpinfo *lp, double *val, double *x)                   */
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
#include "logging-private.h"

#include "allocrus.h"
#include "eg_lpnum.h"
#include "eg_io.h"
#include "except.h"
#include "zeit.h"

#include "priority_EGLPNUM_TYPENAME.h"
#include "sortrus_EGLPNUM_TYPENAME.h"
#include "lpdata_EGLPNUM_TYPENAME.h"
#include "lpdefs_EGLPNUM_TYPENAME.h"
#include "simplex_EGLPNUM_TYPENAME.h"
#include "binary_EGLPNUM_TYPENAME.h"
#include "price_EGLPNUM_TYPENAME.h"
#include "lib_EGLPNUM_TYPENAME.h"
#include "qstruct_EGLPNUM_TYPENAME.h"
#include "qsopt_EGLPNUM_TYPENAME.h"

/*#define  ILL_INTTOL (0.000001)*/
#define ILL_INTTOL EGLPNUM_TYPENAME_PFEAS_TOLER

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
	EGLPNUM_TYPE bound;
	char *cstat;
	char *rstat;
	EGLPNUM_TYPE *rownorms;
	int rownorms_size;
	int bound_cnt;
	int *bound_indx;
	char *lu;
	EGLPNUM_TYPE *bounds;
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
	EGLPNUM_TYPE objectivebound;
	EGLPNUM_TYPE value;
	EGLPNUM_TYPE *downpen;
	EGLPNUM_TYPE *uppen;
	EGLPNUM_TYPE *x;
	EGLPNUM_TYPE *bestx;
	EGLPNUM_TYPE *orig_lower;
	EGLPNUM_TYPE *orig_upper;
	EGLPNUM_TYPE *lower;
	EGLPNUM_TYPE *upper;
	int nstruct;									/* size of all EGLPNUM_TYPE arrays */
	EGLPNUM_TYPENAME_lpinfo *lp;
	EGLPNUM_TYPENAME_price_info *pinf;
	bbnode head_bbnode;
	EGLPNUM_TYPENAME_ILLpriority *que;
	ILLptrworld ptrworld;
}
mipinfo;


ILL_PTRWORLD_ROUTINES (bbnode, bbnodealloc, bbnode_bulkalloc, bbnodefree)
ILL_PTRWORLD_LISTFREE_ROUTINE (bbnode, bbnode_listfree, bbnodefree)
ILL_PTRWORLD_LEAKS_ROUTINE (bbnode, bbnode_check_leaks, depth, int)
static void cleanup_mip ( mipinfo * minf), 
		choose_initial_price ( EGLPNUM_TYPENAME_price_info * pinf), 
		best_bbnode ( mipinfo * minf, bbnode ** best),
		put_bbnode ( mipinfo * minf, bbnode * b),
		remove_bbnode ( bbnode * b),
		find_first_branch ( EGLPNUM_TYPENAME_lpinfo * lp, EGLPNUM_TYPE * x, int *bvar),
		find_middle_branch ( EGLPNUM_TYPENAME_lpinfo * lp, EGLPNUM_TYPE * x, int *bvar),
		check_integral ( EGLPNUM_TYPENAME_lpinfo * lp, EGLPNUM_TYPE * x, int *yesno),
		copy_x ( int nstruct, EGLPNUM_TYPE * from_x, EGLPNUM_TYPE * to_x),
		init_mipinfo ( mipinfo * minf),
		free_mipinfo ( mipinfo * minf),
		init_bbnode ( bbnode * b),
		free_bbnode ( bbnode * b);

static int startup_mip ( mipinfo * minf, EGLPNUM_TYPENAME_lpinfo * lp, EGLPNUM_TYPENAME_price_info * pinf,
			EGLPNUM_TYPE * lpval, itcnt_t*itcnt),
		run_bfs ( mipinfo * minf, itcnt_t*itcnt),
		process_bfs_bbnode ( mipinfo * minf, bbnode * b, itcnt_t*itcnt),
		child_work ( mipinfo * minf, bbnode * active, int bvar, int bdir,
			EGLPNUM_TYPE * cval, int *cp, itcnt_t*itcnt),
		fix_variables ( EGLPNUM_TYPENAME_lpinfo * lp, EGLPNUM_TYPE * bestval, bbnode * b,
			EGLPNUM_TYPE * wupper, EGLPNUM_TYPE * wlower, int *hit),
		find_branch ( mipinfo * minf, EGLPNUM_TYPE * x, EGLPNUM_TYPE * lpval,
			int *bvar, itcnt_t*itcnt),
		find_penalty_branch ( EGLPNUM_TYPENAME_lpinfo * lp, EGLPNUM_TYPENAME_price_info * pinf, EGLPNUM_TYPE * x,
			EGLPNUM_TYPE * downpen, EGLPNUM_TYPE * uppen, EGLPNUM_TYPE * lpval, int *bvar,
			itcnt_t*itcnt),
		find_strong_branch ( EGLPNUM_TYPENAME_lpinfo * lp, EGLPNUM_TYPENAME_price_info * pinf, EGLPNUM_TYPE * x,
			int *bvar, itcnt_t*itcnt),
		plunge ( mipinfo * minf, itcnt_t*itcnt),
		plunge_work ( mipinfo * minf, int depth, itcnt_t*itcnt),
		round_variables ( mipinfo * minf, int *count, EGLPNUM_TYPE * tol);

static void choose_initial_price ( EGLPNUM_TYPENAME_price_info * pinf)
{
	pinf->pI_price = QS_PRICE_PSTEEP;
	pinf->pII_price = QS_PRICE_PSTEEP;
	pinf->dI_price = QS_PRICE_DSTEEP;
	pinf->dII_price = QS_PRICE_DSTEEP;
}

int EGLPNUM_TYPENAME_ILLmip_bfs (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPE * val,
	EGLPNUM_TYPE * x,
	itcnt_t*itcnt)
{
	int tval, rval = 0;
	EGLPNUM_TYPENAME_price_info pinf;
	mipinfo minf;
	bbnode *b;
	EGLPNUM_TYPE lpval;
	double szeit = ILLutil_zeit ();

	EGLPNUM_TYPENAME_EGlpNumInitVar (lpval);
	EGLPNUM_TYPENAME_EGlpNumInitVar (pinf.htrigger);

	EGLPNUM_TYPENAME_ILLprice_init_pricing_info (&pinf);
	init_mipinfo (&minf);

	if (!lp)
	{
		QSlog("EGLPNUM_TYPENAME_ILLmip_bfs called without an LP");
		rval = 1;
		goto CLEANUP;
	}

	rval = startup_mip (&minf, lp, &pinf, &lpval, itcnt);
	ILL_CLEANUP_IF (rval);

	ILL_SAFE_MALLOC (minf.que, 1, EGLPNUM_TYPENAME_ILLpriority);
	rval = EGLPNUM_TYPENAME_ILLutil_priority_init (minf.que, lp->O->nstruct + 1);
	ILL_CLEANUP_IF (rval);

	b = bbnodealloc (&minf.ptrworld);
	init_bbnode (b);
	b->depth = 0;
	b->id = minf.totalnodes++;
	EGLPNUM_TYPENAME_EGlpNumCopy (b->bound, lpval);
	ILL_SAFE_MALLOC (b->cstat, lp->O->nstruct, char);
	ILL_SAFE_MALLOC (b->rstat, lp->nrows, char);

	rval = EGLPNUM_TYPENAME_ILLlib_getbasis (lp, b->cstat, b->rstat);
	ILL_CLEANUP_IF (rval);

	if (pinf.dII_price == QS_PRICE_DSTEEP)
	{
		b->rownorms = EGLPNUM_TYPENAME_EGlpNumAllocArray (lp->nrows);
		tval = EGLPNUM_TYPENAME_ILLlib_getrownorms (lp, &pinf, b->rownorms);
		if (tval)
		{
			QSlog("Row norms not available");
			EGLPNUM_TYPENAME_EGlpNumFreeArray (b->rownorms);
		}
	}

	rval = EGLPNUM_TYPENAME_ILLutil_priority_insert (minf.que, (void *) b, &lpval, &(b->handle));
	ILL_CLEANUP_IF (rval);

	b->prev = &(minf.head_bbnode);
	b->next = 0;
	minf.head_bbnode.next = b;
	minf.activenodes++;

	minf.branching_rule = PENALTYBRANCH;

	rval = run_bfs (&minf, itcnt);
	ILL_CLEANUP_IF (rval);

	QSlog("Total Number of Nodes: %d", minf.totalnodes);
	QSlog("Total Number of Pivots: %d", minf.totalpivots);
	QSlog("BFS MIP Runing Time: %.2f seconds", ILLutil_zeit () - szeit);

	EGLPNUM_TYPENAME_EGlpNumCopy (*val, minf.value);
	if (minf.objsense == EGLPNUM_TYPENAME_ILL_MAX)
		EGLPNUM_TYPENAME_EGlpNumSign (*val);

	if (x && EGLPNUM_TYPENAME_EGlpNumIsNeqq (minf.value, EGLPNUM_TYPENAME_ILL_MAXDOUBLE))
	{
		copy_x (lp->O->nstruct, minf.bestx, x);
	}

CLEANUP:

	if (minf.que)
	{
		EGLPNUM_TYPENAME_ILLutil_priority_free (minf.que);
		ILL_IFFREE (minf.que, EGLPNUM_TYPENAME_ILLpriority);
	}
	cleanup_mip (&minf);
	free_mipinfo (&minf);
	EGLPNUM_TYPENAME_ILLprice_free_pricing_info (&pinf);
	EGLPNUM_TYPENAME_EGlpNumClearVar (lpval);
	EGLPNUM_TYPENAME_EGlpNumClearVar (pinf.htrigger);
	ILL_RETURN (rval, "EGLPNUM_TYPENAME_ILLmip_bfs");
}

static int startup_mip (
	mipinfo * minf,
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_price_info * pinf,
	EGLPNUM_TYPE * lpval,
	itcnt_t*itcnt)
{
	int rval = 0;
	int i, col, status, intcount = 0;
	EGLPNUM_TYPE val;
	EGLPNUM_TYPENAME_ILLlpdata *qlp;

	EGLPNUM_TYPENAME_EGlpNumInitVar (val);

	choose_initial_price (pinf);

	qlp = lp->O;

	rval = EGLPNUM_TYPENAME_ILLlib_optimize (lp, 0, pinf, DUAL_SIMPLEX, &status, 0, itcnt);
	ILL_CLEANUP_IF (rval);

	minf->totalpivots += EGLPNUM_TYPENAME_ILLlib_iter (lp);

	rval = EGLPNUM_TYPENAME_ILLlib_objval (lp, 0, &val);
	ILL_CLEANUP_IF (rval);

	QSlog("LP Value: %.6f", EGLPNUM_TYPENAME_EGlpNumToLf (val));
	if (lpval)
		EGLPNUM_TYPENAME_EGlpNumCopy (*lpval, val);

	if (qlp->intmarker)
	{
		for (i = 0; i < qlp->nstruct; i++)
		{
			if (qlp->intmarker[i])
			{
				col = qlp->structmap[i];
				intcount++;
				if (EGLPNUM_TYPENAME_EGlpNumIsEqqual (qlp->lower[col], EGLPNUM_TYPENAME_ILL_MINDOUBLE)
						|| EGLPNUM_TYPENAME_EGlpNumIsEqqual (qlp->upper[col], EGLPNUM_TYPENAME_ILL_MAXDOUBLE))
				{
					QSlog("Instance has unbounded integer variable");
					rval = 1;
					goto CLEANUP;
				}
			}
		}
	}

	if (intcount == 0)
	{
		QSlog("No integer variables");
		rval = 1;
		goto CLEANUP;
	}
	else
	{
		QSlog("%d integer variables", intcount);
	}

	if (qlp->sinfo)
	{															/* Free the presolve LP and work with orginal */
		EGLPNUM_TYPENAME_ILLlp_sinfo_free (qlp->sinfo);
		ILL_IFFREE (qlp->sinfo, EGLPNUM_TYPENAME_ILLlp_sinfo);
	}


	minf->lp = lp;
	minf->pinf = pinf;
	minf->objsense = qlp->objsense;
	if (qlp->objsense == EGLPNUM_TYPENAME_ILL_MAX)
	{															/* MIP codes work with min */
		for (i = 0; i < lp->ncols; i++)
		{
			EGLPNUM_TYPENAME_EGlpNumCopyNeg (qlp->obj[i], qlp->obj[i]);
		}
		qlp->objsense = EGLPNUM_TYPENAME_ILL_MIN;
	}

	minf->x = EGLPNUM_TYPENAME_EGlpNumAllocArray (qlp->nstruct);
	minf->bestx = EGLPNUM_TYPENAME_EGlpNumAllocArray (qlp->nstruct);
	minf->lower = EGLPNUM_TYPENAME_EGlpNumAllocArray (qlp->nstruct);
	minf->upper = EGLPNUM_TYPENAME_EGlpNumAllocArray (qlp->nstruct);
	minf->orig_lower = EGLPNUM_TYPENAME_EGlpNumAllocArray (qlp->nstruct);
	minf->orig_upper = EGLPNUM_TYPENAME_EGlpNumAllocArray (qlp->nstruct);
	minf->downpen = EGLPNUM_TYPENAME_EGlpNumAllocArray (qlp->nstruct);
	minf->uppen = EGLPNUM_TYPENAME_EGlpNumAllocArray (qlp->nstruct);
	minf->nstruct = qlp->nstruct;

	rval = EGLPNUM_TYPENAME_ILLlib_get_x (lp, 0, minf->x);
	ILL_CLEANUP_IF (rval);

	for (i = 0; i < qlp->nstruct; i++)
	{
		EGLPNUM_TYPENAME_EGlpNumCopy (minf->lower[i], qlp->lower[i]);
		EGLPNUM_TYPENAME_EGlpNumCopy (minf->upper[i], qlp->upper[i]);
		EGLPNUM_TYPENAME_EGlpNumCopy (minf->orig_lower[i], qlp->lower[i]);
		EGLPNUM_TYPENAME_EGlpNumCopy (minf->orig_upper[i], qlp->upper[i]);
		EGLPNUM_TYPENAME_EGlpNumOne (minf->downpen[i]);
		EGLPNUM_TYPENAME_EGlpNumOne (minf->uppen[i]);
		EGLPNUM_TYPENAME_EGlpNumSign (minf->downpen[i]);
		EGLPNUM_TYPENAME_EGlpNumSign (minf->uppen[i]);
	}


CLEANUP:

	EGLPNUM_TYPENAME_EGlpNumClearVar (val);
	ILL_RETURN (rval, "startup_mip");
}

static void cleanup_mip (
	mipinfo * minf)
{
	int i;
	EGLPNUM_TYPENAME_ILLlpdata *qslp;

	if (minf && minf->lp)
	{
		qslp = minf->lp->O;
		if (minf->objsense == EGLPNUM_TYPENAME_ILL_MAX)
		{
			for (i = 0; i < minf->lp->ncols; i++)
			{
				EGLPNUM_TYPENAME_EGlpNumSign (qslp->obj[i]);
			}
			qslp->objsense = EGLPNUM_TYPENAME_ILL_MIN;
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
	EGLPNUM_TYPENAME_lpinfo *lp = minf->lp;
	EGLPNUM_TYPENAME_ILLlp_basis B;
	int status, bvar = 0;
	int i, j, hit, dnp = 0, upp = 0;
	int nstruct = lp->O->nstruct;
	EGLPNUM_TYPE t, lpval, dnval, upval;
	EGLPNUM_TYPE *wupper = 0;
	EGLPNUM_TYPE *wlower = 0;
	int rval = 0;

	EGLPNUM_TYPENAME_EGlpNumInitVar (t);
	EGLPNUM_TYPENAME_EGlpNumInitVar (lpval);
	EGLPNUM_TYPENAME_EGlpNumInitVar (dnval);
	EGLPNUM_TYPENAME_EGlpNumInitVar (upval);

	EGLPNUM_TYPENAME_ILLlp_basis_init (&B);

	if (minf->watch > 1)
	{
		QSlog("Node %4d: %.3f", active->id, EGLPNUM_TYPENAME_EGlpNumToLf (active->bound));
		if (EGLPNUM_TYPENAME_EGlpNumIsNeqq (minf->value, EGLPNUM_TYPENAME_ILL_MAXDOUBLE))
			QSlog(" %.3f", EGLPNUM_TYPENAME_EGlpNumToLf (minf->value));
		else
			QSlog("  None");
		QSlog(", Active %d ", minf->activenodes);
	}
	else if (minf->watch == 1)
	{
		if (minf->lastpivots > 1000)
		{
			minf->lastpivots = 0;
			QSlog("Pivots %d, Active Nodes %d, Bound %.3f, Soln ",
									minf->totalpivots, minf->activenodes,
									EGLPNUM_TYPENAME_EGlpNumToLf (active->bound));
			if (!EGLPNUM_TYPENAME_EGlpNumIsLess (minf->value, EGLPNUM_TYPENAME_ILL_MAXDOUBLE))
				QSlog("%.3f", EGLPNUM_TYPENAME_EGlpNumToLf (minf->value));
			else
				QSlog("None");
		}
	}

	if (EGLPNUM_TYPENAME_EGlpNumIsLeq (minf->objectivebound, active->bound))
	{
		if (minf->watch > 1)
		{
			QSlog("  Node can be purged");
		}
		goto CLEANUP;
	}

	/*  Set the LP bounds for the node. */

	wlower = EGLPNUM_TYPENAME_EGlpNumAllocArray (nstruct);
	wupper = EGLPNUM_TYPENAME_EGlpNumAllocArray (nstruct);

	for (i = 0; i < nstruct; i++)
	{
		EGLPNUM_TYPENAME_EGlpNumCopy (wlower[i], minf->orig_lower[i]);
		EGLPNUM_TYPENAME_EGlpNumCopy (wupper[i], minf->orig_upper[i]);
	}
	for (i = 0; i < active->bound_cnt; i++)
	{
		j = active->bound_indx[i];
		if (active->lu[i] == 'L')
			EGLPNUM_TYPENAME_EGlpNumCopy (wlower[j], active->bounds[i]);
		else
			EGLPNUM_TYPENAME_EGlpNumCopy (wupper[j], active->bounds[i]);
	}

	if (active->bound_cnt > 0)
	{
		rval = EGLPNUM_TYPENAME_ILLlib_chgbnds (lp, active->bound_cnt, active->bound_indx,
													 active->lu, active->bounds);
		ILL_CLEANUP_IF (rval);
	}

	/*  Solve the LP. */

	rval = EGLPNUM_TYPENAME_ILLlib_loadbasis (&B, nstruct, lp->nrows, active->cstat,
													 active->rstat);
	ILL_CLEANUP_IF (rval);
	if (active->rownorms)
	{
		B.rownorms = EGLPNUM_TYPENAME_EGlpNumAllocArray (lp->nrows);
		for (i = 0; i < lp->nrows; i++)
		{
			EGLPNUM_TYPENAME_EGlpNumCopy (B.rownorms[i], active->rownorms[i]);
		}
	}

	rval = EGLPNUM_TYPENAME_ILLlib_optimize (lp, &B, minf->pinf, DUAL_SIMPLEX, &status, 0, itcnt);
	ILL_CLEANUP_IF (rval);

	minf->totalpivots += EGLPNUM_TYPENAME_ILLlib_iter (lp);
	minf->lastpivots += EGLPNUM_TYPENAME_ILLlib_iter (lp);

	if (status == QS_LP_UNSOLVED)
	{
		QSlog("Simplex did not solve the LP");
		rval = 1;
		ILL_CLEANUP;
	}

	if (status == QS_LP_INFEASIBLE)
	{
		QSlog("  Infeasible LP, should have been purged earlier");
		rval = 1;
		ILL_CLEANUP;
	}

	if (active->depth < 0)
	{
		for (i = 0; i < nstruct; i++)
		{
			EGLPNUM_TYPENAME_EGlpNumCopy (minf->lower[i], wlower[i]);
			EGLPNUM_TYPENAME_EGlpNumCopy (minf->upper[i], wupper[i]);
		}
		rval = plunge (minf, itcnt);
		ILL_CLEANUP_IF (rval);
	}

	/*  Fix variables. */

	if (EGLPNUM_TYPENAME_EGlpNumIsLess (minf->value, EGLPNUM_TYPENAME_ILL_MAXDOUBLE))
	{
		rval = fix_variables (lp, &(minf->value), active, wupper, wlower, &hit);
		ILL_CLEANUP_IF (rval);

		if (hit)
		{
			rval = EGLPNUM_TYPENAME_ILLlib_optimize (lp, &B, minf->pinf, DUAL_SIMPLEX, &status, 0, itcnt);
			ILL_CLEANUP_IF (rval);

			minf->totalpivots += EGLPNUM_TYPENAME_ILLlib_iter (lp);
			minf->lastpivots += EGLPNUM_TYPENAME_ILLlib_iter (lp);

			if (status == QS_LP_UNSOLVED)
			{
				QSlog("Simplex did not solve the LP");
				rval = 1;
				ILL_CLEANUP;
			}

			if (status == QS_LP_INFEASIBLE)
			{
				QSlog("  Infeasible LP after fixing");
				rval = 1;
				ILL_CLEANUP;
			}
		}
	}


	/*  Branch. */

	rval = EGLPNUM_TYPENAME_ILLlib_get_x (lp, 0, minf->x);
	ILL_CLEANUP_IF (rval);

	rval = EGLPNUM_TYPENAME_ILLlib_objval (lp, 0, &lpval);
	ILL_CLEANUP_IF (rval);

	rval = find_branch (minf, minf->x, &lpval, &bvar, itcnt);
	ILL_CLEANUP_IF (rval);

	if (bvar == -1)
	{
		QSlog("Found integral solution: %f", EGLPNUM_TYPENAME_EGlpNumToLf (lpval));
		if (EGLPNUM_TYPENAME_EGlpNumIsLess (lpval, minf->value))
		{
			EGLPNUM_TYPENAME_EGlpNumCopy (minf->value, lpval);
			EGLPNUM_TYPENAME_EGlpNumCopy (minf->objectivebound, lpval);
			EGLPNUM_TYPENAME_EGlpNumSubTo (minf->objectivebound, ILL_INTTOL);
			copy_x (nstruct, minf->x, minf->bestx);
		}
	}
	else
	{
		/* Create down child */

		rval = child_work (minf, active, bvar, 'D', &dnval, &dnp, itcnt);
		ILL_CLEANUP_IF (rval);

		/* Restore parent basis */

		rval = EGLPNUM_TYPENAME_ILLlib_loadbasis (&B, nstruct, lp->nrows, active->cstat,
														 active->rstat);
		ILL_CLEANUP_IF (rval);
		if (active->rownorms)
		{
			B.rownorms = EGLPNUM_TYPENAME_EGlpNumAllocArray (lp->nrows);
			for (i = 0; i < lp->nrows; i++)
			{
				EGLPNUM_TYPENAME_EGlpNumCopy (B.rownorms[i], active->rownorms[i]);
			}
		}

		/* Create up child */

		rval = child_work (minf, active, bvar, 'U', &upval, &upp, itcnt);
		ILL_CLEANUP_IF (rval);

		if (minf->watch > 1)
		{
			if (EGLPNUM_TYPENAME_EGlpNumIsEqqual (dnval, EGLPNUM_TYPENAME_ILL_MAXDOUBLE))
			{
				QSlog("DN->XXX");
			}
			else
			{
				QSlog("DN->%.3f%c", EGLPNUM_TYPENAME_EGlpNumToLf (dnval), dnp ? 'X' : ' ');
			}
			if (EGLPNUM_TYPENAME_EGlpNumIsEqqual (upval, EGLPNUM_TYPENAME_ILL_MAXDOUBLE))
			{
				QSlog("UP->XXX");
			}
			else
			{
				QSlog("UP->%.3f%c", EGLPNUM_TYPENAME_EGlpNumToLf (upval), upp ? 'X' : ' ');
			}
		}
	}

	/* Set the LP bounds back to original values */

	for (i = 0; i < active->bound_cnt; i++)
	{
		if (active->lu[i] == 'L')
			EGLPNUM_TYPENAME_EGlpNumCopy (t, minf->orig_lower[active->bound_indx[i]]);
		else
			EGLPNUM_TYPENAME_EGlpNumCopy (t, minf->orig_upper[active->bound_indx[i]]);

		rval = EGLPNUM_TYPENAME_ILLlib_chgbnd (lp, active->bound_indx[i], active->lu[i], t);
		ILL_CLEANUP_IF (rval);
	}

CLEANUP:

	EGLPNUM_TYPENAME_EGlpNumFreeArray (wlower);
	EGLPNUM_TYPENAME_EGlpNumFreeArray (wupper);
	EGLPNUM_TYPENAME_ILLlp_basis_free (&B);
	EGLPNUM_TYPENAME_EGlpNumClearVar (t);
	EGLPNUM_TYPENAME_EGlpNumClearVar (lpval);
	EGLPNUM_TYPENAME_EGlpNumClearVar (dnval);
	EGLPNUM_TYPENAME_EGlpNumClearVar (upval);
	ILL_RETURN (rval, "process_bfs_bbnode");
}

static int child_work (
	mipinfo * minf,
	bbnode * active,
	int bvar,
	int bdir,
	EGLPNUM_TYPE * cval,
	int *cp,
	itcnt_t*itcnt)
{
	int tval, rval = 0;
	int i, status, intsol;
	EGLPNUM_TYPE t, oldt, lpval;
	EGLPNUM_TYPE *xi = &(minf->x[bvar]);
	EGLPNUM_TYPENAME_lpinfo *lp = minf->lp;
	bbnode *b;

	EGLPNUM_TYPENAME_EGlpNumInitVar (t);
	EGLPNUM_TYPENAME_EGlpNumInitVar (lpval);
	EGLPNUM_TYPENAME_EGlpNumInitVar (oldt);

	*cp = 0;

	if (bdir == 'D')
	{
		rval = EGLPNUM_TYPENAME_ILLlib_getbnd (lp, bvar, 'U', &oldt);
		ILL_CLEANUP_IF (rval);
		EGLPNUM_TYPENAME_EGlpNumFloor (t, *xi);
		rval = EGLPNUM_TYPENAME_ILLlib_chgbnd (lp, bvar, 'U', t);
		ILL_CLEANUP_IF (rval);
	}
	else
	{
		rval = EGLPNUM_TYPENAME_ILLlib_getbnd (lp, bvar, 'L', &oldt);
		ILL_CLEANUP_IF (rval);
		EGLPNUM_TYPENAME_EGlpNumCeil (t, *xi);
		rval = EGLPNUM_TYPENAME_ILLlib_chgbnd (lp, bvar, 'L', t);
		ILL_CLEANUP_IF (rval);
	}

	rval = EGLPNUM_TYPENAME_ILLlib_optimize (lp, 0, minf->pinf, DUAL_SIMPLEX, &status, 0, itcnt);
	ILL_CLEANUP_IF (rval);

	minf->totalpivots += EGLPNUM_TYPENAME_ILLlib_iter (lp);
	minf->lastpivots += EGLPNUM_TYPENAME_ILLlib_iter (lp);

	if (status == QS_LP_UNSOLVED)
	{
		QSlog("Simplex did not solve Child LP");
		rval = 1;
		ILL_CLEANUP;
	}

	if (status == QS_LP_INFEASIBLE)
	{
		EGLPNUM_TYPENAME_EGlpNumCopy (*cval, EGLPNUM_TYPENAME_ILL_MAXDOUBLE);
		*cp = 1;
	}
	else
	{
		rval = EGLPNUM_TYPENAME_ILLlib_objval (lp, 0, &lpval);
		ILL_CLEANUP_IF (rval);
		EGLPNUM_TYPENAME_EGlpNumCopy (*cval, lpval);

		/* What about the x vector?  Bico - 020531 */

		check_integral (lp, minf->x, &intsol);
		if (intsol)
		{
			if (EGLPNUM_TYPENAME_EGlpNumIsLess (lpval, minf->value))
			{
				QSlog("Found integral solution: %f", EGLPNUM_TYPENAME_EGlpNumToLf (lpval));
				EGLPNUM_TYPENAME_EGlpNumCopy (minf->value, lpval);
				EGLPNUM_TYPENAME_EGlpNumCopy (minf->objectivebound, lpval);
				EGLPNUM_TYPENAME_EGlpNumSubTo (minf->objectivebound, ILL_INTTOL);
				copy_x (lp->O->nstruct, minf->x, minf->bestx);
			}
		}

		if (EGLPNUM_TYPENAME_EGlpNumIsLeq (minf->objectivebound, lpval))
		{
			*cp = 1;
		}
		else
		{
			b = bbnodealloc (&minf->ptrworld);
			init_bbnode (b);
			b->depth = active->depth + 1;
			b->id = minf->totalnodes;
			EGLPNUM_TYPENAME_EGlpNumCopy (b->bound, lpval);
			ILL_SAFE_MALLOC (b->cstat, lp->O->nstruct, char);
			ILL_SAFE_MALLOC (b->rstat, lp->nrows, char);

			rval = EGLPNUM_TYPENAME_ILLlib_getbasis (lp, b->cstat, b->rstat);
			ILL_CLEANUP_IF (rval);
			if (minf->pinf->dII_price == QS_PRICE_DSTEEP)
			{
				b->rownorms = EGLPNUM_TYPENAME_EGlpNumAllocArray (lp->nrows);
				tval = EGLPNUM_TYPENAME_ILLlib_getrownorms (lp, minf->pinf, b->rownorms);
				if (tval)
				{
					QSlog("Row norms not available");
					QSlog("A");
					exit (1);
					EGLPNUM_TYPENAME_EGlpNumFreeArray (b->rownorms);
				}
			}
			ILL_SAFE_MALLOC (b->bound_indx, active->bound_cnt + 1, int);
			ILL_SAFE_MALLOC (b->lu, active->bound_cnt + 1, char);

			b->bounds = EGLPNUM_TYPENAME_EGlpNumAllocArray (active->bound_cnt + 1);
			for (i = 0; i < active->bound_cnt; i++)
			{
				b->bound_indx[i] = active->bound_indx[i];
				b->lu[i] = active->lu[i];
				EGLPNUM_TYPENAME_EGlpNumCopy (b->bounds[i], active->bounds[i]);
			}
			b->bound_indx[active->bound_cnt] = bvar;
			if (bdir == 'D')
				b->lu[active->bound_cnt] = 'U';
			else
				b->lu[active->bound_cnt] = 'L';
			EGLPNUM_TYPENAME_EGlpNumCopy (b->bounds[active->bound_cnt], t);
			b->bound_cnt = active->bound_cnt + 1;

			rval = EGLPNUM_TYPENAME_ILLutil_priority_insert (minf->que, (void *) b, &lpval,
																			&(b->handle));
			ILL_CLEANUP_IF (rval);

			put_bbnode (minf, b);
			minf->activenodes++;
		}
	}
	minf->totalnodes++;

	if (bdir == 'D')
	{
		rval = EGLPNUM_TYPENAME_ILLlib_chgbnd (lp, bvar, 'U', oldt);
		ILL_CLEANUP_IF (rval);
	}
	else
	{
		rval = EGLPNUM_TYPENAME_ILLlib_chgbnd (lp, bvar, 'L', oldt);
		ILL_CLEANUP_IF (rval);
	}

CLEANUP:

	EGLPNUM_TYPENAME_EGlpNumClearVar (t);
	EGLPNUM_TYPENAME_EGlpNumClearVar (lpval);
	EGLPNUM_TYPENAME_EGlpNumClearVar (oldt);
	return rval;
}

static int fix_variables (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPE * bestval,
	bbnode * b,
	EGLPNUM_TYPE * wupper,
	EGLPNUM_TYPE * wlower,
	int *hit)
{
	int rval = 0;
	int i, nnew = 0;
	int nstruct = lp->O->nstruct;
	EGLPNUM_TYPE delta, lpval;
	int *new_indx = 0;
	char *new_lu = 0;
	EGLPNUM_TYPE *new_bounds = 0;
	EGLPNUM_TYPE *dj = 0;

	EGLPNUM_TYPENAME_EGlpNumInitVar (delta);
	EGLPNUM_TYPENAME_EGlpNumInitVar (lpval);

	*hit = 0;

	if (EGLPNUM_TYPENAME_EGlpNumIsLess (*bestval, EGLPNUM_TYPENAME_ILL_MAXDOUBLE))
	{
		rval = EGLPNUM_TYPENAME_ILLlib_objval (lp, 0, &lpval);
		ILL_CLEANUP_IF (rval);
		//delta = bestval - lpval + ILL_INTTOL;
		EGLPNUM_TYPENAME_EGlpNumCopy (delta, *bestval);
		EGLPNUM_TYPENAME_EGlpNumSubTo (delta, lpval);
		EGLPNUM_TYPENAME_EGlpNumAddTo (delta, ILL_INTTOL);

		ILL_SAFE_MALLOC (new_indx, nstruct, int);
		ILL_SAFE_MALLOC (new_lu, nstruct, char);

		dj = EGLPNUM_TYPENAME_EGlpNumAllocArray (nstruct);
		new_bounds = EGLPNUM_TYPENAME_EGlpNumAllocArray (nstruct);

		rval = EGLPNUM_TYPENAME_ILLlib_solution (lp, 0, 0, 0, 0, 0, dj);
		ILL_CLEANUP_IF (rval);

		for (i = 0; i < nstruct; i++)
		{
			if (lp->O->intmarker[i])
			{
				if (EGLPNUM_TYPENAME_EGlpNumIsNeqq (wlower[i], wupper[i]))
				{
					if (EGLPNUM_TYPENAME_EGlpNumIsLess (delta, dj[i]))
					{
						EGLPNUM_TYPENAME_EGlpNumSubTo (wupper[i], EGLPNUM_TYPENAME_oneLpNum);
						rval = EGLPNUM_TYPENAME_ILLlib_chgbnd (lp, i, 'U', wupper[i]);
						ILL_CLEANUP_IF (rval);
						new_indx[nnew] = i;
						new_lu[nnew] = 'U';
						EGLPNUM_TYPENAME_EGlpNumCopy (new_bounds[nnew], wupper[i]);
						nnew++;
					}
					/*if (-dj[i] > delta) */
					EGLPNUM_TYPENAME_EGlpNumSign (delta);
					if (EGLPNUM_TYPENAME_EGlpNumIsLess (delta, dj[i]))
					{
						EGLPNUM_TYPENAME_EGlpNumAddTo (wlower[i], EGLPNUM_TYPENAME_oneLpNum);
						rval = EGLPNUM_TYPENAME_ILLlib_chgbnd (lp, i, 'L', wlower[i]);
						ILL_CLEANUP_IF (rval);
						new_indx[nnew] = i;
						new_lu[nnew] = 'L';
						EGLPNUM_TYPENAME_EGlpNumCopy (new_bounds[nnew], wlower[i]);
						nnew++;
					}
					EGLPNUM_TYPENAME_EGlpNumSign (delta);
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
			EGLPNUM_TYPENAME_EGlpNumReallocArray (&(b->bounds), b->bound_cnt + nnew);
			for (i = 0; i < nnew; i++)
			{
				b->bound_indx[b->bound_cnt + i] = new_indx[i];
				b->lu[b->bound_cnt + i] = new_lu[i];
				EGLPNUM_TYPENAME_EGlpNumCopy (b->bounds[b->bound_cnt + i], new_bounds[i]);
			}
			b->bound_cnt += nnew;
		}
	}

	*hit = nnew;

CLEANUP:

	ILL_IFFREE (new_indx, int);
	ILL_IFFREE (new_lu, char);

	EGLPNUM_TYPENAME_EGlpNumFreeArray (dj);
	EGLPNUM_TYPENAME_EGlpNumFreeArray (new_bounds);
	EGLPNUM_TYPENAME_EGlpNumClearVar (delta);
	EGLPNUM_TYPENAME_EGlpNumClearVar (lpval);
	return rval;
}

static void best_bbnode (
	mipinfo * minf,
	bbnode ** best)
{
#if 0
	bbnode *b;
	double bestval = EGLPNUM_TYPENAME_ILL_MAXDOUBLE;

	for (b = minf->head_bbnode.next; b; b = b->next)
	{
		if (b->bound < bestval)
		{
			*best = b;
			bestval = b->bound;
		}
	}
#endif

	EGLPNUM_TYPE val;

	EGLPNUM_TYPENAME_EGlpNumInitVar (val);
	EGLPNUM_TYPENAME_ILLutil_priority_deletemin (minf->que, &val, (void **) best);
	EGLPNUM_TYPENAME_EGlpNumClearVar (val);
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
	EGLPNUM_TYPE * x,
	EGLPNUM_TYPE * lpval,
	int *bvar,
	itcnt_t*itcnt)
{
	EGLPNUM_TYPENAME_lpinfo *lp = minf->lp;
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
		QSlog("Unknown branching rule.");
		rval = 1;
		goto CLEANUP;
	}

CLEANUP:

	ILL_RETURN (rval, "find_branch");
}

static void find_first_branch (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPE * x,
	int *bvar)
{
	int i, ibest = -1;
	EGLPNUM_TYPENAME_ILLlpdata *qslp = lp->O;
	EGLPNUM_TYPE t;

	EGLPNUM_TYPENAME_EGlpNumInitVar (t);

	for (i = 0; i < qslp->nstruct; i++)
	{
		if (qslp->intmarker[i])
		{
			/*t = EGLPNUM_TYPENAME_ILLutil_our_frac (x[i]); */
			EGLPNUM_TYPENAME_EGlpNumFloor (t, x[i]);
			EGLPNUM_TYPENAME_EGlpNumSubTo (t, x[i]);
			EGLPNUM_TYPENAME_EGlpNumSign (t);
			if ((EGLPNUM_TYPENAME_EGlpNumIsNeqZero (t, ILL_INTTOL)) &&
					(EGLPNUM_TYPENAME_EGlpNumIsNeq (t, EGLPNUM_TYPENAME_oneLpNum, ILL_INTTOL)))
			{
				ibest = i;
				break;
			}
		}
	}
	*bvar = ibest;
	EGLPNUM_TYPENAME_EGlpNumClearVar (t);
}

static void find_middle_branch (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPE * x,
	int *bvar)
{
	int i, ibest = -1;
	EGLPNUM_TYPE t, tbest;
	EGLPNUM_TYPENAME_ILLlpdata *qlp = lp->O;

	EGLPNUM_TYPENAME_EGlpNumInitVar (t);
	EGLPNUM_TYPENAME_EGlpNumInitVar (tbest);
	EGLPNUM_TYPENAME_EGlpNumSet (tbest, 0.5);

	for (i = 0; i < qlp->nstruct; i++)
	{
		if (qlp->intmarker[i])
		{
			/*t = EGLPNUM_TYPENAME_ILLutil_our_frac (x[i]) - 0.5;
			 * if (t < 0.0)
			 * t = -t; */
			EGLPNUM_TYPENAME_EGlpNumFloor (t, x[i]);
			EGLPNUM_TYPENAME_EGlpNumMultUiTo (t, 2);
			EGLPNUM_TYPENAME_EGlpNumSubTo (t, EGLPNUM_TYPENAME_oneLpNum);
			EGLPNUM_TYPENAME_EGlpNumDivUiTo (t, 2);
			if (EGLPNUM_TYPENAME_EGlpNumIsLessZero (t))
				EGLPNUM_TYPENAME_EGlpNumSign (t);
			/*if (t < tbest) */
			if (EGLPNUM_TYPENAME_EGlpNumIsLess (t, tbest))
			{
				EGLPNUM_TYPENAME_EGlpNumCopy (tbest, t);
				ibest = i;
			}
		}
	}

	/*if (tbest < (0.5 - ILL_INTTOL)) */
	EGLPNUM_TYPENAME_EGlpNumAddTo (tbest, ILL_INTTOL);
	if (EGLPNUM_TYPENAME_EGlpNumIsLessDbl (tbest, 0.5))
	{
		*bvar = ibest;
	}
	else
	{
		*bvar = -1;
	}
	EGLPNUM_TYPENAME_EGlpNumClearVar (t);
	EGLPNUM_TYPENAME_EGlpNumClearVar (tbest);
}

static int find_penalty_branch (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_price_info * pinf,
	EGLPNUM_TYPE * x,
	EGLPNUM_TYPE * downpen,
	EGLPNUM_TYPE * uppen,
	EGLPNUM_TYPE * lpval,
	int *bvar,
	itcnt_t*itcnt)
{
	int rval = 0;
	int i, k, ibest = -1, ncand = 0, nneed = 0;
	EGLPNUM_TYPENAME_ILLlpdata *qslp = lp->O;
	int *candidatelist = 0;
	int *needlist = 0;
	EGLPNUM_TYPE *fval = 0;
	EGLPNUM_TYPE *xlist = 0;
	EGLPNUM_TYPE *newdown = 0;
	EGLPNUM_TYPE *newup = 0;
	EGLPNUM_TYPE a, t, tbest;

	EGLPNUM_TYPENAME_EGlpNumInitVar (a);
	EGLPNUM_TYPENAME_EGlpNumInitVar (t);
	EGLPNUM_TYPENAME_EGlpNumInitVar (tbest);
	EGLPNUM_TYPENAME_EGlpNumCopy (tbest, EGLPNUM_TYPENAME_ILL_MINDOUBLE);

	ILL_SAFE_MALLOC (candidatelist, qslp->nstruct, int);
	ILL_SAFE_MALLOC (needlist, qslp->nstruct, int);

	fval = EGLPNUM_TYPENAME_EGlpNumAllocArray (qslp->nstruct);
	xlist = EGLPNUM_TYPENAME_EGlpNumAllocArray (qslp->nstruct);
	for (i = 0; i < qslp->nstruct; i++)
	{
		if (qslp->intmarker[i])
		{
			/*fval[i] = x[i] - floor(x[i]); */
			EGLPNUM_TYPENAME_EGlpNumFloor (fval[i], x[i]);
			EGLPNUM_TYPENAME_EGlpNumSubTo (fval[i], x[i]);
			EGLPNUM_TYPENAME_EGlpNumSign (fval[i]);
			if ((EGLPNUM_TYPENAME_EGlpNumIsNeqZero (fval[i], ILL_INTTOL)) &&
					(EGLPNUM_TYPENAME_EGlpNumIsNeq (fval[i], EGLPNUM_TYPENAME_oneLpNum, ILL_INTTOL)))
			{
				candidatelist[ncand++] = i;
				/*if (downpen[i] == -1.0) */
				EGLPNUM_TYPENAME_EGlpNumSign (downpen[i]);
				if (EGLPNUM_TYPENAME_EGlpNumIsEqqual (downpen[i], EGLPNUM_TYPENAME_oneLpNum))
				{
					EGLPNUM_TYPENAME_EGlpNumCopy (xlist[nneed], x[i]);
					needlist[nneed++] = i;
				}
				EGLPNUM_TYPENAME_EGlpNumSign (downpen[i]);
			}
		}
	}

	if (nneed > 0)
	{
		newdown = EGLPNUM_TYPENAME_EGlpNumAllocArray (nneed);
		newup = EGLPNUM_TYPENAME_EGlpNumAllocArray (nneed);
		rval = EGLPNUM_TYPENAME_ILLlib_strongbranch (lp, pinf, needlist, nneed,
																0, newdown, newup,
																5 * STRONG_PIVOTS, EGLPNUM_TYPENAME_ILL_MAXDOUBLE, itcnt);
		ILL_CLEANUP_IF (rval);

		for (i = 0; i < nneed; i++)
		{
			k = needlist[i];
			/*uppen[k] = (newup[i] - lpval) / (1.0 - fval[k]); */
			EGLPNUM_TYPENAME_EGlpNumCopyDiff (uppen[k], newup[i], *lpval);
			EGLPNUM_TYPENAME_EGlpNumCopyDiff (downpen[k], EGLPNUM_TYPENAME_oneLpNum, fval[k]);
			EGLPNUM_TYPENAME_EGlpNumDivTo (uppen[k], downpen[k]);
			/*downpen[k] = (newdown[i] - lpval) / fval[k]; */
			EGLPNUM_TYPENAME_EGlpNumCopyDiffRatio (downpen[k], newdown[i], *lpval, fval[k]);

		}
	}

	for (i = 0; i < ncand; i++)
	{
		k = candidatelist[i];
		/*t = ILL_BRANCH_PENALTY_VAL (downpen[k], uppen[k], fval[k]); */
		EGLPNUM_TYPENAME_EGlpNumCopy (t, downpen[k]);
		EGLPNUM_TYPENAME_EGlpNumMultTo (t, fval[k]);
		EGLPNUM_TYPENAME_EGlpNumCopyDiff (a, EGLPNUM_TYPENAME_oneLpNum, fval[k]);
		EGLPNUM_TYPENAME_EGlpNumMultTo (a, uppen[k]);
		if (EGLPNUM_TYPENAME_EGlpNumIsLess (t, a))
		{
			EGLPNUM_TYPENAME_EGlpNumMultUiTo (t, ILL_BRANCH_PENALTY_WEIGHT);
			EGLPNUM_TYPENAME_EGlpNumAddTo (t, a);
		}
		else
		{
			EGLPNUM_TYPENAME_EGlpNumMultUiTo (a, ILL_BRANCH_PENALTY_WEIGHT);
			EGLPNUM_TYPENAME_EGlpNumAddTo (t, a);
		}
		EGLPNUM_TYPENAME_EGlpNumDivUiTo (t, ILL_BRANCH_PENALTY_WEIGHT + 1);

		if (EGLPNUM_TYPENAME_EGlpNumIsLess (tbest, t))
		{
			EGLPNUM_TYPENAME_EGlpNumCopy (tbest, t);
			ibest = k;
		}
	}

	*bvar = ibest;

CLEANUP:

	EGLPNUM_TYPENAME_EGlpNumClearVar (a);
	EGLPNUM_TYPENAME_EGlpNumClearVar (t);
	EGLPNUM_TYPENAME_EGlpNumClearVar (tbest);
	EGLPNUM_TYPENAME_EGlpNumFreeArray (newdown);
	EGLPNUM_TYPENAME_EGlpNumFreeArray (newup);
	EGLPNUM_TYPENAME_EGlpNumFreeArray (fval);
	EGLPNUM_TYPENAME_EGlpNumFreeArray (xlist);
	ILL_IFFREE (candidatelist, int);
	ILL_IFFREE (needlist, int);

	return rval;
}

static int find_strong_branch (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_price_info * pinf,
	EGLPNUM_TYPE * x,
	int *bvar,
	itcnt_t*itcnt)
{
	int rval = 0;
	int i, ibest = -1, ncand = 0;
	int maxtrys = STRONG_CANDIDATES;
	EGLPNUM_TYPE t, tbest;
	EGLPNUM_TYPENAME_ILLlpdata *qlp = lp->O;
	int *candidatelist = 0;
	int *newlist = 0;
	int *perm = 0;
	EGLPNUM_TYPE *tval = 0;
	EGLPNUM_TYPE *xlist = 0;
	EGLPNUM_TYPE *downpen = 0;
	EGLPNUM_TYPE *uppen = 0;
	ILLrandstate rstate;

	EGLPNUM_TYPENAME_EGlpNumInitVar (t);
	EGLPNUM_TYPENAME_EGlpNumInitVar (tbest);
	EGLPNUM_TYPENAME_EGlpNumCopy (tbest, EGLPNUM_TYPENAME_ILL_MINDOUBLE);

	ILLutil_sprand (999, &rstate);
	ILL_SAFE_MALLOC (candidatelist, qlp->nstruct, int);

	tval = EGLPNUM_TYPENAME_EGlpNumAllocArray (qlp->nstruct);

	for (i = 0; i < qlp->nstruct; i++)
	{
		if (qlp->intmarker[i])
		{
			/*t = EGLPNUM_TYPENAME_ILLutil_our_frac (x[i]) - 0.5;
			 * if (t < 0.0)
			 * t = -t; */
			EGLPNUM_TYPENAME_EGlpNumFloor (t, x[i]);
			EGLPNUM_TYPENAME_EGlpNumSubTo (t, x[i]);
			EGLPNUM_TYPENAME_EGlpNumSign (t);
			EGLPNUM_TYPENAME_EGlpNumMultUiTo (t, 2);
			EGLPNUM_TYPENAME_EGlpNumSubTo (t, EGLPNUM_TYPENAME_oneLpNum);
			if (EGLPNUM_TYPENAME_EGlpNumIsLessZero (t))
				EGLPNUM_TYPENAME_EGlpNumSign (t);
			/*if (t < (0.5 - ILL_INTTOL)) */
			if (EGLPNUM_TYPENAME_EGlpNumIsNeq (t, EGLPNUM_TYPENAME_oneLpNum, ILL_INTTOL))
			{
				candidatelist[ncand] = i;
				EGLPNUM_TYPENAME_EGlpNumDivUiTo (t, 2);
				EGLPNUM_TYPENAME_EGlpNumCopy (tval[ncand++], t);
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
			EGLPNUM_TYPENAME_ILLutil_EGlpNum_rselect (perm, 0, ncand - 1, maxtrys, tval, &rstate);

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

		downpen = EGLPNUM_TYPENAME_EGlpNumAllocArray (ncand);
		uppen = EGLPNUM_TYPENAME_EGlpNumAllocArray (ncand);
		xlist = EGLPNUM_TYPENAME_EGlpNumAllocArray (ncand);

		for (i = 0; i < ncand; i++)
		{
			EGLPNUM_TYPENAME_EGlpNumCopy (xlist[i], x[candidatelist[i]]);
		}

		rval = EGLPNUM_TYPENAME_ILLlib_strongbranch (lp, pinf, candidatelist, ncand,
																0, downpen, uppen, STRONG_PIVOTS,
																EGLPNUM_TYPENAME_ILL_MAXDOUBLE, itcnt);
		ILL_CLEANUP_IF (rval);

		for (i = 0; i < ncand; i++)
		{
			/*t = ILL_BRANCH_STRONG_VAL (downpen[i], uppen[i]); */
			if (EGLPNUM_TYPENAME_EGlpNumIsLess (downpen[i], uppen[i]))
			{
				EGLPNUM_TYPENAME_EGlpNumCopy (t, downpen[i]);
				EGLPNUM_TYPENAME_EGlpNumMultUiTo (t, ILL_BRANCH_STRONG_WEIGHT);
				EGLPNUM_TYPENAME_EGlpNumAddTo (t, uppen[i]);
			}
			else
			{
				EGLPNUM_TYPENAME_EGlpNumCopy (t, uppen[i]);
				EGLPNUM_TYPENAME_EGlpNumMultUiTo (t, ILL_BRANCH_STRONG_WEIGHT);
				EGLPNUM_TYPENAME_EGlpNumAddTo (t, downpen[i]);
			}
			EGLPNUM_TYPENAME_EGlpNumDivUiTo (t, ILL_BRANCH_STRONG_WEIGHT + 1);
			if (EGLPNUM_TYPENAME_EGlpNumIsLess (tbest, t))
			{
				EGLPNUM_TYPENAME_EGlpNumCopy (tbest, t);
				ibest = candidatelist[i];
			}
		}
	}

	*bvar = ibest;


CLEANUP:

	EGLPNUM_TYPENAME_EGlpNumClearVar (t);
	EGLPNUM_TYPENAME_EGlpNumClearVar (tbest);
	EGLPNUM_TYPENAME_EGlpNumFreeArray (tval);
	EGLPNUM_TYPENAME_EGlpNumFreeArray (xlist);
	EGLPNUM_TYPENAME_EGlpNumFreeArray (uppen);
	EGLPNUM_TYPENAME_EGlpNumFreeArray (downpen);
	ILL_IFFREE (candidatelist, int);
	ILL_IFFREE (newlist, int);
	ILL_IFFREE (perm, int);

	ILL_RETURN (rval, "find_strong_branch");
}

static void check_integral (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPE * x,
	int *yesno)
{
	int i;
	EGLPNUM_TYPE t;
	EGLPNUM_TYPENAME_ILLlpdata *qlp = lp->O;

	EGLPNUM_TYPENAME_EGlpNumInitVar (t);

	for (i = 0; i < qlp->nstruct; i++)
	{
		if (qlp->intmarker[i])
		{
			/*t = EGLPNUM_TYPENAME_ILLutil_our_frac (x[i]); */
			EGLPNUM_TYPENAME_EGlpNumFloor (t, x[i]);
			EGLPNUM_TYPENAME_EGlpNumSubTo (t, x[i]);
			EGLPNUM_TYPENAME_EGlpNumSign (t);
			/*if (t > ILL_INTTOL && t < 1.0 - ILL_INTTOL) */
			if ((EGLPNUM_TYPENAME_EGlpNumIsNeqZero (t, ILL_INTTOL)) &&
					(EGLPNUM_TYPENAME_EGlpNumIsNeq (t, EGLPNUM_TYPENAME_oneLpNum, ILL_INTTOL)))
			{
				*yesno = 0;
				EGLPNUM_TYPENAME_EGlpNumClearVar (t);
				return;
			}
		}
	}

	*yesno = 1;
	EGLPNUM_TYPENAME_EGlpNumClearVar (t);
}

static int plunge (
	mipinfo * minf,
	itcnt_t*itcnt)
{
	int rval = 0;
	int i, status;
	EGLPNUM_TYPENAME_lpinfo *lp = minf->lp;
	EGLPNUM_TYPENAME_ILLlpdata *qlp = minf->lp->O;
	EGLPNUM_TYPE *oldlower = 0;
	EGLPNUM_TYPE *oldupper = 0;

	if (minf->watch)
	{
		QSlog("Plunging ...");
	}

	oldlower = EGLPNUM_TYPENAME_EGlpNumAllocArray (qlp->nstruct);
	oldupper = EGLPNUM_TYPENAME_EGlpNumAllocArray (qlp->nstruct);

	for (i = 0; i < qlp->nstruct; i++)
	{
		EGLPNUM_TYPENAME_EGlpNumCopy (oldlower[i], minf->lower[i]);
		EGLPNUM_TYPENAME_EGlpNumCopy (oldupper[i], minf->upper[i]);
	}

	rval = plunge_work (minf, 0, itcnt);
	ILL_CLEANUP_IF (rval);

	for (i = 0; i < qlp->nstruct; i++)
	{
		rval = EGLPNUM_TYPENAME_ILLlib_chgbnd (lp, i, 'L', oldlower[i]);
		ILL_CLEANUP_IF (rval);
		rval = EGLPNUM_TYPENAME_ILLlib_chgbnd (lp, i, 'U', oldupper[i]);
		ILL_CLEANUP_IF (rval);
		EGLPNUM_TYPENAME_EGlpNumCopy (minf->lower[i], oldlower[i]);
		EGLPNUM_TYPENAME_EGlpNumCopy (minf->upper[i], oldupper[i]);
	}

	rval = EGLPNUM_TYPENAME_ILLlib_optimize (lp, 0, minf->pinf, DUAL_SIMPLEX, &status, 0, itcnt);
	ILL_CLEANUP_IF (rval);


CLEANUP:

	EGLPNUM_TYPENAME_EGlpNumFreeArray (oldlower);
	EGLPNUM_TYPENAME_EGlpNumFreeArray (oldupper);

	ILL_RETURN (rval, "plunge");
}

static int plunge_work (
	mipinfo * minf,
	int depth,
	itcnt_t*itcnt)
{
	int rval = 0;
	int bvar, status, count;
	EGLPNUM_TYPE lpval, val0, val1, int_tol;
	EGLPNUM_TYPENAME_lpinfo *lp = minf->lp;

	EGLPNUM_TYPENAME_EGlpNumInitVar (lpval);
	EGLPNUM_TYPENAME_EGlpNumInitVar (val0);
	EGLPNUM_TYPENAME_EGlpNumInitVar (val1);
	EGLPNUM_TYPENAME_EGlpNumInitVar (int_tol);
	EGLPNUM_TYPENAME_EGlpNumSet (int_tol, 0.001);

	rval = EGLPNUM_TYPENAME_ILLlib_get_x (lp, 0, minf->x);
	ILL_CLEANUP_IF (rval);

	rval = round_variables (minf, &count, &int_tol /* 0.001 */ );
	ILL_CLEANUP_IF (rval);
	if (count)
	{
		rval = EGLPNUM_TYPENAME_ILLlib_optimize (lp, 0, minf->pinf, DUAL_SIMPLEX, &status, 0, itcnt);
		ILL_CLEANUP_IF (rval);
		if (status != QS_LP_OPTIMAL)
		{
			goto CLEANUP;
		}
		rval = EGLPNUM_TYPENAME_ILLlib_get_x (lp, 0, minf->x);
		ILL_CLEANUP_IF (rval);
	}

	find_middle_branch (lp, minf->x, &bvar);
	if (bvar == -1)
	{
		rval = EGLPNUM_TYPENAME_ILLlib_objval (lp, 0, &lpval);
		ILL_CLEANUP_IF (rval);

		if (EGLPNUM_TYPENAME_EGlpNumIsLess (lpval, minf->value))
		{
			QSlog("Plunge Integral Solution: %.6f (Depth: %d)",
									EGLPNUM_TYPENAME_EGlpNumToLf (lpval), depth);

			EGLPNUM_TYPENAME_EGlpNumCopy (minf->value, lpval);
			EGLPNUM_TYPENAME_EGlpNumCopyDiff (minf->objectivebound, lpval, ILL_INTTOL);
			copy_x (lp->O->nstruct, minf->x, minf->bestx);
		}
		goto CLEANUP;
	}

	EGLPNUM_TYPENAME_EGlpNumOne (minf->lower[bvar]);
	rval = EGLPNUM_TYPENAME_ILLlib_chgbnd (lp, bvar, 'L', EGLPNUM_TYPENAME_oneLpNum);
	ILL_CLEANUP_IF (rval);
	rval = EGLPNUM_TYPENAME_ILLlib_optimize (lp, 0, minf->pinf, DUAL_SIMPLEX, &status, 0, itcnt);
	ILL_CLEANUP_IF (rval);

	if (status == QS_LP_UNSOLVED)
	{
		QSlog("Simplex did not solve the plunge LP");
		rval = 1;
		ILL_CLEANUP;
	}
	else if (status == QS_LP_INFEASIBLE)
	{
		EGLPNUM_TYPENAME_EGlpNumCopy (val1, EGLPNUM_TYPENAME_ILL_MAXDOUBLE);
	}
	else if (status == QS_LP_OPTIMAL)
	{
		rval = EGLPNUM_TYPENAME_ILLlib_objval (lp, 0, &val1);
		ILL_CLEANUP_IF (rval);
	}
	else
	{
		ILL_CLEANUP;
	}

	rval = EGLPNUM_TYPENAME_ILLlib_chgbnd (lp, bvar, 'L', EGLPNUM_TYPENAME_zeroLpNum);
	ILL_CLEANUP_IF (rval);
	EGLPNUM_TYPENAME_EGlpNumZero (minf->lower[bvar]);

	EGLPNUM_TYPENAME_EGlpNumZero (minf->upper[bvar]);
	rval = EGLPNUM_TYPENAME_ILLlib_chgbnd (lp, bvar, 'U', EGLPNUM_TYPENAME_zeroLpNum);
	ILL_CLEANUP_IF (rval);
	rval = EGLPNUM_TYPENAME_ILLlib_optimize (lp, 0, minf->pinf, DUAL_SIMPLEX, &status, 0, itcnt);
	ILL_CLEANUP_IF (rval);

	if (status == QS_LP_UNSOLVED)
	{
		QSlog("Simplex did not solve the plunge LP");
		rval = 1;
		ILL_CLEANUP;
	}
	else if (status == QS_LP_INFEASIBLE)
	{
		EGLPNUM_TYPENAME_EGlpNumCopy (val0, EGLPNUM_TYPENAME_ILL_MAXDOUBLE);
	}
	else if (status == QS_LP_OPTIMAL)
	{
		rval = EGLPNUM_TYPENAME_ILLlib_objval (lp, 0, &val0);
		ILL_CLEANUP_IF (rval);
	}
	else
	{
		ILL_CLEANUP;
	}

	rval = EGLPNUM_TYPENAME_ILLlib_chgbnd (lp, bvar, 'U', EGLPNUM_TYPENAME_oneLpNum);
	ILL_CLEANUP_IF (rval);
	EGLPNUM_TYPENAME_EGlpNumCopy (minf->upper[bvar], EGLPNUM_TYPENAME_oneLpNum);

	if (EGLPNUM_TYPENAME_EGlpNumIsEqqual (val0, EGLPNUM_TYPENAME_ILL_MAXDOUBLE) &&
			EGLPNUM_TYPENAME_EGlpNumIsEqqual (val1, EGLPNUM_TYPENAME_ILL_MAXDOUBLE))
	{
		ILL_CLEANUP;
	}

	if (EGLPNUM_TYPENAME_EGlpNumIsLess (val0, val1))
	{
		EGLPNUM_TYPENAME_EGlpNumZero (minf->upper[bvar]);
		rval = EGLPNUM_TYPENAME_ILLlib_chgbnd (lp, bvar, 'U', EGLPNUM_TYPENAME_zeroLpNum);
		ILL_CLEANUP_IF (rval);
		rval = EGLPNUM_TYPENAME_ILLlib_optimize (lp, 0, minf->pinf, DUAL_SIMPLEX, &status, 0, itcnt);
		ILL_CLEANUP_IF (rval);
		rval = plunge_work (minf, depth + 1, itcnt);
		ILL_CLEANUP_IF (rval);
		rval = EGLPNUM_TYPENAME_ILLlib_chgbnd (lp, bvar, 'U', EGLPNUM_TYPENAME_oneLpNum);
		ILL_CLEANUP_IF (rval);
		EGLPNUM_TYPENAME_EGlpNumOne (minf->upper[bvar]);
	}
	else
	{
		EGLPNUM_TYPENAME_EGlpNumOne (minf->lower[bvar]);
		rval = EGLPNUM_TYPENAME_ILLlib_chgbnd (lp, bvar, 'L', EGLPNUM_TYPENAME_oneLpNum);
		ILL_CLEANUP_IF (rval);
		rval = EGLPNUM_TYPENAME_ILLlib_optimize (lp, 0, minf->pinf, DUAL_SIMPLEX, &status, 0, itcnt);
		ILL_CLEANUP_IF (rval);
		rval = plunge_work (minf, depth + 1, itcnt);
		ILL_CLEANUP_IF (rval);
		rval = EGLPNUM_TYPENAME_ILLlib_chgbnd (lp, bvar, 'L', EGLPNUM_TYPENAME_zeroLpNum);
		ILL_CLEANUP_IF (rval);
		EGLPNUM_TYPENAME_EGlpNumZero (minf->lower[bvar]);
	}

CLEANUP:

	EGLPNUM_TYPENAME_EGlpNumClearVar (lpval);
	EGLPNUM_TYPENAME_EGlpNumClearVar (val0);
	EGLPNUM_TYPENAME_EGlpNumClearVar (val1);
	EGLPNUM_TYPENAME_EGlpNumClearVar (int_tol);
	ILL_RETURN (rval, "plunge_work");
}

static int round_variables (
	mipinfo * minf,
	int *count,
	EGLPNUM_TYPE * tol)
{
	int rval = 0;
	int i, hit = 0;
	EGLPNUM_TYPENAME_lpinfo *lp = minf->lp;
	EGLPNUM_TYPENAME_ILLlpdata *qlp = lp->O;

	*count = 0;

	for (i = 0; i < qlp->nstruct; i++)
	{
		if (qlp->intmarker[i])
		{
			if (EGLPNUM_TYPENAME_EGlpNumIsNeqq (minf->lower[i], minf->upper[i]))
			{
				if (EGLPNUM_TYPENAME_EGlpNumIsLess (minf->x[i], *tol))
				{
					EGLPNUM_TYPENAME_EGlpNumZero (minf->upper[i]);
					rval = EGLPNUM_TYPENAME_ILLlib_chgbnd (lp, i, 'U', EGLPNUM_TYPENAME_zeroLpNum);
					ILL_CLEANUP_IF (rval);
					hit++;
				}
				else if (EGLPNUM_TYPENAME_EGlpNumIsEqual (minf->x[i], EGLPNUM_TYPENAME_oneLpNum, *tol))
				{
					EGLPNUM_TYPENAME_EGlpNumOne (minf->lower[i]);
					rval = EGLPNUM_TYPENAME_ILLlib_chgbnd (lp, i, 'L', EGLPNUM_TYPENAME_oneLpNum);
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
	EGLPNUM_TYPE * from_x,
	EGLPNUM_TYPE * to_x)
{
	int j;

	for (j = 0; j < nstruct; j++)
	{
		EGLPNUM_TYPENAME_EGlpNumCopy (to_x[j], from_x[j]);
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
		EGLPNUM_TYPENAME_EGlpNumInitVar (minf->objectivebound);
		EGLPNUM_TYPENAME_EGlpNumInitVar (minf->value);
		EGLPNUM_TYPENAME_EGlpNumCopy (minf->objectivebound, EGLPNUM_TYPENAME_ILL_MAXDOUBLE);
		EGLPNUM_TYPENAME_EGlpNumCopy (minf->value, EGLPNUM_TYPENAME_ILL_MAXDOUBLE);
		ILLptrworld_init (&minf->ptrworld);
	}
}

static void free_mipinfo (
	mipinfo * minf)
{
	int total, onlist;

	if (minf)
	{
		EGLPNUM_TYPENAME_EGlpNumFreeArray (minf->downpen);
		EGLPNUM_TYPENAME_EGlpNumFreeArray (minf->uppen);
		EGLPNUM_TYPENAME_EGlpNumFreeArray (minf->x);
		EGLPNUM_TYPENAME_EGlpNumFreeArray (minf->bestx);
		EGLPNUM_TYPENAME_EGlpNumFreeArray (minf->lower);
		EGLPNUM_TYPENAME_EGlpNumFreeArray (minf->upper);
		bbnode_listfree (&minf->ptrworld, minf->head_bbnode.next);
		if (bbnode_check_leaks (&minf->ptrworld, &total, &onlist))
		{
			QSlog("WARNING: %d outstanding bbnodes", total - onlist);
		}
		ILLptrworld_delete (&minf->ptrworld);
		EGLPNUM_TYPENAME_EGlpNumClearVar ((minf->objectivebound));
		EGLPNUM_TYPENAME_EGlpNumClearVar ((minf->value));
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
		EGLPNUM_TYPENAME_EGlpNumInitVar ((b->bound));
		EGLPNUM_TYPENAME_EGlpNumCopy (b->bound, EGLPNUM_TYPENAME_ILL_MINDOUBLE);
	}
}

static void free_bbnode (
	bbnode * b)
{
	if (b)
	{
		EGLPNUM_TYPENAME_EGlpNumFreeArray (b->rownorms);
		EGLPNUM_TYPENAME_EGlpNumFreeArray (b->bounds);
		ILL_IFFREE (b->cstat, char);
		ILL_IFFREE (b->rstat, char);
		ILL_IFFREE (b->bound_indx, int);
		ILL_IFFREE (b->lu, char);

		EGLPNUM_TYPENAME_EGlpNumClearVar ((b->bound));
		memset (b, 0, sizeof (bbnode));
	}
}
