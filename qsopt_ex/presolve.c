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

/* RCS_INFO = "$RCSfile: presolve.c,v $ $Revision: 1.2 $ $Date: 2003/11/05 16:49:52 $"; */
static int TRACE = 0;

/****************************************************************************/
/*                                                                          */
/*                 Presolve Routine for Simplex Method                      */
/*                                                                          */
/*  EXPORTED FUNCTIONS                                                      */
/*                                                                          */
/*    int EGLPNUM_TYPENAME_ILLlp_add_logicals (ILLlpata *lp)                                 */
/*    int EGLPNUM_TYPENAME_ILLlp_presolve (EGLPNUM_TYPENAME_ILLlpdata *lp)                                    */
/*    int EGLPNUM_TYPENAME_ILLlp_scale (EGLPNUM_TYPENAME_ILLlpdata *lp)                                       */
/*    void EGLPNUM_TYPENAME_ILLlp_sinfo_init (EGLPNUM_TYPENAME_ILLlp_sinfo *sinfo)                            */
/*    void EGLPNUM_TYPENAME_ILLlp_sinfo_free (EGLPNUM_TYPENAME_ILLlp_sinfo *sinfo)                            */
/*    void EGLPNUM_TYPENAME_ILLlp_predata_init (EGLPNUM_TYPENAME_ILLlp_predata *pre)                          */
/*    void EGLPNUM_TYPENAME_ILLlp_predata_free (EGLPNUM_TYPENAME_ILLlp_predata *pre)                          */
/*                                                                          */
/*  NOTES                                                                   */
/*                                                                          */
/*    presolve will assume that logicals have been added.                   */
/*                                                                          */
/****************************************************************************/

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "presolve_EGLPNUM_TYPENAME.h"

#include <stdlib.h>
#include <string.h>

#include "qs_config.h"
#include "logging-private.h"

#include "allocrus.h"
#include "eg_lpnum.h"
#include "eg_io.h"
#include "except.h"
#include "zeit.h"

#include "lpdata_EGLPNUM_TYPENAME.h"
#include "lpdefs_EGLPNUM_TYPENAME.h"
//extern EGLPNUM_TYPE EGLPNUM_TYPENAME_SZERO_TOLER;

#define ILL_LP_STATUS_OK (0)
#define ILL_PRE_FEAS_TOL EGLPNUM_TYPENAME_PFEAS_TOLER	//(1e-6)
#define ILL_PRE_ZERO_TOL EGLPNUM_TYPENAME_PIVOT_TOLER	//(1e-10)

#define ILL_PRE_DELETE_EMPTY_ROW               (1)
#define ILL_PRE_DELETE_SINGLETON_ROW           (2)
#define ILL_PRE_DELETE_FIXED_VARIABLE          (3)
#define ILL_PRE_DELETE_FORCED_VARIABLE         (4)
#define ILL_PRE_DELETE_SINGLETON_VARIABLE      (5)
#define ILL_PRE_DELETE_FREE_SINGLETON_VARIABLE (6)
#define ILL_PRE_DELETE_EMPTY_COLUMN            (7)

#define ILL_PRE_COL_STRUC                      (0)
#define ILL_PRE_COL_LOGICAL                    (1)

static int debug = 0;

typedef struct
{
	int row;
	int col;
	char coltype;
	char mark;
	char del;
	EGLPNUM_TYPE coef;
}
edge;

typedef struct node
{
	edge **adj;
	EGLPNUM_TYPE obj;
	EGLPNUM_TYPE lower;
	EGLPNUM_TYPE upper;
	EGLPNUM_TYPE rhs;
	int deg;
	char mark;
	char del;
	char coltype;
	char rowsense;
}
node;

typedef struct intptr
{
	int this_val;
	struct intptr *next;
}
intptr;

typedef struct graph
{
	edge *edgelist;
	struct node *rows;
	struct node *cols;
	int ecount;
	int nrows;
	int ncols;
	int nzcount;
	edge **adjspace;
	ILLptrworld intptrworld;
	int objsense;
}
graph;

static void set_fixed_variable (
	graph * G,
	int j,
	EGLPNUM_TYPE val),
  get_implied_rhs_bounds (
	graph * G,
	int i,
	EGLPNUM_TYPE * lb,
	EGLPNUM_TYPE * ub),
  get_implied_variable_bounds (
	graph * G,
	int j,
	edge * a_ij,
	EGLPNUM_TYPE * lb,
	EGLPNUM_TYPE * ub),
  dump_line (
	EGLPNUM_TYPENAME_ILLlp_preline * line),
  init_graph (
	graph * G),
  free_graph (
	graph * G),
  dump_graph (
	graph * G);

static int simple_presolve (
	EGLPNUM_TYPENAME_ILLlpdata * lp,
	EGLPNUM_TYPENAME_ILLlp_predata * pre,
	EGLPNUM_TYPENAME_ILLlp_sinfo * info,
	int pre_types,
	int *status),
  grab_lp_line (
	graph * G,
	int indx,
	EGLPNUM_TYPENAME_ILLlp_preline * line,
	int row_or_col),
  grab_lp_info (
	graph * G,
	char **colnames,
	EGLPNUM_TYPENAME_ILLlp_sinfo * info),
  fixed_variables (
	graph * G,
	EGLPNUM_TYPENAME_ILLlp_predata * pre),
  empty_columns (
	graph * G,
	EGLPNUM_TYPENAME_ILLlp_predata * pre),
  singleton_rows (
	graph * G,
	EGLPNUM_TYPENAME_ILLlp_predata * pre,
	int *hit),
  forcing_constraints (
	graph * G,
	EGLPNUM_TYPENAME_ILLlp_predata * pre,
	int *hit),
  singleton_columns (
	graph * G,
	EGLPNUM_TYPENAME_ILLlp_predata * pre,
	int *hit),
  duplicate_rows (
	graph * G,
	int *hit),
  duplicate_cols (
	graph * G,
	int *hit),
  gather_dup_lists (
	int *s,
	int count,
	int *duptotal,
	int **dupcnt,
	int **dupind),
  get_next_preop (
	EGLPNUM_TYPENAME_ILLlp_predata * pre,
	EGLPNUM_TYPENAME_ILLlp_preop ** op),
  add_to_list (
	ILLptrworld * world,
	intptr ** list,
	int i),
  build_graph (
	EGLPNUM_TYPENAME_ILLlpdata * lp,
	graph * G);


ILL_PTRWORLD_ROUTINES (intptr, intptralloc, intptr_bulkalloc, intptrfree)
ILL_PTRWORLD_LISTFREE_ROUTINE (intptr, intptr_listfree, intptrfree)
ILL_PTRWORLD_LEAKS_ROUTINE (intptr, intptr_check_leaks, this_val, int)
		 int EGLPNUM_TYPENAME_ILLlp_add_logicals (
	EGLPNUM_TYPENAME_ILLlpdata * lp)
{
	int rval = 0;
	int ncols, nrows, nzcount, i, aindex;
	char *sense;
	EGLPNUM_TYPENAME_ILLmatrix *A;

	if (!lp)
	{
		QSlog("EGLPNUM_TYPENAME_ILLlp_add_logicals called with a NULL pointer");
		rval = 1;
		goto CLEANUP;
	}

	QSlog("EGLPNUM_TYPENAME_ILLlp_add_logicals ...");

	A = &lp->A;
	sense = lp->sense;
	ncols = lp->ncols;
	nrows = lp->nrows;
	nzcount = lp->nzcount;

	if (nrows == 0)
		goto CLEANUP;
	EGLPNUM_TYPENAME_EGlpNumReallocArray (&(lp->obj), lp->colsize + nrows);
	EGLPNUM_TYPENAME_EGlpNumReallocArray (&(lp->upper), lp->colsize + nrows);
	EGLPNUM_TYPENAME_EGlpNumReallocArray (&(lp->lower), lp->colsize + nrows);
	lp->colnames =
		EGrealloc (lp->colnames, sizeof (char *) * (lp->colsize + nrows));
	//rval = ILLutil_reallocrus_count ((void **) &(lp->colnames),
	//                                 lp->colsize + nrows, sizeof (char *));
	//ILL_CLEANUP_IF (rval);
	memset (lp->colnames + ncols, 0, sizeof (char *) * nrows);

	ILL_SAFE_MALLOC (lp->rowmap, lp->rowsize, int);


	A->matcnt = EGrealloc (A->matcnt, sizeof (int) * (A->matcolsize + nrows));
	//rval = ILLutil_reallocrus_count ((void **) &(A->matcnt),
	//                                 A->matcolsize + nrows, sizeof (int));
	//ILL_CLEANUP_IF (rval);

	A->matbeg = EGrealloc (A->matbeg, sizeof (int) * (A->matcolsize + nrows));
	//rval = ILLutil_reallocrus_count ((void **) &(A->matbeg),
	//                                 A->matcolsize + nrows, sizeof (int));
	//ILL_CLEANUP_IF (rval);

	A->matind = EGrealloc (A->matind, sizeof (int) * (A->matsize + nrows));
	//rval = ILLutil_reallocrus_count ((void **) &(A->matind),
	//                                 A->matsize + nrows, sizeof (int));
	//ILL_CLEANUP_IF (rval);
	EGLPNUM_TYPENAME_EGlpNumReallocArray (&(A->matval), A->matsize + nrows);

	for (i = 0; i < nrows; i++)
	{
		A->matind[A->matsize + i] = -1;
	}

	aindex = A->matsize - A->matfree;

	for (i = 0; i < nrows; i++)
	{
		lp->rowmap[i] = ncols;
		EGLPNUM_TYPENAME_EGlpNumZero (lp->obj[ncols]);
		A->matcnt[ncols] = 1;
		A->matbeg[ncols] = aindex;
		A->matind[aindex] = i;
		switch (sense[i])
		{
		case 'E':									/* Arificial */
			EGLPNUM_TYPENAME_EGlpNumZero (lp->lower[ncols]);
			EGLPNUM_TYPENAME_EGlpNumZero (lp->upper[ncols]);
			EGLPNUM_TYPENAME_EGlpNumOne (A->matval[aindex]);
			break;
		case 'G':									/* Surplus   */
			EGLPNUM_TYPENAME_EGlpNumZero (lp->lower[ncols]);
			EGLPNUM_TYPENAME_EGlpNumCopy (lp->upper[ncols], EGLPNUM_TYPENAME_ILL_MAXDOUBLE);
			EGLPNUM_TYPENAME_EGlpNumOne (A->matval[aindex]);
			EGLPNUM_TYPENAME_EGlpNumSign (A->matval[aindex]);
			break;
		case 'L':									/* Slack     */
			EGLPNUM_TYPENAME_EGlpNumZero (lp->lower[ncols]);
			EGLPNUM_TYPENAME_EGlpNumCopy (lp->upper[ncols], EGLPNUM_TYPENAME_ILL_MAXDOUBLE);
			EGLPNUM_TYPENAME_EGlpNumOne (A->matval[aindex]);
			break;
		case 'R':									/* Range     */
			EGLPNUM_TYPENAME_EGlpNumZero (lp->lower[ncols]);
			EGLPNUM_TYPENAME_EGlpNumCopy (lp->upper[ncols], lp->rangeval[i]);
			EGLPNUM_TYPENAME_EGlpNumOne (A->matval[aindex]);
			EGLPNUM_TYPENAME_EGlpNumSign (A->matval[aindex]);
			break;
		default:
			QSlog("unknown sense %c in EGLPNUM_TYPENAME_ILLlp_add_logicals", sense[i]);
			rval = 1;
			goto CLEANUP;
		}
		ncols++;
		nzcount++;
		aindex++;
	}

	lp->ncols = ncols;
	lp->nzcount = nzcount;
	A->matcols = ncols;

	lp->colsize += nrows;
	A->matsize += nrows;
	A->matcolsize += nrows;

	if (lp->rA)
	{
		EGLPNUM_TYPENAME_ILLlp_rows_clear (lp->rA);
	}
	else
	{
		ILL_SAFE_MALLOC (lp->rA, 1, EGLPNUM_TYPENAME_ILLlp_rows);
	}

	rval = EGLPNUM_TYPENAME_ILLlp_rows_init (lp->rA, lp, 1);
	ILL_CLEANUP_IF (rval);

CLEANUP:

	ILL_RETURN (rval, "EGLPNUM_TYPENAME_ILLlp_add_logicals");
}

int EGLPNUM_TYPENAME_ILLlp_scale (
	EGLPNUM_TYPENAME_ILLlpdata * lp)
{
	int rval = 0;
	int i, j, k, col, row, nstruct, start, stop;
	EGLPNUM_TYPENAME_ILLmatrix *A;
	EGLPNUM_TYPE rho;
	EGLPNUM_TYPE *gama = 0;

	EGLPNUM_TYPENAME_EGlpNumInitVar (rho);

	/* Columns - divide by largest absolute value */

	if (!lp)
	{
		ILL_ERROR (rval, "EGLPNUM_TYPENAME_ILLlp_scale called with a NULL pointer");
	}

	if (lp->nrows == 0 || lp->ncols == 0)
		goto CLEANUP;

	A = &lp->A;
	nstruct = lp->nstruct;

	for (j = 0; j < nstruct; j++)
	{
		col = lp->structmap[j];
		EGLPNUM_TYPENAME_EGlpNumZero (rho);

		start = A->matbeg[col];
		stop = start + A->matcnt[col];

		for (k = start; k < stop; k++)
		{
			EGLPNUM_TYPENAME_EGlpNumSetToMaxAbs (rho, A->matval[k]);
		}

		if (EGLPNUM_TYPENAME_EGlpNumIsGreatZero (rho))
		{
			for (k = start; k < stop; k++)
			{
				EGLPNUM_TYPENAME_EGlpNumDivTo (A->matval[k], rho);
			}
			EGLPNUM_TYPENAME_EGlpNumDivTo (lp->obj[col], rho);
			if (EGLPNUM_TYPENAME_EGlpNumIsNeqq (lp->lower[col], EGLPNUM_TYPENAME_ILL_MINDOUBLE))
				EGLPNUM_TYPENAME_EGlpNumMultTo (lp->lower[col], rho);
			if (EGLPNUM_TYPENAME_EGlpNumIsNeqq (lp->upper[col], EGLPNUM_TYPENAME_ILL_MAXDOUBLE))
				EGLPNUM_TYPENAME_EGlpNumMultTo (lp->upper[col], rho);
		}
	}

	gama = EGLPNUM_TYPENAME_EGlpNumAllocArray (lp->nrows);
	for (i = 0; i < lp->nrows; i++)
	{
		EGLPNUM_TYPENAME_EGlpNumZero (gama[i]);
	}

	for (j = 0; j < nstruct; j++)
	{
		col = lp->structmap[j];
		start = A->matbeg[col];
		stop = start + A->matcnt[col];

		for (k = start; k < stop; k++)
		{
			row = A->matind[k];
			EGLPNUM_TYPENAME_EGlpNumSetToMaxAbs (gama[row], A->matval[k]);
		}
	}

	for (j = 0; j < nstruct; j++)
	{
		col = lp->structmap[j];
		start = A->matbeg[col];
		stop = start + A->matcnt[col];

		for (k = start; k < stop; k++)
		{
			row = A->matind[k];
			if (EGLPNUM_TYPENAME_EGlpNumIsGreatZero (gama[row]))
			{
				EGLPNUM_TYPENAME_EGlpNumDivTo (A->matval[k], gama[row]);
			}
		}
	}

	for (i = 0; i < lp->nrows; i++)
	{
		if (EGLPNUM_TYPENAME_EGlpNumIsGreatZero ( gama[i]))
		{
			EGLPNUM_TYPENAME_EGlpNumDivTo (lp->rhs[i], gama[i]);
			col = lp->rowmap[i];
			if (EGLPNUM_TYPENAME_EGlpNumIsNeqq (lp->upper[col], EGLPNUM_TYPENAME_ILL_MAXDOUBLE))
			{
				EGLPNUM_TYPENAME_EGlpNumDivTo (lp->upper[col], gama[i]);	/* Ranged row */
			}
		}
	}

	if (lp->rA)
	{															/* Need to clear the row version of data */
		EGLPNUM_TYPENAME_ILLlp_rows_clear (lp->rA);
		ILL_IFFREE (lp->rA, EGLPNUM_TYPENAME_ILLlp_rows);
	}


CLEANUP:

	EGLPNUM_TYPENAME_EGlpNumClearVar (rho);
	EGLPNUM_TYPENAME_EGlpNumFreeArray (gama);
	ILL_RETURN (rval, "EGLPNUM_TYPENAME_ILLlp_scale");
}

int EGLPNUM_TYPENAME_ILLlp_presolve (
	EGLPNUM_TYPENAME_ILLlpdata * lp,
	int pre_types)
{
	int rval = 0;
	int status = ILL_LP_STATUS_OK;
	EGLPNUM_TYPENAME_ILLlp_predata *pre = 0;
	EGLPNUM_TYPENAME_ILLlp_sinfo *info = 0;

	if (!lp)
	{
		QSlog("EGLPNUM_TYPENAME_ILLlp_presolve called with a NULL pointer");
		rval = 1;
		goto CLEANUP;
	}


/*
    ILLlpdata_writelp (lp, 0);
*/

	ILL_SAFE_MALLOC (pre, 1, EGLPNUM_TYPENAME_ILLlp_predata);
	EGLPNUM_TYPENAME_ILLlp_predata_init (pre);

	ILL_SAFE_MALLOC (info, 1, EGLPNUM_TYPENAME_ILLlp_sinfo);
	EGLPNUM_TYPENAME_ILLlp_sinfo_init (info);

	rval = simple_presolve (lp, pre, info, pre_types, &status);
	ILL_CLEANUP_IF (rval);
	if (status != ILL_LP_STATUS_OK)
	{
		QSlog("simple_presolve returned with bad status");
		rval = 1;
		goto CLEANUP;
	}

/*
    rval = EGLPNUM_TYPENAME_ILLlp_sinfo_print (info);
    ILL_CLEANUP_IF (rval);
*/

CLEANUP:

	if (rval)
	{
		if (pre)
		{
			EGLPNUM_TYPENAME_ILLlp_predata_free (pre);
			ILL_IFFREE (pre, EGLPNUM_TYPENAME_ILLlp_predata);
		}

		if (info)
		{
			EGLPNUM_TYPENAME_ILLlp_sinfo_free (info);
			ILL_IFFREE (info, EGLPNUM_TYPENAME_ILLlp_sinfo);
		}
	}
	else
	{
		lp->presolve = pre;
		lp->sinfo = info;
	}

	ILL_RETURN (rval, "EGLPNUM_TYPENAME_ILLlp_presolve");
}


#if 0
int ILLlp_presolve_addrow (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int cnt,
	int *ind,
	double *val,
	double rhs)
{
	int rval = 0;
	EGLPNUM_TYPENAME_ILLlpdata *qslp;
	EGLPNUM_TYPENAME_ILLlp_sinfo *S;
	EGLPNUM_TYPENAME_ILLmatrix *A;

	/* This will need to evolve into a function that handles the task */
	/* of working through the presolve data to determine the new LP   */
	/* created when a row is added to the original LP.                */

	/* The copies of the obj and bound used in the simplex code are   */
	/* also updated in this function.                                 */

	if (!lp)
	{
		QSlog("ILLlp_presolve_addrow is called without an LP");
		rval = 1;
		goto CLEANUP;
	}

	if (lp->presolve != 0)
	{
		QSlog("Not yet set up to handle addrows after presolve");
		rval = 1;
		goto CLEANUP;
	}

	qslp = lp->O;
	S = qslp->sinfo;
	A = S->A;


	rval = ILLlib_matrix_addrow (A, cnt, ind, val, rhs);
	ILL_CLEANUP_IF (rval);


CLEANUP:

	ILL_RETURN (rval, "ILLlp_presolve_addrow");
}
#endif


static int simple_presolve (
	EGLPNUM_TYPENAME_ILLlpdata * lp,
	EGLPNUM_TYPENAME_ILLlp_predata * pre,
	EGLPNUM_TYPENAME_ILLlp_sinfo * info,
	int pre_types,
	int *status)
{
	int rval = 0;
	int i, hit, newhit;
	graph G;

	if (status)
		*status = ILL_LP_STATUS_OK;
	init_graph (&G);

	if (!lp)
	{
		QSlog("simple_presolve called with a NULL pointer");
		rval = 1;
		goto CLEANUP;
	}

	QSlog("Initial Rows = %d, Cols = %d, Nzcount = %d",
							lp->nrows, lp->ncols, lp->nzcount);

	rval = build_graph (lp, &G);
	ILL_CLEANUP_IF (rval);
	if (debug)
		dump_graph (&G);

	if (pre_types & EGLPNUM_TYPENAME_ILL_PRE_FIXED)
	{
		rval = fixed_variables (&G, pre);
		ILL_CLEANUP_IF (rval);
	}

	do
	{
		hit = 0;
		if (pre_types & EGLPNUM_TYPENAME_ILL_PRE_SINGLE_ROW)
		{
			rval = singleton_rows (&G, pre, &newhit);
			ILL_CLEANUP_IF (rval);
			hit += newhit;
		}

		if (pre_types & EGLPNUM_TYPENAME_ILL_PRE_FORCING)
		{
			rval = forcing_constraints (&G, pre, &newhit);
			ILL_CLEANUP_IF (rval);
			hit += newhit;
		}

		if (pre_types & EGLPNUM_TYPENAME_ILL_PRE_SINGLE_COL)
		{
			rval = singleton_columns (&G, pre, &newhit);
			ILL_CLEANUP_IF (rval);
			hit += newhit;
		}

		if (pre_types & EGLPNUM_TYPENAME_ILL_PRE_DUPLICATE_ROW)
		{
			rval = duplicate_rows (&G, &newhit);
			ILL_CLEANUP_IF (rval);
			hit += newhit;
		}

		if (pre_types & EGLPNUM_TYPENAME_ILL_PRE_DUPLICATE_COL)
		{
			rval = duplicate_cols (&G, &newhit);
			ILL_CLEANUP_IF (rval);
			hit += newhit;
		}


/*
        {
            int k, cnt = 0;
            for (i = 0; i < G.ncols; i++) {
                if (G.cols[i].del == 0) {
                    for (k = 0; k < G.cols[i].deg; k++)  {
                        if (G.cols[i].adj[k]->del == 0) {
                            cnt++;
                        }
                    }
                }
            }
            QSlog("Current NZCOUNT = %d", cnt);
        }
*/
	} while (hit);

	if (EGLPNUM_TYPENAME_ILL_PRE_EMPTY_COL)
	{
		rval = empty_columns (&G, pre);
		ILL_CLEANUP_IF (rval);
	}

	if (debug)
	{
		QSlog("Operations");
		for (i = 0; i < pre->opcount; i++)
		{
			switch (pre->oplist[i].ptype)
			{
			case ILL_PRE_DELETE_EMPTY_ROW:
				QSlog("Delete Empty Row: %d", pre->oplist[i].rowindex);
				break;
			case ILL_PRE_DELETE_SINGLETON_ROW:
				QSlog("Delete Singleton Row: %d (col %d)",
										pre->oplist[i].rowindex, pre->oplist[i].colindex);
				dump_line (&pre->oplist[i].line);
				break;
			case ILL_PRE_DELETE_FIXED_VARIABLE:
				QSlog("Delete Fixed Variable: %d", pre->oplist[i].colindex);
				dump_line (&pre->oplist[i].line);
				break;
			case ILL_PRE_DELETE_FORCED_VARIABLE:
				QSlog("Delete Forced Variable: %d", pre->oplist[i].colindex);
				dump_line (&pre->oplist[i].line);
				break;
			case ILL_PRE_DELETE_SINGLETON_VARIABLE:
				QSlog("Delete Singleton Variable: %d", pre->oplist[i].colindex);
				dump_line (&pre->oplist[i].line);
				break;
			case ILL_PRE_DELETE_FREE_SINGLETON_VARIABLE:
				QSlog("Delete Free Singleton Variable: %d",
										pre->oplist[i].colindex);
				dump_line (&pre->oplist[i].line);
				break;
			case ILL_PRE_DELETE_EMPTY_COLUMN:
				QSlog("Delete Empty Column: %d", pre->oplist[i].colindex);
				dump_line (&pre->oplist[i].line);
				break;
			default:
				QSlog("unknon presolve operation");
				rval = 1;
				goto CLEANUP;
			}
		}
	}

	rval = grab_lp_info (&G, lp->colnames, info);
	ILL_CLEANUP_IF (rval);

/*
    QSlog("Final Rows = %d, Cols = %d, Nzcount = %d",
                info->nrows, info->ncols, info->nzcount);
*/


CLEANUP:

	free_graph (&G);
	ILL_RETURN (rval, "simple_presolve");
}

static int grab_lp_line (
	graph * G,
	int indx,
	EGLPNUM_TYPENAME_ILLlp_preline * line,
	int row_or_col)
{
	int rval = 0;
	int k, cnt;
	node *n;

	if (row_or_col == 0)
		n = &G->rows[indx];
	else
		n = &G->cols[indx];

	line->count = 0;

	for (k = 0; k < n->deg; k++)
	{
		if (n->adj[k]->del == 0)
		{
			line->count++;
		}
	}

	if (line->count)
	{
		ILL_SAFE_MALLOC (line->ind, line->count, int);

		line->val = EGLPNUM_TYPENAME_EGlpNumAllocArray (line->count);
		if (!line->ind || !line->val)
		{
			QSlog("out of memory in grab_lp_line");
			rval = 1;
			goto CLEANUP;
		}
		for (k = 0, cnt = 0; k < n->deg; k++)
		{
			if (n->adj[k]->del == 0)
			{
				line->ind[cnt] = n->adj[k]->row;
				EGLPNUM_TYPENAME_EGlpNumCopy (line->val[cnt], n->adj[k]->coef);
				cnt++;
			}
		}
	}

	if (row_or_col == 0)
	{
		EGLPNUM_TYPENAME_EGlpNumCopy (line->rhs, n->rhs);
	}
	else
	{
		EGLPNUM_TYPENAME_EGlpNumCopy (line->obj, n->obj);
		EGLPNUM_TYPENAME_EGlpNumCopy (line->lower, n->lower);
		EGLPNUM_TYPENAME_EGlpNumCopy (line->upper, n->upper);
	}

	line->row_or_col = row_or_col;

CLEANUP:

	ILL_RETURN (rval, "grab_lp_line");
}

static void dump_line (
	EGLPNUM_TYPENAME_ILLlp_preline * line)
{
	int k;

	if (line->row_or_col == 0)
	{
		for (k = 0; k < line->count; k++)
		{
			QSlog(" C%d->%g", line->ind[k], EGLPNUM_TYPENAME_EGlpNumToLf (line->val[k]));
		}
		QSlog(" RHS->%g", EGLPNUM_TYPENAME_EGlpNumToLf (line->rhs));
	}
	else
	{
		for (k = 0; k < line->count; k++)
		{
			QSlog(" R%d->%g", line->ind[k], EGLPNUM_TYPENAME_EGlpNumToLf (line->val[k]));
		}
		QSlog(" Obj->%g  LB->%g  UB->%g", EGLPNUM_TYPENAME_EGlpNumToLf (line->obj),
						EGLPNUM_TYPENAME_EGlpNumToLf (line->lower), EGLPNUM_TYPENAME_EGlpNumToLf (line->upper));
	}
}

static int grab_lp_info (
	graph * G,
	char **colnames,
	EGLPNUM_TYPENAME_ILLlp_sinfo * info)
{
	int rval = 0;
	int ncols = 0, nrows = 0, nzcount = 0;
	int i, j, k, cnt, len;
	node *grows = G->rows;
	node *gcols = G->cols;
	int *tdeg = 0;
	int *map = 0;
	char *buf = 0;
	EGLPNUM_TYPENAME_ILLmatrix *A = &info->A;

	ILL_SAFE_MALLOC (tdeg, G->ncols, int);
	ILL_SAFE_MALLOC (map, G->nrows, int);

	if (!tdeg || !map)
	{
		QSlog("out of memory in grab_lp_info");
		rval = 1;
		goto CLEANUP;
	}

	for (i = 0; i < G->nrows; i++)
	{
		if (grows[i].del == 0)
		{
			map[i] = nrows;
			nrows++;
		}
	}

	for (j = 0; j < G->ncols; j++)
	{
		if (gcols[j].del == 0)
		{
			tdeg[ncols] = 0;
			for (k = 0; k < gcols[j].deg; k++)
			{
				if (gcols[j].adj[k]->del == 0)
				{
					tdeg[ncols]++;
					nzcount++;
				}
			}
			ncols++;
		}
	}

	info->ncols = ncols;
	info->nrows = nrows;
	info->nzcount = nzcount;

	info->rowsize = nrows;
	info->colsize = ncols;

	info->rhs = EGLPNUM_TYPENAME_EGlpNumAllocArray (nrows);
	info->obj = EGLPNUM_TYPENAME_EGlpNumAllocArray (ncols);
	info->upper = EGLPNUM_TYPENAME_EGlpNumAllocArray (ncols);
	info->lower = EGLPNUM_TYPENAME_EGlpNumAllocArray (ncols);
	A->matval = EGLPNUM_TYPENAME_EGlpNumAllocArray (info->nzcount + 1);
	ILL_SAFE_MALLOC (A->matind, info->nzcount + 1, int);
	ILL_SAFE_MALLOC (A->matcnt, info->colsize, int);
	ILL_SAFE_MALLOC (A->matbeg, info->colsize, int);

	if (!info->rhs || !info->obj || !info->lower || !info->upper ||
			!A->matval || !A->matind || !A->matcnt || !A->matbeg)
	{
		QSlog("out of memory in grab_lp");
		rval = 1;
		goto CLEANUP;
	}

	A->matind[info->nzcount] = -1;
	A->matsize = info->nzcount + 1;
	A->matcolsize = info->colsize;
	A->matfree = 1;
	A->matcols = ncols;
	A->matrows = nrows;


	nrows = 0;
	for (i = 0; i < G->nrows; i++)
	{
		if (grows[i].del == 0)
		{
			EGLPNUM_TYPENAME_EGlpNumCopy (info->rhs[nrows], grows[i].rhs);
			nrows++;
		}
	}

	ncols = 0;
	cnt = 0;
	for (j = 0; j < G->ncols; j++)
	{
		if (gcols[j].del == 0)
		{
			EGLPNUM_TYPENAME_EGlpNumCopy (info->obj[ncols], gcols[j].obj);
			EGLPNUM_TYPENAME_EGlpNumCopy (info->lower[ncols], gcols[j].lower);
			EGLPNUM_TYPENAME_EGlpNumCopy (info->upper[ncols], gcols[j].upper);
			A->matcnt[ncols] = tdeg[ncols];
			A->matbeg[ncols] = cnt;
			for (k = 0; k < gcols[j].deg; k++)
			{
				if (gcols[j].adj[k]->del == 0)
				{
					EGLPNUM_TYPENAME_EGlpNumCopy (A->matval[cnt], gcols[j].adj[k]->coef);
					A->matind[cnt] = map[gcols[j].adj[k]->row];
					cnt++;
				}
			}
			ncols++;
		}
	}

	if (colnames)
	{
		ILL_SAFE_MALLOC (info->colnames, info->colsize, char *);

		if (!info->colnames)
		{
			QSlog("out of memory in grab_lp");
			rval = 1;
			goto CLEANUP;
		}
		for (j = 0; j < info->colsize; j++)
		{
			info->colnames[j] = 0;
		}

		ILL_SAFE_MALLOC (buf, ILL_namebufsize, char);

		if (!buf)
		{
			QSlog("out of memory in grab_lp");
			rval = 1;
			goto CLEANUP;
		}
		ncols = 0;
		for (j = 0; j < G->ncols; j++)
		{
			if (gcols[j].del == 0)
			{
				if (gcols[j].coltype == ILL_PRE_COL_STRUC)
				{
					len = strlen (colnames[j]) + 1;
					ILL_SAFE_MALLOC (info->colnames[ncols], len, char);

					if (!info->colnames[ncols])
					{
						QSlog("out of memory in grab_lp");
						rval = 1;
						goto CLEANUP;
					}
					strcpy (info->colnames[ncols], colnames[j]);
				}
				else
				{
					for (k = 0; k < gcols[j].deg; k++)
					{
						if (gcols[j].adj[k]->del == 0)
						{
							i = gcols[j].adj[k]->row;
							break;
						}
					}
					if (k == gcols[j].deg)
					{
						QSlog("problem with graph in grab_lp");
						rval = 1;
						goto CLEANUP;
					}
					sprintf (buf, "s%d", i);
					len = strlen (buf) + 1;
					ILL_SAFE_MALLOC (info->colnames[ncols], len, char);

					if (!info->colnames[ncols])
					{
						QSlog("out of memory in grab_lp");
						rval = 1;
						goto CLEANUP;
					}
					strcpy (info->colnames[ncols], buf);
				}
				ncols++;
			}
		}
	}

/* ADD STRUCT VARIABLE STUFF */


CLEANUP:

	if (rval)
	{
		EGLPNUM_TYPENAME_ILLlp_sinfo_free (info);
	}
	ILL_IFFREE (tdeg, int);
	ILL_IFFREE (map, int);
	ILL_IFFREE (buf, char);

	ILL_RETURN (rval, "grab_lp_info");
}

static int fixed_variables (
	graph * G,
	EGLPNUM_TYPENAME_ILLlp_predata * pre)
{
	int rval = 0;
	int j;
	int ncols = G->ncols;
	node *cols = G->cols;
	EGLPNUM_TYPENAME_ILLlp_preop *op = 0;

	for (j = 0; j < ncols; j++)
	{
		if (cols[j].del == 0)
		{
			if (EGLPNUM_TYPENAME_EGlpNumIsEqqual (cols[j].lower, cols[j].upper))
			{
				rval = get_next_preop (pre, &op);
				ILL_CLEANUP_IF (rval);

				op->colindex = j;
				op->rowindex = -1;
				op->ptype = ILL_PRE_DELETE_FIXED_VARIABLE;

				rval = grab_lp_line (G, op->colindex, &op->line, 1);
				ILL_CLEANUP_IF (rval);
				pre->opcount++;

				set_fixed_variable (G, j, cols[j].lower);
			}
		}
	}

CLEANUP:

	ILL_RETURN (rval, "fixed_variables");
}

static int empty_columns (
	graph * G,
	EGLPNUM_TYPENAME_ILLlp_predata * pre)
{
	int rval = 0;
	int j, k;
	int ncols = G->ncols;
	node *cols = G->cols;
	EGLPNUM_TYPENAME_ILLlp_preop *op = 0;
	EGLPNUM_TYPE objtmp;

	EGLPNUM_TYPENAME_EGlpNumInitVar (objtmp);

	for (j = 0; j < ncols; j++)
	{
		if (cols[j].del == 0)
		{
			for (k = 0; k < cols[j].deg; k++)
			{
				if (cols[j].adj[k]->del == 0)
					break;
			}
			if (k == cols[j].deg)
			{
				rval = get_next_preop (pre, &op);
				ILL_CLEANUP_IF (rval);

				op->colindex = j;
				op->rowindex = -1;
				op->ptype = ILL_PRE_DELETE_EMPTY_COLUMN;

				rval = grab_lp_line (G, op->colindex, &op->line, 1);
				ILL_CLEANUP_IF (rval);
				pre->opcount++;
				EGLPNUM_TYPENAME_EGlpNumCopy (objtmp, cols[j].obj);
				if (G->objsense < 0)
					EGLPNUM_TYPENAME_EGlpNumSign (objtmp);
				if (!EGLPNUM_TYPENAME_EGlpNumIsNeqZero (objtmp, ILL_PRE_FEAS_TOL))
				{
					set_fixed_variable (G, j, cols[j].lower);
				}
				else if (EGLPNUM_TYPENAME_EGlpNumIsGreatZero (objtmp))
				{
					if (EGLPNUM_TYPENAME_EGlpNumIsEqqual (cols[j].lower, EGLPNUM_TYPENAME_ILL_MINDOUBLE))
					{
						QSlog("unbounded prob detected in empty_columns");
						QSlog("col %d, obj %g", j, EGLPNUM_TYPENAME_EGlpNumToLf (cols[j].obj));
						rval = 1;
						goto CLEANUP;
					}
					else
					{
						set_fixed_variable (G, j, cols[j].lower);
					}
				}
				else if (EGLPNUM_TYPENAME_EGlpNumIsLessZero (objtmp))
				{
					if (EGLPNUM_TYPENAME_EGlpNumIsEqqual (cols[j].upper, EGLPNUM_TYPENAME_ILL_MAXDOUBLE))
					{
						QSlog("unbounded prob detected in empty_columns");
						QSlog("col %d, obj %g", j, EGLPNUM_TYPENAME_EGlpNumToLf (cols[j].obj));
						rval = 1;
						goto CLEANUP;
					}
					else
					{
						set_fixed_variable (G, j, cols[j].upper);
					}
				}
				else
				{
					set_fixed_variable (G, j, cols[j].lower);
				}
			}
		}
	}

CLEANUP:

	EGLPNUM_TYPENAME_EGlpNumClearVar (objtmp);
	ILL_RETURN (rval, "empty_columns");
}

static int singleton_rows (
	graph * G,
	EGLPNUM_TYPENAME_ILLlp_predata * pre,
	int *hit)
{
	int rval = 0;
	int rowindex, i, k, h;
	int nrows = G->nrows;
	node *rows = G->rows;
	node *cols = G->cols;
	node *r, *c;
	edge *pivot, *f;
	intptr *next, *list = 0;
	int *tdeg = 0;
	EGLPNUM_TYPE val;
	EGLPNUM_TYPENAME_ILLlp_preop *op = 0;

	EGLPNUM_TYPENAME_EGlpNumInitVar (val);

	*hit = 0;
	if (G->nrows == 0)
		goto CLEANUP;

	ILL_SAFE_MALLOC (tdeg, G->nrows, int);

	if (!tdeg)
	{
		QSlog("out of memory in singleton_rows");
		rval = 1;
		goto CLEANUP;
	}

	for (i = 0; i < nrows; i++)
	{
		if (rows[i].del == 0)
		{
			tdeg[i] = 0;
			for (k = 0; k < rows[i].deg; k++)
			{
				if (rows[i].adj[k]->del == 0)
				{
					tdeg[i]++;
				}
			}
			if (tdeg[i] <= 1)
			{
				rval = add_to_list (&G->intptrworld, &list, i);
				ILL_CLEANUP_IF (rval);
			}
		}
	}

	while (list)
	{
		(*hit)++;
		rowindex = list->this_val;
		next = list->next;
		intptrfree (&G->intptrworld, list);
		list = next;

		rval = get_next_preop (pre, &op);
		ILL_CLEANUP_IF (rval);

		r = &rows[rowindex];

		if (tdeg[rowindex] == 0)
		{
			if (EGLPNUM_TYPENAME_EGlpNumIsNeqZero (r->rhs, ILL_PRE_FEAS_TOL))
			{
				QSlog("infeasible row detected in singleton_row");
				QSlog("empty row with rhs = %g", EGLPNUM_TYPENAME_EGlpNumToLf (r->rhs));
				rval = 1;
				goto CLEANUP;
			}
			op->ptype = ILL_PRE_DELETE_EMPTY_ROW;
			op->rowindex = rowindex;
		}
		else
		{
			/*  Find the "pivot" entry and colum */

			for (k = 0; k < r->deg; k++)
			{
				if (r->adj[k]->del == 0)
					break;
			}
			if (k == r->deg)
			{
				QSlog("lost an edge in singleton_rows");
				rval = 1;
				goto CLEANUP;
			}

			pivot = r->adj[k];
			c = &cols[pivot->col];

			/*  Store data from operation (incluing the col coefs) */

			op->ptype = ILL_PRE_DELETE_SINGLETON_ROW;
			op->rowindex = rowindex;
			op->colindex = c - cols;
			EGLPNUM_TYPENAME_EGlpNumCopy (op->line.rhs, r->rhs);
			rval = grab_lp_line (G, op->colindex, &op->line, 1);
			ILL_CLEANUP_IF (rval);

			/*  Fix the x[c] to its rhs value */
			/*val = r->rhs / pivot->coef; */
			EGLPNUM_TYPENAME_EGlpNumCopyFrac (val, r->rhs, pivot->coef);
			/* if (val < c->lower - ILL_PRE_FEAS_TOL ||
			 * val > c->upper + ILL_PRE_FEAS_TOL) */
			if (EGLPNUM_TYPENAME_EGlpNumIsSumLess (val, ILL_PRE_FEAS_TOL, c->lower) ||
					EGLPNUM_TYPENAME_EGlpNumIsSumLess (c->upper, ILL_PRE_FEAS_TOL, val))
			{
				QSlog("infeasible bounds detected in singleton_row %d", rowindex);
				QSlog("lower->%g  upper->%g  val = %g",
										EGLPNUM_TYPENAME_EGlpNumToLf (c->lower), EGLPNUM_TYPENAME_EGlpNumToLf (c->upper),
										EGLPNUM_TYPENAME_EGlpNumToLf (val));
				rval = 1;
				goto CLEANUP;
			}
			else
			{
				EGLPNUM_TYPENAME_EGlpNumCopy (c->lower, val);
				EGLPNUM_TYPENAME_EGlpNumCopy (c->upper, val);
			}

			/*  Delete x[c] from other rows (and adjust their rhs) */

			c->del = 1;

			for (h = 0; h < c->deg; h++)
			{
				f = c->adj[h];
				if (f->del == 0)
				{
					/*rows[f->row].rhs -= (f->coef * c->lower); */
					EGLPNUM_TYPENAME_EGlpNumSubInnProdTo (rows[f->row].rhs, f->coef, c->lower);
					tdeg[f->row]--;
					if (tdeg[f->row] == 1)
					{
						if (f == pivot)
						{
							QSlog("bad pivot element");
							rval = 1;
							goto CLEANUP;
						}
						rval = add_to_list (&G->intptrworld, &list, f->row);
						ILL_CLEANUP_IF (rval);
					}
					f->del = 1;
				}
			}
		}

		r->del = 1;
		pre->opcount++;
	}

CLEANUP:

	ILL_IFFREE (tdeg, int);

	intptr_listfree (&G->intptrworld, list);
	EGLPNUM_TYPENAME_EGlpNumClearVar (val);
	ILL_RETURN (rval, "singleton_rows");
}

static int forcing_constraints (
	graph * G,
	EGLPNUM_TYPENAME_ILLlp_predata * pre,
	int *hit)
{
	int rval = 0;
	int i, j, k, ts;
	node *rows = G->rows;
	node *cols = G->cols;
	edge *e;
	int nrows = G->nrows;
	EGLPNUM_TYPE ub, lb;
	EGLPNUM_TYPENAME_ILLlp_preop *op = 0;

	EGLPNUM_TYPENAME_EGlpNumInitVar (ub);
	EGLPNUM_TYPENAME_EGlpNumInitVar (lb);

	*hit = 0;

	for (i = 0; i < nrows; i++)
	{
		if (rows[i].del == 0)
		{
			get_implied_rhs_bounds (G, i, &lb, &ub);
			if (EGLPNUM_TYPENAME_EGlpNumIsSumLess (rows[i].rhs, ILL_PRE_FEAS_TOL, lb) ||
					EGLPNUM_TYPENAME_EGlpNumIsSumLess (ub, ILL_PRE_FEAS_TOL, rows[i].rhs))
			{
				QSlog("infeasible row detected in forcing_constraints");
				QSlog("Row %d:  RHS->%g  LBnd->%g  UBnd->%g",
										i, EGLPNUM_TYPENAME_EGlpNumToLf (rows[i].rhs),
										EGLPNUM_TYPENAME_EGlpNumToLf (lb), EGLPNUM_TYPENAME_EGlpNumToLf (ub));
				rval = 1;
				goto CLEANUP;
			}
			else if (EGLPNUM_TYPENAME_EGlpNumIsDiffLess (rows[i].rhs, ILL_PRE_FEAS_TOL, lb) ||
							 EGLPNUM_TYPENAME_EGlpNumIsDiffLess (ub, ILL_PRE_FEAS_TOL, rows[i].rhs))
			{
				(*hit)++;
				ts = (EGLPNUM_TYPENAME_EGlpNumIsDiffLess (rows[i].rhs, ILL_PRE_FEAS_TOL, lb) ? 0 : 1);
				for (k = 0; k < rows[i].deg; k++)
				{
					e = rows[i].adj[k];
					if (e->del == 0)
					{
						j = e->col;

						rval = get_next_preop (pre, &op);
						ILL_CLEANUP_IF (rval);

						op->colindex = j;
						op->rowindex = i;
						op->ptype = ILL_PRE_DELETE_FORCED_VARIABLE;

						rval = grab_lp_line (G, j, &op->line, 1);
						ILL_CLEANUP_IF (rval);
						pre->opcount++;

						if ((ts == 0 && EGLPNUM_TYPENAME_EGlpNumIsLessZero (e->coef)) ||
								(ts == 1 && EGLPNUM_TYPENAME_EGlpNumIsGreatZero (e->coef)))
						{
							set_fixed_variable (G, j, cols[j].upper);
						}
						else
						{
							set_fixed_variable (G, j, cols[j].lower);
						}
					}
				}
			}
		}
	}

CLEANUP:

	EGLPNUM_TYPENAME_EGlpNumClearVar (ub);
	EGLPNUM_TYPENAME_EGlpNumClearVar (lb);
	ILL_RETURN (rval, "forcing_constraints");
}

static int singleton_columns (
	graph * G,
	EGLPNUM_TYPENAME_ILLlp_predata * pre,
	int *hit)
{
	int rval = 0;
	int ncols = G->ncols;
	int j, k, deg, rdeg, single = 0, irow;
	EGLPNUM_TYPE lb, ub, b, eb;
	node *cols = G->cols;
	node *rows = G->rows;
	edge *b_edge;
	EGLPNUM_TYPENAME_ILLlp_preop *op = 0;
	EGLPNUM_TYPE newub, newlb;
	EGLPNUM_TYPE a, c, l, u;

	EGLPNUM_TYPENAME_EGlpNumInitVar (lb);
	EGLPNUM_TYPENAME_EGlpNumInitVar (ub);
	EGLPNUM_TYPENAME_EGlpNumInitVar (eb);
	EGLPNUM_TYPENAME_EGlpNumInitVar (b);
	EGLPNUM_TYPENAME_EGlpNumInitVar (newlb);
	EGLPNUM_TYPENAME_EGlpNumInitVar (newub);
	EGLPNUM_TYPENAME_EGlpNumInitVar (a);
	EGLPNUM_TYPENAME_EGlpNumInitVar (c);
	EGLPNUM_TYPENAME_EGlpNumInitVar (l);
	EGLPNUM_TYPENAME_EGlpNumInitVar (u);

	*hit = 0;
	if (G->ncols == 0)
		goto CLEANUP;

	for (j = 0; j < ncols; j++)
	{
		if (cols[j].del == 0)
		{
			deg = 0;
			for (k = 0; k < cols[j].deg && deg <= 1; k++)
			{
				if (cols[j].adj[k]->del == 0)
				{
					single = k;
					deg++;
				}
			}
			if (deg == 1)
			{
				irow = cols[j].adj[single]->row;
				EGLPNUM_TYPENAME_EGlpNumCopy (b, cols[j].adj[single]->coef);
				b_edge = cols[j].adj[single];

				get_implied_variable_bounds (G, j, b_edge, &lb, &ub);

				/*if (lb >= cols[j].lower && ub <= cols[j].upper) */
				if (EGLPNUM_TYPENAME_EGlpNumIsLeq (cols[j].lower, lb) &&
						EGLPNUM_TYPENAME_EGlpNumIsLeq (ub, cols[j].upper))
				{
					edge *a_edge;

					/*  The jth variable can be substituted out of problem */
					/*        x = (c/b) - (a/b)y                           */


					rval = get_next_preop (pre, &op);
					ILL_CLEANUP_IF (rval);

					op->colindex = j;
					op->rowindex = irow;
					op->ptype = ILL_PRE_DELETE_FREE_SINGLETON_VARIABLE;

					rval = grab_lp_line (G, irow, &op->line, 0);
					ILL_CLEANUP_IF (rval);
					pre->opcount++;

					/*  Adjust the objective function                      */
					/*     dy ==> (d - (e/b))ay   (e is obj coef of y)     */
					/*eb = cols[j].obj / b; */
					EGLPNUM_TYPENAME_EGlpNumCopyFrac (eb, cols[j].obj, b);

					for (k = 0; k < rows[irow].deg; k++)
					{
						a_edge = rows[irow].adj[k];
						if (a_edge->del == 0 && a_edge != b_edge)
						{
							/*cols[a_edge->col].obj -= (eb * a_edge->coef); */
							EGLPNUM_TYPENAME_EGlpNumSubInnProdTo (cols[a_edge->col].obj, eb, a_edge->coef);
						}
					}


					/*  Delete y from graph */

					cols[j].del = 1;

					/*  Delete equation ay + bx = c */

					rows[irow].del = 1;
					for (k = 0; k < rows[irow].deg; k++)
					{
						rows[irow].adj[k]->del = 1;
					}

				}
				else
				{
					rdeg = 0;
					for (k = 0; k < rows[irow].deg && rdeg <= 2; k++)
					{
						if (rows[irow].adj[k]->del == 0)
						{
							rdeg++;
						}
					}
					if (rdeg == 2)
					{
						edge *a_edge = 0;
						int col2 = 0;

						EGLPNUM_TYPENAME_EGlpNumCopy (newub, EGLPNUM_TYPENAME_ILL_MAXDOUBLE);
						EGLPNUM_TYPENAME_EGlpNumCopy (newlb, EGLPNUM_TYPENAME_ILL_MINDOUBLE);
						EGLPNUM_TYPENAME_EGlpNumZero (a);

						/*    ay + bx = c                                */
						/*    l <= x <= u                                */
						/*      x - is column singleton                  */
						/*      derive bounds on y and substitute out x  */

						EGLPNUM_TYPENAME_EGlpNumCopy (c, rows[irow].rhs);
						EGLPNUM_TYPENAME_EGlpNumCopy (l, cols[j].lower);
						EGLPNUM_TYPENAME_EGlpNumCopy (u, cols[j].upper);

						/* Find the ay term */

						for (k = 0; k < rows[irow].deg; k++)
						{
							if (rows[irow].adj[k]->del == 0 && rows[irow].adj[k]->col != j)
							{
								a_edge = rows[irow].adj[k];
								EGLPNUM_TYPENAME_EGlpNumCopy (a, rows[irow].adj[k]->coef);
								col2 = rows[irow].adj[k]->col;
								break;
							}
						}
						if (k == rows[irow].deg)
						{
							QSlog("graph error in singleton_col");
							rval = 1;
							goto CLEANUP;
						}

						/*  Record the operation             */
						/*  x is column j,  y is column col2 */

						rval = get_next_preop (pre, &op);
						ILL_CLEANUP_IF (rval);

						op->colindex = j;
						op->rowindex = irow;
						op->ptype = ILL_PRE_DELETE_SINGLETON_VARIABLE;

						rval = grab_lp_line (G, irow, &op->line, 0);
						ILL_CLEANUP_IF (rval);
						pre->opcount++;

						/*  Adjust the bounds on y           */
						/*  Using x = c/b - (a/b)y            */
						/* we use eb as temporal variable here */
						/*if (a / b > 0) */
						EGLPNUM_TYPENAME_EGlpNumCopyFrac (eb, a, b);
						if (EGLPNUM_TYPENAME_EGlpNumIsGreatZero (eb))
						{
							/*if (l > -EGLPNUM_TYPENAME_ILL_MAXDOUBLE) */
							if (EGLPNUM_TYPENAME_EGlpNumIsLess (EGLPNUM_TYPENAME_ILL_MINDOUBLE, l))
							{
								/*newub = (c / a) - (l * b) / a; */
								EGLPNUM_TYPENAME_EGlpNumCopy (newub, c);
								EGLPNUM_TYPENAME_EGlpNumSubInnProdTo (newub, l, b);
								EGLPNUM_TYPENAME_EGlpNumDivTo (newub, a);
							}
							/*if (u < EGLPNUM_TYPENAME_ILL_MAXDOUBLE) */
							if (EGLPNUM_TYPENAME_EGlpNumIsLess (u, EGLPNUM_TYPENAME_ILL_MAXDOUBLE))
							{
								/*newlb = (c / a) - (u * b) / a; */
								EGLPNUM_TYPENAME_EGlpNumCopy (newlb, c);
								EGLPNUM_TYPENAME_EGlpNumSubInnProdTo (newlb, u, b);
								EGLPNUM_TYPENAME_EGlpNumDivTo (newlb, a);
							}
						}
						else
						{
							/*if (l > -EGLPNUM_TYPENAME_ILL_MAXDOUBLE) */
							if (EGLPNUM_TYPENAME_EGlpNumIsLess (EGLPNUM_TYPENAME_ILL_MINDOUBLE, l))
							{
								/*newlb = (c / a) - (l * b) / a; */
								EGLPNUM_TYPENAME_EGlpNumCopy (newlb, c);
								EGLPNUM_TYPENAME_EGlpNumSubInnProdTo (newlb, l, b);
								EGLPNUM_TYPENAME_EGlpNumDivTo (newlb, a);
							}
							/*if (u < EGLPNUM_TYPENAME_ILL_MAXDOUBLE) */
							if (EGLPNUM_TYPENAME_EGlpNumIsLess (u, EGLPNUM_TYPENAME_ILL_MAXDOUBLE))
							{
								/*newub = (c / a) - (u * b) / a; */
								EGLPNUM_TYPENAME_EGlpNumCopy (newub, c);
								EGLPNUM_TYPENAME_EGlpNumSubInnProdTo (newub, u, b);
								EGLPNUM_TYPENAME_EGlpNumDivTo (newub, a);
							}
						}

						if (EGLPNUM_TYPENAME_EGlpNumIsLess (cols[col2].lower, newlb))
							EGLPNUM_TYPENAME_EGlpNumCopy (cols[col2].lower, newlb);
						if (EGLPNUM_TYPENAME_EGlpNumIsLess (newub, cols[col2].upper))
							EGLPNUM_TYPENAME_EGlpNumCopy (cols[col2].upper, newub);
						EGLPNUM_TYPENAME_EGlpNumSubTo (cols[col2].obj, eb);

						/*  Delete x (and the bx term) from graph */

						cols[j].del = 1;
						b_edge->del = 1;

						/*  Delete equation ay + bx = c (and the ax term) */

						rows[irow].del = 1;
						a_edge->del = 1;
					}
				}
			}
		}
	}


CLEANUP:

	EGLPNUM_TYPENAME_EGlpNumClearVar (lb);
	EGLPNUM_TYPENAME_EGlpNumClearVar (ub);
	EGLPNUM_TYPENAME_EGlpNumClearVar (eb);
	EGLPNUM_TYPENAME_EGlpNumClearVar (b);
	EGLPNUM_TYPENAME_EGlpNumClearVar (newlb);
	EGLPNUM_TYPENAME_EGlpNumClearVar (newub);
	EGLPNUM_TYPENAME_EGlpNumClearVar (a);
	EGLPNUM_TYPENAME_EGlpNumClearVar (c);
	EGLPNUM_TYPENAME_EGlpNumClearVar (l);
	EGLPNUM_TYPENAME_EGlpNumClearVar (u);
	ILL_RETURN (rval, "singleton_columns");
}

static int duplicate_rows (
	graph * G,
	int *hit)
{
	int rval = 0;
	node *cols = G->cols;
	node *rows = G->rows;
	int ncols = G->ncols;
	int nrows = G->nrows;
	int *s = 0;
	EGLPNUM_TYPE *f = 0;
	double szeit = ILLutil_zeit ();
	EGLPNUM_TYPE q;
	int i, j, k, k2, ri, r0 = 0, n, nu = 0, got, t0, t = 1;
	node *c;

	EGLPNUM_TYPENAME_EGlpNumInitVar (q);


	/*  Code follows J. Tomlin and J. S. Welch, OR Letters 5 (1986) 7--11 */

	*hit = 0;
	if (nrows == 0)
		goto CLEANUP;

	ILL_SAFE_MALLOC (s, nrows, int);

	f = EGLPNUM_TYPENAME_EGlpNumAllocArray (nrows);

	for (i = 0; i < nrows; i++)
	{
		if (rows[i].del || rows[i].rowsense != 'E')
		{
			s[i] = EGLPNUM_TYPENAME_ILL_MAXINT;				/* EGLPNUM_TYPENAME_ILL_MAXINT means no longer eligible    */
		}
		else
		{
			s[i] = 0;									/* 0 means eligible, >0 means in a group */
			nu++;											/* Tracks the number of eligible rows    */
		}
	}

	for (j = 0; j < ncols; j++)
	{
		c = &cols[j];
		if (c->del)
			continue;
		if (c->coltype != ILL_PRE_COL_STRUC)
			continue;

		n = 0;
		t0 = t++;

		for (k = 0; k < c->deg; k++)
		{
			if (c->adj[k]->del)
				continue;

			ri = c->adj[k]->row;
			if (s[ri] == 0)
			{
				s[ri] = t0;
				EGLPNUM_TYPENAME_EGlpNumCopy (f[ri], c->adj[k]->coef);
				r0 = ri;
				n++;
			}
			else if (s[ri] < t0)
			{
				got = 0;
				for (k2 = k + 1; k2 < c->deg; k2++)
				{
					if (c->adj[k2]->del)
						continue;

					i = c->adj[k2]->row;
					if (s[i] == s[ri])
					{
						/*q = (c->adj[k]->coef * (f[i])) / (f[ri] * (c->adj[k2]->coef)); */
						EGLPNUM_TYPENAME_EGlpNumCopy (q, c->adj[k]->coef);
						EGLPNUM_TYPENAME_EGlpNumMultTo (q, f[i]);
						EGLPNUM_TYPENAME_EGlpNumDivTo (q, f[ri]);
						EGLPNUM_TYPENAME_EGlpNumDivTo (q, c->adj[k2]->coef);
						if (EGLPNUM_TYPENAME_EGlpNumIsEqual (q, EGLPNUM_TYPENAME_oneLpNum, ILL_PRE_ZERO_TOL))
						{
							s[ri] = t;
							s[i] = t;
							got++;
						}
					}
				}
				if (got)
				{
					t++;
				}
				else
				{
					s[ri] = EGLPNUM_TYPENAME_ILL_MAXINT;
					if (--nu == 0)
						goto DONE;
				}
			}
		}

		if (n == 1)
		{
			s[r0] = EGLPNUM_TYPENAME_ILL_MAXINT;
			if (--nu == 0)
				goto DONE;
		}
	}

DONE:

	{
		int idup = 0;

		for (i = 0; i < nrows; i++)
		{
			if (s[i] > 0 && s[i] < EGLPNUM_TYPENAME_ILL_MAXINT)
			{
				QSlog("Row %d: %d", i, s[i]);
				idup++;
			}
		}
		QSlog("Number of duplicate rows: %d", idup);
	}

	QSlog("Time in duplicate_rows: %.2f (seconds)", ILLutil_zeit () - szeit);

CLEANUP:

	ILL_IFFREE (s, int);

	EGLPNUM_TYPENAME_EGlpNumFreeArray (f);
	EGLPNUM_TYPENAME_EGlpNumClearVar (q);
	ILL_RETURN (rval, "duplicate_rows");
}

static int duplicate_cols (
	graph * G,
	int *hit)
{
	int rval = 0;
	node *cols = G->cols;
	node *rows = G->rows;
	int ncols = G->ncols;
	int nrows = G->nrows;
	int *s = 0;
	EGLPNUM_TYPE *f = 0;
	double szeit = ILLutil_zeit ();
	EGLPNUM_TYPE q;
	int i, j, k, k2, ci, c0 = 0, n, nu = 0, got, t0, t = 1;
	node *r;

	EGLPNUM_TYPENAME_EGlpNumInitVar (q);


	/*  Code follows J. Tomlin and J. S. Welch, OR Letters 5 (1986) 7--11 */

	*hit = 0;
	if (ncols == 0)
		goto CLEANUP;

	ILL_SAFE_MALLOC (s, ncols, int);

	f = EGLPNUM_TYPENAME_EGlpNumAllocArray (ncols);

	for (j = 0; j < ncols; j++)
	{
		if (cols[j].del || cols[j].coltype != ILL_PRE_COL_STRUC)
		{
			s[j] = EGLPNUM_TYPENAME_ILL_MAXINT;				/* EGLPNUM_TYPENAME_ILL_MAXINT means no longer eligible    */
		}
		else
		{
			s[j] = 0;									/* 0 means eligible, >0 means in a group */
			nu++;											/* Tracks the number of eligible rows    */
		}
	}

	for (i = 0; i < nrows; i++)
	{
		r = &rows[i];
		if (r->del)
			continue;

		n = 0;
		t0 = t++;

		for (k = 0; k < r->deg; k++)
		{
			if (r->adj[k]->del)
				continue;

			ci = r->adj[k]->col;
			if (s[ci] == 0)
			{
				s[ci] = t0;
				EGLPNUM_TYPENAME_EGlpNumCopy (f[ci], r->adj[k]->coef);
				c0 = ci;
				n++;
			}
			else if (s[ci] < t0)
			{
				got = 0;
				for (k2 = k + 1; k2 < r->deg; k2++)
				{
					if (r->adj[k2]->del)
						continue;

					j = r->adj[k2]->col;
					if (s[j] == s[ci])
					{
						/*q = (r->adj[k]->coef * (f[j])) / (f[ci] * (r->adj[k2]->coef)); */
						EGLPNUM_TYPENAME_EGlpNumCopy (q, r->adj[k]->coef);
						EGLPNUM_TYPENAME_EGlpNumMultTo (q, f[j]);
						EGLPNUM_TYPENAME_EGlpNumDivTo (q, f[ci]);
						EGLPNUM_TYPENAME_EGlpNumDivTo (q, r->adj[k2]->coef);
						if (EGLPNUM_TYPENAME_EGlpNumIsEqual (q, EGLPNUM_TYPENAME_oneLpNum, ILL_PRE_ZERO_TOL))
						{
							s[ci] = t;
							s[j] = t;
							got++;
						}
					}
				}
				if (got)
				{
					t++;
				}
				else
				{
					s[ci] = EGLPNUM_TYPENAME_ILL_MAXINT;
					if (--nu == 0)
						goto DONE;
				}
			}
		}

		if (n == 1)
		{
			s[c0] = EGLPNUM_TYPENAME_ILL_MAXINT;
			if (--nu == 0)
				goto DONE;
		}
	}

DONE:

	{
		int dcount;
		int *dcnt;
		int *dlist;

		rval = gather_dup_lists (s, ncols, &dcount, &dcnt, &dlist);
		ILL_CLEANUP_IF (rval);
	}

	QSlog("Time in duplicate_cols: %.2f (seconds)", ILLutil_zeit () - szeit);

CLEANUP:

	ILL_IFFREE (s, int);

	EGLPNUM_TYPENAME_EGlpNumFreeArray (f);
	EGLPNUM_TYPENAME_EGlpNumClearVar (q);
	ILL_RETURN (rval, "duplicate_cols");
}

static int gather_dup_lists (
	/* graph *G, */ int *s,
	/* double *f, */
	int count,
	int *duptotal,
	int **dupcnt,
	int **dupind)
{
	int rval = 0;
	int *cnt = 0;
	int *ind = 0;
	int *beg = 0;
	int i, smax = 0, ndup = 0, total = 0;

	*duptotal = 0;
	*dupcnt = 0;
	*dupind = 0;


	for (i = 0; i < count; i++)
	{
		if (s[i] < EGLPNUM_TYPENAME_ILL_MAXINT && s[i] > smax)
			smax = s[i];
	}
	if (smax == 0)
		goto CLEANUP;

	ILL_SAFE_MALLOC (cnt, smax + 1, int);

	ILL_SAFE_MALLOC (ind, smax + 1, int);

	for (i = 0; i < smax + 1; i++)
	{
		cnt[i] = 0;
	}

	for (i = 0; i < count; i++)
	{
		if (s[i] < EGLPNUM_TYPENAME_ILL_MAXINT)
		{
			cnt[s[i]]++;
		}
	}

	if (cnt[0] > 0)
		QSlog("%d Empty Lines", cnt[0]);

	QSlog("Duplicate Classes:");
	for (i = 1; i < smax + 1; i++)
	{
		if (cnt[i] > 1)
		{
			ndup++;
			QSlog(" %d", cnt[i]);
		}
	}
	QSlog("  Number %d\n", ndup);

	if (ndup == 0)
		goto CLEANUP;

	ILL_SAFE_MALLOC (beg, ndup, int);

	for (i = 1, ndup = 0; i < smax + 1; i++)
	{
		if (cnt[i] > 1)
		{
			beg[ndup] = total;
			total += cnt[i];
			ind[i] = ndup;
			ndup++;
		}
	}

	if (total == 0)
		goto CLEANUP;

	ILL_SAFE_MALLOC (*dupcnt, ndup, int);

	ILL_SAFE_MALLOC (*dupind, total, int);

	for (i = 0; i < ndup; i++)
	{
		(*dupcnt)[i] = 0;
	}

	for (i = 0; i < count; i++)
	{
		if (s[i] < EGLPNUM_TYPENAME_ILL_MAXINT && s[i] > 0)
		{
			if (cnt[s[i]] > 1)
			{
				(*dupind)[beg[ind[s[i]]] + (*dupcnt)[ind[s[i]]]] = i;
				(*dupcnt)[ind[s[i]]]++;
			}
		}
	}

	for (i = 0; i < ndup; i++)
	{
		int j;

		for (j = beg[i]; j < beg[i] + (*dupcnt)[i]; j++)
		{
			QSlog(" %d", (*dupind)[j]);
		}
		QSlog(" | ");
	}

	*duptotal = ndup;

CLEANUP:

	ILL_IFFREE (cnt, int);
	ILL_IFFREE (ind, int);
	ILL_IFFREE (beg, int);

	ILL_RETURN (rval, "gather_dup_lists");
}

static void set_fixed_variable (
	graph * G,
	int j,
	EGLPNUM_TYPE val)
{
	int k;
	edge *e;

	G->cols[j].del = 1;
	for (k = 0; k < G->cols[j].deg; k++)
	{
		e = G->cols[j].adj[k];
		if (e->del == 0)
		{
			/*G->rows[e->row].rhs -= (e->coef * val); */
			EGLPNUM_TYPENAME_EGlpNumSubInnProdTo (G->rows[e->row].rhs, e->coef, val);
			e->del = 1;
		}
	}
}

static void get_implied_rhs_bounds (
	graph * G,
	int i,
	EGLPNUM_TYPE * lb,
	EGLPNUM_TYPE * ub)
{
	int k;
	EGLPNUM_TYPE l, u;
	node *cols = G->cols;
	node *rows = G->rows;
	edge *e;

	EGLPNUM_TYPENAME_EGlpNumInitVar (u);
	EGLPNUM_TYPENAME_EGlpNumInitVar (l);

	EGLPNUM_TYPENAME_EGlpNumZero (l);
	for (k = 0; k < rows[i].deg; k++)
	{
		e = rows[i].adj[k];
		if (e->del == 0)
		{
			if (EGLPNUM_TYPENAME_EGlpNumIsLessZero (e->coef))
			{
				if (EGLPNUM_TYPENAME_EGlpNumIsEqqual (cols[e->col].upper, EGLPNUM_TYPENAME_ILL_MAXDOUBLE))
				{
					EGLPNUM_TYPENAME_EGlpNumCopy (l, EGLPNUM_TYPENAME_ILL_MINDOUBLE);
					break;
				}
				else
				{
					/*l += (e->coef * cols[e->col].upper); */
					EGLPNUM_TYPENAME_EGlpNumAddInnProdTo (l, e->coef, cols[e->col].upper);
				}
			}
			else if (EGLPNUM_TYPENAME_EGlpNumIsGreatZero (e->coef))
			{
				if (EGLPNUM_TYPENAME_EGlpNumIsEqqual (cols[e->col].lower, EGLPNUM_TYPENAME_ILL_MINDOUBLE))
				{
					EGLPNUM_TYPENAME_EGlpNumCopy (l, EGLPNUM_TYPENAME_ILL_MINDOUBLE);
					break;
				}
				else
				{
					/*l += (e->coef * cols[e->col].lower); */
					EGLPNUM_TYPENAME_EGlpNumAddInnProdTo (l, e->coef, cols[e->col].lower);
				}
			}
		}
	}

	EGLPNUM_TYPENAME_EGlpNumZero (u);
	for (k = 0; k < rows[i].deg; k++)
	{
		e = rows[i].adj[k];
		if (e->del == 0)
		{
			if (EGLPNUM_TYPENAME_EGlpNumIsLessZero (e->coef ))
			{
				if (EGLPNUM_TYPENAME_EGlpNumIsEqqual (cols[e->col].lower, EGLPNUM_TYPENAME_ILL_MINDOUBLE))
				{
					EGLPNUM_TYPENAME_EGlpNumCopy (u, EGLPNUM_TYPENAME_ILL_MAXDOUBLE);
				}
				else
				{
					/*u += (e->coef * cols[e->col].lower); */
					EGLPNUM_TYPENAME_EGlpNumAddInnProdTo (u, e->coef, cols[e->col].lower);
				}
			}
			else if (EGLPNUM_TYPENAME_EGlpNumIsGreatZero (e->coef))
			{
				if (EGLPNUM_TYPENAME_EGlpNumIsEqqual (cols[e->col].upper, EGLPNUM_TYPENAME_ILL_MAXDOUBLE))
				{
					EGLPNUM_TYPENAME_EGlpNumCopy (u, EGLPNUM_TYPENAME_ILL_MAXDOUBLE);
				}
				else
				{
					/*u += (e->coef * cols[e->col].upper); */
					EGLPNUM_TYPENAME_EGlpNumAddInnProdTo (u, e->coef, cols[e->col].upper);
				}
			}
		}
	}

	EGLPNUM_TYPENAME_EGlpNumCopy (*lb, l);
	EGLPNUM_TYPENAME_EGlpNumCopy (*ub, u);
	EGLPNUM_TYPENAME_EGlpNumClearVar (u);
	EGLPNUM_TYPENAME_EGlpNumClearVar (l);
}

static void get_implied_variable_bounds (
	graph * G,
	int j,
	edge * a_ij,
	EGLPNUM_TYPE * lb,
	EGLPNUM_TYPE * ub)
{
	int i = a_ij->row;
	EGLPNUM_TYPE l, u;

	EGLPNUM_TYPENAME_EGlpNumInitVar (u);
	EGLPNUM_TYPENAME_EGlpNumInitVar (l);

	get_implied_rhs_bounds (G, i, &l, &u);
	EGLPNUM_TYPENAME_EGlpNumCopy (*lb, EGLPNUM_TYPENAME_ILL_MINDOUBLE);
	EGLPNUM_TYPENAME_EGlpNumCopy (*ub, EGLPNUM_TYPENAME_ILL_MAXDOUBLE);

	if (EGLPNUM_TYPENAME_EGlpNumIsLess (ILL_PRE_FEAS_TOL, a_ij->coef))
	{
		if (EGLPNUM_TYPENAME_EGlpNumIsLess (u, EGLPNUM_TYPENAME_ILL_MAXDOUBLE))
		{
			/**lb = (G->rows[i].rhs - u) / a_ij->coef + G->cols[j].upper;*/
			EGLPNUM_TYPENAME_EGlpNumCopyDiffRatio (*lb, G->rows[i].rhs, u, a_ij->coef);
			EGLPNUM_TYPENAME_EGlpNumAddTo (*lb, G->cols[j].upper);
		}
		if (EGLPNUM_TYPENAME_EGlpNumIsLess (EGLPNUM_TYPENAME_ILL_MINDOUBLE, l))
		{
			/**ub = (G->rows[i].rhs - l) / a_ij->coef + G->cols[j].lower;*/
			EGLPNUM_TYPENAME_EGlpNumCopyDiffRatio (*ub, G->rows[i].rhs, l, a_ij->coef);
			EGLPNUM_TYPENAME_EGlpNumAddTo (*ub, G->cols[j].lower);
		}
	}
	else if (EGLPNUM_TYPENAME_EGlpNumIsLess (a_ij->coef, ILL_PRE_FEAS_TOL))
	{
		if (EGLPNUM_TYPENAME_EGlpNumIsLess (EGLPNUM_TYPENAME_ILL_MINDOUBLE, l))
		{
			/**lb = (G->rows[i].rhs - l) / a_ij->coef + G->cols[j].upper;*/
			EGLPNUM_TYPENAME_EGlpNumCopyDiffRatio (*lb, G->rows[i].rhs, l, a_ij->coef);
			EGLPNUM_TYPENAME_EGlpNumAddTo (*lb, G->cols[j].upper);
		}
		if (EGLPNUM_TYPENAME_EGlpNumIsLess (u, EGLPNUM_TYPENAME_ILL_MAXDOUBLE))
		{
			/**ub = (G->rows[i].rhs - u) / a_ij->coef + G->cols[j].lower;*/
			EGLPNUM_TYPENAME_EGlpNumCopyDiffRatio (*ub, G->rows[i].rhs, u, a_ij->coef);
			EGLPNUM_TYPENAME_EGlpNumAddTo (*ub, G->cols[j].lower);
		}
	}
	EGLPNUM_TYPENAME_EGlpNumClearVar (u);
	EGLPNUM_TYPENAME_EGlpNumClearVar (l);
}

static int get_next_preop (
	EGLPNUM_TYPENAME_ILLlp_predata * pre,
	EGLPNUM_TYPENAME_ILLlp_preop ** op)
{
	int rval = 0;

	if (pre->opcount >= pre->opsize)
	{
		pre->opsize *= 1.3;
		pre->opsize += 1000;
		if (pre->opsize < pre->opcount + 1)
			pre->opsize = pre->opcount + 1;
		pre->oplist = EGrealloc (pre->oplist, sizeof (EGLPNUM_TYPENAME_ILLlp_preop) * pre->opsize);
		//rval = ILLutil_reallocrus_scale ((void **) &pre->oplist,
		//                                 &pre->opsize, pre->opcount + 1, 1.3,
		//                                 sizeof (EGLPNUM_TYPENAME_ILLlp_preop));
		//ILL_CLEANUP_IF (rval);
	}
	*op = &pre->oplist[pre->opcount];
	EGLPNUM_TYPENAME_ILLlp_preop_init (*op);

//CLEANUP:

	ILL_RETURN (rval, "get_next_preop");
}

static int add_to_list (
	ILLptrworld * world,
	intptr ** list,
	int i)
{
	int rval = 0;
	intptr *ip;

	ip = intptralloc (world);
	if (!ip)
	{
		rval = 1;
		goto CLEANUP;
	}
	ip->this_val = i;
	ip->next = *list;
	*list = ip;

CLEANUP:

	ILL_RETURN (rval, "add_to_list");
}

static int build_graph (
	EGLPNUM_TYPENAME_ILLlpdata * lp,
	graph * G)
{
	int rval = 0;
	int ncols = lp->ncols;
	int nrows = lp->nrows;
	int nzcount = lp->nzcount;
	int i, j, k, stop, count;
	edge *edgelist;
	node *rows, *cols;
	EGLPNUM_TYPENAME_ILLmatrix *A = &lp->A;

	G->objsense = lp->objsense;

	ILL_SAFE_MALLOC (G->rows, nrows, node);
	if (!G->rows)
	{
		QSlog("out of memory in build_graph");
		rval = 1;
		goto CLEANUP;
	}
	rows = G->rows;

	for (i = 0; i < nrows; i++)
	{
		rows[i].rowsense = lp->sense[i];
		rows[i].deg = 0;
	}

	ILL_SAFE_MALLOC (G->cols, ncols, node);
	ILL_SAFE_MALLOC (G->edgelist, nzcount, edge);
	for (i = nzcount; i--;)
		EGLPNUM_TYPENAME_EGlpNumInitVar ((G->edgelist[i].coef));
	G->nzcount = nzcount;
	ILL_SAFE_MALLOC (G->adjspace, 2 * nzcount, edge *);

	if (!G->cols || !G->edgelist || !G->adjspace)
	{
		QSlog("out of memory in build_graph");
		rval = 1;
		goto CLEANUP;
	}

	cols = G->cols;
	edgelist = G->edgelist;

	for (j = 0; j < ncols; j++)
	{
		stop = A->matbeg[j] + A->matcnt[j];
		for (k = A->matbeg[j]; k < stop; k++)
		{
			rows[A->matind[k]].deg++;
		}
	}

	for (i = 0, count = 0; i < nrows; i++)
	{
		rows[i].adj = G->adjspace + count;
		count += rows[i].deg;
		rows[i].deg = 0;
	}

	for (j = 0; j < ncols; j++)
	{
		cols[j].adj = G->adjspace + count;
		count += A->matcnt[j];
		cols[j].deg = 0;
		cols[j].coltype = ILL_PRE_COL_STRUC;
	}
	for (i = 0; i < nrows; i++)
	{
		cols[lp->rowmap[i]].coltype = ILL_PRE_COL_LOGICAL;
	}

	for (j = 0, count = 0; j < ncols; j++)
	{
		EGLPNUM_TYPENAME_EGlpNumCopy (cols[j].obj, lp->obj[j]);
		EGLPNUM_TYPENAME_EGlpNumCopy (cols[j].lower, lp->lower[j]);
		EGLPNUM_TYPENAME_EGlpNumCopy (cols[j].upper, lp->upper[j]);
		stop = A->matbeg[j] + A->matcnt[j];
		for (k = A->matbeg[j]; k < stop; k++)
		{
			i = A->matind[k];
			rows[i].adj[rows[i].deg++] = &(edgelist[count]);
			cols[j].adj[cols[j].deg++] = &(edgelist[count]);
			edgelist[count].row = i;
			edgelist[count].col = j;
			EGLPNUM_TYPENAME_EGlpNumCopy (edgelist[count].coef, A->matval[k]);
			edgelist[count].mark = 0;
			edgelist[count].del = 0;
			edgelist[count].coltype = cols[j].coltype;
			count++;
		}
	}
	if (count != nzcount)
	{
		QSlog("counts are off in build_graph");
		rval = 1;
		goto CLEANUP;
	}

	G->ecount = count;
	G->nrows = nrows;
	G->ncols = ncols;

	for (i = 0; i < G->nrows; i++)
	{
		G->rows[i].del = 0;
		EGLPNUM_TYPENAME_EGlpNumCopy (G->rows[i].rhs, lp->rhs[i]);
	}
	for (j = 0; j < G->ncols; j++)
	{
		G->cols[j].del = 0;
	}

CLEANUP:

	ILL_RETURN (rval, "build_graph");
}

static void dump_graph (
	graph * G)
{
	int i, j, k;

	QSlog("ecount = %d, nrows = %d, ncols = %d",
							G->ecount, G->nrows, G->ncols);

	for (i = 0; i < G->nrows; i++)
	{
		QSlog("Row %d:", i);
		for (k = 0; k < G->rows[i].deg; k++)
		{
			QSlog(" %d", G->rows[i].adj[k]->col);
			if (G->rows[i].adj[k]->coltype == ILL_PRE_COL_LOGICAL)
				QSlog("S");
			QSlog("(%g)", EGLPNUM_TYPENAME_EGlpNumToLf (G->rows[i].adj[k]->coef));
		}
		QSlog("  rhs: %g", EGLPNUM_TYPENAME_EGlpNumToLf (G->rows[i].rhs));
		if (G->rows[i].del) QSlog(" (deleted)");
	}

	for (j = 0; j < G->ncols; j++)
	{
		if (G->cols[j].coltype == ILL_PRE_COL_LOGICAL)
		{
			QSlog("Slk %d:", j);
		}
		else
		{
			QSlog("Col %d:", j);
		}
		for (k = 0; k < G->cols[j].deg; k++)
		{
			QSlog(" %d", G->cols[j].adj[k]->row);
		}
		QSlog("  obj: %g  bnd: (%g, %g)", EGLPNUM_TYPENAME_EGlpNumToLf (G->cols[j].obj),
								EGLPNUM_TYPENAME_EGlpNumToLf (G->cols[j].lower), EGLPNUM_TYPENAME_EGlpNumToLf (G->cols[j].upper));
		if (G->cols[j].del) QSlog(" (deleted)");
	}
}

static void init_graph (
	graph * G)
{
	if (G)
	{
		G->edgelist = 0;
		G->rows = 0;
		G->cols = 0;
		G->ecount = 0;
		G->nrows = 0;
		G->ncols = 0;
		G->adjspace = 0;
		ILLptrworld_init (&G->intptrworld);
	}
}

static void free_graph (
	graph * G)
{
	register int i;

	if (G)
	{
		int total, onlist;

		for (i = G->nzcount; i--;)
			EGLPNUM_TYPENAME_EGlpNumClearVar ((G->edgelist[i].coef));
		ILL_IFFREE (G->edgelist, edge);
		ILL_IFFREE (G->rows, node);
		ILL_IFFREE (G->cols, node);
		ILL_IFFREE (G->adjspace, edge *);
		if (intptr_check_leaks (&G->intptrworld, &total, &onlist))
		{
			QSlog("WARNING: %d outstanding intptrs", total - onlist);
		}
		ILLptrworld_delete (&G->intptrworld);
		init_graph (G);
	}
}

int EGLPNUM_TYPENAME_ILLlp_sinfo_print (
	EGLPNUM_TYPENAME_ILLlp_sinfo * s)
{
	int rval = 0;
	int i;
	EGLPNUM_TYPENAME_ILLlpdata lp;
	char *sense = 0;

	EGLPNUM_TYPENAME_ILLlpdata_init (&lp);

	lp.nrows = s->nrows;
	lp.ncols = s->ncols;
	lp.nzcount = s->nzcount;
	lp.objsense = s->objsense;
	lp.obj = s->obj;
	lp.rhs = s->rhs;
	lp.lower = s->lower;
	lp.upper = s->upper;
	lp.A.matval = s->A.matval;
	lp.A.matcnt = s->A.matcnt;
	lp.A.matbeg = s->A.matbeg;
	lp.A.matind = s->A.matind;
	lp.rownames = 0;
	lp.colnames = s->colnames;
	lp.objname = 0;
	lp.probname = 0;
	lp.intmarker = 0;

	ILL_SAFE_MALLOC (sense, s->nrows, char);

	if (!sense)
	{
		QSlog("out of memory in EGLPNUM_TYPENAME_ILLlp_sinfo_print");
		rval = 1;
		goto CLEANUP;
	}
	for (i = 0; i < s->nrows; i++)
	{
		sense[i] = 'E';
	}
	lp.sense = sense;

/*
    rval = ILLlpdata_writelp (&lp, 0);
    ILL_CLEANUP_IF (rval);
*/

CLEANUP:

	ILL_IFFREE (sense, char);

	ILL_RETURN (rval, "EGLPNUM_TYPENAME_ILLlp_sinfo_print");
}

void EGLPNUM_TYPENAME_ILLlp_sinfo_init (
	EGLPNUM_TYPENAME_ILLlp_sinfo * sinfo)
{
	if (sinfo)
	{
		sinfo->ncols = 0;
		sinfo->nrows = 0;
		sinfo->nzcount = 0;
		sinfo->rowsize = 0;
		sinfo->colsize = 0;
		sinfo->obj = 0;
		sinfo->rhs = 0;
		sinfo->lower = 0;
		sinfo->upper = 0;
		sinfo->colnames = 0;
		sinfo->objsense = EGLPNUM_TYPENAME_ILL_MIN;
		EGLPNUM_TYPENAME_ILLmatrix_init (&sinfo->A);
	}
}

void EGLPNUM_TYPENAME_ILLlp_sinfo_free (
	EGLPNUM_TYPENAME_ILLlp_sinfo * sinfo)
{
	if (sinfo)
	{
		EGLPNUM_TYPENAME_EGlpNumFreeArray (sinfo->obj);
		EGLPNUM_TYPENAME_EGlpNumFreeArray (sinfo->lower);
		EGLPNUM_TYPENAME_EGlpNumFreeArray (sinfo->upper);
		EGLPNUM_TYPENAME_EGlpNumFreeArray (sinfo->rhs);
		EGLPNUM_TYPENAME_ILLmatrix_free (&sinfo->A);
		if (sinfo->colnames)
		{
			int i;

			for (i = 0; i < sinfo->ncols; i++)
			{
				ILL_IFFREE (sinfo->colnames[i], char);
			}
			ILL_IFFREE (sinfo->colnames, char *);
		}
		EGLPNUM_TYPENAME_ILLlp_sinfo_init (sinfo);
	}
}

void EGLPNUM_TYPENAME_ILLlp_predata_init (
	EGLPNUM_TYPENAME_ILLlp_predata * pre)
{
	if (pre)
	{
		pre->opcount = 0;
		pre->opsize = 0;
		pre->oplist = 0;
		pre->r_nrows = 0;
		pre->r_ncols = 0;
		pre->colmap = 0;
		pre->rowmap = 0;
		pre->colscale = 0;
		pre->rowscale = 0;
		pre->colfixval = 0;
		pre->rowfixval = 0;
	}
}

void EGLPNUM_TYPENAME_ILLlp_predata_free (
	EGLPNUM_TYPENAME_ILLlp_predata * pre)
{
	if (pre)
	{
		int i;

		for (i = 0; i < pre->opcount; i++)
		{
			EGLPNUM_TYPENAME_ILLlp_preop_free (&pre->oplist[i]);
		}
		ILL_IFFREE (pre->oplist, EGLPNUM_TYPENAME_ILLlp_preop);
		ILL_IFFREE (pre->colmap, int);
		ILL_IFFREE (pre->rowmap, int);

		ILL_IFFREE (pre->colscale, EGLPNUM_TYPE);
		ILL_IFFREE (pre->rowscale, EGLPNUM_TYPE);
		ILL_IFFREE (pre->colfixval, EGLPNUM_TYPE);
		ILL_IFFREE (pre->rowfixval, EGLPNUM_TYPE);
		EGLPNUM_TYPENAME_ILLlp_predata_init (pre);
	}
}

void EGLPNUM_TYPENAME_ILLlp_preop_init (
	EGLPNUM_TYPENAME_ILLlp_preop * op)
{
	if (op)
	{
		op->ptype = 0;
		op->rowindex = -1;
		op->colindex = -1;
		EGLPNUM_TYPENAME_ILLlp_preline_init (&op->line);
	}
}

void EGLPNUM_TYPENAME_ILLlp_preop_free (
	EGLPNUM_TYPENAME_ILLlp_preop * op)
{
	if (op)
	{
		EGLPNUM_TYPENAME_ILLlp_preline_free (&op->line);
		EGLPNUM_TYPENAME_ILLlp_preop_init (op);
	}
}

void EGLPNUM_TYPENAME_ILLlp_preline_init (
	EGLPNUM_TYPENAME_ILLlp_preline * line)
{
	if (line)
	{
		EGLPNUM_TYPENAME_EGlpNumInitVar (line->rhs);
		EGLPNUM_TYPENAME_EGlpNumInitVar (line->obj);
		EGLPNUM_TYPENAME_EGlpNumInitVar (line->upper);
		EGLPNUM_TYPENAME_EGlpNumInitVar (line->lower);
		EGLPNUM_TYPENAME_EGlpNumZero (line->rhs);
		EGLPNUM_TYPENAME_EGlpNumZero (line->obj);
		EGLPNUM_TYPENAME_EGlpNumZero (line->upper);
		EGLPNUM_TYPENAME_EGlpNumZero (line->lower);
		line->count = 0;
		line->ind = 0;
		line->val = 0;
	}
}

void EGLPNUM_TYPENAME_ILLlp_preline_free (
	EGLPNUM_TYPENAME_ILLlp_preline * line)
{
	if (line)
	{
		EGLPNUM_TYPENAME_EGlpNumClearVar (line->rhs);
		EGLPNUM_TYPENAME_EGlpNumClearVar (line->obj);
		EGLPNUM_TYPENAME_EGlpNumClearVar (line->upper);
		EGLPNUM_TYPENAME_EGlpNumClearVar (line->lower);
		ILL_IFFREE (line->ind, int);

		EGLPNUM_TYPENAME_EGlpNumFreeArray (line->val);
		//EGLPNUM_TYPENAME_ILLlp_preline_init (line);
	}
}
