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
/*    int ILLlp_add_logicals (ILLlpata *lp)                                 */
/*    int ILLlp_presolve (ILLlpdata *lp)                                    */
/*    int ILLlp_scale (ILLlpdata *lp)                                       */
/*    void ILLlp_sinfo_init (ILLlp_sinfo *sinfo)                            */
/*    void ILLlp_sinfo_free (ILLlp_sinfo *sinfo)                            */
/*    void ILLlp_predata_init (ILLlp_predata *pre)                          */
/*    void ILLlp_predata_free (ILLlp_predata *pre)                          */
/*                                                                          */
/*  NOTES                                                                   */
/*                                                                          */
/*    presolve will assume that logicals have been added.                   */
/*                                                                          */
/****************************************************************************/

#include "presolve.h"

#include <stdlib.h>
#include <string.h>

#include "qs_config.h"

#include "eg_lpnum.h"
#include "eg_io.h"

#include "iqsutil.h"
#include "lpdata.h"
#include "lpdefs.h"
//extern EGlpNum_t SZERO_TOLER;

#define ILL_LP_STATUS_OK (0)
#define ILL_PRE_FEAS_TOL PFEAS_TOLER	//(1e-6)
#define ILL_PRE_ZERO_TOL PIVOT_TOLER	//(1e-10)

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
	EGlpNum_t coef;
}
edge;

typedef struct node
{
	edge **adj;
	EGlpNum_t obj;
	EGlpNum_t lower;
	EGlpNum_t upper;
	EGlpNum_t rhs;
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
	EGlpNum_t val),
  get_implied_rhs_bounds (
	graph * G,
	int i,
	EGlpNum_t * lb,
	EGlpNum_t * ub),
  get_implied_variable_bounds (
	graph * G,
	int j,
	edge * a_ij,
	EGlpNum_t * lb,
	EGlpNum_t * ub),
  dump_line (
	ILLlp_preline * line),
  init_graph (
	graph * G),
  free_graph (
	graph * G),
  dump_graph (
	graph * G);

static int simple_presolve (
	ILLlpdata * lp,
	ILLlp_predata * pre,
	ILLlp_sinfo * info,
	int pre_types,
	int *status),
  grab_lp_line (
	graph * G,
	int indx,
	ILLlp_preline * line,
	int row_or_col),
  grab_lp_info (
	graph * G,
	char **colnames,
	ILLlp_sinfo * info),
  fixed_variables (
	graph * G,
	ILLlp_predata * pre),
  empty_columns (
	graph * G,
	ILLlp_predata * pre),
  singleton_rows (
	graph * G,
	ILLlp_predata * pre,
	int *hit),
  forcing_constraints (
	graph * G,
	ILLlp_predata * pre,
	int *hit),
  singleton_columns (
	graph * G,
	ILLlp_predata * pre,
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
	ILLlp_predata * pre,
	ILLlp_preop ** op),
  add_to_list (
	ILLptrworld * world,
	intptr ** list,
	int i),
  build_graph (
	ILLlpdata * lp,
	graph * G);


ILL_PTRWORLD_ROUTINES (intptr, intptralloc, intptr_bulkalloc, intptrfree)
ILL_PTRWORLD_LISTFREE_ROUTINE (intptr, intptr_listfree, intptrfree)
ILL_PTRWORLD_LEAKS_ROUTINE (intptr, intptr_check_leaks, this_val, int)
		 int ILLlp_add_logicals (
	ILLlpdata * lp)
{
	int rval = 0;
	int ncols, nrows, nzcount, i, aindex;
	char *sense;
	ILLmatrix *A;

	if (!lp)
	{
		fprintf (stderr, "ILLlp_add_logicals called with a NULL pointer\n");
		rval = 1;
		goto CLEANUP;
	}

	printf ("ILLlp_add_logicals ...\n");
	fflush (stdout);

	A = &lp->A;
	sense = lp->sense;
	ncols = lp->ncols;
	nrows = lp->nrows;
	nzcount = lp->nzcount;

	if (nrows == 0)
		goto CLEANUP;
	EGlpNumReallocArray (&(lp->obj), lp->colsize + nrows);
	EGlpNumReallocArray (&(lp->upper), lp->colsize + nrows);
	EGlpNumReallocArray (&(lp->lower), lp->colsize + nrows);
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
	EGlpNumReallocArray (&(A->matval), A->matsize + nrows);

	for (i = 0; i < nrows; i++)
	{
		A->matind[A->matsize + i] = -1;
	}

	aindex = A->matsize - A->matfree;

	for (i = 0; i < nrows; i++)
	{
		lp->rowmap[i] = ncols;
		EGlpNumZero (lp->obj[ncols]);
		A->matcnt[ncols] = 1;
		A->matbeg[ncols] = aindex;
		A->matind[aindex] = i;
		switch (sense[i])
		{
		case 'E':									/* Arificial */
			EGlpNumZero (lp->lower[ncols]);
			EGlpNumZero (lp->upper[ncols]);
			EGlpNumOne (A->matval[aindex]);
			break;
		case 'G':									/* Surplus   */
			EGlpNumZero (lp->lower[ncols]);
			EGlpNumCopy (lp->upper[ncols], ILL_MAXDOUBLE);
			EGlpNumOne (A->matval[aindex]);
			EGlpNumSign (A->matval[aindex]);
			break;
		case 'L':									/* Slack     */
			EGlpNumZero (lp->lower[ncols]);
			EGlpNumCopy (lp->upper[ncols], ILL_MAXDOUBLE);
			EGlpNumOne (A->matval[aindex]);
			break;
		case 'R':									/* Range     */
			EGlpNumZero (lp->lower[ncols]);
			EGlpNumCopy (lp->upper[ncols], lp->rangeval[i]);
			EGlpNumOne (A->matval[aindex]);
			EGlpNumSign (A->matval[aindex]);
			break;
		default:
			fprintf (stderr, "unknown sense %c in ILLlp_add_logicals\n", sense[i]);
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
		ILLlp_rows_clear (lp->rA);
	}
	else
	{
		ILL_SAFE_MALLOC (lp->rA, 1, ILLlp_rows);
	}

	rval = ILLlp_rows_init (lp->rA, lp, 1);
	ILL_CLEANUP_IF (rval);

CLEANUP:

	ILL_RETURN (rval, "ILLlp_add_logicals");
}

int ILLlp_scale (
	ILLlpdata * lp)
{
	int rval = 0;
	int i, j, k, col, row, nstruct, start, stop;
	ILLmatrix *A;
	EGlpNum_t rho;
	EGlpNum_t *gama = 0;

	EGlpNumInitVar (rho);

	/* Columns - divide by largest absolute value */

	if (!lp)
	{
		ILL_ERROR (rval, "ILLlp_scale called with a NULL pointer");
	}

	if (lp->nrows == 0 || lp->ncols == 0)
		goto CLEANUP;

	A = &lp->A;
	nstruct = lp->nstruct;

	for (j = 0; j < nstruct; j++)
	{
		col = lp->structmap[j];
		EGlpNumZero (rho);

		start = A->matbeg[col];
		stop = start + A->matcnt[col];

		for (k = start; k < stop; k++)
		{
			EGlpNumSetToMaxAbs (rho, A->matval[k]);
		}

		if (EGlpNumIsGreatZero (rho))
		{
			for (k = start; k < stop; k++)
			{
				EGlpNumDivTo (A->matval[k], rho);
			}
			EGlpNumDivTo (lp->obj[col], rho);
			if (EGlpNumIsNeqq (lp->lower[col], ILL_MINDOUBLE))
				EGlpNumMultTo (lp->lower[col], rho);
			if (EGlpNumIsNeqq (lp->upper[col], ILL_MAXDOUBLE))
				EGlpNumMultTo (lp->upper[col], rho);
		}
	}

	gama = EGlpNumAllocArray (lp->nrows);
	for (i = 0; i < lp->nrows; i++)
	{
		EGlpNumZero (gama[i]);
	}

	for (j = 0; j < nstruct; j++)
	{
		col = lp->structmap[j];
		start = A->matbeg[col];
		stop = start + A->matcnt[col];

		for (k = start; k < stop; k++)
		{
			row = A->matind[k];
			EGlpNumSetToMaxAbs (gama[row], A->matval[k]);
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
			if (EGlpNumIsGreatZero (gama[row]))
			{
				EGlpNumDivTo (A->matval[k], gama[row]);
			}
		}
	}

	for (i = 0; i < lp->nrows; i++)
	{
		if (EGlpNumIsGreatZero ( gama[i]))
		{
			EGlpNumDivTo (lp->rhs[i], gama[i]);
			col = lp->rowmap[i];
			if (EGlpNumIsNeqq (lp->upper[col], ILL_MAXDOUBLE))
			{
				EGlpNumDivTo (lp->upper[col], gama[i]);	/* Ranged row */
			}
		}
	}

	if (lp->rA)
	{															/* Need to clear the row version of data */
		ILLlp_rows_clear (lp->rA);
		ILL_IFFREE (lp->rA, ILLlp_rows);
	}


CLEANUP:

	EGlpNumClearVar (rho);
	EGlpNumFreeArray (gama);
	ILL_RETURN (rval, "ILLlp_scale");
}

int ILLlp_presolve (
	ILLlpdata * lp,
	int pre_types)
{
	int rval = 0;
	int status = ILL_LP_STATUS_OK;
	ILLlp_predata *pre = 0;
	ILLlp_sinfo *info = 0;

	if (!lp)
	{
		fprintf (stderr, "ILLlp_presolve called with a NULL pointer\n");
		rval = 1;
		goto CLEANUP;
	}


/*
    ILLlpdata_writelp (lp, 0);
    printf ("\n"); fflush (stdout);
*/

	ILL_SAFE_MALLOC (pre, 1, ILLlp_predata);
	ILLlp_predata_init (pre);

	ILL_SAFE_MALLOC (info, 1, ILLlp_sinfo);
	ILLlp_sinfo_init (info);

	rval = simple_presolve (lp, pre, info, pre_types, &status);
	ILL_CLEANUP_IF (rval);
	if (status != ILL_LP_STATUS_OK)
	{
		printf ("simple_presolve returned with bad status\n");
		rval = 1;
		goto CLEANUP;
	}

/*
    rval = ILLlp_sinfo_print (info);
    ILL_CLEANUP_IF (rval);
*/

CLEANUP:

	if (rval)
	{
		if (pre)
		{
			ILLlp_predata_free (pre);
			ILL_IFFREE (pre, ILLlp_predata);
		}

		if (info)
		{
			ILLlp_sinfo_free (info);
			ILL_IFFREE (info, ILLlp_sinfo);
		}
	}
	else
	{
		lp->presolve = pre;
		lp->sinfo = info;
	}

	ILL_RETURN (rval, "ILLlp_presolve");
}


#if 0
int ILLlp_presolve_addrow (
	lpinfo * lp,
	int cnt,
	int *ind,
	double *val,
	double rhs)
{
	int rval = 0;
	ILLlpdata *qslp;
	ILLlp_sinfo *S;
	ILLmatrix *A;

	/* This will need to evolve into a function that handles the task */
	/* of working through the presolve data to determine the new LP   */
	/* created when a row is added to the original LP.                */

	/* The copies of the obj and bound used in the simplex code are   */
	/* also updated in this function.                                 */

	if (!lp)
	{
		fprintf (stderr, "ILLlp_presolve_addrow is called without an LP\n");
		rval = 1;
		goto CLEANUP;
	}

	if (lp->presolve != 0)
	{
		fprintf (stderr, "Not yet set up to handle addrows after presolve\n");
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
	ILLlpdata * lp,
	ILLlp_predata * pre,
	ILLlp_sinfo * info,
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
		fprintf (stderr, "simple_presolve called with a NULL pointer\n");
		rval = 1;
		goto CLEANUP;
	}

	printf ("Initial Rows = %d, Cols = %d, Nzcount = %d\n",
					lp->nrows, lp->ncols, lp->nzcount);
	fflush (stdout);

	rval = build_graph (lp, &G);
	ILL_CLEANUP_IF (rval);
	if (debug)
		dump_graph (&G);

	if (pre_types & ILL_PRE_FIXED)
	{
		rval = fixed_variables (&G, pre);
		ILL_CLEANUP_IF (rval);
	}

	do
	{
		hit = 0;
		if (pre_types & ILL_PRE_SINGLE_ROW)
		{
			rval = singleton_rows (&G, pre, &newhit);
			ILL_CLEANUP_IF (rval);
			hit += newhit;
		}

		if (pre_types & ILL_PRE_FORCING)
		{
			rval = forcing_constraints (&G, pre, &newhit);
			ILL_CLEANUP_IF (rval);
			hit += newhit;
		}

		if (pre_types & ILL_PRE_SINGLE_COL)
		{
			rval = singleton_columns (&G, pre, &newhit);
			ILL_CLEANUP_IF (rval);
			hit += newhit;
		}

		if (pre_types & ILL_PRE_DUPLICATE_ROW)
		{
			rval = duplicate_rows (&G, &newhit);
			ILL_CLEANUP_IF (rval);
			hit += newhit;
		}

		if (pre_types & ILL_PRE_DUPLICATE_COL)
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
            printf ("Current NZCOUNT = %d\n", cnt); fflush (stdout);
        }
*/
	} while (hit);

	if (ILL_PRE_EMPTY_COL)
	{
		rval = empty_columns (&G, pre);
		ILL_CLEANUP_IF (rval);
	}

	if (debug)
	{
		printf ("Operations\n");
		for (i = 0; i < pre->opcount; i++)
		{
			switch (pre->oplist[i].ptype)
			{
			case ILL_PRE_DELETE_EMPTY_ROW:
				printf ("Delete Empty Row: %d\n", pre->oplist[i].rowindex);
				fflush (stdout);
				break;
			case ILL_PRE_DELETE_SINGLETON_ROW:
				printf ("Delete Singleton Row: %d (col %d)\n",
								pre->oplist[i].rowindex, pre->oplist[i].colindex);
				fflush (stdout);
				dump_line (&pre->oplist[i].line);
				break;
			case ILL_PRE_DELETE_FIXED_VARIABLE:
				printf ("Delete Fixed Variable: %d\n", pre->oplist[i].colindex);
				fflush (stdout);
				dump_line (&pre->oplist[i].line);
				break;
			case ILL_PRE_DELETE_FORCED_VARIABLE:
				printf ("Delete Forced Variable: %d\n", pre->oplist[i].colindex);
				fflush (stdout);
				dump_line (&pre->oplist[i].line);
				break;
			case ILL_PRE_DELETE_SINGLETON_VARIABLE:
				printf ("Delete Singleton Variable: %d\n", pre->oplist[i].colindex);
				fflush (stdout);
				dump_line (&pre->oplist[i].line);
				break;
			case ILL_PRE_DELETE_FREE_SINGLETON_VARIABLE:
				printf ("Delete Free Singleton Variable: %d\n",
								pre->oplist[i].colindex);
				fflush (stdout);
				dump_line (&pre->oplist[i].line);
				break;
			case ILL_PRE_DELETE_EMPTY_COLUMN:
				printf ("Delete Empty Column: %d\n", pre->oplist[i].colindex);
				fflush (stdout);
				dump_line (&pre->oplist[i].line);
				break;
			default:
				fprintf (stderr, "unknon presolve operation\n");
				rval = 1;
				goto CLEANUP;
			}
		}
		printf ("\n");
	}

	rval = grab_lp_info (&G, lp->colnames, info);
	ILL_CLEANUP_IF (rval);

/*
    printf ("Final Rows = %d, Cols = %d, Nzcount = %d\n",
               info->nrows, info->ncols, info->nzcount);
    fflush (stdout);
*/


CLEANUP:

	free_graph (&G);
	ILL_RETURN (rval, "simple_presolve");
}

static int grab_lp_line (
	graph * G,
	int indx,
	ILLlp_preline * line,
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

		line->val = EGlpNumAllocArray (line->count);
		if (!line->ind || !line->val)
		{
			fprintf (stderr, "out of memory in grab_lp_line\n");
			rval = 1;
			goto CLEANUP;
		}
		for (k = 0, cnt = 0; k < n->deg; k++)
		{
			if (n->adj[k]->del == 0)
			{
				line->ind[cnt] = n->adj[k]->row;
				EGlpNumCopy (line->val[cnt], n->adj[k]->coef);
				cnt++;
			}
		}
	}

	if (row_or_col == 0)
	{
		EGlpNumCopy (line->rhs, n->rhs);
	}
	else
	{
		EGlpNumCopy (line->obj, n->obj);
		EGlpNumCopy (line->lower, n->lower);
		EGlpNumCopy (line->upper, n->upper);
	}

	line->row_or_col = row_or_col;

CLEANUP:

	ILL_RETURN (rval, "grab_lp_line");
}

static void dump_line (
	ILLlp_preline * line)
{
	int k;

	printf (" ");
	if (line->row_or_col == 0)
	{
		for (k = 0; k < line->count; k++)
		{
			printf (" C%d->%g", line->ind[k], EGlpNumToLf (line->val[k]));
		}
		printf (" RHS->%g\n", EGlpNumToLf (line->rhs));
	}
	else
	{
		for (k = 0; k < line->count; k++)
		{
			printf (" R%d->%g", line->ind[k], EGlpNumToLf (line->val[k]));
		}
		printf (" Obj->%g  LB->%g  UB->%g\n", EGlpNumToLf (line->obj),
						EGlpNumToLf (line->lower), EGlpNumToLf (line->upper));
	}
	fflush (stdout);
}

static int grab_lp_info (
	graph * G,
	char **colnames,
	ILLlp_sinfo * info)
{
	int rval = 0;
	int ncols = 0, nrows = 0, nzcount = 0;
	int i, j, k, cnt, len;
	node *grows = G->rows;
	node *gcols = G->cols;
	int *tdeg = 0;
	int *map = 0;
	char *buf = 0;
	ILLmatrix *A = &info->A;

	ILL_SAFE_MALLOC (tdeg, G->ncols, int);
	ILL_SAFE_MALLOC (map, G->nrows, int);

	if (!tdeg || !map)
	{
		fprintf (stderr, "out of memory in grab_lp_info\n");
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

	info->rhs = EGlpNumAllocArray (nrows);
	info->obj = EGlpNumAllocArray (ncols);
	info->upper = EGlpNumAllocArray (ncols);
	info->lower = EGlpNumAllocArray (ncols);
	A->matval = EGlpNumAllocArray (info->nzcount + 1);
	ILL_SAFE_MALLOC (A->matind, info->nzcount + 1, int);
	ILL_SAFE_MALLOC (A->matcnt, info->colsize, int);
	ILL_SAFE_MALLOC (A->matbeg, info->colsize, int);

	if (!info->rhs || !info->obj || !info->lower || !info->upper ||
			!A->matval || !A->matind || !A->matcnt || !A->matbeg)
	{
		fprintf (stderr, "out of memory in grab_lp\n");
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
			EGlpNumCopy (info->rhs[nrows], grows[i].rhs);
			nrows++;
		}
	}

	ncols = 0;
	cnt = 0;
	for (j = 0; j < G->ncols; j++)
	{
		if (gcols[j].del == 0)
		{
			EGlpNumCopy (info->obj[ncols], gcols[j].obj);
			EGlpNumCopy (info->lower[ncols], gcols[j].lower);
			EGlpNumCopy (info->upper[ncols], gcols[j].upper);
			A->matcnt[ncols] = tdeg[ncols];
			A->matbeg[ncols] = cnt;
			for (k = 0; k < gcols[j].deg; k++)
			{
				if (gcols[j].adj[k]->del == 0)
				{
					EGlpNumCopy (A->matval[cnt], gcols[j].adj[k]->coef);
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
			fprintf (stderr, "out of memory in grab_lp\n");
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
			fprintf (stderr, "out of memory in grab_lp\n");
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
						fprintf (stderr, "out of memory in grab_lp\n");
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
						fprintf (stderr, "problem with graph in grab_lp\n");
						rval = 1;
						goto CLEANUP;
					}
					sprintf (buf, "s%d", i);
					len = strlen (buf) + 1;
					ILL_SAFE_MALLOC (info->colnames[ncols], len, char);

					if (!info->colnames[ncols])
					{
						fprintf (stderr, "out of memory in grab_lp\n");
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
		ILLlp_sinfo_free (info);
	}
	ILL_IFFREE (tdeg, int);
	ILL_IFFREE (map, int);
	ILL_IFFREE (buf, char);

	ILL_RETURN (rval, "grab_lp_info");
}

static int fixed_variables (
	graph * G,
	ILLlp_predata * pre)
{
	int rval = 0;
	int j;
	int ncols = G->ncols;
	node *cols = G->cols;
	ILLlp_preop *op = 0;

	for (j = 0; j < ncols; j++)
	{
		if (cols[j].del == 0)
		{
			if (EGlpNumIsEqqual (cols[j].lower, cols[j].upper))
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
	ILLlp_predata * pre)
{
	int rval = 0;
	int j, k;
	int ncols = G->ncols;
	node *cols = G->cols;
	ILLlp_preop *op = 0;
	EGlpNum_t objtmp;

	EGlpNumInitVar (objtmp);

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
				EGlpNumCopy (objtmp, cols[j].obj);
				if (G->objsense < 0)
					EGlpNumSign (objtmp);
				if (!EGlpNumIsNeqZero (objtmp, ILL_PRE_FEAS_TOL))
				{
					set_fixed_variable (G, j, cols[j].lower);
				}
				else if (EGlpNumIsGreatZero (objtmp))
				{
					if (EGlpNumIsEqqual (cols[j].lower, ILL_MINDOUBLE))
					{
						printf ("unbounded prob detected in empty_columns\n");
						printf ("col %d, obj %g\n", j, EGlpNumToLf (cols[j].obj));
						fflush (stdout);
						rval = 1;
						goto CLEANUP;
					}
					else
					{
						set_fixed_variable (G, j, cols[j].lower);
					}
				}
				else if (EGlpNumIsLessZero (objtmp))
				{
					if (EGlpNumIsEqqual (cols[j].upper, ILL_MAXDOUBLE))
					{
						printf ("unbounded prob detected in empty_columns\n");
						printf ("col %d, obj %g\n", j, EGlpNumToLf (cols[j].obj));
						fflush (stdout);
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

	EGlpNumClearVar (objtmp);
	ILL_RETURN (rval, "empty_columns");
}

static int singleton_rows (
	graph * G,
	ILLlp_predata * pre,
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
	EGlpNum_t val;
	ILLlp_preop *op = 0;

	EGlpNumInitVar (val);

	*hit = 0;
	if (G->nrows == 0)
		goto CLEANUP;

	ILL_SAFE_MALLOC (tdeg, G->nrows, int);

	if (!tdeg)
	{
		fprintf (stderr, "out of memory in singleton_rows\n");
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
			if (EGlpNumIsNeqZero (r->rhs, ILL_PRE_FEAS_TOL))
			{
				printf ("infeasible row detected in singleton_row\n");
				printf ("empty row with rhs = %g\n", EGlpNumToLf (r->rhs));
				fflush (stdout);
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
				fprintf (stderr, "lost an edge in singleton_rows\n");
				rval = 1;
				goto CLEANUP;
			}

			pivot = r->adj[k];
			c = &cols[pivot->col];

			/*  Store data from operation (incluing the col coefs) */

			op->ptype = ILL_PRE_DELETE_SINGLETON_ROW;
			op->rowindex = rowindex;
			op->colindex = c - cols;
			EGlpNumCopy (op->line.rhs, r->rhs);
			rval = grab_lp_line (G, op->colindex, &op->line, 1);
			ILL_CLEANUP_IF (rval);

			/*  Fix the x[c] to its rhs value */
			/*val = r->rhs / pivot->coef; */
			EGlpNumCopyFrac (val, r->rhs, pivot->coef);
			/* if (val < c->lower - ILL_PRE_FEAS_TOL ||
			 * val > c->upper + ILL_PRE_FEAS_TOL) */
			if (EGlpNumIsSumLess (val, ILL_PRE_FEAS_TOL, c->lower) ||
					EGlpNumIsSumLess (c->upper, ILL_PRE_FEAS_TOL, val))
			{
				printf ("infeasible bounds detected in singleton_row %d\n", rowindex);
				printf ("lower->%g  upper->%g  val = %g\n",
								EGlpNumToLf (c->lower), EGlpNumToLf (c->upper),
								EGlpNumToLf (val));
				fflush (stdout);
				rval = 1;
				goto CLEANUP;
			}
			else
			{
				EGlpNumCopy (c->lower, val);
				EGlpNumCopy (c->upper, val);
			}

			/*  Delete x[c] from other rows (and adjust their rhs) */

			c->del = 1;

			for (h = 0; h < c->deg; h++)
			{
				f = c->adj[h];
				if (f->del == 0)
				{
					/*rows[f->row].rhs -= (f->coef * c->lower); */
					EGlpNumSubInnProdTo (rows[f->row].rhs, f->coef, c->lower);
					tdeg[f->row]--;
					if (tdeg[f->row] == 1)
					{
						if (f == pivot)
						{
							fprintf (stderr, "bad pivot element\n");
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
	EGlpNumClearVar (val);
	ILL_RETURN (rval, "singleton_rows");
}

static int forcing_constraints (
	graph * G,
	ILLlp_predata * pre,
	int *hit)
{
	int rval = 0;
	int i, j, k, ts;
	node *rows = G->rows;
	node *cols = G->cols;
	edge *e;
	int nrows = G->nrows;
	EGlpNum_t ub, lb;
	ILLlp_preop *op = 0;

	EGlpNumInitVar (ub);
	EGlpNumInitVar (lb);

	*hit = 0;

	for (i = 0; i < nrows; i++)
	{
		if (rows[i].del == 0)
		{
			get_implied_rhs_bounds (G, i, &lb, &ub);
			if (EGlpNumIsSumLess (rows[i].rhs, ILL_PRE_FEAS_TOL, lb) ||
					EGlpNumIsSumLess (ub, ILL_PRE_FEAS_TOL, rows[i].rhs))
			{
				printf ("infeasible row detected in forcing_constraints\n");
				printf ("Row %d:  RHS->%g  LBnd->%g  UBnd->%g\n",
								i, EGlpNumToLf (rows[i].rhs),
								EGlpNumToLf (lb), EGlpNumToLf (ub));
				fflush (stdout);
				rval = 1;
				goto CLEANUP;
			}
			else if (EGlpNumIsDiffLess (rows[i].rhs, ILL_PRE_FEAS_TOL, lb) ||
							 EGlpNumIsDiffLess (ub, ILL_PRE_FEAS_TOL, rows[i].rhs))
			{
				(*hit)++;
				ts = (EGlpNumIsDiffLess (rows[i].rhs, ILL_PRE_FEAS_TOL, lb) ? 0 : 1);
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

						if ((ts == 0 && EGlpNumIsLessZero (e->coef)) ||
								(ts == 1 && EGlpNumIsGreatZero (e->coef)))
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

	EGlpNumClearVar (ub);
	EGlpNumClearVar (lb);
	ILL_RETURN (rval, "forcing_constraints");
}

static int singleton_columns (
	graph * G,
	ILLlp_predata * pre,
	int *hit)
{
	int rval = 0;
	int ncols = G->ncols;
	int j, k, deg, rdeg, single = 0, irow;
	EGlpNum_t lb, ub, b, eb;
	node *cols = G->cols;
	node *rows = G->rows;
	edge *b_edge;
	ILLlp_preop *op = 0;
	EGlpNum_t newub, newlb;
	EGlpNum_t a, c, l, u;

	EGlpNumInitVar (lb);
	EGlpNumInitVar (ub);
	EGlpNumInitVar (eb);
	EGlpNumInitVar (b);
	EGlpNumInitVar (newlb);
	EGlpNumInitVar (newub);
	EGlpNumInitVar (a);
	EGlpNumInitVar (c);
	EGlpNumInitVar (l);
	EGlpNumInitVar (u);

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
				EGlpNumCopy (b, cols[j].adj[single]->coef);
				b_edge = cols[j].adj[single];

				get_implied_variable_bounds (G, j, b_edge, &lb, &ub);

				/*if (lb >= cols[j].lower && ub <= cols[j].upper) */
				if (EGlpNumIsLeq (cols[j].lower, lb) &&
						EGlpNumIsLeq (ub, cols[j].upper))
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
					EGlpNumCopyFrac (eb, cols[j].obj, b);

					for (k = 0; k < rows[irow].deg; k++)
					{
						a_edge = rows[irow].adj[k];
						if (a_edge->del == 0 && a_edge != b_edge)
						{
							/*cols[a_edge->col].obj -= (eb * a_edge->coef); */
							EGlpNumSubInnProdTo (cols[a_edge->col].obj, eb, a_edge->coef);
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

						EGlpNumCopy (newub, ILL_MAXDOUBLE);
						EGlpNumCopy (newlb, ILL_MINDOUBLE);
						EGlpNumZero (a);

						/*    ay + bx = c                                */
						/*    l <= x <= u                                */
						/*      x - is column singleton                  */
						/*      derive bounds on y and substitute out x  */

						EGlpNumCopy (c, rows[irow].rhs);
						EGlpNumCopy (l, cols[j].lower);
						EGlpNumCopy (u, cols[j].upper);

						/* Find the ay term */

						for (k = 0; k < rows[irow].deg; k++)
						{
							if (rows[irow].adj[k]->del == 0 && rows[irow].adj[k]->col != j)
							{
								a_edge = rows[irow].adj[k];
								EGlpNumCopy (a, rows[irow].adj[k]->coef);
								col2 = rows[irow].adj[k]->col;
								break;
							}
						}
						if (k == rows[irow].deg)
						{
							fprintf (stderr, "graph error in singleton_col\n");
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
						EGlpNumCopyFrac (eb, a, b);
						if (EGlpNumIsGreatZero (eb))
						{
							/*if (l > -ILL_MAXDOUBLE) */
							if (EGlpNumIsLess (ILL_MINDOUBLE, l))
							{
								/*newub = (c / a) - (l * b) / a; */
								EGlpNumCopy (newub, c);
								EGlpNumSubInnProdTo (newub, l, b);
								EGlpNumDivTo (newub, a);
							}
							/*if (u < ILL_MAXDOUBLE) */
							if (EGlpNumIsLess (u, ILL_MAXDOUBLE))
							{
								/*newlb = (c / a) - (u * b) / a; */
								EGlpNumCopy (newlb, c);
								EGlpNumSubInnProdTo (newlb, u, b);
								EGlpNumDivTo (newlb, a);
							}
						}
						else
						{
							/*if (l > -ILL_MAXDOUBLE) */
							if (EGlpNumIsLess (ILL_MINDOUBLE, l))
							{
								/*newlb = (c / a) - (l * b) / a; */
								EGlpNumCopy (newlb, c);
								EGlpNumSubInnProdTo (newlb, l, b);
								EGlpNumDivTo (newlb, a);
							}
							/*if (u < ILL_MAXDOUBLE) */
							if (EGlpNumIsLess (u, ILL_MAXDOUBLE))
							{
								/*newub = (c / a) - (u * b) / a; */
								EGlpNumCopy (newub, c);
								EGlpNumSubInnProdTo (newub, u, b);
								EGlpNumDivTo (newub, a);
							}
						}

						if (EGlpNumIsLess (cols[col2].lower, newlb))
							EGlpNumCopy (cols[col2].lower, newlb);
						if (EGlpNumIsLess (newub, cols[col2].upper))
							EGlpNumCopy (cols[col2].upper, newub);
						EGlpNumSubTo (cols[col2].obj, eb);

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

	EGlpNumClearVar (lb);
	EGlpNumClearVar (ub);
	EGlpNumClearVar (eb);
	EGlpNumClearVar (b);
	EGlpNumClearVar (newlb);
	EGlpNumClearVar (newub);
	EGlpNumClearVar (a);
	EGlpNumClearVar (c);
	EGlpNumClearVar (l);
	EGlpNumClearVar (u);
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
	EGlpNum_t *f = 0;
	double szeit = ILLutil_zeit ();
	EGlpNum_t q;
	int i, j, k, k2, ri, r0 = 0, n, nu = 0, got, t0, t = 1;
	node *c;

	EGlpNumInitVar (q);


	/*  Code follows J. Tomlin and J. S. Welch, OR Letters 5 (1986) 7--11 */

	*hit = 0;
	if (nrows == 0)
		goto CLEANUP;

	ILL_SAFE_MALLOC (s, nrows, int);

	f = EGlpNumAllocArray (nrows);

	for (i = 0; i < nrows; i++)
	{
		if (rows[i].del || rows[i].rowsense != 'E')
		{
			s[i] = ILL_MAXINT;				/* ILL_MAXINT means no longer eligible    */
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
				EGlpNumCopy (f[ri], c->adj[k]->coef);
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
						EGlpNumCopy (q, c->adj[k]->coef);
						EGlpNumMultTo (q, f[i]);
						EGlpNumDivTo (q, f[ri]);
						EGlpNumDivTo (q, c->adj[k2]->coef);
						if (EGlpNumIsEqual (q, oneLpNum, ILL_PRE_ZERO_TOL))
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
					s[ri] = ILL_MAXINT;
					if (--nu == 0)
						goto DONE;
				}
			}
		}

		if (n == 1)
		{
			s[r0] = ILL_MAXINT;
			if (--nu == 0)
				goto DONE;
		}
	}

DONE:

	{
		int idup = 0;

		for (i = 0; i < nrows; i++)
		{
			if (s[i] > 0 && s[i] < ILL_MAXINT)
			{
				printf ("Row %d: %d\n", i, s[i]);
				idup++;
			}
		}
		printf ("Number of duplicate rows: %d\n", idup);
	}

	printf ("Time in duplicate_rows: %.2f (seconds)\n", ILLutil_zeit () - szeit);
	fflush (stdout);

CLEANUP:

	ILL_IFFREE (s, int);

	EGlpNumFreeArray (f);
	EGlpNumClearVar (q);
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
	EGlpNum_t *f = 0;
	double szeit = ILLutil_zeit ();
	EGlpNum_t q;
	int i, j, k, k2, ci, c0 = 0, n, nu = 0, got, t0, t = 1;
	node *r;

	EGlpNumInitVar (q);


	/*  Code follows J. Tomlin and J. S. Welch, OR Letters 5 (1986) 7--11 */

	*hit = 0;
	if (ncols == 0)
		goto CLEANUP;

	ILL_SAFE_MALLOC (s, ncols, int);

	f = EGlpNumAllocArray (ncols);

	for (j = 0; j < ncols; j++)
	{
		if (cols[j].del || cols[j].coltype != ILL_PRE_COL_STRUC)
		{
			s[j] = ILL_MAXINT;				/* ILL_MAXINT means no longer eligible    */
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
				EGlpNumCopy (f[ci], r->adj[k]->coef);
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
						EGlpNumCopy (q, r->adj[k]->coef);
						EGlpNumMultTo (q, f[j]);
						EGlpNumDivTo (q, f[ci]);
						EGlpNumDivTo (q, r->adj[k2]->coef);
						if (EGlpNumIsEqual (q, oneLpNum, ILL_PRE_ZERO_TOL))
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
					s[ci] = ILL_MAXINT;
					if (--nu == 0)
						goto DONE;
				}
			}
		}

		if (n == 1)
		{
			s[c0] = ILL_MAXINT;
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

	printf ("Time in duplicate_cols: %.2f (seconds)\n", ILLutil_zeit () - szeit);
	fflush (stdout);

CLEANUP:

	ILL_IFFREE (s, int);

	EGlpNumFreeArray (f);
	EGlpNumClearVar (q);
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
		if (s[i] < ILL_MAXINT && s[i] > smax)
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
		if (s[i] < ILL_MAXINT)
		{
			cnt[s[i]]++;
		}
	}

	if (cnt[0] > 0)
		printf ("%d Empty Lines\n", cnt[0]);

	printf ("Duplicate Classes:");
	fflush (stdout);
	for (i = 1; i < smax + 1; i++)
	{
		if (cnt[i] > 1)
		{
			ndup++;
			printf (" %d", cnt[i]);
		}
	}
	printf ("  Number %d\n", ndup);
	fflush (stdout);

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
		if (s[i] < ILL_MAXINT && s[i] > 0)
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
			printf (" %d", (*dupind)[j]);
		}
		printf (" | ");
		fflush (stdout);
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
	EGlpNum_t val)
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
			EGlpNumSubInnProdTo (G->rows[e->row].rhs, e->coef, val);
			e->del = 1;
		}
	}
}

static void get_implied_rhs_bounds (
	graph * G,
	int i,
	EGlpNum_t * lb,
	EGlpNum_t * ub)
{
	int k;
	EGlpNum_t l, u;
	node *cols = G->cols;
	node *rows = G->rows;
	edge *e;

	EGlpNumInitVar (u);
	EGlpNumInitVar (l);

	EGlpNumZero (l);
	for (k = 0; k < rows[i].deg; k++)
	{
		e = rows[i].adj[k];
		if (e->del == 0)
		{
			if (EGlpNumIsLessZero (e->coef))
			{
				if (EGlpNumIsEqqual (cols[e->col].upper, ILL_MAXDOUBLE))
				{
					EGlpNumCopy (l, ILL_MINDOUBLE);
					break;
				}
				else
				{
					/*l += (e->coef * cols[e->col].upper); */
					EGlpNumAddInnProdTo (l, e->coef, cols[e->col].upper);
				}
			}
			else if (EGlpNumIsGreatZero (e->coef))
			{
				if (EGlpNumIsEqqual (cols[e->col].lower, ILL_MINDOUBLE))
				{
					EGlpNumCopy (l, ILL_MINDOUBLE);
					break;
				}
				else
				{
					/*l += (e->coef * cols[e->col].lower); */
					EGlpNumAddInnProdTo (l, e->coef, cols[e->col].lower);
				}
			}
		}
	}

	EGlpNumZero (u);
	for (k = 0; k < rows[i].deg; k++)
	{
		e = rows[i].adj[k];
		if (e->del == 0)
		{
			if (EGlpNumIsLessZero (e->coef ))
			{
				if (EGlpNumIsEqqual (cols[e->col].lower, ILL_MINDOUBLE))
				{
					EGlpNumCopy (u, ILL_MAXDOUBLE);
				}
				else
				{
					/*u += (e->coef * cols[e->col].lower); */
					EGlpNumAddInnProdTo (u, e->coef, cols[e->col].lower);
				}
			}
			else if (EGlpNumIsGreatZero (e->coef))
			{
				if (EGlpNumIsEqqual (cols[e->col].upper, ILL_MAXDOUBLE))
				{
					EGlpNumCopy (u, ILL_MAXDOUBLE);
				}
				else
				{
					/*u += (e->coef * cols[e->col].upper); */
					EGlpNumAddInnProdTo (u, e->coef, cols[e->col].upper);
				}
			}
		}
	}

	EGlpNumCopy (*lb, l);
	EGlpNumCopy (*ub, u);
	EGlpNumClearVar (u);
	EGlpNumClearVar (l);
}

static void get_implied_variable_bounds (
	graph * G,
	int j,
	edge * a_ij,
	EGlpNum_t * lb,
	EGlpNum_t * ub)
{
	int i = a_ij->row;
	EGlpNum_t l, u;

	EGlpNumInitVar (u);
	EGlpNumInitVar (l);

	get_implied_rhs_bounds (G, i, &l, &u);
	EGlpNumCopy (*lb, ILL_MINDOUBLE);
	EGlpNumCopy (*ub, ILL_MAXDOUBLE);

	if (EGlpNumIsLess (ILL_PRE_FEAS_TOL, a_ij->coef))
	{
		if (EGlpNumIsLess (u, ILL_MAXDOUBLE))
		{
			/**lb = (G->rows[i].rhs - u) / a_ij->coef + G->cols[j].upper;*/
			EGlpNumCopyDiffRatio (*lb, G->rows[i].rhs, u, a_ij->coef);
			EGlpNumAddTo (*lb, G->cols[j].upper);
		}
		if (EGlpNumIsLess (ILL_MINDOUBLE, l))
		{
			/**ub = (G->rows[i].rhs - l) / a_ij->coef + G->cols[j].lower;*/
			EGlpNumCopyDiffRatio (*ub, G->rows[i].rhs, l, a_ij->coef);
			EGlpNumAddTo (*ub, G->cols[j].lower);
		}
	}
	else if (EGlpNumIsLess (a_ij->coef, ILL_PRE_FEAS_TOL))
	{
		if (EGlpNumIsLess (ILL_MINDOUBLE, l))
		{
			/**lb = (G->rows[i].rhs - l) / a_ij->coef + G->cols[j].upper;*/
			EGlpNumCopyDiffRatio (*lb, G->rows[i].rhs, l, a_ij->coef);
			EGlpNumAddTo (*lb, G->cols[j].upper);
		}
		if (EGlpNumIsLess (u, ILL_MAXDOUBLE))
		{
			/**ub = (G->rows[i].rhs - u) / a_ij->coef + G->cols[j].lower;*/
			EGlpNumCopyDiffRatio (*ub, G->rows[i].rhs, u, a_ij->coef);
			EGlpNumAddTo (*ub, G->cols[j].lower);
		}
	}
	EGlpNumClearVar (u);
	EGlpNumClearVar (l);
}

static int get_next_preop (
	ILLlp_predata * pre,
	ILLlp_preop ** op)
{
	int rval = 0;

	if (pre->opcount >= pre->opsize)
	{
		pre->opsize *= 1.3;
		pre->opsize += 1000;
		if (pre->opsize < pre->opcount + 1)
			pre->opsize = pre->opcount + 1;
		pre->oplist = EGrealloc (pre->oplist, sizeof (ILLlp_preop) * pre->opsize);
		//rval = ILLutil_reallocrus_scale ((void **) &pre->oplist,
		//                                 &pre->opsize, pre->opcount + 1, 1.3,
		//                                 sizeof (ILLlp_preop));
		//ILL_CLEANUP_IF (rval);
	}
	*op = &pre->oplist[pre->opcount];
	ILLlp_preop_init (*op);

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
	ILLlpdata * lp,
	graph * G)
{
	int rval = 0;
	int ncols = lp->ncols;
	int nrows = lp->nrows;
	int nzcount = lp->nzcount;
	int i, j, k, stop, count;
	edge *edgelist;
	node *rows, *cols;
	ILLmatrix *A = &lp->A;

	G->objsense = lp->objsense;

	ILL_SAFE_MALLOC (G->rows, nrows, node);
	if (!G->rows)
	{
		fprintf (stderr, "out of memory in build_graph\n");
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
		EGlpNumInitVar ((G->edgelist[i].coef));
	G->nzcount = nzcount;
	ILL_SAFE_MALLOC (G->adjspace, 2 * nzcount, edge *);

	if (!G->cols || !G->edgelist || !G->adjspace)
	{
		fprintf (stderr, "out of memory in build_graph\n");
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
		EGlpNumCopy (cols[j].obj, lp->obj[j]);
		EGlpNumCopy (cols[j].lower, lp->lower[j]);
		EGlpNumCopy (cols[j].upper, lp->upper[j]);
		stop = A->matbeg[j] + A->matcnt[j];
		for (k = A->matbeg[j]; k < stop; k++)
		{
			i = A->matind[k];
			rows[i].adj[rows[i].deg++] = &(edgelist[count]);
			cols[j].adj[cols[j].deg++] = &(edgelist[count]);
			edgelist[count].row = i;
			edgelist[count].col = j;
			EGlpNumCopy (edgelist[count].coef, A->matval[k]);
			edgelist[count].mark = 0;
			edgelist[count].del = 0;
			edgelist[count].coltype = cols[j].coltype;
			count++;
		}
	}
	if (count != nzcount)
	{
		fprintf (stderr, "counts are off in build_graph\n");
		rval = 1;
		goto CLEANUP;
	}

	G->ecount = count;
	G->nrows = nrows;
	G->ncols = ncols;

	for (i = 0; i < G->nrows; i++)
	{
		G->rows[i].del = 0;
		EGlpNumCopy (G->rows[i].rhs, lp->rhs[i]);
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

	printf ("ecount = %d, nrows = %d, ncols = %d\n",
					G->ecount, G->nrows, G->ncols);
	fflush (stdout);

	for (i = 0; i < G->nrows; i++)
	{
		printf ("Row %d:", i);
		for (k = 0; k < G->rows[i].deg; k++)
		{
			printf (" %d", G->rows[i].adj[k]->col);
			if (G->rows[i].adj[k]->coltype == ILL_PRE_COL_LOGICAL)
				printf ("S");
			printf ("(%g)", EGlpNumToLf (G->rows[i].adj[k]->coef));
		}
		printf ("  rhs: %g", EGlpNumToLf (G->rows[i].rhs));
		if (G->rows[i].del)
		{
			printf (" (deleted)\n");
		}
		else
		{
			printf ("\n");
		}
	}

	for (j = 0; j < G->ncols; j++)
	{
		if (G->cols[j].coltype == ILL_PRE_COL_LOGICAL)
		{
			printf ("Slk %d:", j);
		}
		else
		{
			printf ("Col %d:", j);
		}
		for (k = 0; k < G->cols[j].deg; k++)
		{
			printf (" %d", G->cols[j].adj[k]->row);
		}
		printf ("  obj: %g  bnd: (%g, %g)", EGlpNumToLf (G->cols[j].obj),
						EGlpNumToLf (G->cols[j].lower), EGlpNumToLf (G->cols[j].upper));
		if (G->cols[j].del)
		{
			printf (" (deleted)\n");
		}
		else
		{
			printf ("\n");
		}
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
			EGlpNumClearVar ((G->edgelist[i].coef));
		ILL_IFFREE (G->edgelist, edge);
		ILL_IFFREE (G->rows, node);
		ILL_IFFREE (G->cols, node);
		ILL_IFFREE (G->adjspace, edge *);
		if (intptr_check_leaks (&G->intptrworld, &total, &onlist))
		{
			fprintf (stderr, "WARNING: %d outstanding intptrs\n", total - onlist);
		}
		ILLptrworld_delete (&G->intptrworld);
		init_graph (G);
	}
}

int ILLlp_sinfo_print (
	ILLlp_sinfo * s)
{
	int rval = 0;
	int i;
	ILLlpdata lp;
	char *sense = 0;

	ILLlpdata_init (&lp);

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
		fprintf (stderr, "out of memory in ILLlp_sinfo_print\n");
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

	ILL_RETURN (rval, "ILLlp_sinfo_print");
}

void ILLlp_sinfo_init (
	ILLlp_sinfo * sinfo)
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
		sinfo->objsense = ILL_MIN;
		ILLmatrix_init (&sinfo->A);
	}
}

void ILLlp_sinfo_free (
	ILLlp_sinfo * sinfo)
{
	if (sinfo)
	{
		EGlpNumFreeArray (sinfo->obj);
		EGlpNumFreeArray (sinfo->lower);
		EGlpNumFreeArray (sinfo->upper);
		EGlpNumFreeArray (sinfo->rhs);
		ILLmatrix_free (&sinfo->A);
		if (sinfo->colnames)
		{
			int i;

			for (i = 0; i < sinfo->ncols; i++)
			{
				ILL_IFFREE (sinfo->colnames[i], char);
			}
			ILL_IFFREE (sinfo->colnames, char *);
		}
		ILLlp_sinfo_init (sinfo);
	}
}

void ILLlp_predata_init (
	ILLlp_predata * pre)
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

void ILLlp_predata_free (
	ILLlp_predata * pre)
{
	if (pre)
	{
		int i;

		for (i = 0; i < pre->opcount; i++)
		{
			ILLlp_preop_free (&pre->oplist[i]);
		}
		ILL_IFFREE (pre->oplist, ILLlp_preop);
		ILL_IFFREE (pre->colmap, int);
		ILL_IFFREE (pre->rowmap, int);

		ILL_IFFREE (pre->colscale, EGlpNum_t);
		ILL_IFFREE (pre->rowscale, EGlpNum_t);
		ILL_IFFREE (pre->colfixval, EGlpNum_t);
		ILL_IFFREE (pre->rowfixval, EGlpNum_t);
		ILLlp_predata_init (pre);
	}
}

void ILLlp_preop_init (
	ILLlp_preop * op)
{
	if (op)
	{
		op->ptype = 0;
		op->rowindex = -1;
		op->colindex = -1;
		ILLlp_preline_init (&op->line);
	}
}

void ILLlp_preop_free (
	ILLlp_preop * op)
{
	if (op)
	{
		ILLlp_preline_free (&op->line);
		ILLlp_preop_init (op);
	}
}

void ILLlp_preline_init (
	ILLlp_preline * line)
{
	if (line)
	{
		EGlpNumInitVar (line->rhs);
		EGlpNumInitVar (line->obj);
		EGlpNumInitVar (line->upper);
		EGlpNumInitVar (line->lower);
		EGlpNumZero (line->rhs);
		EGlpNumZero (line->obj);
		EGlpNumZero (line->upper);
		EGlpNumZero (line->lower);
		line->count = 0;
		line->ind = 0;
		line->val = 0;
	}
}

void ILLlp_preline_free (
	ILLlp_preline * line)
{
	if (line)
	{
		EGlpNumClearVar (line->rhs);
		EGlpNumClearVar (line->obj);
		EGlpNumClearVar (line->upper);
		EGlpNumClearVar (line->lower);
		ILL_IFFREE (line->ind, int);

		EGlpNumFreeArray (line->val);
		//ILLlp_preline_init (line);
	}
}
