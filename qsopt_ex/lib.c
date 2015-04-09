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

/* RCS_INFO = "$RCSfile: lib.c,v $ $Revision: 1.2 $ $Date: 2003/11/05 16:49:52 $"; */
//static int TRACE = 0;

/****************************************************************************/
/*                                                                          */
/*               Interface Routines to Core LP Solver                       */
/*                                                                          */
/*  EXPORTED FUNCTIONS                                                      */
/*                                                                          */
/*    int EGLPNUM_TYPENAME_ILLlib_optimize (EGLPNUM_TYPENAME_lpinfo *lp, EGLPNUM_TYPENAME_ILLlp_basis *B, EGLPNUM_TYPENAME_price_info *pinf,    */
/*            int algo, int *status, int simplex_display)                   */
/*    int EGLPNUM_TYPENAME_ILLlib_cache_solution (EGLPNUM_TYPENAME_lpinfo *lp, EGLPNUM_TYPENAME_ILLlp_cache *C)                */
/*    int EGLPNUM_TYPENAME_ILLlib_solution (EGLPNUM_TYPENAME_lpinfo *lp, EGLPNUM_TYPENAME_ILLlp_cache *C, double *val,         */
/*            double *x, double *pi, double *slack, double *rc)             */
/*    int EGLPNUM_TYPENAME_ILLlib_get_x (EGLPNUM_TYPENAME_lpinfo *lp, EGLPNUM_TYPENAME_ILLlp_cache *C, double *x)              */
/*    int EGLPNUM_TYPENAME_ILLlib_get_slack (EGLPNUM_TYPENAME_lpinfo *lp, EGLPNUM_TYPENAME_ILLlp_cache *C, double *slack)      */
/*    int EGLPNUM_TYPENAME_ILLlib_objval (EGLPNUM_TYPENAME_lpinfo *lp, EGLPNUM_TYPENAME_ILLlp_cache *C, double *val)           */
/*    int EGLPNUM_TYPENAME_ILLlib_newrow (EGLPNUM_TYPENAME_lpinfo *lp, EGLPNUM_TYPENAME_ILLlp_basis *B, double rhs,            */
/*            char sense, double range, const char *name)                   */
/*        -range can specify a rangeval for the row (if sense is not 'R',   */
/*         then range is ignored); it should be 0 if no range is needed;    */
/*         if sense is 'R' but no rangeval array exists for the LP, the     */
/*         array will be allocated and initialized.                         */
/*    int EGLPNUM_TYPENAME_ILLlib_newrows (EGLPNUM_TYPENAME_lpinfo *lp, EGLPNUM_TYPENAME_ILLlp_basis *B, int num, double *rhs, */
/*            char *sense, double *range, const char **names)               */
/*        -range is an array specifying the rangevals for the rows; range   */
/*         should be NULL if no rangevals are needed.                       */
/*    int EGLPNUM_TYPENAME_ILLlib_addrow (EGLPNUM_TYPENAME_lpinfo *lp, EGLPNUM_TYPENAME_ILLlp_basis *B, int cnt, int *ind,     */
/*            double *val, double rhs, char sense, double range,            */
/*            const char *name)                                             */
/*    int EGLPNUM_TYPENAME_ILLlib_addrows (EGLPNUM_TYPENAME_lpinfo *lp, EGLPNUM_TYPENAME_ILLlp_basis *B, int num,              */
/*            int *rmatcnt, int *rmatbeg, int *rmatind, double *rmatval,    */
/*            double *rhs, char *sense, double *range, const char **names,  */
/*            int *factorok)                                                */
/*    int EGLPNUM_TYPENAME_ILLlib_delrows (EGLPNUM_TYPENAME_lpinfo *lp, EGLPNUM_TYPENAME_ILLlp_basis *B,                       */
/*            int num, int *dellist, int *basis_ok)                         */
/*    int EGLPNUM_TYPENAME_ILLlib_newcol (EGLPNUM_TYPENAME_lpinfo *lp, EGLPNUM_TYPENAME_ILLlp_basis *B,                        */
/*            double obj, double lower, double upper, const char *name,     */
/*            int factorok)                                                 */
/*    int EGLPNUM_TYPENAME_ILLlib_newcols (EGLPNUM_TYPENAME_lpinfo *lp, EGLPNUM_TYPENAME_ILLlp_basis *B,                       */
/*            int num, double *obj, double *lower, double *upper,           */
/*            const char **names, int factorok)                             */
/*    int EGLPNUM_TYPENAME_ILLlib_addcol (EGLPNUM_TYPENAME_lpinfo *lp, EGLPNUM_TYPENAME_ILLlp_basis *B,                        */
/*            int cnt, int *ind, double *val, double obj, double lower,     */
/*            double upper, const char *name, int factorok)                 */
/*    int EGLPNUM_TYPENAME_ILLlib_addcols (EGLPNUM_TYPENAME_lpinfo *lp, EGLPNUM_TYPENAME_ILLlp_basis *B,                       */
/*            int num, int *cmatcnt, int *cmatbeg, int *cmatind,            */
/*            double *cmatval, double *obj, double *lower, double *upper,   */
/*            const char **names, int factorok)                             */
/*    int EGLPNUM_TYPENAME_ILLlib_delcols (EGLPNUM_TYPENAME_lpinfo *lp, EGLPNUM_TYPENAME_ILLlp_basis *B, int num, int *dellist */
/*            int *basis_ok)                                                */
/*    int EGLPNUM_TYPENAME_ILLlib_chgcoef (EGLPNUM_TYPENAME_lpinfo *lp, int rowindex, int colindex,           */
/*            double coef)                                                  */
/*    int EGLPNUM_TYPENAME_ILLlib_chgsense (EGLPNUM_TYPENAME_lpinfo *lp, int num, int *rowlist, char *sense)  */
/*    int EGLPNUM_TYPENAME_ILLlib_getrows (EGLPNUM_TYPENAME_lpinfo *lp, int num, int *rowlist, int **rowcnt,  */
/*            int **rowbeg, int **rowind, double **rowval, double **rhs,    */
/*            char **sense, char ***names)                                  */
/*    int EGLPNUM_TYPENAME_ILLlib_getcols (EGLPNUM_TYPENAME_lpinfo *lp, int num, int *collist, int **colcnt,  */
/*            int **colbeg, int **colind, double **colval, double **obj,    */
/*            double **lower, double **upper, char ***names)                */
/*    int EGLPNUM_TYPENAME_ILLlib_getobj (EGLPNUM_TYPENAME_lpinfo *lp, double *obj)                           */
/*    int EGLPNUM_TYPENAME_ILLlib_chgobj (EGLPNUM_TYPENAME_lpinfo *lp, int indx, double coef)                 */
/*    int EGLPNUM_TYPENAME_ILLlib_getrhs (EGLPNUM_TYPENAME_lpinfo *lp, double *rhs)                           */
/*    int EGLPNUM_TYPENAME_ILLlib_chgrhs (EGLPNUM_TYPENAME_lpinfo *lp, int indx, double coef)                 */
/*    int EGLPNUM_TYPENAME_ILLlib_getintflags (EGLPNUM_TYPENAME_lpinfo *lp, int *intflags)                    */
/*    int EGLPNUM_TYPENAME_ILLlib_rownames (EGLPNUM_TYPENAME_lpinfo *lp, char **rownames)                     */
/*    int EGLPNUM_TYPENAME_ILLlib_colnames (EGLPNUM_TYPENAME_lpinfo *lp, char **colnames)                     */
/*    int EGLPNUM_TYPENAME_ILLlib_colindex (EGLPNUM_TYPENAME_lpinfo *lp, char *name, int *colindex)           */
/*    int EGLPNUM_TYPENAME_ILLlib_rowindex (EGLPNUM_TYPENAME_lpinfo *lp, char *name, int *rowindex)           */
/*    int EGLPNUM_TYPENAME_ILLlib_chgbnd  (EGLPNUM_TYPENAME_lpinfo *lp, int indx, char lu, double bnd)        */
/*    int EGLPNUM_TYPENAME_ILLlib_chgbnds (EGLPNUM_TYPENAME_lpinfo *lp, int cnt, int *indx, char *lu,         */
/*            double *bnd)                                                  */
/*    int EGLPNUM_TYPENAME_ILLlib_getbnd (EGLPNUM_TYPENAME_lpinfo *lp, int indx, char lu, double *bnd)        */
/*    int EGLPNUM_TYPENAME_ILLlib_getbnds (EGLPNUM_TYPENAME_lpinfo *lp, double *lower, double *upper)         */
/*    int EGLPNUM_TYPENAME_ILLlib_strongbranch (EGLPNUM_TYPENAME_lpinfo *lp, EGLPNUM_TYPENAME_price_info *pinf,                */
/*      int *candidatelist, int ncand, double *xlist, double *downpen,      */
/*      double *uppen, int iterations, double objbound)                     */
/*    int EGLPNUM_TYPENAME_ILLlib_getbasis (EGLPNUM_TYPENAME_lpinfo *lp, char *cstat, char *rstat)            */
/*    int EGLPNUM_TYPENAME_ILLlib_loadbasis (EGLPNUM_TYPENAME_ILLlp_basis *B, int nstruct, int nrows,         */
/*      char *cstat, char *rstat)                                           */
/*    int EGLPNUM_TYPENAME_ILLlib_readbasis (EGLPNUM_TYPENAME_lpinfo *lp, EGLPNUM_TYPENAME_ILLlp_basis *B, char *fname)        */
/*    int EGLPNUM_TYPENAME_ILLlib_writebasis (EGLPNUM_TYPENAME_lpinfo *lp, const char *fname)                 */
/*    int EGLPNUM_TYPENAME_ILLlib_getrownorms (EGLPNUM_TYPENAME_lpinfo *lp, EGLPNUM_TYPENAME_price_info *pinf,                 */
/*            double *rownorms)                                             */
/*    int EGLPNUM_TYPENAME_ILLlib_loadrownorms (EGLPNUM_TYPENAME_lpinfo *lp, EGLPNUM_TYPENAME_price_info *pinf,                */
/*            double *rownorms)                                             */
/*    int EGLPNUM_TYPENAME_ILLlib_recompute_rownorms (EGLPNUM_TYPENAME_lpinfo *lp, EGLPNUM_TYPENAME_price_info *pinf)          */
/*    int EGLPNUM_TYPENAME_ILLlib_print_x (EGioFile_t *fd, EGLPNUM_TYPENAME_lpinfo *lp, EGLPNUM_TYPENAME_ILLlp_cache *C, double *x,  */
/*            int nonZerosOnly)                                             */
/*    int EGLPNUM_TYPENAME_ILLlib_print_x (EGLPNUM_TYPENAME_lpinfo *lp, EGLPNUM_TYPENAME_ILLlp_cache *C)                       */
/*    int EGLPNUM_TYPENAME_ILLlib_iter (EGLPNUM_TYPENAME_lpinfo *lp)                                          */
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

#include "qs_config.h"
#include "logging-private.h"

#include "eg_lpnum.h"
#include "eg_io.h"
#include "except.h"
#include "util.h"

#include "lpdata_EGLPNUM_TYPENAME.h"
#include "lpdefs_EGLPNUM_TYPENAME.h"
#include "simplex_EGLPNUM_TYPENAME.h"
#include "price_EGLPNUM_TYPENAME.h"
#include "basis_EGLPNUM_TYPENAME.h"
#include "lib_EGLPNUM_TYPENAME.h"
#include "qstruct_EGLPNUM_TYPENAME.h"
#include "qsopt_EGLPNUM_TYPENAME.h"
#include "lp_EGLPNUM_TYPENAME.h"
#include "mps_EGLPNUM_TYPENAME.h"


static void check_pinf (
	EGLPNUM_TYPENAME_price_info * pinf,
	int *it_exists);

static int matrix_addrow (
	EGLPNUM_TYPENAME_ILLmatrix * A,
	int rowcnt,
	int *rowind,
	const EGLPNUM_TYPE * rowval),
  matrix_addrow_end (
	EGLPNUM_TYPENAME_ILLmatrix * A,
	int row,
	int rowcnt,
	int *rowind,
	const EGLPNUM_TYPE * rowval),
  matrix_addcoef (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_ILLmatrix * A,
	int row,
	int col,
	EGLPNUM_TYPE val),
  matrix_addcol (
	EGLPNUM_TYPENAME_ILLmatrix * A,
	int colcnt,
	int *colind,
	EGLPNUM_TYPE * colval),
  delcols_work (
	EGLPNUM_TYPENAME_lpinfo * lp,
	char *colmark),
  reset_colindex (
	EGLPNUM_TYPENAME_lpinfo * lp),
  reset_rowindex (
	EGLPNUM_TYPENAME_lpinfo * lp);

int EGLPNUM_TYPENAME_ILLlib_optimize (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_ILLlp_basis * B,
	EGLPNUM_TYPENAME_price_info * pinf,
	int algo,
	int *status,
	int simplex_display,
	itcnt_t*itcnt)
{
	int rval = 0;
	int sol_status;

	if (status)
		*status = QS_LP_UNSOLVED;

	/* EGLPNUM_TYPENAME_ILLprice_free_pricing_info (pinf); *//* Should be removed later */

	rval = EGLPNUM_TYPENAME_ILLsimplex (lp, algo, B, pinf, &sol_status, simplex_display, itcnt);
	CHECKRVALG (rval, CLEANUP);

	if (status)
		*status = sol_status;

CLEANUP:

	if (rval == E_SIMPLEX_ERROR)
	{
		EGioFile_t *eout = 0;
		int tval;

		QSlog("write bad lp to error.lp");
#ifdef HAVE_LIBZ
		eout = EGioOpen ("error.lp.gz", "w");
#else
#ifdef HAVE_LIBBZ2
		eout = EGioOpen ("error.lp.bz2", "w");
#else
		eout = EGioOpen ("error.lp", "w");
#endif
#endif
		if (!eout)
		{
			QSlog("could not open file to write bad lp");
		}
		else
		{
			tval = EGLPNUM_TYPENAME_ILLwrite_lp (lp->O, NULL);
			if (tval)
			{
				QSlog("error while writing bad lp");
			}
			EGioClose (eout);
		}

		QSlog("write bad basis to error.bas");
		tval = EGLPNUM_TYPENAME_ILLlib_writebasis (lp, 0, "error.bas");
		if (tval)
		{
			QSlog("error while writing bad basis");
		}
	}
	if (rval == QS_LP_CHANGE_PREC)
	{
		MESSAGE (__QS_SB_VERB, "Changing precision");
		return rval;
	}
	MESSAGE (rval ? 0 : 1000, "Error code %d", rval);
	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLlib_cache_solution (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_ILLlp_cache * C)
{
	int rval = 0;

	if (C)
	{
		if (C->nstruct != lp->O->nstruct || C->nrows != lp->O->nrows)
		{
			QSlog("lp_cache does not match size of lp");
			rval = 1;
			ILL_CLEANUP;
		}
		rval = EGLPNUM_TYPENAME_ILLlib_solution (lp, 0, &(C->val), C->x, C->pi, C->slack, C->rc);
		CHECKRVALG (rval, CLEANUP);
	}

CLEANUP:

	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLlib_solution (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_ILLlp_cache * C,
	EGLPNUM_TYPE * val,
	EGLPNUM_TYPE * x,
	EGLPNUM_TYPE * pi,
	EGLPNUM_TYPE * slack,
	EGLPNUM_TYPE * rc)
{
	int i, rval = 0;
	EGLPNUM_TYPE *tempx = 0;
	EGLPNUM_TYPE *temprc = 0;
	int ncols = lp->O->ncols;
	int nrows = lp->O->nrows;
	int nstruct = lp->O->nstruct;
	EGLPNUM_TYPENAME_ILLlpdata *qslp = lp->O;

	if (C)
	{
		if (C->nrows != nrows || C->nstruct != nstruct)
		{
			QSlog("cache mismatch in EGLPNUM_TYPENAME_ILLlib_solution");
			rval = 0;
			ILL_CLEANUP;
		}
		if (val)
		{
			EGLPNUM_TYPENAME_EGlpNumCopy (*val, C->val);
		}
		if (x)
		{
			for (i = 0; i < nstruct; i++)
			{
				EGLPNUM_TYPENAME_EGlpNumCopy (x[i], C->x[i]);
			}
		}
		if (pi)
		{
			for (i = 0; i < nrows; i++)
			{
				EGLPNUM_TYPENAME_EGlpNumCopy (pi[i], C->pi[i]);
			}
		}
		if (slack)
		{
			for (i = 0; i < nrows; i++)
			{
				EGLPNUM_TYPENAME_EGlpNumCopy (slack[i], C->slack[i]);
			}
		}
		if (rc)
		{
			for (i = 0; i < nstruct; i++)
			{
				EGLPNUM_TYPENAME_EGlpNumCopy (rc[i], C->rc[i]);
			}
		}
	}
	else
	{
		if (x || slack)
			tempx = EGLPNUM_TYPENAME_EGlpNumAllocArray (ncols);

		if (rc)
			temprc = EGLPNUM_TYPENAME_EGlpNumAllocArray (ncols);

		rval = EGLPNUM_TYPENAME_ILLsimplex_solution (lp, tempx, pi, temprc, val);
		CHECKRVALG (rval, CLEANUP);

		if (x)
		{
			for (i = 0; i < nstruct; i++)
			{
				EGLPNUM_TYPENAME_EGlpNumCopy (x[i], tempx[qslp->structmap[i]]);
			}
		}
		if (slack)
		{
			for (i = 0; i < nrows; i++)
			{
				EGLPNUM_TYPENAME_EGlpNumCopy (slack[i], tempx[qslp->rowmap[i]]);
			}
		}

		if (rc)
		{
			for (i = 0; i < nstruct; i++)
			{
				EGLPNUM_TYPENAME_EGlpNumCopy (rc[i], temprc[qslp->structmap[i]]);
			}
		}


		if (lp->O->objsense == EGLPNUM_TYPENAME_ILL_MAX)
		{														/* Reverse signs for max prob */
			if (val)
			{
				EGLPNUM_TYPENAME_EGlpNumSign (*val);
			}
			if (pi)
			{
				for (i = 0; i < nrows; i++)
				{
					EGLPNUM_TYPENAME_EGlpNumSign (pi[i]);
				}
			}
			if (rc)
			{
				for (i = 0; i < nstruct; i++)
				{
					EGLPNUM_TYPENAME_EGlpNumSign (rc[i]);
				}
			}
		}
	}

CLEANUP:

	EGLPNUM_TYPENAME_EGlpNumFreeArray (tempx);
	EGLPNUM_TYPENAME_EGlpNumFreeArray (temprc);
	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLlib_get_x (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_ILLlp_cache * C,
	EGLPNUM_TYPE * x)
{
	int rval = 0;

	rval = EGLPNUM_TYPENAME_ILLlib_solution (lp, C, 0, x, 0, 0, 0);
	CHECKRVALG (rval, CLEANUP);

CLEANUP:

	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLlib_get_slack (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_ILLlp_cache * C,
	EGLPNUM_TYPE * slack)
{
	int rval = 0;

	rval = EGLPNUM_TYPENAME_ILLlib_solution (lp, C, 0, 0, 0, slack, 0);
	CHECKRVALG (rval, CLEANUP);

CLEANUP:

	EG_RETURN (rval);
}


int EGLPNUM_TYPENAME_ILLlib_objval (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_ILLlp_cache * C,
	EGLPNUM_TYPE * val)
{
	int rval = 0;

	if (lp->basisstat.optimal)
	{
		rval = EGLPNUM_TYPENAME_ILLlib_solution (lp, C, val, 0, 0, 0, 0);
		CHECKRVALG (rval, CLEANUP);
	}
	else
	{
		EGLPNUM_TYPENAME_EGlpNumCopy (*val, lp->dobjval);	/* Ask Sanjeeb */
	}

CLEANUP:

	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLlib_tableau (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int row,
	EGLPNUM_TYPE * binv,
	EGLPNUM_TYPE * tabrow)
{
	int rval = 0;
	int i;
	int ncols = lp->O->ncols;
	int nrows = lp->O->nrows;
	int nstruct = lp->O->nstruct;
	EGLPNUM_TYPE *brow = 0;
	EGLPNUM_TYPE *trow = 0;
	EGLPNUM_TYPENAME_ILLlpdata *qslp = lp->O;

	if (row < 0 || row >= qslp->nrows)
	{
		QSlog("EGLPNUM_TYPENAME_ILLlib_tableau called with bad row: %d", row);
		rval = 1;
		ILL_CLEANUP;
	}
	brow = EGLPNUM_TYPENAME_EGlpNumAllocArray (nrows);

	if (tabrow)
		trow = EGLPNUM_TYPENAME_EGlpNumAllocArray (ncols);

	rval = EGLPNUM_TYPENAME_ILLbasis_tableau_row (lp, row, brow, trow, 0, 0);
	CHECKRVALG (rval, CLEANUP);

	if (binv)
	{
		for (i = 0; i < nrows; i++)
		{
			EGLPNUM_TYPENAME_EGlpNumCopy (binv[i], brow[i]);
		}
	}

	if (tabrow)
	{
		for (i = 0; i < nstruct; i++)
		{
			EGLPNUM_TYPENAME_EGlpNumCopy (tabrow[i], trow[qslp->structmap[i]]);
		}
		for (i = 0; i < nrows; i++)
		{
			EGLPNUM_TYPENAME_EGlpNumCopy (tabrow[nstruct + i], trow[qslp->rowmap[i]]);
		}
	}

CLEANUP:

	EGLPNUM_TYPENAME_EGlpNumFreeArray (brow);
	EGLPNUM_TYPENAME_EGlpNumFreeArray (trow);
	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLlib_basis_order (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int *header)
{
	int rval = 0;
	int i, j;
	int ncols = lp->O->ncols;
	int nrows = lp->O->nrows;
	int nstruct = lp->O->nstruct;
	EGLPNUM_TYPENAME_ILLlpdata *qslp = lp->O;
	int *invmap = 0;

	ILL_SAFE_MALLOC (invmap, ncols, int);

	for (j = 0; j < nstruct; j++)
	{
		invmap[qslp->structmap[j]] = j;
	}
	for (i = 0; i < nrows; i++)
	{
		invmap[qslp->rowmap[i]] = nstruct + i;
	}

	for (i = 0; i < nrows; i++)
	{
		header[i] = invmap[lp->baz[i]];
	}

CLEANUP:

	ILL_IFFREE (invmap, int);

	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLlib_chgbnd (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int indx,
	int lu,
	const EGLPNUM_TYPE bnd)
{
	int rval = 0;
	int col;

	if (!lp)
	{
		QSlog("EGLPNUM_TYPENAME_ILLlib_chgbnd called without an lp");
		rval = 1;
		ILL_CLEANUP;
	}

	if (indx < 0 || indx > lp->O->nstruct)
	{
		QSlog("EGLPNUM_TYPENAME_ILLlib_chgbnd called with bad indx: %d", indx);
		rval = 1;
		ILL_CLEANUP;
	}

	if (lp->O->sinfo)
	{															/* Presolve LP is no longer valid, free the data */
		EGLPNUM_TYPENAME_ILLlp_sinfo_free (lp->O->sinfo);
		ILL_IFFREE (lp->O->sinfo, EGLPNUM_TYPENAME_ILLlp_sinfo);
	}

	col = lp->O->structmap[indx];

	switch (lu)
	{
	case 'L':
		EGLPNUM_TYPENAME_EGlpNumCopy (lp->O->lower[col], bnd);
		break;
	case 'U':
		EGLPNUM_TYPENAME_EGlpNumCopy (lp->O->upper[col], bnd);
		break;
	case 'B':
		EGLPNUM_TYPENAME_EGlpNumCopy (lp->O->lower[col], bnd);
		EGLPNUM_TYPENAME_EGlpNumCopy (lp->O->upper[col], bnd);
		break;
	default:
		QSlog("EGLPNUM_TYPENAME_ILLlib_chgbnd called with lu: %c", lu);
		rval = 1;
		ILL_CLEANUP;
	}

CLEANUP:

	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLlib_chgbnds (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int cnt,
	int *indx,
	char *lu,
	const EGLPNUM_TYPE * bnd)
{
	int rval = 0;
	int i;

	for (i = 0; i < cnt; i++)
	{
		rval = EGLPNUM_TYPENAME_ILLlib_chgbnd (lp, indx[i], lu[i], bnd[i]);
		if (rval)
			ILL_CLEANUP;
	}

CLEANUP:

	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLlib_getbnd (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int indx,
	int lu,
	EGLPNUM_TYPE * bnd)
{
	int rval = 0;
	int col;

	if (!lp)
	{
		QSlog("EGLPNUM_TYPENAME_ILLlib_getbnd called without an lp");
		rval = 1;
		ILL_CLEANUP;
	}

	if (indx < 0 || indx > lp->O->nstruct)
	{
		QSlog("EGLPNUM_TYPENAME_ILLlib_getbnd called with bad indx: %d", indx);
		rval = 1;
		ILL_CLEANUP;
	}

	col = lp->O->structmap[indx];

	switch (lu)
	{
	case 'L':
		EGLPNUM_TYPENAME_EGlpNumCopy (*bnd, lp->O->lower[col]);
		break;
	case 'U':
		EGLPNUM_TYPENAME_EGlpNumCopy (*bnd, lp->O->upper[col]);
		break;
	default:
		QSlog("EGLPNUM_TYPENAME_ILLlib_getbnd called with lu: %c", lu);
		rval = 1;
		ILL_CLEANUP;
	}

CLEANUP:

	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLlib_getbnds_list ( 
	EGLPNUM_TYPENAME_lpinfo *lp,
	int num,
	int*collist, 
	EGLPNUM_TYPE *lower,
	EGLPNUM_TYPE *upper)
{
    int rval = 0;
    EGLPNUM_TYPENAME_ILLlpdata *qslp;
    int nstruct;
    int j, col;

    if (!lp) {
        QSlog("EGLPNUM_TYPENAME_ILLlib_getbnds_list called without an lp");
        rval = 1; ILL_CLEANUP;
    }

    qslp = lp->O;
    nstruct = qslp->nstruct;
		for (j = 0; j < num ; j++) {
		if(collist[j]<0|| collist[j] >= nstruct)
			{
				QSlog("EGLPNUM_TYPENAME_ILLlib_getbnds_list collist[%d] = %d out "
										"of range", j, collist[j]);
			}
			col = qslp->structmap[collist[j]];
			if (lower)
				EGLPNUM_TYPENAME_EGlpNumCopy(lower[j], qslp->lower[col]);
			if (upper)
				EGLPNUM_TYPENAME_EGlpNumCopy(upper[j], qslp->upper[col]);
		}

CLEANUP:

	EG_RETURN(rval);		
}


int EGLPNUM_TYPENAME_ILLlib_getbnds (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPE * lower,
	EGLPNUM_TYPE * upper)
{
	int rval = 0;
	EGLPNUM_TYPENAME_ILLlpdata *qslp;
	int nstruct;
	int j, col;

	if (!lp)
	{
		QSlog("EGLPNUM_TYPENAME_ILLlib_getbnd called without an lp");
		rval = 1;
		ILL_CLEANUP;
	}

	qslp = lp->O;
	nstruct = qslp->nstruct;

	for (j = 0; j < nstruct; j++)
	{
		col = qslp->structmap[j];
		if (lower)
			EGLPNUM_TYPENAME_EGlpNumCopy (lower[j], qslp->lower[col]);
		if (upper)
			EGLPNUM_TYPENAME_EGlpNumCopy (upper[j], qslp->upper[col]);
	}

CLEANUP:

	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLlib_strongbranch (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_price_info * pinf,
	int *candidatelist,
	int ncand,
	EGLPNUM_TYPE * xlist,
	EGLPNUM_TYPE * downpen,
	EGLPNUM_TYPE * uppen,
	int iterations,
	EGLPNUM_TYPE objbound,
	itcnt_t*itcnt)
{
	int rval = 0;
	int i, k, status, have_norms;
	int olditer = lp->maxiter;
	int nstruct = lp->O->nstruct;
	int nrows = lp->O->nrows;
	EGLPNUM_TYPE *myx = 0;
	EGLPNUM_TYPE xi, t, oldbnd;
	EGLPNUM_TYPENAME_price_info lpinf;
	EGLPNUM_TYPENAME_ILLlp_basis B, origB;

	EGLPNUM_TYPENAME_EGlpNumInitVar (lpinf.htrigger);
	EGLPNUM_TYPENAME_EGlpNumInitVar (xi);
	EGLPNUM_TYPENAME_EGlpNumInitVar (t);
	EGLPNUM_TYPENAME_EGlpNumInitVar (oldbnd);
	EGLPNUM_TYPENAME_EGlpNumZero (oldbnd);
	EGLPNUM_TYPENAME_ILLlp_basis_init (&B);
	EGLPNUM_TYPENAME_ILLlp_basis_init (&origB);
	EGLPNUM_TYPENAME_ILLprice_init_pricing_info (&lpinf);
	lpinf.dI_price = QS_PRICE_DSTEEP;
	lpinf.dII_price = QS_PRICE_DSTEEP;

	if (xlist == 0)
	{
		myx = EGLPNUM_TYPENAME_EGlpNumAllocArray (nstruct);
		rval = EGLPNUM_TYPENAME_ILLlib_get_x (lp, 0, myx);
		CHECKRVALG (rval, CLEANUP);
	}

	rval = EGLPNUM_TYPENAME_ILLlp_basis_alloc (&origB, nstruct, nrows);
	CHECKRVALG (rval, CLEANUP);

	rval = EGLPNUM_TYPENAME_ILLlib_getbasis (lp, origB.cstat, origB.rstat);
	CHECKRVALG (rval, CLEANUP);

	check_pinf (pinf, &have_norms);
	if (have_norms == 0)
	{
		origB.rownorms = EGLPNUM_TYPENAME_EGlpNumAllocArray (nrows);
		rval = EGLPNUM_TYPENAME_ILLlib_getrownorms (lp, pinf, origB.rownorms);
		CHECKRVALG (rval, CLEANUP);
	}
	else
	{
		lp->basisid = -1;
		rval = EGLPNUM_TYPENAME_ILLlib_optimize (lp, 0, &lpinf, DUAL_SIMPLEX, &status, 0, itcnt);
		CHECKRVALG (rval, CLEANUP);
	}

	rval = EGLPNUM_TYPENAME_ILLlp_basis_alloc (&B, nstruct, nrows);	/* Note: B and orgiB may */
	/* differ.               */
	CHECKRVALG (rval, CLEANUP);

	rval = EGLPNUM_TYPENAME_ILLlib_getbasis (lp, B.cstat, B.rstat);
	CHECKRVALG (rval, CLEANUP);
	B.rownorms = EGLPNUM_TYPENAME_EGlpNumAllocArray (nrows);

	if (have_norms == 0)
	{
		rval = EGLPNUM_TYPENAME_ILLlib_getrownorms (lp, pinf, B.rownorms);
		CHECKRVALG (rval, CLEANUP);
	}
	else
	{
		rval = EGLPNUM_TYPENAME_ILLlib_getrownorms (lp, &lpinf, B.rownorms);
		CHECKRVALG (rval, CLEANUP);
	}

	lp->maxiter = iterations;

	for (i = 0; i < ncand; i++)
	{
		k = candidatelist[i];
		rval = EGLPNUM_TYPENAME_ILLlib_getbnd (lp, k, 'U', &oldbnd);
		CHECKRVALG (rval, CLEANUP);
		if (xlist)
			EGLPNUM_TYPENAME_EGlpNumCopy (xi, xlist[i]);
		else
			EGLPNUM_TYPENAME_EGlpNumCopy (xi, myx[k]);
		EGLPNUM_TYPENAME_EGlpNumFloor (t, xi);
		if (EGLPNUM_TYPENAME_EGlpNumIsLessDbl (t, 0.1) && EGLPNUM_TYPENAME_EGlpNumIsGreaDbl (t, -0.1))
			EGLPNUM_TYPENAME_EGlpNumZero (t);

		rval = EGLPNUM_TYPENAME_ILLlib_chgbnd (lp, k, 'U', t);
		CHECKRVALG (rval, CLEANUP);

		rval = EGLPNUM_TYPENAME_ILLlib_optimize (lp, &B, &lpinf, DUAL_SIMPLEX, &status, 0, itcnt);
		CHECKRVALG (rval, CLEANUP);

		EGLPNUM_TYPENAME_EGlpNumCopy (downpen[i], lp->dobjval);
		rval = EGLPNUM_TYPENAME_ILLlib_chgbnd (lp, k, 'U', oldbnd);
		CHECKRVALG (rval, CLEANUP);

		rval = EGLPNUM_TYPENAME_ILLlib_getbnd (lp, k, 'L', &oldbnd);
		CHECKRVALG (rval, CLEANUP);
		EGLPNUM_TYPENAME_EGlpNumCeil (t, xi);
		if (EGLPNUM_TYPENAME_EGlpNumIsLessDbl (t, 1.1) && EGLPNUM_TYPENAME_EGlpNumIsGreaDbl (t, 0.9))
			EGLPNUM_TYPENAME_EGlpNumOne (t);
		rval = EGLPNUM_TYPENAME_ILLlib_chgbnd (lp, k, 'L', t);
		CHECKRVALG (rval, CLEANUP);

		rval = EGLPNUM_TYPENAME_ILLlib_optimize (lp, &B, &lpinf, DUAL_SIMPLEX, &status, 0, itcnt);
		CHECKRVALG (rval, CLEANUP);

		EGLPNUM_TYPENAME_EGlpNumCopy (uppen[i], lp->dobjval);
		rval = EGLPNUM_TYPENAME_ILLlib_chgbnd (lp, k, 'L', oldbnd);
		CHECKRVALG (rval, CLEANUP);
	}

	if (lp->O->objsense == EGLPNUM_TYPENAME_ILL_MAX)
	{

	}
	else
	{
		for (i = 0; i < ncand; i++)
		{
			if (EGLPNUM_TYPENAME_EGlpNumIsLess (objbound, downpen[i]))
				EGLPNUM_TYPENAME_EGlpNumCopy (downpen[i], objbound);
			if (EGLPNUM_TYPENAME_EGlpNumIsLess (objbound, uppen[i]))
				EGLPNUM_TYPENAME_EGlpNumCopy (uppen[i], objbound);
		}
	}

	/* Restore the old optimal solution */

	lp->maxiter = olditer;
	rval = EGLPNUM_TYPENAME_ILLlib_optimize (lp, &origB, pinf, DUAL_SIMPLEX, &status, 0, itcnt);
	CHECKRVALG (rval, CLEANUP);

CLEANUP:

	EGLPNUM_TYPENAME_EGlpNumClearVar (xi);
	EGLPNUM_TYPENAME_EGlpNumClearVar (t);
	EGLPNUM_TYPENAME_EGlpNumClearVar (oldbnd);
	lp->maxiter = olditer;
	EGLPNUM_TYPENAME_ILLprice_free_pricing_info (&lpinf);
	EGLPNUM_TYPENAME_ILLlp_basis_free (&B);
	EGLPNUM_TYPENAME_ILLlp_basis_free (&origB);
	if (xlist == 0)
		EGLPNUM_TYPENAME_EGlpNumFreeArray (myx);
	EGLPNUM_TYPENAME_EGlpNumClearVar (lpinf.htrigger);
	EG_RETURN (rval);
}

#define EXTRA_ROWS (100)
#define EXTRA_COLS (100)
#define EXTRA_MAT  (1000)

int EGLPNUM_TYPENAME_ILLlib_newrow (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_ILLlp_basis * B,
	const EGLPNUM_TYPE rhs,
	int sense,
	const EGLPNUM_TYPE range,
	const char *name)
{
	int rval = 0;

	rval = EGLPNUM_TYPENAME_ILLlib_addrow (lp, B, 0, 0, 0, rhs, sense, range, name);
	CHECKRVALG (rval, CLEANUP);

CLEANUP:

	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLlib_newrows (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_ILLlp_basis * B,
	int num,
	const EGLPNUM_TYPE * rhs,
	char *sense,
	const EGLPNUM_TYPE * range,
	const char **names)
{
	int rval = 0;
	int *rmatcnt = 0;
	int *rmatbeg = 0;
	int i;

	if (!num)
		ILL_CLEANUP;

	ILL_SAFE_MALLOC (rmatcnt, num, int);

	ILL_SAFE_MALLOC (rmatbeg, num, int);

	for (i = 0; i < num; i++)
	{
		rmatcnt[i] = 0;
		rmatbeg[i] = 0;
	}

	rval = EGLPNUM_TYPENAME_ILLlib_addrows (lp, B, num, rmatcnt, rmatbeg, 0, 0, rhs, sense,
												 range, names, 0);
	CHECKRVALG (rval, CLEANUP);

CLEANUP:

	ILL_IFFREE (rmatcnt, int);
	ILL_IFFREE (rmatbeg, int);

	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLlib_addrows (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_ILLlp_basis * B,
	int num,
	int *rmatcnt,
	int *rmatbeg,
	int *rmatind,
	const EGLPNUM_TYPE * rmatval,
	const EGLPNUM_TYPE * rhs,
	char *sense,
	const EGLPNUM_TYPE * range,
	const char **names,
	int *factorok)
{
	int rval = 0;
	int i, j, total, bsing;
	int *imap = 0;
	int *bbeg = 0;
	int *bcnt = 0;
	int *bindi = 0;
	int *rindi = 0;
	int *jstat = 0;
	EGLPNUM_TYPE *bval = 0;
	EGLPNUM_TYPE rng;
	int badfactor = 0;

	EGLPNUM_TYPENAME_EGlpNumInitVar (rng);

	if (B == 0 || B->rownorms == 0)
	{
		if (factorok)
			*factorok = 0;
	}

	if (B)
		EGLPNUM_TYPENAME_EGlpNumFreeArray (B->colnorms);

	if (B && B->rownorms && factorok && *factorok == 1)
	{
		int *structmap = lp->O->structmap;

		lp->matbeg = lp->O->A.matbeg;
		lp->matcnt = lp->O->A.matcnt;
		lp->matind = lp->O->A.matind;
		lp->matval = lp->O->A.matval;

		lp->nrows = lp->O->nrows;
		lp->ncols = lp->O->ncols;
		if (B->rownorms_size < lp->O->nrows + num)
			EGLPNUM_TYPENAME_EGlpNumReallocArray (&(B->rownorms), lp->O->nrows + num);

		ILL_SAFE_MALLOC (bcnt, num, int);
		ILL_SAFE_MALLOC (bbeg, num, int);
		ILL_SAFE_MALLOC (imap, lp->O->nstruct, int);

		ILL_SAFE_MALLOC (jstat, lp->ncols, int);

		for (i = 0; i < lp->ncols; i++)
		{
			jstat[i] = -1;
		}
		for (i = 0; i < lp->O->nstruct; i++)
		{
			jstat[structmap[i]] = i;
		}

		for (i = 0; i < lp->O->nstruct; i++)
		{
			imap[i] = -1;
		}
		for (i = 0; i < lp->O->nrows; i++)
		{
			if (jstat[lp->baz[i]] != -1)
			{
				imap[jstat[lp->baz[i]]] = i;
			}
		}

		for (i = 0, total = 0; i < num; i++)
		{
			bcnt[i] = 0;
			bbeg[i] = total;
			for (j = 0; j < rmatcnt[i]; j++)
			{
				if (imap[rmatind[rmatbeg[i] + j]] != -1)
				{
					bcnt[i]++;
					total++;
				}
			}
		}
		if (total)
		{
			ILL_SAFE_MALLOC (bindi, total, int);

			bval = EGLPNUM_TYPENAME_EGlpNumAllocArray (total);
		}
		for (i = 0, total = 0; i < num; i++)
		{
			for (j = 0; j < rmatcnt[i]; j++)
			{
				if (imap[rmatind[rmatbeg[i] + j]] != -1)
				{
					EGLPNUM_TYPENAME_EGlpNumCopy (bval[total], rmatval[rmatbeg[i] + j]);
					bindi[total] = imap[rmatind[rmatbeg[i] + j]];
					total++;
				}
			}
		}

		rval = EGLPNUM_TYPENAME_ILLprice_get_new_rownorms (lp, num, B->rownorms + lp->O->nrows,
																			bcnt, bbeg, bindi, bval);
		CHECKRVALG (rval, CLEANUP);

		ILL_IFFREE (bcnt, int);
		ILL_IFFREE (bbeg, int);
		ILL_IFFREE (bindi, int);

		EGLPNUM_TYPENAME_EGlpNumFreeArray (bval);
		ILL_IFFREE (imap, int);

		badfactor = 1;
	}

	for (i = 0; i < num; i++)
	{
		if (range)
			EGLPNUM_TYPENAME_EGlpNumCopy (rng, range[i]);
		else
			EGLPNUM_TYPENAME_EGlpNumZero (rng);
		if (names)
		{
			rval = EGLPNUM_TYPENAME_ILLlib_addrow (lp, B, rmatcnt[i], rmatind + rmatbeg[i],
														rmatval + rmatbeg[i], rhs[i], sense[i], rng,
														names[i]);
		}
		else
		{
			rval = EGLPNUM_TYPENAME_ILLlib_addrow (lp, B, rmatcnt[i], rmatind + rmatbeg[i],
														rmatval + rmatbeg[i], rhs[i], sense[i], rng, 0);
		}
		CHECKRVALG (rval, CLEANUP);
	}


	if (B && B->rownorms && (factorok && *factorok == 0))
	{
		lp->matbeg = lp->O->A.matbeg;
		lp->matcnt = lp->O->A.matcnt;
		lp->matind = lp->O->A.matind;
		lp->matval = lp->O->A.matval;
		lp->nrows = lp->O->nrows;
		lp->ncols = lp->O->ncols;
		lp->bz = lp->O->rhs;
		lp->nnbasic = lp->ncols - lp->nrows;

		rval = EGLPNUM_TYPENAME_ILLbasis_load (lp, B);
		CHECKRVALG (rval, CLEANUP);

		if (lp->f)
			EGLPNUM_TYPENAME_ILLfactor_free_factor_work (lp->f);

		rval = EGLPNUM_TYPENAME_ILLbasis_factor (lp, &bsing);
		CHECKRVALG (rval, CLEANUP);
		if (bsing)
			MESSAGE (__QS_SB_VERB, "Singular Basis found!");
		*factorok = 1;

		if (B->rownorms_size < lp->O->nrows)
			EGLPNUM_TYPENAME_EGlpNumReallocArray (&(B->rownorms), lp->O->nrows);

		ILL_SAFE_MALLOC (rindi, lp->O->nrows /* num */ , int);

		for (i = 0; i < num; i++)
		{
			rindi[i] = lp->O->nrows - num + i;
		}

		rval = EGLPNUM_TYPENAME_ILLprice_get_dsteep_norms (lp, num, rindi,
																			B->rownorms + lp->O->nrows - num);
		CHECKRVALG (rval, CLEANUP);
	}

	if (factorok != 0 && badfactor == 1)
	{
		*factorok = 0;
	}


CLEANUP:

	ILL_IFFREE (bcnt, int);
	ILL_IFFREE (bbeg, int);
	ILL_IFFREE (bindi, int);

	EGLPNUM_TYPENAME_EGlpNumFreeArray (bval);
	ILL_IFFREE (imap, int);
	ILL_IFFREE (jstat, int);
	ILL_IFFREE (rindi, int);

	EGLPNUM_TYPENAME_EGlpNumClearVar (rng);
	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLlib_addrow (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_ILLlp_basis * B,
	int cnt,
	int *ind,
	const EGLPNUM_TYPE * val,
	const EGLPNUM_TYPE rhs,
	int sense,
	const EGLPNUM_TYPE range,
	const char *name)
{
	int rval = 0;
	EGLPNUM_TYPENAME_ILLlpdata *qslp;
	EGLPNUM_TYPENAME_ILLmatrix *A;
	int i, nrows, ncols;
	char buf[ILL_namebufsize];
	int tind[1];
	EGLPNUM_TYPE tval[1];
	int *tempind = 0;
	int pind, hit;

	EGLPNUM_TYPENAME_EGlpNumInitVar (tval[0]);

	if (!lp)
	{
		QSlog("EGLPNUM_TYPENAME_ILLlib_addrow called without an lp");
		rval = 1;
		ILL_CLEANUP;
	}

	qslp = lp->O;
	A = &qslp->A;

	if (qslp->rA)
	{															/* After an addrow call, needs to be updated */
		EGLPNUM_TYPENAME_ILLlp_rows_clear (qslp->rA);
		ILL_IFFREE (qslp->rA, EGLPNUM_TYPENAME_ILLlp_rows);
	}

	if (qslp->sinfo)
	{															/* Presolve LP is no longer valid, free the data */
		EGLPNUM_TYPENAME_ILLlp_sinfo_free (qslp->sinfo);
		ILL_IFFREE (qslp->sinfo, EGLPNUM_TYPENAME_ILLlp_sinfo);
	}

	nrows = qslp->nrows;
	ncols = qslp->ncols;

	/* If the row has a range, create the rangeval array if needed  */

	if (sense == 'R' && !(qslp->rangeval) && qslp->rowsize > 0)
	{
		qslp->rangeval = EGLPNUM_TYPENAME_EGlpNumAllocArray (qslp->rowsize);
		for (i = 0; i < qslp->nrows; i++)
		{
			EGLPNUM_TYPENAME_EGlpNumZero (qslp->rangeval[i]);
		}
	}

	/* Add the row to the row structures */

	if (qslp->rowsize < nrows + 1)
	{
		EGLPNUM_TYPENAME_EGlpNumReallocArray (&(qslp->rhs), qslp->rowsize + EXTRA_ROWS);
		qslp->sense = EGrealloc (qslp->sense,
														 sizeof (char) * (qslp->rowsize + EXTRA_ROWS));
		//rval = ILLutil_reallocrus_count ((void **) &(qslp->sense),
		//                                 qslp->rowsize + EXTRA_ROWS, sizeof (char));
		//CHECKRVALG(rval,CLEANUP);

		qslp->rowmap = EGrealloc (qslp->rowmap,
															sizeof (int) * (qslp->rowsize + EXTRA_ROWS));
		//rval = ILLutil_reallocrus_count ((void **) &(qslp->rowmap),
		//                                 qslp->rowsize + EXTRA_ROWS, sizeof (int));
		//CHECKRVALG(rval,CLEANUP);

		if (qslp->rangeval || sense == 'R')
			EGLPNUM_TYPENAME_EGlpNumReallocArray (&(qslp->rangeval), qslp->rowsize + EXTRA_ROWS);

		qslp->rownames = EGrealloc (qslp->rownames,
																sizeof (char *) * (qslp->rowsize + EXTRA_ROWS));
		//rval = ILLutil_reallocrus_count ((void **) &(qslp->rownames),
		//                                 qslp->rowsize + EXTRA_ROWS,
		//                                 sizeof (char *));
		//CHECKRVALG(rval,CLEANUP);
		qslp->rowsize += EXTRA_ROWS;
	}

	EGLPNUM_TYPENAME_EGlpNumCopy (qslp->rhs[nrows], rhs);
	qslp->sense[nrows] = sense;
	qslp->rowmap[nrows] = ncols;	/* this will be the new logical */
	if (qslp->rangeval)
	{
		if (sense == 'R')
			EGLPNUM_TYPENAME_EGlpNumCopy (qslp->rangeval[nrows], range);
		else
			EGLPNUM_TYPENAME_EGlpNumZero (qslp->rangeval[nrows]);
	}
	ILL_FAILtrue (qslp->rownames == NULL, "must always be non NULL");
	EGLPNUM_TYPENAME_ILLlib_findName (qslp, 1 /*row */ , name, nrows, buf);
	ILL_UTIL_STR (qslp->rownames[nrows], buf);
	ILLsymboltab_register (&qslp->rowtab, buf, qslp->nrows, &pind, &hit);
	ILL_FAILfalse (hit == 0, "must be new");


	/* Add the logical variable to the column structures */

	if (qslp->colsize < ncols + 1)
	{
		EGLPNUM_TYPENAME_EGlpNumReallocArray (&(qslp->lower), qslp->colsize + EXTRA_COLS);
		EGLPNUM_TYPENAME_EGlpNumReallocArray (&(qslp->upper), qslp->colsize + EXTRA_COLS);
		EGLPNUM_TYPENAME_EGlpNumReallocArray (&(qslp->obj), qslp->colsize + EXTRA_COLS);
		qslp->colsize += EXTRA_COLS;
	}

	EGLPNUM_TYPENAME_EGlpNumZero (qslp->obj[ncols]);
	EGLPNUM_TYPENAME_EGlpNumZero (qslp->lower[ncols]);
	if (sense == 'E')
	{
		EGLPNUM_TYPENAME_EGlpNumZero (qslp->upper[ncols]);	/* Artificial */
	}
	else if (sense == 'R')
	{
		EGLPNUM_TYPENAME_EGlpNumCopy (qslp->upper[ncols], range);	/* Range      */
	}
	else
	{
		EGLPNUM_TYPENAME_EGlpNumCopy (qslp->upper[ncols], EGLPNUM_TYPENAME_ILL_MAXDOUBLE);	/* Slack      */
	}

	/* Add new row and new logical col to matrix */

	/* Need to map the structural indices to their proper place */

	if (cnt)
	{
		ILL_SAFE_MALLOC (tempind, cnt, int);

		for (i = 0; i < cnt; i++)
		{
			tempind[i] = qslp->structmap[ind[i]];
		}
	}

	rval = matrix_addrow (A, cnt, tempind, val);
	CHECKRVALG (rval, CLEANUP);

	tind[0] = nrows;
	EGLPNUM_TYPENAME_EGlpNumOne (*tval);
	if (sense == 'G' || sense == 'R')
		EGLPNUM_TYPENAME_EGlpNumSign (*tval);

	rval = matrix_addcol (A, 1, tind, tval);
	CHECKRVALG (rval, CLEANUP);

	if (B != 0)
	{
		B->rstat = EGrealloc (B->rstat, sizeof (char) * (nrows + 1));
		//rval = ILLutil_reallocrus_count ((void **) &(B->rstat), nrows + 1,
		//                                 sizeof (char));
		//CHECKRVALG(rval,CLEANUP);
		B->rstat[nrows] = QS_ROW_BSTAT_BASIC;
	}

#if 0
	lp->basisid = -1;							/* To get optimizer to reload the basis */
#endif

	qslp->ncols++;
	qslp->nrows++;
	qslp->nzcount += (cnt + 1);

	if (B != 0)
	{
		B->nrows++;
	}

CLEANUP:
	ILL_IFFREE (tempind, int);

	EGLPNUM_TYPENAME_EGlpNumClearVar (tval[0]);
	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLlib_delrows (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_ILLlp_basis * B,
	EGLPNUM_TYPENAME_ILLlp_cache * C,
	int num,
	int *dellist,
	int *basis_ok,
	int *cache_ok)
{
	int rval = 0;
	int i, j, k, nrows, ncols, nstruct, spot, dk, bok = 0, cok = 0;
	EGLPNUM_TYPENAME_ILLlpdata *qslp;
	EGLPNUM_TYPENAME_ILLmatrix *A;
	char *rowmark = 0;
	char *colmark = 0;
	int *newrowindex = 0;
	int *newcolindex = 0;
	int *ind, *beg, *cnt;
	EGLPNUM_TYPE *val;

	if (!lp)
	{
		QSlog("EGLPNUM_TYPENAME_ILLlib_delrows called without an lp");
		rval = 1;
		ILL_CLEANUP;
	}

	if (num <= 0)
	{
		if (basis_ok)
			*basis_ok = 1;
		if (cache_ok)
			*cache_ok = 1;
		ILL_CLEANUP;
	}

	if (basis_ok)
		*basis_ok = 0;
	if (cache_ok)
		*cache_ok = 0;

	qslp = lp->O;
	A = &qslp->A;

	if (qslp->rA)
	{															/* After a delrow call, needs to be updated */
		EGLPNUM_TYPENAME_ILLlp_rows_clear (qslp->rA);
		ILL_IFFREE (qslp->rA, EGLPNUM_TYPENAME_ILLlp_rows);
	}

	nrows = A->matrows;
	ncols = A->matcols;
	ind = A->matind;
	beg = A->matbeg;
	cnt = A->matcnt;
	val = A->matval;
	nstruct = qslp->nstruct;

	ILL_SAFE_MALLOC (rowmark, nrows, char);

	for (i = 0; i < nrows; i++)
	{
		rowmark[i] = 0;
	}
	for (i = 0; i < num; i++)
	{
		rowmark[dellist[i]] = 1;
	}


	/* Try to update the basis */

	if (B)
	{
		bok = 1;
		cok = 1;
		for (i = 0; i < num; i++)
		{
			j = dellist[i];
			if (B->rstat[j] == QS_ROW_BSTAT_LOWER ||
					B->rstat[j] == QS_ROW_BSTAT_UPPER)
			{
				bok = 0;
				break;
			}
			if (C && EGLPNUM_TYPENAME_EGlpNumIsLess (EGLPNUM_TYPENAME_DFEAS_TOLER, C->pi[j]))
			{
/*
                QSlog("XXXX: Postive pi (%f) at basic row", C->pi[j]);
*/
				cok = 0;
			}
		}
		if (bok == 1)
		{
			EGLPNUM_TYPENAME_EGlpNumFreeArray (B->colnorms);
			if (B->rownorms)
			{
				for (i = 0, k = 0; i < nstruct; i++)
				{
					if (B->cstat[i] == QS_COL_BSTAT_BASIC)
						k++;
				}
				for (i = 0, j = k; i < nrows; i++)
				{
					if (B->rstat[i] == QS_ROW_BSTAT_BASIC)
					{
						if (rowmark[i] == 0)
						{
							EGLPNUM_TYPENAME_EGlpNumCopy (B->rownorms[k++], B->rownorms[j]);
						}
						j++;
					}
				}
				if (k != nrows - num)
				{
					QSlog("error in  EGLPNUM_TYPENAME_ILLlib_delrows");
					rval = 1;
					ILL_CLEANUP;
				}
			}

			for (i = 0, j = 0; i < nrows; i++)
			{
				if (rowmark[i] == 0)
				{
					B->rstat[j++] = B->rstat[i];
				}
			}
			B->nrows = j;

			if (C && cok == 1)
			{
				for (i = 0, j = 0; i < nrows; i++)
				{
					if (rowmark[i] == 0)
					{
						EGLPNUM_TYPENAME_EGlpNumCopy (C->pi[j], C->pi[i]);
						EGLPNUM_TYPENAME_EGlpNumCopy (C->slack[j++], C->slack[i]);
					}
				}
				C->nrows = j;
				if (cache_ok)
					*cache_ok = 1;
			}
			if (basis_ok)
				*basis_ok = 1;
		}
	}

	ILL_SAFE_MALLOC (newrowindex, nrows, int);


	/* Delete the marked rows */

	ILL_FAILtrue (qslp->rownames == NULL, "must always be non NULL");
	for (i = 0, j = 0; i < nrows; i++)
	{
		if (rowmark[i] == 0)
		{
			if (i != j)
			{
				EGLPNUM_TYPENAME_EGlpNumCopy (qslp->rhs[j], qslp->rhs[i]);
				qslp->sense[j] = qslp->sense[i];
				if (qslp->rangeval)
					EGLPNUM_TYPENAME_EGlpNumCopy (qslp->rangeval[j], qslp->rangeval[i]);
				if (qslp->rownames)
					qslp->rownames[j] = qslp->rownames[i];
			}
			newrowindex[i] = j++;
		}
		else
		{
			if (qslp->rownames)
			{
				rval = ILLsymboltab_delete (&qslp->rowtab, qslp->rownames[i]);
				CHECKRVALG (rval, CLEANUP);
				ILL_IFFREE (qslp->rownames[i], char);
			}
		}
	}


	/* Delete the logicals */

	ILL_SAFE_MALLOC (colmark, ncols, char);

	for (i = 0; i < ncols; i++)
	{
		colmark[i] = 0;
	}
	for (i = 0; i < num; i++)
	{
		colmark[qslp->rowmap[dellist[i]]] = 1;
	}

	rval = delcols_work (lp, colmark);
	CHECKRVALG (rval, CLEANUP);

	A->matcols -= num;
	qslp->ncols -= num;


	/* Pack the rowmap  */

	for (i = 0, j = 0; i < nrows; i++)
	{
		if (rowmark[i] == 0)
		{
			qslp->rowmap[j++] = qslp->rowmap[i];
		}
	}

	/* Remove the entries to deleted rows, and update the indices */

	for (i = 0; i < ncols - num; i++)
	{
		dk = 0;
		spot = beg[i];
		for (j = 0; j < cnt[i]; j++)
		{
			if (rowmark[ind[beg[i] + j]] == 1)
			{
				dk++;
			}
			else
			{
				EGLPNUM_TYPENAME_EGlpNumCopy (val[spot], val[beg[i] + j]);
				ind[spot] = newrowindex[ind[beg[i] + j]];
				spot++;
			}
		}
		for (; spot < beg[i] + cnt[i]; spot++)
		{
			ind[spot] = -1;
		}

		cnt[i] -= dk;
		if (cnt[i] == 0)
		{
			ind[beg[i]] = 1;					/* we always mark the empty cols */
		}
	}

	A->matrows -= num;
	qslp->nrows -= num;

#if 0
	lp->basisid = -1;							/* To get optimizer to reload the basis */
#endif

	/* if the base is OK, we MUST load the status variables again */
	if(bok)
	{
		rval = EGLPNUM_TYPENAME_ILLbasis_load( lp, B);
		CHECKRVALG (rval, CLEANUP);
	}
CLEANUP:

	ILL_IFFREE (rowmark, char);
	ILL_IFFREE (colmark, char);
	ILL_IFFREE (newcolindex, int);
	ILL_IFFREE (newrowindex, int);

	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLlib_delcols (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_ILLlp_basis * B,
	int num,
	int *dellist,
	int *basis_ok)
{
	int rval = 0;
	int i, j, bok = 0, ncols;
	char *colmark = 0;
	EGLPNUM_TYPENAME_ILLlpdata *qslp;

	if (!lp)
	{
		QSlog("EGLPNUM_TYPENAME_ILLlib_delcols called without an lp");
		rval = 1;
		ILL_CLEANUP;
	}

	if (basis_ok)
		*basis_ok = 0;

	if (num <= 0)
	{
		*basis_ok = 1;
		ILL_CLEANUP;
	}

	qslp = lp->O;
	ncols = qslp->A.matcols;

	if (qslp->rA)
	{															/* After a delcol call, needs to be updated */
		EGLPNUM_TYPENAME_ILLlp_rows_clear (qslp->rA);
		ILL_IFFREE (qslp->rA, EGLPNUM_TYPENAME_ILLlp_rows);
	}

	ILL_SAFE_MALLOC (colmark, ncols, char);

	for (i = 0; i < ncols; i++)
	{
		colmark[i] = 0;
	}
	for (i = 0; i < num; i++)
	{
		colmark[qslp->structmap[dellist[i]]] = 1;
	}

	if (B)
	{
		B->nstruct -= num;
		bok = 1;
		for (i = 0; i < num; i++)
		{
			j = dellist[i];
			if (B->cstat[j] == QS_COL_BSTAT_BASIC)
			{
				bok = 0;
				//QSlog("BONG");
				break;
			}
		}
		if (bok == 1)
		{
			EGLPNUM_TYPENAME_EGlpNumFreeArray (B->colnorms);
			for (i = 0, j = 0; i < qslp->nstruct; i++)
			{
				if (colmark[qslp->structmap[i]] == 0)
				{
					B->cstat[j++] = B->cstat[i];
				}
			}
			if (basis_ok)
				*basis_ok = 1;
		}
	}

	rval = delcols_work (lp, colmark);
	CHECKRVALG (rval, CLEANUP);


	qslp->A.matcols -= num;
	qslp->ncols -= num;
	qslp->nstruct -= num;

	/* if the base is OK, we MUST load the status variables again */
	if(bok)
	{
		rval = EGLPNUM_TYPENAME_ILLbasis_load( lp, B);
		CHECKRVALG (rval, CLEANUP);
	}
#if 0
	lp->basisid = -1;							/* To get optimizer to reload the basis */
#endif

CLEANUP:

	ILL_IFFREE (colmark, char);

	EG_RETURN (rval);
}

static int matrix_getcoef (
	EGLPNUM_TYPENAME_ILLmatrix *A, 
	int row,
	int col,
	EGLPNUM_TYPE*val)
{
	int i;
	int rval = 0;
	if (row >= A->matrows || row < 0)
	{
		QSlog("illegal row index in matrix_getcoef");
		rval= 1;
		ILL_CLEANUP;
	}

	if (col >= A->matcols || col < 0)
	{
		QSlog("illegal col index in matrix_getcoef");
		rval= 1;
		ILL_CLEANUP;
	}

	/* by default value is zero */
	EGLPNUM_TYPENAME_EGlpNumZero(*val);
	for (i = A->matbeg[col]; i < A->matbeg[col] + A->matcnt[col]; i++)
	{
		if (A->matind[i] == row)
		{
			EGLPNUM_TYPENAME_EGlpNumCopy(*val, A->matval[i]);
			ILL_CLEANUP;
		}
	}

CLEANUP:

	EG_RETURN(rval);
}

static int delcols_work (
	EGLPNUM_TYPENAME_lpinfo * lp,
	char *colmark)
{
	int rval = 0;
	int i, j, k, nrows, ncols;
	EGLPNUM_TYPENAME_ILLlpdata *qslp;
	EGLPNUM_TYPENAME_ILLmatrix *A;
	int *newcolindex = 0;
	int *ind, *beg, *cnt;

	/* Allows logicals to be deleted, to handle call from delcols. */

	qslp = lp->O;
	A = &qslp->A;
	nrows = A->matrows;
	ncols = A->matcols;
	ind = A->matind;
	beg = A->matbeg;
	cnt = A->matcnt;

	ILL_SAFE_MALLOC (newcolindex, ncols, int);

	/* Delete the columns */

	for (i = 0, j = 0; i < ncols; i++)
	{
		if (colmark[i] == 0)
		{
			if (i != j)
			{
				beg[j] = beg[i];
				cnt[j] = cnt[i];
				EGLPNUM_TYPENAME_EGlpNumCopy (qslp->obj[j], qslp->obj[i]);
				EGLPNUM_TYPENAME_EGlpNumCopy (qslp->lower[j], qslp->lower[i]);
				EGLPNUM_TYPENAME_EGlpNumCopy (qslp->upper[j], qslp->upper[i]);
			}
			newcolindex[i] = j++;
		}
		else
		{
			for (k = 0; k < cnt[i]; k++)
			{
				ind[beg[i] + k] = -1;
			}
			newcolindex[i] = -1;
		}
	}

	/* Update the struct arrays */

	for (i = 0, j = 0; i < qslp->nstruct; i++)
	{
		k = qslp->structmap[i];
		if (colmark[k] == 0)
		{
			qslp->structmap[j] = newcolindex[k];
			qslp->colnames[j] = qslp->colnames[i];
			if (qslp->intmarker)
				qslp->intmarker[j] = qslp->intmarker[i];
			j++;
		}
		else
		{
			rval = ILLsymboltab_delete (&qslp->coltab, qslp->colnames[i]);
			CHECKRVALG (rval, CLEANUP);
			ILL_IFFREE (qslp->colnames[i], char);
		}
	}

	/* Update the rowmap: note if logicals deleted, map will be -1 */

	for (i = 0; i < nrows; i++)
	{
		qslp->rowmap[i] = newcolindex[qslp->rowmap[i]];
	}

CLEANUP:

	ILL_IFFREE (newcolindex, int);

	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLlib_getcoef (
	EGLPNUM_TYPENAME_lpinfo *lp,
	int rowindex,
	int colindex,
	EGLPNUM_TYPE* coef)
{ 
	int rval = 0;
	EGLPNUM_TYPENAME_ILLlpdata *qslp;
	EGLPNUM_TYPENAME_ILLmatrix *A;
	int nrows, nstruct, j;
	if (!lp)
	{
		QSlog("EGLPNUM_TYPENAME_ILLlib_chgcoef called without an lp");
		rval = 1;
		ILL_CLEANUP;
	}

	qslp = lp->O;
	A = &qslp->A;
	nrows = qslp->nrows;
	nstruct = qslp->nstruct;

	if (rowindex < 0 || rowindex >= nrows || colindex < 0 || colindex >= nstruct)
	{
		QSlog("EGLPNUM_TYPENAME_ILLlib_getcoef called with out-of-range index");
		rval = 1;
		ILL_CLEANUP;
	}
	
	j = qslp->structmap[colindex];
	rval = matrix_getcoef (A, rowindex, j, coef);
	CHECKRVALG(rval, CLEANUP);

CLEANUP:

	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLlib_chgcoef (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int rowindex,
	int colindex,
	EGLPNUM_TYPE coef)
{
	int rval = 0;
	EGLPNUM_TYPENAME_ILLlpdata *qslp;
	EGLPNUM_TYPENAME_ILLmatrix *A;
	int nrows, nstruct, j;

	if (!lp)
	{
		QSlog("EGLPNUM_TYPENAME_ILLlib_chgcoef called without an lp");
		rval = 1;
		ILL_CLEANUP;
	}

	qslp = lp->O;
	A = &qslp->A;

	nrows = qslp->nrows;
	nstruct = qslp->nstruct;

	if (rowindex < 0 || rowindex >= nrows || colindex < 0 || colindex >= nstruct)
	{
		QSlog("EGLPNUM_TYPENAME_ILLlib_chgcoef called with out-of-range index");
		rval = 1;
		ILL_CLEANUP;
	}

	if (qslp->rA)
	{															/* After a chgcoef call, needs to be updated */
		EGLPNUM_TYPENAME_ILLlp_rows_clear (qslp->rA);
		ILL_IFFREE (qslp->rA, EGLPNUM_TYPENAME_ILLlp_rows);
	}

	if (qslp->sinfo)
	{															/* Presolve LP is no longer valid, free the data */
		EGLPNUM_TYPENAME_ILLlp_sinfo_free (qslp->sinfo);
		ILL_IFFREE (qslp->sinfo, EGLPNUM_TYPENAME_ILLlp_sinfo);
	}

	j = qslp->structmap[colindex];

	rval = matrix_addcoef (lp, A, rowindex, j, coef);
	CHECKRVALG (rval, CLEANUP);

CLEANUP:

	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLlib_chgsense (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int num,
	int *rowlist,
	char *sense)
{
	int rval = 0;
	int i, j, k;
	EGLPNUM_TYPENAME_ILLlpdata *qslp = lp->O;
	EGLPNUM_TYPENAME_ILLmatrix *A = &(lp->O->A);

	for (i = 0; i < num; i++)
	{
		j = qslp->rowmap[rowlist[i]];
		if (A->matcnt[j] != 1)
		{
			QSlog("logical variable is not a singleton");
			rval = 1;
			ILL_CLEANUP;
		}
		k = A->matbeg[j];
		switch (sense[i])
		{
		case 'R':									/* Range constraint, we will set its upper bound
																 once we call EGLPNUM_TYPENAME_QSchange_range, by default it 
																 will be zero, i.e. an equation. */
			qslp->sense[rowlist[i]] = 'R';
			EGLPNUM_TYPENAME_EGlpNumZero(qslp->lower[j]);
			EGLPNUM_TYPENAME_EGlpNumZero(qslp->upper[j]);
			EGLPNUM_TYPENAME_EGlpNumOne(A->matval[k]);
			break;
		case 'E':									/* Artificial */
			qslp->sense[rowlist[i]] = 'E';
			EGLPNUM_TYPENAME_EGlpNumZero (qslp->lower[j]);
			EGLPNUM_TYPENAME_EGlpNumZero (qslp->upper[j]);
			EGLPNUM_TYPENAME_EGlpNumOne (A->matval[k]);
			break;
		case 'G':									/* Surplus   */
			qslp->sense[rowlist[i]] = 'G';
			EGLPNUM_TYPENAME_EGlpNumZero (qslp->lower[j]);
			EGLPNUM_TYPENAME_EGlpNumCopy (qslp->upper[j], EGLPNUM_TYPENAME_ILL_MAXDOUBLE);
			EGLPNUM_TYPENAME_EGlpNumOne (A->matval[k]);
			EGLPNUM_TYPENAME_EGlpNumSign (A->matval[k]);
			break;
		case 'L':									/* Slack     */
			qslp->sense[rowlist[i]] = 'L';
			EGLPNUM_TYPENAME_EGlpNumZero (qslp->lower[j]);
			EGLPNUM_TYPENAME_EGlpNumCopy (qslp->upper[j], EGLPNUM_TYPENAME_ILL_MAXDOUBLE);
			EGLPNUM_TYPENAME_EGlpNumOne (A->matval[k]);
			break;
		default:
			QSlog("illegal sense %c in EGLPNUM_TYPENAME_ILLlib_chgsense", sense[i]);
			rval = 1;
			ILL_CLEANUP;
		}
	}

CLEANUP:

	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLlib_getsenses (
	EGLPNUM_TYPENAME_lpinfo *lp,
	char *senses)
{
	EGLPNUM_TYPENAME_ILLlpdata *qslp;
	int nrows, i;
	int rval = 0;

	if (!lp) {
		QSlog("ILLlib_getsense called without an LP");
		rval = 1;
		ILL_CLEANUP;
	}

	qslp = lp->O;
	nrows = qslp->nrows;

	for (i = 0; i < nrows; i++)
	{
		senses[i] = qslp->sense[i];
	} 

CLEANUP:

    EG_RETURN(rval);
}

int EGLPNUM_TYPENAME_ILLlib_newcol (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_ILLlp_basis * B,
	const EGLPNUM_TYPE obj,
	const EGLPNUM_TYPE lower,
	const EGLPNUM_TYPE upper,
	const char *name,
	int factorok)
{
	int rval = 0;

	rval = EGLPNUM_TYPENAME_ILLlib_addcol (lp, B, 0, 0, 0, obj, lower, upper, name, factorok);
	CHECKRVALG (rval, CLEANUP);

CLEANUP:

	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLlib_newcols (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_ILLlp_basis * B,
	int num,
	EGLPNUM_TYPE * obj,
	EGLPNUM_TYPE * lower,
	EGLPNUM_TYPE * upper,
	const char **names,
	int factorok)
{
	int rval = 0;
	int *cmatcnt = 0;
	int *cmatbeg = 0;
	int i;

	ILL_SAFE_MALLOC (cmatcnt, num, int);

	ILL_SAFE_MALLOC (cmatbeg, num, int);

	for (i = 0; i < num; i++)
	{
		cmatcnt[i] = 0;
		cmatbeg[i] = 0;
	}

	rval = EGLPNUM_TYPENAME_ILLlib_addcols (lp, B, num, cmatcnt, cmatbeg, 0,
												 0, obj, lower, upper, names, factorok);
	CHECKRVALG (rval, CLEANUP);


CLEANUP:

	ILL_IFFREE (cmatcnt, int);
	ILL_IFFREE (cmatbeg, int);

	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLlib_addcols (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_ILLlp_basis * B,
	int num,
	int *cmatcnt,
	int *cmatbeg,
	int *cmatind,
	EGLPNUM_TYPE * cmatval,
	EGLPNUM_TYPE * obj,
	EGLPNUM_TYPE * lower,
	EGLPNUM_TYPE * upper,
	const char **names,
	int factorok)
{
	int rval = 0;
	int i;

	for (i = 0; i < num; i++)
	{
		if (names)
		{
			rval = EGLPNUM_TYPENAME_ILLlib_addcol (lp, B, cmatcnt[i], cmatind + cmatbeg[i],
														cmatval + cmatbeg[i], obj[i], lower[i],
														upper[i], names[i], factorok);
		}
		else
		{
			rval = EGLPNUM_TYPENAME_ILLlib_addcol (lp, B, cmatcnt[i], cmatind + cmatbeg[i],
														cmatval + cmatbeg[i], obj[i], lower[i],
														upper[i], 0, factorok);
		}
		CHECKRVALG (rval, CLEANUP);
	}

CLEANUP:

	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLlib_addcol (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_ILLlp_basis * B,
	int cnt,
	int *ind,
	EGLPNUM_TYPE * val,
	const EGLPNUM_TYPE obj,
	const EGLPNUM_TYPE lower,
	const EGLPNUM_TYPE upper,
	const char *name,
	int factorok)
{
	int rval = 0;
	EGLPNUM_TYPENAME_ILLlpdata *qslp;
	EGLPNUM_TYPENAME_ILLmatrix *A;
	int ncols;
	char buf[ILL_namebufsize];
	int pind, hit;
	EGLPNUM_TYPE l, u;

	EGLPNUM_TYPENAME_EGlpNumInitVar (l);
	EGLPNUM_TYPENAME_EGlpNumInitVar (u);

	if (!lp)
	{
		QSlog("EGLPNUM_TYPENAME_ILLlib_addcol called without an lp");
		rval = 1;
		ILL_CLEANUP;
	}

	qslp = lp->O;
	A = &qslp->A;
	ncols = qslp->ncols;

	if (qslp->rA)
	{															/* After an addcol call, needs to be updated */
		EGLPNUM_TYPENAME_ILLlp_rows_clear (qslp->rA);
		ILL_IFFREE (qslp->rA, EGLPNUM_TYPENAME_ILLlp_rows);
	}

	if (qslp->sinfo)
	{															/* Presolve LP is no longer valid, free the data */
		EGLPNUM_TYPENAME_ILLlp_sinfo_free (qslp->sinfo);
		ILL_IFFREE (qslp->sinfo, EGLPNUM_TYPENAME_ILLlp_sinfo);
	}


	/* Add the new variable to the column structures */

	if (qslp->colsize < ncols + 1)
	{
		EGLPNUM_TYPENAME_EGlpNumReallocArray (&(qslp->lower), qslp->colsize + EXTRA_COLS);
		EGLPNUM_TYPENAME_EGlpNumReallocArray (&(qslp->upper), qslp->colsize + EXTRA_COLS);
		EGLPNUM_TYPENAME_EGlpNumReallocArray (&(qslp->obj), qslp->colsize + EXTRA_COLS);
		qslp->colsize += EXTRA_COLS;
	}

	EGLPNUM_TYPENAME_EGlpNumCopy (qslp->obj[ncols], obj);
	EGLPNUM_TYPENAME_EGlpNumCopy (qslp->lower[ncols], lower);
	EGLPNUM_TYPENAME_EGlpNumCopy (qslp->upper[ncols], upper);

	/*  Add the variable to the structural arrays */

	if (qslp->structsize < qslp->nstruct + 1)
	{
		qslp->structmap = EGrealloc (qslp->structmap,
																 sizeof (int) * (qslp->structsize +
																								 EXTRA_COLS));
		//rval = ILLutil_reallocrus_count ((void **) &(qslp->structmap),
		//                                 qslp->structsize + EXTRA_COLS,
		//                                 sizeof (int));
		//CHECKRVALG(rval,CLEANUP);

		qslp->colnames = EGrealloc (qslp->colnames,
																sizeof (char *) * (qslp->structsize +
																									 EXTRA_COLS));
		//rval = ILLutil_reallocrus_count ((void **) &(qslp->colnames),
		//                                 qslp->structsize + EXTRA_COLS,
		//                                 sizeof (char *));
		//CHECKRVALG(rval,CLEANUP);

		if (qslp->intmarker)
		{
			qslp->intmarker = EGrealloc (qslp->intmarker,
																	 sizeof (char) * (qslp->structsize +
																										EXTRA_COLS));
			//rval = ILLutil_reallocrus_count ((void **) &(qslp->intmarker),
			//                                 qslp->structsize + EXTRA_COLS,
			//                                 sizeof (char));
			//CHECKRVALG(rval,CLEANUP);
		}
		qslp->structsize += EXTRA_COLS;
	}

	qslp->structmap[qslp->nstruct] = ncols;
	if (qslp->intmarker)
	{
		/* NOTE: If we want to add integer variables, this is the place. */
		qslp->intmarker[qslp->nstruct] = (char) 0;
	}

	ILL_FAILtrue (qslp->colnames == NULL, "must always be non NULL");
	EGLPNUM_TYPENAME_ILLlib_findName (qslp, 0 /*isRow */ , name, qslp->nstruct, buf);
	ILLsymboltab_register (&qslp->coltab, buf, qslp->nstruct, &pind, &hit);
	ILL_FAILfalse ((pind == qslp->nstruct) && (hit == 0), "must be new");
	ILL_UTIL_STR (qslp->colnames[qslp->nstruct], buf);


	/*  Add col to the matrix */

	rval = matrix_addcol (A, cnt, ind, val);
	CHECKRVALG (rval, CLEANUP);


	if (B)
	{
		B->cstat = EGrealloc (B->cstat, sizeof (char) * (qslp->nstruct + 1));
		//rval = ILLutil_reallocrus_count ((void **) &(B->cstat),
		//                                 qslp->nstruct + 1, sizeof (char));
		//CHECKRVALG(rval,CLEANUP);
		if (EGLPNUM_TYPENAME_EGlpNumIsEqual (lower, EGLPNUM_TYPENAME_ILL_MINDOUBLE, EGLPNUM_TYPENAME_oneLpNum) &&
				EGLPNUM_TYPENAME_EGlpNumIsEqual (upper, EGLPNUM_TYPENAME_ILL_MAXDOUBLE, EGLPNUM_TYPENAME_oneLpNum))
		{
			B->cstat[qslp->nstruct] = QS_COL_BSTAT_FREE;
		}
		else if (EGLPNUM_TYPENAME_EGlpNumIsEqual (upper, EGLPNUM_TYPENAME_ILL_MAXDOUBLE, EGLPNUM_TYPENAME_oneLpNum))
		{
			B->cstat[qslp->nstruct] = QS_COL_BSTAT_LOWER;
		}
		//else if (lower == EGLPNUM_TYPENAME_ILL_MAXDOUBLE)
		else if (EGLPNUM_TYPENAME_EGlpNumIsEqual (lower, EGLPNUM_TYPENAME_ILL_MAXDOUBLE, EGLPNUM_TYPENAME_oneLpNum))
		{
			B->cstat[qslp->nstruct] = QS_COL_BSTAT_UPPER;
		}
		else
		{
			/*l = fabs (lower);
			 * u = fabs (upper); */
			EGLPNUM_TYPENAME_EGlpNumCopyAbs (l, lower);
			EGLPNUM_TYPENAME_EGlpNumCopyAbs (u, upper);
			if (EGLPNUM_TYPENAME_EGlpNumIsLess (l, u))
			{
				B->cstat[qslp->nstruct] = QS_COL_BSTAT_LOWER;
			}
			else
			{
				B->cstat[qslp->nstruct] = QS_COL_BSTAT_UPPER;
			}
		}

		/* UPDATE THE PINF PRIMAL NORMS */
		EGLPNUM_TYPENAME_EGlpNumFreeArray (B->colnorms);
	}

	if (factorok == 0)
	{
#if 0
		lp->basisid = -1;						/* To get optimizer to reload the basis */
#endif
	}
	else
	{
		if (!lp->nbaz || !lp->vindex || !lp->vstat)
		{
			QSlog("ERROR: factorok set without a current basis");
			rval = 1;
			ILL_CLEANUP;
		}

		lp->nbaz = EGrealloc (lp->nbaz, sizeof (int) * (qslp->nstruct + 1));
		//rval = ILLutil_reallocrus_count ((void **) &(lp->nbaz),
		//                                 qslp->nstruct + 1, sizeof (int));
		//CHECKRVALG(rval,CLEANUP);

		lp->vindex = EGrealloc (lp->vindex, sizeof (int) * (qslp->ncols + 1));
		//rval = ILLutil_reallocrus_count ((void **) &(lp->vindex),
		//                                 qslp->ncols + 1, sizeof (int));
		//CHECKRVALG(rval,CLEANUP);

		lp->vstat = EGrealloc (lp->vstat, sizeof (int) * (qslp->ncols + 1));
		//rval = ILLutil_reallocrus_count ((void **) &(lp->vstat),
		//                                 qslp->ncols + 1, sizeof (int));


		lp->nbaz[qslp->nstruct] = qslp->ncols;
		lp->vindex[qslp->ncols] = qslp->nstruct;

		if (EGLPNUM_TYPENAME_EGlpNumIsEqual (lower, EGLPNUM_TYPENAME_ILL_MINDOUBLE, EGLPNUM_TYPENAME_oneLpNum) &&
				EGLPNUM_TYPENAME_EGlpNumIsEqual (upper, EGLPNUM_TYPENAME_ILL_MAXDOUBLE, EGLPNUM_TYPENAME_oneLpNum))
		{
			lp->vstat[qslp->ncols] = STAT_ZERO;
		}
		else if (EGLPNUM_TYPENAME_EGlpNumIsEqual (upper, EGLPNUM_TYPENAME_ILL_MAXDOUBLE, EGLPNUM_TYPENAME_oneLpNum))
		{
			lp->vstat[qslp->ncols] = STAT_LOWER;
		}
		//else if (lower == EGLPNUM_TYPENAME_ILL_MAXDOUBLE)
		else if (EGLPNUM_TYPENAME_EGlpNumIsEqual (lower, EGLPNUM_TYPENAME_ILL_MAXDOUBLE, EGLPNUM_TYPENAME_oneLpNum))
		{
			lp->vstat[qslp->ncols] = STAT_UPPER;
		}
		else
		{
			/*l = fabs (lower);
			 * u = fabs (upper); */
			EGLPNUM_TYPENAME_EGlpNumCopyAbs (l, lower);
			EGLPNUM_TYPENAME_EGlpNumCopyAbs (u, upper);
			if (EGLPNUM_TYPENAME_EGlpNumIsLess (l, u))
			{
				lp->vstat[qslp->ncols] = STAT_LOWER;
			}
			else
			{
				lp->vstat[qslp->ncols] = STAT_UPPER;
			}
		}
	}


	qslp->ncols++;
	qslp->nstruct++;
	(qslp->nzcount) += cnt;

	if (B)
	{
		B->nstruct++;
	}

CLEANUP:
	EGLPNUM_TYPENAME_EGlpNumClearVar (l);
	EGLPNUM_TYPENAME_EGlpNumClearVar (u);
	EG_RETURN (rval);
}

static int matrix_addrow (
	EGLPNUM_TYPENAME_ILLmatrix * A,
	int rowcnt,
	int *rowind,
	const EGLPNUM_TYPE * rowval)
{
	int rval = 0;
	int i, j, k, ind, memo, stop, delta = 0;

	/* matsize will be the length of the array.                   */
	/* matfree will keep track of the free space at end of array. */

	for (i = 0; i < rowcnt; i++)
	{
		if (rowind[i] >= A->matcols || rowind[i] < 0)
		{
			QSlog("illegal col index in matrix_addrow");
			rval = 1;
			ILL_CLEANUP;
		}
	}

	for (i = 0; i < rowcnt; i++)
	{
		j = rowind[i];
		if (A->matcnt[j] > 0 &&
				(A->matbeg[j] + A->matcnt[j] + 1 > A->matsize ||
				 A->matind[A->matbeg[j] + A->matcnt[j]] != -1))
		{
			delta += (A->matcnt[j] + 2);	/* 1 for the new coef and 1 for */
			/* an extra space               */
		}
	}

	if (delta < A->matfree)
	{
		for (i = 0; i < rowcnt; i++)
		{
			j = rowind[i];
			if (A->matcnt[j] == 0)
			{
				A->matind[A->matbeg[j]] = A->matrows;
				EGLPNUM_TYPENAME_EGlpNumCopy (A->matval[A->matbeg[j]], rowval[i]);
				A->matcnt[j] = 1;
			}
			else if (A->matind[A->matbeg[j] + A->matcnt[j]] == -1)
			{
				/* Since A->matfree is positive, we know that we are not */
				/* sitting at the end of the array.                      */
				A->matind[A->matbeg[j] + A->matcnt[j]] = A->matrows;
				EGLPNUM_TYPENAME_EGlpNumCopy (A->matval[A->matbeg[j] + A->matcnt[j]], rowval[i]);
				if ((A->matbeg[j] + A->matcnt[j]) == (A->matsize - A->matfree))
				{
					A->matfree--;					/* at end of used space */
				}
				(A->matcnt[j])++;
			}
			else
			{
				ind = A->matsize - A->matfree + 1;	/* leave space for -1 */
				memo = ind;
				stop = A->matbeg[j] + A->matcnt[j];
				for (k = A->matbeg[j]; k < stop; k++)
				{
					if (ind >= A->matsize)
					{
						QSlog("WHAT: %d, %d", A->matsize, ind);
						exit (1);
					}
					A->matind[ind] = A->matind[k];
					EGLPNUM_TYPENAME_EGlpNumCopy (A->matval[ind], A->matval[k]);
					A->matind[k] = -1;
					ind++;
				}
				A->matind[ind] = A->matrows;
				EGLPNUM_TYPENAME_EGlpNumCopy (A->matval[ind], rowval[i]);
				A->matbeg[j] = memo;
				(A->matcnt[j])++;
				(A->matfree) -= (A->matcnt[j] + 1);
			}
		}
	}
	else
	{
		rval = matrix_addrow_end (A, A->matrows, rowcnt, rowind, rowval);
		CHECKRVALG (rval, CLEANUP);
	}
	A->matrows++;

CLEANUP:

	EG_RETURN (rval);
}

static int matrix_addrow_end (
	EGLPNUM_TYPENAME_ILLmatrix * A,
	int row,
	int rowcnt,
	int *rowind,
	const EGLPNUM_TYPE * rowval)
{
	int rval = 0;
	int i, j, k, start, stop, total;
	int *newbeg = 0;
	int *newind = 0;
	EGLPNUM_TYPE *newval = 0;
	int ncols = A->matcols;

	if (A->matcolsize > 0)
	{
		ILL_SAFE_MALLOC (newbeg, A->matcolsize, int);
	}
	ILL_SAFE_MALLOC (newind, A->matsize + rowcnt + EXTRA_MAT, int);

	newval = EGLPNUM_TYPENAME_EGlpNumAllocArray (A->matsize + rowcnt + EXTRA_MAT);

	A->matsize += (rowcnt + EXTRA_MAT);

	for (i = 0; i < rowcnt; i++)
	{
		A->matcnt[rowind[i]]++;
	}
	for (total = 0, j = 0; j < ncols; j++)
	{
		newbeg[j] = total;
		if (A->matcnt[j] > 0)
			total += A->matcnt[j];
		else
			total += 1;
	}
	for (i = 0; i < rowcnt; i++)
	{
		A->matcnt[rowind[i]]--;
	}
	for (j = total; j < A->matsize; j++)
	{
		newind[j] = -1;
	}
	A->matfree = A->matsize - total;

	for (j = 0; j < ncols; j++)
	{
		if (A->matcnt[j] > 0)
		{
			stop = A->matbeg[j] + A->matcnt[j];
			start = newbeg[j];
			for (k = A->matbeg[j]; k < stop; k++)
			{
				newind[start] = A->matind[k];
				EGLPNUM_TYPENAME_EGlpNumCopy (newval[start], A->matval[k]);
				start++;
			}
		}
		else
		{
			newind[newbeg[j]] = 1;
		}
	}
	for (i = 0; i < rowcnt; i++)
	{
		j = rowind[i];
		newind[newbeg[j] + A->matcnt[j]] = row;
		EGLPNUM_TYPENAME_EGlpNumCopy (newval[newbeg[j] + A->matcnt[j]], rowval[i]);
		(A->matcnt[j])++;
	}

	ILL_IFFREE (A->matbeg, int);
	ILL_IFFREE (A->matind, int);

	EGLPNUM_TYPENAME_EGlpNumFreeArray (A->matval);

	A->matbeg = newbeg;
	A->matind = newind;
	A->matval = newval;

CLEANUP:

	if (rval)
	{
		ILL_IFFREE (newbeg, int);
		ILL_IFFREE (newind, int);

		EGLPNUM_TYPENAME_EGlpNumFreeArray (newval);
	}

	EG_RETURN (rval);
}

static int matrix_addcoef (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_ILLmatrix * A,
	int row,
	int col,
	EGLPNUM_TYPE val)
{
	int i, k, delta, ind, stop, memo;
	int tind[1];
	EGLPNUM_TYPE tval[1];
	int rval = 0;

	EGLPNUM_TYPENAME_EGlpNumInitVar (tval[0]);
	EGLPNUM_TYPENAME_EGlpNumCopy (tval[0], val);

	if (row >= A->matrows || row < 0)
	{
		QSlog("illegal row index in matrix_addcoef");
		rval = 1;
		ILL_CLEANUP;
	}

	if (col >= A->matcols || col < 0)
	{
		QSlog("illegal col index in matrix_addcoef");
		rval = 1;
		ILL_CLEANUP;
	}

	for (i = A->matbeg[col]; i < A->matbeg[col] + A->matcnt[col]; i++)
	{
		if (A->matind[i] == row)
		{
			EGLPNUM_TYPENAME_EGlpNumCopy (A->matval[i], val);
			ILL_CLEANUP;
		}
	}

	/* The coef is new, we need to add it to A */

	lp->O->nzcount++;
	delta = A->matcnt[col] + 2;

	if (A->matcnt[col] == 0)
	{
		/* First entry, always a free space */
		A->matind[A->matbeg[col]] = row;
		EGLPNUM_TYPENAME_EGlpNumCopy (A->matval[A->matbeg[col]], val);
		A->matcnt[col] = 1;
	}
	else if (A->matbeg[col] + A->matcnt[col] < A->matsize &&
					 A->matind[A->matbeg[col] + A->matcnt[col]] == -1)
	{
		/* Free space in the column */
		A->matind[A->matbeg[col] + A->matcnt[col]] = row;
		EGLPNUM_TYPENAME_EGlpNumCopy (A->matval[A->matbeg[col] + A->matcnt[col]], val);
		if ((A->matbeg[col] + A->matcnt[col]) == (A->matsize - A->matfree))
		{
			A->matfree--;
		}
		(A->matcnt[col])++;
	}
	else if (A->matfree > delta)
	{
		/* Enough space to move column to end of array */
		ind = A->matsize - A->matfree + 1;
		memo = ind;
		stop = A->matbeg[col] + A->matcnt[col];
		for (k = A->matbeg[col]; k < stop; k++)
		{
			A->matind[ind] = A->matind[k];
			EGLPNUM_TYPENAME_EGlpNumCopy (A->matval[ind], A->matval[k]);
			A->matind[k] = -1;
			ind++;
		}
		A->matind[ind] = row;
		EGLPNUM_TYPENAME_EGlpNumCopy (A->matval[ind], val);

		A->matbeg[col] = memo;
		(A->matcnt[col])++;
		(A->matfree) -= (A->matcnt[col] + 1);
	}
	else
	{
		/* Need to malloc space to move column to end of array */

		tind[0] = col;

		rval = matrix_addrow_end (A, row, 1, tind, tval);
		CHECKRVALG (rval, CLEANUP);
	}

CLEANUP:

	EGLPNUM_TYPENAME_EGlpNumClearVar (tval[0]);
	EG_RETURN (rval);
}

static int matrix_addcol (
	EGLPNUM_TYPENAME_ILLmatrix * A,
	int colcnt,
	int *colind,
	EGLPNUM_TYPE * colval)
{
	int rval = 0;
	int i, ind;

	for (i = 0; i < colcnt; i++)
	{
		if (colind[i] >= A->matrows || colind[i] < 0)
		{
			QSlog("illegal row index in matrix_addcol");
			rval = 1;
			ILL_CLEANUP;
		}
	}

	if (A->matcolsize < A->matcols + 1)
	{
		A->matbeg =
			EGrealloc (A->matbeg, sizeof (int) * (A->matcolsize + EXTRA_COLS));
		//rval = ILLutil_reallocrus_count ((void **) &(A->matbeg),
		//                                 A->matcolsize + EXTRA_COLS, sizeof (int));
		//CHECKRVALG(rval,CLEANUP);

		A->matcnt =
			EGrealloc (A->matcnt, sizeof (int) * (A->matcolsize + EXTRA_COLS));
		//rval = ILLutil_reallocrus_count ((void **) &(A->matcnt),
		//                                 A->matcolsize + EXTRA_COLS, sizeof (int));
		//CHECKRVALG(rval,CLEANUP);

		(A->matcolsize) += EXTRA_COLS;
	}

	if (A->matfree < colcnt + 1)
	{
		A->matind = EGrealloc (A->matind,
													 sizeof (int) * (A->matsize + colcnt + EXTRA_MAT +
																					 1));
		//rval = ILLutil_reallocrus_count ((void **) &(A->matind),
		//                                 A->matsize + colcnt + EXTRA_MAT + 1,
		//                                 sizeof (int));
		//CHECKRVALG(rval,CLEANUP);
		EGLPNUM_TYPENAME_EGlpNumReallocArray (&(A->matval), A->matsize + colcnt + EXTRA_MAT + 1);

		for (i = 0; i < colcnt + EXTRA_MAT + 1; i++)
		{
			A->matind[A->matsize + i] = -1;
		}
		A->matsize += (colcnt + EXTRA_MAT + 1);
		A->matfree += (colcnt + EXTRA_MAT + 1);
	}

	ind = A->matsize - A->matfree;
	A->matbeg[A->matcols] = ind;
	A->matcnt[A->matcols] = colcnt;
	if (colcnt == 0)
	{
		A->matind[ind] = 1;					/* Dummy value to stop columns from stealing */
		/* this space in addrows.                    */
		A->matfree -= 1;
	}
	else
	{
		for (i = 0; i < colcnt; i++)
		{
			EGLPNUM_TYPENAME_EGlpNumCopy (A->matval[ind], colval[i]);
			A->matind[ind] = colind[i];
			ind++;
		}
		A->matfree -= colcnt;
	}
	A->matcols++;

CLEANUP:

	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLlib_getrows (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int num,
	int *rowlist,
	int **rowcnt,
	int **rowbeg,
	int **rowind,
	EGLPNUM_TYPE ** rowval,
	EGLPNUM_TYPE ** rhs,
	char **sense,
	EGLPNUM_TYPE ** range,
	char ***names)
{
	int rval = 0;
	int *allbeg = 0;
	int *allcnt = 0;
	int *allind = 0;
	EGLPNUM_TYPE *allval = 0;
	int i, row, k, start, stop, len, tcnt, cnt = 0;
	EGLPNUM_TYPENAME_ILLlpdata *qslp;
	EGLPNUM_TYPENAME_ILLlp_rows lprows;

	if (rowcnt) *rowcnt = 0;
	if (rowbeg) *rowbeg = 0;
	if (rowind) *rowind = 0;
	if (rowval) *rowval = 0;
	if (rhs) *rhs = 0;
	if (range) *range = 0;
	if (sense) *sense = 0;
	if (names) *names = 0;

	if (!lp)
	{
		QSlog("EGLPNUM_TYPENAME_ILLlib_getrows called without an LP");
		rval = 1;
		ILL_CLEANUP;
	}

	if (!num)
		ILL_CLEANUP;

	qslp = lp->O;

	rval = EGLPNUM_TYPENAME_ILLlp_rows_init (&lprows, qslp, 0);
	CHECKRVALG (rval, CLEANUP);
	allbeg = lprows.rowbeg;
	allcnt = lprows.rowcnt;
	allind = lprows.rowind;
	allval = lprows.rowval;

	for (i = 0; i < num; i++)
	{
		cnt += allcnt[rowlist[i]];
	}

	if (rowcnt)
	{
		ILL_SAFE_MALLOC (*rowcnt, num, int);

		for (i = 0; i < num; i++)
		{
			(*rowcnt)[i] = allcnt[rowlist[i]];
		}
	}

	if (rowbeg)
	{
		ILL_SAFE_MALLOC (*rowbeg, num, int);

		tcnt = 0;
		for (i = 0; i < num; i++)
		{
			(*rowbeg)[i] = tcnt;
			tcnt += allcnt[rowlist[i]];
		}
	}

	if (cnt && rowind)
	{
		ILL_SAFE_MALLOC (*rowind, cnt, int);

		tcnt = 0;
		for (i = 0; i < num; i++)
		{
			row = rowlist[i];
			start = allbeg[row];
			stop = start + allcnt[row];
			for (k = start; k < stop; k++)
			{
				(*rowind)[tcnt++] = allind[k];
			}
		}
	}

	if (cnt && rowval)
	{
		*rowval = EGLPNUM_TYPENAME_EGlpNumAllocArray (cnt);
		tcnt = 0;
		for (i = 0; i < num; i++)
		{
			row = rowlist[i];
			start = allbeg[row];
			stop = start + allcnt[row];
			for (k = start; k < stop; k++)
			{
				EGLPNUM_TYPENAME_EGlpNumCopy ((*rowval)[tcnt++], allval[k]);
			}
		}
	}

	if (rhs)
	{
		*rhs = EGLPNUM_TYPENAME_EGlpNumAllocArray (num);
		for (i = 0; i < num; i++)
		{
			EGLPNUM_TYPENAME_EGlpNumCopy ((*rhs)[i], qslp->rhs[rowlist[i]]);
		}
	}

	if (range)
	{
		*range = EGLPNUM_TYPENAME_EGlpNumAllocArray(num);
		if(qslp->rangeval)
		{
			for(i = 0; i < num ; i++)
			{
				EGLPNUM_TYPENAME_EGlpNumCopy((*range)[i], qslp->rangeval[rowlist[i]]);
			}
		}
		else
		{
			for(i = 0; i < num ; i++)
			{
				EGLPNUM_TYPENAME_EGlpNumZero((*range)[i]);
			}
		}
	}

	if (sense)
	{
		ILL_SAFE_MALLOC (*sense, num, char);

		for (i = 0; i < num; i++)
		{
			(*sense)[i] = qslp->sense[rowlist[i]];
		}
	}

	if (names)
	{
		if (qslp->rownames == 0)
		{
			QSlog("LP does not have row names");
			rval = 1;
			ILL_CLEANUP;
		}
		ILL_SAFE_MALLOC (*names, num, char *);

		for (i = 0; i < num; i++)
		{
			(*names)[i] = 0;
		}
		for (i = 0; i < num; i++)
		{
			len = strlen (qslp->rownames[rowlist[i]]) + 1;
			ILL_SAFE_MALLOC ((*names)[i], len, char);

			strcpy ((*names)[i], qslp->rownames[rowlist[i]]);
		}
	}

CLEANUP:

	ILL_IFFREE (allbeg, int);
	ILL_IFFREE (allcnt, int);
	ILL_IFFREE (allind, int);

	EGLPNUM_TYPENAME_EGlpNumFreeArray (allval);

	if (rval)
	{
		if (rowcnt)
			ILL_IFFREE (*rowcnt, int);

		if (rowbeg)
			ILL_IFFREE (*rowbeg, int);

		if (rowind)
			ILL_IFFREE (*rowind, int);

		if (rowval)
			EGLPNUM_TYPENAME_EGlpNumFreeArray (*rowval);
		if (rhs)
			EGLPNUM_TYPENAME_EGlpNumFreeArray (*rhs);
		if (sense)
			ILL_IFFREE (*sense, char);

		if (names && (*names))
		{
			for (i = 0; i < num; i++)
			{
				ILL_IFFREE ((*names)[i], char);
			}
			ILL_IFFREE (*names, char *);
		}
	}

	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLlib_getcols (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int num,
	int *collist,
	int **colcnt,
	int **colbeg,
	int **colind,
	EGLPNUM_TYPE ** colval,
	EGLPNUM_TYPE ** obj,
	EGLPNUM_TYPE ** lower,
	EGLPNUM_TYPE ** upper,
	char ***names)
{
	int rval = 0;
	int i, col, k, start, stop, len, tcnt, cnt = 0;
	EGLPNUM_TYPENAME_ILLlpdata *qslp;
	EGLPNUM_TYPENAME_ILLmatrix *A;
	int *tlist = 0;

	if (colcnt)
		*colcnt = 0;
	if (colbeg)
		*colbeg = 0;
	if (colind)
		*colind = 0;
	if (colval)
		*colval = 0;
	if (lower)
		*lower = 0;
	if (upper)
		*upper = 0;
	if (obj)
		*obj = 0;
	if (names)
		*names = 0;

	if (!lp)
	{
		QSlog("EGLPNUM_TYPENAME_ILLlib_getcols called without an LP");
		rval = 1;
		ILL_CLEANUP;
	}

	if (!num)
		ILL_CLEANUP;

	qslp = lp->O;
	A = &(qslp->A);

	ILL_SAFE_MALLOC (tlist, num, int);

	for (i = 0; i < num; i++)
	{
		tlist[i] = qslp->structmap[collist[i]];
	}

	for (i = 0; i < num; i++)
	{
		cnt += A->matcnt[tlist[i]];
	}

	if (colcnt)
	{
		ILL_SAFE_MALLOC (*colcnt, num, int);

		for (i = 0; i < num; i++)
		{
			(*colcnt)[i] = A->matcnt[tlist[i]];
		}
	}

	if (colbeg)
	{
		ILL_SAFE_MALLOC (*colbeg, num, int);

		tcnt = 0;
		for (i = 0; i < num; i++)
		{
			(*colbeg)[i] = tcnt;
			tcnt += A->matcnt[tlist[i]];
		}
	}

	if (cnt && colind)
	{
		ILL_SAFE_MALLOC (*colind, cnt, int);

		tcnt = 0;
		for (i = 0; i < num; i++)
		{
			col = tlist[i];
			start = A->matbeg[col];
			stop = start + A->matcnt[col];
			for (k = start; k < stop; k++)
			{
				(*colind)[tcnt++] = A->matind[k];
			}
		}
	}

	if (cnt && colval)
	{
		*colval = EGLPNUM_TYPENAME_EGlpNumAllocArray (cnt);
		tcnt = 0;
		for (i = 0; i < num; i++)
		{
			col = tlist[i];
			start = A->matbeg[col];
			stop = start + A->matcnt[col];
			for (k = start; k < stop; k++)
			{
				EGLPNUM_TYPENAME_EGlpNumCopy ((*colval)[tcnt++], A->matval[k]);
			}
		}
	}

	if (obj)
	{
		*obj = EGLPNUM_TYPENAME_EGlpNumAllocArray (num);
		for (i = 0; i < num; i++)
		{
			EGLPNUM_TYPENAME_EGlpNumCopy ((*obj)[i], qslp->obj[tlist[i]]);
		}
	}

	if (lower)
	{
		*lower = EGLPNUM_TYPENAME_EGlpNumAllocArray (num);
		for (i = 0; i < num; i++)
		{
			EGLPNUM_TYPENAME_EGlpNumCopy ((*lower)[i], qslp->lower[tlist[i]]);
		}
	}

	if (upper)
	{
		*upper = EGLPNUM_TYPENAME_EGlpNumAllocArray (num);
		for (i = 0; i < num; i++)
		{
			EGLPNUM_TYPENAME_EGlpNumCopy ((*upper)[i], qslp->upper[tlist[i]]);
		}
	}

	if (names)
	{
		if (qslp->colnames == 0)
		{
			QSlog("LP does not have col names");
			rval = 1;
			ILL_CLEANUP;
		}
		ILL_SAFE_MALLOC (*names, num, char *);

		for (i = 0; i < num; i++)
		{
			(*names)[i] = 0;
		}
		for (i = 0; i < num; i++)
		{
			len = strlen (qslp->colnames[collist[i]]) + 1;
			ILL_SAFE_MALLOC ((*names)[i], len, char);

			strcpy ((*names)[i], qslp->colnames[collist[i]]);
		}
	}

CLEANUP:

	if (rval)
	{
		if (colcnt)
			ILL_IFFREE (*colcnt, int);

		if (colbeg)
			ILL_IFFREE (*colbeg, int);

		if (colind)
			ILL_IFFREE (*colind, int);

		if (colval)
			EGLPNUM_TYPENAME_EGlpNumFreeArray (*colval);
		if (obj)
			EGLPNUM_TYPENAME_EGlpNumFreeArray (*obj);
		if (lower)
			EGLPNUM_TYPENAME_EGlpNumFreeArray (*lower);
		if (upper)
			EGLPNUM_TYPENAME_EGlpNumFreeArray (*upper);
		if (names && (*names))
		{
			for (i = 0; i < num; i++)
			{
				ILL_IFFREE ((*names)[i], char);
			}
			ILL_IFFREE (*names, char *);
		}
	}
	ILL_IFFREE (tlist, int);

	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLlib_getobj_list (
	EGLPNUM_TYPENAME_lpinfo *lp,
	int num,
	int* collist,
	EGLPNUM_TYPE* obj)
{
	const int*const structmap = lp->O->structmap;
	EGLPNUM_TYPENAME_ILLlpdata *qslp;
	int nstruct, j, col;
	int rval = 0;
	
	if (!lp)
	{
		QSlog("EGLPNUM_TYPENAME_ILLlib_getobj_list called without an LP");
		rval = 1;
		ILL_CLEANUP;
	}
	
	qslp = lp->O;
	nstruct = qslp->nstruct;
	
	for (j = 0; j < num; j++)
	{
		col = collist[j];
		if(col<0 || col >= nstruct)
		{
			QSlog("EGLPNUM_TYPENAME_ILLlib_getobj_list collist[%d] = %d outside"
									" valid range", j, col);
			rval = 1;
			ILL_CLEANUP;
		}
		EGLPNUM_TYPENAME_EGlpNumCopy(obj[j],qslp->obj[structmap[col]]);
	}  

CLEANUP:

	EG_RETURN (rval);
}


int EGLPNUM_TYPENAME_ILLlib_getobj (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPE * obj)
{
	const int*const structmap = lp->O->structmap;
	EGLPNUM_TYPENAME_ILLlpdata *qslp;
	int nstruct, j;
	int rval = 0;

	if (!lp)
	{
		QSlog("EGLPNUM_TYPENAME_ILLlib_getobj called without an LP");
		rval = 1;
		ILL_CLEANUP;
	}

	qslp = lp->O;
	nstruct = qslp->nstruct;

	for (j = 0; j < nstruct; j++)
	{
		EGLPNUM_TYPENAME_EGlpNumCopy (obj[j], qslp->obj[structmap[j]]);
	}

CLEANUP:

	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLlib_chgobj (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int indx,
	EGLPNUM_TYPE coef)
{
	int rval = 0;
	int col;

	if (!lp)
	{
		QSlog("EGLPNUM_TYPENAME_ILLlib_chgobj called without an lp");
		rval = 1;
		ILL_CLEANUP;
	}

	if (indx < 0 || indx >= lp->O->nstruct)
	{
		QSlog("EGLPNUM_TYPENAME_ILLlib_chgrhs called with bad indx: %d", indx);
		rval = 1;
		ILL_CLEANUP;
	}

	if (lp->O->sinfo)
	{															/* Presolve LP is no longer valid, free the data */
		EGLPNUM_TYPENAME_ILLlp_sinfo_free (lp->O->sinfo);
		ILL_IFFREE (lp->O->sinfo, EGLPNUM_TYPENAME_ILLlp_sinfo);
	}

	col = lp->O->structmap[indx];
	EGLPNUM_TYPENAME_EGlpNumCopy (lp->O->obj[col], coef);

CLEANUP:

	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLlib_getrhs (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPE * rhs)
{
	EGLPNUM_TYPENAME_ILLlpdata *qslp;
	int nrows, i;
	int rval = 0;

	if (!lp)
	{
		QSlog("EGLPNUM_TYPENAME_ILLlib_getrhs called without an LP");
		rval = 1;
		ILL_CLEANUP;
	}

	qslp = lp->O;
	nrows = qslp->nrows;

	for (i = 0; i < nrows; i++)
	{
		EGLPNUM_TYPENAME_EGlpNumCopy (rhs[i], qslp->rhs[i]);
	}

CLEANUP:

	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLlib_chgrange (
	EGLPNUM_TYPENAME_lpinfo *lp,
	int indx,
	EGLPNUM_TYPE coef)
{
	register int i;
	int rval = 0;
	EGLPNUM_TYPENAME_ILLlpdata *qslp;
	
	if (!lp)
	{
		QSlog("EGLPNUM_TYPENAME_ILLlib_chgrhs called without an lp");
		rval = 1;
		ILL_CLEANUP;
	}
	
	if (indx < 0 || indx >= lp->O->nrows)
	{
		QSlog("EGLPNUM_TYPENAME_ILLlib_chgrhs called with bad indx: %d", indx);
		rval = 1;
		ILL_CLEANUP;
	}
	
	if (lp->O->sinfo)
	{/* Presolve LP is no longer valid, free the data */
		EGLPNUM_TYPENAME_ILLlp_sinfo_free (lp->O->sinfo);
		ILL_IFFREE (lp->O->sinfo, EGLPNUM_TYPENAME_ILLlp_sinfo);
	}
	
	qslp = lp->O;
	if(qslp->rangeval == 0)
	{
		qslp->rangeval = EGLPNUM_TYPENAME_EGlpNumAllocArray(qslp->rowsize);
		for( i = qslp->nrows ; i-- ; )
		{
			EGLPNUM_TYPENAME_EGlpNumZero(qslp->rangeval[i]);
		}
	}
	
	if(qslp->sense[indx] != 'R')
	{
		QSlog("setting range for non-range constraint");
		rval = 1;
		ILL_CLEANUP;
	}
	
	EGLPNUM_TYPENAME_EGlpNumCopy(qslp->rangeval[indx], coef);

CLEANUP:

	EG_RETURN (rval);
}


int EGLPNUM_TYPENAME_ILLlib_chgrhs (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int indx,
	EGLPNUM_TYPE coef)
{
	int rval = 0;

	if (!lp)
	{
		QSlog("EGLPNUM_TYPENAME_ILLlib_chgrhs called without an lp");
		rval = 1;
		ILL_CLEANUP;
	}

	if (indx < 0 || indx >= lp->O->nrows)
	{
		QSlog("EGLPNUM_TYPENAME_ILLlib_chgrhs called with bad indx: %d", indx);
		rval = 1;
		ILL_CLEANUP;
	}

	if (lp->O->sinfo)
	{															/* Presolve LP is no longer valid, free the data */
		EGLPNUM_TYPENAME_ILLlp_sinfo_free (lp->O->sinfo);
		ILL_IFFREE (lp->O->sinfo, EGLPNUM_TYPENAME_ILLlp_sinfo);
	}

	EGLPNUM_TYPENAME_EGlpNumCopy (lp->O->rhs[indx], coef);

CLEANUP:

	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLlib_rownames (
	EGLPNUM_TYPENAME_lpinfo * lp,
	char **rownames)
{
	EGLPNUM_TYPENAME_ILLlpdata *qslp;
	int nrows, len, i, rcount = 0;
	int rval = 0;

	if (!lp)
	{
		QSlog("EGLPNUM_TYPENAME_ILLlib_rownames called without an LP");
		rval = 1;
		ILL_CLEANUP;
	}
	if (!rownames)
	{
		QSlog("EGLPNUM_TYPENAME_ILLlib_rownames called with NULL rownames");
		rval = 1;
		ILL_CLEANUP;
	}

	qslp = lp->O;
	nrows = qslp->nrows;

	if (qslp->rownames == 0)
	{
		QSlog("LP does not have rownames assigned");
		rval = 1;
		ILL_CLEANUP;
	}

	for (i = 0; i < nrows; i++)
	{
		len = strlen (qslp->rownames[i]) + 1;
		ILL_SAFE_MALLOC (rownames[i], len, char);

		strcpy (rownames[i], qslp->rownames[i]);
		rcount++;
	}

CLEANUP:

	if (rval)
	{
		for (i = 0; i < rcount; i++)
		{
			ILL_IFFREE (rownames[i], char);
		}
	}
	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLlib_getintflags (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int *intflags)
{
	int j, nstruct, rval = 0;
	EGLPNUM_TYPENAME_ILLlpdata *qslp;

	if (!lp)
	{
		QSlog("EGLPNUM_TYPENAME_ILLlib_getintflags called without an LP");
		rval = 1;
		ILL_CLEANUP;
	}

	qslp = lp->O;
	nstruct = qslp->nstruct;

	if (qslp->intmarker == 0)
	{
		for (j = 0; j < nstruct; j++)
		{
			intflags[j] = 0;
		}
	}
	else
	{
		for (j = 0; j < nstruct; j++)
		{
			if (qslp->intmarker[j])
			{
				intflags[j] = 1;
			}
			else
			{
				intflags[j] = 0;
			}
		}
	}

CLEANUP:

	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLlib_colnames (
	EGLPNUM_TYPENAME_lpinfo * lp,
	char **colnames)
{
	EGLPNUM_TYPENAME_ILLlpdata *qslp;
	int nstruct, len, i, ccount = 0;
	int rval = 0;

	if (!lp)
	{
		QSlog("EGLPNUM_TYPENAME_ILLlib_colnames called without an LP");
		rval = 1;
		ILL_CLEANUP;
	}
	if (!colnames)
	{
		QSlog("EGLPNUM_TYPENAME_ILLlib_colnames called with NULL colnames");
		rval = 1;
		ILL_CLEANUP;
	}

	qslp = lp->O;
	nstruct = qslp->nstruct;

	if (qslp->colnames == 0)
	{
		QSlog("LP does not have colnames assigned");
		rval = 1;
		ILL_CLEANUP;
	}

	for (i = 0; i < nstruct; i++)
	{
		len = strlen (qslp->colnames[i]) + 1;
		ILL_SAFE_MALLOC (colnames[i], len, char);

		strcpy (colnames[i], qslp->colnames[i]);
		ccount++;
	}


CLEANUP:

	if (rval)
	{
		for (i = 0; i < ccount; i++)
		{
			ILL_IFFREE (colnames[i], char);
		}
	}

	EG_RETURN (rval);
}

static int reset_colindex (
	EGLPNUM_TYPENAME_lpinfo * lp)
{
	int rval = 0;
	int test;
	EGLPNUM_TYPENAME_ILLlpdata *qslp = lp->O;

	test = ILLsymboltab_index_ok (&qslp->coltab);
	if (!test)
	{
		rval = ILLsymboltab_index_reset (&qslp->coltab, qslp->nstruct,
																		 qslp->colnames);
		CHECKRVALG (rval, CLEANUP);
	}

CLEANUP:

	EG_RETURN (rval);
}

static int reset_rowindex (
	EGLPNUM_TYPENAME_lpinfo * lp)
{
	int rval = 0;
	int test;
	EGLPNUM_TYPENAME_ILLlpdata *qslp = lp->O;

	test = ILLsymboltab_index_ok (&qslp->rowtab);
	if (!test)
	{
		rval = ILLsymboltab_index_reset (&qslp->rowtab, qslp->nrows,
																		 qslp->rownames);
		CHECKRVALG (rval, CLEANUP);
	}

CLEANUP:

	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLlib_colindex (
	EGLPNUM_TYPENAME_lpinfo * lp,
	const char *name,
	int *colindex)
{
	int rval = 0;
	EGLPNUM_TYPENAME_ILLlpdata *qslp;

	*colindex = -1;

	if (!lp)
	{
		QSlog("EGLPNUM_TYPENAME_ILLlib_colindex called without an LP");
		rval = 1;
		ILL_CLEANUP;
	}

	qslp = lp->O;

	rval = reset_colindex (lp);
	CHECKRVALG (rval, CLEANUP);

	rval = ILLsymboltab_getindex (&qslp->coltab, name, colindex);
	CHECKRVALG (rval, CLEANUP);

CLEANUP:

	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLlib_rowindex (
	EGLPNUM_TYPENAME_lpinfo * lp,
	const char *name,
	int *rowindex)
{
	int rval = 0;
	EGLPNUM_TYPENAME_ILLlpdata *qslp;

	*rowindex = -1;

	if (!lp)
	{
		QSlog("EGLPNUM_TYPENAME_ILLlib_rowindex called without an LP");
		rval = 1;
		ILL_CLEANUP;
	}

	qslp = lp->O;

	rval = reset_rowindex (lp);
	CHECKRVALG (rval, CLEANUP);

	rval = ILLsymboltab_getindex (&qslp->rowtab, name, rowindex);
	CHECKRVALG (rval, CLEANUP);

CLEANUP:

	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLlib_getbasis (
	EGLPNUM_TYPENAME_lpinfo * lp,
	char *cstat,
	char *rstat)
{
	int rval = 0;
	int i, j, nrows;
	EGLPNUM_TYPENAME_ILLlpdata *qslp;

	if (!lp)
	{
		QSlog("EGLPNUM_TYPENAME_ILLlib_getbasis called without an LP");
		rval = 1;
		ILL_CLEANUP;
	}

	if (lp->basisid == -1)
	{
		QSlog("EGLPNUM_TYPENAME_ILLlib_getbasis called with modifed LP");
		rval = 1;
		ILL_CLEANUP;
	}

	nrows = lp->O->nrows;
	qslp = lp->O;

	for (i = 0; i < qslp->nstruct; i++)
	{
		j = qslp->structmap[i];
		switch (lp->vstat[j])
		{
		case STAT_BASIC:
			cstat[i] = QS_COL_BSTAT_BASIC;
			break;
		case STAT_LOWER:
			cstat[i] = QS_COL_BSTAT_LOWER;
			break;
		case STAT_UPPER:
			cstat[i] = QS_COL_BSTAT_UPPER;
			break;
		case STAT_ZERO:
			cstat[i] = QS_COL_BSTAT_FREE;
			break;
		default:
			QSlog("unknown vstat in EGLPNUM_TYPENAME_ILLlib_getbasis: %d", lp->vstat[j]);
			rval = 1;
			ILL_CLEANUP;
		}
	}

	for (i = 0; i < nrows; i++)
	{
		j = qslp->rowmap[i];
		if (qslp->rangeval && EGLPNUM_TYPENAME_EGlpNumIsNeqqZero (qslp->rangeval[i]))
		{
			switch (lp->vstat[j])
			{
			case STAT_BASIC:
				rstat[i] = QS_ROW_BSTAT_BASIC;
				break;
			case STAT_LOWER:
				rstat[i] = QS_ROW_BSTAT_LOWER;
				break;
			case STAT_UPPER:
				rstat[i] = QS_ROW_BSTAT_UPPER;
				break;
			default:
				QSlog("unknown vstat in EGLPNUM_TYPENAME_ILLlib_getbasis 2");
				rval = 1;
				ILL_CLEANUP;
			}
		}
		else
		{
			switch (lp->vstat[j])
			{
			case STAT_BASIC:
				rstat[i] = QS_ROW_BSTAT_BASIC;
				break;
			case STAT_UPPER:
			case STAT_LOWER:
				rstat[i] = QS_ROW_BSTAT_LOWER;
				break;
			default:
				QSlog("unknown vstat in EGLPNUM_TYPENAME_ILLlib_getbasis 3: %d, %d",
										i, lp->vstat[j]);
				rval = 1;
				ILL_CLEANUP;
			}
		}
	}

CLEANUP:

	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLlib_loadbasis (
	EGLPNUM_TYPENAME_ILLlp_basis * B,
	int nstruct,
	int nrows,
	char *cstat,
	char *rstat)
{
	int i;
	int rval = 0;

	EGLPNUM_TYPENAME_ILLlp_basis_init (B);

	if (!cstat || !rstat)
	{
		rval = 1;
		CHECKRVALG (rval, CLEANUP);
	}

	rval = EGLPNUM_TYPENAME_ILLlp_basis_alloc (B, nstruct, nrows);
	CHECKRVALG (rval, CLEANUP);

	for (i = 0; i < nstruct; i++)
	{
		B->cstat[i] = cstat[i];
	}
	for (i = 0; i < nrows; i++)
	{
		B->rstat[i] = rstat[i];
	}

CLEANUP:

	EG_RETURN (rval);
}

#define READ_BASIS_XL 0
#define READ_BASIS_XU 1
#define READ_BASIS_LL 2
#define READ_BASIS_UL 3

int EGLPNUM_TYPENAME_ILLlib_readbasis (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_ILLlp_basis * B,
	const char *fname)
{
	int rval = 0;
	EGLPNUM_TYPENAME_ILLlpdata *qslp = lp->O;
	int nstruct = qslp->nstruct;
	int nrows = qslp->nrows;
	int i, j, end = 0, sec, havename = 0;
	int rowtype, row, col;
	char *bname = 0;
	EGioFile_t *file_in = 0;
	EGLPNUM_TYPENAME_ILLread_mps_state state;
	EGLPNUM_TYPENAME_qsline_reader *in = NULL;

	EGLPNUM_TYPENAME_ILLlp_basis_init (B);

	ILL_SAFE_MALLOC (B->cstat, qslp->nstruct, char);
	ILL_SAFE_MALLOC (B->rstat, qslp->nrows, char);

	B->nstruct = nstruct;
	B->nrows = nrows;

	for (j = 0; j < nstruct; j++)
	{
		B->cstat[j] = QS_COL_BSTAT_LOWER;
	}
	for (i = 0; i < nrows; i++)
	{
		B->rstat[i] = QS_ROW_BSTAT_BASIC;
	}

	file_in = EGioOpen (fname, "r");
	if (file_in == 0)
	{
		QSlog("unable to open %s for reading", fname);
		rval = 1;
		ILL_CLEANUP;
	}

	in = EGLPNUM_TYPENAME_ILLline_reader_new ((EGLPNUM_TYPENAME_qsread_line_fct) EGioGets, file_in);
	rval = EGLPNUM_TYPENAME_ILLmps_state_init (&state, in, fname);
	CHECKRVALG (rval, CLEANUP);

	while (EGLPNUM_TYPENAME_ILLmps_next_line (&state) == 0)
	{
		if (EGLPNUM_TYPENAME_ILLmps_empty_key (&state))
		{

			/* Get the XL XU LL UL line */

			if (!havename)
			{
				rval = EGLPNUM_TYPENAME_ILLmps_error (&state, "BASIS data before NAME\n");
				ILL_CLEANUP;
			}

			if (!strcmp (state.field, "XL"))
			{
				rowtype = READ_BASIS_XL;
			}
			else if (!strcmp (state.field, "XU"))
			{
				rowtype = READ_BASIS_XU;
			}
			else if (!strcmp (state.field, "LL"))
			{
				rowtype = READ_BASIS_LL;
			}
			else if (!strcmp (state.field, "UL"))
			{
				rowtype = READ_BASIS_UL;
			}
			else
			{
				rval = EGLPNUM_TYPENAME_ILLmps_error (&state, "BASIS \"%s\" is invalid\n", state.field);
				ILL_CLEANUP;
			}

			if (EGLPNUM_TYPENAME_ILLmps_next_field (&state) == 0)
			{

				rval = EGLPNUM_TYPENAME_ILLlib_colindex (lp, (const char *) state.field, &col);
				CHECKRVALG (rval, CLEANUP);
				if (col == -1)
				{
					rval = EGLPNUM_TYPENAME_ILLmps_error (&state, "BASIS col not in LP\n");
					ILL_CLEANUP;
				}

				if (rowtype == READ_BASIS_XL || rowtype == READ_BASIS_XU)
				{
					if (EGLPNUM_TYPENAME_ILLmps_next_field (&state) == 0)
					{
						rval = EGLPNUM_TYPENAME_ILLlib_rowindex (lp, (const char *) state.field, &row);
						CHECKRVALG (rval, CLEANUP);
						if (row == -1)
						{
							rval = EGLPNUM_TYPENAME_ILLmps_error (&state, "BASIS row not in LP\n");
							ILL_CLEANUP;
						}
						if (rowtype == READ_BASIS_XL)
						{
							B->cstat[col] = QS_COL_BSTAT_BASIC;
							B->rstat[row] = QS_ROW_BSTAT_LOWER;

						}
						else
						{
							B->cstat[col] = QS_COL_BSTAT_BASIC;
							B->rstat[row] = QS_ROW_BSTAT_UPPER;
						}
					}
					else
					{
						rval = EGLPNUM_TYPENAME_ILLmps_error (&state, "BASIS line needs row and column\n");
						ILL_CLEANUP;
					}
				}
				else
				{
					if (rowtype == READ_BASIS_LL)
					{
						B->cstat[col] = QS_COL_BSTAT_LOWER;
					}
					else
					{
						B->cstat[col] = QS_COL_BSTAT_UPPER;
					}
				}
			}
			else
			{
				rval = EGLPNUM_TYPENAME_ILLmps_error (&state, "BASIS line has no row/column\n");
				ILL_CLEANUP;
			}
		}
		else
		{
			/* found a section indicator in col 1 */
			if (!strcmp (state.key, EGLPNUM_TYPENAME_ILLmps_section_name[ILL_MPS_ENDATA]))
			{
				end = 1;
				break;									/* done reading */
			}

			sec = ILLutil_index (EGLPNUM_TYPENAME_ILLmps_section_name, state.key);
			if (sec < 0 || sec != ILL_MPS_NAME)
			{
				rval = EGLPNUM_TYPENAME_ILLmps_error (&state, "BASIS \"%s\" is not a key\n", state.key);
				ILL_CLEANUP;
			}

			if (havename)
			{
				rval = EGLPNUM_TYPENAME_ILLmps_error (&state, "BASIS two name sections\n");
				ILL_CLEANUP;
			}

			havename = 1;

			if (EGLPNUM_TYPENAME_ILLmps_empty_field (&state))
			{
				EGLPNUM_TYPENAME_ILLmps_warn (&state, "BASIS blank NAME.");
			}
			else
			{
				ILL_UTIL_STR (bname, state.field);
				QSlog("Basis Name: %s", bname);
				if (strcmp (bname, qslp->probname))
				{
					EGLPNUM_TYPENAME_ILLmps_warn (&state, "BASIS name does not match LP.");
				}
			}
		}
	}

	if (!end)
	{
		EGLPNUM_TYPENAME_ILLmps_warn (&state, "Missing ENDATA in basis file.");
	}
	if (!EGLPNUM_TYPENAME_ILLmps_next_line (&state))
	{
		EGLPNUM_TYPENAME_ILLmps_warn (&state, "Ignoring text after ENDATA.");
	}

	if (!havename)
	{
		rval = EGLPNUM_TYPENAME_ILLmps_error (&state, "BASIS no name section\n");
		ILL_CLEANUP;
	}

	/* Correct the free variables */

	for (j = 0; j < nstruct; j++)
	{
		col = lp->O->structmap[j];
		if (EGLPNUM_TYPENAME_EGlpNumIsEqqual (qslp->lower[col], EGLPNUM_TYPENAME_ILL_MINDOUBLE) &&
				EGLPNUM_TYPENAME_EGlpNumIsEqqual (qslp->upper[col], EGLPNUM_TYPENAME_ILL_MAXDOUBLE) &&
				B->cstat[j] == QS_COL_BSTAT_LOWER)
		{
			B->cstat[j] = QS_COL_BSTAT_FREE;
		}
	}

CLEANUP:

	if (file_in)
		EGioClose (file_in);
	EGLPNUM_TYPENAME_ILLline_reader_free (in);

	if (rval)
	{
		EGLPNUM_TYPENAME_ILLlp_basis_free (B);
	}
	ILL_IFFREE (bname, char);

	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLlib_writebasis (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_ILLlp_basis * B,
	const char *fname)
{
	int rval = 0;
	EGioFile_t *out = 0;
	char *cstat = 0;
	char *rstat = 0;
	EGLPNUM_TYPENAME_ILLlpdata *qslp;
	int i, j, nstruct, nrows;

	/* NOTE: non-basic free variables are encoded as non-basic at lower */

	if (!lp)
	{
		QSlog("EGLPNUM_TYPENAME_ILLlib_writebasis called without an LP");
		rval = 1;
		ILL_CLEANUP;
	}
	if (!B && lp->basisid == -1)
	{
		QSlog("EGLPNUM_TYPENAME_ILLlib_writebasis called with unsolved LP");
		rval = 1;
		ILL_CLEANUP;
	}

	qslp = lp->O;
	nstruct = qslp->nstruct;
	nrows = qslp->nrows;

	out = EGioOpen (fname, "w");
	if (out == 0)
	{
		QSlog("unable to open %s for writing", fname);
		rval = 1;
		ILL_CLEANUP;
	}

	if (B)
	{
		cstat = B->cstat;
		rstat = B->rstat;
	}
	else
	{
		ILL_SAFE_MALLOC (cstat, nstruct, char);
		ILL_SAFE_MALLOC (rstat, nrows, char);

		rval = EGLPNUM_TYPENAME_ILLlib_getbasis (lp, cstat, rstat);
		CHECKRVALG (rval, CLEANUP);
	}

	EGioPrintf (out, "NAME    %s\n", qslp->probname);

	/* Pick out the non-basic rows and find a matching basic column */

	i = 0;
	j = 0;
	do
	{
		while (i < nrows && rstat[i] == QS_ROW_BSTAT_BASIC)
		{
			i++;
		}
		if (i < nrows)
		{
			while (j < nstruct && cstat[j] != QS_COL_BSTAT_BASIC)
			{
				j++;
			}
			if (j == nstruct)
			{
				/* No basic column to match the non-basic row */
				QSlog("No basic column to match non-basic row %d", i);
				rval = 1;
				goto CLEANUP;
			}

			if (rstat[i] == QS_ROW_BSTAT_LOWER)
			{
				EGioPrintf (out, " XL %s %s\n", qslp->colnames[j], qslp->rownames[i]);
			}
			else
			{
				EGioPrintf (out, " XU %s %s\n", qslp->colnames[j], qslp->rownames[i]);
			}
			i++;
			j++;
		}
	} while (i < nrows);

	/* Now go through and output the non-basic cols at upper bound */

	for (j = 0; j < nstruct; j++)
	{
		if (cstat[j] == QS_COL_BSTAT_UPPER)
		{
			EGioPrintf (out, " UL %s\n", qslp->colnames[j]);
		}
	}

	EGioPrintf (out, "ENDATA\n");

CLEANUP:

	if (out)
		EGioClose (out);
	if (!B)
	{
		ILL_IFFREE (cstat, char);
		ILL_IFFREE (rstat, char);
	}
	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLlib_getrownorms (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_price_info * pinf,
	EGLPNUM_TYPE * rownorms)
{
	int rval = 0;
	int i, j, basic = 0;
	EGLPNUM_TYPENAME_ILLlpdata *qslp = lp->O;
	int nstruct = lp->O->nstruct;
	int nrows = lp->O->nrows;

	check_pinf (pinf, &rval);
	if (rval)
	{
/*
        QSlog("dual steepest edge norms not available");
*/
		ILL_CLEANUP;
	}

	for (i = 0; i < nstruct; i++)
	{
		j = qslp->structmap[i];
		if (lp->vstat[j] == STAT_BASIC)
		{
			EGLPNUM_TYPENAME_EGlpNumCopy (rownorms[basic++], pinf->dsinfo.norms[lp->vindex[j]]);
		}
	}
	for (i = 0; i < nrows; i++)
	{
		j = qslp->rowmap[i];
		if (lp->vstat[j] == STAT_BASIC)
		{
			EGLPNUM_TYPENAME_EGlpNumCopy (rownorms[basic++], pinf->dsinfo.norms[lp->vindex[j]]);
		}
	}

	if (basic != nrows)
	{
		QSlog("error in EGLPNUM_TYPENAME_ILLlib_getrownorms");
		rval = 1;
		ILL_CLEANUP;
	}

CLEANUP:

/*
    EG_RETURN(rval);
*/
	return rval;									/* Don't want error message */
}

int EGLPNUM_TYPENAME_ILLlib_loadrownorms (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_price_info * pinf,
	EGLPNUM_TYPE * rownorms)
{
	int rval = 0;

	rval = EGLPNUM_TYPENAME_ILLprice_load_rownorms (lp, rownorms, pinf);
	CHECKRVALG (rval, CLEANUP);

CLEANUP:

	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLlib_recompute_rownorms (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_price_info * pinf)
{
	int rval = 0;

	rval = EGLPNUM_TYPENAME_ILLprice_build_pricing_info (lp, pinf, DUAL_PHASEII);
	CHECKRVALG (rval, CLEANUP);

CLEANUP:

	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLlib_iter (
	EGLPNUM_TYPENAME_lpinfo * lp)
{
	int iter = 0;

	if (lp && lp->cnts)
	{
		iter = lp->cnts->pI_iter + lp->cnts->pII_iter +
			lp->cnts->dI_iter + lp->cnts->dII_iter;
	}

	return iter;
}

//#define PRINT_TOL 0.000001
#define PRINT_TOL EGLPNUM_TYPENAME_PFEAS_TOLER

int EGLPNUM_TYPENAME_ILLlib_print_x (
	EGioFile_t * fd,
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_ILLlp_cache * C,
	EGLPNUM_TYPE * x,
	int nonZerosOnly)
{
	int rval = 0;
	int j;
	EGLPNUM_TYPENAME_ILLlpdata *qslp = lp->O;
	EGLPNUM_TYPE *dx, *myx = 0;
	char *strtmp;

	/* If x is not specified, grab the LP solution */

	if (!x)
	{
		myx = EGLPNUM_TYPENAME_EGlpNumAllocArray (lp->ncols);
		rval = EGLPNUM_TYPENAME_ILLlib_get_x (lp, C, myx);
		CHECKRVALG (rval, CLEANUP);
		dx = myx;
	}
	else
	{
		dx = x;
	}

	EGioPrintf (fd, "Solution Values\n");
	for (j = 0; j < qslp->nstruct; j++)
	{
		/*if (!nonZerosOnly || dx[j] > PRINT_TOL || dx[j] < -PRINT_TOL) */
		if (!nonZerosOnly || EGLPNUM_TYPENAME_EGlpNumIsNeqZero (dx[j], PRINT_TOL))
		{
			strtmp = EGLPNUM_TYPENAME_EGlpNumGetStr (dx[j]);
			ILL_FAILfalse (qslp->colnames[j] != NULL, "no NULL names PLEASE!");
			EGioPrintf (fd, "%s = %s\n", qslp->colnames[j], strtmp);
			EGioFlush (fd);
			EGfree (strtmp);
		}
	}

CLEANUP:
	EGLPNUM_TYPENAME_EGlpNumFreeArray (myx);
	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLlib_findName (
	EGLPNUM_TYPENAME_ILLlpdata * qslp,
	int forRow,
	const char *name,
	int id,
	char buf[ILL_namebufsize])
{
	ILLsymboltab *tab;
	const char *mode;
	const char *p1, *p2;
	int sind, rval = 0;

	id++;
	tab = (forRow) ? &qslp->rowtab : &qslp->coltab;
	if (tab->tablesize == 0)
		ILLsymboltab_create (tab, 100);
	p1 = (forRow) ? "c" : "x";
	p2 = (forRow) ? "c_" : "x_";
	mode = (forRow) ? "row" : "column";
	if (name == 0)
	{
		ILLsymboltab_unique_name (tab, id, p1, buf);
		/*
		 * QSlog("Generating %s name \"%s\".", mode, buf);
		 */
	}
	else
	{
		strcpy (buf, name);
	}
	if (!ILLsymboltab_lookup (tab, buf, &sind))
	{
		rval = ILLsymboltab_uname (&qslp->rowtab, buf, p1, p2);
		if (name != NULL)
		{
			QSlog("Changing %s name \"%s\" to \"%s\".", mode, name, buf);
		}
		CHECKRVALG (rval, CLEANUP);
	}
CLEANUP:
	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLwrite_lp_file (
	EGLPNUM_TYPENAME_ILLlpdata * lp,
	EGioFile_t * out,
	EGLPNUM_TYPENAME_qserror_collector * c)
{
	int rval = 0;
	qsstring_reporter rep;

	ILLstring_reporter_copy (&rep, &lp->reporter);
	ILLstring_reporter_init (&lp->reporter, (qsreport_string_fct) EGioWrite, out);
	rval = EGLPNUM_TYPENAME_ILLwrite_lp (lp, c);
	ILLstring_reporter_copy (&lp->reporter, &rep);
	return rval;
}

static void check_pinf (
	EGLPNUM_TYPENAME_price_info * pinf,
	int *it_exists)
{
	if (!pinf || pinf->dI_price != QS_PRICE_DSTEEP ||
			pinf->dII_price != QS_PRICE_DSTEEP || pinf->dsinfo.norms == 0)
	{
		*it_exists = 1;
	}
	else
	{
		*it_exists = 0;
	}
}

#if 0
static int test_matrix (
	EGLPNUM_TYPENAME_ILLmatrix * A)
{
	int rval = 0;
	int i, j, k;
	int ncols = A->matcols;
	int nrows = A->matrows;
	int matsize = A->matsize;
	int *mbeg = A->matbeg;
	int *mcnt = A->matcnt;
	int *mind = A->matind;
	int *tempi = 0;


	if (matsize == 0)
		ILL_CLEANUP;

	ILL_SAFE_MALLOC (tempi, matsize, int);

	for (i = 0; i < matsize; i++)
		tempi[i] = 0;

	for (i = 0; i < ncols; i++)
	{
		for (j = 0; j < mcnt[i]; j++)
		{
			k = mind[mbeg[i] + j];
			if (k < 0 || k >= nrows)
			{
				QSlog("ERROR IN MATRIX: %d", k);
				QSlog("ncols = %d, bad col = %d", ncols, i);
				QSlog("bad cnt = %d, bad index = %d", mcnt[i], mbeg[i] + j);
				QSlog("matcolsize = %d, matsize = %d", A->matcolsize, A->matsize);
				rval = 1;
				ILL_CLEANUP;
			}
			if (tempi[mbeg[i] + j] != 0)
			{
				QSlog("ERROR: over written matrix");
				QSlog("ncols = %d, bad col = %d", ncols, i);
				QSlog("nrows = %d", nrows);
				QSlog("bad cnt = %d, bad index = %d", mcnt[i], mbeg[i] + j);
				rval = 1;
				ILL_CLEANUP;
			}
			else
			{
				tempi[mbeg[i] + j] = 1;
			}
		}
	}

	for (i = A->matsize - A->matfree; i < A->matsize; i++)
	{
		if (tempi[i] != 0)
		{
			QSlog("ERROR: free space is being used");
			rval = 1;
			ILL_CLEANUP;
		}
	}

CLEANUP:

	ILL_IFFREE (tempi, int);

	EG_RETURN (rval);
}
#endif
