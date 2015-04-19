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

/* RCS_INFO = "$RCSfile: qsopt.c,v $ $Revision: 1.2 $ $Date: 2003/11/05 16:49:52 $"; */
static int TRACE = 0;

/****************************************************************************/
/*                                                                          */
/*                     User-level Functions                                 */
/*                                                                          */
/*  EXPORTED FUNCTIONS                                                      */
/*                                                                          */
/*    int EGLPNUM_TYPENAME_QSopt_primal (EGLPNUM_TYPENAME_QSdata *p, int *status)                             */
/*    int EGLPNUM_TYPENAME_QSopt_dual (EGLPNUM_TYPENAME_QSdata *p, int *status)                               */
/*    EGLPNUM_TYPENAME_QSdata *EGLPNUM_TYPENAME_QScreate_prob (const char *name, int objsense)                */
/*    EGLPNUM_TYPENAME_QSdata *EGLPNUM_TYPENAME_QSread_prob (const char *filename, const char *filetype)      */
/*    EGLPNUM_TYPENAME_QSdata *EGLPNUM_TYPENAME_QSload_prob (const char *probname, int ncols, int nrows,      */
/*        int *cmatcnt, int *cmatbeg, int *cmatind, double *cmatval,        */
/*        int objsense, double *obj, double *rhs, char *sense,              */
/*        double *lower, double *upper, const char **colnames,              */
/*        const char **rownames)                                            */
/*    EGLPNUM_TYPENAME_QSdata *EGLPNUM_TYPENAME_QScopy_prob (EGLPNUM_TYPENAME_QSdata *p, const char *newname)                  */
/*    int EGLPNUM_TYPENAME_QSchange_objsense (EGLPNUM_TYPENAME_QSdata *p, int newsense)                       */
/*    int EGLPNUM_TYPENAME_QSget_objsense (EGLPNUM_TYPENAME_QSdata *p, int *objsense)                         */
/*    int EGLPNUM_TYPENAME_QSnew_col (EGLPNUM_TYPENAME_QSdata *p, double obj, double lower, double upper,     */
/*        const char *name)                                                 */
/*    int EGLPNUM_TYPENAME_QSadd_cols (EGLPNUM_TYPENAME_QSdata *p, int num, int *cmatcnt, int *cmatbeg,       */
/*        int *cmatind, double *cmatval, double *obj, double *lower,        */
/*        double *upper, const char **names)                                */
/*    int EGLPNUM_TYPENAME_QSadd_col (EGLPNUM_TYPENAME_QSdata *p, int cnt, int *cmatind, double *cmatval,     */
/*        double obj, double lower, double upper, const char *name)         */
/*    int EGLPNUM_TYPENAME_QSnew_row (EGLPNUM_TYPENAME_QSdata *p, double rhs, const char sense, char *name)   */
/*    int EGLPNUM_TYPENAME_QSadd_rows (EGLPNUM_TYPENAME_QSdata *p, int num, int *rmatcnt, int *rmatbeg,       */
/*        int *rmatind, double *rmatval, double *rhs, char *sense,          */
/*        char **names)                                                     */
/*    int EGLPNUM_TYPENAME_QSadd_row (EGLPNUM_TYPENAME_QSdata *p, int cnt, int *rmatind, double *rmatval,     */
/*        double rhs, char sense, const char *name)                         */
/*    int EGLPNUM_TYPENAME_QSdelete_rows (EGLPNUM_TYPENAME_QSdata *p, int num, int *dellist)                  */
/*    int EGLPNUM_TYPENAME_QSdelete_row (EGLPNUM_TYPENAME_QSdata *p, int rowindex)                            */
/*    int EGLPNUM_TYPENAME_QSdelete_setrows (EGLPNUM_TYPENAME_QSdata *p, int *flags)                          */
/*    int EGLPNUM_TYPENAME_QSdelete_cols (EGLPNUM_TYPENAME_QSdata *p, int num, int *dellist)                  */
/*    int EGLPNUM_TYPENAME_QSdelete_col (EGLPNUM_TYPENAME_QSdata *p, int colindex)                            */
/*    int EGLPNUM_TYPENAME_QSdelete_setcols (EGLPNUM_TYPENAME_QSdata *p, int *flags)                          */
/*    int EGLPNUM_TYPENAME_QSdelete_named_column (EGLPNUM_TYPENAME_QSdata *p, const char *colname)            */
/*    int EGLPNUM_TYPENAME_QSdelete_named_columns_list (EGLPNUM_TYPENAME_QSdata *p, int num,                  */
/*        const char **colnames)                                            */
/*    int EGLPNUM_TYPENAME_QSdelete_named_row (EGLPNUM_TYPENAME_QSdata *p, const char *rowname)               */
/*    int EGLPNUM_TYPENAME_QSdelete_named_rows_list (EGLPNUM_TYPENAME_QSdata *p, int num,                     */
/*        const char **rownames)                                            */
/*    int EGLPNUM_TYPENAME_QSchange_senses (EGLPNUM_TYPENAME_QSdata *p, int num, int *rowlist, char *sense)   */
/*    int EGLPNUM_TYPENAME_QSchange_sense (EGLPNUM_TYPENAME_QSdata *p, int rowindex, char sense)              */
/*    int EGLPNUM_TYPENAME_QSchange_coef (EGLPNUM_TYPENAME_QSdata *p, int rowindex, int colindex,             */
/*        double coef)                                                      */
/*    int EGLPNUM_TYPENAME_QSchange_objcoef (EGLPNUM_TYPENAME_QSdata *p, int indx, double coef)               */
/*    int EGLPNUM_TYPENAME_QSchange_rhscoef (EGLPNUM_TYPENAME_QSdata *p, int indx, double coef)               */
/*    int EGLPNUM_TYPENAME_QSchange_bounds (EGLPNUM_TYPENAME_QSdata *p, int num, int *collist, char *lu,      */
/*        double *bounds)                                                   */
/*    int EGLPNUM_TYPENAME_QSchange_bound (EGLPNUM_TYPENAME_QSdata *p, int indx, char lu, double bound)       */
/*    int EGLPNUM_TYPENAME_QSwrite_basis (EGLPNUM_TYPENAME_QSdata *p, QSbasis *B, const char *filename)       */
/*    QSbasis *EGLPNUM_TYPENAME_QSget_basis (EGLPNUM_TYPENAME_QSdata *p)                                      */
/*    QSbasis *EGLPNUM_TYPENAME_QSread_basis (EGLPNUM_TYPENAME_QSdata *p, const char *filename)               */
/*    int EGLPNUM_TYPENAME_QSload_basis (EGLPNUM_TYPENAME_QSdata *p, QSbasis *B)                              */
/*    int EGLPNUM_TYPENAME_QSread_and_load_basis (EGLPNUM_TYPENAME_QSdata *p, const char *filename)           */
/*    int EGLPNUM_TYPENAME_QSload_basis_array (EGLPNUM_TYPENAME_QSdata *p, char *cstat, char *rstat)          */
/*    int EGLPNUM_TYPENAME_QSload_basis_and_row_norms_array (EGLPNUM_TYPENAME_QSdata *p, char *cstat,         */
/*      char *rstat, double *rownorms)                                      */
/*    int EGLPNUM_TYPENAME_QSget_basis_array (EGLPNUM_TYPENAME_QSdata *p, char *cstat, char *rstat)           */
/*    int EGLPNUM_TYPENAME_QSget_basis_and_row_norms_array (EGLPNUM_TYPENAME_QSdata *p, char *cstat,          */
/*      char *rstat, double *rownorms)                                      */
/*    int EGLPNUM_TYPENAME_QSget_binv_row (EGLPNUM_TYPENAME_QSdata *p, int indx, double *binvrow)             */
/*    int EGLPNUM_TYPENAME_QSget_tableau_row (EGLPNUM_TYPENAME_QSdata *p, int indx, double *tableaurow)       */
/*    int EGLPNUM_TYPENAME_QSget_basis_order (EGLPNUM_TYPENAME_QSdata *p, int *basorder)                      */
/*    int EGLPNUM_TYPENAME_QSget_status (EGLPNUM_TYPENAME_QSdata *p, int *status)                             */
/*    int EGLPNUM_TYPENAME_QSget_solution (EGLPNUM_TYPENAME_QSdata *p, double *value, double *x,              */
/*        double *pi, double *slack, double *rc),                           */
/*    int EGLPNUM_TYPENAME_QSget_objval (EGLPNUM_TYPENAME_QSdata *p, double *value)                           */
/*    int EGLPNUM_TYPENAME_QSget_x_array (EGLPNUM_TYPENAME_QSdata *p, double *x)                              */
/*    int EGLPNUM_TYPENAME_QSget_rc_array (EGLPNUM_TYPENAME_QSdata *p, double *rc)                            */
/*    int EGLPNUM_TYPENAME_QSget_pi_array (EGLPNUM_TYPENAME_QSdata *p, double *pi)                            */
/*    int EGLPNUM_TYPENAME_QSget_slack_array (EGLPNUM_TYPENAME_QSdata *p, double *slack)                      */
/*    int EGLPNUM_TYPENAME_QSget_infeas_array (EGLPNUM_TYPENAME_QSdata *p, double *pi)                        */
/*    int EGLPNUM_TYPENAME_QSget_named_x (EGLPNUM_TYPENAME_QSdata *p, const char *colname, double *val)       */
/*    int EGLPNUM_TYPENAME_QSget_named_rc (EGLPNUM_TYPENAME_QSdata *p, const char *colname, double *val)      */
/*    int EGLPNUM_TYPENAME_QSget_named_pi (EGLPNUM_TYPENAME_QSdata *p, const char *rowname, double *val)      */
/*    int EGLPNUM_TYPENAME_QSget_named_slack (EGLPNUM_TYPENAME_QSdata *p, const char *rowname, double *val)   */
/*    int EGLPNUM_TYPENAME_QSget_colcount (EGLPNUM_TYPENAME_QSdata *p)                                        */
/*    int EGLPNUM_TYPENAME_QSget_rowcount (EGLPNUM_TYPENAME_QSdata *p)                                        */
/*    int EGLPNUM_TYPENAME_QSget_nzcount (EGLPNUM_TYPENAME_QSdata *p)                                         */
/*    int EGLPNUM_TYPENAME_QSget_obj (EGLPNUM_TYPENAME_QSdata *p, double *obj),                               */
/*    int EGLPNUM_TYPENAME_QSget_rhs (EGLPNUM_TYPENAME_QSdata *p, double *rhs)                                */
/*    char* EGLPNUM_TYPENAME_QSget_probname (EGLPNUM_TYPENAME_QSdata *p)                                      */
/*    char* EGLPNUM_TYPENAME_QSget_objname (EGLPNUM_TYPENAME_QSdata *p)                                       */
/*    int EGLPNUM_TYPENAME_QSget_columns (EGLPNUM_TYPENAME_QSdata *p, int **colcnt, int **colbeg,             */
/*        int **colind, double **colval, double **obj, double **lower,      */
/*        double **upper, char ***names)                                    */
/*    int EGLPNUM_TYPENAME_QSget_columns_list (EGLPNUM_TYPENAME_QSdata *p, int num, int *collist,             */
/*        int **colcnt, int **colbeg, int **colind, double **colval,        */
/*        double **obj, double **lower, double **upper, char ***names)      */
/*    int EGLPNUM_TYPENAME_QSget_rows (EGLPNUM_TYPENAME_QSdata *p, int **rowcnt, int **rowbeg, int **rowind,  */
/*        double **rowval, double **rhs, char **sense, char ***names)       */
/*    int EGLPNUM_TYPENAME_QSget_rows_list (EGLPNUM_TYPENAME_QSdata *p, int num, int *rowlist, int **rowcnt,  */
/*        int **rowbeg, int **rowind, double **rowval, double **rhs,        */
/*        char **sense, char ***names)                                      */
/*    int EGLPNUM_TYPENAME_QSget_column_index (EGLPNUM_TYPENAME_QSdata *p, const char *name, int *colindex)   */
/*    int EGLPNUM_TYPENAME_QSget_row_index (EGLPNUM_TYPENAME_QSdata *p, const char *name, int *rowindex)      */
/*    int EGLPNUM_TYPENAME_QSget_rownames (EGLPNUM_TYPENAME_QSdata *p, char **rownames)                       */
/*    int EGLPNUM_TYPENAME_QSget_colnames (EGLPNUM_TYPENAME_QSdata *p, char **colnames)                       */
/*    int EGLPNUM_TYPENAME_QSget_bound (EGLPNUM_TYPENAME_QSdata *p, int colindex, char lu, double *bound)     */
/*    int EGLPNUM_TYPENAME_QSget_bounds (EGLPNUM_TYPENAME_QSdata *p, double *lower, double *upper)            */
/*    int EGLPNUM_TYPENAME_QSget_intcount (EGLPNUM_TYPENAME_QSdata *p, int *count)                            */
/*    int EGLPNUM_TYPENAME_QSget_intflags (EGLPNUM_TYPENAME_QSdata *p, int *intflags)                         */
/*    int EGLPNUM_TYPENAME_QScompute_row_norms (EGLPNUM_TYPENAME_QSdata *p)                                   */
/*    void EGLPNUM_TYPENAME_QSfree_prob (EGLPNUM_TYPENAME_QSdata *p)                                          */
/*    void EGLPNUM_TYPENAME_QSfree_basis (QSbasis *B)                                        */
/*    int EGLPNUM_TYPENAME_QSwrite_prob (EGLPNUM_TYPENAME_QSdata *p, const char *filename,                    */
/*        const char *filetype)                                             */
/*    int EGLPNUM_TYPENAME_QSwrite_prob_file (EGLPNUM_TYPENAME_QSdata *p, FILE *file, const char *filetype)   */
/*    int EGLPNUM_TYPENAME_QSset_param (EGLPNUM_TYPENAME_QSdata *p, int whichparam, int newvalue)             */
/*    int QSset_param_double (EGLPNUM_TYPENAME_QSdata *p, int whichparam, double newvalue)   */
/*    int EGLPNUM_TYPENAME_QSget_param (EGLPNUM_TYPENAME_QSdata *p, int whichparam, int *value)               */
/*    int QSget_param_double (EGLPNUM_TYPENAME_QSdata *p, int whichparam, double *value)     */
/*    int EGLPNUM_TYPENAME_QStest_row_norms (EGLPNUM_TYPENAME_QSdata *p)                                      */
/*    int EGLPNUM_TYPENAME_QSopt_strongbranch (EGLPNUM_TYPENAME_QSdata *p, int ncand, int *candidatelist,     */
/*        double *xlist, double *down_vals, double *up_vals,                */
/*        int iterations, double objbound)                                  */
/*    int EGLPNUM_TYPENAME_QSopt_pivotin_row (EGLPNUM_TYPENAME_QSdata *p, int rcnt, int *rlist)               */
/*    int EGLPNUM_TYPENAME_QSopt_pivotin_col (EGLPNUM_TYPENAME_QSdata *p, int ccnt, int *clist)               */
/*    void EGLPNUM_TYPENAME_QSfree (void *ptr)                                               */
/*    void EGLPNUM_TYPENAME_QSstart (void)                                                   */
/*    void EGLPNUM_TYPENAME_QSend (void)                                                     */
/*    char *EGLPNUM_TYPENAME_QSversion (void))                                               */
/*                                                                          */
/*    NEW FUNCTIONS - Add to Docs                                           */
/*                                                                          */
/*    char *EGLPNUM_TYPENAME_QSversion (void))                                               */
/*    int EGLPNUM_TYPENAME_QSget_objsense (EGLPNUM_TYPENAME_QSdata *p, int *objsense)                         */
/*                                                                          */
/****************************************************************************/

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <stdlib.h>
#include <string.h>

#include "qsopt_EGLPNUM_TYPENAME.h"

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
#include "qstruct_EGLPNUM_TYPENAME.h"
#include "lib_EGLPNUM_TYPENAME.h"
#include "mps_EGLPNUM_TYPENAME.h"
#include "lp_EGLPNUM_TYPENAME.h"


void EGLPNUM_TYPENAME_QSset_precision (
	const unsigned prec)
{
	EGlpNumSetPrecision (prec);
	EGLPNUM_TYPENAME_ILLchange_precision ();
	/* change the numbers */
}

static void init_basis (
	QSbasis * B),
  free_cache (
	EGLPNUM_TYPENAME_QSdata * p);

static int opt_work ( EGLPNUM_TYPENAME_QSdata * p, int *status, int primal_or_dual),
  qsbasis_to_illbasis ( QSbasis * qB, EGLPNUM_TYPENAME_ILLlp_basis * B),
  illbasis_to_qsbasis ( EGLPNUM_TYPENAME_ILLlp_basis * B, QSbasis * qB),
  grab_basis ( EGLPNUM_TYPENAME_QSdata * p),
  check_qsdata_pointer ( EGLPNUM_TYPENAME_QSdata * p);


EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSopt_primal (
	EGLPNUM_TYPENAME_QSdata * p,
	int *status)
{
	int rval = 0;

	if (status)
		*status = QS_LP_UNSOLVED;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	/* If both the basis and the cache exist, then skip the optimization */

	if (!p->basis || !p->cache)
	{
		rval = opt_work (p, status, 0);
		CHECKRVALG (rval, CLEANUP);
	}
	else
	{
		if (status)
			*status = p->cache->status;
	}

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSopt_dual (
	EGLPNUM_TYPENAME_QSdata * p,
	int *status)
{
	int rval = 0;

	if (status)
		*status = QS_LP_UNSOLVED;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	if (!p->basis || !p->cache || !p->factorok)
	{
		rval = opt_work (p, status, 1);
		CHECKRVALG (rval, CLEANUP);
	}
	else
	{
		if (status)
			*status = p->cache->status;
	}

CLEANUP:

	if (rval == QS_LP_CHANGE_PREC)
	{
		MESSAGE (__QS_SB_VERB, "Changing precision");
		return rval;
	}
	EG_RETURN (rval);
}

static int opt_work (
	EGLPNUM_TYPENAME_QSdata * p,
	int *status,
	int primal_or_dual)
{
	int rval = 0;
	int rstatus = QS_LP_UNSOLVED;
	EGLPNUM_TYPENAME_QSdata *p2 = 0;

	if (p->basis)
	{
		if (p->basis->nstruct != p->qslp->nstruct ||
				p->basis->nrows != p->qslp->nrows)
		{
			QSlog("Size of basis does not match LP");
			rval = 1;
			goto CLEANUP;
		}
	}

	if (!p->basis && p->lp->basisid == -1 && p->simplex_scaling == 1)
	{
		/* Try scaling by copying the LP and solving */

		EGLPNUM_TYPENAME_ILLprice_free_pricing_info (p->pricing);	/* Just to be sure  */
		p->factorok = 0;						/* that p is clean. */

		p2 = EGLPNUM_TYPENAME_QScopy_prob (p, "scaled_lp");
		if (p2 == 0)
			goto CLEANUP;

		rval = EGLPNUM_TYPENAME_ILLlp_scale (p2->qslp);
		CHECKRVALG (rval, CLEANUP);

		if (primal_or_dual == 0)
		{
			rval = EGLPNUM_TYPENAME_ILLlib_optimize (p2->lp, p2->basis, p2->pricing,
															PRIMAL_SIMPLEX, 0, p2->simplex_display, 
															&(p->itcnt));
		}
		else
		{
			rval = EGLPNUM_TYPENAME_ILLlib_optimize (p2->lp, p2->basis, p2->pricing,
															DUAL_SIMPLEX, 0, p2->simplex_display,
															&(p->itcnt));
		}
		CHECKRVALG (rval, CLEANUP);

		rval = grab_basis (p2);
		CHECKRVALG (rval, CLEANUP);

		if (p->basis)
		{
			EGLPNUM_TYPENAME_ILLlp_basis_free (p->basis);
			ILL_IFFREE (p->basis, EGLPNUM_TYPENAME_ILLlp_basis);
		}
		p->basis = p2->basis;
		p2->basis = 0;
		EGLPNUM_TYPENAME_QSfree_prob (p2);
		p2 = 0;
	}

	if (primal_or_dual == 0)
	{
		if (p->factorok == 0)
		{
			if (p->basis == 0)
				p->lp->basisid = -1;
			rval = EGLPNUM_TYPENAME_ILLlib_optimize (p->lp, p->basis, p->pricing, PRIMAL_SIMPLEX,
															&rstatus, p->simplex_display, &(p->itcnt));
		}
		else
		{
			EGLPNUM_TYPENAME_ILLprice_free_pricing_info (p->pricing);
			if (p->lp->basisid != -1)
				p->lp->fbasisid = p->lp->basisid;
			rval = EGLPNUM_TYPENAME_ILLlib_optimize (p->lp, 0, p->pricing,
															PRIMAL_SIMPLEX, &rstatus, p->simplex_display, 
															&(p->itcnt));
		}
	}
	else
	{
		if (p->factorok == 0)
		{
			if (p->basis == 0)
				p->lp->basisid = -1;
			rval = EGLPNUM_TYPENAME_ILLlib_optimize (p->lp, p->basis, p->pricing, DUAL_SIMPLEX,
															&rstatus, p->simplex_display, &(p->itcnt));
		}
		else
		{
			/* The factorization and rownorms should be up-to-date */
			if (p->lp->basisid != -1)
			{
				p->lp->fbasisid = p->lp->basisid;
			}
			else
			{
				EGLPNUM_TYPENAME_ILLprice_free_pricing_info (p->pricing);
			}
			rval = EGLPNUM_TYPENAME_ILLlib_optimize (p->lp, 0, p->pricing,
															DUAL_SIMPLEX, &rstatus, p->simplex_display, 
															&(p->itcnt));
		}
	}
	CHECKRVALG (rval, CLEANUP);

	rval = grab_basis (p);
	CHECKRVALG (rval, CLEANUP);

	if (rstatus == QS_LP_OPTIMAL)
	{
		rval = EGLPNUM_TYPENAME_QSgrab_cache (p, rstatus);
		CHECKRVALG (rval, CLEANUP);
	}
	else
	{
		free_cache (p);
	}

	p->factorok = 1;

#if 0
	p->lp->basisid = -1;					/* This will cause the basis to be reloaded at the */
	/* next optimization - it could be moved into the  */
	/* add/del routines if we want to cache the        */
	/* factored basis.                                 */
	-switched to having qs_simplex load a basis whenever it is passed;
	the trouble with keeping basisid == -1 is that the QS - level routines
		will not know if the
		current lp has been optimized (so, for example,
																	 getbasis will not work)
		.
#endif
	CLEANUP:p->qstatus = rstatus;

	if (status)
		*status = rstatus;

	if (p2)
		EGLPNUM_TYPENAME_QSfree_prob (p2);
	if (rval == QS_LP_CHANGE_PREC)
	{
		MESSAGE (__QS_SB_VERB, "Changing precision");
		return rval;
	}
	MESSAGE (rval ? 0 : 1000, "Error code %d", rval);
	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSopt_pivotin_row (
	EGLPNUM_TYPENAME_QSdata * p,
	int rcnt,
	int *rlist)
{
	int basismod = 0;
	int rval = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	if (p->pricing == 0)
	{
		ILL_ERROR (rval, "pricing info not available in EGLPNUM_TYPENAME_QSopt_pivotin_row\n");
	}

	rval = EGLPNUM_TYPENAME_ILLsimplex_pivotin (p->lp, p->pricing, rcnt, rlist,
														 SIMPLEX_PIVOTINROW, &basismod);
	CHECKRVALG (rval, CLEANUP);

	rval = grab_basis (p);
	CHECKRVALG (rval, CLEANUP);

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSopt_pivotin_col (
	EGLPNUM_TYPENAME_QSdata * p,
	int ccnt,
	int *clist)
{
	int basismod = 0;
	int rval = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	if (p->pricing == 0)
	{
		ILL_ERROR (rval, "pricing info not available in QSopt_pivotin\n");
	}

	rval = EGLPNUM_TYPENAME_ILLsimplex_pivotin (p->lp, p->pricing, ccnt, clist,
														 SIMPLEX_PIVOTINCOL, &basismod);
	CHECKRVALG (rval, CLEANUP);

	rval = grab_basis (p);
	CHECKRVALG (rval, CLEANUP);

CLEANUP:

	EG_RETURN (rval);
}


EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSopt_strongbranch (
	EGLPNUM_TYPENAME_QSdata * p,
	int ncand,
	int *candidatelist,
	EGLPNUM_TYPE * xlist,
	EGLPNUM_TYPE * down_vals,
	EGLPNUM_TYPE * up_vals,
	int iterations,
	EGLPNUM_TYPE objbound)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	if (p->pricing == 0)
	{
		rval = 1;
		CHECKRVALG (rval, CLEANUP);
	}

	rval = EGLPNUM_TYPENAME_ILLlib_strongbranch (p->lp, p->pricing, candidatelist, ncand,
															xlist, down_vals, up_vals, iterations, objbound,
															&(p->itcnt));
	CHECKRVALG (rval, CLEANUP);

	p->factorok = 0;
	free_cache (p);
	p->qstatus = QS_LP_UNSOLVED;	/* Was set to MODIFIED in free_cache () */

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE EGLPNUM_TYPENAME_QSdata *EGLPNUM_TYPENAME_QScreate_prob (
	const char *name,
	int objsense)
{
	int rval = 0;
	EGLPNUM_TYPENAME_QSdata *p = 0;
	int len;

	ILL_SAFE_MALLOC (p, 1, EGLPNUM_TYPENAME_QSdata);
	if (!p)
	{
		QSlog("out of memory in EGLPNUM_TYPENAME_QScreate_prob");
		rval = 1;
		goto CLEANUP;
	}

	p->qslp = 0;
	p->lp = 0;
	p->pricing = 0;
	p->basis = 0;
	p->cache = 0;
	p->qstatus = QS_LP_UNSOLVED;
	p->factorok = 0;

	p->itcnt.pI_iter = 0;
	p->itcnt.pII_iter = 0;
	p->itcnt.dI_iter = 0;
	p->itcnt.dII_iter = 0;
	p->itcnt.tot_iter = 0;
	EGLPNUM_TYPENAME_EGlpNumInitVar(p->uobjlim);
	EGLPNUM_TYPENAME_EGlpNumInitVar(p->lobjlim);
	EGLPNUM_TYPENAME_EGlpNumCopy(p->uobjlim, EGLPNUM_TYPENAME_ILL_MAXDOUBLE);
	EGLPNUM_TYPENAME_EGlpNumCopy(p->lobjlim, EGLPNUM_TYPENAME_ILL_MINDOUBLE);

	p->simplex_display = 0;
	p->simplex_scaling = 1;

	ILL_SAFE_MALLOC (p->qslp, 1, EGLPNUM_TYPENAME_ILLlpdata);
	if (!p->qslp)
	{
		QSlog("out of memory in EGLPNUM_TYPENAME_QScreate_prob");
		rval = 1;
		goto CLEANUP;
	}
	EGLPNUM_TYPENAME_ILLlpdata_init (p->qslp);

	ILL_SAFE_MALLOC (p->lp, 1, EGLPNUM_TYPENAME_lpinfo);
	if (!p->lp)
	{
		QSlog("out of memory in EGLPNUM_TYPENAME_QScreate_prob");
		rval = 1;
		goto CLEANUP;
	}
	EGLPNUM_TYPENAME_EGlpNumInitVar (p->lp->objval);
	EGLPNUM_TYPENAME_EGlpNumInitVar (p->lp->pobjval);
	EGLPNUM_TYPENAME_EGlpNumInitVar (p->lp->dobjval);
	EGLPNUM_TYPENAME_EGlpNumInitVar (p->lp->pinfeas);
	EGLPNUM_TYPENAME_EGlpNumInitVar (p->lp->dinfeas);
	EGLPNUM_TYPENAME_EGlpNumInitVar (p->lp->objbound);
	EGLPNUM_TYPENAME_EGlpNumInitVar (p->lp->upd.piv);
	EGLPNUM_TYPENAME_EGlpNumInitVar (p->lp->upd.dty);
	EGLPNUM_TYPENAME_EGlpNumInitVar (p->lp->upd.c_obj);
	EGLPNUM_TYPENAME_EGlpNumInitVar (p->lp->upd.tz);
	EGLPNUM_TYPENAME_ILLsimplex_init_lpinfo (p->lp);
	EGLPNUM_TYPENAME_ILLsimplex_load_lpinfo (p->qslp, p->lp);

	ILL_SAFE_MALLOC (p->pricing, 1, EGLPNUM_TYPENAME_price_info);
	if (!p->pricing)
	{
		QSlog("out of memory in EGLPNUM_TYPENAME_QScreate_prob");
		rval = 1;
		goto CLEANUP;
	}
	EGLPNUM_TYPENAME_EGlpNumInitVar (p->pricing->htrigger);
	EGLPNUM_TYPENAME_ILLprice_init_pricing_info (p->pricing);
	p->pricing->pI_price = QS_DEFAULT_PRICE_PI;
	p->pricing->pII_price = QS_DEFAULT_PRICE_PII;
	p->pricing->dI_price = QS_DEFAULT_PRICE_DI;
	p->pricing->dII_price = QS_DEFAULT_PRICE_DII;

	if (name)
	{
		len = strlen (name) + 1;
		ILL_SAFE_MALLOC (p->name, len, char);

		strcpy (p->name, name);
	}
	else
	{
		ILL_SAFE_MALLOC (p->name, 7, char);

		sprintf (p->name, "noname");
	}

	len = strlen (p->name) + 1;
	ILL_SAFE_MALLOC (p->qslp->probname, len, char);

	strcpy (p->qslp->probname, p->name);

	if (objsense == QS_MAX)
	{
		p->qslp->objsense = QS_MAX;
	}

CLEANUP:

	if (rval)
	{
		EGLPNUM_TYPENAME_QSfree_prob (p);
		p = 0;
	}

	return p;
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE EGLPNUM_TYPENAME_QSdata *EGLPNUM_TYPENAME_QSread_prob (
	const char *filename,
	const char *filetype)
{
	EGLPNUM_TYPENAME_QSdata *p = 0;
	EGioFile_t *file = 0;
	EGLPNUM_TYPENAME_QSline_reader reader;

	if ((file = EGioOpen (filename, "r")) == 0)
	{
		perror (filename);
		QSlog("Unable to open \"%s\" for input.", filename);
	}
	if (file == NULL)
		goto CLEANUP;

	reader = EGLPNUM_TYPENAME_ILLline_reader_new ((EGLPNUM_TYPENAME_qsread_line_fct) EGioGets, file);
	p = EGLPNUM_TYPENAME_QSget_prob (reader, filename, filetype);
	EGLPNUM_TYPENAME_QSline_reader_free (reader);	/* Bico - 040723 */

CLEANUP:
	if (file != NULL)
	{
		EGioClose (file);
	}
	return p;
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE EGLPNUM_TYPENAME_QSdata *EGLPNUM_TYPENAME_QSload_prob (
	const char *probname,
	int ncols,
	int nrows,
	int *cmatcnt,
	int *cmatbeg,
	int *cmatind,
	EGLPNUM_TYPE * cmatval,
	int objsense,
	EGLPNUM_TYPE * obj,
	EGLPNUM_TYPE * rhs,
	char *sense,
	EGLPNUM_TYPE * lower,
	EGLPNUM_TYPE * upper,
	const char **colnames,
	const char **rownames)
{
	int rval = 0;
	EGLPNUM_TYPENAME_QSdata *p = 0;

	p = EGLPNUM_TYPENAME_QScreate_prob (probname, objsense);
	if (p == 0)
		goto CLEANUP;

	rval = EGLPNUM_TYPENAME_ILLlib_newrows (p->lp, 0, nrows, rhs, sense, 0, rownames);
	CHECKRVALG (rval, CLEANUP);

	rval = EGLPNUM_TYPENAME_ILLlib_addcols (p->lp, 0, ncols, cmatcnt, cmatbeg, cmatind,
												 cmatval, obj, lower, upper, colnames, 0);
	CHECKRVALG (rval, CLEANUP);

	p->factorok = 0;

CLEANUP:

	if (rval)
	{
		EGLPNUM_TYPENAME_QSfree_prob (p);
		p = 0;
	}

	return p;
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE EGLPNUM_TYPENAME_QSdata *EGLPNUM_TYPENAME_QScopy_prob (
	EGLPNUM_TYPENAME_QSdata * p,
	const char *newname)
{
	int rval = 0;
	int j, col, beg, pindex, hit;
	EGLPNUM_TYPENAME_QSdata *p2 = 0;
	char *coln;
	char buf[ILL_namebufsize];

	/* QSlog("EGLPNUM_TYPENAME_QScopy_prob ..."); */

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	p2 = EGLPNUM_TYPENAME_QScreate_prob (newname, p->qslp->objsense);
	if (p2 == 0)
		goto CLEANUP;

	rval = EGLPNUM_TYPENAME_ILLlib_newrows (p2->lp, 0, p->qslp->nrows,
												 p->qslp->rhs, p->qslp->sense, p->qslp->rangeval,
												 (const char **) p->qslp->rownames);
	CHECKRVALG (rval, CLEANUP);

	for (j = 0; j < p->qslp->nstruct; j++)
	{
		col = p->qslp->structmap[j];
		if (p->qslp->colnames)
			coln = p->qslp->colnames[j];
		else
			coln = 0;
		beg = p->qslp->A.matbeg[col];

		/* Monika: Note that Java will need to handle these arrays */
		/*         without using the beg offset.  The easiest way  */
		/*         may be to copy the arrays, as in the getcols()  */
		/*         code in lib.c.                                  */

		rval = EGLPNUM_TYPENAME_ILLlib_addcol (p2->lp, 0,
													p->qslp->A.matcnt[col],
													p->qslp->A.matind + beg, p->qslp->A.matval + beg,
													p->qslp->obj[col], p->qslp->lower[col],
													p->qslp->upper[col], coln, 0);
		CHECKRVALG (rval, CLEANUP);
	}

	p2->qslp->objsense = p->qslp->objsense;

	p2->factorok = 0;
	p2->simplex_display = p->simplex_display;
	p2->simplex_scaling = p->simplex_scaling;
	EGLPNUM_TYPENAME_EGlpNumClearVar (p2->pricing->htrigger);
	*(p2->pricing) = *(p->pricing);
	/* I added this line because copying the EGLPNUM_TYPENAME_heap (as a pointer) doesn't make any
	 * sense ! */
	EGLPNUM_TYPENAME_ILLheap_init (&(p2->pricing->h));
	EGLPNUM_TYPENAME_EGlpNumInitVar (p2->pricing->htrigger);
	EGLPNUM_TYPENAME_EGlpNumCopy (p2->pricing->htrigger, p->pricing->htrigger);

	if (p->qslp->intmarker != 0)
	{
		ILL_SAFE_MALLOC (p2->qslp->intmarker, p->qslp->nstruct, char);

		for (j = 0; j < p->qslp->nstruct; j++)
		{
			p2->qslp->intmarker[j] = p->qslp->intmarker[j];
		}
	}

	if (p->qslp->objname != 0)
	{
		ILL_UTIL_STR (p2->qslp->objname, p->qslp->objname);
	}
	else
	{
		strcpy (buf, "obj");
		rval = ILLsymboltab_uname (&p2->qslp->rowtab, buf, "", NULL);
		CHECKRVALG (rval, CLEANUP);
		ILL_UTIL_STR (p2->qslp->objname, buf);
	}
	rval = ILLsymboltab_register (&p2->qslp->rowtab, p2->qslp->objname,
																-1, &pindex, &hit);
	rval = rval || hit;
	CHECKRVALG (rval, CLEANUP);

	ILLstring_reporter_copy (&p2->qslp->reporter, &p->qslp->reporter);

CLEANUP:

	if (rval)
	{
		EGLPNUM_TYPENAME_QSfree_prob (p2);
		p2 = 0;
	}

	return p2;
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSchange_objsense (
	EGLPNUM_TYPENAME_QSdata * p,
	int newsense)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	if (newsense != QS_MIN && newsense != QS_MAX)
	{
		QSlog("Illegal objective sense %d", newsense);
		rval = 1;
		goto CLEANUP;
	}

	if (p->qslp->objsense != newsense)
	{
		if(newsense == QS_MAX) EGLPNUM_TYPENAME_ILLsimplex_set_bound(p->lp,(const EGLPNUM_TYPE *)(&(p->lobjlim)), newsense);
		else EGLPNUM_TYPENAME_ILLsimplex_set_bound(p->lp,(const EGLPNUM_TYPE*)(&(p->uobjlim)), newsense);
		p->qslp->objsense = newsense;
		free_cache (p);
	}

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int
EGLPNUM_TYPENAME_QSget_itcnt(
	EGLPNUM_TYPENAME_QSdata* p,
	int *pI_iter,
	int *pII_iter,
	int *dI_iter, 
	int *dII_iter,
	int *tot_iter)
{
	int rval = 0;
	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);
	if(pI_iter) *pI_iter = p->itcnt.pI_iter;
	if(pII_iter) *pII_iter = p->itcnt.pII_iter;
	if(dI_iter) *dI_iter = p->itcnt.dI_iter;
	if(dII_iter) *dII_iter = p->itcnt.dII_iter;
	if(tot_iter) *tot_iter = p->itcnt.tot_iter;

CLEANUP:

 	EG_RETURN(rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSget_objsense (
	EGLPNUM_TYPENAME_QSdata * p,
	int *objsense)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	if (objsense)
		*objsense = p->qslp->objsense;

CLEANUP:

	EG_RETURN (rval);
}


EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSnew_col (
	EGLPNUM_TYPENAME_QSdata * p,
	const EGLPNUM_TYPE obj,
	const EGLPNUM_TYPE lower,
	const EGLPNUM_TYPE upper,
	const char *name)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	rval = EGLPNUM_TYPENAME_ILLlib_newcol (p->lp, p->basis, obj, lower, upper, name, p->factorok);
	CHECKRVALG (rval, CLEANUP);

	free_cache (p);

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSadd_cols (
	EGLPNUM_TYPENAME_QSdata * p,
	int num,
	int *cmatcnt,
	int *cmatbeg,
	int *cmatind,
	EGLPNUM_TYPE * cmatval,
	EGLPNUM_TYPE * obj,
	EGLPNUM_TYPE * lower,
	EGLPNUM_TYPE * upper,
	const char **names)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	rval = EGLPNUM_TYPENAME_ILLlib_addcols (p->lp, p->basis, num, cmatcnt, cmatbeg,
												 cmatind, cmatval, obj, lower, upper, names,
												 p->factorok);
	CHECKRVALG (rval, CLEANUP);

	free_cache (p);

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSadd_col (
	EGLPNUM_TYPENAME_QSdata * p,
	int cnt,
	int *cmatind,
	EGLPNUM_TYPE * cmatval,
	EGLPNUM_TYPE obj,
	EGLPNUM_TYPE lower,
	EGLPNUM_TYPE upper,
	const char *name)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	rval = EGLPNUM_TYPENAME_ILLlib_addcol (p->lp, p->basis, cnt, cmatind, cmatval,
												obj, lower, upper, name, p->factorok);
	CHECKRVALG (rval, CLEANUP);

	free_cache (p);

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSnew_row (
	EGLPNUM_TYPENAME_QSdata * p,
	const EGLPNUM_TYPE rhs,
	int sense,
	const char *name)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	rval = EGLPNUM_TYPENAME_ILLlib_newrow (p->lp, p->basis, rhs, sense, EGLPNUM_TYPENAME_zeroLpNum, name);
	CHECKRVALG (rval, CLEANUP);

	p->factorok = 0;
	free_cache (p);

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSadd_ranged_rows (
	EGLPNUM_TYPENAME_QSdata * p,
	int num,
	int *rmatcnt,
	int *rmatbeg,
	int *rmatind,
	const EGLPNUM_TYPE * rmatval,
	const EGLPNUM_TYPE * rhs,
	char *sense,
	const EGLPNUM_TYPE * range,
	const char **names)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	rval = EGLPNUM_TYPENAME_ILLlib_addrows (p->lp, p->basis, num, rmatcnt, rmatbeg,
												 rmatind, rmatval, rhs, sense, range,
												 names, &(p->factorok));
	CHECKRVALG (rval, CLEANUP);

	if (p->factorok == 1 && p->basis->rownorms)
	{
		rval = EGLPNUM_TYPENAME_ILLlib_loadrownorms (p->lp, p->pricing, p->basis->rownorms);
		CHECKRVALG (rval, CLEANUP);
		/* This really should go inside of EGLPNUM_TYPENAME_ILLlib_addrows, once pinf is  */
		/* is moved into the lp struct.                                  */
	}

	free_cache (p);

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSadd_ranged_row (
	EGLPNUM_TYPENAME_QSdata * p,
	int cnt,
	int *rmatind,
	const EGLPNUM_TYPE * rmatval,
	const EGLPNUM_TYPE * rhs,
	int sense,
	const EGLPNUM_TYPE * range,
	const char *name)
{
	int rval = 0;
	int vmatcnt[1];
	int vmatbeg[1];
	char vsense[1];
	const char *vnames[1];

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	vmatcnt[0] = cnt;
	vmatbeg[0] = 0;
	vsense[0] = sense;
	vnames[0] = name;

	rval = EGLPNUM_TYPENAME_QSadd_ranged_rows (p, 1, vmatcnt, vmatbeg, rmatind, rmatval, rhs,
														vsense, range, vnames);
	CHECKRVALG (rval, CLEANUP);

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSadd_rows (
	EGLPNUM_TYPENAME_QSdata * p,
	int num,
	int *rmatcnt,
	int *rmatbeg,
	int *rmatind,
	const EGLPNUM_TYPE * rmatval,
	const EGLPNUM_TYPE * rhs,
	char *sense,
	const char **names)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	rval = EGLPNUM_TYPENAME_ILLlib_addrows (p->lp, p->basis, num, rmatcnt, rmatbeg,
												 rmatind, rmatval, rhs, sense, 0, names,
												 &(p->factorok));
	CHECKRVALG (rval, CLEANUP);

	if (p->factorok == 1 && p->basis->rownorms)
	{
		rval = EGLPNUM_TYPENAME_ILLlib_loadrownorms (p->lp, p->pricing, p->basis->rownorms);
		CHECKRVALG (rval, CLEANUP);
		/* This really should go inside of EGLPNUM_TYPENAME_ILLlib_addrows, once pinf is  */
		/* is moved into the lp struct.                                  */
	}

	free_cache (p);

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSadd_row (
	EGLPNUM_TYPENAME_QSdata * p,
	int cnt,
	int *rmatind,
	const EGLPNUM_TYPE * rmatval,
	const EGLPNUM_TYPE * rhs,
	int sense,
	const char *name)
{
	int rval = 0;
	int vmatcnt[1];
	int vmatbeg[1];
	char vsense[1];
	const char *vnames[1];

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	vmatcnt[0] = cnt;
	vmatbeg[0] = 0;
	vsense[0] = sense;
	vnames[0] = name;

	rval = EGLPNUM_TYPENAME_QSadd_rows (p, 1, vmatcnt, vmatbeg, rmatind, rmatval, rhs, vsense,
										 vnames);
	CHECKRVALG (rval, CLEANUP);

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSdelete_rows (
	EGLPNUM_TYPENAME_QSdata * p,
	int num,
	int *dellist)
{
	int rval = 0;
	int basis_ok = 0;
	int cache_ok = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	rval = EGLPNUM_TYPENAME_ILLlib_delrows (p->lp, p->basis, p->cache, num, dellist, &basis_ok,
												 &cache_ok);
	CHECKRVALG (rval, CLEANUP);

	/* For now, just remove the basis - wait for pivotin */

	if (p->basis && !basis_ok)
	{
		EGLPNUM_TYPENAME_ILLlp_basis_free (p->basis);
		ILL_IFFREE (p->basis, EGLPNUM_TYPENAME_ILLlp_basis);
	}

	p->factorok = 0;

	if (!p->basis || !basis_ok || !cache_ok)
	{
		/* Note: If we only delete basic rows then cached soln is valid */
		free_cache (p);
	}

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSdelete_row (
	EGLPNUM_TYPENAME_QSdata * p,
	int rowindex)
{
	int rval = 0;
	int vdellist[1];

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	vdellist[0] = rowindex;

	rval = EGLPNUM_TYPENAME_QSdelete_rows (p, 1, vdellist);
	CHECKRVALG (rval, CLEANUP);

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSdelete_setrows (
	EGLPNUM_TYPENAME_QSdata * p,
	int *flags)
{
	int rval = 0;
	int j, num = 0;
	int *dellist = 0;
	int nrows;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	nrows = p->qslp->nrows;

	for (j = 0; j < nrows; j++)
	{
		if (flags[j] == 1)
			num++;
	}

	if (num > 0)
	{
		ILL_SAFE_MALLOC (dellist, num, int);

		for (j = 0, num = 0; j < nrows; j++)
		{
			if (flags[j] == 1)
			{
				dellist[num++] = j;
			}
		}

		rval = EGLPNUM_TYPENAME_QSdelete_rows (p, num, dellist);
		CHECKRVALG (rval, CLEANUP);
	}

CLEANUP:

	ILL_IFFREE (dellist, int);

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSdelete_named_row (
	EGLPNUM_TYPENAME_QSdata * p,
	const char *rowname)
{
	int rval = 0;
	int i, vdellist[1];

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	rval = EGLPNUM_TYPENAME_QSget_row_index (p, rowname, &i);
	CHECKRVALG (rval, CLEANUP);

	vdellist[0] = i;

	rval = EGLPNUM_TYPENAME_QSdelete_rows (p, 1, vdellist);
	CHECKRVALG (rval, CLEANUP);

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSdelete_named_rows_list (
	EGLPNUM_TYPENAME_QSdata * p,
	int num,
	const char **rownames)
{
	int rval = 0;
	int i, k;
	int *vdellist = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	if (num > 0)
	{
		ILL_SAFE_MALLOC (vdellist, num, int);

		for (k = 0; k < num; k++)
		{
			rval = EGLPNUM_TYPENAME_QSget_row_index (p, rownames[k], &i);
			CHECKRVALG (rval, CLEANUP);
			vdellist[k] = i;
		}

		rval = EGLPNUM_TYPENAME_QSdelete_rows (p, num, vdellist);
		CHECKRVALG (rval, CLEANUP);
	}

CLEANUP:

	ILL_IFFREE (vdellist, int);

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSdelete_cols (
	EGLPNUM_TYPENAME_QSdata * p,
	int num,
	int *dellist)
{
	int rval = 0;
	int basis_ok;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	rval = EGLPNUM_TYPENAME_ILLlib_delcols (p->lp, p->basis, num, dellist, &basis_ok);
	CHECKRVALG (rval, CLEANUP);

	/* For now, just remove the basis - wait for pivotout */

	if (p->basis && !basis_ok)
	{
		EGLPNUM_TYPENAME_ILLlp_basis_free (p->basis);
		ILL_IFFREE (p->basis, EGLPNUM_TYPENAME_ILLlp_basis);
	}

	p->factorok = 0;
	free_cache (p);

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSdelete_col (
	EGLPNUM_TYPENAME_QSdata * p,
	int colindex)
{
	int rval = 0;
	int vdellist[1];

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	vdellist[0] = colindex;

	rval = EGLPNUM_TYPENAME_QSdelete_cols (p, 1, vdellist);
	CHECKRVALG (rval, CLEANUP);

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSdelete_setcols (
	EGLPNUM_TYPENAME_QSdata * p,
	int *flags)
{
	int rval = 0;
	int j, num = 0;
	int *dellist = 0;
	int ncols;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	ncols = p->qslp->nstruct;

	for (j = 0; j < ncols; j++)
	{
		if (flags[j] == 1)
			num++;
	}

	if (num > 0)
	{
		ILL_SAFE_MALLOC (dellist, num, int);

		for (j = 0, num = 0; j < ncols; j++)
		{
			if (flags[j] == 1)
			{
				dellist[num++] = j;
			}
		}

		rval = EGLPNUM_TYPENAME_QSdelete_cols (p, num, dellist);
		CHECKRVALG (rval, CLEANUP);
	}

CLEANUP:

	ILL_IFFREE (dellist, int);

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSdelete_named_column (
	EGLPNUM_TYPENAME_QSdata * p,
	const char *colname)
{
	int rval = 0;
	int j, vdellist[1];

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	rval = EGLPNUM_TYPENAME_QSget_column_index (p, colname, &j);
	CHECKRVALG (rval, CLEANUP);

	vdellist[0] = j;

	rval = EGLPNUM_TYPENAME_QSdelete_cols (p, 1, vdellist);
	CHECKRVALG (rval, CLEANUP);

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSdelete_named_columns_list (
	EGLPNUM_TYPENAME_QSdata * p,
	int num,
	const char **colnames)
{
	int rval = 0;
	int i, j;
	int *vdellist = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	if (num > 0)
	{
		ILL_SAFE_MALLOC (vdellist, num, int);

		for (i = 0; i < num; i++)
		{
			rval = EGLPNUM_TYPENAME_QSget_column_index (p, colnames[i], &j);
			CHECKRVALG (rval, CLEANUP);
			vdellist[i] = j;
		}

		rval = EGLPNUM_TYPENAME_QSdelete_cols (p, num, vdellist);
		CHECKRVALG (rval, CLEANUP);
	}

CLEANUP:

	ILL_IFFREE (vdellist, int);

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSchange_senses (
	EGLPNUM_TYPENAME_QSdata * p,
	int num,
	int *rowlist,
	char *sense)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	rval = EGLPNUM_TYPENAME_ILLlib_chgsense (p->lp, num, rowlist, sense);
	CHECKRVALG (rval, CLEANUP);

	free_cache (p);

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSchange_range (
	EGLPNUM_TYPENAME_QSdata*p,
	int rowindex,
	EGLPNUM_TYPE range)
{
	int rval = 0; 
	rval = check_qsdata_pointer (p);
	CHECKRVALG(rval, CLEANUP);
	
	rval = EGLPNUM_TYPENAME_ILLlib_chgrange (p->lp, rowindex, range);
	CHECKRVALG(rval, CLEANUP);
	
	p->factorok = 0; /* Sanjeeb -- 050911  */
	free_cache (p);
 
CLEANUP:

	EG_RETURN(rval);

}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSchange_sense (
	EGLPNUM_TYPENAME_QSdata * p,
	int rowindex,
	int sense)
{
	int rval = 0;
	int vrowlist[1];
	char vsenselist[1];

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	vrowlist[0] = rowindex;
	vsenselist[0] = sense;

	rval = EGLPNUM_TYPENAME_QSchange_senses (p, 1, vrowlist, vsenselist);
	CHECKRVALG (rval, CLEANUP);

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE  int EGLPNUM_TYPENAME_QSget_coef (
	EGLPNUM_TYPENAME_QSdata *p,
	int rowindex,
	int colindex,
	EGLPNUM_TYPE* coef)
{
	int rval = 0;
	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);
	rval = EGLPNUM_TYPENAME_ILLlib_getcoef (p->lp, rowindex, colindex, coef);
	CHECKRVALG (rval, CLEANUP);
	
CLEANUP:
	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSchange_coef (
	EGLPNUM_TYPENAME_QSdata * p,
	int rowindex,
	int colindex,
	EGLPNUM_TYPE coef)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	rval = EGLPNUM_TYPENAME_ILLlib_chgcoef (p->lp, rowindex, colindex, coef);
	CHECKRVALG (rval, CLEANUP);

	free_cache (p);

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSchange_objcoef (
	EGLPNUM_TYPENAME_QSdata * p,
	int indx,
	EGLPNUM_TYPE coef)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	rval = EGLPNUM_TYPENAME_ILLlib_chgobj (p->lp, indx, coef);
	CHECKRVALG (rval, CLEANUP);

	free_cache (p);

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSchange_rhscoef (
	EGLPNUM_TYPENAME_QSdata * p,
	int indx,
	EGLPNUM_TYPE coef)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	rval = EGLPNUM_TYPENAME_ILLlib_chgrhs (p->lp, indx, coef);
	CHECKRVALG (rval, CLEANUP);

	free_cache (p);

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSchange_bounds (
	EGLPNUM_TYPENAME_QSdata * p,
	int num,
	int *collist,
	char *lu,
	const EGLPNUM_TYPE * bounds)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	rval = EGLPNUM_TYPENAME_ILLlib_chgbnds (p->lp, num, collist, lu, bounds);
	CHECKRVALG (rval, CLEANUP);

	free_cache (p);

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSchange_bound (
	EGLPNUM_TYPENAME_QSdata * p,
	int indx,
	int lu,
	const EGLPNUM_TYPE bound)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	rval = EGLPNUM_TYPENAME_ILLlib_chgbnd (p->lp, indx, lu, bound);
	CHECKRVALG (rval, CLEANUP);

	free_cache (p);

CLEANUP:

	EG_RETURN (rval);
}

#if 0
		/*
		 * Bico - I removed this on 02.04.22.  I don't think we need to support
		 * this type of interface (the loading via arrays can do the job)
		 */
EGLPNUM_TYPENAME_QSLIB_INTERFACE QSbasis *QScreate_basis (
	int nstruct,
	int nrows)
{
	int rval = 0;
	int i;
	QSbasis *B = 0;

	ILL_SAFE_MALLOC (B, 1, QSbasis);

	B->nstruct = nstruct;
	B->nrows = nrows;
	B->cstat = 0;
	B->rstat = 0;

	if (nstruct)
	{
		ILL_SAFE_MALLOC (B->cstat, nstruct, char);
	}

	if (nrows)
	{
		ILL_SAFE_MALLOC (B->rstat, nrows, char);
	}

	for (i = 0; i < nstruct; i++)
		B->cstat[i] = 0;
	for (i = 0; i < nrows; i++)
		B->rstat[i] = 0;

CLEANUP:

	if (rval)
	{
		EGLPNUM_TYPENAME_QSfree_basis (B);
		B = 0;
	}

	return B;
}
#endif

EGLPNUM_TYPENAME_QSLIB_INTERFACE QSbasis *EGLPNUM_TYPENAME_QSread_basis (
	EGLPNUM_TYPENAME_QSdata * p,
	const char *filename)
{
	int rval = 0;
	QSbasis *qB = 0;
	EGLPNUM_TYPENAME_ILLlp_basis B;

	EGLPNUM_TYPENAME_ILLlp_basis_init (&B);

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	ILL_NEW (qB, QSbasis);
	init_basis (qB);

	rval = EGLPNUM_TYPENAME_ILLlib_readbasis (p->lp, &B, filename);
	CHECKRVALG (rval, CLEANUP);

	rval = illbasis_to_qsbasis (&B, qB);
	CHECKRVALG (rval, CLEANUP);

CLEANUP:

	if (rval && qB)
	{
		EGLPNUM_TYPENAME_QSfree_basis (qB);
		qB = 0;
	}
	EGLPNUM_TYPENAME_ILLlp_basis_free (&B);

	return qB;
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSload_basis (
	EGLPNUM_TYPENAME_QSdata * p,
	QSbasis * B)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	if (B->nstruct != p->qslp->nstruct || B->nrows != p->qslp->nrows)
	{
		QSlog("size of basis does not match lp");
		rval = 1;
		goto CLEANUP;
	}

	if (p->basis == 0)
	{
		ILL_SAFE_MALLOC (p->basis, 1, EGLPNUM_TYPENAME_ILLlp_basis);
		EGLPNUM_TYPENAME_ILLlp_basis_init (p->basis);
	}
	else
	{
		EGLPNUM_TYPENAME_ILLlp_basis_free (p->basis);
	}

	rval = qsbasis_to_illbasis (B, p->basis);
	CHECKRVALG (rval, CLEANUP);

	p->factorok = 0;

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSread_and_load_basis (
	EGLPNUM_TYPENAME_QSdata * p,
	const char *filename)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	if (p->basis == 0)
	{
		ILL_SAFE_MALLOC (p->basis, 1, EGLPNUM_TYPENAME_ILLlp_basis);
		EGLPNUM_TYPENAME_ILLlp_basis_init (p->basis);
	}
	else
	{
		EGLPNUM_TYPENAME_ILLlp_basis_free (p->basis);
	}

	rval = EGLPNUM_TYPENAME_ILLlib_readbasis (p->lp, p->basis, filename);
	CHECKRVALG (rval, CLEANUP);

CLEANUP:

	return rval;
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSload_basis_array (
	EGLPNUM_TYPENAME_QSdata * p,
	char *cstat,
	char *rstat)
{
	int rval = 0;
	int i;
	EGLPNUM_TYPENAME_ILLlp_basis *B;
	EGLPNUM_TYPENAME_ILLlpdata *qslp;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	qslp = p->qslp;

	if (qslp->nstruct > 0 && cstat == 0)
	{
		QSlog("EGLPNUM_TYPENAME_QSload_basis_array called without cstat");
		rval = 1;
		goto CLEANUP;
	}

	if (qslp->nrows > 0 && rstat == 0)
	{
		QSlog("EGLPNUM_TYPENAME_QSload_basis_array called without rstat");
		rval = 1;
		goto CLEANUP;
	}

	if (p->basis == 0)
	{
		ILL_SAFE_MALLOC (p->basis, 1, EGLPNUM_TYPENAME_ILLlp_basis);
		EGLPNUM_TYPENAME_ILLlp_basis_init (p->basis);
	}
	else
	{
		EGLPNUM_TYPENAME_ILLlp_basis_free (p->basis);
	}

	B = p->basis;

	B->nstruct = qslp->nstruct;
	B->nrows = qslp->nrows;
	ILL_SAFE_MALLOC (B->cstat, qslp->nstruct, char);
	ILL_SAFE_MALLOC (B->rstat, qslp->nrows, char);

	for (i = 0; i < qslp->nstruct; i++)
	{
		B->cstat[i] = cstat[i];
	}

	for (i = 0; i < qslp->nrows; i++)
	{
		B->rstat[i] = rstat[i];
	}

	p->factorok = 0;

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSload_basis_and_row_norms_array (
	EGLPNUM_TYPENAME_QSdata * p,
	char *cstat,
	char *rstat,
	EGLPNUM_TYPE * rownorms)
{
	int rval = 0;
	int i, nrows;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	nrows = p->qslp->nrows;

	rval = EGLPNUM_TYPENAME_QSload_basis_array (p, cstat, rstat);
	CHECKRVALG (rval, CLEANUP);
	p->basis->rownorms = EGLPNUM_TYPENAME_EGlpNumAllocArray (nrows);

	for (i = 0; i < nrows; i++)
	{
		EGLPNUM_TYPENAME_EGlpNumCopy (p->basis->rownorms[i], rownorms[i]);
	}

	p->factorok = 0;

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSwrite_basis (
	EGLPNUM_TYPENAME_QSdata * p,
	QSbasis * B,
	const char *filename)
{
	int rval = 0;
	EGLPNUM_TYPENAME_ILLlp_basis iB, *basis = 0;

	EGLPNUM_TYPENAME_ILLlp_basis_init (&iB);

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	if (B)
	{
		rval = qsbasis_to_illbasis (B, &iB);
		CHECKRVALG (rval, CLEANUP);
		basis = &iB;
	}
	else
	{
		if (p->basis == 0)
		{
			QSlog("no basis available in EGLPNUM_TYPENAME_QSwrite_basis");
			rval = 1;
			goto CLEANUP;
		}
		basis = p->basis;
	}

	rval = EGLPNUM_TYPENAME_ILLlib_writebasis (p->lp, basis, filename);
	CHECKRVALG (rval, CLEANUP);


CLEANUP:

	EGLPNUM_TYPENAME_ILLlp_basis_free (basis);
	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE QSbasis *EGLPNUM_TYPENAME_QSget_basis (
	EGLPNUM_TYPENAME_QSdata * p)
{
	int rval = 0;
	QSbasis *B = 0;

	if (p->basis == 0)
	{
		QSlog("no basis available in EGLPNUM_TYPENAME_QSget_basis");
		rval = 1;
		goto CLEANUP;
	}

	ILL_SAFE_MALLOC (B, 1, QSbasis);
	init_basis (B);
	rval = illbasis_to_qsbasis (p->basis, B);
	CHECKRVALG (rval, CLEANUP);

CLEANUP:

	if (rval)
	{
		EGLPNUM_TYPENAME_QSfree_basis (B);
		B = 0;
	}

	return B;
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSget_basis_array (
	EGLPNUM_TYPENAME_QSdata * p,
	char *cstat,
	char *rstat)
{
	int rval = 0;
	int i;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	if (p->basis == 0)
	{
		QSlog("no basis available in EGLPNUM_TYPENAME_QSget_basis_array");
		rval = 1;
		goto CLEANUP;
	}

	for (i = 0; i < p->basis->nstruct; i++)
		cstat[i] = p->basis->cstat[i];
	for (i = 0; i < p->basis->nrows; i++)
		rstat[i] = p->basis->rstat[i];

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSget_basis_and_row_norms_array (
	EGLPNUM_TYPENAME_QSdata * p,
	char *cstat,
	char *rstat,
	EGLPNUM_TYPE * rownorms)
{
	int rval = 0;
	int i;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	if (p->basis == 0)
	{
		QSlog("no basis available");
		rval = 1;
		goto CLEANUP;
	}

	for (i = 0; i < p->basis->nstruct; i++)
		cstat[i] = p->basis->cstat[i];
	for (i = 0; i < p->basis->nrows; i++)
		rstat[i] = p->basis->rstat[i];

	if (p->basis->rownorms == 0)
	{
		QSlog("no row norms available");
		rval = 1;
		goto CLEANUP;
	}
	else
	{
		for (i = 0; i < p->basis->nrows; i++)
		{
			EGLPNUM_TYPENAME_EGlpNumCopy (rownorms[i], p->basis->rownorms[i]);
		}
	}

CLEANUP:

	EG_RETURN (rval);
}

static int illbasis_to_qsbasis (
	EGLPNUM_TYPENAME_ILLlp_basis * B,
	QSbasis * qB)
{
	int rval = 0;
	int i;

	qB->nstruct = B->nstruct;
	qB->nrows = B->nrows;
	ILL_SAFE_MALLOC (qB->cstat, B->nstruct, char);
	ILL_SAFE_MALLOC (qB->rstat, B->nrows, char);

	for (i = 0; i < B->nstruct; i++)
	{
		qB->cstat[i] = B->cstat[i];
	}

	for (i = 0; i < B->nrows; i++)
	{
		qB->rstat[i] = B->rstat[i];
	}

CLEANUP:

	EG_RETURN (rval);
}

static int qsbasis_to_illbasis (
	QSbasis * qB,
	EGLPNUM_TYPENAME_ILLlp_basis * B)
{
	int rval = 0;
	int i;
	int nbas = 0;

	B->nstruct = qB->nstruct;
	B->nrows = qB->nrows;
	ILL_SAFE_MALLOC (B->cstat, qB->nstruct, char);
	ILL_SAFE_MALLOC (B->rstat, qB->nrows, char);

	for (i = 0; i < qB->nstruct; i++)
	{
		if(qB->cstat[i] == QS_COL_BSTAT_BASIC) nbas++;
		B->cstat[i] = qB->cstat[i];
	}

	for (i = 0; i < qB->nrows; i++)
	{
		if(qB->rstat[i] == QS_ROW_BSTAT_BASIC) nbas++;
		B->rstat[i] = qB->rstat[i];
	}

	if(nbas != qB->nrows)
	{
		QSlog("Received basis is not valid, in qsbasis_to_illbasis");
		rval = 1;
		ILL_CLEANUP;
	}

CLEANUP:

	EG_RETURN (rval);
}

int grab_basis (
	EGLPNUM_TYPENAME_QSdata * p)
{
	int rval = 0;
	EGLPNUM_TYPENAME_ILLlp_basis *B = p->basis;
	int nstruct = p->qslp->nstruct;
	int nrows = p->qslp->nrows;

	if (!B)
	{
		ILL_SAFE_MALLOC (p->basis, 1, EGLPNUM_TYPENAME_ILLlp_basis);
		EGLPNUM_TYPENAME_ILLlp_basis_init (p->basis);
		B = p->basis;
	}

	if (nstruct != B->nstruct)
	{
		ILL_IFFREE (B->cstat, char);
		ILL_SAFE_MALLOC (B->cstat, nstruct, char);

		B->nstruct = nstruct;
	}

	if (nrows != B->nrows)
	{
		ILL_IFFREE (B->rstat, char);
		ILL_SAFE_MALLOC (B->rstat, nrows, char);

		B->nrows = nrows;
	}

	rval = EGLPNUM_TYPENAME_ILLlib_getbasis (p->lp, B->cstat, B->rstat);
	CHECKRVALG (rval, CLEANUP);

	EGLPNUM_TYPENAME_EGlpNumFreeArray (B->rownorms);
	EGLPNUM_TYPENAME_EGlpNumFreeArray (B->colnorms);

	if (p->pricing->dII_price == QS_PRICE_DSTEEP)
	{
		B->rownorms = EGLPNUM_TYPENAME_EGlpNumAllocArray (nrows);
		rval = EGLPNUM_TYPENAME_ILLlib_getrownorms (p->lp, p->pricing, B->rownorms);
		if (rval)
		{
/*
            QSlog("no edge norms, continue anyway");
*/
			EGLPNUM_TYPENAME_EGlpNumFreeArray (B->rownorms);
			rval = 0;
		}
	}

CLEANUP:

	if (rval)
	{
		if (B)
		{
			EGLPNUM_TYPENAME_ILLlp_basis_free (B);
			ILL_IFFREE (p->basis, EGLPNUM_TYPENAME_ILLlp_basis);
		}
	}

	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_QSgrab_cache (
	EGLPNUM_TYPENAME_QSdata * p,
	int status)
{
	int rval = 0;
	EGLPNUM_TYPENAME_ILLlp_cache *C = p->cache;
	int nstruct = p->qslp->nstruct;
	int nrows = p->qslp->nrows;
	#if 0
	/* we may need to fix basic status for fixed variables */
	register int i;
	char *const cstat = p->basis->cstat;
	const int *const structmap = p->lp->O->structmap;
	const int sense = p->lp->O->objsense;
	const int *const vtype = p->lp->vtype;
	int *const vstat = p->lp->vstat;
	/* end extra variables needed */
	#endif

	if (C == 0)
	{
		ILL_SAFE_MALLOC (p->cache, 1, EGLPNUM_TYPENAME_ILLlp_cache);
		EGLPNUM_TYPENAME_EGlpNumInitVar (p->cache->val);
		EGLPNUM_TYPENAME_ILLlp_cache_init (p->cache);
		C = p->cache;
	}

	if (nstruct != C->nstruct || nrows != C->nrows)
	{
		EGLPNUM_TYPENAME_ILLlp_cache_free (C);
		rval = EGLPNUM_TYPENAME_ILLlp_cache_alloc (C, nstruct, nrows);
		CHECKRVALG (rval, CLEANUP);
	}

	rval = EGLPNUM_TYPENAME_ILLlib_cache_solution (p->lp, C);
	CHECKRVALG (rval, CLEANUP);

	#if 0
	/* we fix basis status and vstat */
	for( i = nstruct ; i-- ; )
	{
		if(vtype[structmap[i]] != VFIXED) continue;
		if(cstat[i] == QS_COL_BSTAT_BASIC) continue;
		if(( sense > 0 && EGLPNUM_TYPENAME_EGlpNumIsLess(EGLPNUM_TYPENAME_SZERO_TOLER,C->rc[i])) ||
			 ( sense < 0 && !EGLPNUM_TYPENAME_EGlpNumIsLess(EGLPNUM_TYPENAME_SZERO_TOLER,C->rc[i])))
		{
			vstat[structmap[i]] = STAT_LOWER;
			cstat[i] = QS_COL_BSTAT_LOWER;
		}
		else
		{
			vstat[structmap[i]] = STAT_UPPER;
			cstat[i] = QS_COL_BSTAT_UPPER;
		}
	}/* end fix */
	#endif
	C->status = status;

CLEANUP:

	if (rval)
	{
		if (C)
		{
			EGLPNUM_TYPENAME_ILLlp_cache_free (C);
			EGLPNUM_TYPENAME_EGlpNumClearVar (p->cache->val);
			ILL_IFFREE (p->cache, EGLPNUM_TYPENAME_ILLlp_cache);
		}
	}

	EG_RETURN (rval);
}

void free_cache (
	EGLPNUM_TYPENAME_QSdata * p)
{
	if (p->cache)
	{
		EGLPNUM_TYPENAME_ILLlp_cache_free (p->cache);
		EGLPNUM_TYPENAME_EGlpNumClearVar (p->cache->val);
		ILL_IFFREE (p->cache, EGLPNUM_TYPENAME_ILLlp_cache);
	}
	p->qstatus = QS_LP_MODIFIED;
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSget_binv_row (
	EGLPNUM_TYPENAME_QSdata * p,
	int indx,
	EGLPNUM_TYPE * binvrow)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	if(!p->basis)
	{
		QSlog("no active basis in store");
		rval = 1;
		goto CLEANUP;
	}
	if(0>indx || indx >= EGLPNUM_TYPENAME_QSget_rowcount(p))
	{
		QSlog("row index %d outside valid bounds [%d:%d]",
								indx, 0, EGLPNUM_TYPENAME_QSget_rowcount(p)-1);
		rval = 1;
		goto CLEANUP;
	}
	if (p->cache == 0)
	{
		QSlog("LP has not been optimized in EGLPNUM_TYPENAME_QSget_binv_row");
		rval = 1;
		goto CLEANUP;
	}

	rval = EGLPNUM_TYPENAME_ILLlib_tableau (p->lp, indx, binvrow, 0);
	CHECKRVALG (rval, CLEANUP);

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSget_tableau_row (
	EGLPNUM_TYPENAME_QSdata * p,
	int indx,
	EGLPNUM_TYPE * tableaurow)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	if (p->cache == 0)
	{
		QSlog("LP has not been optimized in EGLPNUM_TYPENAME_QSget_tableau_row");
		rval = 1;
		goto CLEANUP;
	}

	rval = EGLPNUM_TYPENAME_ILLlib_tableau (p->lp, indx, 0, tableaurow);
	CHECKRVALG (rval, CLEANUP);

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSget_basis_order (
	EGLPNUM_TYPENAME_QSdata * p,
	int *basorder)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	if (p->cache == 0)
	{
		QSlog("LP has not been optimized in EGLPNUM_TYPENAME_QSget_basis_order");
		rval = 1;
		goto CLEANUP;
	}

	rval = EGLPNUM_TYPENAME_ILLlib_basis_order (p->lp, basorder);
	CHECKRVALG (rval, CLEANUP);

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QScompute_row_norms (
	EGLPNUM_TYPENAME_QSdata * p)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	if (p->pricing->dII_price != QS_PRICE_DSTEEP)
	{
		QSlog("not using dual steepest edge");
		rval = 1;
		goto CLEANUP;
	}

	rval = EGLPNUM_TYPENAME_ILLlib_recompute_rownorms (p->lp, p->pricing);
	CHECKRVALG (rval, CLEANUP);

	rval = grab_basis (p);
	CHECKRVALG (rval, CLEANUP);

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE void EGLPNUM_TYPENAME_QSfree_prob (
	EGLPNUM_TYPENAME_QSdata * p)
{
	if (p)
	{
		EGLPNUM_TYPENAME_EGlpNumClearVar(p->uobjlim);
		EGLPNUM_TYPENAME_EGlpNumClearVar(p->lobjlim);
		if (p->qslp)
		{
			EGLPNUM_TYPENAME_ILLlpdata_free (p->qslp);
			ILL_IFFREE (p->qslp, EGLPNUM_TYPENAME_ILLlpdata);
		}
		if (p->lp)
		{
			EGLPNUM_TYPENAME_ILLsimplex_free_lpinfo (p->lp);
			EGLPNUM_TYPENAME_EGlpNumClearVar (p->lp->objval);
			EGLPNUM_TYPENAME_EGlpNumClearVar (p->lp->pobjval);
			EGLPNUM_TYPENAME_EGlpNumClearVar (p->lp->dobjval);
			EGLPNUM_TYPENAME_EGlpNumClearVar (p->lp->pinfeas);
			EGLPNUM_TYPENAME_EGlpNumClearVar (p->lp->dinfeas);
			EGLPNUM_TYPENAME_EGlpNumClearVar (p->lp->objbound);
			EGLPNUM_TYPENAME_EGlpNumClearVar (p->lp->upd.piv);
			EGLPNUM_TYPENAME_EGlpNumClearVar (p->lp->upd.dty);
			EGLPNUM_TYPENAME_EGlpNumClearVar (p->lp->upd.c_obj);
			EGLPNUM_TYPENAME_EGlpNumClearVar (p->lp->upd.tz);
			ILL_IFFREE (p->lp, EGLPNUM_TYPENAME_lpinfo);
		}
		if (p->basis)
		{
			EGLPNUM_TYPENAME_ILLlp_basis_free (p->basis);
			ILL_IFFREE (p->basis, EGLPNUM_TYPENAME_ILLlp_basis);
		}
		if (p->cache)
		{
			EGLPNUM_TYPENAME_ILLlp_cache_free (p->cache);
			EGLPNUM_TYPENAME_EGlpNumClearVar (p->cache->val);
			ILL_IFFREE (p->cache, EGLPNUM_TYPENAME_ILLlp_cache);
		}
		if (p->pricing)
		{
			EGLPNUM_TYPENAME_EGlpNumClearVar (p->pricing->htrigger);
			EGLPNUM_TYPENAME_ILLprice_free_pricing_info (p->pricing);
			ILL_IFFREE (p->pricing, EGLPNUM_TYPENAME_price_info);
		}
		ILL_IFFREE (p->name, char);

		ILL_IFFREE (p, EGLPNUM_TYPENAME_QSdata);
	}
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE void EGLPNUM_TYPENAME_QSfree_basis (
	QSbasis * B)
{
	if (B)
	{
		ILL_IFFREE (B->rstat, char);
		ILL_IFFREE (B->cstat, char);

		ILL_IFFREE (B, QSbasis);
	}
}

static void init_basis (
	QSbasis * B)
{
	if (B)
	{
		B->nstruct = 0;
		B->nrows = 0;
		B->cstat = 0;
		B->rstat = 0;
	}
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSget_status (
	EGLPNUM_TYPENAME_QSdata * p,
	int *status)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	if (status)
		*status = p->qstatus;

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSget_solution (
	EGLPNUM_TYPENAME_QSdata * p,
	EGLPNUM_TYPE * value,
	EGLPNUM_TYPE * x,
	EGLPNUM_TYPE * pi,
	EGLPNUM_TYPE * slack,
	EGLPNUM_TYPE * rc)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	if (p->cache == 0)
	{
		QSlog("no solution available in EGLPNUM_TYPENAME_QSget_solution");
		rval = 1;
		goto CLEANUP;
	}

	rval = EGLPNUM_TYPENAME_ILLlib_solution (p->lp, p->cache, value, x, pi, slack, rc);
	CHECKRVALG (rval, CLEANUP);

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSget_objval (
	EGLPNUM_TYPENAME_QSdata * p,
	EGLPNUM_TYPE * value)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

/* Want to get objval after limited number of pivots. */

	if (p->qstatus == QS_LP_MODIFIED)
	{
		QSlog("QSmsg: LP has been modified since last solve.");
		rval = 1;
		ILL_CLEANUP;
	}

	rval = EGLPNUM_TYPENAME_ILLlib_objval (p->lp, p->cache, value);
	CHECKRVALG (rval, CLEANUP);

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSget_x_array (
	EGLPNUM_TYPENAME_QSdata * p,
	EGLPNUM_TYPE * x)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	if (p->cache == 0)
	{
		QSlog("no solution available in EGLPNUM_TYPENAME_QSget_x_array");
		rval = 1;
		goto CLEANUP;
	}

	rval = EGLPNUM_TYPENAME_ILLlib_get_x (p->lp, p->cache, x);
	CHECKRVALG (rval, CLEANUP);

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSget_slack_array (
	EGLPNUM_TYPENAME_QSdata * p,
	EGLPNUM_TYPE * slack)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	if (p->cache == 0)
	{
		QSlog("no solution available in EGLPNUM_TYPENAME_QSget_slack_array");
		rval = 1;
		goto CLEANUP;
	}

	rval = EGLPNUM_TYPENAME_ILLlib_get_slack (p->lp, p->cache, slack);
	CHECKRVALG (rval, CLEANUP);

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSget_rc_array (
	EGLPNUM_TYPENAME_QSdata * p,
	EGLPNUM_TYPE * rc)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	if (p->cache == 0)
	{
		QSlog("no solution available in EGLPNUM_TYPENAME_QSget_rc_array");
		rval = 1;
		goto CLEANUP;
	}

	rval = EGLPNUM_TYPENAME_ILLlib_solution (p->lp, p->cache, 0, 0, 0, 0, rc);
	CHECKRVALG (rval, CLEANUP);

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSget_pi_array (
	EGLPNUM_TYPENAME_QSdata * p,
	EGLPNUM_TYPE * pi)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	if (p->cache == 0)
	{
		QSlog("no solution available in EGLPNUM_TYPENAME_QSget_pi_array");
		rval = 1;
		goto CLEANUP;
	}

	rval = EGLPNUM_TYPENAME_ILLlib_solution (p->lp, p->cache, 0, 0, pi, 0, 0);
	CHECKRVALG (rval, CLEANUP);

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSget_infeas_array (
	EGLPNUM_TYPENAME_QSdata * p,
	EGLPNUM_TYPE * pi)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	if (pi == 0)
	{
		ILL_ERROR (rval, "QS_get_infeas_array called with NULL pi vector\n");
	}

	rval = EGLPNUM_TYPENAME_ILLsimplex_infcertificate (p->lp, pi);
	CHECKRVALG (rval, CLEANUP);

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSget_named_x (
	EGLPNUM_TYPENAME_QSdata * p,
	const char *colname,
	EGLPNUM_TYPE * val)
{
	int rval = 0;
	int j;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	if (p->cache == 0)
	{
		QSlog("no solution available in EGLPNUM_TYPENAME_QSget_named_x");
		rval = 1;
		goto CLEANUP;
	}

	rval = EGLPNUM_TYPENAME_QSget_column_index (p, colname, &j);
	CHECKRVALG (rval, CLEANUP);

	if (j != -1)
	{
		EGLPNUM_TYPENAME_EGlpNumCopy (*val, p->cache->x[j]);
	}
	else
	{
		rval = 1;
	}

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSget_named_rc (
	EGLPNUM_TYPENAME_QSdata * p,
	const char *colname,
	EGLPNUM_TYPE * val)
{
	int rval = 0;
	int j;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	if (p->cache == 0)
	{
		QSlog("no solution available in EGLPNUM_TYPENAME_QSget_named_rc");
		rval = 1;
		goto CLEANUP;
	}

	rval = EGLPNUM_TYPENAME_QSget_column_index (p, colname, &j);
	CHECKRVALG (rval, CLEANUP);

	if (j != -1)
	{
		EGLPNUM_TYPENAME_EGlpNumCopy (*val, p->cache->rc[j]);
	}
	else
	{
		rval = 1;
	}

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSget_named_pi (
	EGLPNUM_TYPENAME_QSdata * p,
	const char *rowname,
	EGLPNUM_TYPE * val)
{
	int rval = 0;
	int i;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	if (p->cache == 0)
	{
		QSlog("no solution available in EGLPNUM_TYPENAME_QSget_named_pi");
		rval = 1;
		goto CLEANUP;
	}

	rval = EGLPNUM_TYPENAME_QSget_row_index (p, rowname, &i);
	CHECKRVALG (rval, CLEANUP);

	if (i != -1)
	{
		EGLPNUM_TYPENAME_EGlpNumCopy (*val, p->cache->pi[i]);
	}
	else
	{
		rval = 1;
	}

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSget_named_slack (
	EGLPNUM_TYPENAME_QSdata * p,
	const char *rowname,
	EGLPNUM_TYPE * val)
{
	int rval = 0;
	int i;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	if (p->cache == 0)
	{
		QSlog("no solution available in EGLPNUM_TYPENAME_QSget_named_slack");
		rval = 1;
		goto CLEANUP;
	}

	rval = EGLPNUM_TYPENAME_QSget_row_index (p, rowname, &i);
	CHECKRVALG (rval, CLEANUP);

	if (i != -1)
	{
		EGLPNUM_TYPENAME_EGlpNumCopy (*val, p->cache->slack[i]);
	}
	else
	{
		rval = 1;
	}

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSget_colcount (
	EGLPNUM_TYPENAME_QSdata * p)
{
	int cnt;

	if (check_qsdata_pointer (p))
		cnt = 0;
	else
		cnt = p->qslp->nstruct;

	return cnt;
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSget_rowcount (
	EGLPNUM_TYPENAME_QSdata * p)
{
	int cnt;

	if (check_qsdata_pointer (p))
		cnt = 0;
	else
		cnt = p->qslp->nrows;

	return cnt;
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSget_nzcount (
	EGLPNUM_TYPENAME_QSdata * p)
{
	int cnt;

	if (check_qsdata_pointer (p))
		cnt = 0;
	else
		cnt = p->qslp->nzcount - p->qslp->nrows;

	return cnt;
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QStest_row_norms (
	EGLPNUM_TYPENAME_QSdata * p)
{
	int yesno;

	if (check_qsdata_pointer (p))
	{
		yesno = 0;
	}
	else
	{
		if (p->basis && p->basis->rownorms)
		{
			yesno = 1;
		}
		else
		{
			yesno = 0;
		}
	}

	return yesno;
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSget_obj_list(EGLPNUM_TYPENAME_QSprob p,
	int num,
	int*collist,
	EGLPNUM_TYPE*obj)
{
	int rval = 0;
	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);
	rval = EGLPNUM_TYPENAME_ILLlib_getobj_list (p->lp, num, collist, obj);
	CHECKRVALG (rval, CLEANUP);

CLEANUP:
	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSget_obj (
	EGLPNUM_TYPENAME_QSdata * p,
	EGLPNUM_TYPE * obj)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	rval = EGLPNUM_TYPENAME_ILLlib_getobj (p->lp, obj);
	CHECKRVALG (rval, CLEANUP);

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSget_rhs (
	EGLPNUM_TYPENAME_QSdata * p,
	EGLPNUM_TYPE * rhs)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	rval = EGLPNUM_TYPENAME_ILLlib_getrhs (p->lp, rhs);
	CHECKRVALG (rval, CLEANUP);

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSget_ranged_rows_list (
	EGLPNUM_TYPENAME_QSdata * p,
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
	int i, nrows;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	nrows = EGLPNUM_TYPENAME_QSget_rowcount (p);
	for (i = 0; i < num; i++)
	{
		if (rowlist[i] < 0 || rowlist[i] >= nrows)
		{
			QSlog("entry %d in rowlist out of range", i);
			rval = 1;
			goto CLEANUP;
		}
	}

	rval = EGLPNUM_TYPENAME_ILLlib_getrows (p->lp, num, rowlist, rowcnt, rowbeg, rowind,
												 rowval, rhs, sense, range, names);
	CHECKRVALG (rval, CLEANUP);


CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSget_ranged_rows (
	EGLPNUM_TYPENAME_QSdata * p,
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
	int i, nrows;
	int *rowlist = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	nrows = EGLPNUM_TYPENAME_QSget_rowcount (p);
	if (nrows > 0)
	{
		ILL_SAFE_MALLOC (rowlist, nrows, int);

		for (i = 0; i < nrows; i++)
		{
			rowlist[i] = i;
		}
		rval = EGLPNUM_TYPENAME_ILLlib_getrows (p->lp, nrows, rowlist, rowcnt, rowbeg, rowind,
													 rowval, rhs, sense, range, names);
		CHECKRVALG (rval, CLEANUP);
	}

CLEANUP:

	ILL_IFFREE (rowlist, int);

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSget_senses (
	EGLPNUM_TYPENAME_QSdata *p,
	char *senses)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	rval = EGLPNUM_TYPENAME_ILLlib_getsenses (p->lp, senses);
	CHECKRVALG (rval, CLEANUP);

CLEANUP:

	EG_RETURN (rval);
}



EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSget_rows_list (
	EGLPNUM_TYPENAME_QSdata * p,
	int num,
	int *rowlist,
	int **rowcnt,
	int **rowbeg,
	int **rowind,
	EGLPNUM_TYPE ** rowval,
	EGLPNUM_TYPE ** rhs,
	char **sense,
	char ***names)
{
	int rval = 0;
	int i, nrows;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	nrows = EGLPNUM_TYPENAME_QSget_rowcount (p);
	for (i = 0; i < num; i++)
	{
		if (rowlist[i] < 0 || rowlist[i] >= nrows)
		{
			QSlog("entry %d in rowlist out of range", i);
			rval = 1;
			goto CLEANUP;
		}
	}

	rval = EGLPNUM_TYPENAME_ILLlib_getrows (p->lp, num, rowlist, rowcnt, rowbeg, rowind,
												 rowval, rhs, sense, 0, names);
	CHECKRVALG (rval, CLEANUP);


CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSget_rows (
	EGLPNUM_TYPENAME_QSdata * p,
	int **rowcnt,
	int **rowbeg,
	int **rowind,
	EGLPNUM_TYPE ** rowval,
	EGLPNUM_TYPE ** rhs,
	char **sense,
	char ***names)
{
	int rval = 0;
	int i, nrows;
	int *rowlist = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	nrows = EGLPNUM_TYPENAME_QSget_rowcount (p);
	if (nrows > 0)
	{
		ILL_SAFE_MALLOC (rowlist, nrows, int);

		for (i = 0; i < nrows; i++)
		{
			rowlist[i] = i;
		}
		rval = EGLPNUM_TYPENAME_ILLlib_getrows (p->lp, nrows, rowlist, rowcnt, rowbeg, rowind,
													 rowval, rhs, sense, 0, names);
		CHECKRVALG (rval, CLEANUP);
	}

CLEANUP:

	ILL_IFFREE (rowlist, int);

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSget_columns_list (
	EGLPNUM_TYPENAME_QSdata * p,
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
	int j, ncols;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	ncols = EGLPNUM_TYPENAME_QSget_colcount (p);
	for (j = 0; j < num; j++)
	{
		if (collist[j] < 0 || collist[j] >= ncols)
		{
			QSlog("entry %d in collist out of range", j);
			rval = 1;
			goto CLEANUP;
		}
	}

	rval = EGLPNUM_TYPENAME_ILLlib_getcols (p->lp, num, collist, colcnt, colbeg, colind,
												 colval, obj, lower, upper, names);
	CHECKRVALG (rval, CLEANUP);

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSget_columns (
	EGLPNUM_TYPENAME_QSdata * p,
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
	int j, ncols;
	int *collist = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	ncols = EGLPNUM_TYPENAME_QSget_colcount (p);
	if (ncols > 0)
	{
		ILL_SAFE_MALLOC (collist, ncols, int);

		for (j = 0; j < ncols; j++)
		{
			collist[j] = j;
		}
		rval = EGLPNUM_TYPENAME_ILLlib_getcols (p->lp, ncols, collist, colcnt, colbeg, colind,
													 colval, obj, lower, upper, names);
		CHECKRVALG (rval, CLEANUP);
	}

CLEANUP:

	ILL_IFFREE (collist, int);

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE char *EGLPNUM_TYPENAME_QSget_probname (
	EGLPNUM_TYPENAME_QSdata * p)
{
	int rval = 0;
	char *name = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	ILL_UTIL_STR (name, p->name);

CLEANUP:
	ILL_RETURN_PTR (name, "EGLPNUM_TYPENAME_QSget_probname");
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE char *EGLPNUM_TYPENAME_QSget_objname (
	EGLPNUM_TYPENAME_QSdata * p)
{
	int rval = 0;
	char *name = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	if (p->qslp->objname != 0)
	{
		ILL_UTIL_STR (name, p->qslp->objname);
	}

CLEANUP:
	ILL_RETURN_PTR (name, "EGLPNUM_TYPENAME_QSget_objname");
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSget_rownames (
	EGLPNUM_TYPENAME_QSdata * p,
	char **rownames)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	rval = EGLPNUM_TYPENAME_ILLlib_rownames (p->lp, rownames);
	CHECKRVALG (rval, CLEANUP);

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSget_colnames (
	EGLPNUM_TYPENAME_QSdata * p,
	char **colnames)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	rval = EGLPNUM_TYPENAME_ILLlib_colnames (p->lp, colnames);
	CHECKRVALG (rval, CLEANUP);

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSget_bound (
	EGLPNUM_TYPENAME_QSdata * p,
	int colindex,
	int lu,
	EGLPNUM_TYPE * bound)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	rval = EGLPNUM_TYPENAME_ILLlib_getbnd (p->lp, colindex, lu, bound);
	CHECKRVALG (rval, CLEANUP);

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSget_bounds_list(
	EGLPNUM_TYPENAME_QSdata* p,
	int num,
	int*collist,
	EGLPNUM_TYPE*lb,
	EGLPNUM_TYPE*ub)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	rval = EGLPNUM_TYPENAME_ILLlib_getbnds_list (p->lp, num, collist, lb, ub);
	CHECKRVALG (rval, CLEANUP);
	
CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSget_bounds (
	EGLPNUM_TYPENAME_QSdata * p,
	EGLPNUM_TYPE * lower,
	EGLPNUM_TYPE * upper)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	rval = EGLPNUM_TYPENAME_ILLlib_getbnds (p->lp, lower, upper);
	CHECKRVALG (rval, CLEANUP);

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSget_intflags (
	EGLPNUM_TYPENAME_QSdata * p,
	int *intflags)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	if (intflags == 0)
	{
		rval = 1;
		ILL_CLEANUP;
	}

	rval = EGLPNUM_TYPENAME_ILLlib_getintflags (p->lp, intflags);
	CHECKRVALG (rval, CLEANUP);

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSget_intcount (
	EGLPNUM_TYPENAME_QSdata * p,
	int *count)
{
	int j, ncols, cnt = 0, rval = 0;
	int *intflags = 0;

	*count = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	ncols = EGLPNUM_TYPENAME_QSget_colcount (p);

	if (ncols > 0)
	{
		ILL_SAFE_MALLOC (intflags, ncols, int);

		rval = EGLPNUM_TYPENAME_ILLlib_getintflags (p->lp, intflags);
		CHECKRVALG (rval, CLEANUP);

		for (j = 0; j < ncols; j++)
		{
			if (intflags[j] > 0)
				cnt++;
		}
	}

CLEANUP:

	*count = cnt;
	ILL_IFFREE (intflags, int);

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSget_column_index (
	EGLPNUM_TYPENAME_QSdata * p,
	const char *name,
	int *colindex)
{
	int rval = 0;

	*colindex = -1;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	rval = EGLPNUM_TYPENAME_ILLlib_colindex (p->lp, name, colindex);
	CHECKRVALG (rval, CLEANUP);

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSget_row_index (
	EGLPNUM_TYPENAME_QSdata * p,
	const char *name,
	int *rowindex)
{
	int rval = 0;

	*rowindex = -1;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	rval = EGLPNUM_TYPENAME_ILLlib_rowindex (p->lp, name, rowindex);
	CHECKRVALG (rval, CLEANUP);

CLEANUP:

	EG_RETURN (rval);
}

static int QSwrite_prob_EGioFile(
		EGLPNUM_TYPENAME_QSdata*p,
		EGioFile_t*out,
		const char*filetype)
{
    int rval = 0; 
    qsstring_reporter rep; 

    ILLstring_reporter_copy(&rep, &p->qslp->reporter); 
    ILLstring_reporter_init(&p->qslp->reporter, 
                            (qsreport_string_fct) EGioWrite, out); 
    rval = EGLPNUM_TYPENAME_QSreport_prob(p, filetype, NULL);
    ILLstring_reporter_copy(&p->qslp->reporter, &rep); 
    ILL_RESULT(rval, "QSwrite_prob_EGioFile");
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSwrite_prob (
	EGLPNUM_TYPENAME_QSdata * p,
	const char *filename,
	const char *filetype)
{
	EGioFile_t *file = NULL;
	int rval = 0;

	if (filename == NULL)
	{
		file = EGioOpenFILE(stdout);
	}
	else
	{
		if ((file = EGioOpen (filename, "w")) == 0)
		{
			QSlog("Unable to open \"%s\" for output: %s.", filename,
									strerror(errno));
		}
	}
	ILL_CHECKnull (file, NULL);
	rval = QSwrite_prob_EGioFile (p, file, filetype);
	if (file)
	{
		EGioClose (file);
	}
CLEANUP:
	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSwrite_prob_file (
	EGLPNUM_TYPENAME_QSdata * p,
	FILE * out,
	const char *filetype)
{
	int rval = 0;
	EGioFile_t*lout = EGioOpenFILE(out);
	rval = QSwrite_prob_EGioFile(p,lout,filetype);
	CHECKRVALG(rval,CLEANUP);
	CLEANUP:
	free(lout);
	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE void EGLPNUM_TYPENAME_QSfree (
	void *ptr)
{
	ILL_IFFREE (ptr, void);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSset_param (
	EGLPNUM_TYPENAME_QSdata * p,
	int whichparam,
	int newvalue)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	switch (whichparam)
	{
	case QS_PARAM_PRIMAL_PRICING:
		if (newvalue == QS_PRICE_PDANTZIG ||
				newvalue == QS_PRICE_PDEVEX ||
				newvalue == QS_PRICE_PSTEEP || newvalue == QS_PRICE_PMULTPARTIAL)
		{
			p->pricing->pI_price = newvalue;
			p->pricing->pII_price = newvalue;
		}
		else
		{
			QSlog("illegal value for QS_PARAM_PRIMAL_PRICING");
			rval = 1;
			goto CLEANUP;
		}
		break;
	case QS_PARAM_DUAL_PRICING:
		if (newvalue == QS_PRICE_DDANTZIG ||
				newvalue == QS_PRICE_DSTEEP ||
				newvalue == QS_PRICE_DMULTPARTIAL || newvalue == QS_PRICE_DDEVEX)
		{
			p->pricing->dI_price = newvalue;
			p->pricing->dII_price = newvalue;
		}
		else
		{
			QSlog("illegal value for QS_PARAM_DUAL_PRICING");
			rval = 1;
			goto CLEANUP;
		}
		break;
	case QS_PARAM_SIMPLEX_DISPLAY:
		if (newvalue == 0 || (newvalue > 0 && newvalue < 4))
		{
			p->simplex_display = newvalue;
		}
		else
		{
			QSlog("illegal value for QS_PARAM_SIMPLEX_DISPLAY");
			rval = 1;
			goto CLEANUP;
		}
		break;
	case QS_PARAM_SIMPLEX_MAX_ITERATIONS:
		if (newvalue > 0)
		{
			p->lp->maxiter = newvalue;
		}
		else
		{
			QSlog("illegal value for QS_PARAM_SIMPLEX_MAX_ITERATIONS");
			rval = 1;
			goto CLEANUP;
		}
		break;
	case QS_PARAM_SIMPLEX_SCALING:
		if (newvalue == 0 || newvalue == 1)
		{
			p->simplex_scaling = newvalue;
		}
		else
		{
			QSlog("illegal value for QS_PARAM_SIMPLEX_SCALING");
			rval = 1;
			goto CLEANUP;
		}
		break;
	default:
		QSlog("unknown parameter: %d", whichparam);
		rval = 1;
		goto CLEANUP;
	}

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSset_param_EGlpNum (
	EGLPNUM_TYPENAME_QSdata * p,
	int whichparam,
	EGLPNUM_TYPE newvalue)
{
	int rval = 0;
	int sense;
	EGLPNUM_TYPE lvar;

	EGLPNUM_TYPENAME_EGlpNumInitVar(lvar);
	EGLPNUM_TYPENAME_EGlpNumCopy(lvar,newvalue);
	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	switch (whichparam)
	{
	case QS_PARAM_SIMPLEX_MAX_TIME:
		if (EGLPNUM_TYPENAME_EGlpNumIsGreatZero (lvar))
		{
			p->lp->maxtime = EGLPNUM_TYPENAME_EGlpNumToLf (lvar);
		}
		else
		{
			QSlog("illegal value for QS_PARAM_SIMPLEX_MAX_TIME");
			rval = 1;
			goto CLEANUP;
		}
		break;
	case QS_PARAM_OBJULIM:
		EGLPNUM_TYPENAME_QSget_objsense(p,&sense);
		if(EGLPNUM_TYPENAME_EGlpNumIsLeq(EGLPNUM_TYPENAME_ILL_MAXDOUBLE,lvar)) EGLPNUM_TYPENAME_EGlpNumCopy(lvar,EGLPNUM_TYPENAME_ILL_MAXDOUBLE);
		EGLPNUM_TYPENAME_EGlpNumCopy(p->uobjlim,lvar);
		if(sense == QS_MIN) EGLPNUM_TYPENAME_ILLsimplex_set_bound(p->lp,(const EGLPNUM_TYPE*)(&lvar), sense);
		break;
	case QS_PARAM_OBJLLIM:
		EGLPNUM_TYPENAME_QSget_objsense(p,&sense);
		if(EGLPNUM_TYPENAME_EGlpNumIsLeq(newvalue,EGLPNUM_TYPENAME_ILL_MINDOUBLE)) EGLPNUM_TYPENAME_EGlpNumCopy(lvar,EGLPNUM_TYPENAME_ILL_MINDOUBLE);
		EGLPNUM_TYPENAME_EGlpNumCopy(p->lobjlim,lvar);
		if(sense == QS_MAX) EGLPNUM_TYPENAME_ILLsimplex_set_bound(p->lp,(const EGLPNUM_TYPE*)(&lvar), sense);
		break;
	default:
		QSlog("unknown parameter: %d", whichparam);
		rval = 1;
		goto CLEANUP;
	}

CLEANUP:

	EGLPNUM_TYPENAME_EGlpNumClearVar(lvar);
	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSget_param (
	EGLPNUM_TYPENAME_QSdata * p,
	int whichparam,
	int *value)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	if (!value)
	{
		QSlog("EGLPNUM_TYPENAME_QSget_param call without a value pointer");
		rval = 1;
		goto CLEANUP;
	}

	switch (whichparam)
	{
	case QS_PARAM_PRIMAL_PRICING:
		*value = p->pricing->pII_price;
		break;
	case QS_PARAM_DUAL_PRICING:
		*value = p->pricing->dII_price;
		break;
	case QS_PARAM_SIMPLEX_DISPLAY:
		*value = p->simplex_display;
		break;
	case QS_PARAM_SIMPLEX_MAX_ITERATIONS:
		*value = p->lp->maxiter;
		break;
	case QS_PARAM_SIMPLEX_SCALING:
		*value = p->simplex_scaling;
		break;
	default:
		QSlog("unknown parameter: %d", whichparam);
		rval = 1;
		goto CLEANUP;
	}

CLEANUP:

	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSget_param_EGlpNum (
	EGLPNUM_TYPENAME_QSdata * p,
	int whichparam,
	EGLPNUM_TYPE * value)
{
	int rval = 0;

	rval = check_qsdata_pointer (p);
	CHECKRVALG (rval, CLEANUP);

	if (!value)
	{
		QSlog("QSget_param_double call without a value pointer");
		rval = 1;
		goto CLEANUP;
	}

	switch (whichparam)
	{
	case QS_PARAM_SIMPLEX_MAX_TIME:
		EGLPNUM_TYPENAME_EGlpNumSet (*value, p->lp->maxtime);
		break;
	case QS_PARAM_OBJULIM:
		EGLPNUM_TYPENAME_EGlpNumCopy(*value,p->uobjlim);
		break;
	case QS_PARAM_OBJLLIM:
		EGLPNUM_TYPENAME_EGlpNumCopy(*value,p->lobjlim);
		break;
	default:
		QSlog("unknown parameter: %d", whichparam);
		rval = 1;
		goto CLEANUP;
	}

CLEANUP:

	EG_RETURN (rval);
}

static int check_qsdata_pointer (
	EGLPNUM_TYPENAME_QSdata * p)
{
	if (p == NULL)
	{
		QSlog("NULL EGLPNUM_TYPENAME_QSprob pointer");
		return 1;
	}
	else
	{
		return 0;
	}
}

static int formatIsMps (
	const char *filetype,
	int *isMps)
{
	int rval = 0;

	if (!strcasecmp (filetype, "MPS"))
	{
		*isMps = 1;
	}
	else if (!strcasecmp (filetype, "LP"))
	{
		*isMps = 0;
	}
	else
	{
		QSlog("Unknown prob-file type: %s", filetype);
		rval = 1;
		ILL_CLEANUP;
	}
CLEANUP:
	return rval;
}


/****************************************************************************/
/* 
 * undocumentyed functions 
 */

static void check_pointer (
	void *p,
	const char *fct,
	const char *param)
{
	if (p == NULL)
		QSlog("NULL %s argument to %s", param, fct);
}

/* EGLPNUM_TYPENAME_QSline_reader: 
 *    used by mps/lp reader to get input lines 
 *    by default input is read froma FILE* via fgets 
 */
EGLPNUM_TYPENAME_QSLIB_INTERFACE EGLPNUM_TYPENAME_QSline_reader EGLPNUM_TYPENAME_QSline_reader_new (
	void *fct,
	void *data_src)
{
	check_pointer (fct, "EGLPNUM_TYPENAME_QSline_reader_new", "fct");
	check_pointer (data_src, "EGLPNUM_TYPENAME_QSline_reader_new", "data_src");
	return EGLPNUM_TYPENAME_ILLline_reader_new ((EGLPNUM_TYPENAME_qsread_line_fct) fct, data_src);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE void EGLPNUM_TYPENAME_QSline_reader_set_error_collector (
	EGLPNUM_TYPENAME_QSline_reader reader,
	EGLPNUM_TYPENAME_QSerror_collector collector)
{
	check_pointer (reader, "EGLPNUM_TYPENAME_QSline_reader_set_error_collector", "reader");
	check_pointer (collector, "EGLPNUM_TYPENAME_QSline_reader_set_error_collector", "collector");
	reader->error_collector = collector;
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE void EGLPNUM_TYPENAME_QSline_reader_free (
	EGLPNUM_TYPENAME_QSline_reader reader)
{
	EGLPNUM_TYPENAME_ILLline_reader_free (reader);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE char *EGLPNUM_TYPENAME_QSline_reader_get (
	EGLPNUM_TYPENAME_QSline_reader reader,
	char *s,
	int size)
{
	check_pointer (reader, "EGLPNUM_TYPENAME_QSline_reader_get", "reader");
	check_pointer (s, "EGLPNUM_TYPENAME_QSline_reader_get", "s");
	return EGLPNUM_TYPENAME_ILLline_reader_get (s, size, reader);
}


EGLPNUM_TYPENAME_QSLIB_INTERFACE EGLPNUM_TYPENAME_QSerror_collector EGLPNUM_TYPENAME_QSerror_collector_new (
	void *fct,
	void *dest)
{
	check_pointer (fct, "EGLPNUM_TYPENAME_QSerror_collector_new", "fct");
	return EGLPNUM_TYPENAME_ILLerror_collector_new ((EGLPNUM_TYPENAME_qsadd_error_fct) fct, dest);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE
	EGLPNUM_TYPENAME_QSerror_collector EGLPNUM_TYPENAME_QSerror_memory_collector_new (EGLPNUM_TYPENAME_QSerror_memory mem)
{
	check_pointer (mem, "EGLPNUM_TYPENAME_QSerror_memory_collector_new", "mem");
	return EGLPNUM_TYPENAME_ILLerror_memory_collector_new (mem);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE void EGLPNUM_TYPENAME_QSerror_collector_free (
	EGLPNUM_TYPENAME_QSerror_collector c)
{
	EGLPNUM_TYPENAME_ILLerror_collector_free (c);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE EGLPNUM_TYPENAME_QSdata *EGLPNUM_TYPENAME_QSget_prob (
	EGLPNUM_TYPENAME_QSline_reader reader,
	const char *probname,
	const char *filetype)
{
	int isMps, rval = 0;
	EGLPNUM_TYPENAME_QSdata *p = 0;

	if ((filetype != NULL) && !strcasecmp (filetype, "MPS"))
	{
		isMps = 1;
	}
	else if ((filetype != NULL) && !strcasecmp (filetype, "LP"))
	{
		isMps = 0;
	}
	else
	{
		QSlog("Unknown prob-file type: %s",
								(filetype != NULL) ? filetype : "NULL");
		rval = 1;
		ILL_CLEANUP;
	}

	p = EGLPNUM_TYPENAME_ILLread (reader, probname, isMps);
	ILL_CHECKnull (p, NULL);

	ILL_FAILfalse (p->qslp != NULL, "If there's a p there must be a p-qslp");
	ILL_IFFREE (p->name, char);

	ILL_UTIL_STR (p->name, p->qslp->probname);
	EGLPNUM_TYPENAME_ILLsimplex_load_lpinfo (p->qslp, p->lp);

CLEANUP:

	if (rval != 0)
	{
		EGLPNUM_TYPENAME_QSfree_prob (p);
		p = 0;
	}
	return p;
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSreport_prob (
	EGLPNUM_TYPENAME_QSdata * p,
	const char *filetype,
	EGLPNUM_TYPENAME_qserror_collector * c)
{
	int isMps, rval = 0;

	rval = formatIsMps (filetype, &isMps);
	CHECKRVALG (rval, CLEANUP);
	if (isMps)
	{
		rval = EGLPNUM_TYPENAME_ILLwrite_mps (p->qslp, c);
	}
	else
	{
		rval = EGLPNUM_TYPENAME_ILLwrite_lp (p->qslp, c);
	}
CLEANUP:
	EG_RETURN (rval);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE char *EGLPNUM_TYPENAME_QSversion (
	void)
{
	char *name = 0;
	name = EGsMalloc (char, 256);

	snprintf (name, (size_t) 255, "%s (build %s-%s)",PACKAGE_STRING, __DATE__,
						__TIME__);
	return name;
}

/* QSstring_reporter: 
 *    used by solver code to report feedback 
 *    by default feedback is sent to logging
 */
EGLPNUM_TYPENAME_QSLIB_INTERFACE void EGLPNUM_TYPENAME_QSset_reporter (
	EGLPNUM_TYPENAME_QSprob prob,
	int skip,
	void *fct,
	void *dest)
{
	int rval = 0;

	rval = check_qsdata_pointer (prob);
	if (rval != 0)
		return;

	check_pointer (fct, "EGLPNUM_TYPENAME_QSset_reporter", "fct");

	ILL_FAILtrue (prob->lp == NULL, "EGLPNUM_TYPENAME_QSprob internal error: prob->lp == NULL");
	ILLstring_reporter_init (&prob->qslp->reporter,
													 (qsreport_string_fct) fct, dest);

	prob->lp->iterskip = skip;
CLEANUP:
	return;
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE const char *EGLPNUM_TYPENAME_QSformat_error_type_string (
	int tp)
{
	const char *type = "Error";

	if (tp == QS_DATA_ERROR)
	{
		type = "Data Error";
	}
	if (tp == QS_DATA_WARN)
	{
		type = "Data Warning";
	}
	if (tp == QS_MPS_FORMAT_ERROR)
	{
		type = "MPS Error";
	}
	if (tp == QS_MPS_FORMAT_WARN)
	{
		type = "MPS Warning";
	}
	if (tp == QS_LP_FORMAT_ERROR)
	{
		type = "LP Error";
	}
	if (tp == QS_LP_FORMAT_WARN)
	{
		type = "LP Warning";
	}
	return type;
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSerror_get_type (
	EGLPNUM_TYPENAME_QSformat_error error)
{
	check_pointer (error, "EGLPNUM_TYPENAME_QSerror_get_type", "error");
	return error->type;
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE const char *EGLPNUM_TYPENAME_QSerror_get_desc (
	EGLPNUM_TYPENAME_QSformat_error error)
{
	check_pointer (error, "EGLPNUM_TYPENAME_QSerror_get_desc", "error");
	return error->desc;
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSerror_get_line_number (
	EGLPNUM_TYPENAME_QSformat_error error)
{
	check_pointer (error, "EGLPNUM_TYPENAME_QSerror_get_line_number", "error");
	return error->lineNumber;
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSerror_get_pos (
	EGLPNUM_TYPENAME_QSformat_error error)
{
	check_pointer (error, "EGLPNUM_TYPENAME_QSerror_get_pos", "error");
	return error->at;
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE const char *EGLPNUM_TYPENAME_QSerror_get_line (
	EGLPNUM_TYPENAME_QSformat_error error)
{
	check_pointer (error, "EGLPNUM_TYPENAME_QSerror_get_line", "error");
	return error->theLine;
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE void EGLPNUM_TYPENAME_QSerror_print (
	FILE * f,
	EGLPNUM_TYPENAME_QSformat_error error)
{
	check_pointer (f, "EGLPNUM_TYPENAME_QSerror_print", "f");
	if (error == NULL)
	{
		QSlog("0");
	}
	else
	{
		EGioFile_t*out = EGioOpenFILE(f);
		EGLPNUM_TYPENAME_ILLformat_error_print (out, error);
		EGioClose(out);
	}
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE EGLPNUM_TYPENAME_QSerror_memory EGLPNUM_TYPENAME_QSerror_memory_create (
	int takeErrorLines)
{
	return EGLPNUM_TYPENAME_ILLerror_memory_create (takeErrorLines);
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE void EGLPNUM_TYPENAME_QSerror_memory_free (
	EGLPNUM_TYPENAME_QSerror_memory mem)
{
	if (mem != NULL)
	{
		EGLPNUM_TYPENAME_ILLerror_memory_free (mem);
	}
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSerror_memory_get_nerrors (
	EGLPNUM_TYPENAME_QSerror_memory mem)
{
	check_pointer (mem, "EGLPNUM_TYPENAME_QSerror_memory_get_nerrors", "mem");
	return mem->nerror;
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSerror_memory_get_nof (
	EGLPNUM_TYPENAME_QSerror_memory mem,
	int type)
{
	check_pointer (mem, "EGLPNUM_TYPENAME_QSerror_memory_get_nerrors", "mem");
	if (0 <= type && type < QS_INPUT_NERROR)
	{
		return mem->has_error[type];
	}
	else
	{
		ILL_REPRT ("bad error type");
		return 0;
	}
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE EGLPNUM_TYPENAME_QSformat_error EGLPNUM_TYPENAME_QSerror_memory_get_last_error (
	EGLPNUM_TYPENAME_QSerror_memory mem)
{
	check_pointer (mem, "QSerror_memory_get_last_errors", "mem");
	return mem->error_list;
}

EGLPNUM_TYPENAME_QSLIB_INTERFACE EGLPNUM_TYPENAME_QSformat_error EGLPNUM_TYPENAME_QSerror_memory_get_prev_error (
	EGLPNUM_TYPENAME_QSformat_error e)
{
	check_pointer (e, "QSerror_memory_get_prev_errors", "e");
	if (e != NULL)
		e = e->next;
	return e;
}
