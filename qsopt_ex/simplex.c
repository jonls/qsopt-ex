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

/* RCS_INFO = "$RCSfile: simplex.c,v $ $Revision: 1.2 $ $Date: 2003/11/05 16:49:52 $"; */
static int TRACE = 0;

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <stdlib.h>
#include <string.h>

#include "logging-private.h"

#include "eg_lpnum.h"
#include "eg_io.h"
#include "except.h"
#include "zeit.h"
#include "util.h"

#define QSOPT_CURRENT_PRECICION
#include "basicdefs.h"
#include "lpdata_EGLPNUM_TYPENAME.h"
#include "lpdefs_EGLPNUM_TYPENAME.h"

#include "stddefs.h"
#include "fct_EGLPNUM_TYPENAME.h"
#include "ratio_EGLPNUM_TYPENAME.h"
#include "price_EGLPNUM_TYPENAME.h"
#include "basis_EGLPNUM_TYPENAME.h"
#include "simplex_EGLPNUM_TYPENAME.h"
#include "dstruct_EGLPNUM_TYPENAME.h"
#include "qstruct_EGLPNUM_TYPENAME.h"
#include "qsopt_EGLPNUM_TYPENAME.h"
#include "lib_EGLPNUM_TYPENAME.h"								/* for EGLPNUM_TYPENAME_ILLlib_writebasis */
#include "lp_EGLPNUM_TYPENAME.h"									/* for EGLPNUM_TYPENAME_ILLwrite_lp */


static void init_lp_status_info (
	EGLPNUM_TYPENAME_lp_status_info * ls),
  init_simplex_tols (
	EGLPNUM_TYPENAME_lpinfo * lp),
  monitor_iter (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_price_info * p,
	EGLPNUM_TYPENAME_iter_info * it,
	int cphase),
  get_current_stat (
	EGLPNUM_TYPENAME_lp_status_info * p,
	int algorithm,
	int *bstat);

static int terminate_simplex (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int phase,
	EGLPNUM_TYPENAME_iter_info * it),
  primal_phaseI_step (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_price_info * pinf,
	EGLPNUM_TYPENAME_svector * updz,
	EGLPNUM_TYPENAME_svector * wz,
	EGLPNUM_TYPENAME_iter_info * it),
  primal_phaseII_step (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_price_info * pinf,
	EGLPNUM_TYPENAME_svector * updz,
	EGLPNUM_TYPENAME_svector * wz,
	EGLPNUM_TYPENAME_iter_info * it),
  dual_phaseI_step (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_price_info * pinf,
	EGLPNUM_TYPENAME_svector * updz,
	EGLPNUM_TYPENAME_svector * wz,
	EGLPNUM_TYPENAME_iter_info * it),
  dual_phaseII_step (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_price_info * pinf,
	EGLPNUM_TYPENAME_svector * updz,
	EGLPNUM_TYPENAME_svector * wz,
	EGLPNUM_TYPENAME_iter_info * it),
  report_value (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_iter_info * it,
	const char *value_name,
	EGLPNUM_TYPE value);


void EGLPNUM_TYPENAME_ILLsimplex_init_lpinfo (
	EGLPNUM_TYPENAME_lpinfo * lp)
{
	EGLPNUM_TYPENAME_ILLbasis_init_basisinfo (lp);
	EGLPNUM_TYPENAME_init_internal_lpinfo (lp);
}

void EGLPNUM_TYPENAME_ILLsimplex_free_lpinfo (
	EGLPNUM_TYPENAME_lpinfo * lp)
{
	if (lp)
	{
		EGLPNUM_TYPENAME_EGlpNumFreeArray (lp->lz);
		EGLPNUM_TYPENAME_EGlpNumFreeArray (lp->uz);
		EGLPNUM_TYPENAME_EGlpNumFreeArray (lp->cz);
		EGLPNUM_TYPENAME_ILLbasis_free_basisinfo (lp);
		EGLPNUM_TYPENAME_free_internal_lpinfo (lp);
	}
}

void EGLPNUM_TYPENAME_ILLsimplex_load_lpinfo (
	EGLPNUM_TYPENAME_ILLlpdata * qslp,
	EGLPNUM_TYPENAME_lpinfo * lp)
{
	lp->basisid = -1;
	lp->maxiter = 500000;
	lp->maxtime = 300000;
	//lp->iterskip = 10;
	lp->iterskip = 100;
	EGLPNUM_TYPENAME_EGlpNumCopy (lp->objbound, EGLPNUM_TYPENAME_INFTY);
	lp->O = qslp;
}

void EGLPNUM_TYPENAME_ILLsimplex_set_bound (
	EGLPNUM_TYPENAME_lpinfo * lp,
	const EGLPNUM_TYPE * objbound,
	int sense)
{
	EGLPNUM_TYPENAME_EGlpNumCopy (lp->objbound, *objbound);
	if (sense == EGLPNUM_TYPENAME_ILL_MAX)
		EGLPNUM_TYPENAME_EGlpNumSign (lp->objbound);
}

static void init_lp_status_info (
	EGLPNUM_TYPENAME_lp_status_info * ls)
{
	ls->optimal = 0;
	ls->primal_feasible = 0;
	ls->primal_infeasible = 0;
	ls->primal_unbounded = 0;
	ls->dual_feasible = 0;
	ls->dual_infeasible = 0;
	ls->dual_unbounded = 0;
}

static void init_simplex_tols (
	EGLPNUM_TYPENAME_lpinfo * lp)
{
	EGLPNUM_TYPENAME_EGlpNumCopy (lp->tol->pfeas_tol, EGLPNUM_TYPENAME_PFEAS_TOLER);
	EGLPNUM_TYPENAME_EGlpNumCopy (lp->tol->dfeas_tol, EGLPNUM_TYPENAME_DFEAS_TOLER);
	EGLPNUM_TYPENAME_EGlpNumCopy (lp->tol->pivot_tol, EGLPNUM_TYPENAME_PIVOT_TOLER);
	EGLPNUM_TYPENAME_EGlpNumCopy (lp->tol->szero_tol, EGLPNUM_TYPENAME_SZERO_TOLER);
	EGLPNUM_TYPENAME_EGlpNumCopy (lp->tol->ip_tol, lp->tol->pfeas_tol);
	EGLPNUM_TYPENAME_EGlpNumCopy (lp->tol->id_tol, lp->tol->dfeas_tol);
	if (EGLPNUM_TYPENAME_EGlpNumIsNeqqZero (lp->tol->ip_tol))
	{
#if VERBOSE_LEVEL <= DEBUG
		MESSAGE (VERBOSE_LEVEL, "ip_tol %lg", EGLPNUM_TYPENAME_EGlpNumToLf (lp->tol->ip_tol));
		MESSAGE (VERBOSE_LEVEL, "eps %lg", EGLPNUM_TYPENAME_EGlpNumToLf (EGLPNUM_TYPENAME_epsLpNum));
		MESSAGE (VERBOSE_LEVEL, "EGLPNUM_TYPENAME_PFEAS_TOLER %lg", EGLPNUM_TYPENAME_EGlpNumToLf (EGLPNUM_TYPENAME_PFEAS_TOLER));
#endif
		EGLPNUM_TYPENAME_EGlpNumDivUiTo (lp->tol->ip_tol, 2UL);
	}
	if (EGLPNUM_TYPENAME_EGlpNumIsNeqqZero (lp->tol->id_tol))
	{
#if VERBOSE_LEVEL <= DEBUG
		MESSAGE (VERBOSE_LEVEL, "id_tol %lg", EGLPNUM_TYPENAME_EGlpNumToLf (lp->tol->id_tol));
#endif
		EGLPNUM_TYPENAME_EGlpNumDivUiTo (lp->tol->id_tol, 2UL);
	}
}

void EGLPNUM_TYPENAME_init_internal_lpinfo (
	EGLPNUM_TYPENAME_lpinfo * lp)
{
	int rval = 0;

	lp->nrows = 0;
	lp->nnbasic = 0;
	lp->localrows = 0;
	lp->rowcnt = 0;
	lp->rowbeg = 0;
	lp->rowind = 0;
	lp->rowval = 0;
	lp->cz = 0;
	lp->lz = 0;
	lp->uz = 0;
	lp->xbz = 0;
	lp->piz = 0;
	lp->dz = 0;
	lp->pIxbz = 0;
	lp->pIpiz = 0;
	lp->pIdz = 0;
	lp->vtype = 0;
	lp->vclass = 0;
	lp->iwork = 0;
	lp->upd.perm = 0;
	lp->upd.ix = 0;
	lp->upd.t = 0;
	lp->bfeas = 0;
	lp->dfeas = 0;
	lp->tol = 0;
	lp->cnts = 0;
	lp->bchanges = 0;
	lp->cchanges = 0;
	EGLPNUM_TYPENAME_ILLsvector_init (&(lp->zz));
	EGLPNUM_TYPENAME_ILLsvector_init (&(lp->yjz));
	EGLPNUM_TYPENAME_ILLsvector_init (&(lp->zA));
	EGLPNUM_TYPENAME_ILLsvector_init (&(lp->work));
	EGLPNUM_TYPENAME_ILLsvector_init (&(lp->srhs));
	EGLPNUM_TYPENAME_ILLsvector_init (&(lp->ssoln));
	ILL_SAFE_MALLOC (lp->tol, 1, EGLPNUM_TYPENAME_tol_struct);
	EGLPNUM_TYPENAME_EGlpNumInitVar (lp->tol->pfeas_tol);
	EGLPNUM_TYPENAME_EGlpNumInitVar (lp->tol->dfeas_tol);
	EGLPNUM_TYPENAME_EGlpNumInitVar (lp->tol->pivot_tol);
	EGLPNUM_TYPENAME_EGlpNumInitVar (lp->tol->szero_tol);
	EGLPNUM_TYPENAME_EGlpNumInitVar (lp->tol->ip_tol);
	EGLPNUM_TYPENAME_EGlpNumInitVar (lp->tol->id_tol);
	ILL_SAFE_MALLOC (lp->cnts, 1, EGLPNUM_TYPENAME_count_struct);
	EGLPNUM_TYPENAME_EGlpNumInitVar (lp->cnts->y_ravg);
	EGLPNUM_TYPENAME_EGlpNumInitVar (lp->cnts->z_ravg);
	EGLPNUM_TYPENAME_EGlpNumInitVar (lp->cnts->za_ravg);
CLEANUP:
	if (rval)
	{
		QSlog("no memory, in %s, exit", __func__);
		exit (1);
	}
}

void EGLPNUM_TYPENAME_free_internal_lpinfo (
	EGLPNUM_TYPENAME_lpinfo * lp)
{
	EGLPNUM_TYPENAME_bndinfo *binfo = 0;
	EGLPNUM_TYPENAME_coefinfo *cinfo = 0;

	if (lp->localrows)
	{
		ILL_IFFREE (lp->rowcnt, int);
		ILL_IFFREE (lp->rowbeg, int);
		ILL_IFFREE (lp->rowind, int);

		EGLPNUM_TYPENAME_EGlpNumFreeArray (lp->rowval);
		lp->localrows = 0;
	}
	EGLPNUM_TYPENAME_EGlpNumFreeArray (lp->lz);
	EGLPNUM_TYPENAME_EGlpNumFreeArray (lp->uz);
	EGLPNUM_TYPENAME_EGlpNumFreeArray (lp->cz);
	EGLPNUM_TYPENAME_EGlpNumFreeArray (lp->xbz);
	EGLPNUM_TYPENAME_EGlpNumFreeArray (lp->piz);
	EGLPNUM_TYPENAME_EGlpNumFreeArray (lp->pIpiz);
	EGLPNUM_TYPENAME_EGlpNumFreeArray (lp->dz);
	EGLPNUM_TYPENAME_EGlpNumFreeArray (lp->pIdz);
	EGLPNUM_TYPENAME_EGlpNumFreeArray (lp->pIxbz);

	ILL_IFFREE (lp->vtype, int);
	ILL_IFFREE (lp->vclass, char);

	EGLPNUM_TYPENAME_ILLsvector_free (&(lp->zz));
	EGLPNUM_TYPENAME_ILLsvector_free (&(lp->yjz));
	EGLPNUM_TYPENAME_ILLsvector_free (&(lp->zA));
	EGLPNUM_TYPENAME_ILLsvector_free (&(lp->work));
	EGLPNUM_TYPENAME_ILLsvector_free (&(lp->srhs));
	EGLPNUM_TYPENAME_ILLsvector_free (&(lp->ssoln));
	ILL_IFFREE (lp->iwork, int);
	ILL_IFFREE (lp->upd.perm, int);
	ILL_IFFREE (lp->upd.ix, int);

	EGLPNUM_TYPENAME_EGlpNumFreeArray (lp->upd.t);

	ILL_IFFREE (lp->bfeas, int);
	ILL_IFFREE (lp->dfeas, int);

	if (lp->tol)
	{
		EGLPNUM_TYPENAME_EGlpNumClearVar (lp->tol->pfeas_tol);
		EGLPNUM_TYPENAME_EGlpNumClearVar (lp->tol->dfeas_tol);
		EGLPNUM_TYPENAME_EGlpNumClearVar (lp->tol->pivot_tol);
		EGLPNUM_TYPENAME_EGlpNumClearVar (lp->tol->szero_tol);
		EGLPNUM_TYPENAME_EGlpNumClearVar (lp->tol->ip_tol);
		EGLPNUM_TYPENAME_EGlpNumClearVar (lp->tol->id_tol);
		ILL_IFFREE (lp->tol, EGLPNUM_TYPENAME_tol_struct);
	}
	if (lp->cnts)
	{
		EGLPNUM_TYPENAME_EGlpNumClearVar (lp->cnts->y_ravg);
		EGLPNUM_TYPENAME_EGlpNumClearVar (lp->cnts->z_ravg);
		EGLPNUM_TYPENAME_EGlpNumClearVar (lp->cnts->za_ravg);
		ILL_IFFREE (lp->cnts, EGLPNUM_TYPENAME_count_struct);
	}

	while (lp->bchanges)
	{
		binfo = lp->bchanges;
		EGLPNUM_TYPENAME_EGlpNumClearVar (binfo->pbound);
		EGLPNUM_TYPENAME_EGlpNumClearVar (binfo->cbound);
		lp->bchanges = binfo->next;
		ILL_IFFREE (binfo, EGLPNUM_TYPENAME_bndinfo);
	}

	while (lp->cchanges)
	{
		cinfo = lp->cchanges;
		EGLPNUM_TYPENAME_EGlpNumClearVar (cinfo->pcoef);
		EGLPNUM_TYPENAME_EGlpNumClearVar (cinfo->ccoef);
		lp->cchanges = cinfo->next;
		ILL_IFFREE (cinfo, EGLPNUM_TYPENAME_coefinfo);
	}
}

int EGLPNUM_TYPENAME_build_internal_lpinfo (
	EGLPNUM_TYPENAME_lpinfo * lp)
{
	int rval = 0;
	int i, n;
	EGLPNUM_TYPENAME_ILLlpdata *qslp = lp->O;
	EGLPNUM_TYPENAME_ILLlp_sinfo *S = lp->O->sinfo;
	EGLPNUM_TYPE *lower, *upper, *obj;
	EGLPNUM_TYPENAME_ILLlp_rows lprows;
	EGLPNUM_TYPENAME_ILLmatrix *A;

	init_lp_status_info (&(lp->probstat));
	init_lp_status_info (&(lp->basisstat));

	if (S != 0)
	{
		lp->nrows = S->nrows;
		lp->ncols = S->ncols;
		lp->bz = S->rhs;
		lower = S->lower;
		upper = S->upper;
		obj = S->obj;
		A = &(S->A);
	}
	else
	{
		lp->nrows = qslp->nrows;
		lp->ncols = qslp->ncols;
		lp->bz = qslp->rhs;
		lower = qslp->lower;
		upper = qslp->upper;
		obj = qslp->obj;
		A = &(qslp->A);
	}

	lp->matbeg = A->matbeg;
	lp->matcnt = A->matcnt;
	lp->matind = A->matind;
	lp->matval = A->matval;

	lp->nnbasic = lp->ncols - lp->nrows;

	lp->lz = EGLPNUM_TYPENAME_EGlpNumAllocArray (lp->ncols);
	lp->cz = EGLPNUM_TYPENAME_EGlpNumAllocArray (lp->ncols);
	lp->uz = EGLPNUM_TYPENAME_EGlpNumAllocArray (lp->ncols);
	if (!lp->lz || !lp->uz || !lp->cz)
	{
		QSlog("EGLPNUM_TYPENAME_build_internal_lpinfo");
		rval = 1;
		goto CLEANUP;
	}
	for (i = 0; i < lp->ncols; i++)
	{
		EGLPNUM_TYPENAME_EGlpNumCopy (lp->lz[i], lower[i]);
		EGLPNUM_TYPENAME_EGlpNumCopy (lp->uz[i], upper[i]);
		EGLPNUM_TYPENAME_EGlpNumCopy (lp->cz[i], obj[i]);
		if (qslp->objsense == EGLPNUM_TYPENAME_ILL_MAX)
		{
			EGLPNUM_TYPENAME_EGlpNumSign (lp->cz[i]);
		}
	}

	if (!lp->O->rA)
	{
		rval = EGLPNUM_TYPENAME_ILLlp_rows_init (&lprows, lp->O, 1);
		CHECKRVALG (rval, CLEANUP);
		lp->rowbeg = lprows.rowbeg;
		lp->rowcnt = lprows.rowcnt;
		lp->rowind = lprows.rowind;
		lp->rowval = lprows.rowval;
		lp->localrows = 1;
	}
	else
	{
		/* row format exists, just use pointers */
		lp->rowbeg = lp->O->rA->rowbeg;
		lp->rowcnt = lp->O->rA->rowcnt;
		lp->rowind = lp->O->rA->rowind;
		lp->rowval = lp->O->rA->rowval;
		lp->localrows = 0;
	}

	lp->xbz = EGLPNUM_TYPENAME_EGlpNumAllocArray (lp->nrows);
	lp->piz = EGLPNUM_TYPENAME_EGlpNumAllocArray (lp->nrows);
	lp->dz = EGLPNUM_TYPENAME_EGlpNumAllocArray (lp->nnbasic);
	lp->final_phase = -1;
	lp->infub_ix = -1;

	ILL_SAFE_MALLOC (lp->vtype, lp->ncols, int);
	ILL_SAFE_MALLOC (lp->vclass, lp->ncols, char);

	rval = EGLPNUM_TYPENAME_ILLsvector_alloc (&(lp->zz), lp->nrows);
	CHECKRVALG (rval, CLEANUP);
	rval = EGLPNUM_TYPENAME_ILLsvector_alloc (&(lp->yjz), lp->nrows);
	CHECKRVALG (rval, CLEANUP);
	rval = EGLPNUM_TYPENAME_ILLsvector_alloc (&(lp->zA), lp->nnbasic);
	CHECKRVALG (rval, CLEANUP);
	rval = EGLPNUM_TYPENAME_ILLsvector_alloc (&(lp->work), lp->ncols);
	CHECKRVALG (rval, CLEANUP);
	rval = EGLPNUM_TYPENAME_ILLsvector_alloc (&(lp->srhs), lp->nrows);
	CHECKRVALG (rval, CLEANUP);
	rval = EGLPNUM_TYPENAME_ILLsvector_alloc (&(lp->ssoln), lp->nrows);
	CHECKRVALG (rval, CLEANUP);
	ILL_SAFE_MALLOC (lp->iwork, lp->ncols, int);

	for (i = 0; i < lp->ncols; i++)
	{
		lp->work.indx[i] = 0;
		EGLPNUM_TYPENAME_EGlpNumZero (lp->work.coef[i]);
		lp->iwork[i] = 0;
	}
	n = lp->nrows > lp->ncols ? 2 * (lp->nrows) + 1 : 2 * (lp->ncols) + 1;
	lp->upd.t = EGLPNUM_TYPENAME_EGlpNumAllocArray (n);
	ILL_SAFE_MALLOC (lp->upd.perm, n, int);
	ILL_SAFE_MALLOC (lp->upd.ix, n, int);


	ILL_SAFE_MALLOC (lp->bfeas, lp->nrows, int);
	ILL_SAFE_MALLOC (lp->dfeas, lp->nnbasic, int);

	init_simplex_tols (lp);
	EGLPNUM_TYPENAME_ILLfct_init_counts (lp);

	lp->nbchange = 0;
	lp->ncchange = 0;

	lp->pIratio = RATIOTEST_HARRIS;
	lp->pIIratio = RATIOTEST_HARRIS;
	lp->dIratio = RATIOTEST_HARRIS;
	lp->dIIratio = RATIOTEST_HARRIS;
	lp->starttime = ILLutil_zeit ();
	ILLutil_sprand (1, &(lp->rstate));

CLEANUP:
	if (rval)
		EGLPNUM_TYPENAME_free_internal_lpinfo (lp);
	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLsimplex_retest_psolution (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_price_info * p,
	int phase,
	EGLPNUM_TYPENAME_feas_info * fi)
{
	int rval = 0;
	int fbid = lp->fbasisid;
	int bid = lp->basisid;
	EGLPNUM_TYPE *ptol = &(lp->tol->pfeas_tol);
	EGLPNUM_TYPE *dtol = &(lp->tol->dfeas_tol);
	EGLPNUM_TYPE *iptol = &(lp->tol->ip_tol);
	EGLPNUM_TYPE *idtol = &(lp->tol->id_tol);

	fi->pstatus = -1;
	fi->dstatus = -1;
	if (fbid < bid - PARAM_PRIMAL_REFACTORGAP)
	{
		rval = EGLPNUM_TYPENAME_ILLbasis_refactor (lp);
		CHECKRVALG (rval, CLEANUP);
	}
	if (fbid < bid - PARAM_PRIMAL_RESOLVEGAP)
		EGLPNUM_TYPENAME_ILLfct_compute_xbz (lp);

	if (phase == PRIMAL_PHASEII)
	{
		if (fbid < bid - PARAM_PRIMAL_RESOLVEGAP)
		{
			EGLPNUM_TYPENAME_ILLfct_compute_piz (lp);
			EGLPNUM_TYPENAME_ILLfct_compute_dz (lp);
			if (p != NULL && p->p_strategy == COMPLETE_PRICING)
				EGLPNUM_TYPENAME_ILLprice_compute_dual_inf (lp, p, NULL, 0, PRIMAL_PHASEII);
		}
		EGLPNUM_TYPENAME_ILLfct_compute_pobj (lp);
		EGLPNUM_TYPENAME_ILLfct_check_pfeasible (lp, fi, *ptol);
		EGLPNUM_TYPENAME_ILLfct_check_dfeasible (lp, fi, *dtol);
	}
	else if (phase == PRIMAL_PHASEI)
	{
		EGLPNUM_TYPENAME_ILLfct_check_pfeasible (lp, fi, *iptol);
		if (fi->pstatus != PRIMAL_FEASIBLE)
		{
			if (lp->pIpiz)
			{
				EGLPNUM_TYPENAME_ILLfct_compute_phaseI_piz (lp);
				EGLPNUM_TYPENAME_ILLfct_compute_phaseI_dz (lp);
				EGLPNUM_TYPENAME_ILLfct_check_pIdfeasible (lp, fi, *idtol);
				if (p != NULL && p->p_strategy == COMPLETE_PRICING)
					EGLPNUM_TYPENAME_ILLprice_compute_dual_inf (lp, p, NULL, 0, PRIMAL_PHASEI);
			}
		}
	}
CLEANUP:
	if (rval == QS_LP_CHANGE_PREC)
	{
		MESSAGE (__QS_SB_VERB, "Changing precision");
		return rval;
	}
	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLsimplex_retest_dsolution (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_price_info * p,
	int phase,
	EGLPNUM_TYPENAME_feas_info * fi)
{
	int rval = 0;
	int fbid = lp->fbasisid;
	int bid = lp->basisid;
	EGLPNUM_TYPE *ptol = &(lp->tol->pfeas_tol);
	EGLPNUM_TYPE *dtol = &(lp->tol->dfeas_tol);
	EGLPNUM_TYPE *iptol = &(lp->tol->ip_tol);
	EGLPNUM_TYPE *idtol = &(lp->tol->id_tol);

	//ILL_IFTRACE("\t%s:%s:%d\n",__func__,__FILE__,__LINE__);

	fi->pstatus = -1;
	fi->dstatus = -1;
	if (fbid < bid - PARAM_DUAL_REFACTORGAP)
	{
		//ILL_IFTRACE("Refactor: %s:%s:%d\n",__func__,__FILE__,__LINE__);
		rval = EGLPNUM_TYPENAME_ILLbasis_refactor (lp);
		CHECKRVALG (rval, CLEANUP);
	}
	if (fbid < bid - PARAM_DUAL_RESOLVEGAP)
	{
		//ILL_IFTRACE("\t%s:%s:%d\n",__func__,__FILE__,__LINE__);
		EGLPNUM_TYPENAME_ILLfct_compute_piz (lp);
		EGLPNUM_TYPENAME_ILLfct_compute_dz (lp);
	}

	if (phase == DUAL_PHASEII)
	{
		//ILL_IFTRACE("\t%s:%s:%d\n",__func__,__FILE__,__LINE__);
		if (fbid < bid - PARAM_DUAL_RESOLVEGAP)
		{
			//ILL_IFTRACE("\t%s:%s:%d\n",__func__,__FILE__,__LINE__);
			EGLPNUM_TYPENAME_ILLfct_compute_xbz (lp);
			CHECKRVALG (rval, CLEANUP);
			if (p != NULL)
			{
				//ILL_IFTRACE("\t%s:%s:%d\n",__func__,__FILE__,__LINE__);
				if (p->d_strategy == COMPLETE_PRICING)
				{
					//ILL_IFTRACE("\t%s:%s:%d\n",__func__,__FILE__,__LINE__);
					EGLPNUM_TYPENAME_ILLprice_compute_primal_inf (lp, p, NULL, 0, DUAL_PHASEII);
				}
				else
				{
					//ILL_IFTRACE("\t%s:%s:%d\n",__func__,__FILE__,__LINE__);
					EGLPNUM_TYPENAME_ILLprice_update_mpartial_price (lp, p, DUAL_PHASEII, ROW_PRICING);
				}
			}
		}
		//ILL_IFTRACE("\t%s:%s:%d\n",__func__,__FILE__,__LINE__);
		EGLPNUM_TYPENAME_ILLfct_compute_dobj (lp);
		EGLPNUM_TYPENAME_ILLfct_check_dfeasible (lp, fi, *dtol);
		EGLPNUM_TYPENAME_ILLfct_check_pfeasible (lp, fi, *ptol);
	}
	else if (phase == DUAL_PHASEI)
	{
		EGLPNUM_TYPENAME_ILLfct_check_dfeasible (lp, fi, *idtol);
		if (fi->dstatus != DUAL_FEASIBLE)
		{
			EGLPNUM_TYPENAME_ILLfct_compute_phaseI_xbz (lp);
			EGLPNUM_TYPENAME_ILLfct_check_pIpfeasible (lp, fi, *iptol);
			if (p != NULL)
			{
				if (p->d_strategy == COMPLETE_PRICING)
					EGLPNUM_TYPENAME_ILLprice_compute_primal_inf (lp, p, NULL, 0, DUAL_PHASEI);
				else
					EGLPNUM_TYPENAME_ILLprice_update_mpartial_price (lp, p, DUAL_PHASEI, ROW_PRICING);
			}
		}
	}
CLEANUP:
	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLsimplex_solution (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPE * xz,
	EGLPNUM_TYPE * piz,
	EGLPNUM_TYPE * dz,
	EGLPNUM_TYPE * objval)
{
	int i, j;
	int col;

	if (xz != NULL)
	{
		if (lp->basisstat.optimal == 0)
		{
			EG_RETURN (1);
		}
		for (i = 0; i < lp->nrows; i++)
			EGLPNUM_TYPENAME_EGlpNumCopy (xz[lp->baz[i]], lp->xbz[i]);
		for (j = 0; j < lp->nnbasic; j++)
		{
			col = lp->nbaz[j];
			if (lp->vstat[col] == STAT_UPPER)
				EGLPNUM_TYPENAME_EGlpNumCopy (xz[col], lp->uz[col]);
			else if (lp->vstat[col] == STAT_LOWER)
				EGLPNUM_TYPENAME_EGlpNumCopy (xz[col], lp->lz[col]);
			else
				EGLPNUM_TYPENAME_EGlpNumZero (xz[col]);
		}
	}
	if (piz != NULL)
	{
		if (lp->basisstat.optimal == 0)
		{
			EG_RETURN (1);
		}
		for (i = 0; i < lp->nrows; i++)
			EGLPNUM_TYPENAME_EGlpNumCopy (piz[i], lp->piz[i]);
	}
	if (dz != NULL)
	{
		if (lp->basisstat.optimal == 0)
		{
			EG_RETURN (1);
		}
		for (i = 0; i < lp->nrows; i++)
			EGLPNUM_TYPENAME_EGlpNumZero (dz[lp->baz[i]]);
		for (j = 0; j < lp->nnbasic; j++)
			EGLPNUM_TYPENAME_EGlpNumCopy (dz[lp->nbaz[j]], lp->dz[j]);
	}
	if (objval != NULL)
		EGLPNUM_TYPENAME_EGlpNumCopy (*objval, lp->objval);
	return 0;
}

int EGLPNUM_TYPENAME_ILLsimplex_infcertificate (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPE * pi)
{
	int i, col, nz;
	char *sense;
	EGLPNUM_TYPE *x, *l, *u;
	EGLPNUM_TYPENAME_lp_status_info *ls;

	if (pi == NULL)
		return 0;

	ls = &(lp->basisstat);
	if (ls->primal_infeasible == 0 && ls->dual_unbounded == 0)
	{
		EG_RETURN (1);
	}

	if (lp->final_phase == PRIMAL_PHASEI && lp->pIpiz != NULL)
	{
		for (i = 0; i < lp->nrows; i++)
			EGLPNUM_TYPENAME_EGlpNumCopy (pi[i], lp->pIpiz[i]);
	}
	else if (lp->final_phase == DUAL_PHASEII && lp->infub_ix != -1)
	{
		col = lp->baz[lp->infub_ix];
		x = &(lp->xbz[lp->infub_ix]);
		l = &(lp->lz[col]);
		u = &(lp->uz[col]);

		for (i = 0; i < lp->nrows; i++)
			EGLPNUM_TYPENAME_EGlpNumZero (pi[i]);

		if (EGLPNUM_TYPENAME_EGlpNumIsNeqq (*l, EGLPNUM_TYPENAME_NINFTY) && EGLPNUM_TYPENAME_EGlpNumIsLess (*x, *l))
		{
			for (i = 0, nz = lp->zz.nzcnt; i < nz; i++)
				EGLPNUM_TYPENAME_EGlpNumCopyNeg (pi[lp->zz.indx[i]], lp->zz.coef[i]);
		}
		else
		{
			for (i = 0, nz = lp->zz.nzcnt; i < nz; i++)
				EGLPNUM_TYPENAME_EGlpNumCopy (pi[lp->zz.indx[i]], lp->zz.coef[i]);
		}
	}
	else
	{
		QSlog("Invalid call to inf. certificate routine");
		EG_RETURN (1);
	}

	sense = lp->O->sense;
	for (i = 0; i < lp->nrows; i++)
	{
		if (sense[i] == 'G' && EGLPNUM_TYPENAME_EGlpNumIsLessZero (pi[i]))
			EGLPNUM_TYPENAME_EGlpNumZero (pi[i]);
		if (sense[i] == 'L' && EGLPNUM_TYPENAME_EGlpNumIsGreatZero (pi[i]))
			EGLPNUM_TYPENAME_EGlpNumZero (pi[i]);
	}
	return 0;
}

#if SIMPLEX_DEBUG > 1
static void test_cert (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPE * pi)
{
	int i, j;
	int mcnt, mbeg;
	EGLPNUM_TYPE fsum, sum;

	EGLPNUM_TYPENAME_EGlpNumInitVar (fsum);
	EGLPNUM_TYPENAME_EGlpNumInitVar (sum);
	EGLPNUM_TYPENAME_EGlpNumZero (fsum);

	for (i = 0; i < lp->nrows; i++)
	{
		if (lp->O->sense[i] == 'G' && EGLPNUM_TYPENAME_EGlpNumIsLessZero (pi[i]))
			QSlog("compl");
		if (lp->O->sense[i] == 'L' && EGLPNUM_TYPENAME_EGlpNumIsGreatZero (pi[i]))
			QSlog("compll");
	}

	for (i = 0; i < lp->nrows; i++)
		EGLPNUM_TYPENAME_EGlpNumAddInnProdTo (fsum, pi[i], lp->bz[i]);

	for (j = 0; j < lp->nnbasic; j++)
	{
		EGLPNUM_TYPENAME_EGlpNumZero (sum);
		mcnt = lp->matcnt[j];
		mbeg = lp->matbeg[j];
		for (i = 0; i < mcnt; i++)
			EGLPNUM_TYPENAME_EGlpNumAddInnProdTo (sum, pi[lp->matind[mbeg + i]], lp->matval[mbeg + i]);

		if (EGLPNUM_TYPENAME_EGlpNumIsLess (EGLPNUM_TYPENAME_PFEAS_TOLER, sum) &&
				(lp->vtype[j] == VLOWER || lp->vtype[j] == VFREE))
			QSlog("compl2");
		else
		{
			EGLPNUM_TYPENAME_EGlpNumSign (sum);
			if (EGLPNUM_TYPENAME_EGlpNumIsLess (EGLPNUM_TYPENAME_PFEAS_TOLER, sum) &&
					(lp->vtype[j] == VUPPER || lp->vtype[j] == VFREE))
				QSlog("compl1");
			EGLPNUM_TYPENAME_EGlpNumSign (sum);
		}

		if (EGLPNUM_TYPENAME_EGlpNumIsLessZero (sum)
				&& (lp->vtype[j] & (VFREE | VUPPER)) == 0)
			EGLPNUM_TYPENAME_EGlpNumSubInnProdTo (fsum, sum, lp->lz[j]);
		else if (EGLPNUM_TYPENAME_EGlpNumIsGreatZero (sum)
						 && (lp->vtype[j] & (VFREE | VLOWER)) == 0)
			EGLPNUM_TYPENAME_EGlpNumSubInnProdTo (fsum, sum, lp->uz[j]);
	}
	QSlog("fsum = %.8f", EGLPNUM_TYPENAME_EGlpNumToLf (fsum));
	EGLPNUM_TYPENAME_EGlpNumClearVar (fsum);
	EGLPNUM_TYPENAME_EGlpNumClearVar (sum);
}
#endif

static void save_paraminfo (
	EGLPNUM_TYPENAME_price_info * pinf,
	EGLPNUM_TYPENAME_iter_info * it)
{
	EGLPNUM_TYPENAME_param_info *pr = &(it->oldinfo);

	pr->origalgo = it->algorithm;
	pr->pphaseI = pinf->pI_price;
	pr->pphaseII = pinf->pII_price;
	pr->dphaseI = pinf->dI_price;
	pr->dphaseII = pinf->dII_price;
	pr->p_strategy = pinf->p_strategy;
	pr->d_strategy = pinf->d_strategy;
}

static void restore_paraminfo (
	EGLPNUM_TYPENAME_iter_info * it,
	EGLPNUM_TYPENAME_price_info * pinf)
{
	EGLPNUM_TYPENAME_param_info *pr = &(it->oldinfo);

	it->algorithm = pr->origalgo;
	pinf->pI_price = pr->pphaseI;
	pinf->pII_price = pr->pphaseII;
	pinf->dI_price = pr->dphaseI;
	pinf->dII_price = pr->dphaseII;
	pinf->p_strategy = pr->p_strategy;
	pinf->d_strategy = pr->d_strategy;
}

int EGLPNUM_TYPENAME_ILLsimplex (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int algorithm,
	EGLPNUM_TYPENAME_ILLlp_basis * B,
	EGLPNUM_TYPENAME_price_info * pinf,
	int *status,
	int sdisplay,
	itcnt_t*itcnt)
{
	int phase = -1;
	int singular = -1;
	int rval = 0;
	int new_price = -1;
	EGLPNUM_TYPENAME_svector wz;
	EGLPNUM_TYPENAME_svector updz;
	EGLPNUM_TYPENAME_feas_info fi;
	EGLPNUM_TYPENAME_iter_info it;

	EGLPNUM_TYPENAME_EGlpNumInitVar (fi.totinfeas);
	EGLPNUM_TYPENAME_EGlpNumInitVar (it.prevobj);
	EGLPNUM_TYPENAME_EGlpNumInitVar (it.objtol);

	it.newphase = -1;
	it.nextphase = -1;
	it.nextstep = -1;
	it.sdisplay = sdisplay;
	it.n_pivot_fail = 0;
	it.itercnt = 0;
	it.n_restart = 0;
	it.solstatus = ILL_LP_UNSOLVED;
	it.curtime = 0;
	it.rounds = 0;
	EGLPNUM_TYPENAME_EGlpNumCopy (it.prevobj, EGLPNUM_TYPENAME_INFTY);
	it.nosolve = 0;
	it.noprog = 0;
	EGLPNUM_TYPENAME_EGlpNumCopy (it.objtol, EGLPNUM_TYPENAME_OBJBND_TOLER);
	it.chkobj = PARAM_MAX_NOPROG;
	it.inner = 0;
	it.algorithm = algorithm;
	it.pricetype = -1;
	it.resumeid = -1;
	save_paraminfo (pinf, &it);

#if SIMPLEX_DEBUG > 0
	if (lp->O->nrows > 1000)
		it.sdisplay = 1;
#endif
	if (status)
		*status = QS_LP_UNSOLVED;

	EGLPNUM_TYPENAME_free_internal_lpinfo (lp);
	EGLPNUM_TYPENAME_init_internal_lpinfo (lp);
	rval = EGLPNUM_TYPENAME_build_internal_lpinfo (lp);
	CHECKRVALG (rval, CLEANUP);

	EGLPNUM_TYPENAME_ILLsvector_init (&wz);
	rval = EGLPNUM_TYPENAME_ILLsvector_alloc (&wz, lp->nrows);
	CHECKRVALG (rval, CLEANUP);
	EGLPNUM_TYPENAME_ILLsvector_init (&updz);
	rval = EGLPNUM_TYPENAME_ILLsvector_alloc (&updz, lp->nrows);
	CHECKRVALG (rval, CLEANUP);

	if (it.sdisplay)
	{
		char buffer[256];
		int nonzero = 0;
		register int i = lp->ncols;

		while (i--)
			nonzero += lp->matcnt[i];
		sprintf (buffer, "starting EGLPNUM_TYPENAME_ILLsimplex on %s...", lp->O->probname);
		/* depending on LP's reporter 
		 * string is printed to stdout 
		 * or handed to GUI */
		rval = rval || ILLstring_report (buffer, &lp->O->reporter);
		QSlog("Problem has %d rows and %d cols and %d nonzeros", lp->nrows,
								lp->ncols, nonzero);
	}
	EGLPNUM_TYPENAME_ILLfct_set_variable_type (lp);

	if (B != 0)
	{
		rval = EGLPNUM_TYPENAME_ILLbasis_load (lp, B);
		CHECKRVALG (rval, CLEANUP);
		if (it.algorithm == DUAL_SIMPLEX)
		{
			if (B->rownorms)
			{
				rval = EGLPNUM_TYPENAME_ILLprice_load_rownorms (lp, B->rownorms, pinf);
				CHECKRVALG (rval, CLEANUP);
			}
			else
				EGLPNUM_TYPENAME_EGlpNumFreeArray (pinf->dsinfo.norms);
		}
		else if (it.algorithm == PRIMAL_SIMPLEX)
		{
			if (B->colnorms)
			{
				rval = EGLPNUM_TYPENAME_ILLprice_load_colnorms (lp, B->colnorms, pinf);
				CHECKRVALG (rval, CLEANUP);
			}
			else
				EGLPNUM_TYPENAME_EGlpNumFreeArray (pinf->psinfo.norms);
		}
		else if (it.algorithm != PRIMAL_OR_DUAL)
		{
			QSlog("Unknown algorithm %d in EGLPNUM_TYPENAME_ILLsimplex", it.algorithm);
			rval = 1;
			ILL_CLEANUP;
		}
	}
	else if (lp->basisid == -1)
	{
		if (lp->nrows < 200 && lp->ncols < 400)
			rval = EGLPNUM_TYPENAME_ILLbasis_get_initial (lp, it.algorithm);
		else
			rval = EGLPNUM_TYPENAME_ILLbasis_get_cinitial (lp, it.algorithm);
		CHECKRVALG (rval, CLEANUP);
		EGLPNUM_TYPENAME_ILLprice_free_pricing_info (pinf);
	}

	if (lp->fbasisid != lp->basisid)
	{
		rval = EGLPNUM_TYPENAME_ILLbasis_factor (lp, &singular);
		CHECKRVALG (rval, CLEANUP);
		if (singular)
		{
			MESSAGE (__QS_SB_VERB, "Singular basis found!");
			EGLPNUM_TYPENAME_ILLprice_free_pricing_info (pinf);
		}
	}

START:
#if 0
	if (it.resumeid == SIMPLEX_RESUME_UNSHIFT)
		QSlog("Resuming Unshift");
	else if (it.resumeid == SIMPLEX_RESUME_SING)
		QSlog("Resuming Singular");
	else if (it.resumeid == SIMPLEX_RESUME_NUMER)
		QSlog("Resuming Numer");
	else if (it.resumeid != -1)
		QSlog("Resuming for other reason... %d", it.resumeid);
#endif
	it.solstatus = ILL_LP_UNSOLVED;
	init_lp_status_info (&(lp->basisstat));

	EGLPNUM_TYPENAME_ILLfct_compute_piz (lp);
	EGLPNUM_TYPENAME_ILLfct_compute_dz (lp);
	if (it.algorithm == DUAL_SIMPLEX)
	{
		if (B != NULL || it.resumeid == SIMPLEX_RESUME_UNSHIFT)
			EGLPNUM_TYPENAME_ILLfct_dual_adjust (lp, lp->tol->dfeas_tol);
		else
			EGLPNUM_TYPENAME_ILLfct_dual_adjust (lp, EGLPNUM_TYPENAME_zeroLpNum);
	}
	EGLPNUM_TYPENAME_ILLfct_compute_xbz (lp);

	EGLPNUM_TYPENAME_ILLfct_check_pfeasible (lp, &fi, lp->tol->pfeas_tol);
	EGLPNUM_TYPENAME_ILLfct_check_dfeasible (lp, &fi, lp->tol->dfeas_tol);
	EGLPNUM_TYPENAME_ILLfct_set_status_values (lp, fi.pstatus, fi.dstatus, PHASEII, PHASEII);

	if (fi.dstatus == DUAL_FEASIBLE && EGLPNUM_TYPENAME_EGlpNumIsNeqq (lp->objbound, EGLPNUM_TYPENAME_INFTY))
	{
		EGLPNUM_TYPENAME_ILLfct_compute_dobj (lp);
		if (EGLPNUM_TYPENAME_EGlpNumIsLess (lp->objbound, lp->dobjval))
		{
			it.solstatus = ILL_BND_REACHED;
			QSlog("solstatus = ILL_BND_REACHED = 5 %lf %lf",
									EGLPNUM_TYPENAME_EGlpNumToLf (lp->objbound),
									EGLPNUM_TYPENAME_EGlpNumToLf (lp->dobjval));
			goto TERMINATE;
		}
	}
	if (fi.pstatus == PRIMAL_FEASIBLE && fi.dstatus == DUAL_FEASIBLE)
	{
		it.solstatus = ILL_LP_SOLVED;
		EGLPNUM_TYPENAME_ILLfct_compute_pobj (lp);
		goto TERMINATE;
	}

	if (it.algorithm == PRIMAL_OR_DUAL)
	{
		if (fi.pstatus == PRIMAL_FEASIBLE)
			it.algorithm = PRIMAL_SIMPLEX;
		else if (fi.dstatus == DUAL_FEASIBLE)
			it.algorithm = DUAL_SIMPLEX;
		else if (EGLPNUM_TYPENAME_EGlpNumToLf (lp->pinfeas) < 10 * EGLPNUM_TYPENAME_EGlpNumToLf (lp->dinfeas))
			it.algorithm = PRIMAL_SIMPLEX;
		else
			it.algorithm = DUAL_SIMPLEX;
	}

	if (it.algorithm == PRIMAL_SIMPLEX)
	{
		if (fi.pstatus == PRIMAL_FEASIBLE)
			phase = PRIMAL_PHASEII;
		else
			phase = PRIMAL_PHASEI;
	}
	else if (it.algorithm == DUAL_SIMPLEX)
	{
		if (fi.dstatus == DUAL_FEASIBLE)
			phase = DUAL_PHASEII;
		else
			phase = DUAL_PHASEI;
	}

	rval = EGLPNUM_TYPENAME_ILLprice_build_pricing_info (lp, pinf, phase);
	CHECKRVALG (rval, CLEANUP);

	it.newphase = SIMPLEX_PHASE_NEW;
	it.nextstep = SIMPLEX_CONTINUE;

	while (it.nextstep == SIMPLEX_CONTINUE)
	{
		//ILL_IFTRACE("\t%s:%s:%d\n",__func__,__FILE__,__LINE__);

		if (phase == PRIMAL_PHASEI)
		{
			rval = primal_phaseI_step (lp, pinf, &updz, &wz, &it);
			CHECKRVALG (rval, CLEANUP);
		}

		else if (phase == PRIMAL_PHASEII)
		{
			rval = primal_phaseII_step (lp, pinf, &updz, &wz, &it);
			CHECKRVALG (rval, CLEANUP);
		}

		else if (phase == DUAL_PHASEI)
		{
			rval = dual_phaseI_step (lp, pinf, &updz, &wz, &it);
			CHECKRVALG (rval, CLEANUP);
		}

		else if (phase == DUAL_PHASEII)
		{
			rval = dual_phaseII_step (lp, pinf, &updz, &wz, &it);
			CHECKRVALG (rval, CLEANUP);
		}
		if (it.nextstep == SIMPLEX_RESUME)
		{
			//ILL_IFTRACE("\t%s:%s:%d\n",__func__,__FILE__,__LINE__);
			EGLPNUM_TYPENAME_ILLprice_free_pricing_info (pinf);
			if (it.resumeid == SIMPLEX_RESUME_UNSHIFT)
			{
				if (it.pricetype == QS_PRICE_PDEVEX)
				{
					pinf->pI_price = QS_PRICE_PDEVEX;
					pinf->pII_price = QS_PRICE_PDEVEX;
				}
				else if (it.pricetype == QS_PRICE_DDEVEX)
				{
					pinf->dI_price = QS_PRICE_DDEVEX;
					pinf->dII_price = QS_PRICE_DDEVEX;
				}
			}
			else if (it.resumeid == SIMPLEX_RESUME_NUMER)
			{
				EGLPNUM_TYPENAME_ILLfct_unroll_bound_change (lp);
				EGLPNUM_TYPENAME_ILLfct_unroll_coef_change (lp);
				/* we are disabling re-do under this circunstances ! */
				rval = EGLPNUM_TYPENAME_ILLbasis_get_initial (lp, it.algorithm);
				CHECKRVALG (rval, CLEANUP);
				rval = EGLPNUM_TYPENAME_ILLbasis_factor (lp, &singular);
				if (singular)
					MESSAGE (__QS_SB_VERB, "Singular basis found!");
				CHECKRVALG (rval, CLEANUP);
			}
			it.pricetype = -1;
			if (it.n_restart > SIMPLEX_MAX_RESTART)
			{
				it.solstatus = ILL_MAX_ITER;
				goto LIMIT_TERMINATE;
			}
			goto START;
		}
		else if (it.nextstep == SIMPLEX_CONTINUE)
		{
			//ILL_IFTRACE("\t%s:%s:%d\n",__func__,__FILE__,__LINE__);
			it.itercnt++;

			if (it.nextphase != phase)
			{
				//ILL_IFTRACE("\t%s:%s:%d\n",__func__,__FILE__,__LINE__);
				it.newphase = SIMPLEX_PHASE_NEW;
				phase = it.nextphase;
				new_price = EGLPNUM_TYPENAME_ILLprice_get_price (pinf, phase);

				if (pinf->cur_price != new_price)
				{
					//ILL_IFTRACE("\t%s:%s:%d\n",__func__,__FILE__,__LINE__);
					EGLPNUM_TYPENAME_ILLprice_free_pricing_info (pinf);
					rval = EGLPNUM_TYPENAME_ILLprice_build_pricing_info (lp, pinf, phase);
					CHECKRVALG (rval, CLEANUP);
				}
			}
		}
	}

#if SIMPLEX_DEBUG > 0
	EGLPNUM_TYPENAME_ILLfct_print_counts (lp);
#endif
LIMIT_TERMINATE:
	rval = terminate_simplex (lp, phase, &it);
	CHECKRVALG (rval, CLEANUP);

TERMINATE:
	restore_paraminfo (&it, pinf);

	if (it.sdisplay)
	{
		QSlog("completed EGLPNUM_TYPENAME_ILLsimplex");
		QSlog("%s: ", lp->O->probname);
	}

	if (status)
	{
		if (it.solstatus == ILL_MAX_ITER)
		{
			*status = QS_LP_ITER_LIMIT;
		}
		else if (it.solstatus == ILL_MAX_TIME)
		{
			*status = QS_LP_TIME_LIMIT;
		}
		else if (it.solstatus == ILL_LP_ABORTED)
		{
			*status = QS_LP_ABORTED;
		}
		else if (it.solstatus == ILL_PPHASEI_ERROR ||
						 it.solstatus == ILL_PPHASEII_ERROR ||
						 it.solstatus == ILL_DPHASEI_ERROR ||
						 it.solstatus == ILL_DPHASEII_ERROR)
		{
			*status = QS_LP_NUMERR;
		}
		else if(it.solstatus == ILL_LP_UNSOLVED)
		{
			*status = QS_LP_UNSOLVED;
		}
		else if (it.solstatus == ILL_BND_REACHED)
		{
			*status = QS_LP_OBJ_LIMIT;
		}
		else if (it.solstatus == ILL_LP_SOLVED)
		{
			if (lp->basisstat.optimal)
			{
				*status = QS_LP_OPTIMAL;
			}
			else if (lp->basisstat.primal_infeasible || lp->basisstat.dual_unbounded)
			{
				*status = QS_LP_INFEASIBLE;
				if (it.sdisplay)
				{
					if (lp->basisstat.primal_infeasible)
						QSlog("Primal Infeasible");
					else
						QSlog("Dual Unbounded");
				}
			}
			else if (lp->basisstat.primal_unbounded)
			{
				*status = QS_LP_UNBOUNDED;
			}
		}
		else
		{
			QSlog("unknown solution status in EGLPNUM_TYPENAME_ILLsimplex %d",
									it.solstatus);
			rval = 1;
			CHECKRVALG (rval, CLEANUP);
		}
	}

#if SIMPLEX_DEBUG > 1
	{
		int rva = 0;
		EGLPNUM_TYPE *pi = NULL;

		pi = EGLPNUM_TYPENAME_EGlpNumAllocArray (lp->nrows);
		rva = EGLPNUM_TYPENAME_ILLsimplex_infcertificate (lp, pi);
		QSlog("rva = %d", rva);
		if (!rva)
		{
			test_cert (lp, pi);
		}
		EGLPNUM_TYPENAME_EGlpNumFreeArray (pi);
	}
#endif
	 /* update counter */
	 itcnt->pI_iter += lp->cnts->pI_iter;
	 itcnt->pII_iter += lp->cnts->pII_iter;
	 itcnt->dI_iter += lp->cnts->dI_iter;
	 itcnt->dII_iter += lp->cnts->dII_iter;
	 itcnt->tot_iter = itcnt->pI_iter + itcnt->pII_iter + itcnt->dI_iter +
	 										itcnt->dII_iter;
   /* end update */
	if (it.sdisplay)
	{
		int bstat = 0;

		QSlog("time = %.3f, pI = %d, pII = %d, dI = %d, dII = %d,",
								ILLutil_zeit () - lp->starttime, lp->cnts->pI_iter,
								lp->cnts->pII_iter, lp->cnts->dI_iter, lp->cnts->dII_iter);
		get_current_stat (&(lp->basisstat), it.algorithm, &bstat);
		switch (bstat)
		{
		case OPTIMAL:
			QSlog("opt = %f", EGLPNUM_TYPENAME_EGlpNumToLf (lp->objval));
			break;
		case PRIMAL_INFEASIBLE:
			QSlog("no primal soln");
			break;
		case PRIMAL_UNBOUNDED:
			QSlog("primal unbounded");
			break;
		case PRIMAL_FEASIBLE:
			QSlog("primal obj = %f", EGLPNUM_TYPENAME_EGlpNumToLf (lp->pobjval));
			break;
		case DUAL_INFEASIBLE:
			QSlog("no dual soln");
			break;
		case DUAL_UNBOUNDED:
			QSlog("dual unbounded");
			break;
		case DUAL_FEASIBLE:
			QSlog("dual obj = %f", EGLPNUM_TYPENAME_EGlpNumToLf (lp->dobjval));
			break;
		}

		if (it.sdisplay > 1)
		{
			if (it.algorithm == PRIMAL_SIMPLEX && pinf->pI_price == QS_PRICE_PDEVEX)
				QSlog("Devex norms initialised %d times", pinf->pdinfo.ninit);
		}
	}

CLEANUP:
	EGLPNUM_TYPENAME_ILLsvector_free (&wz);
	EGLPNUM_TYPENAME_ILLsvector_free (&updz);
	EGLPNUM_TYPENAME_EGlpNumClearVar (it.prevobj);
	EGLPNUM_TYPENAME_EGlpNumClearVar (it.objtol);
	EGLPNUM_TYPENAME_EGlpNumClearVar (fi.totinfeas);
	if (rval == QS_LP_CHANGE_PREC)
	{
		MESSAGE (__QS_SB_VERB, "Changing precision");
		return rval;
	}
	else
	{
		MESSAGE (rval ? 0 : 1000, "Error code %d", rval);
		EG_RETURN (rval);
	}
}

static int terminate_simplex (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int phase,
	EGLPNUM_TYPENAME_iter_info * it)
{
	int rval = 0;
	int sphase;
	EGLPNUM_TYPENAME_feas_info fi;

	EGLPNUM_TYPENAME_EGlpNumInitVar (fi.totinfeas);

	if (it->solstatus != ILL_MAX_TIME && it->solstatus != ILL_MAX_ITER)
		ILL_CLEANUP;

	if (it->algorithm == PRIMAL_SIMPLEX)
	{
		if (lp->nbchange != 0)
		{
			if (it->sdisplay > 1)
			{
				QSlog("unrolling %d bound shifts", lp->nbchange);
			}
			EGLPNUM_TYPENAME_ILLfct_unroll_bound_change (lp);
		}
		rval = EGLPNUM_TYPENAME_ILLsimplex_retest_psolution (lp, NULL, phase, &fi);
		CHECKRVALG (rval, CLEANUP);

		sphase = (phase == PRIMAL_PHASEI) ? PHASEI : PHASEII;
		EGLPNUM_TYPENAME_ILLfct_set_status_values (lp, fi.pstatus, fi.dstatus, PHASEII, sphase);
	}
	else if (it->algorithm == DUAL_SIMPLEX)
	{
		if (lp->ncchange != 0)
		{
			if (it->sdisplay > 1)
			{
				QSlog("unrolling %d coef shifts", lp->ncchange);
			}
			EGLPNUM_TYPENAME_ILLfct_unroll_coef_change (lp);
		}
		rval = EGLPNUM_TYPENAME_ILLsimplex_retest_dsolution (lp, NULL, phase, &fi);
		CHECKRVALG (rval, CLEANUP);

		sphase = (phase == DUAL_PHASEI) ? PHASEI : PHASEII;
		EGLPNUM_TYPENAME_ILLfct_set_status_values (lp, fi.pstatus, fi.dstatus, sphase, PHASEII);
	}

CLEANUP:
	EGLPNUM_TYPENAME_EGlpNumClearVar (fi.totinfeas);
	EG_RETURN (rval);
}

static int test_progress (
	EGLPNUM_TYPE objval,
	EGLPNUM_TYPE prevobj)
{
	EGLPNUM_TYPE denom;

	EGLPNUM_TYPENAME_EGlpNumInitVar (denom);
	EGLPNUM_TYPENAME_EGlpNumCopyDiff (denom, objval, prevobj);
	if (EGLPNUM_TYPENAME_EGlpNumIsNeqZero (objval, EGLPNUM_TYPENAME_PROGRESS_ZERO))
		EGLPNUM_TYPENAME_EGlpNumDivTo (denom, objval);
	if (!EGLPNUM_TYPENAME_EGlpNumIsNeqZero (denom, EGLPNUM_TYPENAME_PROGRESS_THRESH))
	{
		EGLPNUM_TYPENAME_EGlpNumClearVar (denom);
		return 0;
	}
	else
	{
		EGLPNUM_TYPENAME_EGlpNumClearVar (denom);
		return 1;
	}
}

static void monitor_iter (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_price_info * p,
	EGLPNUM_TYPENAME_iter_info * it,
	int phase)
{
	EGLPNUM_TYPE print_val;
	double tottime = ILLutil_zeit () - lp->starttime;
	int curtime = ILLutil_our_floor (tottime);	/* MONIKA */
	char print_str[20];
	EGLPNUM_TYPENAME_feas_info fi;
	int aborted = 0;

	EGLPNUM_TYPENAME_EGlpNumInitVar (print_val);
	EGLPNUM_TYPENAME_EGlpNumInitVar (fi.totinfeas);
	EGLPNUM_TYPENAME_EGlpNumZero (fi.totinfeas);
	EGLPNUM_TYPENAME_EGlpNumZero (print_val);

	/* one of the following two time display mechanisms */
	switch (phase)
	{
	case PRIMAL_PHASEI:
		EGLPNUM_TYPENAME_EGlpNumCopy (print_val, lp->pinfeas);
		EGLPNUM_TYPENAME_EGlpNumAddTo (fi.totinfeas, lp->pinfeas);
		strcpy (print_str, "primal infeas");
		if (EGLPNUM_TYPENAME_EGlpNumIsLessZero (lp->pinfeas) &&
				(EGLPNUM_TYPENAME_EGlpNumIsNeqZero (lp->pinfeas, EGLPNUM_TYPENAME_oneLpNum)))
		{
			/*QSlog("Negative Infeasibility! Imposible %lg %la, iter %d",
			 * EGLPNUM_TYPENAME_EGlpNumToLf (print_val), EGLPNUM_TYPENAME_EGlpNumToLf (print_val), it->itercnt);
			 */
			//exit(1);
		}
		break;
	case PRIMAL_PHASEII:
		EGLPNUM_TYPENAME_EGlpNumCopy (print_val, lp->pobjval);
		strcpy (print_str, "primal objval");
		break;
	case DUAL_PHASEI:
		EGLPNUM_TYPENAME_EGlpNumCopy (print_val, lp->dinfeas);
		EGLPNUM_TYPENAME_EGlpNumAddTo (fi.totinfeas, lp->dinfeas);
		strcpy (print_str, "dual infeas");
		break;
	case DUAL_PHASEII:
		EGLPNUM_TYPENAME_EGlpNumCopy (print_val, lp->dobjval);
		strcpy (print_str, "dual objval");
		break;
	}

	aborted = report_value (lp, it, print_str, print_val);
	/*if (it->sdisplay && it->itercnt % lp->iterskip == 0) {
	 * // QSlog("(%d): %s = %f", it->itercnt, print_str, print_val);
	 * } */
	if (curtime != it->curtime)
	{
		it->curtime = curtime;
		/*
		 * if (it->sdisplay){
		 * QSlog("time = %d.0, ", curtime);
		 * QSlog("(%d): %s = %f", it->itercnt, print_str, print_val);
		 * }
		 */
	}

	EGLPNUM_TYPENAME_EGlpNumAddUiTo (fi.totinfeas, 1000);
	if (EGLPNUM_TYPENAME_EGlpNumIsLessZero (fi.totinfeas))
	{
		it->nextstep = SIMPLEX_TERMINATE;
		it->solstatus = ILL_MAX_ITER;
		MESSAGE (it->sdisplay ? 0 : __QS_SB_VERB,
						 "early finish by excess infeasibility");
		ILL_CLEANUP;
	}

	if (phase == DUAL_PHASEII && EGLPNUM_TYPENAME_EGlpNumIsNeqq (lp->objbound, EGLPNUM_TYPENAME_INFTY))
	{
		/*if (lp->dobjval > lp->objbound + it->objtol) */
		EGLPNUM_TYPENAME_EGlpNumCopyDiff (print_val, lp->dobjval, lp->objbound);
		if (EGLPNUM_TYPENAME_EGlpNumIsLess (it->objtol, print_val))
		{
			EGLPNUM_TYPENAME_ILLfct_unroll_coef_change (lp);
			EGLPNUM_TYPENAME_ILLfct_check_dfeasible (lp, &fi, lp->tol->dfeas_tol);
			EGLPNUM_TYPENAME_ILLfct_set_status_values (lp, -1, fi.dstatus, -1, PHASEII);

			if (fi.dstatus == DUAL_FEASIBLE)
			{
				EGLPNUM_TYPENAME_ILLfct_compute_dobj (lp);
				if (EGLPNUM_TYPENAME_EGlpNumIsLess (lp->objbound, lp->dobjval))
				{
					it->solstatus = ILL_BND_REACHED;
					it->nextstep = SIMPLEX_TERMINATE;
					/*if (it->sdisplay) */
					{
						QSlog("bound reached %lf %lf", EGLPNUM_TYPENAME_EGlpNumToLf (lp->objbound),
										EGLPNUM_TYPENAME_EGlpNumToLf (lp->dobjval));
					}
				}
				else
					EGLPNUM_TYPENAME_EGlpNumMultUiTo (it->objtol, 10);
			}
			else
			{
				it->nextphase = DUAL_PHASEI;
				it->newphase = SIMPLEX_PHASE_NEW;
				EGLPNUM_TYPENAME_EGlpNumMultUiTo (it->objtol, 5);
			}
		}
	}
	if (it->itercnt >= lp->maxiter)
	{
		it->solstatus = ILL_MAX_ITER;
		it->nextstep = SIMPLEX_TERMINATE;
		if (it->sdisplay)
		{
			QSlog("iter limit reached");
		}
		ILL_CLEANUP;
	}
	else if (tottime >= lp->maxtime)
	{
		it->solstatus = ILL_MAX_TIME;
		it->nextstep = SIMPLEX_TERMINATE;
		if (it->sdisplay)
		{
			QSlog("time limit reached");
		}
		ILL_CLEANUP;
	}
	else if (aborted)
	{
		it->solstatus = ILL_LP_ABORTED;
		it->nextstep = SIMPLEX_TERMINATE;
		if (it->sdisplay)
		{
			QSlog("aborted");
		}
		ILL_CLEANUP;
	}
	/* why is this commented out? */
	if(0){
		if (it->rounds && it->inner){
			it->inner --;
			if (it->inner == 0){
				QSlog("restoring ..");
				restore_paraminfo (it, p);
				it->newphase   = SIMPLEX_PHASE_NEW;
				it->nextstep   = SIMPLEX_RESUME;
				/*it->resumeid   = SIMPLEX_RESUME_OUTER;*/
				ILL_CLEANUP;
			}
		}
	}
	if (phase == DUAL_PHASEII)
	{
		if (it->noprog > it->chkobj)
		{
			EGLPNUM_TYPENAME_ILLfct_perturb_coefs (lp);
			it->noprog = 0;
			EGLPNUM_TYPENAME_EGlpNumCopy (it->prevobj, lp->dobjval);
		}
	}
	else if (phase == PRIMAL_PHASEII)
	{
		if (it->noprog > it->chkobj)
		{
			EGLPNUM_TYPENAME_ILLfct_perturb_bounds (lp);
			it->noprog = 0;
			EGLPNUM_TYPENAME_EGlpNumCopy (it->prevobj, lp->pobjval);
		}
	}
	else if (phase == PRIMAL_PHASEI)
	{
		if (it->noprog > it->chkobj)
		{
			it->algorithm = DUAL_SIMPLEX;
			it->nextstep = SIMPLEX_RESUME;
			it->resumeid = SIMPLEX_RESUME_NUMER;
			/* this is to force to exit in the case of bad basis */
			EGLPNUM_TYPENAME_EGlpNumMultUiTo (lp->tol->pfeas_tol, SIMPLEX_FACTOR);
			EGLPNUM_TYPENAME_EGlpNumMultUiTo (lp->tol->dfeas_tol, SIMPLEX_FACTOR);
			it->n_restart++;
			//QSlog("Resume Numerical %s:%s:%d",__func__,__FILE__,__LINE__);
		}
	}
	else if (phase == DUAL_PHASEI)
	{
		if (it->noprog > it->chkobj)
		{
			it->algorithm = PRIMAL_SIMPLEX;
			it->nextstep = SIMPLEX_RESUME;
			it->resumeid = SIMPLEX_RESUME_NUMER;
			/* this is to force to exit in the case of bad basis */
			EGLPNUM_TYPENAME_EGlpNumMultUiTo (lp->tol->pfeas_tol, SIMPLEX_FACTOR);
			EGLPNUM_TYPENAME_EGlpNumMultUiTo (lp->tol->dfeas_tol, SIMPLEX_FACTOR);
			it->n_restart++;
			//QSlog("Resume Numerical %s:%s:%d",__func__,__FILE__,__LINE__);
		}
	}
CLEANUP:
	EGLPNUM_TYPENAME_EGlpNumClearVar (fi.totinfeas);
	EGLPNUM_TYPENAME_EGlpNumClearVar (print_val);
	return;
}

static int primal_phaseI_step (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_price_info * pinf,
	EGLPNUM_TYPENAME_svector * updz,
	EGLPNUM_TYPENAME_svector * wz,
	EGLPNUM_TYPENAME_iter_info * it)
{
	int rval = 0;
	int singular = 0;
	int refactor = 0;
	int cphase = PRIMAL_PHASEI;
	EGLPNUM_TYPE alpha;
	EGLPNUM_TYPENAME_feas_info fi;
	EGLPNUM_TYPENAME_ratio_res rs;
	EGLPNUM_TYPENAME_price_res pr;

	EGLPNUM_TYPENAME_EGlpNumInitVar (alpha);
	EGLPNUM_TYPENAME_EGlpNumInitVar (fi.totinfeas);
	EGLPNUM_TYPENAME_EGlpNumInitVar (pr.dinfeas);
	EGLPNUM_TYPENAME_EGlpNumInitVar (pr.pinfeas);
	EGLPNUM_TYPENAME_EGlpNumInitVar (rs.tz);
	EGLPNUM_TYPENAME_EGlpNumInitVar (rs.lbound);
	EGLPNUM_TYPENAME_EGlpNumInitVar (rs.ecoeff);
	EGLPNUM_TYPENAME_EGlpNumInitVar (rs.pivotval);
	EGLPNUM_TYPENAME_EGlpNumZero (alpha);

	EGLPNUM_TYPENAME_ILLfct_update_counts (lp, CNT_PPHASE1ITER, 0, EGLPNUM_TYPENAME_zeroLpNum);
	it->nextstep = SIMPLEX_CONTINUE;
	it->nextphase = PRIMAL_PHASEI;
	lp->final_phase = PRIMAL_PHASEI;
	it->nosolve++;

	if (it->newphase != 0)
	{
		EGLPNUM_TYPENAME_ILLfct_check_pfeasible (lp, &fi, lp->tol->ip_tol);
		if (it->newphase == SIMPLEX_PHASE_NEW)
		{
			it->noprog = 0;
			if (it->sdisplay)
			{
				QSlog("starting primal phase I, nosolve %d", it->nosolve);
			}
		}
		it->newphase = 0;
		it->nosolve = 0;
		EGLPNUM_TYPENAME_EGlpNumCopy (it->prevobj, lp->pinfeas);
		lp->pIpiz = EGLPNUM_TYPENAME_EGlpNumAllocArray (lp->nrows);
		lp->pIdz = EGLPNUM_TYPENAME_EGlpNumAllocArray (lp->nnbasic);

		EGLPNUM_TYPENAME_ILLfct_compute_phaseI_piz (lp);
		if (pinf->p_strategy == COMPLETE_PRICING)
		{
			EGLPNUM_TYPENAME_ILLfct_compute_phaseI_dz (lp);
#if USEHEAP > 0
			EGLPNUM_TYPENAME_ILLprice_free_heap (pinf);
#endif
			EGLPNUM_TYPENAME_ILLprice_compute_dual_inf (lp, pinf, NULL, 0, PRIMAL_PHASEI);
#if USEHEAP > 0
			rval = EGLPNUM_TYPENAME_ILLprice_test_for_heap (lp, pinf, lp->nnbasic, pinf->d_scaleinf,
																		 PRIMAL_SIMPLEX, 0);
			CHECKRVALG (rval, CLEANUP);
#endif
		}
		else if (pinf->p_strategy == MULTI_PART_PRICING)
			EGLPNUM_TYPENAME_ILLprice_init_mpartial_price (lp, pinf, cphase, COL_PRICING);
	}

	monitor_iter (lp, pinf, it, cphase);
	if (it->nextstep == SIMPLEX_TERMINATE || it->nextstep == SIMPLEX_RESUME ||
			it->newphase != 0)
		ILL_CLEANUP;

	EGLPNUM_TYPENAME_ILLprice_primal (lp, pinf, &pr, cphase);
	ILL_IFTRACE2 ("%s:after_price\n", __func__);

	if (pr.price_stat == PRICE_OPTIMAL)
	{
		if (it->sdisplay > 1)
		{
			QSlog("primal phase I seemingly done");
			QSlog("retesting soln");
		}
		rval = EGLPNUM_TYPENAME_ILLsimplex_retest_psolution (lp, pinf, cphase, &fi);

		CHECKRVALG (rval, CLEANUP);
		EGLPNUM_TYPENAME_ILLfct_set_status_values (lp, fi.pstatus, fi.dstatus, PHASEII, PHASEI);

		if (fi.pstatus == PRIMAL_FEASIBLE)
		{
			it->nextphase = PRIMAL_PHASEII;
		}
		else if (fi.dstatus == DUAL_FEASIBLE)
		{
			it->solstatus = ILL_LP_SOLVED;
			it->nextstep = SIMPLEX_TERMINATE;
		}
		ILL_CLEANUP;
	}

	EGLPNUM_TYPENAME_ILLfct_compute_yz (lp, &(lp->yjz), updz, lp->nbaz[pr.eindex]);
	EGLPNUM_TYPENAME_ILLfct_update_counts (lp, CNT_YNZ, lp->yjz.nzcnt, EGLPNUM_TYPENAME_zeroLpNum);
	EGLPNUM_TYPENAME_ILLfct_update_counts (lp, CNT_UPNZ, updz->nzcnt, EGLPNUM_TYPENAME_zeroLpNum);

	EGLPNUM_TYPENAME_ILLratio_pI_test (lp, pr.eindex, pr.dir, &rs);
	//ILL_IFTRACE(":%d",rs.lindex);

	if (rs.ratio_stat == RATIO_FAILED)
	{
		/*
		 * rval = E_SIMPLEX_ERROR;
		 * it->solstatus = ILL_PPHASEI_ERROR;
		 */
		//ILL_IFTRACE("ratio_failed\n");
		it->algorithm = DUAL_SIMPLEX;
		it->nextstep = SIMPLEX_RESUME;
		it->resumeid = SIMPLEX_RESUME_NUMER;
		/* this is to force to exit in the case of bad basis */
		it->n_restart++;
		EGLPNUM_TYPENAME_EGlpNumMultUiTo (lp->tol->pfeas_tol, SIMPLEX_FACTOR);
		EGLPNUM_TYPENAME_EGlpNumMultUiTo (lp->tol->dfeas_tol, SIMPLEX_FACTOR);
		//QSlog("Resume Numerical %s:%s:%d",__func__,__FILE__,__LINE__);
		ILL_CLEANUP;
	}
	else if (rs.ratio_stat == RATIO_NEGATIVE)
	{
		EGLPNUM_TYPE itol;

		EGLPNUM_TYPENAME_EGlpNumInitVar (itol);
		//ILL_IFTRACE("ratio_negative\n");
		EGLPNUM_TYPENAME_EGlpNumCopy (itol, lp->tol->ip_tol);
		EGLPNUM_TYPENAME_EGlpNumZero (lp->tol->ip_tol);
		EGLPNUM_TYPENAME_EGlpNumAddTo (lp->pinfeas, lp->upd.c_obj);
		if (!test_progress (lp->pinfeas, it->prevobj))
			it->noprog++;
		else
		{
			EGLPNUM_TYPENAME_EGlpNumCopy (it->prevobj, lp->pinfeas);
			it->noprog = 0;
		}
		EGLPNUM_TYPENAME_ILLfct_update_pfeas (lp, rs.lindex, &(lp->srhs));
		EGLPNUM_TYPENAME_EGlpNumCopy (lp->tol->ip_tol, itol);
		EGLPNUM_TYPENAME_ILLfct_compute_ppIzz (lp, &(lp->srhs), &(lp->ssoln));
		EGLPNUM_TYPENAME_ILLfct_update_ppI_prices (lp, pinf, &(lp->srhs), &(lp->ssoln), pr.eindex,
															rs.lindex, EGLPNUM_TYPENAME_zeroLpNum);
		EGLPNUM_TYPENAME_EGlpNumClearVar (itol);
	}
	else if (rs.ratio_stat == RATIO_NOBCHANGE)
	{
		//ILL_IFTRACE("ratio_nobchange\n");
		EGLPNUM_TYPENAME_EGlpNumAddTo (lp->pinfeas, lp->upd.c_obj);
		if (!test_progress (lp->pinfeas, it->prevobj))
			it->noprog++;
		else
		{
			EGLPNUM_TYPENAME_EGlpNumCopy (it->prevobj, lp->pinfeas);
			it->noprog = 0;
		}

		//ILL_IFTRACE("%s:a\n",__func__);
		EGLPNUM_TYPENAME_ILLfct_update_xz (lp, rs.tz, pr.eindex, rs.lindex);
		EGLPNUM_TYPENAME_ILLfct_update_pfeas (lp, rs.lindex, &(lp->srhs));
		EGLPNUM_TYPENAME_ILLfct_compute_ppIzz (lp, &(lp->srhs), &(lp->ssoln));
		EGLPNUM_TYPENAME_ILLfct_update_basis_info (lp, pr.eindex, rs.lindex, rs.lvstat);
#if DENSE_PI > 0
		EGLPNUM_TYPENAME_fct_test_workvector (lp);
		EGLPNUM_TYPENAME_fct_test_pfeasible (lp);
#endif
		EGLPNUM_TYPENAME_ILLfct_update_ppI_prices (lp, pinf, &(lp->srhs), &(lp->ssoln), pr.eindex,
															rs.lindex, EGLPNUM_TYPENAME_zeroLpNum);
	}
	else if (rs.ratio_stat == RATIO_BCHANGE)
	{
		//ILL_IFTRACE("ratio_bchange\n");
		EGLPNUM_TYPENAME_EGlpNumCopyFrac (alpha, lp->pIdz[pr.eindex], rs.pivotval);
		EGLPNUM_TYPENAME_EGlpNumAddTo (lp->pinfeas, lp->upd.c_obj);

		if (!test_progress (lp->pinfeas, it->prevobj))
		{
			if (lp->vtype[lp->nbaz[pr.eindex]] == VFREE ||
					lp->vtype[lp->baz[rs.lindex]] == VARTIFICIAL)
			{
				if (it->noprog > 0)
					it->noprog--;
			}
			else
				it->noprog++;
		}
		else
		{
			EGLPNUM_TYPENAME_EGlpNumCopy (it->prevobj, lp->pinfeas);
			it->noprog = 0;
		}

		EGLPNUM_TYPENAME_ILLfct_compute_zz (lp, &(lp->zz), rs.lindex);
		EGLPNUM_TYPENAME_ILLfct_update_counts (lp, CNT_ZNZ, lp->zz.nzcnt, EGLPNUM_TYPENAME_zeroLpNum);
		if (pinf->p_strategy == COMPLETE_PRICING)
		{
			//ILL_IFTRACE("%s:a\n",__func__);
			EGLPNUM_TYPENAME_ILLfct_compute_zA (lp, &(lp->zz), &(lp->zA));
			EGLPNUM_TYPENAME_ILLfct_update_counts (lp, CNT_ZANZ, lp->zA.nzcnt, EGLPNUM_TYPENAME_zeroLpNum);

			if (pinf->pI_price == QS_PRICE_PSTEEP)
			{
				EGLPNUM_TYPENAME_ILLfct_compute_psteep_upv (lp, wz);
			}
		}

		rval =
			EGLPNUM_TYPENAME_ILLprice_update_pricing_info (lp, pinf, cphase, wz, pr.eindex, rs.lindex,
																		rs.pivotval);
		CHECKRVALG (rval, CLEANUP);

		//ILL_IFTRACE("%s:b:%d\n",__func__,rs.lindex);
		EGLPNUM_TYPENAME_ILLfct_update_xz (lp, rs.tz, pr.eindex, rs.lindex);
		EGLPNUM_TYPENAME_ILLfct_update_pfeas (lp, rs.lindex, &(lp->srhs));
		//ILL_IFTRACE("%s:%d:%d\n",__func__,rs.lindex,lp->srhs.nzcnt);
		EGLPNUM_TYPENAME_ILLfct_compute_ppIzz (lp, &(lp->srhs), &(lp->ssoln));
		EGLPNUM_TYPENAME_ILLfct_update_basis_info (lp, pr.eindex, rs.lindex, rs.lvstat);
#if DENSE_PI > 0
		EGLPNUM_TYPENAME_fct_test_workvector (lp);
		EGLPNUM_TYPENAME_fct_test_pfeasible (lp);
#endif
		rval = EGLPNUM_TYPENAME_ILLbasis_update (lp, updz, rs.lindex, &refactor, &singular);
		CHECKRVALG (rval, CLEANUP);

		if (singular)
		{
			it->nextstep = SIMPLEX_RESUME;
			it->resumeid = SIMPLEX_RESUME_SING;
			/* this is to force to exit in the case of bad basis */
			it->n_restart++;
			EGLPNUM_TYPENAME_EGlpNumMultUiTo (lp->tol->pfeas_tol, SIMPLEX_FACTOR);
			EGLPNUM_TYPENAME_EGlpNumMultUiTo (lp->tol->dfeas_tol, SIMPLEX_FACTOR);
			//QSlog("Resume Singular %s:%s:%d",__func__,__FILE__,__LINE__);
			ILL_CLEANUP;
		}
		if (!refactor)
		{
			EGLPNUM_TYPENAME_ILLfct_update_ppI_prices (lp, pinf, &(lp->srhs), &(lp->ssoln), pr.eindex,
																rs.lindex, alpha);
		}
		if (refactor != 0 || it->nosolve > PARAM_MAX_NOSOLVE)
		{
			EGLPNUM_TYPENAME_ILLfct_compute_xbz (lp);
			EGLPNUM_TYPENAME_ILLfct_check_pfeasible (lp, &fi, lp->tol->ip_tol);
			EGLPNUM_TYPENAME_ILLfct_set_status_values (lp, fi.pstatus, -1, PHASEII, -1);
			if (fi.pstatus == PRIMAL_FEASIBLE)
				it->nextphase = PRIMAL_PHASEII;

			it->newphase = SIMPLEX_PHASE_RECOMP;
			ILL_CLEANUP;
		}
	}

#if DENSE_PI > 1
	EGLPNUM_TYPENAME_fct_test_workvector (lp);
	fct_test_pi_dz (lp, pinf);
#endif

CLEANUP:
	if (it->nextphase != PRIMAL_PHASEI || it->nextstep == SIMPLEX_RESUME ||
			it->newphase != 0 || rval != 0)
	{
		EGLPNUM_TYPENAME_EGlpNumFreeArray (lp->pIpiz);
		EGLPNUM_TYPENAME_EGlpNumFreeArray (lp->pIdz);
	}
	EGLPNUM_TYPENAME_EGlpNumClearVar (alpha);
	EGLPNUM_TYPENAME_EGlpNumClearVar (fi.totinfeas);
	EGLPNUM_TYPENAME_EGlpNumClearVar (pr.dinfeas);
	EGLPNUM_TYPENAME_EGlpNumClearVar (pr.pinfeas);
	EGLPNUM_TYPENAME_EGlpNumClearVar (rs.tz);
	EGLPNUM_TYPENAME_EGlpNumClearVar (rs.lbound);
	EGLPNUM_TYPENAME_EGlpNumClearVar (rs.ecoeff);
	EGLPNUM_TYPENAME_EGlpNumClearVar (rs.pivotval);
	//EG_RETURN(rval);
	return rval;
}

static int primal_phaseII_step (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_price_info * pinf,
	EGLPNUM_TYPENAME_svector * updz,
	EGLPNUM_TYPENAME_svector * wz,
	EGLPNUM_TYPENAME_iter_info * it)
{
	int boundch;
	int rval = 0;
	int bndtype = 0;
	int singular = 0;
	int refactor = 0;
	int ratio_iter = 0;
	int cphase = PRIMAL_PHASEII;
	EGLPNUM_TYPE lbound;
	EGLPNUM_TYPE alpha;
	EGLPNUM_TYPENAME_feas_info fi;
	EGLPNUM_TYPENAME_ratio_res rs;
	EGLPNUM_TYPENAME_price_res pr;

	EGLPNUM_TYPENAME_EGlpNumInitVar (alpha);
	EGLPNUM_TYPENAME_EGlpNumInitVar (lbound);
	EGLPNUM_TYPENAME_EGlpNumInitVar (fi.totinfeas);
	EGLPNUM_TYPENAME_EGlpNumInitVar (pr.dinfeas);
	EGLPNUM_TYPENAME_EGlpNumInitVar (pr.pinfeas);
	EGLPNUM_TYPENAME_EGlpNumInitVar (rs.tz);
	EGLPNUM_TYPENAME_EGlpNumInitVar (rs.lbound);
	EGLPNUM_TYPENAME_EGlpNumInitVar (rs.ecoeff);
	EGLPNUM_TYPENAME_EGlpNumInitVar (rs.pivotval);

	EGLPNUM_TYPENAME_ILLfct_update_counts (lp, CNT_PPHASE2ITER, 0, EGLPNUM_TYPENAME_zeroLpNum);
	it->nextstep = SIMPLEX_CONTINUE;
	it->nextphase = PRIMAL_PHASEII;
	lp->final_phase = PRIMAL_PHASEII;
	it->nosolve++;

	if (it->newphase != 0)
	{
		EGLPNUM_TYPENAME_ILLfct_compute_pobj (lp);
		if (it->newphase == SIMPLEX_PHASE_NEW)
		{
			it->noprog = 0;
			if (it->sdisplay)
			{
				QSlog("starting primal phase II, nosolve %d", it->nosolve);
			}
		}
		it->newphase = 0;
		it->nosolve = 0;
		EGLPNUM_TYPENAME_EGlpNumCopy (it->prevobj, lp->pobjval);
		EGLPNUM_TYPENAME_ILLfct_compute_piz (lp);
		if (pinf->p_strategy == COMPLETE_PRICING)
		{
			EGLPNUM_TYPENAME_ILLfct_compute_dz (lp);
#if USEHEAP > 0
			EGLPNUM_TYPENAME_ILLprice_free_heap (pinf);
#endif
			EGLPNUM_TYPENAME_ILLprice_compute_dual_inf (lp, pinf, NULL, 0, PRIMAL_PHASEII);
#if USEHEAP > 0
			rval = EGLPNUM_TYPENAME_ILLprice_test_for_heap (lp, pinf, lp->nnbasic, pinf->d_scaleinf,
																		 PRIMAL_SIMPLEX, 0);
			CHECKRVALG (rval, CLEANUP);
#endif
		}
		else if (pinf->p_strategy == MULTI_PART_PRICING)
		{
			EGLPNUM_TYPENAME_ILLprice_init_mpartial_price (lp, pinf, cphase, COL_PRICING);
		}
	}

	monitor_iter (lp, pinf, it, cphase);
	if (it->nextstep == SIMPLEX_TERMINATE || it->nextstep == SIMPLEX_RESUME ||
			it->newphase != 0)
		ILL_CLEANUP;

	EGLPNUM_TYPENAME_ILLprice_primal (lp, pinf, &pr, cphase);

	if (pr.price_stat == PRICE_OPTIMAL)
	{
		//ILL_IFTRACE("%s:PRICE_OPTIMAL\n",__func__);
		if (lp->nbchange != 0)
		{
			if (it->sdisplay > 1)
			{
				QSlog("unrolling %d bound shifts", lp->nbchange);
			}
			EGLPNUM_TYPENAME_ILLfct_unroll_bound_change (lp);
			EGLPNUM_TYPENAME_ILLfct_check_pfeasible (lp, &fi, lp->tol->pfeas_tol);
			EGLPNUM_TYPENAME_ILLfct_set_status_values (lp, fi.pstatus, -1, PHASEII, -1);

			 /*HHH*/ EGLPNUM_TYPENAME_ILLfct_check_dfeasible (lp, &fi, lp->tol->dfeas_tol);
			/*HHH* QSlog("primal (opt) infeas %.6f", lp->pinfeas);
			 *HHH* QSlog("dual (opt) infeas %.6f", lp->dinfeas);*/

			if (fi.pstatus != PRIMAL_FEASIBLE)
			{
				it->algorithm = DUAL_SIMPLEX;
				it->nextstep = SIMPLEX_RESUME;
				it->resumeid = SIMPLEX_RESUME_UNSHIFT;
				it->pricetype = QS_PRICE_DDEVEX;
				/* this is to force to exit in the case of bad basis */
				//QSlog("Resume Unshift %s:%s:%d",__func__,__FILE__,__LINE__);
				EGLPNUM_TYPENAME_EGlpNumMultUiTo (lp->tol->pfeas_tol, SIMPLEX_FACTOR);
				EGLPNUM_TYPENAME_EGlpNumMultUiTo (lp->tol->dfeas_tol, SIMPLEX_FACTOR);
				it->n_restart++;
				ILL_CLEANUP;
				/*
				 * it->nextphase = PRIMAL_PHASEI;
				 * lp->tol->ip_tol /= 5.0;
				 * lp->tol->id_tol /= 5.0;
				 * ILL_CLEANUP;
				 */
			}
		}

		if (it->sdisplay > 1)
		{
			QSlog("problem seemingly solved");
			QSlog("seemingly opt = %f",
									EGLPNUM_TYPENAME_EGlpNumToLf (lp->pobjval));
			QSlog("retesting soln");
		}
		rval = EGLPNUM_TYPENAME_ILLsimplex_retest_psolution (lp, pinf, cphase, &fi);
		CHECKRVALG (rval, CLEANUP);
		EGLPNUM_TYPENAME_ILLfct_set_status_values (lp, fi.pstatus, fi.dstatus, PHASEII, PHASEII);

		if (fi.pstatus == PRIMAL_INFEASIBLE)
		{
			it->nextphase = PRIMAL_PHASEI;
			EGLPNUM_TYPENAME_EGlpNumDivUiTo (lp->tol->ip_tol, 5);
			EGLPNUM_TYPENAME_EGlpNumDivUiTo (lp->tol->id_tol, 5);
			ILL_IFTRACE ("%s:PINF:%lg\n", __func__, EGLPNUM_TYPENAME_EGlpNumToLf (lp->tol->ip_tol));
		}
		else if (fi.dstatus == DUAL_FEASIBLE)
		{
			//ILL_IFTRACE("%s:PFEAS_DFEAS\n",__func__);
			it->solstatus = ILL_LP_SOLVED;
			EGLPNUM_TYPENAME_EGlpNumCopy (lp->objval, lp->pobjval);
			it->nextstep = SIMPLEX_TERMINATE;
		}
		else
			ILL_IFTRACE ("%s:DINF:%la:%lf\n", __func__, EGLPNUM_TYPENAME_EGlpNumToLf (lp->dinfeas),
									 EGLPNUM_TYPENAME_EGlpNumToLf (lp->dinfeas));
		ILL_CLEANUP;
	}

	EGLPNUM_TYPENAME_ILLfct_compute_yz (lp, &(lp->yjz), updz, lp->nbaz[pr.eindex]);
	EGLPNUM_TYPENAME_ILLfct_update_counts (lp, CNT_YNZ, lp->yjz.nzcnt, EGLPNUM_TYPENAME_zeroLpNum);
	EGLPNUM_TYPENAME_ILLfct_update_counts (lp, CNT_UPNZ, updz->nzcnt, EGLPNUM_TYPENAME_zeroLpNum);
	ratio_iter = 0;
	do
	{
		EGLPNUM_TYPENAME_ILLratio_pII_test (lp, pr.eindex, pr.dir, &rs);
		//ILL_IFTRACE("all:%d",rs.lindex);
		EGLPNUM_TYPENAME_EGlpNumCopy (lbound, rs.lbound);
		boundch = rs.boundch;
		ratio_iter++;

		if (boundch)
		{
			/*
			 * if (ratio_iter > PARAM_PRATIOTESTS){
			 * lbound = lp->xbz[rs.lindex];
			 * boundch = 0;
			 * }
			 */
			boundch = 0;
			bndtype = (rs.lvstat == STAT_UPPER) ? BOUND_UPPER : BOUND_LOWER;
			rval = EGLPNUM_TYPENAME_ILLfct_bound_shift (lp, lp->baz[rs.lindex], bndtype, lbound);
			CHECKRVALG (rval, CLEANUP);
		}
	} while (boundch);

	if (rs.ratio_stat == RATIO_FAILED)
	{
		//ILL_IFTRACE(":%d",rs.lindex);
		/*
		 * rval = E_SIMPLEX_ERROR;
		 * it->solstatus = ILL_PPHASEII_ERROR;
		 */
		it->algorithm = DUAL_SIMPLEX;
		it->nextstep = SIMPLEX_RESUME;
		it->resumeid = SIMPLEX_RESUME_NUMER;
		/* this is to force to exit in the case of bad basis */
		it->n_restart++;
		EGLPNUM_TYPENAME_EGlpNumMultUiTo (lp->tol->pfeas_tol, SIMPLEX_FACTOR);
		EGLPNUM_TYPENAME_EGlpNumMultUiTo (lp->tol->dfeas_tol, SIMPLEX_FACTOR);
		//QSlog("Resume Numerical %s:%s:%d",__func__,__FILE__,__LINE__);
		ILL_CLEANUP;
	}
	else if (rs.ratio_stat == RATIO_UNBOUNDED)
	{
		//ILL_IFTRACE(":%d",rs.lindex);
		if (lp->nbchange != 0)
		{
			if (it->sdisplay > 1)
			{
				QSlog("unrolling %d bound shifts", lp->nbchange);
			}
			EGLPNUM_TYPENAME_ILLfct_unroll_bound_change (lp);
		}
		EGLPNUM_TYPENAME_ILLfct_set_status_values (lp, PRIMAL_UNBOUNDED, -1, PHASEII, -1);
		it->solstatus = ILL_LP_SOLVED;
		it->nextstep = SIMPLEX_TERMINATE;
		ILL_CLEANUP;
	}
	else if (rs.ratio_stat == RATIO_NOBCHANGE)
	{
		//ILL_IFTRACE(":%d",rs.lindex);
		EGLPNUM_TYPENAME_EGlpNumAddInnProdTo (lp->pobjval, rs.tz, lp->dz[pr.eindex]);
		EGLPNUM_TYPENAME_EGlpNumCopy (lp->objval, lp->pobjval);
		if (!test_progress (lp->pobjval, it->prevobj))
			it->noprog++;
		else
		{
			EGLPNUM_TYPENAME_EGlpNumCopy (it->prevobj, lp->pobjval);
			it->noprog = 0;
		}

		//ILL_IFTRACE("%s:c:%d\n",__func__,rs.lindex);
		EGLPNUM_TYPENAME_ILLfct_update_xz (lp, rs.tz, pr.eindex, rs.lindex);
		EGLPNUM_TYPENAME_ILLfct_update_basis_info (lp, pr.eindex, rs.lindex, rs.lvstat);
		if (pinf->p_strategy == COMPLETE_PRICING)
			EGLPNUM_TYPENAME_ILLprice_compute_dual_inf (lp, pinf, &pr.eindex, 1, PRIMAL_PHASEII);
		else if (pinf->p_strategy == MULTI_PART_PRICING)
			EGLPNUM_TYPENAME_ILLprice_update_mpartial_price (lp, pinf, cphase, COL_PRICING);
	}
	else if (rs.ratio_stat == RATIO_BCHANGE)
	{
		EGLPNUM_TYPENAME_EGlpNumCopyFrac (alpha, lp->dz[pr.eindex], rs.pivotval);
		EGLPNUM_TYPENAME_EGlpNumAddInnProdTo (lp->pobjval, rs.tz, lp->dz[pr.eindex]);
		EGLPNUM_TYPENAME_EGlpNumCopy (lp->objval, lp->pobjval);

		if (!test_progress (lp->pobjval, it->prevobj))
		{
			//ILL_IFTRACE(":%d",rs.lindex);
			if (lp->vtype[lp->nbaz[pr.eindex]] == VFREE ||
					lp->vtype[lp->baz[rs.lindex]] == VARTIFICIAL)
			{
				if (it->noprog > 0)
					it->noprog--;
			}
			else
				it->noprog++;
		}
		else
		{
			EGLPNUM_TYPENAME_EGlpNumCopy (it->prevobj, lp->pobjval);
			it->noprog = 0;
		}

		//ILL_IFTRACE(":%d",rs.lindex);
		EGLPNUM_TYPENAME_ILLfct_compute_zz (lp, &(lp->zz), rs.lindex);
		EGLPNUM_TYPENAME_ILLfct_update_counts (lp, CNT_ZNZ, lp->zz.nzcnt, EGLPNUM_TYPENAME_zeroLpNum);
		if (pinf->p_strategy == COMPLETE_PRICING)
		{
			//ILL_IFTRACE("%s:b\n",__func__);
			EGLPNUM_TYPENAME_ILLfct_compute_zA (lp, &(lp->zz), &(lp->zA));
			EGLPNUM_TYPENAME_ILLfct_update_counts (lp, CNT_ZANZ, lp->zA.nzcnt, EGLPNUM_TYPENAME_zeroLpNum);
			if (pinf->pII_price == QS_PRICE_PSTEEP)
				EGLPNUM_TYPENAME_ILLfct_compute_psteep_upv (lp, wz);
		}
		rval =
			EGLPNUM_TYPENAME_ILLprice_update_pricing_info (lp, pinf, cphase, wz, pr.eindex, rs.lindex,
																		rs.pivotval);
		CHECKRVALG (rval, CLEANUP);

		//ILL_IFTRACE("%s:d:%d\n",__func__,rs.lindex);
		EGLPNUM_TYPENAME_ILLfct_update_xz (lp, rs.tz, pr.eindex, rs.lindex);
		EGLPNUM_TYPENAME_ILLfct_update_basis_info (lp, pr.eindex, rs.lindex, rs.lvstat);
		rval = EGLPNUM_TYPENAME_ILLbasis_update (lp, updz, rs.lindex, &refactor, &singular);
		CHECKRVALG (rval, CLEANUP);

		if (singular)
		{
			it->nextstep = SIMPLEX_RESUME;
			it->resumeid = SIMPLEX_RESUME_SING;
			/* this is to force to exit in the case of bad basis */
			it->n_restart++;
			//QSlog("Resume Singular %s:%s:%d",__func__,__FILE__,__LINE__);
			EGLPNUM_TYPENAME_EGlpNumMultUiTo (lp->tol->pfeas_tol, SIMPLEX_FACTOR);
			EGLPNUM_TYPENAME_EGlpNumMultUiTo (lp->tol->dfeas_tol, SIMPLEX_FACTOR);
			ILL_CLEANUP;
		}
		if (!refactor)
		{
			EGLPNUM_TYPENAME_ILLfct_update_piz (lp, alpha);

			if (pinf->p_strategy == COMPLETE_PRICING)
			{
				EGLPNUM_TYPENAME_ILLfct_update_dz (lp, pr.eindex, alpha);
				EGLPNUM_TYPENAME_ILLprice_compute_dual_inf (lp, pinf, lp->zA.indx, lp->zA.nzcnt,
																	 PRIMAL_PHASEII);
				EGLPNUM_TYPENAME_ILLfct_update_counts (lp, CNT_ZARAVG, lp->zA.nzcnt, EGLPNUM_TYPENAME_zeroLpNum);
			}
			else if (pinf->p_strategy == MULTI_PART_PRICING)
			{
				EGLPNUM_TYPENAME_ILLprice_update_mpartial_price (lp, pinf, cphase, COL_PRICING);
			}
		}
		if (refactor != 0 || it->nosolve > PARAM_MAX_NOSOLVE)
		{
			EGLPNUM_TYPENAME_ILLfct_compute_xbz (lp);
			it->newphase = SIMPLEX_PHASE_RECOMP;
		}
	}

CLEANUP:
	EGLPNUM_TYPENAME_EGlpNumClearVar (alpha);
	EGLPNUM_TYPENAME_EGlpNumClearVar (lbound);
	EGLPNUM_TYPENAME_EGlpNumClearVar (fi.totinfeas);
	EGLPNUM_TYPENAME_EGlpNumClearVar (pr.dinfeas);
	EGLPNUM_TYPENAME_EGlpNumClearVar (pr.pinfeas);
	EGLPNUM_TYPENAME_EGlpNumClearVar (rs.tz);
	EGLPNUM_TYPENAME_EGlpNumClearVar (rs.lbound);
	EGLPNUM_TYPENAME_EGlpNumClearVar (rs.ecoeff);
	EGLPNUM_TYPENAME_EGlpNumClearVar (rs.pivotval);
	//EG_RETURN(rval);
	return rval;
}

static int dual_phaseI_step (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_price_info * pinf,
	EGLPNUM_TYPENAME_svector * updz,
	EGLPNUM_TYPENAME_svector * wz,
	EGLPNUM_TYPENAME_iter_info * it)
{
	int rval = 0;
	int singular = 0;
	int refactor = 0;
	int cphase = DUAL_PHASEI;
	EGLPNUM_TYPE alpha;
	EGLPNUM_TYPE alpha1;
	EGLPNUM_TYPENAME_feas_info fi;
	EGLPNUM_TYPENAME_ratio_res rs;
	EGLPNUM_TYPENAME_price_res pr;

	EGLPNUM_TYPENAME_EGlpNumInitVar (alpha);
	EGLPNUM_TYPENAME_EGlpNumInitVar (alpha1);
	EGLPNUM_TYPENAME_EGlpNumInitVar (fi.totinfeas);
	EGLPNUM_TYPENAME_EGlpNumInitVar (pr.dinfeas);
	EGLPNUM_TYPENAME_EGlpNumInitVar (pr.pinfeas);
	EGLPNUM_TYPENAME_EGlpNumInitVar (rs.tz);
	EGLPNUM_TYPENAME_EGlpNumInitVar (rs.lbound);
	EGLPNUM_TYPENAME_EGlpNumInitVar (rs.ecoeff);
	EGLPNUM_TYPENAME_EGlpNumInitVar (rs.pivotval);
	EGLPNUM_TYPENAME_EGlpNumZero (alpha1);

	EGLPNUM_TYPENAME_ILLfct_update_counts (lp, CNT_DPHASE1ITER, 0, EGLPNUM_TYPENAME_zeroLpNum);
	it->nextstep = SIMPLEX_CONTINUE;
	it->nextphase = DUAL_PHASEI;
	lp->final_phase = DUAL_PHASEI;
	it->nosolve++;

	if (it->newphase != 0)
	{
		EGLPNUM_TYPENAME_ILLfct_check_dfeasible (lp, &fi, lp->tol->id_tol);
		if (it->newphase == SIMPLEX_PHASE_NEW)
		{
			it->noprog = 0;
			if (it->sdisplay)
			{
				QSlog("starting dual phase I, nosolve %d", it->nosolve);
			}
		}
		it->newphase = 0;
		it->nosolve = 0;
		EGLPNUM_TYPENAME_EGlpNumCopy (it->prevobj, lp->dinfeas);

		EGLPNUM_TYPENAME_ILLfct_compute_phaseI_xbz (lp);
		if (pinf->d_strategy == COMPLETE_PRICING)
		{
#if USEHEAP > 0
			EGLPNUM_TYPENAME_ILLprice_free_heap (pinf);
#endif
			EGLPNUM_TYPENAME_ILLprice_compute_primal_inf (lp, pinf, NULL, 0, DUAL_PHASEI);
#if USEHEAP > 0
			rval = EGLPNUM_TYPENAME_ILLprice_test_for_heap (lp, pinf, lp->nrows, pinf->p_scaleinf,
																		 DUAL_SIMPLEX, 0);
			CHECKRVALG (rval, CLEANUP);
#endif
		}
		else if (pinf->d_strategy == MULTI_PART_PRICING)
		{
			EGLPNUM_TYPENAME_ILLprice_init_mpartial_price (lp, pinf, cphase, ROW_PRICING);
		}
	}

	monitor_iter (lp, pinf, it, cphase);
	if (it->nextstep == SIMPLEX_TERMINATE || it->nextstep == SIMPLEX_RESUME ||
			it->newphase != 0)
		ILL_CLEANUP;

	EGLPNUM_TYPENAME_ILLprice_dual (lp, pinf, cphase, &pr);

	if (pr.price_stat == PRICE_OPTIMAL)
	{
		if (it->sdisplay > 1)
		{
			QSlog("dual phase I seemingly done");
			QSlog("retesting soln");
		}

		rval = EGLPNUM_TYPENAME_ILLsimplex_retest_dsolution (lp, pinf, cphase, &fi);
		CHECKRVALG (rval, CLEANUP);
		EGLPNUM_TYPENAME_ILLfct_set_status_values (lp, fi.pstatus, fi.dstatus, PHASEI, PHASEII);

		if (fi.dstatus == DUAL_FEASIBLE)
		{
			it->nextphase = DUAL_PHASEII;
		}
		else if (fi.pstatus == PRIMAL_FEASIBLE)
		{
			it->solstatus = ILL_LP_SOLVED;
			it->nextstep = SIMPLEX_TERMINATE;
		}
		it->newphase = SIMPLEX_PHASE_NEW;
		ILL_CLEANUP;
	}

	EGLPNUM_TYPENAME_ILLfct_compute_zz (lp, &(lp->zz), pr.lindex);
	//ILL_IFTRACE("%s:c\n",__func__);
	EGLPNUM_TYPENAME_ILLfct_compute_zA (lp, &(lp->zz), &(lp->zA));
	EGLPNUM_TYPENAME_ILLfct_update_counts (lp, CNT_ZNZ, lp->zz.nzcnt, EGLPNUM_TYPENAME_zeroLpNum);
	EGLPNUM_TYPENAME_ILLfct_update_counts (lp, CNT_ZANZ, lp->zA.nzcnt, EGLPNUM_TYPENAME_zeroLpNum);

	EGLPNUM_TYPENAME_ILLratio_dI_test (lp, pr.lindex, pr.lvstat, &rs);

	if (rs.ratio_stat == RATIO_FAILED)
	{
		/*
		 * rval = E_SIMPLEX_ERROR;
		 * it->solstatus = ILL_DPHASEI_ERROR;
		 */
		it->algorithm = PRIMAL_SIMPLEX;
		it->nextstep = SIMPLEX_RESUME;
		it->resumeid = SIMPLEX_RESUME_NUMER;
		/* this is to force to exit in the case of bad basis */
		it->n_restart++;
		EGLPNUM_TYPENAME_EGlpNumMultUiTo (lp->tol->pfeas_tol, SIMPLEX_FACTOR);
		EGLPNUM_TYPENAME_EGlpNumMultUiTo (lp->tol->dfeas_tol, SIMPLEX_FACTOR);
		//QSlog("Resume Numerical %s:%s:%d",__func__,__FILE__,__LINE__);
		ILL_CLEANUP;
	}
	else if (rs.ratio_stat == RATIO_BCHANGE)
	{
		EGLPNUM_TYPENAME_ILLfct_compute_yz (lp, &(lp->yjz), updz, lp->nbaz[rs.eindex]);
		rval = EGLPNUM_TYPENAME_ILLfct_test_pivot (lp, pr.lindex, ROW_PIVOT, rs.pivotval);
		if (rval)
		{
			it->n_pivot_fail++;
			if (it->n_pivot_fail > SIMPLEX_MAX_PIVOT_FAIL)
			{
				it->n_pivot_fail = 0;
				/* this is to force to exit in the case of bad basis */
				it->n_restart++;
				it->algorithm = PRIMAL_SIMPLEX;
				it->nextstep = SIMPLEX_RESUME;
				it->resumeid = SIMPLEX_RESUME_NUMER;
				EGLPNUM_TYPENAME_EGlpNumMultUiTo (lp->tol->pfeas_tol, SIMPLEX_FACTOR);
				EGLPNUM_TYPENAME_EGlpNumMultUiTo (lp->tol->dfeas_tol, SIMPLEX_FACTOR);
				//QSlog("Resume Pivot %s:%s:%d",__func__,__FILE__,__LINE__);
				rval = 0;
				ILL_CLEANUP;
			}
			rval = EGLPNUM_TYPENAME_ILLbasis_factor (lp, &singular);
			if (singular)
				MESSAGE (__QS_SB_VERB, "Singular basis found!");
			CHECKRVALG (rval, CLEANUP);
			if (singular == 0)
				refactor = 1;
			goto END;
		}
		EGLPNUM_TYPENAME_ILLfct_update_counts (lp, CNT_YNZ, lp->yjz.nzcnt, EGLPNUM_TYPENAME_zeroLpNum);
		EGLPNUM_TYPENAME_ILLfct_update_counts (lp, CNT_UPNZ, updz->nzcnt, EGLPNUM_TYPENAME_zeroLpNum);

		if (pinf->dI_price == QS_PRICE_DSTEEP)
			EGLPNUM_TYPENAME_ILLfct_compute_dsteep_upv (lp, wz);
		rval =
			EGLPNUM_TYPENAME_ILLprice_update_pricing_info (lp, pinf, cphase, wz, rs.eindex, pr.lindex,
																		rs.pivotval);
		CHECKRVALG (rval, CLEANUP);

		EGLPNUM_TYPENAME_EGlpNumSubTo (lp->dinfeas, lp->upd.c_obj);

		if (!test_progress (lp->dinfeas, it->prevobj))
		{
			if (lp->vtype[lp->baz[pr.lindex]] == VARTIFICIAL ||
					lp->vtype[lp->nbaz[rs.eindex]] == VFREE)
			{
				if (it->noprog > 0)
					it->noprog--;
			}
			else
				it->noprog++;
		}
		else
		{
			EGLPNUM_TYPENAME_EGlpNumCopy (it->prevobj, lp->dinfeas);
			it->noprog = 0;
		}

		EGLPNUM_TYPENAME_EGlpNumCopyFrac (alpha, lp->dz[rs.eindex], rs.pivotval);
		EGLPNUM_TYPENAME_EGlpNumCopyFrac (alpha1, lp->xbz[pr.lindex], rs.pivotval);

		EGLPNUM_TYPENAME_ILLfct_update_piz (lp, alpha);
		EGLPNUM_TYPENAME_ILLfct_update_dz (lp, rs.eindex, alpha);
		EGLPNUM_TYPENAME_ILLfct_update_dfeas (lp, rs.eindex, &(lp->srhs));
		EGLPNUM_TYPENAME_ILLfct_compute_dpIy (lp, &(lp->srhs), &(lp->ssoln));
		EGLPNUM_TYPENAME_ILLfct_update_basis_info (lp, rs.eindex, pr.lindex, pr.lvstat);

#if DENSE_PI > 0
		EGLPNUM_TYPENAME_fct_test_workvector (lp);
		EGLPNUM_TYPENAME_fct_test_dfeasible (lp);
#endif
		rval = EGLPNUM_TYPENAME_ILLbasis_update (lp, updz, pr.lindex, &refactor, &singular);
		CHECKRVALG (rval, CLEANUP);

#if DENSE_NORM > 0
		EGLPNUM_TYPENAME_test_dsteep_norms (lp, pinf);
#endif

		EGLPNUM_TYPENAME_ILLfct_update_dpI_prices (lp, pinf, &(lp->srhs), &(lp->ssoln), pr.lindex,
															alpha1);

	END:
		if (singular)
		{
			it->nextstep = SIMPLEX_RESUME;
			it->resumeid = SIMPLEX_RESUME_SING;
			/* this is to force to exit in the case of bad basis */
			it->n_restart++;
			//QSlog("Resume Singular %s:%s:%d",__func__,__FILE__,__LINE__);
			EGLPNUM_TYPENAME_EGlpNumMultUiTo (lp->tol->pfeas_tol, SIMPLEX_FACTOR);
			EGLPNUM_TYPENAME_EGlpNumMultUiTo (lp->tol->dfeas_tol, SIMPLEX_FACTOR);
			ILL_CLEANUP;
		}
		if (refactor != 0 || it->nosolve > PARAM_MAX_NOSOLVE)
		{
			EGLPNUM_TYPENAME_ILLfct_compute_piz (lp);
			EGLPNUM_TYPENAME_ILLfct_compute_dz (lp);
			EGLPNUM_TYPENAME_ILLfct_dual_adjust (lp, EGLPNUM_TYPENAME_zeroLpNum);
			EGLPNUM_TYPENAME_ILLfct_check_dfeasible (lp, &fi, lp->tol->id_tol);
			EGLPNUM_TYPENAME_ILLfct_set_status_values (lp, -1, fi.dstatus, -1, PHASEII);
			if (fi.dstatus == DUAL_FEASIBLE)
				it->nextphase = DUAL_PHASEII;

			it->newphase = SIMPLEX_PHASE_RECOMP;
			ILL_CLEANUP;
		}
	}

#if DENSE_PI > 1
	EGLPNUM_TYPENAME_fct_test_workvector (lp);
	EGLPNUM_TYPENAME_fct_test_pI_x (lp, pinf);
#endif

CLEANUP:
	EGLPNUM_TYPENAME_EGlpNumClearVar (alpha);
	EGLPNUM_TYPENAME_EGlpNumClearVar (alpha1);
	EGLPNUM_TYPENAME_EGlpNumClearVar (fi.totinfeas);
	EGLPNUM_TYPENAME_EGlpNumClearVar (pr.dinfeas);
	EGLPNUM_TYPENAME_EGlpNumClearVar (pr.pinfeas);
	EGLPNUM_TYPENAME_EGlpNumClearVar (rs.tz);
	EGLPNUM_TYPENAME_EGlpNumClearVar (rs.lbound);
	EGLPNUM_TYPENAME_EGlpNumClearVar (rs.ecoeff);
	EGLPNUM_TYPENAME_EGlpNumClearVar (rs.pivotval);
	//EG_RETURN(rval);
	return rval;
}

static int dual_phaseII_step (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_price_info * pinf,
	EGLPNUM_TYPENAME_svector * updz,
	EGLPNUM_TYPENAME_svector * wz,
	EGLPNUM_TYPENAME_iter_info * it)
{
	int coeffch;
	int rval = 0;
	int singular = 0;
	int refactor = 0;
	int ratio_iter = 0;
	int cphase = DUAL_PHASEII;
	int lcol, ecol;
	int estat, newphase;
	EGLPNUM_TYPE x_bi, v_l, eval;
	EGLPNUM_TYPE ecoeff;
	EGLPNUM_TYPE alpha;
	EGLPNUM_TYPE alpha1;
	EGLPNUM_TYPENAME_feas_info fi;
	EGLPNUM_TYPENAME_ratio_res rs;
	EGLPNUM_TYPENAME_price_res pr;

	EGLPNUM_TYPENAME_EGlpNumInitVar (x_bi);
	EGLPNUM_TYPENAME_EGlpNumInitVar (v_l);
	EGLPNUM_TYPENAME_EGlpNumInitVar (eval);
	EGLPNUM_TYPENAME_EGlpNumInitVar (ecoeff);
	EGLPNUM_TYPENAME_EGlpNumInitVar (alpha);
	EGLPNUM_TYPENAME_EGlpNumInitVar (alpha1);
	EGLPNUM_TYPENAME_EGlpNumInitVar (fi.totinfeas);
	EGLPNUM_TYPENAME_EGlpNumInitVar (pr.dinfeas);
	EGLPNUM_TYPENAME_EGlpNumInitVar (pr.pinfeas);
	EGLPNUM_TYPENAME_EGlpNumInitVar (rs.tz);
	EGLPNUM_TYPENAME_EGlpNumInitVar (rs.lbound);
	EGLPNUM_TYPENAME_EGlpNumInitVar (rs.ecoeff);
	EGLPNUM_TYPENAME_EGlpNumInitVar (rs.pivotval);
	EGLPNUM_TYPENAME_EGlpNumZero (rs.ecoeff);
	EGLPNUM_TYPENAME_EGlpNumZero (alpha1);

	EGLPNUM_TYPENAME_ILLfct_update_counts (lp, CNT_DPHASE2ITER, 0, EGLPNUM_TYPENAME_zeroLpNum);
	it->nextstep = SIMPLEX_CONTINUE;
	it->nextphase = DUAL_PHASEII;
	lp->final_phase = DUAL_PHASEII;
	newphase = it->newphase;
	it->nosolve++;

	if (it->newphase != 0)
	{
		EGLPNUM_TYPENAME_ILLfct_compute_dobj (lp);
		if (it->newphase == SIMPLEX_PHASE_NEW)
		{
			it->noprog = 0;
			if (it->sdisplay)
			{
				QSlog("starting dual phase II, nosolve %d", it->nosolve);
			}
		}
		it->newphase = 0;
		it->nosolve = 0;
		EGLPNUM_TYPENAME_EGlpNumCopy (it->prevobj, lp->dobjval);
		EGLPNUM_TYPENAME_ILLfct_compute_xbz (lp);

		if (pinf->d_strategy == COMPLETE_PRICING)
		{
#if USEHEAP > 0
			EGLPNUM_TYPENAME_ILLprice_free_heap (pinf);
#endif
			EGLPNUM_TYPENAME_ILLprice_compute_primal_inf (lp, pinf, NULL, 0, DUAL_PHASEII);
#if USEHEAP > 0
			rval = EGLPNUM_TYPENAME_ILLprice_test_for_heap (lp, pinf, lp->nrows, pinf->p_scaleinf,
																		 DUAL_SIMPLEX, 0);
			CHECKRVALG (rval, CLEANUP);
#endif
		}
		else if (pinf->d_strategy == MULTI_PART_PRICING)
		{
			EGLPNUM_TYPENAME_ILLprice_init_mpartial_price (lp, pinf, cphase, ROW_PRICING);
		}
	}

	monitor_iter (lp, pinf, it, cphase);
	if (it->nextstep == SIMPLEX_TERMINATE || it->nextstep == SIMPLEX_RESUME ||
			it->newphase != 0)
		ILL_CLEANUP;

	EGLPNUM_TYPENAME_ILLprice_dual (lp, pinf, cphase, &pr);

	if (pr.price_stat == PRICE_OPTIMAL)
	{
		if (lp->ncchange != 0)
		{
			if (it->sdisplay > 1)
			{
				QSlog("unrolling %d coef shifts", lp->ncchange);
			}
			EGLPNUM_TYPENAME_ILLfct_unroll_coef_change (lp);
			EGLPNUM_TYPENAME_ILLfct_check_dfeasible (lp, &fi, lp->tol->dfeas_tol);
			EGLPNUM_TYPENAME_ILLfct_set_status_values (lp, -1, fi.dstatus, -1, PHASEII);

			 /*HHH*/ EGLPNUM_TYPENAME_ILLfct_check_pfeasible (lp, &fi, lp->tol->pfeas_tol);
			/*HHH* QSlog("dual (opt) infeas %.6f", lp->dinfeas);
			 *HHH* QSlog("primal (opt) infeas %.6f", lp->pinfeas);*/

			if (fi.dstatus != DUAL_FEASIBLE)
			{
				it->algorithm = PRIMAL_SIMPLEX;
				it->nextstep = SIMPLEX_RESUME;
				it->resumeid = SIMPLEX_RESUME_UNSHIFT;
				it->pricetype = QS_PRICE_PDEVEX;
				/* this is to force to exit in the case of bad basis */
				it->n_restart++;
				EGLPNUM_TYPENAME_EGlpNumMultUiTo (lp->tol->pfeas_tol, SIMPLEX_FACTOR);
				EGLPNUM_TYPENAME_EGlpNumMultUiTo (lp->tol->dfeas_tol, SIMPLEX_FACTOR);
				//QSlog("Resume Unshift %s:%s:%d",__func__,__FILE__,__LINE__);
				ILL_CLEANUP;
				/*
				 * it->nextphase = DUAL_PHASEI;
				 * lp->tol->ip_tol /= 5.0;
				 * lp->tol->id_tol /= 5.0;
				 * ILL_CLEANUP;
				 */
			}
		}
		if (it->sdisplay > 1)
		{
			//ILL_IFTRACE("\t%s:%s:%d\n",__func__,__FILE__,__LINE__);
			QSlog("problem seemingly solved");
			QSlog("seemingly dual opt = %f", EGLPNUM_TYPENAME_EGlpNumToLf (lp->dobjval));
			QSlog("retesting soln");
		}

		rval = EGLPNUM_TYPENAME_ILLsimplex_retest_dsolution (lp, pinf, cphase, &fi);
		CHECKRVALG (rval, CLEANUP);
		EGLPNUM_TYPENAME_ILLfct_set_status_values (lp, fi.pstatus, fi.dstatus, PHASEII, PHASEII);

		if (fi.dstatus == DUAL_INFEASIBLE)
		{
			ILL_IFTRACE ("DUAL_INFEAS: %s\n", __func__);
			it->nextphase = DUAL_PHASEI;
			EGLPNUM_TYPENAME_EGlpNumDivUiTo (lp->tol->ip_tol, 5);
			EGLPNUM_TYPENAME_EGlpNumDivUiTo (lp->tol->id_tol, 5);
		}
		else if (fi.pstatus == PRIMAL_FEASIBLE)
		{
			ILL_IFTRACE ("PRIM_FEAS: %s\n", __func__);
			EGLPNUM_TYPENAME_EGlpNumCopy (lp->objval, lp->dobjval);
			it->solstatus = ILL_LP_SOLVED;
			it->nextstep = SIMPLEX_TERMINATE;
		}
		else
			ILL_IFTRACE ("PRIM_INFEAS: %s\n", __func__);
		ILL_CLEANUP;
	}

	EGLPNUM_TYPENAME_ILLfct_compute_zz (lp, &(lp->zz), pr.lindex);
	//ILL_IFTRACE("%s:d\n",__func__);
	EGLPNUM_TYPENAME_ILLfct_compute_zA (lp, &(lp->zz), &(lp->zA));
	EGLPNUM_TYPENAME_ILLfct_update_counts (lp, CNT_ZNZ, lp->zz.nzcnt, EGLPNUM_TYPENAME_zeroLpNum);
	EGLPNUM_TYPENAME_ILLfct_update_counts (lp, CNT_ZANZ, lp->zA.nzcnt, EGLPNUM_TYPENAME_zeroLpNum);

	ratio_iter = 0;
	do
	{
		EGLPNUM_TYPENAME_ILLratio_longdII_test (lp, pr.lindex, pr.lvstat, &rs);
		if (rs.ratio_stat == RATIO_NEGATIVE)
		{
			if (it->sdisplay > 1)
			{
				QSlog("adjust coefs to remove negative ratio tests");
			}
			EGLPNUM_TYPENAME_ILLfct_adjust_viol_coefs (lp);
			EGLPNUM_TYPENAME_ILLratio_longdII_test (lp, pr.lindex, pr.lvstat, &rs);
			if (rs.ratio_stat == RATIO_NEGATIVE)
			{
				MESSAGE (__QS_SB_VERB, "internal error: bad ratio test");
				rs.ratio_stat = RATIO_FAILED;
				break;
			}
		}

		coeffch = rs.coeffch;
		EGLPNUM_TYPENAME_EGlpNumCopy (ecoeff, rs.ecoeff);
		ratio_iter++;

		if (coeffch)
		{
			/*
			 * if (ratio_iter > PARAM_DRATIOTESTS){
			 * ecoeff = lp->cz[lp->nbaz[rs.eindex]] - lp->dz[rs.eindex];
			 * coeffch = 0;
			 * }
			 */
			coeffch = 0;
			rval = EGLPNUM_TYPENAME_ILLfct_coef_shift (lp, lp->nbaz[rs.eindex], ecoeff);
			CHECKRVALG (rval, CLEANUP);
		}
		if (rs.ratio_stat == RATIO_BCHANGE)
			if (lp->vstat[lp->nbaz[rs.eindex]] == STAT_ZERO)
				break;

	} while (coeffch);

	if (rs.ratio_stat == RATIO_FAILED)
	{
		/*
		 * rval = E_SIMPLEX_ERROR;
		 * it->solstatus = ILL_DPHASEII_ERROR;
		 */
		it->algorithm = PRIMAL_SIMPLEX;
		it->nextstep = SIMPLEX_RESUME;
		it->resumeid = SIMPLEX_RESUME_NUMER;
		/* this is to force to exit in the case of bad basis */
		it->n_restart++;
		EGLPNUM_TYPENAME_EGlpNumMultUiTo (lp->tol->pfeas_tol, SIMPLEX_FACTOR);
		EGLPNUM_TYPENAME_EGlpNumMultUiTo (lp->tol->dfeas_tol, SIMPLEX_FACTOR);
		//QSlog("Resume Numerical %s:%s:%d",__func__,__FILE__,__LINE__);
		ILL_CLEANUP;
	}
	else if (rs.ratio_stat == RATIO_UNBOUNDED)
	{
		lp->infub_ix = pr.lindex;
		if (lp->ncchange != 0)
		{
			if (it->sdisplay > 1)
			{
				QSlog("unrolling %d coef shifts", lp->ncchange);
			}
			EGLPNUM_TYPENAME_ILLfct_unroll_coef_change (lp);
		}
		EGLPNUM_TYPENAME_ILLfct_set_status_values (lp, -1, DUAL_UNBOUNDED, -1, PHASEII);
		it->solstatus = ILL_LP_SOLVED;
		it->nextstep = SIMPLEX_TERMINATE;
	}
	else if (rs.ratio_stat == RATIO_BCHANGE)
	{
		lcol = lp->baz[pr.lindex];
		ecol = lp->nbaz[rs.eindex];

		EGLPNUM_TYPENAME_ILLfct_compute_yz (lp, &(lp->yjz), updz, ecol);
		EGLPNUM_TYPENAME_ILLfct_update_counts (lp, CNT_YNZ, lp->yjz.nzcnt, EGLPNUM_TYPENAME_zeroLpNum);
		EGLPNUM_TYPENAME_ILLfct_update_counts (lp, CNT_UPNZ, updz->nzcnt, EGLPNUM_TYPENAME_zeroLpNum);
		rval = EGLPNUM_TYPENAME_ILLfct_test_pivot (lp, pr.lindex, ROW_PIVOT, rs.pivotval);
		if (rval != 0)
		{
			it->n_pivot_fail++;
			if (it->n_pivot_fail > SIMPLEX_MAX_PIVOT_FAIL)
			{
				it->n_pivot_fail = 0;
				/* this is to force to exit in the case of bad basis */
				it->algorithm = PRIMAL_SIMPLEX;
				it->nextstep = SIMPLEX_RESUME;
				it->resumeid = SIMPLEX_RESUME_NUMER;
				it->n_restart++;
				EGLPNUM_TYPENAME_EGlpNumMultUiTo (lp->tol->pfeas_tol, SIMPLEX_FACTOR);
				EGLPNUM_TYPENAME_EGlpNumMultUiTo (lp->tol->dfeas_tol, SIMPLEX_FACTOR);
				//QSlog("Resume Pivot %s:%s:%d",__func__,__FILE__,__LINE__);
				rval = 0;
				ILL_CLEANUP;
			}
			if (newphase == 0)
			{
				rval = EGLPNUM_TYPENAME_ILLbasis_factor (lp, &singular);
				CHECKRVALG (rval, CLEANUP);
				if (singular)
					MESSAGE (__QS_SB_VERB, "Singular basis found!");
#ifdef dbl_QSOPT_CURRENT_PRECICION
				if (singular)
				{
					MESSAGE (__QS_SB_VERB, "Forcing fail!");
					rval = QS_LP_CHANGE_PREC;
				}
#endif
				if (singular == 0)
					refactor = 1;
				goto END;
			}
			else
			{
				if (it->sdisplay > 1)
				{
					QSlog("warning: bad step");
				}
			}
		}

		EGLPNUM_TYPENAME_EGlpNumAddTo (lp->dobjval, lp->upd.c_obj);
		EGLPNUM_TYPENAME_EGlpNumCopy (lp->objval, lp->dobjval);

		if (!test_progress (lp->dobjval, it->prevobj))
		{
			if (lp->vtype[lcol] == VARTIFICIAL || lp->vtype[ecol] == VFREE)
			{
				if (it->noprog > 0)
					it->noprog--;
			}
			else
				it->noprog++;
		}
		else
		{
			EGLPNUM_TYPENAME_EGlpNumCopy (it->prevobj, lp->dobjval);
			it->noprog = 0;
		}

		if (pinf->dII_price == QS_PRICE_DSTEEP)
			EGLPNUM_TYPENAME_ILLfct_compute_dsteep_upv (lp, wz);
		rval =
			EGLPNUM_TYPENAME_ILLprice_update_pricing_info (lp, pinf, cphase, wz, rs.eindex, pr.lindex,
																		rs.pivotval);
		CHECKRVALG (rval, CLEANUP);

		EGLPNUM_TYPENAME_EGlpNumCopy (x_bi, lp->xbz[pr.lindex]);
		if (pr.lvstat == STAT_LOWER)
			EGLPNUM_TYPENAME_EGlpNumCopy (v_l, lp->lz[lcol]);
		else
			EGLPNUM_TYPENAME_EGlpNumCopy (v_l, lp->uz[lcol]);
		EGLPNUM_TYPENAME_EGlpNumCopy (alpha, rs.tz);
		if (pr.lvstat == STAT_LOWER)
			EGLPNUM_TYPENAME_EGlpNumSign (alpha);
		estat = lp->vstat[ecol];
		if (estat == STAT_LOWER)
			EGLPNUM_TYPENAME_EGlpNumCopy (eval, lp->lz[ecol]);
		else if (estat == STAT_ZERO)
			EGLPNUM_TYPENAME_EGlpNumZero (eval);
		else
			EGLPNUM_TYPENAME_EGlpNumCopy (eval, lp->uz[ecol]);

		EGLPNUM_TYPENAME_ILLfct_update_piz (lp, alpha);
		EGLPNUM_TYPENAME_ILLfct_update_dz (lp, rs.eindex, alpha);
		EGLPNUM_TYPENAME_ILLfct_update_dIIfeas (lp, rs.eindex, &(lp->srhs));
		EGLPNUM_TYPENAME_ILLfct_compute_dpIIy (lp, &(lp->srhs), &(lp->ssoln));
		EGLPNUM_TYPENAME_EGlpNumCopyDiff (alpha1, x_bi, v_l);
		EGLPNUM_TYPENAME_EGlpNumSubTo (alpha1, lp->upd.dty);
		EGLPNUM_TYPENAME_EGlpNumDivTo (alpha1, rs.pivotval);
		EGLPNUM_TYPENAME_ILLfct_update_basis_info (lp, rs.eindex, pr.lindex, pr.lvstat);
		rval = EGLPNUM_TYPENAME_ILLbasis_update (lp, updz, pr.lindex, &refactor, &singular);
		CHECKRVALG (rval, CLEANUP);

		EGLPNUM_TYPENAME_ILLfct_update_dpII_prices (lp, pinf, &(lp->srhs), &(lp->ssoln), /*rs.eindex,*/
															 pr.lindex, eval, alpha1);

#if DENSE_NORM > 0
		EGLPNUM_TYPENAME_test_dsteep_norms (lp, pinf);
#endif

	END:
		if (singular)
		{
			it->nextstep = SIMPLEX_RESUME;
			it->resumeid = SIMPLEX_RESUME_SING;
			/* this is to force to exit in the case of bad basis */
			it->n_restart++;
			EGLPNUM_TYPENAME_EGlpNumMultUiTo (lp->tol->pfeas_tol, SIMPLEX_FACTOR);
			EGLPNUM_TYPENAME_EGlpNumMultUiTo (lp->tol->dfeas_tol, SIMPLEX_FACTOR);
			//QSlog("Resume Singular %s:%s:%d",__func__,__FILE__,__LINE__);
			ILL_CLEANUP;
		}
		if (refactor != 0 || it->nosolve > PARAM_MAX_NOSOLVE)
		{
			EGLPNUM_TYPENAME_ILLfct_compute_piz (lp);
			EGLPNUM_TYPENAME_ILLfct_compute_dz (lp);
			EGLPNUM_TYPENAME_ILLfct_dual_adjust (lp, EGLPNUM_TYPENAME_zeroLpNum);
			it->newphase = SIMPLEX_PHASE_RECOMP;
		}
	}

#if DENSE_PIIPI > 0
	EGLPNUM_TYPENAME_fct_test_workvector (lp);
	if (!refactor)
	{
		EGLPNUM_TYPENAME_fct_test_pII_x (lp, pinf);
		EGLPNUM_TYPENAME_fct_test_pII_pi_dz (lp, pinf);
	}
#endif

CLEANUP:
	EGLPNUM_TYPENAME_EGlpNumClearVar (x_bi);
	EGLPNUM_TYPENAME_EGlpNumClearVar (v_l);
	EGLPNUM_TYPENAME_EGlpNumClearVar (eval);
	EGLPNUM_TYPENAME_EGlpNumClearVar (ecoeff);
	EGLPNUM_TYPENAME_EGlpNumClearVar (alpha);
	EGLPNUM_TYPENAME_EGlpNumClearVar (alpha1);
	EGLPNUM_TYPENAME_EGlpNumClearVar (fi.totinfeas);
	EGLPNUM_TYPENAME_EGlpNumClearVar (pr.dinfeas);
	EGLPNUM_TYPENAME_EGlpNumClearVar (pr.pinfeas);
	EGLPNUM_TYPENAME_EGlpNumClearVar (rs.tz);
	EGLPNUM_TYPENAME_EGlpNumClearVar (rs.lbound);
	EGLPNUM_TYPENAME_EGlpNumClearVar (rs.ecoeff);
	EGLPNUM_TYPENAME_EGlpNumClearVar (rs.pivotval);
	//EG_RETURN(rval);
	return rval;
}

static void get_current_stat (
	EGLPNUM_TYPENAME_lp_status_info * p,
	int algorithm,
	int *bstat)
{
	if (p->optimal)
		*bstat = OPTIMAL;
	else if (algorithm == PRIMAL_SIMPLEX)
	{
		if (p->primal_feasible)
			*bstat = PRIMAL_FEASIBLE;
		else if (p->primal_infeasible)
			*bstat = PRIMAL_INFEASIBLE;
		else if (p->primal_unbounded)
			*bstat = PRIMAL_UNBOUNDED;
		else
			*bstat = NONOPTIMAL;
	}
	else if (algorithm == DUAL_SIMPLEX)
	{
		if (p->dual_feasible)
			*bstat = DUAL_FEASIBLE;
		else if (p->dual_infeasible)
			*bstat = DUAL_INFEASIBLE;
		else if (p->dual_unbounded)
			*bstat = DUAL_UNBOUNDED;
		else
			*bstat = NONOPTIMAL;
	}
}

int EGLPNUM_TYPENAME_ILLsimplex_pivotin (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_price_info * pinf,
	int rcnt,
	int *rlist,
	int pivot_opt,
	int *basis_mod)
{
	int i, npiv = 0;
	int eindex;
	int rval = 0;
	int singular = 0;
	int refactor = 0;
	int *rowmap = lp->O->rowmap;
	int *clist = NULL;
	EGLPNUM_TYPENAME_svector wz;
	EGLPNUM_TYPENAME_svector updz;
	EGLPNUM_TYPE alpha;
	EGLPNUM_TYPENAME_ratio_res rs;
	EGLPNUM_TYPENAME_feas_info fi;

	EGLPNUM_TYPENAME_EGlpNumInitVar (alpha);
	EGLPNUM_TYPENAME_EGlpNumInitVar (fi.totinfeas);
	EGLPNUM_TYPENAME_EGlpNumInitVar (rs.tz);
	EGLPNUM_TYPENAME_EGlpNumInitVar (rs.lbound);
	EGLPNUM_TYPENAME_EGlpNumInitVar (rs.ecoeff);
	EGLPNUM_TYPENAME_EGlpNumInitVar (rs.pivotval);
	EGLPNUM_TYPENAME_EGlpNumZero (alpha);

	*basis_mod = 0;
	if (rcnt <= 0)
	{
		EG_RETURN (rval);
	}

	if (pivot_opt == SIMPLEX_PIVOTINROW)
	{
		ILL_SAFE_MALLOC (clist, rcnt, int);

		for (i = 0; i < rcnt; i++)
			clist[i] = rowmap[rlist[i]];
	}
	else
		clist = rlist;

	for (i = 0; i < rcnt; i++)
	{
		if (lp->vstat[clist[i]] != STAT_BASIC)
		{
			*basis_mod = 1;
			break;
		}
	}
	if (*basis_mod == 0)
	{
		if (pivot_opt == SIMPLEX_PIVOTINROW)
		{
			ILL_IFFREE (clist, int);
		}
		EG_RETURN (rval);
	}

	/* QSlog("Forcing vars into basis in EGLPNUM_TYPENAME_ILLsimplex_pivotin"); */
	EGLPNUM_TYPENAME_ILLsvector_init (&wz);
	rval = EGLPNUM_TYPENAME_ILLsvector_alloc (&wz, lp->nrows);
	CHECKRVALG (rval, CLEANUP);
	EGLPNUM_TYPENAME_ILLsvector_init (&updz);
	rval = EGLPNUM_TYPENAME_ILLsvector_alloc (&updz, lp->nrows);
	CHECKRVALG (rval, CLEANUP);

	EGLPNUM_TYPENAME_EGlpNumCopy (lp->pobjval, lp->dobjval);
	for (i = 0; i < rcnt; i++)
	{
		if (lp->vstat[clist[i]] == STAT_BASIC)
			continue;
		npiv++;

		eindex = lp->vindex[clist[i]];
		EGLPNUM_TYPENAME_ILLfct_compute_yz (lp, &(lp->yjz), &updz, lp->nbaz[eindex]);
		EGLPNUM_TYPENAME_ILLfct_update_counts (lp, CNT_YNZ, lp->yjz.nzcnt, EGLPNUM_TYPENAME_zeroLpNum);
		EGLPNUM_TYPENAME_ILLfct_update_counts (lp, CNT_UPNZ, updz.nzcnt, EGLPNUM_TYPENAME_zeroLpNum);

		EGLPNUM_TYPENAME_ILLratio_pivotin_test (lp, clist, rcnt, &rs);

		if (rs.ratio_stat == RATIO_UNBOUNDED || rs.ratio_stat == RATIO_FAILED)
		{
			QSlog("Pivot_in failed");
			rval = E_SIMPLEX_ERROR;
			ILL_CLEANUP;
		}
		else if (rs.ratio_stat == RATIO_BCHANGE)
		{
			if (rs.lvstat == STAT_LOWER)
			{
				EGLPNUM_TYPENAME_EGlpNumCopyDiff (alpha, lp->lz[lp->baz[rs.lindex]], lp->xbz[rs.lindex]);
				EGLPNUM_TYPENAME_EGlpNumAddInnProdTo (lp->dobjval, rs.tz, alpha);
			}
			else
			{
				EGLPNUM_TYPENAME_EGlpNumCopyDiff (alpha, lp->xbz[rs.lindex], lp->uz[lp->baz[rs.lindex]]);
				EGLPNUM_TYPENAME_EGlpNumAddInnProdTo (lp->dobjval, rs.tz, alpha);
			}
			EGLPNUM_TYPENAME_EGlpNumCopyFrac (alpha, lp->dz[eindex], rs.pivotval);
			EGLPNUM_TYPENAME_EGlpNumCopy (lp->objval, lp->dobjval);

			EGLPNUM_TYPENAME_ILLfct_compute_zz (lp, &(lp->zz), rs.lindex);
			//ILL_IFTRACE("%s:e\n",__func__);
			EGLPNUM_TYPENAME_ILLfct_compute_zA (lp, &(lp->zz), &(lp->zA));
			EGLPNUM_TYPENAME_ILLfct_update_counts (lp, CNT_ZNZ, lp->zz.nzcnt, EGLPNUM_TYPENAME_zeroLpNum);
			EGLPNUM_TYPENAME_ILLfct_update_counts (lp, CNT_ZANZ, lp->zA.nzcnt, EGLPNUM_TYPENAME_zeroLpNum);

			if (pinf->dsinfo.norms && pinf->dII_price == QS_PRICE_DSTEEP)
			{
				EGLPNUM_TYPENAME_ILLfct_compute_dsteep_upv (lp, &wz);
				rval = EGLPNUM_TYPENAME_ILLprice_update_pricing_info (lp, pinf, DUAL_PHASEII, &wz,
																						 eindex, rs.lindex, rs.pivotval);
				CHECKRVALG (rval, CLEANUP);
			}
			else if (pinf->psinfo.norms && pinf->pII_price == QS_PRICE_PSTEEP)
			{
				EGLPNUM_TYPENAME_ILLfct_compute_psteep_upv (lp, &wz);
				rval = EGLPNUM_TYPENAME_ILLprice_update_pricing_info (lp, pinf, PRIMAL_PHASEII, &wz,
																						 eindex, rs.lindex, rs.pivotval);
				CHECKRVALG (rval, CLEANUP);
			}

			//ILL_IFTRACE("%s:e\n",__func__);
			EGLPNUM_TYPENAME_ILLfct_update_xz (lp, rs.tz, eindex, rs.lindex);
			EGLPNUM_TYPENAME_ILLfct_update_basis_info (lp, eindex, rs.lindex, rs.lvstat);
			rval = EGLPNUM_TYPENAME_ILLbasis_update (lp, &updz, rs.lindex, &refactor, &singular);
			CHECKRVALG (rval, CLEANUP);

			if (singular)
			{
				QSlog("singular matrix in pivot_in");
				rval = E_SIMPLEX_ERROR;
				ILL_CLEANUP;
			}
			if (!refactor)
			{
				EGLPNUM_TYPENAME_ILLfct_update_piz (lp, alpha);
				EGLPNUM_TYPENAME_ILLfct_update_dz (lp, eindex, alpha);
			}
			else
			{
				EGLPNUM_TYPENAME_ILLfct_compute_xbz (lp);
				EGLPNUM_TYPENAME_ILLfct_compute_piz (lp);
				EGLPNUM_TYPENAME_ILLfct_compute_dz (lp);
				EGLPNUM_TYPENAME_ILLfct_compute_dobj (lp);
			}
		}
	}
	/*
	 * EGLPNUM_TYPENAME_ILLfct_dphaseI_simple_update (lp, lp->tol->dfeas_tol);
	 * EGLPNUM_TYPENAME_ILLfct_compute_xbz (lp);
	 * EGLPNUM_TYPENAME_ILLfct_compute_dobj (lp);
	 */

	EGLPNUM_TYPENAME_ILLfct_check_pfeasible (lp, &fi, lp->tol->pfeas_tol);
	EGLPNUM_TYPENAME_ILLfct_check_dfeasible (lp, &fi, lp->tol->dfeas_tol);
	EGLPNUM_TYPENAME_ILLfct_set_status_values (lp, fi.pstatus, fi.dstatus, PHASEII, PHASEII);

CLEANUP:
	if (pivot_opt == SIMPLEX_PIVOTINROW)
		ILL_IFFREE (clist, int);

	EGLPNUM_TYPENAME_ILLsvector_free (&wz);
	EGLPNUM_TYPENAME_ILLsvector_free (&updz);
	EGLPNUM_TYPENAME_EGlpNumClearVar (alpha);
	EGLPNUM_TYPENAME_EGlpNumClearVar (fi.totinfeas);
	EGLPNUM_TYPENAME_EGlpNumClearVar (rs.tz);
	EGLPNUM_TYPENAME_EGlpNumClearVar (rs.lbound);
	EGLPNUM_TYPENAME_EGlpNumClearVar (rs.ecoeff);
	EGLPNUM_TYPENAME_EGlpNumClearVar (rs.pivotval);
	EG_RETURN (rval);
}

static int report_value (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_iter_info * it,
	const char *value_name,
	EGLPNUM_TYPE value)
{
	int rval = 0;

	if (it->sdisplay && it->itercnt % lp->iterskip == 0)
	{
		char buffer[1024];

		snprintf (buffer, (size_t) 1023, "(%d): %s = %10.7lf", it->itercnt,
							value_name, EGLPNUM_TYPENAME_EGlpNumToLf (value));
		rval = ILLstring_report (buffer, &lp->O->reporter);
	}
	else
	{
		/* make sure ILLstring_report is called at least every 10 iterations */
		if (it->itercnt % (lp->iterskip / 10))
		{
			rval = ILLstring_report (NULL, &lp->O->reporter);
		}
	}
	if (rval != 0)
	{															/* ILLstring_report was called and failed, which means we should abort */
		it->solstatus = QS_LP_ABORTED;
	}
	return rval;
}

#undef QSOPT_CURRENT_PRECICION
