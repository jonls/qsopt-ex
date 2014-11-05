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

#include "eg_lpnum.h"
#include "eg_io.h"

#define QSOPT_CURRENT_PRECICION
#include "basicdefs.h"
#include "config.h"
#include "iqsutil.h"
#include "lpdata.h"
#include "lpdefs.h"

#include "stddefs.h"
#include "fct.h"
#include "ratio.h"
#include "price.h"
#include "basis.h"
#include "simplex.h"
#include "dstruct.h"
#include "qstruct.h"
#include "qsopt.h"
#include "lib.h"								/* for ILLlib_writebasis */
#include "lp.h"									/* for ILLwrite_lp */

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

static void init_lp_status_info (
	lp_status_info * ls),
  init_simplex_tols (
	lpinfo * lp),
  monitor_iter (
	lpinfo * lp,
	price_info * p,
	iter_info * it,
	int cphase),
  get_current_stat (
	lp_status_info * p,
	int algorithm,
	int *bstat);

static int terminate_simplex (
	lpinfo * lp,
	int phase,
	iter_info * it),
  primal_phaseI_step (
	lpinfo * lp,
	price_info * pinf,
	svector * updz,
	svector * wz,
	iter_info * it),
  primal_phaseII_step (
	lpinfo * lp,
	price_info * pinf,
	svector * updz,
	svector * wz,
	iter_info * it),
  dual_phaseI_step (
	lpinfo * lp,
	price_info * pinf,
	svector * updz,
	svector * wz,
	iter_info * it),
  dual_phaseII_step (
	lpinfo * lp,
	price_info * pinf,
	svector * updz,
	svector * wz,
	iter_info * it),
  report_value (
	lpinfo * lp,
	iter_info * it,
	const char *value_name,
	EGlpNum_t value);


void ILLsimplex_init_lpinfo (
	lpinfo * lp)
{
	ILLbasis_init_basisinfo (lp);
	init_internal_lpinfo (lp);
}

void ILLsimplex_free_lpinfo (
	lpinfo * lp)
{
	if (lp)
	{
		EGlpNumFreeArray (lp->lz);
		EGlpNumFreeArray (lp->uz);
		EGlpNumFreeArray (lp->cz);
		ILLbasis_free_basisinfo (lp);
		free_internal_lpinfo (lp);
	}
}

void ILLsimplex_load_lpinfo (
	ILLlpdata * qslp,
	lpinfo * lp)
{
	lp->basisid = -1;
	lp->maxiter = 500000;
	lp->maxtime = 300000;
	//lp->iterskip = 10;
	lp->iterskip = 100;
	EGlpNumCopy (lp->objbound, INFTY);
	lp->O = qslp;
}

void ILLsimplex_set_bound (
	lpinfo * lp,
	const EGlpNum_t * objbound,
	int sense)
{
	EGlpNumCopy (lp->objbound, *objbound);
	if (sense == ILL_MAX)
		EGlpNumSign (lp->objbound);
}

static void init_lp_status_info (
	lp_status_info * ls)
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
	lpinfo * lp)
{
	EGlpNumCopy (lp->tol->pfeas_tol, PFEAS_TOLER);
	EGlpNumCopy (lp->tol->dfeas_tol, DFEAS_TOLER);
	EGlpNumCopy (lp->tol->pivot_tol, PIVOT_TOLER);
	EGlpNumCopy (lp->tol->szero_tol, SZERO_TOLER);
	EGlpNumCopy (lp->tol->ip_tol, lp->tol->pfeas_tol);
	EGlpNumCopy (lp->tol->id_tol, lp->tol->dfeas_tol);
	if (EGlpNumIsNeqqZero (lp->tol->ip_tol))
	{
#if VERBOSE_LEVEL <= DEBUG
		MESSAGE (VERBOSE_LEVEL, "ip_tol %lg", EGlpNumToLf (lp->tol->ip_tol));
		MESSAGE (VERBOSE_LEVEL, "eps %lg", EGlpNumToLf (epsLpNum));
		MESSAGE (VERBOSE_LEVEL, "PFEAS_TOLER %lg", EGlpNumToLf (PFEAS_TOLER));
#endif
		EGlpNumDivUiTo (lp->tol->ip_tol, 2UL);
	}
	if (EGlpNumIsNeqqZero (lp->tol->id_tol))
	{
#if VERBOSE_LEVEL <= DEBUG
		MESSAGE (VERBOSE_LEVEL, "id_tol %lg", EGlpNumToLf (lp->tol->id_tol));
#endif
		EGlpNumDivUiTo (lp->tol->id_tol, 2UL);
	}
}

void init_internal_lpinfo (
	lpinfo * lp)
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
	ILLsvector_init (&(lp->zz));
	ILLsvector_init (&(lp->yjz));
	ILLsvector_init (&(lp->zA));
	ILLsvector_init (&(lp->work));
	ILLsvector_init (&(lp->srhs));
	ILLsvector_init (&(lp->ssoln));
	ILL_SAFE_MALLOC (lp->tol, 1, tol_struct);
	EGlpNumInitVar (lp->tol->pfeas_tol);
	EGlpNumInitVar (lp->tol->dfeas_tol);
	EGlpNumInitVar (lp->tol->pivot_tol);
	EGlpNumInitVar (lp->tol->szero_tol);
	EGlpNumInitVar (lp->tol->ip_tol);
	EGlpNumInitVar (lp->tol->id_tol);
	ILL_SAFE_MALLOC (lp->cnts, 1, count_struct);
	EGlpNumInitVar (lp->cnts->y_ravg);
	EGlpNumInitVar (lp->cnts->z_ravg);
	EGlpNumInitVar (lp->cnts->za_ravg);
CLEANUP:
	if (rval)
	{
		fprintf (stderr, "\nno memory, in %s, exit\n", __func__);
		exit (1);
	}
}

void free_internal_lpinfo (
	lpinfo * lp)
{
	bndinfo *binfo = 0;
	coefinfo *cinfo = 0;

	if (lp->localrows)
	{
		ILL_IFFREE (lp->rowcnt, int);
		ILL_IFFREE (lp->rowbeg, int);
		ILL_IFFREE (lp->rowind, int);

		EGlpNumFreeArray (lp->rowval);
		lp->localrows = 0;
	}
	EGlpNumFreeArray (lp->lz);
	EGlpNumFreeArray (lp->uz);
	EGlpNumFreeArray (lp->cz);
	EGlpNumFreeArray (lp->xbz);
	EGlpNumFreeArray (lp->piz);
	EGlpNumFreeArray (lp->pIpiz);
	EGlpNumFreeArray (lp->dz);
	EGlpNumFreeArray (lp->pIdz);
	EGlpNumFreeArray (lp->pIxbz);

	ILL_IFFREE (lp->vtype, int);
	ILL_IFFREE (lp->vclass, char);

	ILLsvector_free (&(lp->zz));
	ILLsvector_free (&(lp->yjz));
	ILLsvector_free (&(lp->zA));
	ILLsvector_free (&(lp->work));
	ILLsvector_free (&(lp->srhs));
	ILLsvector_free (&(lp->ssoln));
	ILL_IFFREE (lp->iwork, int);
	ILL_IFFREE (lp->upd.perm, int);
	ILL_IFFREE (lp->upd.ix, int);

	EGlpNumFreeArray (lp->upd.t);

	ILL_IFFREE (lp->bfeas, int);
	ILL_IFFREE (lp->dfeas, int);

	if (lp->tol)
	{
		EGlpNumClearVar (lp->tol->pfeas_tol);
		EGlpNumClearVar (lp->tol->dfeas_tol);
		EGlpNumClearVar (lp->tol->pivot_tol);
		EGlpNumClearVar (lp->tol->szero_tol);
		EGlpNumClearVar (lp->tol->ip_tol);
		EGlpNumClearVar (lp->tol->id_tol);
		ILL_IFFREE (lp->tol, tol_struct);
	}
	if (lp->cnts)
	{
		EGlpNumClearVar (lp->cnts->y_ravg);
		EGlpNumClearVar (lp->cnts->z_ravg);
		EGlpNumClearVar (lp->cnts->za_ravg);
		ILL_IFFREE (lp->cnts, count_struct);
	}

	while (lp->bchanges)
	{
		binfo = lp->bchanges;
		EGlpNumClearVar (binfo->pbound);
		EGlpNumClearVar (binfo->cbound);
		lp->bchanges = binfo->next;
		ILL_IFFREE (binfo, bndinfo);
	}

	while (lp->cchanges)
	{
		cinfo = lp->cchanges;
		EGlpNumClearVar (cinfo->pcoef);
		EGlpNumClearVar (cinfo->ccoef);
		lp->cchanges = cinfo->next;
		ILL_IFFREE (cinfo, coefinfo);
	}
}

int build_internal_lpinfo (
	lpinfo * lp)
{
	int rval = 0;
	int i, n;
	ILLlpdata *qslp = lp->O;
	ILLlp_sinfo *S = lp->O->sinfo;
	EGlpNum_t *lower, *upper, *obj;
	ILLlp_rows lprows;
	ILLmatrix *A;

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

	lp->lz = EGlpNumAllocArray (lp->ncols);
	lp->cz = EGlpNumAllocArray (lp->ncols);
	lp->uz = EGlpNumAllocArray (lp->ncols);
	if (!lp->lz || !lp->uz || !lp->cz)
	{
		fprintf (stderr, "build_internal_lpinfo\n");
		rval = 1;
		goto CLEANUP;
	}
	for (i = 0; i < lp->ncols; i++)
	{
		EGlpNumCopy (lp->lz[i], lower[i]);
		EGlpNumCopy (lp->uz[i], upper[i]);
		EGlpNumCopy (lp->cz[i], obj[i]);
		if (qslp->objsense == ILL_MAX)
		{
			EGlpNumSign (lp->cz[i]);
		}
	}

	if (!lp->O->rA)
	{
		rval = ILLlp_rows_init (&lprows, lp->O, 1);
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

	lp->xbz = EGlpNumAllocArray (lp->nrows);
	lp->piz = EGlpNumAllocArray (lp->nrows);
	lp->dz = EGlpNumAllocArray (lp->nnbasic);
	lp->final_phase = -1;
	lp->infub_ix = -1;

	ILL_SAFE_MALLOC (lp->vtype, lp->ncols, int);
	ILL_SAFE_MALLOC (lp->vclass, lp->ncols, char);

	rval = ILLsvector_alloc (&(lp->zz), lp->nrows);
	CHECKRVALG (rval, CLEANUP);
	rval = ILLsvector_alloc (&(lp->yjz), lp->nrows);
	CHECKRVALG (rval, CLEANUP);
	rval = ILLsvector_alloc (&(lp->zA), lp->nnbasic);
	CHECKRVALG (rval, CLEANUP);
	rval = ILLsvector_alloc (&(lp->work), lp->ncols);
	CHECKRVALG (rval, CLEANUP);
	rval = ILLsvector_alloc (&(lp->srhs), lp->nrows);
	CHECKRVALG (rval, CLEANUP);
	rval = ILLsvector_alloc (&(lp->ssoln), lp->nrows);
	CHECKRVALG (rval, CLEANUP);
	ILL_SAFE_MALLOC (lp->iwork, lp->ncols, int);

	for (i = 0; i < lp->ncols; i++)
	{
		lp->work.indx[i] = 0;
		EGlpNumZero (lp->work.coef[i]);
		lp->iwork[i] = 0;
	}
	n = lp->nrows > lp->ncols ? 2 * (lp->nrows) + 1 : 2 * (lp->ncols) + 1;
	lp->upd.t = EGlpNumAllocArray (n);
	ILL_SAFE_MALLOC (lp->upd.perm, n, int);
	ILL_SAFE_MALLOC (lp->upd.ix, n, int);


	ILL_SAFE_MALLOC (lp->bfeas, lp->nrows, int);
	ILL_SAFE_MALLOC (lp->dfeas, lp->nnbasic, int);

	init_simplex_tols (lp);
	ILLfct_init_counts (lp);

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
		free_internal_lpinfo (lp);
	EG_RETURN (rval);
}

int ILLsimplex_retest_psolution (
	lpinfo * lp,
	price_info * p,
	int phase,
	feas_info * fi)
{
	int rval = 0;
	int fbid = lp->fbasisid;
	int bid = lp->basisid;
	EGlpNum_t *ptol = &(lp->tol->pfeas_tol);
	EGlpNum_t *dtol = &(lp->tol->dfeas_tol);
	EGlpNum_t *iptol = &(lp->tol->ip_tol);
	EGlpNum_t *idtol = &(lp->tol->id_tol);

	fi->pstatus = -1;
	fi->dstatus = -1;
	if (fbid < bid - PARAM_PRIMAL_REFACTORGAP)
	{
		rval = ILLbasis_refactor (lp);
		CHECKRVALG (rval, CLEANUP);
	}
	if (fbid < bid - PARAM_PRIMAL_RESOLVEGAP)
		ILLfct_compute_xbz (lp);

	if (phase == PRIMAL_PHASEII)
	{
		if (fbid < bid - PARAM_PRIMAL_RESOLVEGAP)
		{
			ILLfct_compute_piz (lp);
			ILLfct_compute_dz (lp);
			if (p != NULL && p->p_strategy == COMPLETE_PRICING)
				ILLprice_compute_dual_inf (lp, p, NULL, 0, PRIMAL_PHASEII);
		}
		ILLfct_compute_pobj (lp);
		ILLfct_check_pfeasible (lp, fi, *ptol);
		ILLfct_check_dfeasible (lp, fi, *dtol);
	}
	else if (phase == PRIMAL_PHASEI)
	{
		ILLfct_check_pfeasible (lp, fi, *iptol);
		if (fi->pstatus != PRIMAL_FEASIBLE)
		{
			if (lp->pIpiz)
			{
				ILLfct_compute_phaseI_piz (lp);
				ILLfct_compute_phaseI_dz (lp);
				ILLfct_check_pIdfeasible (lp, fi, *idtol);
				if (p != NULL && p->p_strategy == COMPLETE_PRICING)
					ILLprice_compute_dual_inf (lp, p, NULL, 0, PRIMAL_PHASEI);
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

int ILLsimplex_retest_dsolution (
	lpinfo * lp,
	price_info * p,
	int phase,
	feas_info * fi)
{
	int rval = 0;
	int fbid = lp->fbasisid;
	int bid = lp->basisid;
	EGlpNum_t *ptol = &(lp->tol->pfeas_tol);
	EGlpNum_t *dtol = &(lp->tol->dfeas_tol);
	EGlpNum_t *iptol = &(lp->tol->ip_tol);
	EGlpNum_t *idtol = &(lp->tol->id_tol);

	//ILL_IFTRACE("\t%s:%s:%d\n",__func__,__FILE__,__LINE__);

	fi->pstatus = -1;
	fi->dstatus = -1;
	if (fbid < bid - PARAM_DUAL_REFACTORGAP)
	{
		//ILL_IFTRACE("Refactor: %s:%s:%d\n",__func__,__FILE__,__LINE__);
		rval = ILLbasis_refactor (lp);
		CHECKRVALG (rval, CLEANUP);
	}
	if (fbid < bid - PARAM_DUAL_RESOLVEGAP)
	{
		//ILL_IFTRACE("\t%s:%s:%d\n",__func__,__FILE__,__LINE__);
		ILLfct_compute_piz (lp);
		ILLfct_compute_dz (lp);
	}

	if (phase == DUAL_PHASEII)
	{
		//ILL_IFTRACE("\t%s:%s:%d\n",__func__,__FILE__,__LINE__);
		if (fbid < bid - PARAM_DUAL_RESOLVEGAP)
		{
			//ILL_IFTRACE("\t%s:%s:%d\n",__func__,__FILE__,__LINE__);
			ILLfct_compute_xbz (lp);
			CHECKRVALG (rval, CLEANUP);
			if (p != NULL)
			{
				//ILL_IFTRACE("\t%s:%s:%d\n",__func__,__FILE__,__LINE__);
				if (p->d_strategy == COMPLETE_PRICING)
				{
					//ILL_IFTRACE("\t%s:%s:%d\n",__func__,__FILE__,__LINE__);
					ILLprice_compute_primal_inf (lp, p, NULL, 0, DUAL_PHASEII);
				}
				else
				{
					//ILL_IFTRACE("\t%s:%s:%d\n",__func__,__FILE__,__LINE__);
					ILLprice_update_mpartial_price (lp, p, DUAL_PHASEII, ROW_PRICING);
				}
			}
		}
		//ILL_IFTRACE("\t%s:%s:%d\n",__func__,__FILE__,__LINE__);
		ILLfct_compute_dobj (lp);
		ILLfct_check_dfeasible (lp, fi, *dtol);
		ILLfct_check_pfeasible (lp, fi, *ptol);
	}
	else if (phase == DUAL_PHASEI)
	{
		ILLfct_check_dfeasible (lp, fi, *idtol);
		if (fi->dstatus != DUAL_FEASIBLE)
		{
			ILLfct_compute_phaseI_xbz (lp);
			ILLfct_check_pIpfeasible (lp, fi, *iptol);
			if (p != NULL)
			{
				if (p->d_strategy == COMPLETE_PRICING)
					ILLprice_compute_primal_inf (lp, p, NULL, 0, DUAL_PHASEI);
				else
					ILLprice_update_mpartial_price (lp, p, DUAL_PHASEI, ROW_PRICING);
			}
		}
	}
CLEANUP:
	EG_RETURN (rval);
}

int ILLsimplex_solution (
	lpinfo * lp,
	EGlpNum_t * xz,
	EGlpNum_t * piz,
	EGlpNum_t * dz,
	EGlpNum_t * objval)
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
			EGlpNumCopy (xz[lp->baz[i]], lp->xbz[i]);
		for (j = 0; j < lp->nnbasic; j++)
		{
			col = lp->nbaz[j];
			if (lp->vstat[col] == STAT_UPPER)
				EGlpNumCopy (xz[col], lp->uz[col]);
			else if (lp->vstat[col] == STAT_LOWER)
				EGlpNumCopy (xz[col], lp->lz[col]);
			else
				EGlpNumZero (xz[col]);
		}
	}
	if (piz != NULL)
	{
		if (lp->basisstat.optimal == 0)
		{
			EG_RETURN (1);
		}
		for (i = 0; i < lp->nrows; i++)
			EGlpNumCopy (piz[i], lp->piz[i]);
	}
	if (dz != NULL)
	{
		if (lp->basisstat.optimal == 0)
		{
			EG_RETURN (1);
		}
		for (i = 0; i < lp->nrows; i++)
			EGlpNumZero (dz[lp->baz[i]]);
		for (j = 0; j < lp->nnbasic; j++)
			EGlpNumCopy (dz[lp->nbaz[j]], lp->dz[j]);
	}
	if (objval != NULL)
		EGlpNumCopy (*objval, lp->objval);
	return 0;
}

int ILLsimplex_infcertificate (
	lpinfo * lp,
	EGlpNum_t * pi)
{
	int i, col, nz;
	char *sense;
	EGlpNum_t *x, *l, *u;
	lp_status_info *ls;

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
			EGlpNumCopy (pi[i], lp->pIpiz[i]);
	}
	else if (lp->final_phase == DUAL_PHASEII && lp->infub_ix != -1)
	{
		col = lp->baz[lp->infub_ix];
		x = &(lp->xbz[lp->infub_ix]);
		l = &(lp->lz[col]);
		u = &(lp->uz[col]);

		for (i = 0; i < lp->nrows; i++)
			EGlpNumZero (pi[i]);

		if (EGlpNumIsNeqq (*l, NINFTY) && EGlpNumIsLess (*x, *l))
		{
			for (i = 0, nz = lp->zz.nzcnt; i < nz; i++)
				EGlpNumCopyNeg (pi[lp->zz.indx[i]], lp->zz.coef[i]);
		}
		else
		{
			for (i = 0, nz = lp->zz.nzcnt; i < nz; i++)
				EGlpNumCopy (pi[lp->zz.indx[i]], lp->zz.coef[i]);
		}
	}
	else
	{
		fprintf (stderr, "Invalid call to inf. certificate routine\n");
		EG_RETURN (1);
	}

	sense = lp->O->sense;
	for (i = 0; i < lp->nrows; i++)
	{
		if (sense[i] == 'G' && EGlpNumIsLessZero (pi[i]))
			EGlpNumZero (pi[i]);
		if (sense[i] == 'L' && EGlpNumIsGreatZero (pi[i]))
			EGlpNumZero (pi[i]);
	}
	return 0;
}

#if SIMPLEX_DEBUG > 1
static void test_cert (
	lpinfo * lp,
	EGlpNum_t * pi)
{
	int i, j;
	int mcnt, mbeg;
	EGlpNum_t fsum, sum;

	EGlpNumInitVar (fsum);
	EGlpNumInitVar (sum);
	EGlpNumZero (fsum);

	for (i = 0; i < lp->nrows; i++)
	{
		if (lp->O->sense[i] == 'G' && EGlpNumIsLessZero (pi[i]))
			printf ("compl \n");
		if (lp->O->sense[i] == 'L' && EGlpNumIsGreatZero (pi[i]))
			printf ("compll \n");
	}

	for (i = 0; i < lp->nrows; i++)
		EGlpNumAddInnProdTo (fsum, pi[i], lp->bz[i]);

	for (j = 0; j < lp->nnbasic; j++)
	{
		EGlpNumZero (sum);
		mcnt = lp->matcnt[j];
		mbeg = lp->matbeg[j];
		for (i = 0; i < mcnt; i++)
			EGlpNumAddInnProdTo (sum, pi[lp->matind[mbeg + i]], lp->matval[mbeg + i]);

		if (EGlpNumIsLess (PFEAS_TOLER, sum) &&
				(lp->vtype[j] == VLOWER || lp->vtype[j] == VFREE))
			printf ("compl2\n");
		else
		{
			EGlpNumSign (sum);
			if (EGlpNumIsLess (PFEAS_TOLER, sum) &&
					(lp->vtype[j] == VUPPER || lp->vtype[j] == VFREE))
				printf ("compl1\n");
			EGlpNumSign (sum);
		}

		if (EGlpNumIsLessZero (sum)
				&& (lp->vtype[j] & (VFREE | VUPPER)) == 0)
			EGlpNumSubInnProdTo (fsum, sum, lp->lz[j]);
		else if (EGlpNumIsGreatZero (sum)
						 && (lp->vtype[j] & (VFREE | VLOWER)) == 0)
			EGlpNumSubInnProdTo (fsum, sum, lp->uz[j]);
	}
	printf ("fsum = %.8f\n", EGlpNumToLf (fsum));
	EGlpNumClearVar (fsum);
	EGlpNumClearVar (sum);
}
#endif

static void save_paraminfo (
	price_info * pinf,
	iter_info * it)
{
	param_info *pr = &(it->oldinfo);

	pr->origalgo = it->algorithm;
	pr->pphaseI = pinf->pI_price;
	pr->pphaseII = pinf->pII_price;
	pr->dphaseI = pinf->dI_price;
	pr->dphaseII = pinf->dII_price;
	pr->p_strategy = pinf->p_strategy;
	pr->d_strategy = pinf->d_strategy;
}

static void restore_paraminfo (
	iter_info * it,
	price_info * pinf)
{
	param_info *pr = &(it->oldinfo);

	it->algorithm = pr->origalgo;
	pinf->pI_price = pr->pphaseI;
	pinf->pII_price = pr->pphaseII;
	pinf->dI_price = pr->dphaseI;
	pinf->dII_price = pr->dphaseII;
	pinf->p_strategy = pr->p_strategy;
	pinf->d_strategy = pr->d_strategy;
}

int ILLsimplex (
	lpinfo * lp,
	int algorithm,
	ILLlp_basis * B,
	price_info * pinf,
	int *status,
	int sdisplay,
	itcnt_t*itcnt)
{
	int phase = -1;
	int singular = -1;
	int rval = 0;
	int new_price = -1;
	svector wz;
	svector updz;
	feas_info fi;
	iter_info it;

	EGlpNumInitVar (fi.totinfeas);
	EGlpNumInitVar (it.prevobj);
	EGlpNumInitVar (it.objtol);

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
	EGlpNumCopy (it.prevobj, INFTY);
	it.nosolve = 0;
	it.noprog = 0;
	EGlpNumCopy (it.objtol, OBJBND_TOLER);
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

	free_internal_lpinfo (lp);
	init_internal_lpinfo (lp);
	rval = build_internal_lpinfo (lp);
	CHECKRVALG (rval, CLEANUP);

	ILLsvector_init (&wz);
	rval = ILLsvector_alloc (&wz, lp->nrows);
	CHECKRVALG (rval, CLEANUP);
	ILLsvector_init (&updz);
	rval = ILLsvector_alloc (&updz, lp->nrows);
	CHECKRVALG (rval, CLEANUP);

	if (it.sdisplay)
	{
		char buffer[256];
		int nonzero = 0;
		register int i = lp->ncols;

		while (i--)
			nonzero += lp->matcnt[i];
		sprintf (buffer, "starting ILLsimplex on %s...\n", lp->O->probname);
		/* depending on LP's reporter 
		 * string is printed to stdout 
		 * or handed to GUI */
		rval = rval || ILLstring_report (buffer, &lp->O->reporter);
		printf ("Problem has %d rows and %d cols and %d nonzeros\n", lp->nrows,
						lp->ncols, nonzero);
		fflush (stdout);
	}
	ILLfct_set_variable_type (lp);

	if (B != 0)
	{
		rval = ILLbasis_load (lp, B);
		CHECKRVALG (rval, CLEANUP);
		if (it.algorithm == DUAL_SIMPLEX)
		{
			if (B->rownorms)
			{
				rval = ILLprice_load_rownorms (lp, B->rownorms, pinf);
				CHECKRVALG (rval, CLEANUP);
			}
			else
				EGlpNumFreeArray (pinf->dsinfo.norms);
		}
		else if (it.algorithm == PRIMAL_SIMPLEX)
		{
			if (B->colnorms)
			{
				rval = ILLprice_load_colnorms (lp, B->colnorms, pinf);
				CHECKRVALG (rval, CLEANUP);
			}
			else
				EGlpNumFreeArray (pinf->psinfo.norms);
		}
		else if (it.algorithm != PRIMAL_OR_DUAL)
		{
			fprintf (stderr, "Unknown algorithm %d in ILLsimplex\n", it.algorithm);
			rval = 1;
			ILL_CLEANUP;
		}
	}
	else if (lp->basisid == -1)
	{
		if (lp->nrows < 200 && lp->ncols < 400)
			rval = ILLbasis_get_initial (lp, it.algorithm);
		else
			rval = ILLbasis_get_cinitial (lp, it.algorithm);
		CHECKRVALG (rval, CLEANUP);
		ILLprice_free_pricing_info (pinf);
	}

	if (lp->fbasisid != lp->basisid)
	{
		rval = ILLbasis_factor (lp, &singular);
		CHECKRVALG (rval, CLEANUP);
		if (singular)
		{
			MESSAGE (__QS_SB_VERB, "Singular basis found!");
			ILLprice_free_pricing_info (pinf);
		}
	}

START:
#if 0
	if (it.resumeid == SIMPLEX_RESUME_UNSHIFT)
		fprintf (stderr, "Resuming Unshift\n");
	else if (it.resumeid == SIMPLEX_RESUME_SING)
		fprintf (stderr, "Resuming Singular\n");
	else if (it.resumeid == SIMPLEX_RESUME_NUMER)
		fprintf (stderr, "Resuming Numer\n");
	else if (it.resumeid != -1)
		fprintf (stderr, "Resuming for other reason... %d\n", it.resumeid);
#endif
	it.solstatus = ILL_LP_UNSOLVED;
	init_lp_status_info (&(lp->basisstat));

	ILLfct_compute_piz (lp);
	ILLfct_compute_dz (lp);
	if (it.algorithm == DUAL_SIMPLEX)
	{
		if (B != NULL || it.resumeid == SIMPLEX_RESUME_UNSHIFT)
			ILLfct_dual_adjust (lp, lp->tol->dfeas_tol);
		else
			ILLfct_dual_adjust (lp, zeroLpNum);
	}
	ILLfct_compute_xbz (lp);

	ILLfct_check_pfeasible (lp, &fi, lp->tol->pfeas_tol);
	ILLfct_check_dfeasible (lp, &fi, lp->tol->dfeas_tol);
	ILLfct_set_status_values (lp, fi.pstatus, fi.dstatus, PHASEII, PHASEII);

	if (fi.dstatus == DUAL_FEASIBLE && EGlpNumIsNeqq (lp->objbound, INFTY))
	{
		ILLfct_compute_dobj (lp);
		if (EGlpNumIsLess (lp->objbound, lp->dobjval))
		{
			it.solstatus = ILL_BND_REACHED;
			printf ("solstatus = ILL_BND_REACHED = 5 %lf %lf\n",
							EGlpNumToLf (lp->objbound), EGlpNumToLf (lp->dobjval));
			goto TERMINATE;
		}
	}
	if (fi.pstatus == PRIMAL_FEASIBLE && fi.dstatus == DUAL_FEASIBLE)
	{
		it.solstatus = ILL_LP_SOLVED;
		ILLfct_compute_pobj (lp);
		goto TERMINATE;
	}

	if (it.algorithm == PRIMAL_OR_DUAL)
	{
		if (fi.pstatus == PRIMAL_FEASIBLE)
			it.algorithm = PRIMAL_SIMPLEX;
		else if (fi.dstatus == DUAL_FEASIBLE)
			it.algorithm = DUAL_SIMPLEX;
		else if (EGlpNumToLf (lp->pinfeas) < 10 * EGlpNumToLf (lp->dinfeas))
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

	rval = ILLprice_build_pricing_info (lp, pinf, phase);
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
			ILLprice_free_pricing_info (pinf);
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
				ILLfct_unroll_bound_change (lp);
				ILLfct_unroll_coef_change (lp);
				/* we are disabling re-do under this circunstances ! */
				rval = ILLbasis_get_initial (lp, it.algorithm);
				CHECKRVALG (rval, CLEANUP);
				rval = ILLbasis_factor (lp, &singular);
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
				new_price = ILLprice_get_price (pinf, phase);

				if (pinf->cur_price != new_price)
				{
					//ILL_IFTRACE("\t%s:%s:%d\n",__func__,__FILE__,__LINE__);
					ILLprice_free_pricing_info (pinf);
					rval = ILLprice_build_pricing_info (lp, pinf, phase);
					CHECKRVALG (rval, CLEANUP);
				}
			}
		}
	}

#if SIMPLEX_DEBUG > 0
	ILLfct_print_counts (lp);
#endif
LIMIT_TERMINATE:
	rval = terminate_simplex (lp, phase, &it);
	CHECKRVALG (rval, CLEANUP);

TERMINATE:
	restore_paraminfo (&it, pinf);

	if (it.sdisplay)
	{
		printf ("completed ILLsimplex\n");
		printf ("%s: ", lp->O->probname);
		fflush (stdout);
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
						fprintf (stdout, "Primal Infeasible\n");
					else
						fprintf (stdout, "Dual Unbounded\n");
				}
			}
			else if (lp->basisstat.primal_unbounded)
			{
				*status = QS_LP_UNBOUNDED;
			}
		}
		else
		{
			fprintf (stderr, "unknown solution status in ILLsimplex %d\n",
							 it.solstatus);
			rval = 1;
			CHECKRVALG (rval, CLEANUP);
		}
	}

#if SIMPLEX_DEBUG > 1
	{
		int rva = 0;
		EGlpNum_t *pi = NULL;

		pi = EGlpNumAllocArray (lp->nrows);
		rva = ILLsimplex_infcertificate (lp, pi);
		printf ("rva = %d\n", rva);
		if (!rva)
		{
			test_cert (lp, pi);
		}
		EGlpNumFreeArray (pi);
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

		printf ("time = %.3f, pI = %d, pII = %d, dI = %d, dII = %d, ",
						ILLutil_zeit () - lp->starttime, lp->cnts->pI_iter,
						lp->cnts->pII_iter, lp->cnts->dI_iter, lp->cnts->dII_iter);
		fflush (stdout);
		get_current_stat (&(lp->basisstat), it.algorithm, &bstat);
		switch (bstat)
		{
		case OPTIMAL:
			printf ("opt = %f\n", EGlpNumToLf (lp->objval));
			break;
		case PRIMAL_INFEASIBLE:
			printf ("no primal soln\n");
			break;
		case PRIMAL_UNBOUNDED:
			printf ("primal unbounded\n");
			break;
		case PRIMAL_FEASIBLE:
			printf ("primal obj = %f\n", EGlpNumToLf (lp->pobjval));
			break;
		case DUAL_INFEASIBLE:
			printf ("no dual soln\n");
			break;
		case DUAL_UNBOUNDED:
			printf ("dual unbounded\n");
			break;
		case DUAL_FEASIBLE:
			printf ("dual obj = %f\n", EGlpNumToLf (lp->dobjval));
			break;
		}
		fflush (stdout);

		if (it.sdisplay > 1)
		{
			if (it.algorithm == PRIMAL_SIMPLEX && pinf->pI_price == QS_PRICE_PDEVEX)
				printf ("Devex norms initialised %d times\n", pinf->pdinfo.ninit);
			fflush (stdout);
		}
	}

CLEANUP:
	ILLsvector_free (&wz);
	ILLsvector_free (&updz);
	EGlpNumClearVar (it.prevobj);
	EGlpNumClearVar (it.objtol);
	EGlpNumClearVar (fi.totinfeas);
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
	lpinfo * lp,
	int phase,
	iter_info * it)
{
	int rval = 0;
	int sphase;
	feas_info fi;

	EGlpNumInitVar (fi.totinfeas);

	if (it->solstatus != ILL_MAX_TIME && it->solstatus != ILL_MAX_ITER)
		ILL_CLEANUP;

	if (it->algorithm == PRIMAL_SIMPLEX)
	{
		if (lp->nbchange != 0)
		{
			if (it->sdisplay > 1)
			{
				printf ("unrolling %d bound shifts\n", lp->nbchange);
				fflush (stdout);
			}
			ILLfct_unroll_bound_change (lp);
		}
		rval = ILLsimplex_retest_psolution (lp, NULL, phase, &fi);
		CHECKRVALG (rval, CLEANUP);

		sphase = (phase == PRIMAL_PHASEI) ? PHASEI : PHASEII;
		ILLfct_set_status_values (lp, fi.pstatus, fi.dstatus, PHASEII, sphase);
	}
	else if (it->algorithm == DUAL_SIMPLEX)
	{
		if (lp->ncchange != 0)
		{
			if (it->sdisplay > 1)
			{
				printf ("unrolling %d coef shifts\n", lp->ncchange);
				fflush (stdout);
			}
			ILLfct_unroll_coef_change (lp);
		}
		rval = ILLsimplex_retest_dsolution (lp, NULL, phase, &fi);
		CHECKRVALG (rval, CLEANUP);

		sphase = (phase == DUAL_PHASEI) ? PHASEI : PHASEII;
		ILLfct_set_status_values (lp, fi.pstatus, fi.dstatus, sphase, PHASEII);
	}

CLEANUP:
	EGlpNumClearVar (fi.totinfeas);
	EG_RETURN (rval);
}

static int test_progress (
	EGlpNum_t objval,
	EGlpNum_t prevobj)
{
	EGlpNum_t denom;

	EGlpNumInitVar (denom);
	EGlpNumCopyDiff (denom, objval, prevobj);
	if (EGlpNumIsNeqZero (objval, PROGRESS_ZERO))
		EGlpNumDivTo (denom, objval);
	if (!EGlpNumIsNeqZero (denom, PROGRESS_THRESH))
	{
		EGlpNumClearVar (denom);
		return 0;
	}
	else
	{
		EGlpNumClearVar (denom);
		return 1;
	}
}

static void monitor_iter (
	lpinfo * lp,
	price_info * p,
	iter_info * it,
	int phase)
{
	EGlpNum_t print_val;
	double tottime = ILLutil_zeit () - lp->starttime;
	int curtime = ILLutil_our_floor (tottime);	/* MONIKA */
	char print_str[20];
	feas_info fi;
	int aborted = 0;

	EGlpNumInitVar (print_val);
	EGlpNumInitVar (fi.totinfeas);
	EGlpNumZero (fi.totinfeas);
	EGlpNumZero (print_val);

	/* one of the following two time display mechanisms */
	switch (phase)
	{
	case PRIMAL_PHASEI:
		EGlpNumCopy (print_val, lp->pinfeas);
		EGlpNumAddTo (fi.totinfeas, lp->pinfeas);
		strcpy (print_str, "primal infeas");
		if (EGlpNumIsLessZero (lp->pinfeas) &&
				(EGlpNumIsNeqZero (lp->pinfeas, oneLpNum)))
		{
			/*printf ("Negative Infeasibility! Imposible %lg %la, iter %d\n",
			 * EGlpNumToLf (print_val), EGlpNumToLf (print_val), it->itercnt);
			 */
			//exit(1);
		}
		break;
	case PRIMAL_PHASEII:
		EGlpNumCopy (print_val, lp->pobjval);
		strcpy (print_str, "primal objval");
		break;
	case DUAL_PHASEI:
		EGlpNumCopy (print_val, lp->dinfeas);
		EGlpNumAddTo (fi.totinfeas, lp->dinfeas);
		strcpy (print_str, "dual infeas");
		break;
	case DUAL_PHASEII:
		EGlpNumCopy (print_val, lp->dobjval);
		strcpy (print_str, "dual objval");
		break;
	}

	aborted = report_value (lp, it, print_str, print_val);
	/*if (it->sdisplay && it->itercnt % lp->iterskip == 0) {
	 * // printf ("(%d): %s = %f\n", it->itercnt, print_str, print_val);
	 * // fflush (stdout);
	 * } */
	if (curtime != it->curtime)
	{
		it->curtime = curtime;
		/*
		 * if (it->sdisplay){
		 * printf ("time = %d.0, ", curtime);
		 * printf ("(%d): %s = %f\n", it->itercnt, print_str, print_val);
		 * fflush (stdout);
		 * }
		 */
	}

	EGlpNumAddUiTo (fi.totinfeas, 1000);
	if (EGlpNumIsLessZero (fi.totinfeas))
	{
		it->nextstep = SIMPLEX_TERMINATE;
		it->solstatus = ILL_MAX_ITER;
		MESSAGE (it->sdisplay ? 0 : __QS_SB_VERB,
						 "early finish by excess infeasibility");
		ILL_CLEANUP;
	}

	if (phase == DUAL_PHASEII && EGlpNumIsNeqq (lp->objbound, INFTY))
	{
		/*if (lp->dobjval > lp->objbound + it->objtol) */
		EGlpNumCopyDiff (print_val, lp->dobjval, lp->objbound);
		if (EGlpNumIsLess (it->objtol, print_val))
		{
			ILLfct_unroll_coef_change (lp);
			ILLfct_check_dfeasible (lp, &fi, lp->tol->dfeas_tol);
			ILLfct_set_status_values (lp, -1, fi.dstatus, -1, PHASEII);

			if (fi.dstatus == DUAL_FEASIBLE)
			{
				ILLfct_compute_dobj (lp);
				if (EGlpNumIsLess (lp->objbound, lp->dobjval))
				{
					it->solstatus = ILL_BND_REACHED;
					it->nextstep = SIMPLEX_TERMINATE;
					/*if (it->sdisplay) */
					{
						printf ("bound reached %lf %lf\n", EGlpNumToLf (lp->objbound),
										EGlpNumToLf (lp->dobjval));
						fflush (stdout);
					}
				}
				else
					EGlpNumMultUiTo (it->objtol, 10);
			}
			else
			{
				it->nextphase = DUAL_PHASEI;
				it->newphase = SIMPLEX_PHASE_NEW;
				EGlpNumMultUiTo (it->objtol, 5);
			}
		}
	}
	if (it->itercnt >= lp->maxiter)
	{
		it->solstatus = ILL_MAX_ITER;
		it->nextstep = SIMPLEX_TERMINATE;
		if (it->sdisplay)
		{
			printf ("iter limit reached\n");
			fflush (stdout);
		}
		ILL_CLEANUP;
	}
	else if (tottime >= lp->maxtime)
	{
		it->solstatus = ILL_MAX_TIME;
		it->nextstep = SIMPLEX_TERMINATE;
		if (it->sdisplay)
		{
			printf ("time limit reached\n");
			fflush (stdout);
		}
		ILL_CLEANUP;
	}
	else if (aborted)
	{
		it->solstatus = ILL_LP_ABORTED;
		it->nextstep = SIMPLEX_TERMINATE;
		if (it->sdisplay)
		{
			printf ("aborted\n");
			fflush (stdout);
		}
		ILL_CLEANUP;
	}
	/* why is this commented out? */
	if(0){
		if (it->rounds && it->inner){
			it->inner --;
			if (it->inner == 0){
				printf ("restoring ..\n");
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
			ILLfct_perturb_coefs (lp);
			it->noprog = 0;
			EGlpNumCopy (it->prevobj, lp->dobjval);
		}
	}
	else if (phase == PRIMAL_PHASEII)
	{
		if (it->noprog > it->chkobj)
		{
			ILLfct_perturb_bounds (lp);
			it->noprog = 0;
			EGlpNumCopy (it->prevobj, lp->pobjval);
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
			EGlpNumMultUiTo (lp->tol->pfeas_tol, SIMPLEX_FACTOR);
			EGlpNumMultUiTo (lp->tol->dfeas_tol, SIMPLEX_FACTOR);
			it->n_restart++;
			//fprintf(stderr,"Resume Numerical %s:%s:%d\n",__func__,__FILE__,__LINE__);
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
			EGlpNumMultUiTo (lp->tol->pfeas_tol, SIMPLEX_FACTOR);
			EGlpNumMultUiTo (lp->tol->dfeas_tol, SIMPLEX_FACTOR);
			it->n_restart++;
			//fprintf(stderr,"Resume Numerical %s:%s:%d\n",__func__,__FILE__,__LINE__);
		}
	}
CLEANUP:
	EGlpNumClearVar (fi.totinfeas);
	EGlpNumClearVar (print_val);
	return;
}

static int primal_phaseI_step (
	lpinfo * lp,
	price_info * pinf,
	svector * updz,
	svector * wz,
	iter_info * it)
{
	int rval = 0;
	int singular = 0;
	int refactor = 0;
	int cphase = PRIMAL_PHASEI;
	EGlpNum_t alpha;
	feas_info fi;
	ratio_res rs;
	price_res pr;

	EGlpNumInitVar (alpha);
	EGlpNumInitVar (fi.totinfeas);
	EGlpNumInitVar (pr.dinfeas);
	EGlpNumInitVar (pr.pinfeas);
	EGlpNumInitVar (rs.tz);
	EGlpNumInitVar (rs.lbound);
	EGlpNumInitVar (rs.ecoeff);
	EGlpNumInitVar (rs.pivotval);
	EGlpNumZero (alpha);

	ILLfct_update_counts (lp, CNT_PPHASE1ITER, 0, zeroLpNum);
	it->nextstep = SIMPLEX_CONTINUE;
	it->nextphase = PRIMAL_PHASEI;
	lp->final_phase = PRIMAL_PHASEI;
	it->nosolve++;

	if (it->newphase != 0)
	{
		ILLfct_check_pfeasible (lp, &fi, lp->tol->ip_tol);
		if (it->newphase == SIMPLEX_PHASE_NEW)
		{
			it->noprog = 0;
			if (it->sdisplay)
			{
				printf ("starting primal phase I, nosolve %d\n", it->nosolve);
				fflush (stdout);
			}
		}
		it->newphase = 0;
		it->nosolve = 0;
		EGlpNumCopy (it->prevobj, lp->pinfeas);
		lp->pIpiz = EGlpNumAllocArray (lp->nrows);
		lp->pIdz = EGlpNumAllocArray (lp->nnbasic);

		ILLfct_compute_phaseI_piz (lp);
		if (pinf->p_strategy == COMPLETE_PRICING)
		{
			ILLfct_compute_phaseI_dz (lp);
#if USEHEAP > 0
			ILLprice_free_heap (pinf);
#endif
			ILLprice_compute_dual_inf (lp, pinf, NULL, 0, PRIMAL_PHASEI);
#if USEHEAP > 0
			rval = ILLprice_test_for_heap (lp, pinf, lp->nnbasic, pinf->d_scaleinf,
																		 PRIMAL_SIMPLEX, 0);
			CHECKRVALG (rval, CLEANUP);
#endif
		}
		else if (pinf->p_strategy == MULTI_PART_PRICING)
			ILLprice_init_mpartial_price (lp, pinf, cphase, COL_PRICING);
	}

	monitor_iter (lp, pinf, it, cphase);
	if (it->nextstep == SIMPLEX_TERMINATE || it->nextstep == SIMPLEX_RESUME ||
			it->newphase != 0)
		ILL_CLEANUP;

	ILLprice_primal (lp, pinf, &pr, cphase);
	ILL_IFTRACE2 ("%s:after_price\n", __func__);

	if (pr.price_stat == PRICE_OPTIMAL)
	{
		if (it->sdisplay > 1)
		{
			printf ("primal phase I seemingly done\n");
			printf ("retesting soln\n");
			fflush (stdout);
		}
		rval = ILLsimplex_retest_psolution (lp, pinf, cphase, &fi);

		CHECKRVALG (rval, CLEANUP);
		ILLfct_set_status_values (lp, fi.pstatus, fi.dstatus, PHASEII, PHASEI);

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

	ILLfct_compute_yz (lp, &(lp->yjz), updz, lp->nbaz[pr.eindex]);
	ILLfct_update_counts (lp, CNT_YNZ, lp->yjz.nzcnt, zeroLpNum);
	ILLfct_update_counts (lp, CNT_UPNZ, updz->nzcnt, zeroLpNum);

	ILLratio_pI_test (lp, pr.eindex, pr.dir, &rs);
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
		EGlpNumMultUiTo (lp->tol->pfeas_tol, SIMPLEX_FACTOR);
		EGlpNumMultUiTo (lp->tol->dfeas_tol, SIMPLEX_FACTOR);
		//fprintf(stderr,"Resume Numerical %s:%s:%d\n",__func__,__FILE__,__LINE__);
		ILL_CLEANUP;
	}
	else if (rs.ratio_stat == RATIO_NEGATIVE)
	{
		EGlpNum_t itol;

		EGlpNumInitVar (itol);
		//ILL_IFTRACE("ratio_negative\n");
		EGlpNumCopy (itol, lp->tol->ip_tol);
		EGlpNumZero (lp->tol->ip_tol);
		EGlpNumAddTo (lp->pinfeas, lp->upd.c_obj);
		if (!test_progress (lp->pinfeas, it->prevobj))
			it->noprog++;
		else
		{
			EGlpNumCopy (it->prevobj, lp->pinfeas);
			it->noprog = 0;
		}
		ILLfct_update_pfeas (lp, rs.lindex, &(lp->srhs));
		EGlpNumCopy (lp->tol->ip_tol, itol);
		ILLfct_compute_ppIzz (lp, &(lp->srhs), &(lp->ssoln));
		ILLfct_update_ppI_prices (lp, pinf, &(lp->srhs), &(lp->ssoln), pr.eindex,
															rs.lindex, zeroLpNum);
		EGlpNumClearVar (itol);
	}
	else if (rs.ratio_stat == RATIO_NOBCHANGE)
	{
		//ILL_IFTRACE("ratio_nobchange\n");
		EGlpNumAddTo (lp->pinfeas, lp->upd.c_obj);
		if (!test_progress (lp->pinfeas, it->prevobj))
			it->noprog++;
		else
		{
			EGlpNumCopy (it->prevobj, lp->pinfeas);
			it->noprog = 0;
		}

		//ILL_IFTRACE("%s:a\n",__func__);
		ILLfct_update_xz (lp, rs.tz, pr.eindex, rs.lindex);
		ILLfct_update_pfeas (lp, rs.lindex, &(lp->srhs));
		ILLfct_compute_ppIzz (lp, &(lp->srhs), &(lp->ssoln));
		ILLfct_update_basis_info (lp, pr.eindex, rs.lindex, rs.lvstat);
#if DENSE_PI > 0
		fct_test_workvector (lp);
		fct_test_pfeasible (lp);
#endif
		ILLfct_update_ppI_prices (lp, pinf, &(lp->srhs), &(lp->ssoln), pr.eindex,
															rs.lindex, zeroLpNum);
	}
	else if (rs.ratio_stat == RATIO_BCHANGE)
	{
		//ILL_IFTRACE("ratio_bchange\n");
		EGlpNumCopyFrac (alpha, lp->pIdz[pr.eindex], rs.pivotval);
		EGlpNumAddTo (lp->pinfeas, lp->upd.c_obj);

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
			EGlpNumCopy (it->prevobj, lp->pinfeas);
			it->noprog = 0;
		}

		ILLfct_compute_zz (lp, &(lp->zz), rs.lindex);
		ILLfct_update_counts (lp, CNT_ZNZ, lp->zz.nzcnt, zeroLpNum);
		if (pinf->p_strategy == COMPLETE_PRICING)
		{
			//ILL_IFTRACE("%s:a\n",__func__);
			ILLfct_compute_zA (lp, &(lp->zz), &(lp->zA));
			ILLfct_update_counts (lp, CNT_ZANZ, lp->zA.nzcnt, zeroLpNum);

			if (pinf->pI_price == QS_PRICE_PSTEEP)
			{
				ILLfct_compute_psteep_upv (lp, wz);
			}
		}

		rval =
			ILLprice_update_pricing_info (lp, pinf, cphase, wz, pr.eindex, rs.lindex,
																		rs.pivotval);
		CHECKRVALG (rval, CLEANUP);

		//ILL_IFTRACE("%s:b:%d\n",__func__,rs.lindex);
		ILLfct_update_xz (lp, rs.tz, pr.eindex, rs.lindex);
		ILLfct_update_pfeas (lp, rs.lindex, &(lp->srhs));
		//ILL_IFTRACE("%s:%d:%d\n",__func__,rs.lindex,lp->srhs.nzcnt);
		ILLfct_compute_ppIzz (lp, &(lp->srhs), &(lp->ssoln));
		ILLfct_update_basis_info (lp, pr.eindex, rs.lindex, rs.lvstat);
#if DENSE_PI > 0
		fct_test_workvector (lp);
		fct_test_pfeasible (lp);
#endif
		rval = ILLbasis_update (lp, updz, rs.lindex, &refactor, &singular);
		CHECKRVALG (rval, CLEANUP);

		if (singular)
		{
			it->nextstep = SIMPLEX_RESUME;
			it->resumeid = SIMPLEX_RESUME_SING;
			/* this is to force to exit in the case of bad basis */
			it->n_restart++;
			EGlpNumMultUiTo (lp->tol->pfeas_tol, SIMPLEX_FACTOR);
			EGlpNumMultUiTo (lp->tol->dfeas_tol, SIMPLEX_FACTOR);
			//fprintf(stderr,"Resume Singular %s:%s:%d\n",__func__,__FILE__,__LINE__);
			ILL_CLEANUP;
		}
		if (!refactor)
		{
			ILLfct_update_ppI_prices (lp, pinf, &(lp->srhs), &(lp->ssoln), pr.eindex,
																rs.lindex, alpha);
		}
		if (refactor != 0 || it->nosolve > PARAM_MAX_NOSOLVE)
		{
			ILLfct_compute_xbz (lp);
			ILLfct_check_pfeasible (lp, &fi, lp->tol->ip_tol);
			ILLfct_set_status_values (lp, fi.pstatus, -1, PHASEII, -1);
			if (fi.pstatus == PRIMAL_FEASIBLE)
				it->nextphase = PRIMAL_PHASEII;

			it->newphase = SIMPLEX_PHASE_RECOMP;
			ILL_CLEANUP;
		}
	}

#if DENSE_PI > 1
	fct_test_workvector (lp);
	fct_test_pi_dz (lp, pinf);
#endif

CLEANUP:
	if (it->nextphase != PRIMAL_PHASEI || it->nextstep == SIMPLEX_RESUME ||
			it->newphase != 0 || rval != 0)
	{
		EGlpNumFreeArray (lp->pIpiz);
		EGlpNumFreeArray (lp->pIdz);
	}
	EGlpNumClearVar (alpha);
	EGlpNumClearVar (fi.totinfeas);
	EGlpNumClearVar (pr.dinfeas);
	EGlpNumClearVar (pr.pinfeas);
	EGlpNumClearVar (rs.tz);
	EGlpNumClearVar (rs.lbound);
	EGlpNumClearVar (rs.ecoeff);
	EGlpNumClearVar (rs.pivotval);
	//EG_RETURN(rval);
	return rval;
}

static int primal_phaseII_step (
	lpinfo * lp,
	price_info * pinf,
	svector * updz,
	svector * wz,
	iter_info * it)
{
	int boundch;
	int rval = 0;
	int bndtype = 0;
	int singular = 0;
	int refactor = 0;
	int ratio_iter = 0;
	int cphase = PRIMAL_PHASEII;
	EGlpNum_t lbound;
	EGlpNum_t alpha;
	feas_info fi;
	ratio_res rs;
	price_res pr;

	EGlpNumInitVar (alpha);
	EGlpNumInitVar (lbound);
	EGlpNumInitVar (fi.totinfeas);
	EGlpNumInitVar (pr.dinfeas);
	EGlpNumInitVar (pr.pinfeas);
	EGlpNumInitVar (rs.tz);
	EGlpNumInitVar (rs.lbound);
	EGlpNumInitVar (rs.ecoeff);
	EGlpNumInitVar (rs.pivotval);

	ILLfct_update_counts (lp, CNT_PPHASE2ITER, 0, zeroLpNum);
	it->nextstep = SIMPLEX_CONTINUE;
	it->nextphase = PRIMAL_PHASEII;
	lp->final_phase = PRIMAL_PHASEII;
	it->nosolve++;

	if (it->newphase != 0)
	{
		ILLfct_compute_pobj (lp);
		if (it->newphase == SIMPLEX_PHASE_NEW)
		{
			it->noprog = 0;
			if (it->sdisplay)
			{
				printf ("starting primal phase II, nosolve %d\n", it->nosolve);
				fflush (stdout);
			}
		}
		it->newphase = 0;
		it->nosolve = 0;
		EGlpNumCopy (it->prevobj, lp->pobjval);
		ILLfct_compute_piz (lp);
		if (pinf->p_strategy == COMPLETE_PRICING)
		{
			ILLfct_compute_dz (lp);
#if USEHEAP > 0
			ILLprice_free_heap (pinf);
#endif
			ILLprice_compute_dual_inf (lp, pinf, NULL, 0, PRIMAL_PHASEII);
#if USEHEAP > 0
			rval = ILLprice_test_for_heap (lp, pinf, lp->nnbasic, pinf->d_scaleinf,
																		 PRIMAL_SIMPLEX, 0);
			CHECKRVALG (rval, CLEANUP);
#endif
		}
		else if (pinf->p_strategy == MULTI_PART_PRICING)
		{
			ILLprice_init_mpartial_price (lp, pinf, cphase, COL_PRICING);
		}
	}

	monitor_iter (lp, pinf, it, cphase);
	if (it->nextstep == SIMPLEX_TERMINATE || it->nextstep == SIMPLEX_RESUME ||
			it->newphase != 0)
		ILL_CLEANUP;

	ILLprice_primal (lp, pinf, &pr, cphase);

	if (pr.price_stat == PRICE_OPTIMAL)
	{
		//ILL_IFTRACE("%s:PRICE_OPTIMAL\n",__func__);
		if (lp->nbchange != 0)
		{
			if (it->sdisplay > 1)
			{
				printf ("unrolling %d bound shifts\n", lp->nbchange);
				fflush (stdout);
			}
			ILLfct_unroll_bound_change (lp);
			ILLfct_check_pfeasible (lp, &fi, lp->tol->pfeas_tol);
			ILLfct_set_status_values (lp, fi.pstatus, -1, PHASEII, -1);

			 /*HHH*/ ILLfct_check_dfeasible (lp, &fi, lp->tol->dfeas_tol);
			/*HHH* printf ("primal (opt) infeas %.6f\n", lp->pinfeas); fflush (stdout); 
			 *HHH* printf ("dual (opt) infeas %.6f\n", lp->dinfeas); fflush (stdout);*/

			if (fi.pstatus != PRIMAL_FEASIBLE)
			{
				it->algorithm = DUAL_SIMPLEX;
				it->nextstep = SIMPLEX_RESUME;
				it->resumeid = SIMPLEX_RESUME_UNSHIFT;
				it->pricetype = QS_PRICE_DDEVEX;
				/* this is to force to exit in the case of bad basis */
				//fprintf(stderr,"Resume Unshift %s:%s:%d\n",__func__,__FILE__,__LINE__);
				EGlpNumMultUiTo (lp->tol->pfeas_tol, SIMPLEX_FACTOR);
				EGlpNumMultUiTo (lp->tol->dfeas_tol, SIMPLEX_FACTOR);
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
			printf ("problem seemingly solved\n");
			printf ("seemingly opt = %f\nretesting soln\n",
							EGlpNumToLf (lp->pobjval));
			fflush (stdout);
		}
		rval = ILLsimplex_retest_psolution (lp, pinf, cphase, &fi);
		CHECKRVALG (rval, CLEANUP);
		ILLfct_set_status_values (lp, fi.pstatus, fi.dstatus, PHASEII, PHASEII);

		if (fi.pstatus == PRIMAL_INFEASIBLE)
		{
			it->nextphase = PRIMAL_PHASEI;
			EGlpNumDivUiTo (lp->tol->ip_tol, 5);
			EGlpNumDivUiTo (lp->tol->id_tol, 5);
			ILL_IFTRACE ("%s:PINF:%lg\n", __func__, EGlpNumToLf (lp->tol->ip_tol));
		}
		else if (fi.dstatus == DUAL_FEASIBLE)
		{
			//ILL_IFTRACE("%s:PFEAS_DFEAS\n",__func__);
			it->solstatus = ILL_LP_SOLVED;
			EGlpNumCopy (lp->objval, lp->pobjval);
			it->nextstep = SIMPLEX_TERMINATE;
		}
		else
			ILL_IFTRACE ("%s:DINF:%la:%lf\n", __func__, EGlpNumToLf (lp->dinfeas),
									 EGlpNumToLf (lp->dinfeas));
		ILL_CLEANUP;
	}

	ILLfct_compute_yz (lp, &(lp->yjz), updz, lp->nbaz[pr.eindex]);
	ILLfct_update_counts (lp, CNT_YNZ, lp->yjz.nzcnt, zeroLpNum);
	ILLfct_update_counts (lp, CNT_UPNZ, updz->nzcnt, zeroLpNum);
	ratio_iter = 0;
	do
	{
		ILLratio_pII_test (lp, pr.eindex, pr.dir, &rs);
		//ILL_IFTRACE("all:%d",rs.lindex);
		EGlpNumCopy (lbound, rs.lbound);
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
			rval = ILLfct_bound_shift (lp, lp->baz[rs.lindex], bndtype, lbound);
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
		EGlpNumMultUiTo (lp->tol->pfeas_tol, SIMPLEX_FACTOR);
		EGlpNumMultUiTo (lp->tol->dfeas_tol, SIMPLEX_FACTOR);
		//fprintf(stderr,"Resume Numerical %s:%s:%d\n",__func__,__FILE__,__LINE__);
		ILL_CLEANUP;
	}
	else if (rs.ratio_stat == RATIO_UNBOUNDED)
	{
		//ILL_IFTRACE(":%d",rs.lindex);
		if (lp->nbchange != 0)
		{
			if (it->sdisplay > 1)
			{
				printf ("unrolling %d bound shifts\n", lp->nbchange);
				fflush (stdout);
			}
			ILLfct_unroll_bound_change (lp);
		}
		ILLfct_set_status_values (lp, PRIMAL_UNBOUNDED, -1, PHASEII, -1);
		it->solstatus = ILL_LP_SOLVED;
		it->nextstep = SIMPLEX_TERMINATE;
		ILL_CLEANUP;
	}
	else if (rs.ratio_stat == RATIO_NOBCHANGE)
	{
		//ILL_IFTRACE(":%d",rs.lindex);
		EGlpNumAddInnProdTo (lp->pobjval, rs.tz, lp->dz[pr.eindex]);
		EGlpNumCopy (lp->objval, lp->pobjval);
		if (!test_progress (lp->pobjval, it->prevobj))
			it->noprog++;
		else
		{
			EGlpNumCopy (it->prevobj, lp->pobjval);
			it->noprog = 0;
		}

		//ILL_IFTRACE("%s:c:%d\n",__func__,rs.lindex);
		ILLfct_update_xz (lp, rs.tz, pr.eindex, rs.lindex);
		ILLfct_update_basis_info (lp, pr.eindex, rs.lindex, rs.lvstat);
		if (pinf->p_strategy == COMPLETE_PRICING)
			ILLprice_compute_dual_inf (lp, pinf, &pr.eindex, 1, PRIMAL_PHASEII);
		else if (pinf->p_strategy == MULTI_PART_PRICING)
			ILLprice_update_mpartial_price (lp, pinf, cphase, COL_PRICING);
	}
	else if (rs.ratio_stat == RATIO_BCHANGE)
	{
		EGlpNumCopyFrac (alpha, lp->dz[pr.eindex], rs.pivotval);
		EGlpNumAddInnProdTo (lp->pobjval, rs.tz, lp->dz[pr.eindex]);
		EGlpNumCopy (lp->objval, lp->pobjval);

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
			EGlpNumCopy (it->prevobj, lp->pobjval);
			it->noprog = 0;
		}

		//ILL_IFTRACE(":%d",rs.lindex);
		ILLfct_compute_zz (lp, &(lp->zz), rs.lindex);
		ILLfct_update_counts (lp, CNT_ZNZ, lp->zz.nzcnt, zeroLpNum);
		if (pinf->p_strategy == COMPLETE_PRICING)
		{
			//ILL_IFTRACE("%s:b\n",__func__);
			ILLfct_compute_zA (lp, &(lp->zz), &(lp->zA));
			ILLfct_update_counts (lp, CNT_ZANZ, lp->zA.nzcnt, zeroLpNum);
			if (pinf->pII_price == QS_PRICE_PSTEEP)
				ILLfct_compute_psteep_upv (lp, wz);
		}
		rval =
			ILLprice_update_pricing_info (lp, pinf, cphase, wz, pr.eindex, rs.lindex,
																		rs.pivotval);
		CHECKRVALG (rval, CLEANUP);

		//ILL_IFTRACE("%s:d:%d\n",__func__,rs.lindex);
		ILLfct_update_xz (lp, rs.tz, pr.eindex, rs.lindex);
		ILLfct_update_basis_info (lp, pr.eindex, rs.lindex, rs.lvstat);
		rval = ILLbasis_update (lp, updz, rs.lindex, &refactor, &singular);
		CHECKRVALG (rval, CLEANUP);

		if (singular)
		{
			it->nextstep = SIMPLEX_RESUME;
			it->resumeid = SIMPLEX_RESUME_SING;
			/* this is to force to exit in the case of bad basis */
			it->n_restart++;
			//fprintf(stderr,"Resume Singular %s:%s:%d\n",__func__,__FILE__,__LINE__);
			EGlpNumMultUiTo (lp->tol->pfeas_tol, SIMPLEX_FACTOR);
			EGlpNumMultUiTo (lp->tol->dfeas_tol, SIMPLEX_FACTOR);
			ILL_CLEANUP;
		}
		if (!refactor)
		{
			ILLfct_update_piz (lp, alpha);

			if (pinf->p_strategy == COMPLETE_PRICING)
			{
				ILLfct_update_dz (lp, pr.eindex, alpha);
				ILLprice_compute_dual_inf (lp, pinf, lp->zA.indx, lp->zA.nzcnt,
																	 PRIMAL_PHASEII);
				ILLfct_update_counts (lp, CNT_ZARAVG, lp->zA.nzcnt, zeroLpNum);
			}
			else if (pinf->p_strategy == MULTI_PART_PRICING)
			{
				ILLprice_update_mpartial_price (lp, pinf, cphase, COL_PRICING);
			}
		}
		if (refactor != 0 || it->nosolve > PARAM_MAX_NOSOLVE)
		{
			ILLfct_compute_xbz (lp);
			it->newphase = SIMPLEX_PHASE_RECOMP;
		}
	}

CLEANUP:
	EGlpNumClearVar (alpha);
	EGlpNumClearVar (lbound);
	EGlpNumClearVar (fi.totinfeas);
	EGlpNumClearVar (pr.dinfeas);
	EGlpNumClearVar (pr.pinfeas);
	EGlpNumClearVar (rs.tz);
	EGlpNumClearVar (rs.lbound);
	EGlpNumClearVar (rs.ecoeff);
	EGlpNumClearVar (rs.pivotval);
	//EG_RETURN(rval);
	return rval;
}

static int dual_phaseI_step (
	lpinfo * lp,
	price_info * pinf,
	svector * updz,
	svector * wz,
	iter_info * it)
{
	int rval = 0;
	int singular = 0;
	int refactor = 0;
	int cphase = DUAL_PHASEI;
	EGlpNum_t alpha;
	EGlpNum_t alpha1;
	feas_info fi;
	ratio_res rs;
	price_res pr;

	EGlpNumInitVar (alpha);
	EGlpNumInitVar (alpha1);
	EGlpNumInitVar (fi.totinfeas);
	EGlpNumInitVar (pr.dinfeas);
	EGlpNumInitVar (pr.pinfeas);
	EGlpNumInitVar (rs.tz);
	EGlpNumInitVar (rs.lbound);
	EGlpNumInitVar (rs.ecoeff);
	EGlpNumInitVar (rs.pivotval);
	EGlpNumZero (alpha1);

	ILLfct_update_counts (lp, CNT_DPHASE1ITER, 0, zeroLpNum);
	it->nextstep = SIMPLEX_CONTINUE;
	it->nextphase = DUAL_PHASEI;
	lp->final_phase = DUAL_PHASEI;
	it->nosolve++;

	if (it->newphase != 0)
	{
		ILLfct_check_dfeasible (lp, &fi, lp->tol->id_tol);
		if (it->newphase == SIMPLEX_PHASE_NEW)
		{
			it->noprog = 0;
			if (it->sdisplay)
			{
				printf ("starting dual phase I, nosolve %d\n", it->nosolve);
				fflush (stdout);
			}
		}
		it->newphase = 0;
		it->nosolve = 0;
		EGlpNumCopy (it->prevobj, lp->dinfeas);

		ILLfct_compute_phaseI_xbz (lp);
		if (pinf->d_strategy == COMPLETE_PRICING)
		{
#if USEHEAP > 0
			ILLprice_free_heap (pinf);
#endif
			ILLprice_compute_primal_inf (lp, pinf, NULL, 0, DUAL_PHASEI);
#if USEHEAP > 0
			rval = ILLprice_test_for_heap (lp, pinf, lp->nrows, pinf->p_scaleinf,
																		 DUAL_SIMPLEX, 0);
			CHECKRVALG (rval, CLEANUP);
#endif
		}
		else if (pinf->d_strategy == MULTI_PART_PRICING)
		{
			ILLprice_init_mpartial_price (lp, pinf, cphase, ROW_PRICING);
		}
	}

	monitor_iter (lp, pinf, it, cphase);
	if (it->nextstep == SIMPLEX_TERMINATE || it->nextstep == SIMPLEX_RESUME ||
			it->newphase != 0)
		ILL_CLEANUP;

	ILLprice_dual (lp, pinf, cphase, &pr);

	if (pr.price_stat == PRICE_OPTIMAL)
	{
		if (it->sdisplay > 1)
		{
			printf ("dual phase I seemingly done\n");
			printf ("retesting soln\n");
			fflush (stdout);
		}

		rval = ILLsimplex_retest_dsolution (lp, pinf, cphase, &fi);
		CHECKRVALG (rval, CLEANUP);
		ILLfct_set_status_values (lp, fi.pstatus, fi.dstatus, PHASEI, PHASEII);

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

	ILLfct_compute_zz (lp, &(lp->zz), pr.lindex);
	//ILL_IFTRACE("%s:c\n",__func__);
	ILLfct_compute_zA (lp, &(lp->zz), &(lp->zA));
	ILLfct_update_counts (lp, CNT_ZNZ, lp->zz.nzcnt, zeroLpNum);
	ILLfct_update_counts (lp, CNT_ZANZ, lp->zA.nzcnt, zeroLpNum);

	ILLratio_dI_test (lp, pr.lindex, pr.lvstat, &rs);

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
		EGlpNumMultUiTo (lp->tol->pfeas_tol, SIMPLEX_FACTOR);
		EGlpNumMultUiTo (lp->tol->dfeas_tol, SIMPLEX_FACTOR);
		//fprintf(stderr,"Resume Numerical %s:%s:%d\n",__func__,__FILE__,__LINE__);
		ILL_CLEANUP;
	}
	else if (rs.ratio_stat == RATIO_BCHANGE)
	{
		ILLfct_compute_yz (lp, &(lp->yjz), updz, lp->nbaz[rs.eindex]);
		rval = ILLfct_test_pivot (lp, pr.lindex, ROW_PIVOT, rs.pivotval);
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
				EGlpNumMultUiTo (lp->tol->pfeas_tol, SIMPLEX_FACTOR);
				EGlpNumMultUiTo (lp->tol->dfeas_tol, SIMPLEX_FACTOR);
				//fprintf(stderr,"Resume Pivot %s:%s:%d\n",__func__,__FILE__,__LINE__);
				rval = 0;
				ILL_CLEANUP;
			}
			rval = ILLbasis_factor (lp, &singular);
			if (singular)
				MESSAGE (__QS_SB_VERB, "Singular basis found!");
			CHECKRVALG (rval, CLEANUP);
			if (singular == 0)
				refactor = 1;
			goto END;
		}
		ILLfct_update_counts (lp, CNT_YNZ, lp->yjz.nzcnt, zeroLpNum);
		ILLfct_update_counts (lp, CNT_UPNZ, updz->nzcnt, zeroLpNum);

		if (pinf->dI_price == QS_PRICE_DSTEEP)
			ILLfct_compute_dsteep_upv (lp, wz);
		rval =
			ILLprice_update_pricing_info (lp, pinf, cphase, wz, rs.eindex, pr.lindex,
																		rs.pivotval);
		CHECKRVALG (rval, CLEANUP);

		EGlpNumSubTo (lp->dinfeas, lp->upd.c_obj);

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
			EGlpNumCopy (it->prevobj, lp->dinfeas);
			it->noprog = 0;
		}

		EGlpNumCopyFrac (alpha, lp->dz[rs.eindex], rs.pivotval);
		EGlpNumCopyFrac (alpha1, lp->xbz[pr.lindex], rs.pivotval);

		ILLfct_update_piz (lp, alpha);
		ILLfct_update_dz (lp, rs.eindex, alpha);
		ILLfct_update_dfeas (lp, rs.eindex, &(lp->srhs));
		ILLfct_compute_dpIy (lp, &(lp->srhs), &(lp->ssoln));
		ILLfct_update_basis_info (lp, rs.eindex, pr.lindex, pr.lvstat);

#if DENSE_PI > 0
		fct_test_workvector (lp);
		fct_test_dfeasible (lp);
#endif
		rval = ILLbasis_update (lp, updz, pr.lindex, &refactor, &singular);
		CHECKRVALG (rval, CLEANUP);

#if DENSE_NORM > 0
		test_dsteep_norms (lp, pinf);
#endif

		ILLfct_update_dpI_prices (lp, pinf, &(lp->srhs), &(lp->ssoln), pr.lindex,
															alpha1);

	END:
		if (singular)
		{
			it->nextstep = SIMPLEX_RESUME;
			it->resumeid = SIMPLEX_RESUME_SING;
			/* this is to force to exit in the case of bad basis */
			it->n_restart++;
			//fprintf(stderr,"Resume Singular %s:%s:%d\n",__func__,__FILE__,__LINE__);
			EGlpNumMultUiTo (lp->tol->pfeas_tol, SIMPLEX_FACTOR);
			EGlpNumMultUiTo (lp->tol->dfeas_tol, SIMPLEX_FACTOR);
			ILL_CLEANUP;
		}
		if (refactor != 0 || it->nosolve > PARAM_MAX_NOSOLVE)
		{
			ILLfct_compute_piz (lp);
			ILLfct_compute_dz (lp);
			ILLfct_dual_adjust (lp, zeroLpNum);
			ILLfct_check_dfeasible (lp, &fi, lp->tol->id_tol);
			ILLfct_set_status_values (lp, -1, fi.dstatus, -1, PHASEII);
			if (fi.dstatus == DUAL_FEASIBLE)
				it->nextphase = DUAL_PHASEII;

			it->newphase = SIMPLEX_PHASE_RECOMP;
			ILL_CLEANUP;
		}
	}

#if DENSE_PI > 1
	fct_test_workvector (lp);
	fct_test_pI_x (lp, pinf);
#endif

CLEANUP:
	EGlpNumClearVar (alpha);
	EGlpNumClearVar (alpha1);
	EGlpNumClearVar (fi.totinfeas);
	EGlpNumClearVar (pr.dinfeas);
	EGlpNumClearVar (pr.pinfeas);
	EGlpNumClearVar (rs.tz);
	EGlpNumClearVar (rs.lbound);
	EGlpNumClearVar (rs.ecoeff);
	EGlpNumClearVar (rs.pivotval);
	//EG_RETURN(rval);
	return rval;
}

static int dual_phaseII_step (
	lpinfo * lp,
	price_info * pinf,
	svector * updz,
	svector * wz,
	iter_info * it)
{
	int coeffch;
	int rval = 0;
	int singular = 0;
	int refactor = 0;
	int ratio_iter = 0;
	int cphase = DUAL_PHASEII;
	int lcol, ecol;
	int estat, newphase;
	EGlpNum_t x_bi, v_l, eval;
	EGlpNum_t ecoeff;
	EGlpNum_t alpha;
	EGlpNum_t alpha1;
	feas_info fi;
	ratio_res rs;
	price_res pr;

	EGlpNumInitVar (x_bi);
	EGlpNumInitVar (v_l);
	EGlpNumInitVar (eval);
	EGlpNumInitVar (ecoeff);
	EGlpNumInitVar (alpha);
	EGlpNumInitVar (alpha1);
	EGlpNumInitVar (fi.totinfeas);
	EGlpNumInitVar (pr.dinfeas);
	EGlpNumInitVar (pr.pinfeas);
	EGlpNumInitVar (rs.tz);
	EGlpNumInitVar (rs.lbound);
	EGlpNumInitVar (rs.ecoeff);
	EGlpNumInitVar (rs.pivotval);
	EGlpNumZero (rs.ecoeff);
	EGlpNumZero (alpha1);

	ILLfct_update_counts (lp, CNT_DPHASE2ITER, 0, zeroLpNum);
	it->nextstep = SIMPLEX_CONTINUE;
	it->nextphase = DUAL_PHASEII;
	lp->final_phase = DUAL_PHASEII;
	newphase = it->newphase;
	it->nosolve++;

	if (it->newphase != 0)
	{
		ILLfct_compute_dobj (lp);
		if (it->newphase == SIMPLEX_PHASE_NEW)
		{
			it->noprog = 0;
			if (it->sdisplay)
			{
				printf ("starting dual phase II, nosolve %d\n", it->nosolve);
				fflush (stdout);
			}
		}
		it->newphase = 0;
		it->nosolve = 0;
		EGlpNumCopy (it->prevobj, lp->dobjval);
		ILLfct_compute_xbz (lp);

		if (pinf->d_strategy == COMPLETE_PRICING)
		{
#if USEHEAP > 0
			ILLprice_free_heap (pinf);
#endif
			ILLprice_compute_primal_inf (lp, pinf, NULL, 0, DUAL_PHASEII);
#if USEHEAP > 0
			rval = ILLprice_test_for_heap (lp, pinf, lp->nrows, pinf->p_scaleinf,
																		 DUAL_SIMPLEX, 0);
			CHECKRVALG (rval, CLEANUP);
#endif
		}
		else if (pinf->d_strategy == MULTI_PART_PRICING)
		{
			ILLprice_init_mpartial_price (lp, pinf, cphase, ROW_PRICING);
		}
	}

	monitor_iter (lp, pinf, it, cphase);
	if (it->nextstep == SIMPLEX_TERMINATE || it->nextstep == SIMPLEX_RESUME ||
			it->newphase != 0)
		ILL_CLEANUP;

	ILLprice_dual (lp, pinf, cphase, &pr);

	if (pr.price_stat == PRICE_OPTIMAL)
	{
		if (lp->ncchange != 0)
		{
			if (it->sdisplay > 1)
			{
				printf ("unrolling %d coef shifts\n", lp->ncchange);
				fflush (stdout);
			}
			ILLfct_unroll_coef_change (lp);
			ILLfct_check_dfeasible (lp, &fi, lp->tol->dfeas_tol);
			ILLfct_set_status_values (lp, -1, fi.dstatus, -1, PHASEII);

			 /*HHH*/ ILLfct_check_pfeasible (lp, &fi, lp->tol->pfeas_tol);
			/*HHH* printf ("dual (opt) infeas %.6f\n", lp->dinfeas); fflush (stdout);
			 *HHH* printf ("primal (opt) infeas %.6f\n", lp->pinfeas); fflush (stdout);*/

			if (fi.dstatus != DUAL_FEASIBLE)
			{
				it->algorithm = PRIMAL_SIMPLEX;
				it->nextstep = SIMPLEX_RESUME;
				it->resumeid = SIMPLEX_RESUME_UNSHIFT;
				it->pricetype = QS_PRICE_PDEVEX;
				/* this is to force to exit in the case of bad basis */
				it->n_restart++;
				EGlpNumMultUiTo (lp->tol->pfeas_tol, SIMPLEX_FACTOR);
				EGlpNumMultUiTo (lp->tol->dfeas_tol, SIMPLEX_FACTOR);
				//fprintf(stderr,"Resume Unshift %s:%s:%d\n",__func__,__FILE__,__LINE__);
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
			printf ("problem seemingly solved\n");
			printf ("seemingly dual opt = %f\n", EGlpNumToLf (lp->dobjval));
			printf ("retesting soln\n");
			fflush (stdout);
		}

		rval = ILLsimplex_retest_dsolution (lp, pinf, cphase, &fi);
		CHECKRVALG (rval, CLEANUP);
		ILLfct_set_status_values (lp, fi.pstatus, fi.dstatus, PHASEII, PHASEII);

		if (fi.dstatus == DUAL_INFEASIBLE)
		{
			ILL_IFTRACE ("DUAL_INFEAS: %s\n", __func__);
			it->nextphase = DUAL_PHASEI;
			EGlpNumDivUiTo (lp->tol->ip_tol, 5);
			EGlpNumDivUiTo (lp->tol->id_tol, 5);
		}
		else if (fi.pstatus == PRIMAL_FEASIBLE)
		{
			ILL_IFTRACE ("PRIM_FEAS: %s\n", __func__);
			EGlpNumCopy (lp->objval, lp->dobjval);
			it->solstatus = ILL_LP_SOLVED;
			it->nextstep = SIMPLEX_TERMINATE;
		}
		else
			ILL_IFTRACE ("PRIM_INFEAS: %s\n", __func__);
		ILL_CLEANUP;
	}

	ILLfct_compute_zz (lp, &(lp->zz), pr.lindex);
	//ILL_IFTRACE("%s:d\n",__func__);
	ILLfct_compute_zA (lp, &(lp->zz), &(lp->zA));
	ILLfct_update_counts (lp, CNT_ZNZ, lp->zz.nzcnt, zeroLpNum);
	ILLfct_update_counts (lp, CNT_ZANZ, lp->zA.nzcnt, zeroLpNum);

	ratio_iter = 0;
	do
	{
		ILLratio_longdII_test (lp, pr.lindex, pr.lvstat, &rs);
		if (rs.ratio_stat == RATIO_NEGATIVE)
		{
			if (it->sdisplay > 1)
			{
				printf ("adjust coefs to remove negative ratio tests\n");
				fflush (stdout);
			}
			ILLfct_adjust_viol_coefs (lp);
			ILLratio_longdII_test (lp, pr.lindex, pr.lvstat, &rs);
			if (rs.ratio_stat == RATIO_NEGATIVE)
			{
				MESSAGE (__QS_SB_VERB, "internal error: bad ratio test");
				fflush (stdout);
				rs.ratio_stat = RATIO_FAILED;
				break;
			}
		}

		coeffch = rs.coeffch;
		EGlpNumCopy (ecoeff, rs.ecoeff);
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
			rval = ILLfct_coef_shift (lp, lp->nbaz[rs.eindex], ecoeff);
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
		EGlpNumMultUiTo (lp->tol->pfeas_tol, SIMPLEX_FACTOR);
		EGlpNumMultUiTo (lp->tol->dfeas_tol, SIMPLEX_FACTOR);
		//fprintf(stderr,"Resume Numerical %s:%s:%d\n",__func__,__FILE__,__LINE__);
		ILL_CLEANUP;
	}
	else if (rs.ratio_stat == RATIO_UNBOUNDED)
	{
		lp->infub_ix = pr.lindex;
		if (lp->ncchange != 0)
		{
			if (it->sdisplay > 1)
			{
				printf ("unrolling %d coef shifts\n", lp->ncchange);
				fflush (stdout);
			}
			ILLfct_unroll_coef_change (lp);
		}
		ILLfct_set_status_values (lp, -1, DUAL_UNBOUNDED, -1, PHASEII);
		it->solstatus = ILL_LP_SOLVED;
		it->nextstep = SIMPLEX_TERMINATE;
	}
	else if (rs.ratio_stat == RATIO_BCHANGE)
	{
		lcol = lp->baz[pr.lindex];
		ecol = lp->nbaz[rs.eindex];

		ILLfct_compute_yz (lp, &(lp->yjz), updz, ecol);
		ILLfct_update_counts (lp, CNT_YNZ, lp->yjz.nzcnt, zeroLpNum);
		ILLfct_update_counts (lp, CNT_UPNZ, updz->nzcnt, zeroLpNum);
		rval = ILLfct_test_pivot (lp, pr.lindex, ROW_PIVOT, rs.pivotval);
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
				EGlpNumMultUiTo (lp->tol->pfeas_tol, SIMPLEX_FACTOR);
				EGlpNumMultUiTo (lp->tol->dfeas_tol, SIMPLEX_FACTOR);
				//fprintf(stderr,"Resume Pivot %s:%s:%d\n",__func__,__FILE__,__LINE__);
				rval = 0;
				ILL_CLEANUP;
			}
			if (newphase == 0)
			{
				rval = ILLbasis_factor (lp, &singular);
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
					printf ("warning: bad step\n");
					fflush (stdout);
				}
			}
		}

		EGlpNumAddTo (lp->dobjval, lp->upd.c_obj);
		EGlpNumCopy (lp->objval, lp->dobjval);

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
			EGlpNumCopy (it->prevobj, lp->dobjval);
			it->noprog = 0;
		}

		if (pinf->dII_price == QS_PRICE_DSTEEP)
			ILLfct_compute_dsteep_upv (lp, wz);
		rval =
			ILLprice_update_pricing_info (lp, pinf, cphase, wz, rs.eindex, pr.lindex,
																		rs.pivotval);
		CHECKRVALG (rval, CLEANUP);

		EGlpNumCopy (x_bi, lp->xbz[pr.lindex]);
		if (pr.lvstat == STAT_LOWER)
			EGlpNumCopy (v_l, lp->lz[lcol]);
		else
			EGlpNumCopy (v_l, lp->uz[lcol]);
		EGlpNumCopy (alpha, rs.tz);
		if (pr.lvstat == STAT_LOWER)
			EGlpNumSign (alpha);
		estat = lp->vstat[ecol];
		if (estat == STAT_LOWER)
			EGlpNumCopy (eval, lp->lz[ecol]);
		else if (estat == STAT_ZERO)
			EGlpNumZero (eval);
		else
			EGlpNumCopy (eval, lp->uz[ecol]);

		ILLfct_update_piz (lp, alpha);
		ILLfct_update_dz (lp, rs.eindex, alpha);
		ILLfct_update_dIIfeas (lp, rs.eindex, &(lp->srhs));
		ILLfct_compute_dpIIy (lp, &(lp->srhs), &(lp->ssoln));
		EGlpNumCopyDiff (alpha1, x_bi, v_l);
		EGlpNumSubTo (alpha1, lp->upd.dty);
		EGlpNumDivTo (alpha1, rs.pivotval);
		ILLfct_update_basis_info (lp, rs.eindex, pr.lindex, pr.lvstat);
		rval = ILLbasis_update (lp, updz, pr.lindex, &refactor, &singular);
		CHECKRVALG (rval, CLEANUP);

		ILLfct_update_dpII_prices (lp, pinf, &(lp->srhs), &(lp->ssoln), /*rs.eindex,*/
															 pr.lindex, eval, alpha1);

#if DENSE_NORM > 0
		test_dsteep_norms (lp, pinf);
#endif

	END:
		if (singular)
		{
			it->nextstep = SIMPLEX_RESUME;
			it->resumeid = SIMPLEX_RESUME_SING;
			/* this is to force to exit in the case of bad basis */
			it->n_restart++;
			EGlpNumMultUiTo (lp->tol->pfeas_tol, SIMPLEX_FACTOR);
			EGlpNumMultUiTo (lp->tol->dfeas_tol, SIMPLEX_FACTOR);
			//fprintf(stderr,"Resume Singular %s:%s:%d\n",__func__,__FILE__,__LINE__);
			ILL_CLEANUP;
		}
		if (refactor != 0 || it->nosolve > PARAM_MAX_NOSOLVE)
		{
			ILLfct_compute_piz (lp);
			ILLfct_compute_dz (lp);
			ILLfct_dual_adjust (lp, zeroLpNum);
			it->newphase = SIMPLEX_PHASE_RECOMP;
		}
	}

#if DENSE_PIIPI > 0
	fct_test_workvector (lp);
	if (!refactor)
	{
		fct_test_pII_x (lp, pinf);
		fct_test_pII_pi_dz (lp, pinf);
	}
#endif

CLEANUP:
	EGlpNumClearVar (x_bi);
	EGlpNumClearVar (v_l);
	EGlpNumClearVar (eval);
	EGlpNumClearVar (ecoeff);
	EGlpNumClearVar (alpha);
	EGlpNumClearVar (alpha1);
	EGlpNumClearVar (fi.totinfeas);
	EGlpNumClearVar (pr.dinfeas);
	EGlpNumClearVar (pr.pinfeas);
	EGlpNumClearVar (rs.tz);
	EGlpNumClearVar (rs.lbound);
	EGlpNumClearVar (rs.ecoeff);
	EGlpNumClearVar (rs.pivotval);
	//EG_RETURN(rval);
	return rval;
}

static void get_current_stat (
	lp_status_info * p,
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

int ILLsimplex_pivotin (
	lpinfo * lp,
	price_info * pinf,
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
	svector wz;
	svector updz;
	EGlpNum_t alpha;
	ratio_res rs;
	feas_info fi;

	EGlpNumInitVar (alpha);
	EGlpNumInitVar (fi.totinfeas);
	EGlpNumInitVar (rs.tz);
	EGlpNumInitVar (rs.lbound);
	EGlpNumInitVar (rs.ecoeff);
	EGlpNumInitVar (rs.pivotval);
	EGlpNumZero (alpha);

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

	/* printf ("Forcing vars into basis in ILLsimplex_pivotin \n"); */
	ILLsvector_init (&wz);
	rval = ILLsvector_alloc (&wz, lp->nrows);
	CHECKRVALG (rval, CLEANUP);
	ILLsvector_init (&updz);
	rval = ILLsvector_alloc (&updz, lp->nrows);
	CHECKRVALG (rval, CLEANUP);

	EGlpNumCopy (lp->pobjval, lp->dobjval);
	for (i = 0; i < rcnt; i++)
	{
		if (lp->vstat[clist[i]] == STAT_BASIC)
			continue;
		npiv++;

		eindex = lp->vindex[clist[i]];
		ILLfct_compute_yz (lp, &(lp->yjz), &updz, lp->nbaz[eindex]);
		ILLfct_update_counts (lp, CNT_YNZ, lp->yjz.nzcnt, zeroLpNum);
		ILLfct_update_counts (lp, CNT_UPNZ, updz.nzcnt, zeroLpNum);

		ILLratio_pivotin_test (lp, clist, rcnt, &rs);

		if (rs.ratio_stat == RATIO_UNBOUNDED || rs.ratio_stat == RATIO_FAILED)
		{
			fprintf (stderr, "Pivot_in failed\n");
			rval = E_SIMPLEX_ERROR;
			ILL_CLEANUP;
		}
		else if (rs.ratio_stat == RATIO_BCHANGE)
		{
			if (rs.lvstat == STAT_LOWER)
			{
				EGlpNumCopyDiff (alpha, lp->lz[lp->baz[rs.lindex]], lp->xbz[rs.lindex]);
				EGlpNumAddInnProdTo (lp->dobjval, rs.tz, alpha);
			}
			else
			{
				EGlpNumCopyDiff (alpha, lp->xbz[rs.lindex], lp->uz[lp->baz[rs.lindex]]);
				EGlpNumAddInnProdTo (lp->dobjval, rs.tz, alpha);
			}
			EGlpNumCopyFrac (alpha, lp->dz[eindex], rs.pivotval);
			EGlpNumCopy (lp->objval, lp->dobjval);

			ILLfct_compute_zz (lp, &(lp->zz), rs.lindex);
			//ILL_IFTRACE("%s:e\n",__func__);
			ILLfct_compute_zA (lp, &(lp->zz), &(lp->zA));
			ILLfct_update_counts (lp, CNT_ZNZ, lp->zz.nzcnt, zeroLpNum);
			ILLfct_update_counts (lp, CNT_ZANZ, lp->zA.nzcnt, zeroLpNum);

			if (pinf->dsinfo.norms && pinf->dII_price == QS_PRICE_DSTEEP)
			{
				ILLfct_compute_dsteep_upv (lp, &wz);
				rval = ILLprice_update_pricing_info (lp, pinf, DUAL_PHASEII, &wz,
																						 eindex, rs.lindex, rs.pivotval);
				CHECKRVALG (rval, CLEANUP);
			}
			else if (pinf->psinfo.norms && pinf->pII_price == QS_PRICE_PSTEEP)
			{
				ILLfct_compute_psteep_upv (lp, &wz);
				rval = ILLprice_update_pricing_info (lp, pinf, PRIMAL_PHASEII, &wz,
																						 eindex, rs.lindex, rs.pivotval);
				CHECKRVALG (rval, CLEANUP);
			}

			//ILL_IFTRACE("%s:e\n",__func__);
			ILLfct_update_xz (lp, rs.tz, eindex, rs.lindex);
			ILLfct_update_basis_info (lp, eindex, rs.lindex, rs.lvstat);
			rval = ILLbasis_update (lp, &updz, rs.lindex, &refactor, &singular);
			CHECKRVALG (rval, CLEANUP);

			if (singular)
			{
				fprintf (stderr, "singular matrix in pivot_in\n");
				rval = E_SIMPLEX_ERROR;
				ILL_CLEANUP;
			}
			if (!refactor)
			{
				ILLfct_update_piz (lp, alpha);
				ILLfct_update_dz (lp, eindex, alpha);
			}
			else
			{
				ILLfct_compute_xbz (lp);
				ILLfct_compute_piz (lp);
				ILLfct_compute_dz (lp);
				ILLfct_compute_dobj (lp);
			}
		}
	}
	/*
	 * ILLfct_dphaseI_simple_update (lp, lp->tol->dfeas_tol);
	 * ILLfct_compute_xbz (lp);
	 * ILLfct_compute_dobj (lp);
	 */

	ILLfct_check_pfeasible (lp, &fi, lp->tol->pfeas_tol);
	ILLfct_check_dfeasible (lp, &fi, lp->tol->dfeas_tol);
	ILLfct_set_status_values (lp, fi.pstatus, fi.dstatus, PHASEII, PHASEII);

CLEANUP:
	if (pivot_opt == SIMPLEX_PIVOTINROW)
		ILL_IFFREE (clist, int);

	ILLsvector_free (&wz);
	ILLsvector_free (&updz);
	EGlpNumClearVar (alpha);
	EGlpNumClearVar (fi.totinfeas);
	EGlpNumClearVar (rs.tz);
	EGlpNumClearVar (rs.lbound);
	EGlpNumClearVar (rs.ecoeff);
	EGlpNumClearVar (rs.pivotval);
	EG_RETURN (rval);
}

static int report_value (
	lpinfo * lp,
	iter_info * it,
	const char *value_name,
	EGlpNum_t value)
{
	int rval = 0;

	if (it->sdisplay && it->itercnt % lp->iterskip == 0)
	{
		char buffer[1024];

		snprintf (buffer, (size_t) 1023, "(%d): %s = %10.7lf\n", it->itercnt,
							value_name, EGlpNumToLf (value));
		buffer[1022] = '\n';
		buffer[1023] = '\0';
		rval = ILLstring_report (buffer, &lp->O->reporter);
		fflush (stdout);
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
