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

/* RCS_INFO = "$RCSfile: price.c,v $ $Revision: 1.2 $ $Date: 2003/11/05 16:49:52 $"; */
//static int TRACE = 0;

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "qs_config.h"
#include "logging-private.h"

#include "allocrus.h"
#include "eg_lpnum.h"
#include "eg_numutil_EGLPNUM_TYPENAME.h"
#include "eg_io.h"
#include "except.h"
#include "util.h"

#include "stddefs.h"
#include "qsopt_EGLPNUM_TYPENAME.h"
#include "lpdefs_EGLPNUM_TYPENAME.h"
#include "fct_EGLPNUM_TYPENAME.h"
#include "price_EGLPNUM_TYPENAME.h"
#include "basis_EGLPNUM_TYPENAME.h"
#include "dstruct_EGLPNUM_TYPENAME.h"

#define  MULTIP 1
#define  PRICE_DEBUG 0


static void update_d_scaleinf (
	EGLPNUM_TYPENAME_price_info * const p,
	EGLPNUM_TYPENAME_heap * const h,
	int const j,
	EGLPNUM_TYPE inf,
	int const prule),
  update_p_scaleinf (
	EGLPNUM_TYPENAME_price_info * const p,
	EGLPNUM_TYPENAME_heap * const h,
	int const i,
	EGLPNUM_TYPE inf,
	int const prule);

static void compute_dualI_inf (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	int const j,
	EGLPNUM_TYPE * const inf),
  compute_dualII_inf (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	int const j,
	EGLPNUM_TYPE * const inf),
  compute_primalI_inf (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	int const i,
	EGLPNUM_TYPE * const inf),
  compute_primalII_inf (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	int const i,
	EGLPNUM_TYPE * const inf);

void EGLPNUM_TYPENAME_ILLprice_free_heap (
	EGLPNUM_TYPENAME_price_info * const pinf)
{
	EGLPNUM_TYPENAME_ILLheap_free (&(pinf->h));
}

int EGLPNUM_TYPENAME_ILLprice_build_heap (
	EGLPNUM_TYPENAME_price_info * const pinf,
	int const nkeys,
	EGLPNUM_TYPE * keylist)
{
	EGLPNUM_TYPENAME_ILLheap_init (&(pinf->h));
	EGLPNUM_TYPENAME_EGlpNumSet (pinf->htrigger,
							1.0 +
							(double) nkeys / (PARAM_HEAP_RATIO * ILLutil_our_log2 (nkeys)));
	return EGLPNUM_TYPENAME_ILLheap_build (&(pinf->h), nkeys, keylist);
}

int EGLPNUM_TYPENAME_ILLprice_test_for_heap (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	EGLPNUM_TYPENAME_price_info * const pinf,
	int const nkeys,
	EGLPNUM_TYPE * keylist,
	int const algo,
	int const upd)
{
	EGLPNUM_TYPENAME_heap *const h = &(pinf->h);
	int rval = 0;
	EGLPNUM_TYPE ravg;

	if (upd != 0)
	{
		EGLPNUM_TYPENAME_EGlpNumInitVar (ravg);
		if (algo == PRIMAL_SIMPLEX)
			EGLPNUM_TYPENAME_EGlpNumCopy (ravg, lp->cnts->za_ravg);
		else
			EGLPNUM_TYPENAME_EGlpNumCopy (ravg, lp->cnts->y_ravg);
		if (EGLPNUM_TYPENAME_EGlpNumIsLeq (ravg, pinf->htrigger))
			pinf->hineff--;
		else
		{
			EGLPNUM_TYPENAME_EGlpNumDivUiTo (ravg, 2U);
			if (EGLPNUM_TYPENAME_EGlpNumIsLess (pinf->htrigger, ravg))
				pinf->hineff++;
		}
		EGLPNUM_TYPENAME_EGlpNumClearVar (ravg);
	}
	if (h->hexist == 0 && pinf->hineff <= 0)
	{
		rval = EGLPNUM_TYPENAME_ILLprice_build_heap (pinf, nkeys, keylist);
		CHECKRVALG (rval, CLEANUP);
	}
	else if (h->hexist != 0 && pinf->hineff >= PARAM_HEAP_UTRIGGER)
	{
		EGLPNUM_TYPENAME_ILLprice_free_heap (pinf);
		/*
		 * QSlog("freeing EGLPNUM_TYPENAME_heap ..");
		 * QSlog("iter = %d, ravg = %.2f, trigger = %.2f",
		 * lp->cnts->tot_iter, ravg, pinf->htrigger);
		 */
	}

CLEANUP:
	if (rval)
		EGLPNUM_TYPENAME_ILLprice_free_heap (pinf);
	return rval;
}

void EGLPNUM_TYPENAME_ILLprice_init_pricing_info (
	EGLPNUM_TYPENAME_price_info * const pinf)
{
	pinf->p_strategy = -1;
	pinf->d_strategy = -1;
	pinf->pI_price = -1;
	pinf->pII_price = -1;
	pinf->dI_price = -1;
	pinf->dII_price = -1;
	pinf->cur_price = -1;
	pinf->p_scaleinf = 0;
	pinf->d_scaleinf = 0;
	pinf->pdinfo.norms = 0;
	pinf->pdinfo.refframe = 0;
	pinf->psinfo.norms = 0;
	pinf->ddinfo.norms = 0;
	pinf->ddinfo.refframe = 0;
	pinf->dsinfo.norms = 0;
	pinf->dmpinfo.gstart = pinf->pmpinfo.gstart = 0;
	pinf->dmpinfo.gshift = pinf->pmpinfo.gshift = 0;
	pinf->dmpinfo.gsize = pinf->pmpinfo.gsize = 0;
	pinf->dmpinfo.bucket = pinf->pmpinfo.bucket = 0;
	pinf->dmpinfo.perm = pinf->pmpinfo.perm = 0;
	pinf->dmpinfo.infeas = pinf->pmpinfo.infeas = 0;
	EGLPNUM_TYPENAME_ILLheap_init (&(pinf->h));
	EGLPNUM_TYPENAME_EGlpNumZero (pinf->htrigger);
	pinf->hineff = 0;
}

void EGLPNUM_TYPENAME_ILLprice_free_pricing_info (
	EGLPNUM_TYPENAME_price_info * const pinf)
{
	EGLPNUM_TYPENAME_EGlpNumFreeArray (pinf->p_scaleinf);
	EGLPNUM_TYPENAME_EGlpNumFreeArray (pinf->d_scaleinf);
	EGLPNUM_TYPENAME_EGlpNumFreeArray (pinf->pdinfo.norms);
	ILL_IFFREE (pinf->pdinfo.refframe, int);
	EGLPNUM_TYPENAME_EGlpNumFreeArray (pinf->psinfo.norms);
	EGLPNUM_TYPENAME_EGlpNumFreeArray (pinf->ddinfo.norms);
	ILL_IFFREE (pinf->ddinfo.refframe, int);
	EGLPNUM_TYPENAME_EGlpNumFreeArray (pinf->dsinfo.norms);

	EGLPNUM_TYPENAME_ILLprice_free_mpartial_info (&(pinf->pmpinfo));
	EGLPNUM_TYPENAME_ILLprice_free_mpartial_info (&(pinf->dmpinfo));
	EGLPNUM_TYPENAME_ILLprice_free_heap (pinf);
}

int EGLPNUM_TYPENAME_ILLprice_build_pricing_info (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	EGLPNUM_TYPENAME_price_info * const pinf,
	int const phase)
{
	int rval = 0;
	int p_price = -1;
	int d_price = -1;

	switch (phase)
	{
	case PRIMAL_PHASEI:
		p_price = pinf->pI_price;
		break;
	case PRIMAL_PHASEII:
		p_price = pinf->pII_price;
		break;
	case DUAL_PHASEI:
		d_price = pinf->dI_price;
		break;
	case DUAL_PHASEII:
		d_price = pinf->dII_price;
		break;
	}

	if (p_price != -1)
	{
		pinf->cur_price = p_price;

		if (p_price == QS_PRICE_PDANTZIG || p_price == QS_PRICE_PDEVEX ||
				p_price == QS_PRICE_PSTEEP)
		{
			pinf->p_strategy = COMPLETE_PRICING;
			EGLPNUM_TYPENAME_EGlpNumFreeArray (pinf->d_scaleinf);
			pinf->d_scaleinf = EGLPNUM_TYPENAME_EGlpNumAllocArray (lp->nnbasic);
		}
		else if (p_price == QS_PRICE_PMULTPARTIAL)
			pinf->p_strategy = MULTI_PART_PRICING;

		switch (p_price)
		{
		case QS_PRICE_PDEVEX:
			if (pinf->pdinfo.norms)
				return rval;
			rval = EGLPNUM_TYPENAME_ILLprice_build_pdevex_norms (lp, &(pinf->pdinfo), 0);
			CHECKRVALG(rval,CLEANUP);
			break;
		case QS_PRICE_PSTEEP:
			if (pinf->psinfo.norms)
				return rval;
			rval = EGLPNUM_TYPENAME_ILLprice_build_psteep_norms (lp, &(pinf->psinfo));
			CHECKRVALG(rval,CLEANUP);
			break;
		case QS_PRICE_PMULTPARTIAL:
			rval = EGLPNUM_TYPENAME_ILLprice_build_mpartial_info (lp, pinf, COL_PRICING);
			CHECKRVALG(rval,CLEANUP);
			break;
		}
	}
	else if (d_price != -1)
	{
		pinf->cur_price = d_price;

		if (d_price == QS_PRICE_DDANTZIG || d_price == QS_PRICE_DSTEEP ||
				d_price == QS_PRICE_DDEVEX)
		{
			pinf->d_strategy = COMPLETE_PRICING;
			EGLPNUM_TYPENAME_EGlpNumFreeArray (pinf->p_scaleinf);
			pinf->p_scaleinf = EGLPNUM_TYPENAME_EGlpNumAllocArray (lp->nrows);
		}
		else if (d_price == QS_PRICE_DMULTPARTIAL)
			pinf->d_strategy = MULTI_PART_PRICING;

		switch (d_price)
		{
		case QS_PRICE_DSTEEP:
			if (pinf->dsinfo.norms)
				return rval;
			rval = EGLPNUM_TYPENAME_ILLprice_build_dsteep_norms (lp, &(pinf->dsinfo));
			CHECKRVALG(rval,CLEANUP);
			break;
		case QS_PRICE_DMULTPARTIAL:
			rval = EGLPNUM_TYPENAME_ILLprice_build_mpartial_info (lp, pinf, ROW_PRICING);
			CHECKRVALG(rval,CLEANUP);
			break;
		case QS_PRICE_DDEVEX:
			if (pinf->ddinfo.norms)
				return rval;
			rval = EGLPNUM_TYPENAME_ILLprice_build_ddevex_norms (lp, &(pinf->ddinfo), 0);
			CHECKRVALG(rval,CLEANUP);
			break;
		}
	}

CLEANUP:
	if (rval)
		EGLPNUM_TYPENAME_ILLprice_free_pricing_info (pinf);
	EG_RETURN(rval);
}

int EGLPNUM_TYPENAME_ILLprice_update_pricing_info (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	EGLPNUM_TYPENAME_price_info * const pinf,
	int const phase,
	EGLPNUM_TYPENAME_svector * const wz,
	int const eindex,
	int const lindex,
	EGLPNUM_TYPE y)
{
	int rval = 0;
	int p_price = -1;
	int d_price = -1;

	switch (phase)
	{
	case PRIMAL_PHASEI:
		p_price = pinf->pI_price;
		break;
	case PRIMAL_PHASEII:
		p_price = pinf->pII_price;
		break;
	case DUAL_PHASEI:
		d_price = pinf->dI_price;
		break;
	case DUAL_PHASEII:
		d_price = pinf->dII_price;
		break;
	}

	if (p_price != -1)
	{
		if (p_price == QS_PRICE_PDEVEX)
		{
			rval = EGLPNUM_TYPENAME_ILLprice_update_pdevex_norms (lp, &(pinf->pdinfo), eindex, y);
			CHECKRVALG(rval,CLEANUP);
		}
		else if (p_price == QS_PRICE_PSTEEP)
			EGLPNUM_TYPENAME_ILLprice_update_psteep_norms (lp, &(pinf->psinfo), wz, eindex, y);
	}
	else if (d_price != -1)
	{
		if (d_price == QS_PRICE_DSTEEP)
			EGLPNUM_TYPENAME_ILLprice_update_dsteep_norms (lp, &(pinf->dsinfo), wz, lindex, y);
		else if (d_price == QS_PRICE_DDEVEX)
		{
			rval = EGLPNUM_TYPENAME_ILLprice_update_ddevex_norms (lp, &(pinf->ddinfo), lindex, y);
			CHECKRVALG(rval,CLEANUP);
		}
	}
CLEANUP:
	EG_RETURN(rval);
}

int EGLPNUM_TYPENAME_ILLprice_get_price (
	EGLPNUM_TYPENAME_price_info * const p,
	int const phase)
{
	int pri = -1;

	switch (phase)
	{
	case PRIMAL_PHASEI:
		return p->pI_price;
	case PRIMAL_PHASEII:
		return p->pII_price;
	case DUAL_PHASEI:
		return p->dI_price;
	case DUAL_PHASEII:
		return p->dII_price;
	}
	return pri;
}

void EGLPNUM_TYPENAME_ILLprice_free_mpartial_info (
	EGLPNUM_TYPENAME_mpart_info * p)
{
	ILL_IFFREE (p->gstart, int);
	ILL_IFFREE (p->gshift, int);
	ILL_IFFREE (p->gsize, int);
	ILL_IFFREE (p->bucket, int);
	EGLPNUM_TYPENAME_EGlpNumFreeArray (p->infeas);
	ILL_IFFREE (p->perm, int);
}

int EGLPNUM_TYPENAME_ILLprice_build_mpartial_info (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	EGLPNUM_TYPENAME_price_info * const pinf,
	int const pricetype)
{
	int i = 0;
	int rval = 0;
	int extra = 0;
	int nelems;
	EGLPNUM_TYPENAME_mpart_info *p;

	p = (pricetype == COL_PRICING) ? &(pinf->pmpinfo) : &(pinf->dmpinfo);
	p->k = 50;
	p->cgroup = 0;
	nelems = (pricetype == COL_PRICING) ? lp->nnbasic : lp->nrows;

	if (nelems % p->k)
		extra = nelems - p->k * (nelems / p->k);
	p->ngroups = nelems / p->k;
	if (extra != 0)
		p->ngroups++;

	ILL_SAFE_MALLOC (p->gstart, p->ngroups, int);
	ILL_SAFE_MALLOC (p->gshift, p->ngroups, int);
	ILL_SAFE_MALLOC (p->gsize, p->ngroups, int);
	ILL_SAFE_MALLOC (p->bucket, 2 * p->k, int);
	p->infeas = EGLPNUM_TYPENAME_EGlpNumAllocArray (2 * p->k);
	ILL_SAFE_MALLOC (p->perm, 2 * p->k, int);

	p->bsize = 0;

	if (extra != 0)
	{
		p->gstart[0] = 0;
		p->gshift[0] = 1;
		p->gsize[0] = extra;
		for (i = 1; i < p->ngroups; i++)
		{
			p->gstart[i] = extra + i - 1;
			p->gshift[i] = p->ngroups - 1;
			p->gsize[i] = p->k;
		}
	}
	else
	{
		for (i = 0; i < p->ngroups; i++)
		{
			p->gstart[i] = i;
			p->gshift[i] = p->ngroups;
			p->gsize[i] = p->k;
		}
	}

CLEANUP:
	if (rval)
		EGLPNUM_TYPENAME_ILLprice_free_mpartial_info (p);
	EG_RETURN(rval);
}

void EGLPNUM_TYPENAME_ILLprice_init_mpartial_price (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	EGLPNUM_TYPENAME_price_info * const pinf,
	int const phase,
	int const pricetype)
{
	int i;
	EGLPNUM_TYPENAME_mpart_info *p;

	p = (pricetype == COL_PRICING) ? &(pinf->pmpinfo) : &(pinf->dmpinfo);
	p->bsize = 0;
	i = p->cgroup;
	do
	{
		EGLPNUM_TYPENAME_ILLprice_mpartial_group (lp, p, phase, i, pricetype);
		i = (i + 1) % p->ngroups;
	} while (i != p->cgroup && p->bsize <= p->k);
	p->cgroup = i;
}

void EGLPNUM_TYPENAME_ILLprice_update_mpartial_price (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	EGLPNUM_TYPENAME_price_info * const pinf,
	int const phase,
	int const pricetype)
{
	int i = 0;
	int csize = 0;
	EGLPNUM_TYPE infeas;
	EGLPNUM_TYPENAME_mpart_info *p;
	EGLPNUM_TYPENAME_price_res pr;

	p = (pricetype == COL_PRICING) ? &(pinf->pmpinfo) : &(pinf->dmpinfo);
	EGLPNUM_TYPENAME_EGlpNumInitVar (pr.dinfeas);
	EGLPNUM_TYPENAME_EGlpNumInitVar (pr.pinfeas);
	EGLPNUM_TYPENAME_EGlpNumInitVar (infeas);

#ifdef MULTIP
	i = 0;
	while (i < p->bsize)
	{
		if (pricetype == COL_PRICING)
		{
			EGLPNUM_TYPENAME_ILLprice_column (lp, p->bucket[i], phase, &pr);
			EGLPNUM_TYPENAME_EGlpNumCopy (infeas, pr.dinfeas);
		}
		else
		{
			EGLPNUM_TYPENAME_ILLprice_row (lp, p->bucket[i], phase, &pr);
			EGLPNUM_TYPENAME_EGlpNumCopy (infeas, pr.pinfeas);
		}
		if (!EGLPNUM_TYPENAME_EGlpNumIsNeqqZero (infeas))
		{
			p->bucket[i] = p->bucket[p->bsize - 1];
			p->bsize--;
		}
		else
		{
			EGLPNUM_TYPENAME_EGlpNumCopy (p->infeas[i], infeas);
			i++;
		}
	}
	if (p->bsize > 0)
	{
		for (i = 0; i < p->bsize; i++)
			p->perm[i] = i;
		EGLPNUM_TYPENAME_EGutilPermSort ((size_t) (p->bsize), p->perm,
										(const EGLPNUM_TYPE * const) p->infeas);

		csize = QSMIN (p->bsize, p->k);
		for (i = csize - 1; i >= 0; i--)
			lp->iwork[p->bucket[p->perm[i]]] = 1;

		for (i = 0, csize = 0; i < p->bsize; i++)
			if (lp->iwork[p->bucket[i]] == 1)
			{
				EGLPNUM_TYPENAME_EGlpNumCopy (p->infeas[csize], p->infeas[i]);
				p->bucket[csize] = p->bucket[i];
				csize++;
			}
		p->bsize = csize;
	}
#else
	p->bsize = 0;
#endif

	i = p->cgroup;
	do
	{
		EGLPNUM_TYPENAME_ILLprice_mpartial_group (lp, p, phase, i, pricetype);
		i = (i + 1) % p->ngroups;
	} while (i != p->cgroup && p->bsize <= p->k);
	p->cgroup = i;

#ifdef MULTIP
	for (i = 0; i < csize; i++)
		lp->iwork[p->bucket[i]] = 0;
#endif
	EGLPNUM_TYPENAME_EGlpNumClearVar (infeas);
	EGLPNUM_TYPENAME_EGlpNumClearVar (pr.pinfeas);
	EGLPNUM_TYPENAME_EGlpNumClearVar (pr.dinfeas);
}

void EGLPNUM_TYPENAME_ILLprice_delete_onempart_price (
	/*EGLPNUM_TYPENAME_lpinfo * const lp,*/
	EGLPNUM_TYPENAME_price_info * const pinf,
	int const indx,
	int const pricetype)
{
	int i = 0;
	EGLPNUM_TYPENAME_mpart_info *p;

	p = (pricetype == COL_PRICING) ? &(pinf->pmpinfo) : &(pinf->dmpinfo);

	for (i = 0; i < p->bsize; i++)
		if (p->bucket[i] == indx)
		{
			p->bucket[i] = p->bucket[p->bsize - 1];
			EGLPNUM_TYPENAME_EGlpNumCopy (p->infeas[i], p->infeas[p->bsize - 1]);
			p->bsize--;
			break;
		}
}

void EGLPNUM_TYPENAME_ILLprice_mpartial_group (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	EGLPNUM_TYPENAME_mpart_info * const p,
	int const phase,
	int const g,
	int const pricetype)
{
	int i, ix;
	int gstart = p->gstart[g];
	int gsize = p->gsize[g];
	int gshift = p->gshift[g];
	EGLPNUM_TYPE infeas;
	EGLPNUM_TYPENAME_price_res pr;

	EGLPNUM_TYPENAME_EGlpNumInitVar (pr.dinfeas);
	EGLPNUM_TYPENAME_EGlpNumInitVar (pr.pinfeas);
	EGLPNUM_TYPENAME_EGlpNumInitVar (infeas);

	for (i = 0, ix = gstart; i < gsize; i++, ix += gshift)
	{
#ifdef MULTIP
		if (lp->iwork[ix])
			continue;
#endif
		if (pricetype == COL_PRICING)
		{
			EGLPNUM_TYPENAME_ILLprice_column (lp, ix, phase, &pr);
			EGLPNUM_TYPENAME_EGlpNumCopy (infeas, pr.dinfeas);
		}
		else
		{
			EGLPNUM_TYPENAME_ILLprice_row (lp, ix, phase, &pr);
			EGLPNUM_TYPENAME_EGlpNumCopy (infeas, pr.pinfeas);
		}
		if (EGLPNUM_TYPENAME_EGlpNumIsNeqqZero (infeas))
		{
			EGLPNUM_TYPENAME_EGlpNumCopy (p->infeas[p->bsize], infeas);
			p->bucket[p->bsize] = ix;
			p->bsize++;
		}
	}
	EGLPNUM_TYPENAME_EGlpNumClearVar (infeas);
	EGLPNUM_TYPENAME_EGlpNumClearVar (pr.dinfeas);
	EGLPNUM_TYPENAME_EGlpNumClearVar (pr.pinfeas);
}

void EGLPNUM_TYPENAME_ILLprice_column (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	int const ix,
	int const phase,
	EGLPNUM_TYPENAME_price_res * const pr)
{
	int i;
	int col;
	int mcnt;
	int mbeg;
	EGLPNUM_TYPE sum;

	EGLPNUM_TYPENAME_EGlpNumZero (pr->dinfeas);
	col = lp->nbaz[ix];
	if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED)
		return;
	EGLPNUM_TYPENAME_EGlpNumInitVar (sum);
	EGLPNUM_TYPENAME_EGlpNumZero (sum);
	mcnt = lp->matcnt[col];
	mbeg = lp->matbeg[col];

	if (phase == PRIMAL_PHASEII)
	{
		for (i = 0; i < mcnt; i++)
			EGLPNUM_TYPENAME_EGlpNumAddInnProdTo (sum, lp->piz[lp->matind[mbeg + i]], lp->matval[mbeg + i]);
		EGLPNUM_TYPENAME_EGlpNumCopyDiff (lp->dz[ix], lp->cz[col], sum);
		compute_dualII_inf (lp, ix, &(pr->dinfeas));
	}
	else
	{
		for (i = 0; i < mcnt; i++)
			EGLPNUM_TYPENAME_EGlpNumAddInnProdTo (sum, lp->pIpiz[lp->matind[mbeg + i]], lp->matval[mbeg + i]);
		EGLPNUM_TYPENAME_EGlpNumCopyNeg (lp->pIdz[ix], sum);
		compute_dualI_inf (lp, ix, &(pr->dinfeas));
	}
	EGLPNUM_TYPENAME_EGlpNumClearVar (sum);
}

void EGLPNUM_TYPENAME_ILLprice_row (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	int const ix,
	int const phase,
	EGLPNUM_TYPENAME_price_res * const pr)
{
	if (phase == DUAL_PHASEII)
		compute_primalII_inf (lp, ix, &(pr->pinfeas));
	else
		compute_primalI_inf (lp, ix, &(pr->pinfeas));
}

int EGLPNUM_TYPENAME_ILLprice_build_pdevex_norms (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	EGLPNUM_TYPENAME_p_devex_info * const pdinfo,
	int const reinit)
{
	int j;
	int rval = 0;

	if (reinit == 0)
	{
		pdinfo->ninit = 0;
		pdinfo->norms = EGLPNUM_TYPENAME_EGlpNumAllocArray (lp->nnbasic);
		ILL_SAFE_MALLOC (pdinfo->refframe, lp->ncols, int);
	}

	if (reinit != 0)
		pdinfo->ninit++;

	for (j = 0; j < lp->ncols; j++)
	{
		if (lp->vstat[j] == STAT_BASIC)
			pdinfo->refframe[j] = 0;
		else
		{
			EGLPNUM_TYPENAME_EGlpNumOne (pdinfo->norms[lp->vindex[j]]);
			pdinfo->refframe[j] = 1;
		}
	}

CLEANUP:
	if (rval)
	{
		EGLPNUM_TYPENAME_EGlpNumFreeArray (pdinfo->norms);
		ILL_IFFREE (pdinfo->refframe, int);
	}
	EG_RETURN(rval);
}

int EGLPNUM_TYPENAME_ILLprice_update_pdevex_norms (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	EGLPNUM_TYPENAME_p_devex_info * const pdinfo,
	int const eindex,
	EGLPNUM_TYPE yl)
{
	int i, j;
	EGLPNUM_TYPE normj;
	EGLPNUM_TYPE zAj;
	EGLPNUM_TYPE ntmp, ntmp2;

	EGLPNUM_TYPENAME_EGlpNumInitVar (normj);
	EGLPNUM_TYPENAME_EGlpNumInitVar (zAj);
	EGLPNUM_TYPENAME_EGlpNumInitVar (ntmp);
	EGLPNUM_TYPENAME_EGlpNumInitVar (ntmp2);
	EGLPNUM_TYPENAME_EGlpNumZero (normj);

	for (i = 0; i < lp->yjz.nzcnt; i++)
		if (pdinfo->refframe[lp->baz[lp->yjz.indx[i]]])
			EGLPNUM_TYPENAME_EGlpNumAddInnProdTo (normj, lp->yjz.coef[i], lp->yjz.coef[i]);

	if (pdinfo->refframe[lp->nbaz[eindex]])
		EGLPNUM_TYPENAME_EGlpNumAddTo (normj, EGLPNUM_TYPENAME_oneLpNum);

	EGLPNUM_TYPENAME_EGlpNumSet(ntmp,1000.0);
	EGLPNUM_TYPENAME_EGlpNumSet(ntmp2,0.001);
	EGLPNUM_TYPENAME_EGlpNumMultTo(ntmp,pdinfo->norms[eindex]);
	EGLPNUM_TYPENAME_EGlpNumMultTo(ntmp2,pdinfo->norms[eindex]);
	if (EGLPNUM_TYPENAME_EGlpNumIsLess (normj, ntmp2) || EGLPNUM_TYPENAME_EGlpNumIsLess (ntmp, normj))
	{
		EGLPNUM_TYPENAME_EGlpNumClearVar (zAj);
		EGLPNUM_TYPENAME_EGlpNumClearVar (normj);
		EGLPNUM_TYPENAME_EGlpNumClearVar (ntmp);
		EGLPNUM_TYPENAME_EGlpNumClearVar (ntmp2);
		return EGLPNUM_TYPENAME_ILLprice_build_pdevex_norms (lp, pdinfo, 1);
	}

	for (i = 0; i < lp->zA.nzcnt; i++)
	{
		j = lp->zA.indx[i];
		EGLPNUM_TYPENAME_EGlpNumCopyFrac (zAj, lp->zA.coef[i], yl);
		EGLPNUM_TYPENAME_EGlpNumMultTo (zAj, zAj);
		EGLPNUM_TYPENAME_EGlpNumMultTo (zAj, normj);
		if (EGLPNUM_TYPENAME_EGlpNumIsLess (pdinfo->norms[j], zAj))
			EGLPNUM_TYPENAME_EGlpNumCopy (pdinfo->norms[j], zAj);
	}
	EGLPNUM_TYPENAME_EGlpNumDivTo (normj, yl);
	EGLPNUM_TYPENAME_EGlpNumDivTo (normj, yl);
	if (EGLPNUM_TYPENAME_EGlpNumIsLess (normj, EGLPNUM_TYPENAME_oneLpNum))
		EGLPNUM_TYPENAME_EGlpNumCopy (pdinfo->norms[eindex], EGLPNUM_TYPENAME_oneLpNum);
	else
		EGLPNUM_TYPENAME_EGlpNumCopy (pdinfo->norms[eindex], normj);
	EGLPNUM_TYPENAME_EGlpNumClearVar (zAj);
	EGLPNUM_TYPENAME_EGlpNumClearVar (normj);
	EGLPNUM_TYPENAME_EGlpNumClearVar (ntmp);
	EGLPNUM_TYPENAME_EGlpNumClearVar (ntmp2);
	return 0;
}

int EGLPNUM_TYPENAME_ILLprice_build_psteep_norms (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	EGLPNUM_TYPENAME_p_steep_info * const psinfo)
{
	int j;
	int rval = 0;
	EGLPNUM_TYPENAME_svector yz;

	EGLPNUM_TYPENAME_ILLsvector_init (&yz);
	rval = EGLPNUM_TYPENAME_ILLsvector_alloc (&yz, lp->nrows);
	CHECKRVALG(rval,CLEANUP);
	psinfo->norms = EGLPNUM_TYPENAME_EGlpNumAllocArray (lp->nnbasic);

	for (j = 0; j < lp->nnbasic; j++)
	{
		rval = ILLstring_report (NULL, &lp->O->reporter);
		CHECKRVALG(rval,CLEANUP);
		EGLPNUM_TYPENAME_ILLfct_compute_yz (lp, &yz, 0, lp->nbaz[j]);
		EGLPNUM_TYPENAME_EGlpNumInnProd (psinfo->norms[j], yz.coef, yz.coef, (size_t) yz.nzcnt);
		EGLPNUM_TYPENAME_EGlpNumAddTo (psinfo->norms[j], EGLPNUM_TYPENAME_oneLpNum);
	}

CLEANUP:
	EGLPNUM_TYPENAME_ILLsvector_free (&yz);
	if (rval)
		EGLPNUM_TYPENAME_EGlpNumFreeArray (psinfo->norms);

	EG_RETURN(rval);
}

void EGLPNUM_TYPENAME_ILLprice_update_psteep_norms (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	EGLPNUM_TYPENAME_p_steep_info * const psinfo,
	EGLPNUM_TYPENAME_svector * const wz,
	int const eindex,
	EGLPNUM_TYPE yl)
{
	int i, j, k;
	int mcnt, mbeg;
	EGLPNUM_TYPE normj,ntmp;
	EGLPNUM_TYPE zAj, wAj;
	EGLPNUM_TYPE *v = 0;

	EGLPNUM_TYPENAME_EGlpNumInitVar (normj);
	EGLPNUM_TYPENAME_EGlpNumInitVar (zAj);
	EGLPNUM_TYPENAME_EGlpNumInitVar (ntmp);
	EGLPNUM_TYPENAME_EGlpNumInitVar (wAj);
	EGLPNUM_TYPENAME_EGlpNumInnProd (normj, lp->yjz.coef, lp->yjz.coef, (size_t) (lp->yjz.nzcnt));
	EGLPNUM_TYPENAME_EGlpNumAddTo (normj, EGLPNUM_TYPENAME_oneLpNum);

#if 0
	Bico - remove warnings for dist
		if (fabs ((normj - psinfo->norms[eindex]) / normj) > 1000.0 /* 0.01 */ )
		{
			QSlog("warning: incorrect norm values");
			QSlog("anorm = %.6f, pnorm = %.6f", normj, psinfo->norms[eindex]);
		}
#endif

	EGLPNUM_TYPENAME_ILLfct_load_workvector (lp, wz);
	v = lp->work.coef;

	for (k = 0; k < lp->zA.nzcnt; k++)
	{
		j = lp->zA.indx[k];
		EGLPNUM_TYPENAME_EGlpNumCopy (zAj, lp->zA.coef[k]);
		EGLPNUM_TYPENAME_EGlpNumZero (wAj);
		mcnt = lp->matcnt[lp->nbaz[j]];
		mbeg = lp->matbeg[lp->nbaz[j]];
		for (i = 0; i < mcnt; i++)
			EGLPNUM_TYPENAME_EGlpNumAddInnProdTo (wAj, lp->matval[mbeg + i], v[lp->matind[mbeg + i]]);

		/* compute ntmp = (zAj * ((zAj * normj / yl) - (2.0 * wAj))) / yl; */ 
		EGLPNUM_TYPENAME_EGlpNumCopy(ntmp,zAj);
		EGLPNUM_TYPENAME_EGlpNumMultTo(ntmp,normj);
		EGLPNUM_TYPENAME_EGlpNumDivTo(ntmp,yl);
		EGLPNUM_TYPENAME_EGlpNumSubTo(ntmp,wAj);
		EGLPNUM_TYPENAME_EGlpNumSubTo(ntmp,wAj);
		EGLPNUM_TYPENAME_EGlpNumMultTo(ntmp,zAj);
		EGLPNUM_TYPENAME_EGlpNumDivTo(ntmp,yl);
		/* set psinfo->norms[j] += (zAj * ((zAj * normj / yl) - (2.0 * wAj))) / yl; */
		EGLPNUM_TYPENAME_EGlpNumAddTo(psinfo->norms[j],ntmp);
		if (EGLPNUM_TYPENAME_EGlpNumIsLess (psinfo->norms[j], EGLPNUM_TYPENAME_oneLpNum))
			EGLPNUM_TYPENAME_EGlpNumOne (psinfo->norms[j]);
	}

	EGLPNUM_TYPENAME_EGlpNumCopyFrac (psinfo->norms[eindex], normj, yl);
	EGLPNUM_TYPENAME_EGlpNumDivTo (psinfo->norms[eindex], yl);
	if (EGLPNUM_TYPENAME_EGlpNumIsLess (psinfo->norms[eindex], EGLPNUM_TYPENAME_oneLpNum))
		EGLPNUM_TYPENAME_EGlpNumOne (psinfo->norms[eindex]);

	EGLPNUM_TYPENAME_ILLfct_zero_workvector (lp);
	EGLPNUM_TYPENAME_EGlpNumClearVar (wAj);
	EGLPNUM_TYPENAME_EGlpNumClearVar (zAj);
	EGLPNUM_TYPENAME_EGlpNumClearVar (normj);
	EGLPNUM_TYPENAME_EGlpNumClearVar (ntmp);
}

int EGLPNUM_TYPENAME_ILLprice_build_ddevex_norms (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	EGLPNUM_TYPENAME_d_devex_info * const ddinfo,
	int const reinit)
{
	int i;
	int rval = 0;

	if (reinit == 0)
	{
		ddinfo->ninit = 0;
		ddinfo->norms = EGLPNUM_TYPENAME_EGlpNumAllocArray (lp->nrows);
		ILL_SAFE_MALLOC (ddinfo->refframe, lp->ncols, int);
	}
	if (reinit != 0)
		ddinfo->ninit++;

	for (i = 0; i < lp->ncols; i++)
		ddinfo->refframe[i] = (lp->vstat[i] == STAT_BASIC) ? 1 : 0;

	for (i = 0; i < lp->nrows; i++)
		EGLPNUM_TYPENAME_EGlpNumOne (ddinfo->norms[i]);

CLEANUP:
	if (rval)
	{
		EGLPNUM_TYPENAME_EGlpNumFreeArray (ddinfo->norms);
		ILL_IFFREE (ddinfo->refframe, int);
	}
	EG_RETURN(rval);
}

int EGLPNUM_TYPENAME_ILLprice_update_ddevex_norms (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	EGLPNUM_TYPENAME_d_devex_info * const ddinfo,
	int const lindex,
	EGLPNUM_TYPE yl)
{
	int i, r;
	EGLPNUM_TYPE normi;
	EGLPNUM_TYPE yr;
	EGLPNUM_TYPE ntmp,ntmp2;

	EGLPNUM_TYPENAME_EGlpNumInitVar (ntmp);
	EGLPNUM_TYPENAME_EGlpNumInitVar (ntmp2);
	EGLPNUM_TYPENAME_EGlpNumInitVar (normi);
	EGLPNUM_TYPENAME_EGlpNumInitVar (yr);
	EGLPNUM_TYPENAME_EGlpNumZero (normi);

	for (i = 0; i < lp->zA.nzcnt; i++)
		if (ddinfo->refframe[lp->nbaz[lp->zA.indx[i]]])
			EGLPNUM_TYPENAME_EGlpNumAddInnProdTo (normi, lp->zA.coef[i], lp->zA.coef[i]);

	if (ddinfo->refframe[lp->baz[lindex]])
		EGLPNUM_TYPENAME_EGlpNumAddTo (normi, EGLPNUM_TYPENAME_oneLpNum);

	EGLPNUM_TYPENAME_EGlpNumSet(ntmp,1000.0);
	EGLPNUM_TYPENAME_EGlpNumSet(ntmp2,0.001);
	EGLPNUM_TYPENAME_EGlpNumMultTo(ntmp,ddinfo->norms[lindex]);
	EGLPNUM_TYPENAME_EGlpNumMultTo(ntmp2,ddinfo->norms[lindex]);
	if (EGLPNUM_TYPENAME_EGlpNumIsLess(normi, ntmp2) || EGLPNUM_TYPENAME_EGlpNumIsLess(ntmp, normi))
	{
		EGLPNUM_TYPENAME_EGlpNumClearVar (normi);
		EGLPNUM_TYPENAME_EGlpNumClearVar (yr);
		EGLPNUM_TYPENAME_EGlpNumClearVar (ntmp);
		EGLPNUM_TYPENAME_EGlpNumClearVar (ntmp2);
		return EGLPNUM_TYPENAME_ILLprice_build_ddevex_norms (lp, ddinfo, 1);
	}

	for (i = 0; i < lp->yjz.nzcnt; i++)
	{
		r = lp->yjz.indx[i];
		EGLPNUM_TYPENAME_EGlpNumCopy(yr, lp->yjz.coef[i]);
		EGLPNUM_TYPENAME_EGlpNumCopy(ntmp,yr);
		EGLPNUM_TYPENAME_EGlpNumMultTo(ntmp,yr);
		EGLPNUM_TYPENAME_EGlpNumMultTo(ntmp,normi);
		EGLPNUM_TYPENAME_EGlpNumDivTo(ntmp,yl);
		EGLPNUM_TYPENAME_EGlpNumDivTo(ntmp,yl);
		if (EGLPNUM_TYPENAME_EGlpNumIsLess (ddinfo->norms[r], ntmp))
			EGLPNUM_TYPENAME_EGlpNumCopy (ddinfo->norms[r], ntmp);
	}
	EGLPNUM_TYPENAME_EGlpNumCopy (ddinfo->norms[lindex], normi);
	EGLPNUM_TYPENAME_EGlpNumDivTo(ddinfo->norms[lindex], yl);
	EGLPNUM_TYPENAME_EGlpNumDivTo(ddinfo->norms[lindex], yl);
	if (EGLPNUM_TYPENAME_EGlpNumIsLess (ddinfo->norms[lindex], EGLPNUM_TYPENAME_oneLpNum))
		EGLPNUM_TYPENAME_EGlpNumOne (ddinfo->norms[lindex]);
	EGLPNUM_TYPENAME_EGlpNumClearVar (normi);
	EGLPNUM_TYPENAME_EGlpNumClearVar (yr);
	EGLPNUM_TYPENAME_EGlpNumClearVar (ntmp);
	EGLPNUM_TYPENAME_EGlpNumClearVar (ntmp2);
	return 0;
}

int EGLPNUM_TYPENAME_ILLprice_build_dsteep_norms (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	EGLPNUM_TYPENAME_d_steep_info * const dsinfo)
{
	int i;
	int rval = 0;
	EGLPNUM_TYPENAME_svector z;

	EGLPNUM_TYPENAME_ILLsvector_init (&z);
	rval = EGLPNUM_TYPENAME_ILLsvector_alloc (&z, lp->nrows);
	CHECKRVALG(rval,CLEANUP);
	dsinfo->norms = EGLPNUM_TYPENAME_EGlpNumAllocArray (lp->nrows);

	for (i = 0; i < lp->nrows; i++)
	{
		rval = ILLstring_report (NULL, &lp->O->reporter);
		CHECKRVALG(rval,CLEANUP);

		EGLPNUM_TYPENAME_ILLfct_compute_zz (lp, &z, i);

		EGLPNUM_TYPENAME_EGlpNumInnProd (dsinfo->norms[i], z.coef, z.coef, (size_t) z.nzcnt);
		if (EGLPNUM_TYPENAME_EGlpNumIsLess (dsinfo->norms[i], EGLPNUM_TYPENAME_PARAM_MIN_DNORM))
			EGLPNUM_TYPENAME_EGlpNumCopy (dsinfo->norms[i], EGLPNUM_TYPENAME_PARAM_MIN_DNORM);
	}

CLEANUP:
	EGLPNUM_TYPENAME_ILLsvector_free (&z);
	if (rval)
		EGLPNUM_TYPENAME_EGlpNumFreeArray (dsinfo->norms);

	EG_RETURN(rval);
}

int EGLPNUM_TYPENAME_ILLprice_get_dsteep_norms (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	int const count,
	int *const rowind,
	EGLPNUM_TYPE * const norms)
{
	int i;
	int rval = 0;
	EGLPNUM_TYPENAME_svector z;

	EGLPNUM_TYPENAME_ILLsvector_init (&z);
	rval = EGLPNUM_TYPENAME_ILLsvector_alloc (&z, lp->nrows);
	CHECKRVALG(rval,CLEANUP);

	for (i = 0; i < count; i++)
	{
		EGLPNUM_TYPENAME_ILLfct_compute_zz (lp, &z, rowind[i]);
		EGLPNUM_TYPENAME_EGlpNumInnProd (norms[i], z.coef, z.coef, (size_t) z.nzcnt);
	}

CLEANUP:
	EGLPNUM_TYPENAME_ILLsvector_free (&z);
	EG_RETURN(rval);
}

void EGLPNUM_TYPENAME_ILLprice_update_dsteep_norms (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	EGLPNUM_TYPENAME_d_steep_info * const dsinfo,
	EGLPNUM_TYPENAME_svector * const wz,
	int const lindex,
	EGLPNUM_TYPE yl)
{
	int i, k;
	EGLPNUM_TYPE yij;
	EGLPNUM_TYPE norml;
	EGLPNUM_TYPE *v = 0;
	EGLPNUM_TYPE ntmp;

	EGLPNUM_TYPENAME_EGlpNumInitVar (ntmp);
	EGLPNUM_TYPENAME_EGlpNumInitVar (norml);
	EGLPNUM_TYPENAME_EGlpNumInitVar (yij);
	EGLPNUM_TYPENAME_EGlpNumInnProd (norml, lp->zz.coef, lp->zz.coef, (size_t) (lp->zz.nzcnt));

#if 0
	Bico - remove warnings for dist
		if (fabs ((norml - dsinfo->norms[lindex]) / norml) > 1000.0 /*0.01 */ )
		{
			QSlog("warning: incorrect dnorm values");
			QSlog("anorm = %.6f, pnorm = %.6f", norml, dsinfo->norms[lindex]);
		}
#endif

	EGLPNUM_TYPENAME_ILLfct_load_workvector (lp, wz);
	v = lp->work.coef;

	for (k = 0; k < lp->yjz.nzcnt; k++)
	{
		i = lp->yjz.indx[k];
		EGLPNUM_TYPENAME_EGlpNumCopy (yij, lp->yjz.coef[k]);
		/* compute in ntmp (yij * ((yij * norml / yl) - (2.0 * v[i]))) / yl; */
		EGLPNUM_TYPENAME_EGlpNumCopy(ntmp,yij);
		EGLPNUM_TYPENAME_EGlpNumMultTo(ntmp,norml);
		EGLPNUM_TYPENAME_EGlpNumDivTo(ntmp,yl);
		EGLPNUM_TYPENAME_EGlpNumSubTo(ntmp,v[i]);
		EGLPNUM_TYPENAME_EGlpNumSubTo(ntmp,v[i]);
		EGLPNUM_TYPENAME_EGlpNumMultTo (ntmp, yij);
		EGLPNUM_TYPENAME_EGlpNumDivTo (ntmp, yl);
		/* set dsinfo->norms[i] += (yij * ((yij * norml / yl) - (2.0 * v[i]))) / yl;*/
		EGLPNUM_TYPENAME_EGlpNumAddTo(dsinfo->norms[i], ntmp);
		if (EGLPNUM_TYPENAME_EGlpNumIsLess (dsinfo->norms[i], EGLPNUM_TYPENAME_PARAM_MIN_DNORM))
			EGLPNUM_TYPENAME_EGlpNumCopy (dsinfo->norms[i], EGLPNUM_TYPENAME_PARAM_MIN_DNORM);
	}
	EGLPNUM_TYPENAME_EGlpNumCopyFrac (dsinfo->norms[lindex], norml, yl);
	EGLPNUM_TYPENAME_EGlpNumDivTo (dsinfo->norms[lindex], yl);
	if (EGLPNUM_TYPENAME_EGlpNumIsLess (dsinfo->norms[lindex], EGLPNUM_TYPENAME_PARAM_MIN_DNORM))
		EGLPNUM_TYPENAME_EGlpNumCopy (dsinfo->norms[lindex], EGLPNUM_TYPENAME_PARAM_MIN_DNORM);

	EGLPNUM_TYPENAME_ILLfct_zero_workvector (lp);
	EGLPNUM_TYPENAME_EGlpNumClearVar (norml);
	EGLPNUM_TYPENAME_EGlpNumClearVar (ntmp);
	EGLPNUM_TYPENAME_EGlpNumClearVar (yij);
}

static void update_d_scaleinf (
	EGLPNUM_TYPENAME_price_info * const p,
	EGLPNUM_TYPENAME_heap * const h,
	int const j,
	EGLPNUM_TYPE inf,
	int const prule)
{
	if (!EGLPNUM_TYPENAME_EGlpNumIsNeqqZero (inf))
	{
		EGLPNUM_TYPENAME_EGlpNumZero (p->d_scaleinf[j]);
		if (h->hexist != 0 && h->loc[j] != -1)
			EGLPNUM_TYPENAME_ILLheap_delete (h, j);
	}
	else
	{
		if (prule == QS_PRICE_PDANTZIG)
			EGLPNUM_TYPENAME_EGlpNumCopy (p->d_scaleinf[j], inf);
		else if (prule == QS_PRICE_PDEVEX)
			EGLPNUM_TYPENAME_EGlpNumCopySqrOver (p->d_scaleinf[j], inf, p->pdinfo.norms[j]);
		else if (prule == QS_PRICE_PSTEEP)
			EGLPNUM_TYPENAME_EGlpNumCopySqrOver (p->d_scaleinf[j], inf, p->psinfo.norms[j]);

		if (h->hexist != 0)
		{
			if (h->loc[j] == -1)
				EGLPNUM_TYPENAME_ILLheap_insert (h, j);
			else
				EGLPNUM_TYPENAME_ILLheap_modify (h, j);
		}
	}
}

static void compute_dualI_inf (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	const int j,
	EGLPNUM_TYPE * const inf)
{
	int col = lp->nbaz[j];
	int vt = lp->vtype[col];
	int vs = lp->vstat[col];
	EGLPNUM_TYPE*dj = &(lp->pIdz[j]);
	EGLPNUM_TYPE*ftol = &(lp->tol->id_tol);
	EGLPNUM_TYPENAME_EGlpNumZero (*inf);
	if (vt != VARTIFICIAL && vt != VFIXED)
	{
		if( EGLPNUM_TYPENAME_EGlpNumIsSumLess(*dj,*ftol,EGLPNUM_TYPENAME_zeroLpNum) && (vs == STAT_LOWER || vs == STAT_ZERO))
			EGLPNUM_TYPENAME_EGlpNumCopyNeg(*inf,*dj);
		else if (EGLPNUM_TYPENAME_EGlpNumIsLess(*ftol, *dj) && (vs == STAT_UPPER || vs == STAT_ZERO))
			EGLPNUM_TYPENAME_EGlpNumCopy (*inf, *dj);
	}
}

static void compute_dualII_inf (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	int const j,
	EGLPNUM_TYPE * const inf)
{
	int col = lp->nbaz[j];
	int vt = lp->vtype[col];
	int vs = lp->vstat[col];
	EGLPNUM_TYPE*dj = &(lp->dz[j]);
	EGLPNUM_TYPE*ftol = &(lp->tol->dfeas_tol);
	EGLPNUM_TYPENAME_EGlpNumZero (*inf);
	if (vt != VARTIFICIAL && vt != VFIXED)
	{
		if( EGLPNUM_TYPENAME_EGlpNumIsSumLess(*dj,*ftol,EGLPNUM_TYPENAME_zeroLpNum) && (vs == STAT_LOWER || vs == STAT_ZERO))
			EGLPNUM_TYPENAME_EGlpNumCopyNeg(*inf,*dj);
		else if (EGLPNUM_TYPENAME_EGlpNumIsLess(*ftol,*dj) && (vs == STAT_UPPER || vs == STAT_ZERO))
			EGLPNUM_TYPENAME_EGlpNumCopy (*inf, *dj);
	}
}

void EGLPNUM_TYPENAME_ILLprice_compute_dual_inf (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	EGLPNUM_TYPENAME_price_info * const p,
	int *const ix,
	int const icnt,
	int const phase)
{
	int i;
	int price;
	EGLPNUM_TYPE inf;
	EGLPNUM_TYPENAME_heap *h = &(p->h);

	price = (phase == PRIMAL_PHASEI) ? p->pI_price : p->pII_price;
	EGLPNUM_TYPENAME_EGlpNumInitVar (inf);
	EGLPNUM_TYPENAME_EGlpNumZero (inf);

	if (phase == PRIMAL_PHASEI)
	{
		if (ix == NULL)
			for (i = 0; i < lp->nnbasic; i++)
			{
				compute_dualI_inf (lp, i, &(inf));
				update_d_scaleinf (p, h, i, inf, price);
			}
		else
			for (i = 0; i < icnt; i++)
			{
				compute_dualI_inf (lp, ix[i], &(inf));
				update_d_scaleinf (p, h, ix[i], inf, price);
			}
	}
	else if (phase == PRIMAL_PHASEII)
	{
		if (ix == NULL)
			for (i = 0; i < lp->nnbasic; i++)
			{
				compute_dualII_inf (lp, i, &inf);
				update_d_scaleinf (p, h, i, inf, price);
			}
		else
			for (i = 0; i < icnt; i++)
			{
				compute_dualII_inf (lp, ix[i], &inf);
				update_d_scaleinf (p, h, ix[i], inf, price);
			}
	}
	EGLPNUM_TYPENAME_EGlpNumClearVar (inf);
}

void EGLPNUM_TYPENAME_ILLprice_primal (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	EGLPNUM_TYPENAME_price_info * const pinf,
	EGLPNUM_TYPENAME_price_res * const pr,
	int const phase)
{
	int j, vs;
	EGLPNUM_TYPE d_e, d_max;
	EGLPNUM_TYPE *ftol = &(lp->tol->dfeas_tol);
	EGLPNUM_TYPENAME_heap *const h = &(pinf->h);

	EGLPNUM_TYPENAME_EGlpNumInitVar (d_e);
	EGLPNUM_TYPENAME_EGlpNumInitVar (d_max);
	pr->eindex = -1;
	EGLPNUM_TYPENAME_EGlpNumZero(d_max);

#if USEHEAP > 0
	EGLPNUM_TYPENAME_ILLprice_test_for_heap (lp, pinf, lp->nnbasic, pinf->d_scaleinf,
													PRIMAL_SIMPLEX, 1);
#endif

	if (pinf->p_strategy == COMPLETE_PRICING)
	{
		if (h->hexist)
		{
			pr->eindex = EGLPNUM_TYPENAME_ILLheap_findmin (h);
			if (pr->eindex != -1)
				EGLPNUM_TYPENAME_ILLheap_delete (h, pr->eindex);
		}
		else
		{
			for (j = 0; j < lp->nnbasic; j++)
			{
				if (EGLPNUM_TYPENAME_EGlpNumIsLess (d_max, pinf->d_scaleinf[j]))
				{
					EGLPNUM_TYPENAME_EGlpNumCopy (d_max, pinf->d_scaleinf[j]);
					pr->eindex = j;
				}
			}
		}
	}
	else if (pinf->p_strategy == MULTI_PART_PRICING)
	{
		for (j = 0; j < pinf->pmpinfo.bsize; j++)
		{
			if (EGLPNUM_TYPENAME_EGlpNumIsLess (d_max, pinf->pmpinfo.infeas[j]))
			{
				EGLPNUM_TYPENAME_EGlpNumCopy (d_max, pinf->pmpinfo.infeas[j]);
				pr->eindex = pinf->pmpinfo.bucket[j];
			}
		}
	}

	if (pr->eindex < 0)
		pr->price_stat = PRICE_OPTIMAL;
	else
	{
		if (phase == PRIMAL_PHASEI)
			EGLPNUM_TYPENAME_EGlpNumCopy (d_e, lp->pIdz[pr->eindex]);
		else
			EGLPNUM_TYPENAME_EGlpNumCopy (d_e, lp->dz[pr->eindex]);
		vs = lp->vstat[lp->nbaz[pr->eindex]];

		pr->price_stat = PRICE_NONOPTIMAL;
		if (vs == STAT_UPPER || (vs == STAT_ZERO && EGLPNUM_TYPENAME_EGlpNumIsLess (*ftol, d_e)))
			pr->dir = VDECREASE;
		else
			pr->dir = VINCREASE;
	}
	EGLPNUM_TYPENAME_EGlpNumClearVar (d_e);
	EGLPNUM_TYPENAME_EGlpNumClearVar (d_max);
}

static void update_p_scaleinf (
	EGLPNUM_TYPENAME_price_info * const p,
	EGLPNUM_TYPENAME_heap * const h,
	int const i,
	EGLPNUM_TYPE inf,
	int const prule)
{
	if (!EGLPNUM_TYPENAME_EGlpNumIsNeqqZero (inf))
	{
		EGLPNUM_TYPENAME_EGlpNumZero (p->p_scaleinf[i]);
		if (h->hexist != 0 && h->loc[i] != -1)
			EGLPNUM_TYPENAME_ILLheap_delete (h, i);
	}
	else
	{
		if (prule == QS_PRICE_DDANTZIG)
			EGLPNUM_TYPENAME_EGlpNumCopy (p->p_scaleinf[i], inf);
		else if (prule == QS_PRICE_DSTEEP)
			EGLPNUM_TYPENAME_EGlpNumCopySqrOver (p->p_scaleinf[i], inf, p->dsinfo.norms[i]);
		else if (prule == QS_PRICE_DDEVEX)
			EGLPNUM_TYPENAME_EGlpNumCopySqrOver (p->p_scaleinf[i], inf, p->ddinfo.norms[i]);

		if (h->hexist != 0)
		{
			if (h->loc[i] == -1)
				EGLPNUM_TYPENAME_ILLheap_insert (h, i);
			else
				EGLPNUM_TYPENAME_ILLheap_modify (h, i);
		}
	}
}

static void compute_primalI_inf (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	int const i,
	EGLPNUM_TYPE * const inf)
{
	int const col = lp->baz[i];
	EGLPNUM_TYPE*x = &(lp->xbz[i]);
	EGLPNUM_TYPE*l = &(lp->lz[col]);
	EGLPNUM_TYPE*u = &(lp->uz[col]);
	EGLPNUM_TYPE*ftol = &(lp->tol->ip_tol);
	EGLPNUM_TYPENAME_EGlpNumZero (*inf);

	if (EGLPNUM_TYPENAME_EGlpNumIsLess (*ftol, *x) && EGLPNUM_TYPENAME_EGlpNumIsNeqq (*u, EGLPNUM_TYPENAME_INFTY))
		EGLPNUM_TYPENAME_EGlpNumCopy (*inf, *x);
	else if (EGLPNUM_TYPENAME_EGlpNumIsNeqq (*l, EGLPNUM_TYPENAME_NINFTY) && EGLPNUM_TYPENAME_EGlpNumIsSumLess (*x, *ftol,EGLPNUM_TYPENAME_zeroLpNum))
		EGLPNUM_TYPENAME_EGlpNumCopy (*inf, *x);
}

static void compute_primalII_inf (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	int const i,
	EGLPNUM_TYPE * const inf)
{
	int const col = lp->baz[i];
	EGLPNUM_TYPE*x = &(lp->xbz[i]);
	EGLPNUM_TYPE*l = &(lp->lz[col]);
	EGLPNUM_TYPE*u = &(lp->uz[col]);
	EGLPNUM_TYPE*ftol = &(lp->tol->pfeas_tol);
	EGLPNUM_TYPENAME_EGlpNumZero (*inf);

	if (EGLPNUM_TYPENAME_EGlpNumIsNeqq (*u, EGLPNUM_TYPENAME_INFTY) && EGLPNUM_TYPENAME_EGlpNumIsSumLess (*u, *ftol, *x))
		EGLPNUM_TYPENAME_EGlpNumCopyDiff (*inf, *x, *u);
	else if (EGLPNUM_TYPENAME_EGlpNumIsNeqq (*l, EGLPNUM_TYPENAME_NINFTY) && EGLPNUM_TYPENAME_EGlpNumIsSumLess (*x, *ftol, *l))
		EGLPNUM_TYPENAME_EGlpNumCopyDiff (*inf, *l, *x);
}

void EGLPNUM_TYPENAME_ILLprice_compute_primal_inf (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	EGLPNUM_TYPENAME_price_info * const p,
	int *const ix,
	int const icnt,
	int const phase)
{
	int i;
	int price;
	EGLPNUM_TYPE inf;
	EGLPNUM_TYPENAME_heap *h = &(p->h);

	price = (phase == DUAL_PHASEI) ? p->dI_price : p->dII_price;
	EGLPNUM_TYPENAME_EGlpNumInitVar (inf);
	EGLPNUM_TYPENAME_EGlpNumZero (inf);

	if (phase == DUAL_PHASEI)
	{
		if (ix == NULL)
			for (i = 0; i < lp->nrows; i++)
			{
				compute_primalI_inf (lp, i, &inf);
				update_p_scaleinf (p, h, i, inf, price);
			}
		else
			for (i = 0; i < icnt; i++)
			{
				compute_primalI_inf (lp, ix[i], &inf);
				update_p_scaleinf (p, h, ix[i], inf, price);
			}
	}
	else if (phase == DUAL_PHASEII)
	{
		if (ix == NULL)
			for (i = 0; i < lp->nrows; i++)
			{
				compute_primalII_inf (lp, i, &inf);
				update_p_scaleinf (p, h, i, inf, price);
			}
		else
			for (i = 0; i < icnt; i++)
			{
				compute_primalII_inf (lp, ix[i], &inf);
				update_p_scaleinf (p, h, ix[i], inf, price);
			}
	}
	EGLPNUM_TYPENAME_EGlpNumClearVar (inf);
}

void EGLPNUM_TYPENAME_ILLprice_dual (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	EGLPNUM_TYPENAME_price_info * const pinf,
	int const phase,
	EGLPNUM_TYPENAME_price_res * const pr)
{
	int i;
	EGLPNUM_TYPE p_max;
	EGLPNUM_TYPE ubound;
	EGLPNUM_TYPE*ftol = &(lp->tol->pfeas_tol);
	EGLPNUM_TYPENAME_heap *const h = &(pinf->h);

	EGLPNUM_TYPENAME_EGlpNumInitVar (p_max);
	EGLPNUM_TYPENAME_EGlpNumInitVar (ubound);
	pr->lindex = -1;
	EGLPNUM_TYPENAME_EGlpNumZero(p_max);

#if USEHEAP > 0
	EGLPNUM_TYPENAME_ILLprice_test_for_heap (lp, pinf, lp->nrows, pinf->p_scaleinf, DUAL_SIMPLEX,
													1);
#endif

	if (pinf->d_strategy == COMPLETE_PRICING)
	{
		if (h->hexist)
		{
			pr->lindex = EGLPNUM_TYPENAME_ILLheap_findmin (h);
			if (pr->lindex != -1)
				EGLPNUM_TYPENAME_ILLheap_delete (h, pr->lindex);
		}
		else
		{
			for (i = 0; i < lp->nrows; i++)
			{
				if (EGLPNUM_TYPENAME_EGlpNumIsLess (p_max, pinf->p_scaleinf[i]))
				{
					EGLPNUM_TYPENAME_EGlpNumCopy (p_max, pinf->p_scaleinf[i]);
					pr->lindex = i;
				}
			}
		}
	}
	else if (pinf->d_strategy == MULTI_PART_PRICING)
	{
		for (i = 0; i < pinf->dmpinfo.bsize; i++)
		{
			if (EGLPNUM_TYPENAME_EGlpNumIsLess (p_max, pinf->dmpinfo.infeas[i]))
			{
				EGLPNUM_TYPENAME_EGlpNumCopy (p_max, pinf->dmpinfo.infeas[i]);
				pr->lindex = pinf->dmpinfo.bucket[i];
			}
		}
	}

	if (pr->lindex < 0)
		pr->price_stat = PRICE_OPTIMAL;
	else
	{
		pr->price_stat = NONOPTIMAL;

		if (EGLPNUM_TYPENAME_EGlpNumIsNeqq (lp->uz[lp->baz[pr->lindex]], EGLPNUM_TYPENAME_INFTY))
		{
			if (phase == DUAL_PHASEI)
				EGLPNUM_TYPENAME_EGlpNumZero(ubound);
			else
				EGLPNUM_TYPENAME_EGlpNumCopy(ubound,lp->uz[lp->baz[pr->lindex]]);
			if (EGLPNUM_TYPENAME_EGlpNumIsSumLess (*ftol, ubound, lp->xbz[pr->lindex]))
				pr->lvstat = STAT_UPPER;
			else
				pr->lvstat = STAT_LOWER;
		}
		else
			pr->lvstat = STAT_LOWER;
	}
	EGLPNUM_TYPENAME_EGlpNumClearVar (p_max);
	EGLPNUM_TYPENAME_EGlpNumClearVar (ubound);
}

int EGLPNUM_TYPENAME_ILLprice_get_rownorms (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	EGLPNUM_TYPENAME_price_info * const pinf,
	EGLPNUM_TYPE * const rnorms)
{
	int rval = 0;
	int i;

	if (pinf->dsinfo.norms == NULL)
	{
		rval = EGLPNUM_TYPENAME_ILLprice_build_dsteep_norms (lp, &(pinf->dsinfo));
		CHECKRVALG(rval,CLEANUP);
	}
	for (i = 0; i < lp->nrows; i++)
		EGLPNUM_TYPENAME_EGlpNumCopy (rnorms[i], pinf->dsinfo.norms[i]);

CLEANUP:
	if (rval)
		EGLPNUM_TYPENAME_EGlpNumFreeArray (pinf->dsinfo.norms);

	return rval;
}

int EGLPNUM_TYPENAME_ILLprice_get_colnorms (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	EGLPNUM_TYPENAME_price_info * const pinf,
	EGLPNUM_TYPE * const cnorms)
{
	int rval = 0;
	int i, j;

	if (pinf->psinfo.norms == NULL)
	{
		rval = EGLPNUM_TYPENAME_ILLprice_build_psteep_norms (lp, &(pinf->psinfo));
		CHECKRVALG(rval,CLEANUP);
	}
	for (i = 0; i < lp->nrows; i++)
		EGLPNUM_TYPENAME_EGlpNumZero (cnorms[lp->baz[i]]);
	for (j = 0; j < lp->nnbasic; j++)
		EGLPNUM_TYPENAME_EGlpNumCopy (cnorms[lp->nbaz[j]], pinf->psinfo.norms[j]);

CLEANUP:
	if (rval)
		EGLPNUM_TYPENAME_EGlpNumFreeArray (pinf->psinfo.norms);

	return rval;
}

int EGLPNUM_TYPENAME_ILLprice_get_newnorms (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	int const nelems,
	EGLPNUM_TYPE * const norms,
	int *const matcnt,
	int *const matbeg,
	int *const matind,
	EGLPNUM_TYPE * const matval,
	int const option)
{
	int i, j;
	int rval = 0;
	EGLPNUM_TYPENAME_svector a;
	EGLPNUM_TYPENAME_svector y;

	EGLPNUM_TYPENAME_ILLsvector_init (&y);
	rval = EGLPNUM_TYPENAME_ILLsvector_alloc (&y, lp->nrows);
	CHECKRVALG(rval,CLEANUP);

	for (j = 0; j < nelems; j++)
	{
		a.nzcnt = matcnt[j];
		a.indx = &(matind[matbeg[j]]);
		a.coef = &(matval[matbeg[j]]);

		if (option == COLUMN_SOLVE)
			EGLPNUM_TYPENAME_ILLbasis_column_solve (lp, &a, &y);
		else
			EGLPNUM_TYPENAME_ILLbasis_row_solve (lp, &a, &y);

		EGLPNUM_TYPENAME_EGlpNumOne (norms[j]);
		for (i = 0; i < y.nzcnt; i++)
			EGLPNUM_TYPENAME_EGlpNumAddInnProdTo (norms[j], y.coef[i], y.coef[i]);
	}

CLEANUP:
	EGLPNUM_TYPENAME_ILLsvector_free (&y);
	EG_RETURN(rval);
}

int EGLPNUM_TYPENAME_ILLprice_get_new_rownorms (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	int const newrows,
	EGLPNUM_TYPE * const rnorms,
	int *const rmatcnt,
	int *const rmatbeg,
	int *const rmatind,
	EGLPNUM_TYPE * const rmatval)
{
	return EGLPNUM_TYPENAME_ILLprice_get_newnorms (lp, newrows, rnorms, rmatcnt, rmatbeg, rmatind,
																rmatval, ROW_SOLVE);
}

int EGLPNUM_TYPENAME_ILLprice_get_new_colnorms (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	int const newrows,
	EGLPNUM_TYPE * const rnorms,
	int *const matcnt,
	int *const matbeg,
	int *const matind,
	EGLPNUM_TYPE * const matval)
{
	return EGLPNUM_TYPENAME_ILLprice_get_newnorms (lp, newrows, rnorms, matcnt, matbeg, matind,
																matval, COLUMN_SOLVE);
}

int EGLPNUM_TYPENAME_ILLprice_load_rownorms (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	EGLPNUM_TYPE * const rnorms,
	EGLPNUM_TYPENAME_price_info * const pinf)
{
	int i;
	int rval = 0;

	EGLPNUM_TYPENAME_EGlpNumFreeArray (pinf->dsinfo.norms);
	pinf->dsinfo.norms = EGLPNUM_TYPENAME_EGlpNumAllocArray (lp->nrows);

	for (i = 0; i < lp->nrows; i++)
	{
		EGLPNUM_TYPENAME_EGlpNumCopy (pinf->dsinfo.norms[i], rnorms[i]);
		if (EGLPNUM_TYPENAME_EGlpNumIsLess (pinf->dsinfo.norms[i], EGLPNUM_TYPENAME_PARAM_MIN_DNORM))
			EGLPNUM_TYPENAME_EGlpNumCopy (pinf->dsinfo.norms[i], EGLPNUM_TYPENAME_PARAM_MIN_DNORM);
	}

	EG_RETURN(rval);
}

int EGLPNUM_TYPENAME_ILLprice_load_colnorms (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	EGLPNUM_TYPE * const cnorms,
	EGLPNUM_TYPENAME_price_info * const pinf)
{
	int j;
	int rval = 0;

	EGLPNUM_TYPENAME_EGlpNumFreeArray (pinf->psinfo.norms);
	pinf->psinfo.norms = EGLPNUM_TYPENAME_EGlpNumAllocArray (lp->nnbasic);

	for (j = 0; j < lp->nnbasic; j++)
	{
		EGLPNUM_TYPENAME_EGlpNumCopy (pinf->psinfo.norms[j], cnorms[lp->nbaz[j]]);
		if (EGLPNUM_TYPENAME_EGlpNumIsLess (pinf->psinfo.norms[j], EGLPNUM_TYPENAME_oneLpNum))
			EGLPNUM_TYPENAME_EGlpNumOne (pinf->psinfo.norms[j]);
	}

	EG_RETURN(rval);
}

#if PRICE_DEBUG > 0
void EGLPNUM_TYPENAME_test_dsteep_norms (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_price_info * p)
{
	int i, errn = 0;
	EGLPNUM_TYPE *pn = EGLPNUM_TYPENAME_EGlpNumAllocArray(lp->nrows);
	EGLPNUM_TYPE err, diff;
	EGLPNUM_TYPENAME_EGlpNumZero (err);

	EGLPNUM_TYPENAME_EGlpNumInitVar (err);
	EGLPNUM_TYPENAME_EGlpNumInitVar (diff);

	EGLPNUM_TYPENAME_ILLprice_get_dsteep_norms (lp, lp->yjz.nzcnt, lp->yjz.indx, pn);
	for (i = 0; i < lp->yjz.nzcnt; i++)
	{
		EGLPNUM_TYPENAME_EGlpNumCopyDiff (diff, pn[i], p->dsinfo.norms[lp->yjz.indx[i]]);
		EGLPNUM_TYPENAME_EGlpNumCopyAbs(diff,diff);
		if (EGLPNUM_TYPENAME_EGlpNumIsLess (EGLPNUM_TYPENAME_PFEAS_TOLER, diff))
		{
			errn++;
			EGLPNUM_TYPENAME_EGlpNumAddTo (err, diff);
			EGLPNUM_TYPENAME_EGlpNumCopy (p->dsinfo.norms[lp->yjz.indx[i]], pn[i]);
		}
	}
	if (errn)
		QSlog("%d: dnorm errn = %d, err = %.6f", lp->cnts->tot_iter, errn,
								EGLPNUM_TYPENAME_EGlpNumToLf (err));
	EGLPNUM_TYPENAME_EGlpNumFreeArray (pn);
	EGLPNUM_TYPENAME_EGlpNumClearVar (diff);
	EGLPNUM_TYPENAME_EGlpNumClearVar (err);
}
#endif
