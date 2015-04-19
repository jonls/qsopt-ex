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

/* RCS_INFO = "$RCSfile: fct.c,v $ $Revision: 1.2 $ $Date: 2003/11/05 16:49:52 $"; */
static int TRACE = 0;

//#define FCT_DEBUG 10
#define FCT_DEBUG 0

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "qs_config.h"
#include "logging-private.h"

#include "allocrus.h"
#include "eg_lpnum.h"
#include "eg_io.h"
#include "except.h"
#include "trace.h"
#include "util.h"

#include "lpdefs_EGLPNUM_TYPENAME.h"
#include "stddefs.h"
#include "basis_EGLPNUM_TYPENAME.h"
#include "fct_EGLPNUM_TYPENAME.h"
#include "price_EGLPNUM_TYPENAME.h"
#include "ratio_EGLPNUM_TYPENAME.h"
#include "dstruct_EGLPNUM_TYPENAME.h"


EGLPNUM_TYPENAME_bndinfo *EGLPNUM_TYPENAME_ILLfct_new_bndinfo (
	void)
{
	EGLPNUM_TYPENAME_bndinfo *nbnd = (EGLPNUM_TYPENAME_bndinfo *) malloc (sizeof (EGLPNUM_TYPENAME_bndinfo));

	if (!nbnd)
	{
		QSlog("not enough memory, in %s", __func__);
		exit (1);
	}
	EGLPNUM_TYPENAME_EGlpNumInitVar ((nbnd->pbound));
	EGLPNUM_TYPENAME_EGlpNumInitVar ((nbnd->cbound));
	return nbnd;
}

void EGLPNUM_TYPENAME_ILLfct_free_bndinfo (
	EGLPNUM_TYPENAME_bndinfo * binfo)
{
	EGLPNUM_TYPENAME_EGlpNumClearVar ((binfo->pbound));
	EGLPNUM_TYPENAME_EGlpNumClearVar ((binfo->cbound));
	ILL_IFFREE (binfo, EGLPNUM_TYPENAME_bndinfo);
	return;
}

static int compute_zA1 (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_svector * z,
	EGLPNUM_TYPENAME_svector * zA,
	EGLPNUM_TYPE ztoler),
/*
  compute_zA2 (EGLPNUM_TYPENAME_lpinfo * lp,
							 EGLPNUM_TYPENAME_svector * z,
							 EGLPNUM_TYPENAME_svector * zA,
							 const EGLPNUM_TYPE* ztoler), */
  compute_zA3 (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_svector * z,
	EGLPNUM_TYPENAME_svector * zA,
	EGLPNUM_TYPE ztoler),
  expand_var_bounds (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPE ftol,
	int *chgb),
  expand_var_coefs (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPE ftol,
	int *chgc);

static void update_piv_values (
	EGLPNUM_TYPENAME_count_struct * c,
	int phase,
	const EGLPNUM_TYPE piv),
/*  copy_vectors (EGLPNUM_TYPENAME_svector * a,
								EGLPNUM_TYPENAME_svector * b),*/
  add_vectors (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_svector * a,
	EGLPNUM_TYPENAME_svector * b,
	EGLPNUM_TYPENAME_svector * c,
	const EGLPNUM_TYPE t);

static double my_rand (
	int bound,
	ILLrandstate * r);


void EGLPNUM_TYPENAME_ILLfct_load_workvector (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_svector * s)
{
	int i;

	for (i = 0; i < s->nzcnt; i++)
	{
		lp->work.indx[i] = s->indx[i];
		EGLPNUM_TYPENAME_EGlpNumCopy (lp->work.coef[s->indx[i]], s->coef[i]);
	}
	lp->work.nzcnt = s->nzcnt;
}

void EGLPNUM_TYPENAME_ILLfct_zero_workvector (
	EGLPNUM_TYPENAME_lpinfo * lp)
{
	int i;

	for (i = 0; i < lp->work.nzcnt; i++)
		EGLPNUM_TYPENAME_EGlpNumZero (lp->work.coef[lp->work.indx[i]]);
	lp->work.nzcnt = 0;
}

void EGLPNUM_TYPENAME_ILLfct_set_variable_type (
	EGLPNUM_TYPENAME_lpinfo * lp)
{
	int j;

	for (j = 0; j < lp->ncols; j++)
	{

		if (lp->matcnt[j] == 1 && lp->O->rowmap[lp->matind[lp->matbeg[j]]] == j)
			lp->vclass[j] = CLASS_LOGICAL;
		else
			lp->vclass[j] = CLASS_STRUCT;
		switch ((EGLPNUM_TYPENAME_EGlpNumIsEqqual (lp->uz[j], EGLPNUM_TYPENAME_INFTY) ? 1U : 0U) |
						(EGLPNUM_TYPENAME_EGlpNumIsEqqual (lp->lz[j], EGLPNUM_TYPENAME_NINFTY) ? 2U : 0U))
		{
		case 0:
			if (EGLPNUM_TYPENAME_EGlpNumIsLess (lp->lz[j], lp->uz[j]))
				lp->vtype[j] = VBOUNDED;
			else if (!EGLPNUM_TYPENAME_EGlpNumIsNeqqZero (lp->lz[j]) &&
							 (lp->vclass[j] == CLASS_LOGICAL))
				lp->vtype[j] = VARTIFICIAL;
			else
				lp->vtype[j] = VFIXED;
			break;
		case 3:
			lp->vtype[j] = VFREE;
			break;
		case 1:
			lp->vtype[j] = VLOWER;
			break;
		case 2:
			lp->vtype[j] = VUPPER;
			break;
		}
	}
}

/* compute various vectors */

void EGLPNUM_TYPENAME_ILLfct_compute_pobj (
	EGLPNUM_TYPENAME_lpinfo * lp)
{
	int i, j;
	int col;
	EGLPNUM_TYPE sum;

	EGLPNUM_TYPENAME_EGlpNumInitVar (sum);
	EGLPNUM_TYPENAME_EGlpNumZero (sum);

	for (i = 0; i < lp->nrows; i++)
		EGLPNUM_TYPENAME_EGlpNumAddInnProdTo (sum, lp->cz[lp->baz[i]], lp->xbz[i]);

	for (j = 0; j < lp->nnbasic; j++)
	{
		col = lp->nbaz[j];
		if (lp->vstat[col] == STAT_UPPER)
			EGLPNUM_TYPENAME_EGlpNumAddInnProdTo (sum, lp->cz[col], lp->uz[col]);
		else if (lp->vstat[col] == STAT_LOWER)
			EGLPNUM_TYPENAME_EGlpNumAddInnProdTo (sum, lp->cz[col], lp->lz[col]);
	}
	EGLPNUM_TYPENAME_EGlpNumCopy (lp->pobjval, sum);
	EGLPNUM_TYPENAME_EGlpNumCopy (lp->objval, sum);
	EGLPNUM_TYPENAME_EGlpNumClearVar (sum);
}

void EGLPNUM_TYPENAME_ILLfct_compute_dobj (
	EGLPNUM_TYPENAME_lpinfo * lp)
{
	int i, j;
	int col;
	EGLPNUM_TYPE sum;

	EGLPNUM_TYPENAME_EGlpNumInitVar (sum);
	EGLPNUM_TYPENAME_EGlpNumZero (sum);

	for (i = 0; i < lp->nrows; i++)
		EGLPNUM_TYPENAME_EGlpNumAddInnProdTo (sum, lp->piz[i], lp->bz[i]);

	for (j = 0; j < lp->nnbasic; j++)
	{
		col = lp->nbaz[j];
		if (lp->vstat[col] == STAT_UPPER)
			EGLPNUM_TYPENAME_EGlpNumAddInnProdTo (sum, lp->dz[j], lp->uz[col]);
		else if (lp->vstat[col] == STAT_LOWER)
			EGLPNUM_TYPENAME_EGlpNumAddInnProdTo (sum, lp->dz[j], lp->lz[col]);
	}
	EGLPNUM_TYPENAME_EGlpNumCopy (lp->dobjval, sum);
	EGLPNUM_TYPENAME_EGlpNumCopy (lp->objval, sum);
	EGLPNUM_TYPENAME_EGlpNumClearVar (sum);
}

void EGLPNUM_TYPENAME_ILLfct_compute_xbz (
	EGLPNUM_TYPENAME_lpinfo * lp)
{
	int i, j, r;
	int col, mcnt, mbeg;
	EGLPNUM_TYPENAME_svector *srhs = &(lp->srhs);
	EGLPNUM_TYPENAME_svector *ssoln = &(lp->ssoln);
	EGLPNUM_TYPE xval;

	EGLPNUM_TYPENAME_EGlpNumInitVar (xval);

	for (i = 0; i < lp->nrows; i++)
	{
		EGLPNUM_TYPENAME_EGlpNumZero (lp->xbz[i]);
		EGLPNUM_TYPENAME_EGlpNumCopy (srhs->coef[i], lp->bz[i]);
	}
	for (j = 0; j < lp->nnbasic; j++)
	{
		col = lp->nbaz[j];
		EGLPNUM_TYPENAME_EGlpNumZero (xval);
		if (lp->vstat[col] == STAT_UPPER && EGLPNUM_TYPENAME_EGlpNumIsNeqqZero (lp->uz[col]))
			EGLPNUM_TYPENAME_EGlpNumCopy (xval, lp->uz[col]);
		else if (lp->vstat[col] == STAT_LOWER && EGLPNUM_TYPENAME_EGlpNumIsNeqqZero (lp->lz[col]))
			EGLPNUM_TYPENAME_EGlpNumCopy (xval, lp->lz[col]);

		if (EGLPNUM_TYPENAME_EGlpNumIsNeqqZero (xval))
		{
			mcnt = lp->matcnt[col];
			mbeg = lp->matbeg[col];
			for (i = 0; i < mcnt; i++)
				EGLPNUM_TYPENAME_EGlpNumSubInnProdTo (srhs->coef[lp->matind[mbeg + i]], xval,
														 lp->matval[mbeg + i]);
		}
	}
	for (i = 0, r = 0; i < lp->nrows; i++)
		if (EGLPNUM_TYPENAME_EGlpNumIsNeqqZero (srhs->coef[i]))
		{
			EGLPNUM_TYPENAME_EGlpNumCopy (srhs->coef[r], srhs->coef[i]);
			srhs->indx[r] = i;
			r++;
		}
	srhs->nzcnt = r;

	EGLPNUM_TYPENAME_ILLbasis_column_solve (lp, srhs, ssoln);
	for (i = 0; i < ssoln->nzcnt; i++)
		EGLPNUM_TYPENAME_EGlpNumCopy (lp->xbz[ssoln->indx[i]], ssoln->coef[i]);
	EGLPNUM_TYPENAME_EGlpNumClearVar (xval);
}

void EGLPNUM_TYPENAME_ILLfct_compute_piz (
	EGLPNUM_TYPENAME_lpinfo * lp)
{
	int i, r;
	EGLPNUM_TYPENAME_svector *srhs = &(lp->srhs);
	EGLPNUM_TYPENAME_svector *ssoln = &(lp->ssoln);

	for (i = 0, r = 0; i < lp->nrows; i++)
	{
		EGLPNUM_TYPENAME_EGlpNumZero (lp->piz[i]);
		if (EGLPNUM_TYPENAME_EGlpNumIsNeqqZero (lp->cz[lp->baz[i]]))
		{
			srhs->indx[r] = i;
			EGLPNUM_TYPENAME_EGlpNumCopy (srhs->coef[r], lp->cz[lp->baz[i]]);
			r++;
		}
	}
	srhs->nzcnt = r;

	EGLPNUM_TYPENAME_ILLbasis_row_solve (lp, srhs, ssoln);
	for (i = 0; i < ssoln->nzcnt; i++)
		EGLPNUM_TYPENAME_EGlpNumCopy (lp->piz[ssoln->indx[i]], ssoln->coef[i]);
}

void EGLPNUM_TYPENAME_ILLfct_compute_dz (
	EGLPNUM_TYPENAME_lpinfo * lp)
{
	int i, j;
	int col;
	int mcnt, mbeg;
	EGLPNUM_TYPE sum;

	EGLPNUM_TYPENAME_EGlpNumInitVar (sum);

	for (j = 0; j < lp->nnbasic; j++)
	{
		EGLPNUM_TYPENAME_EGlpNumZero (sum);
		col = lp->nbaz[j];
		mcnt = lp->matcnt[col];
		mbeg = lp->matbeg[col];
		for (i = 0; i < mcnt; i++)
			EGLPNUM_TYPENAME_EGlpNumAddInnProdTo (sum, lp->piz[lp->matind[mbeg + i]],
													 lp->matval[mbeg + i]);
		EGLPNUM_TYPENAME_EGlpNumCopyDiff (lp->dz[j], lp->cz[col], sum);
	}
	EGLPNUM_TYPENAME_EGlpNumClearVar (sum);
}

void EGLPNUM_TYPENAME_ILLfct_compute_phaseI_xbz (
	EGLPNUM_TYPENAME_lpinfo * lp)
{
	int i, j, r;
	int col, mcnt, mbeg;
	EGLPNUM_TYPENAME_svector *srhs = &(lp->srhs);
	EGLPNUM_TYPENAME_svector *ssoln = &(lp->ssoln);

	for (i = 0; i < lp->nrows; i++)
	{
		EGLPNUM_TYPENAME_EGlpNumZero (lp->xbz[i]);
		EGLPNUM_TYPENAME_EGlpNumZero (srhs->coef[i]);
	}
	for (j = 0; j < lp->nnbasic; j++)
	{
		col = lp->nbaz[j];

		if (lp->dfeas[j])
		{
			mcnt = lp->matcnt[col];
			mbeg = lp->matbeg[col];
			if (lp->dfeas[j] == -1)
				for (i = 0; i < mcnt; i++)
					EGLPNUM_TYPENAME_EGlpNumSubTo (srhs->coef[lp->matind[mbeg + i]], lp->matval[mbeg + i]);
			else
				for (i = 0; i < mcnt; i++)
					EGLPNUM_TYPENAME_EGlpNumAddTo (srhs->coef[lp->matind[mbeg + i]], lp->matval[mbeg + i]);
		}
	}
	for (i = 0, r = 0; i < lp->nrows; i++)
		if (EGLPNUM_TYPENAME_EGlpNumIsNeqqZero (srhs->coef[i]))
		{
			EGLPNUM_TYPENAME_EGlpNumCopy (srhs->coef[r], srhs->coef[i]);
			srhs->indx[r] = i;
			r++;
		}
	srhs->nzcnt = r;

	EGLPNUM_TYPENAME_ILLbasis_column_solve (lp, srhs, ssoln);
	for (i = 0; i < ssoln->nzcnt; i++)
		EGLPNUM_TYPENAME_EGlpNumCopy (lp->xbz[ssoln->indx[i]], ssoln->coef[i]);
}

void EGLPNUM_TYPENAME_ILLfct_compute_phaseI_piz (
	EGLPNUM_TYPENAME_lpinfo * lp)
{
	int i, r;
	EGLPNUM_TYPENAME_svector *srhs = &(lp->srhs);
	EGLPNUM_TYPENAME_svector *ssoln = &(lp->ssoln);

	for (i = 0, r = 0; i < lp->nrows; i++)
	{
		EGLPNUM_TYPENAME_EGlpNumZero (lp->pIpiz[i]);
		if (lp->bfeas[i] != 0)
		{
			srhs->indx[r] = i;
			EGLPNUM_TYPENAME_EGlpNumSet (srhs->coef[r], (double) lp->bfeas[i]);
			r++;
		}
	}
	srhs->nzcnt = r;

	EGLPNUM_TYPENAME_ILLbasis_row_solve (lp, srhs, ssoln);
	for (i = 0; i < ssoln->nzcnt; i++)
		EGLPNUM_TYPENAME_EGlpNumCopy (lp->pIpiz[ssoln->indx[i]], ssoln->coef[i]);
	EGLPNUM_TYPENAME_ILLfct_update_counts (lp, CNT_P1PINZ, ssoln->nzcnt, EGLPNUM_TYPENAME_zeroLpNum);
}

void EGLPNUM_TYPENAME_ILLfct_compute_phaseI_dz (
	EGLPNUM_TYPENAME_lpinfo * lp)
{
	int i, j;
	int col;
	int mcnt, mbeg;
	EGLPNUM_TYPE sum;

	EGLPNUM_TYPENAME_EGlpNumInitVar (sum);
	ILL_IFTRACE ("%s\n", __func__);

	for (j = 0; j < lp->nnbasic; j++)
	{
		EGLPNUM_TYPENAME_EGlpNumZero (sum);
		col = lp->nbaz[j];
		mcnt = lp->matcnt[col];
		mbeg = lp->matbeg[col];
		for (i = 0; i < mcnt; i++)
			EGLPNUM_TYPENAME_EGlpNumAddInnProdTo (sum, lp->pIpiz[lp->matind[mbeg + i]],
													 lp->matval[mbeg + i]);
		EGLPNUM_TYPENAME_EGlpNumCopyNeg (lp->pIdz[j], sum);
		ILL_IFTRACE ("%d:%d:%lf:%la\n", j, col, EGLPNUM_TYPENAME_EGlpNumToLf (sum),
								 EGLPNUM_TYPENAME_EGlpNumToLf (sum));
	}
	EGLPNUM_TYPENAME_EGlpNumClearVar (sum);
}

void EGLPNUM_TYPENAME_ILLfct_compute_yz (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_svector * yz,
	EGLPNUM_TYPENAME_svector * updz,
	int col)
{
	EGLPNUM_TYPENAME_svector a;

	a.nzcnt = lp->matcnt[col];
	a.indx = &(lp->matind[lp->matbeg[col]]);
	a.coef = &(lp->matval[lp->matbeg[col]]);

	EGLPNUM_TYPENAME_ILLfactor_set_factor_dparam (lp->f, QS_FACTOR_SZERO_TOL, EGLPNUM_TYPENAME_PIVZ_TOLER);
	if (updz)
		EGLPNUM_TYPENAME_ILLbasis_column_solve_update (lp, &a, updz, yz);
	else
		EGLPNUM_TYPENAME_ILLbasis_column_solve (lp, &a, yz);
	EGLPNUM_TYPENAME_ILLfactor_set_factor_dparam (lp->f, QS_FACTOR_SZERO_TOL, EGLPNUM_TYPENAME_SZERO_TOLER);
}

void EGLPNUM_TYPENAME_ILLfct_compute_zz (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_svector * zz,
	int row)
{
	EGLPNUM_TYPENAME_ILLfct_compute_binvrow (lp, zz, row, EGLPNUM_TYPENAME_PIVZ_TOLER);
}

void EGLPNUM_TYPENAME_ILLfct_compute_binvrow (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_svector * zz,
	int row,
	EGLPNUM_TYPE ztoler)
{
	EGLPNUM_TYPENAME_svector a;
	EGLPNUM_TYPE e;

	EGLPNUM_TYPENAME_EGlpNumInitVar (e);
	EGLPNUM_TYPENAME_EGlpNumOne (e);

	a.nzcnt = 1;
	a.coef = &e;
	a.indx = &row;

	if (EGLPNUM_TYPENAME_EGlpNumIsGreatZero (ztoler))
		EGLPNUM_TYPENAME_ILLfactor_set_factor_dparam (lp->f, QS_FACTOR_SZERO_TOL, ztoler);
	EGLPNUM_TYPENAME_ILLbasis_row_solve (lp, &a, zz);
	if (EGLPNUM_TYPENAME_EGlpNumIsGreatZero (ztoler))
		EGLPNUM_TYPENAME_ILLfactor_set_factor_dparam (lp->f, QS_FACTOR_SZERO_TOL, EGLPNUM_TYPENAME_SZERO_TOLER);
	EGLPNUM_TYPENAME_EGlpNumClearVar (e);
}

void EGLPNUM_TYPENAME_ILLfct_compute_psteep_upv (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_svector * swz)
{
	EGLPNUM_TYPENAME_ILLbasis_row_solve (lp, &(lp->yjz), swz);
}

void EGLPNUM_TYPENAME_ILLfct_compute_dsteep_upv (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_svector * swz)
{
	EGLPNUM_TYPENAME_ILLbasis_column_solve (lp, &(lp->zz), swz);
}

static int compute_zA1 (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_svector * z,
	EGLPNUM_TYPENAME_svector * zA,
	EGLPNUM_TYPE ztoler)
{
	int rval = 0;
	int i, j, nz = 0;
	int col, mcnt, mbeg;
	EGLPNUM_TYPE sum;
	EGLPNUM_TYPE *v = 0;

	EGLPNUM_TYPENAME_EGlpNumInitVar (sum);
	v = EGLPNUM_TYPENAME_EGlpNumAllocArray (lp->nrows);

	for (i = 0; i < lp->nrows; i++)
		EGLPNUM_TYPENAME_EGlpNumZero (v[i]);
	for (i = 0; i < z->nzcnt; i++)
		EGLPNUM_TYPENAME_EGlpNumCopy (v[z->indx[i]], z->coef[i]);

	for (j = 0; j < lp->nnbasic; j++)
	{
		EGLPNUM_TYPENAME_EGlpNumZero (sum);
		col = lp->nbaz[j];
		mcnt = lp->matcnt[col];
		mbeg = lp->matbeg[col];
		for (i = 0; i < mcnt; i++)
			EGLPNUM_TYPENAME_EGlpNumAddInnProdTo (sum, v[lp->matind[mbeg + i]], lp->matval[mbeg + i]);

		if (EGLPNUM_TYPENAME_EGlpNumIsNeqZero (sum, ztoler))
		{
			EGLPNUM_TYPENAME_EGlpNumCopy (zA->coef[nz], sum);
			zA->indx[nz] = j;
			nz++;
		}
	}
	zA->nzcnt = nz;

	EGLPNUM_TYPENAME_EGlpNumClearVar (sum);
	EGLPNUM_TYPENAME_EGlpNumFreeArray (v);
	EG_RETURN (rval);
}


static int compute_zA3 (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_svector * z,
	EGLPNUM_TYPENAME_svector * zA,
	EGLPNUM_TYPE ztoler)
{
	int rval = 0;
	int i, j, k, ix;
	int nz = 0;
	int row, col;
	int rcnt, rbeg;
	EGLPNUM_TYPE val;

	EGLPNUM_TYPENAME_EGlpNumInitVar (val);
	k = 0;
	for (i = 0; i < z->nzcnt; i++)
	{
		row = z->indx[i];
		EGLPNUM_TYPENAME_EGlpNumCopy (val, z->coef[i]);
		rcnt = lp->rowcnt[row];
		rbeg = lp->rowbeg[row];
		for (j = 0; j < rcnt; j++)
		{
			col = lp->rowind[rbeg + j];
			if (lp->vstat[col] != STAT_BASIC)
			{
				ix = lp->vindex[col];
				if (lp->iwork[ix] == 0)
				{
					lp->iwork[ix] = 1;
					lp->work.indx[k++] = ix;
				}
				EGLPNUM_TYPENAME_EGlpNumAddInnProdTo (lp->work.coef[ix], val, lp->rowval[rbeg + j]);
			}
		}
	}
	for (j = 0; j < k; j++)
	{
		ix = lp->work.indx[j];
		EGLPNUM_TYPENAME_EGlpNumCopy (val, lp->work.coef[ix]);
		EGLPNUM_TYPENAME_EGlpNumZero (lp->work.coef[ix]);
		lp->iwork[ix] = 0;
		if (EGLPNUM_TYPENAME_EGlpNumIsNeqZero (val, ztoler))
		{
			EGLPNUM_TYPENAME_EGlpNumCopy (zA->coef[nz], val);
			zA->indx[nz] = ix;
			nz++;
		}
	}
	zA->nzcnt = nz;
	EGLPNUM_TYPENAME_EGlpNumClearVar (val);
	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLfct_compute_zA (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_svector * z,
	EGLPNUM_TYPENAME_svector * zA)
{
	if (z->nzcnt < lp->nrows / 2)
		return compute_zA3 (lp, z, zA, EGLPNUM_TYPENAME_PIVZ_TOLER);
	else
		return compute_zA1 (lp, z, zA, EGLPNUM_TYPENAME_PIVZ_TOLER);
}

/* compute v^T A */
void EGLPNUM_TYPENAME_ILLfct_compute_vA (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_svector * v,
	EGLPNUM_TYPE * vA)
{
	int i, j;
	int row, col;
	int rcnt, rbeg;
	EGLPNUM_TYPE val;

	EGLPNUM_TYPENAME_EGlpNumInitVar (val);

	for (j = 0; j < lp->ncols; j++)
		EGLPNUM_TYPENAME_EGlpNumZero (vA[j]);

	for (i = 0; i < v->nzcnt; i++)
	{
		row = v->indx[i];
		EGLPNUM_TYPENAME_EGlpNumCopy (val, v->coef[i]);
		rcnt = lp->rowcnt[row];
		rbeg = lp->rowbeg[row];
		for (j = 0; j < rcnt; j++)
		{
			col = lp->rowind[rbeg + j];
			EGLPNUM_TYPENAME_EGlpNumAddInnProdTo (vA[col], val, lp->rowval[rbeg + j]);
		}
	}

	for (j = 0; j < lp->ncols; j++)
		if (!EGLPNUM_TYPENAME_EGlpNumIsNeqZero (vA[j], EGLPNUM_TYPENAME_SZERO_TOLER))
			EGLPNUM_TYPENAME_EGlpNumZero (vA[j]);

	EGLPNUM_TYPENAME_EGlpNumClearVar (val);
	return;
}

/* update information */

/*
1) lvstat - new status of leaving var.
*/
void EGLPNUM_TYPENAME_ILLfct_update_basis_info (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int eindex,
	int lindex,
	int lvstat)
{
	int evar;
	int lvar;

	evar = lp->nbaz[eindex];

	if (lindex >= 0)
	{															/* variable leaves basis */
		lvar = lp->baz[lindex];
		lp->vstat[evar] = STAT_BASIC;
		lp->vstat[lvar] = lvstat;
		lp->vindex[evar] = lindex;
		lp->vindex[lvar] = eindex;
		lp->baz[lindex] = evar;
		lp->nbaz[eindex] = lvar;
		(lp->basisid)++;
	}
	else
	{
		lp->vstat[evar] = (lp->vstat[evar] == STAT_LOWER) ? STAT_UPPER : STAT_LOWER;
	}
}

void EGLPNUM_TYPENAME_ILLfct_update_xz (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPE tz,
	int eindex,
	int lindex)
{
	int i, evar, estat;

	ILL_IFTRACE ("%s:%la:%d:%d:%d\n", __func__, EGLPNUM_TYPENAME_EGlpNumToLf (tz), eindex,
							 lindex, lp->yjz.nzcnt);

	if (EGLPNUM_TYPENAME_EGlpNumIsNeqqZero (tz))
		for (i = 0; i < lp->yjz.nzcnt; i++)
			EGLPNUM_TYPENAME_EGlpNumSubInnProdTo (lp->xbz[lp->yjz.indx[i]], tz, lp->yjz.coef[i]);

	if (lindex >= 0)
	{															/* variable leaves basis */
		evar = lp->nbaz[eindex];
		estat = lp->vstat[evar];
		if (estat == STAT_LOWER)
			EGLPNUM_TYPENAME_EGlpNumCopySum (lp->xbz[lindex], lp->lz[evar], tz);
		else if (estat == STAT_UPPER)
			EGLPNUM_TYPENAME_EGlpNumCopySum (lp->xbz[lindex], lp->uz[evar], tz);
		else if (estat == STAT_ZERO)
			EGLPNUM_TYPENAME_EGlpNumCopy (lp->xbz[lindex], tz);
	}
}

void EGLPNUM_TYPENAME_ILLfct_update_piz (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPE alpha)
{
	int i;

	for (i = 0; i < lp->zz.nzcnt; i++)
		EGLPNUM_TYPENAME_EGlpNumAddInnProdTo (lp->piz[lp->zz.indx[i]], alpha, lp->zz.coef[i]);
}

void EGLPNUM_TYPENAME_ILLfct_update_pIpiz (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_svector * z,
	const EGLPNUM_TYPE alpha)
{
	int i;

	if (!EGLPNUM_TYPENAME_EGlpNumIsNeqqZero (alpha))
		return;
	if (EGLPNUM_TYPENAME_EGlpNumIsEqqual (alpha, EGLPNUM_TYPENAME_oneLpNum))
	{
		for (i = 0; i < z->nzcnt; i++)
			EGLPNUM_TYPENAME_EGlpNumAddTo (lp->pIpiz[z->indx[i]], z->coef[i]);
	}
	else
	{
		for (i = 0; i < z->nzcnt; i++)
			EGLPNUM_TYPENAME_EGlpNumAddInnProdTo (lp->pIpiz[z->indx[i]], alpha, z->coef[i]);
	}
}

void EGLPNUM_TYPENAME_ILLfct_update_dz (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int eindex,
	EGLPNUM_TYPE alpha)
{
	int i;

	for (i = 0; i < lp->zA.nzcnt; i++)
		EGLPNUM_TYPENAME_EGlpNumSubInnProdTo (lp->dz[lp->zA.indx[i]], alpha, lp->zA.coef[i]);
	EGLPNUM_TYPENAME_EGlpNumCopyNeg (lp->dz[eindex], alpha);
}

void EGLPNUM_TYPENAME_ILLfct_update_pIdz (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_svector * zA,
	int eindex,
	const EGLPNUM_TYPE alpha)
{
	int i;

	if (!EGLPNUM_TYPENAME_EGlpNumIsNeqqZero (alpha))
		return;

	if (EGLPNUM_TYPENAME_EGlpNumIsEqqual (alpha, EGLPNUM_TYPENAME_oneLpNum))
	{
		for (i = 0; i < zA->nzcnt; i++)
			EGLPNUM_TYPENAME_EGlpNumSubTo (lp->pIdz[zA->indx[i]], zA->coef[i]);
	}
	else
	{
		for (i = 0; i < zA->nzcnt; i++)
			EGLPNUM_TYPENAME_EGlpNumSubInnProdTo (lp->pIdz[zA->indx[i]], alpha, zA->coef[i]);
	}
	if (eindex > -1)
		EGLPNUM_TYPENAME_EGlpNumCopyNeg (lp->pIdz[eindex], alpha);
}

/* bound and coef shift routines */

/* scale bound in my_rand to get more random digits, unless bound is large */
static double my_rand (
	int bound,
	ILLrandstate * r)
{
	int k = bound, scale = 1;
	double v = 0.0;

	if (bound < 100000)
	{
		k = 20000 * bound;
		scale = 20000;
	}
	v = 1 + (ILLutil_lprand (r) % (k));
	return v / (double) scale;
}

static int expand_var_bounds (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPE ftol,
	int *chgb)
{
	int rval = 0;
	int i, col, nchg = 0;
	EGLPNUM_TYPE newb, cftol;
	EGLPNUM_TYPE *x, *l, *u;
	ILLrandstate r;

	EGLPNUM_TYPENAME_EGlpNumInitVar (newb);
	EGLPNUM_TYPENAME_EGlpNumInitVar (cftol);
	EGLPNUM_TYPENAME_EGlpNumCopyAbs (cftol, ftol);
	EGLPNUM_TYPENAME_EGlpNumDivUiTo (cftol, 10);

	ILLutil_sprand (1, &r);

	for (i = 0; i < lp->nrows; i++)
	{
		col = lp->baz[i];
		if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFREE)
			continue;
		x = &(lp->xbz[i]);
		l = &(lp->lz[col]);
		u = &(lp->uz[col]);
		/* we use newb as temporal variable outside the if's scope */
		EGLPNUM_TYPENAME_EGlpNumCopyDiff (newb, *x, ftol);
		if (EGLPNUM_TYPENAME_EGlpNumIsNeqq (*l, EGLPNUM_TYPENAME_NINFTY) && EGLPNUM_TYPENAME_EGlpNumIsLess (newb, *l))
		{
			EGLPNUM_TYPENAME_EGlpNumSet (newb, -1.0 * (my_rand (50, &(lp->rstate)) + 1.0));
			EGLPNUM_TYPENAME_EGlpNumMultTo (newb, cftol);
			if (EGLPNUM_TYPENAME_EGlpNumIsLess (*x, *l))
				EGLPNUM_TYPENAME_EGlpNumAddTo (newb, *x);
			else
				EGLPNUM_TYPENAME_EGlpNumAddTo (newb, *l);
			rval = EGLPNUM_TYPENAME_ILLfct_bound_shift (lp, col, BOUND_LOWER, newb);
			CHECKRVALG (rval, CLEANUP);
			nchg++;
		}
		EGLPNUM_TYPENAME_EGlpNumCopySum (newb, *x, ftol);
		if (EGLPNUM_TYPENAME_EGlpNumIsNeqq (*u, EGLPNUM_TYPENAME_INFTY) && EGLPNUM_TYPENAME_EGlpNumIsLess (*u, newb))
		{
			EGLPNUM_TYPENAME_EGlpNumSet (newb, my_rand (50, &(lp->rstate)) + 1.0);
			EGLPNUM_TYPENAME_EGlpNumMultTo (newb, cftol);
			if (EGLPNUM_TYPENAME_EGlpNumIsLess (*x, *u))
				EGLPNUM_TYPENAME_EGlpNumAddTo (newb, *u);
			else
				EGLPNUM_TYPENAME_EGlpNumAddTo (newb, *x);
			rval = EGLPNUM_TYPENAME_ILLfct_bound_shift (lp, col, BOUND_UPPER, newb);
			CHECKRVALG (rval, CLEANUP);
			nchg++;
		}
	}
	*chgb = nchg;

CLEANUP:
	EGLPNUM_TYPENAME_EGlpNumClearVar (newb);
	EGLPNUM_TYPENAME_EGlpNumClearVar (cftol);
	EG_RETURN (rval);
}

static int expand_phaseI_bounds (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int *chgb)
{
	int rval = 0;
	int i, col, nchg = 0;
	EGLPNUM_TYPE newb, cftol;
	EGLPNUM_TYPE *u, *l, *x;
	ILLrandstate r;

	EGLPNUM_TYPENAME_EGlpNumInitVar (newb);
	EGLPNUM_TYPENAME_EGlpNumInitVar (cftol);
	EGLPNUM_TYPENAME_EGlpNumCopyAbs (cftol, lp->tol->ip_tol);
	EGLPNUM_TYPENAME_EGlpNumDivUiTo (cftol, 10);
	ILLutil_sprand (1, &r);

	for (i = 0; i < lp->nrows; i++)
	{
		col = lp->baz[i];
		if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFREE)
			continue;
		x = &(lp->xbz[i]);
		l = &(lp->lz[col]);
		u = &(lp->uz[col]);

		if (EGLPNUM_TYPENAME_EGlpNumIsNeqq (*l, EGLPNUM_TYPENAME_NINFTY) && EGLPNUM_TYPENAME_EGlpNumIsEqual (*x, *l, cftol))
		{
			EGLPNUM_TYPENAME_EGlpNumSet (newb, my_rand (50, &(lp->rstate)) + 1.0);
			EGLPNUM_TYPENAME_EGlpNumMultTo (newb, cftol);
			EGLPNUM_TYPENAME_EGlpNumSign (newb);
			EGLPNUM_TYPENAME_EGlpNumAddTo (newb, *l);
			rval = EGLPNUM_TYPENAME_ILLfct_bound_shift (lp, col, BOUND_LOWER, newb);
			CHECKRVALG (rval, CLEANUP);
			nchg++;
		}
		if (EGLPNUM_TYPENAME_EGlpNumIsNeqq (*u, EGLPNUM_TYPENAME_INFTY) && EGLPNUM_TYPENAME_EGlpNumIsEqual (*x, *u, cftol))
		{
			EGLPNUM_TYPENAME_EGlpNumSet (newb, my_rand (50, &(lp->rstate)) + 1.0);
			EGLPNUM_TYPENAME_EGlpNumMultTo (newb, cftol);
			EGLPNUM_TYPENAME_EGlpNumAddTo (newb, *u);
			rval = EGLPNUM_TYPENAME_ILLfct_bound_shift (lp, col, BOUND_UPPER, newb);
			CHECKRVALG (rval, CLEANUP);
			nchg++;
		}
	}
	*chgb = nchg;

CLEANUP:
	EGLPNUM_TYPENAME_EGlpNumClearVar (newb);
	EGLPNUM_TYPENAME_EGlpNumClearVar (cftol);
	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLfct_adjust_viol_bounds (
	EGLPNUM_TYPENAME_lpinfo * lp)
{
	int rval = 0;
	int chgb = 0;
	EGLPNUM_TYPE tol;

	EGLPNUM_TYPENAME_EGlpNumInitVar (tol);
	EGLPNUM_TYPENAME_EGlpNumCopyNeg (tol, lp->tol->pfeas_tol);
	rval = expand_var_bounds (lp, tol, &chgb);
#if FCT_DEBUG > 0
	if (rval == 0)
		QSlog("adjusting %d bounds", chgb);
#endif
	EGLPNUM_TYPENAME_EGlpNumClearVar (tol);
	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLfct_perturb_bounds (
	EGLPNUM_TYPENAME_lpinfo * lp)
{
	int rval = 0;
	int chgb = 0;

	rval = expand_var_bounds (lp, lp->tol->ip_tol, &chgb);
#if FCT_DEBUG > 0
	if (rval == 0)
		QSlog("perturbing %d bounds", chgb);
#endif
	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLfct_perturb_phaseI_bounds (
	EGLPNUM_TYPENAME_lpinfo * lp)
{
	int rval = 0;
	int chgb = 0;

	rval = expand_phaseI_bounds (lp, &chgb);
#if FCT_DEBUG > 0
	if (rval == 0)
		QSlog("perturbing %d phase I bounds", chgb);
#endif
	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLfct_bound_shift (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int col,
	int bndtype,
	EGLPNUM_TYPE newbnd)
{
	int rval = 0;
	EGLPNUM_TYPENAME_bndinfo *nbnd = 0;

	ILL_IFTRACE ("\n%s:%d:%d:%la", __func__, col, bndtype, EGLPNUM_TYPENAME_EGlpNumToLf (newbnd));
	nbnd = EGLPNUM_TYPENAME_ILLfct_new_bndinfo ();

	nbnd->varnum = col;
	nbnd->btype = bndtype;
	if (bndtype == BOUND_LOWER)
	{
		EGLPNUM_TYPENAME_EGlpNumCopy (nbnd->pbound, lp->lz[col]);
		EGLPNUM_TYPENAME_EGlpNumCopy (nbnd->cbound, newbnd);
		EGLPNUM_TYPENAME_EGlpNumCopy (lp->lz[col], newbnd);
	}
	else
	{
		EGLPNUM_TYPENAME_EGlpNumCopy (nbnd->pbound, lp->uz[col]);
		EGLPNUM_TYPENAME_EGlpNumCopy (nbnd->cbound, newbnd);
		EGLPNUM_TYPENAME_EGlpNumCopy (lp->uz[col], newbnd);
	}
	ILL_IFTRACE (":%la", EGLPNUM_TYPENAME_EGlpNumToLf (nbnd->pbound));
	if (lp->vtype[col] == VFIXED || lp->vtype[col] == VARTIFICIAL)
	{
		/* QSlog("changing f/a bound"); */
		if (EGLPNUM_TYPENAME_EGlpNumIsLess (lp->lz[col], lp->uz[col]))
			lp->vtype[col] = VBOUNDED;
	}

	nbnd->next = lp->bchanges;
	lp->bchanges = nbnd;
	lp->nbchange++;

//CLEANUP:
	if (rval)
		EGLPNUM_TYPENAME_ILLfct_free_bndinfo (nbnd);
	ILL_IFTRACE ("\n");
	EG_RETURN (rval);
}

void EGLPNUM_TYPENAME_ILLfct_unroll_bound_change (
	EGLPNUM_TYPENAME_lpinfo * lp)
{
	int col;
	int changex = 0;
	EGLPNUM_TYPENAME_bndinfo *bptr = lp->bchanges;
	EGLPNUM_TYPENAME_bndinfo *nptr = 0;

	ILL_IFTRACE ("%s:", __func__);

	while (lp->nbchange != 0)
	{
		col = bptr->varnum;
		ILL_IFTRACE (":%d", col);

		if (bptr->btype == BOUND_UPPER)
			EGLPNUM_TYPENAME_EGlpNumCopy (lp->uz[col], bptr->pbound);
		else
			EGLPNUM_TYPENAME_EGlpNumCopy (lp->lz[col], bptr->pbound);

		if (lp->vtype[col] == VBOUNDED)
		{
			if (EGLPNUM_TYPENAME_EGlpNumIsEqqual (lp->lz[col], lp->uz[col]))
				lp->vtype[col] = (!EGLPNUM_TYPENAME_EGlpNumIsNeqqZero (lp->lz[col])) ?
					VARTIFICIAL : VFIXED;
		}

		if (lp->vstat[col] != STAT_BASIC)
		{
			if ((bptr->btype == BOUND_UPPER && lp->vstat[col] == STAT_UPPER) ||
					(bptr->btype == BOUND_LOWER && lp->vstat[col] == STAT_LOWER))
				changex++;
		}
		nptr = bptr->next;
		EGLPNUM_TYPENAME_EGlpNumClearVar ((bptr->cbound));
		EGLPNUM_TYPENAME_EGlpNumClearVar ((bptr->pbound));
		ILL_IFFREE (bptr, EGLPNUM_TYPENAME_bndinfo);
		bptr = nptr;
		lp->nbchange--;
	}
	lp->bchanges = bptr;
	ILL_IFTRACE ("\n");
	if (changex)
		EGLPNUM_TYPENAME_ILLfct_compute_xbz (lp);
}

static int expand_var_coefs (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPE ftol,
	int *chgc)
{
	int rval = 0;
	int i, col, vs, vt;
	int nchg = 0;
	EGLPNUM_TYPE newc, cftol, mftol[1];
	EGLPNUM_TYPE *c, *dj;
	ILLrandstate r;

	EGLPNUM_TYPENAME_EGlpNumInitVar (newc);
	EGLPNUM_TYPENAME_EGlpNumInitVar (cftol);
	EGLPNUM_TYPENAME_EGlpNumInitVar (mftol[0]);
	EGLPNUM_TYPENAME_EGlpNumCopyAbs (cftol, ftol);
	EGLPNUM_TYPENAME_EGlpNumDivUiTo (cftol, 10);
	EGLPNUM_TYPENAME_EGlpNumCopyNeg (mftol[0], ftol);
	ILLutil_sprand (1, &r);

	for (i = 0; i < lp->nnbasic; i++)
	{
		dj = &(lp->dz[i]);
		col = lp->nbaz[i];
		c = &(lp->cz[col]);
		vs = lp->vstat[col];
		vt = lp->vtype[col];

		if (vt == VARTIFICIAL || vt == VFIXED)
			continue;
		switch (vs)
		{
		case STAT_ZERO:
			EGLPNUM_TYPENAME_EGlpNumCopyDiff (newc, *c, *dj);
			rval = EGLPNUM_TYPENAME_ILLfct_coef_shift (lp, col, newc);
			CHECKRVALG (rval, CLEANUP);
			nchg++;
			break;
		case STAT_LOWER:
			if (EGLPNUM_TYPENAME_EGlpNumIsLess (*dj, ftol))
			{
				EGLPNUM_TYPENAME_EGlpNumSet (newc, my_rand (50, &(lp->rstate)) + 1.0);
				EGLPNUM_TYPENAME_EGlpNumMultTo (newc, cftol);
				EGLPNUM_TYPENAME_EGlpNumAddTo (newc, *c);
				if (EGLPNUM_TYPENAME_EGlpNumIsLessZero (*dj))
					EGLPNUM_TYPENAME_EGlpNumSubTo (newc, *dj);
				rval = EGLPNUM_TYPENAME_ILLfct_coef_shift (lp, col, newc);
				CHECKRVALG (rval, CLEANUP);
				nchg++;
			}
			break;
		case STAT_UPPER:
			if (EGLPNUM_TYPENAME_EGlpNumIsLess (mftol[0], *dj))
			{
				EGLPNUM_TYPENAME_EGlpNumSet (newc, my_rand (50, &(lp->rstate)) + 1.0);
				EGLPNUM_TYPENAME_EGlpNumMultTo (newc, cftol);
				EGLPNUM_TYPENAME_EGlpNumSign (newc);
				EGLPNUM_TYPENAME_EGlpNumAddTo (newc, *c);
				if (EGLPNUM_TYPENAME_EGlpNumIsGreatZero (*dj))
					EGLPNUM_TYPENAME_EGlpNumSubTo (newc, *dj);
				rval = EGLPNUM_TYPENAME_ILLfct_coef_shift (lp, col, newc);
				CHECKRVALG (rval, CLEANUP);
				nchg++;
			}
			break;
		default:
			break;
		}
	}
	*chgc = nchg;

CLEANUP:
	EGLPNUM_TYPENAME_EGlpNumClearVar (mftol[0]);
	EGLPNUM_TYPENAME_EGlpNumClearVar (newc);
	EGLPNUM_TYPENAME_EGlpNumClearVar (cftol);
	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLfct_adjust_viol_coefs (
	EGLPNUM_TYPENAME_lpinfo * lp)
{
	int rval = 0;
	int chgc = 0;
	EGLPNUM_TYPE tol;

	EGLPNUM_TYPENAME_EGlpNumInitVar (tol);
	EGLPNUM_TYPENAME_EGlpNumCopyNeg (tol, lp->tol->dfeas_tol);

	rval = expand_var_coefs (lp, tol, &chgc);
#if FCT_DEBUG > 0
	if (rval == 0)
		QSlog("perturbing %d coefs", chgc);
#endif
	EGLPNUM_TYPENAME_EGlpNumClearVar (tol);
	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLfct_perturb_coefs (
	EGLPNUM_TYPENAME_lpinfo * lp)
{
	int rval = 0;
	int chgc = 0;

	rval = expand_var_coefs (lp, lp->tol->id_tol, &chgc);
#if FCT_DEBUG > 0
	if (rval == 0)
		QSlog("perturbing %d coefs", chgc);
#endif
	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLfct_coef_shift (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int col,
	EGLPNUM_TYPE newcoef)
{
	int rval = 0;
	EGLPNUM_TYPENAME_coefinfo *ncoef = 0;

	ILL_SAFE_MALLOC (ncoef, 1, EGLPNUM_TYPENAME_coefinfo);
	EGLPNUM_TYPENAME_EGlpNumInitVar ((ncoef->pcoef));
	EGLPNUM_TYPENAME_EGlpNumInitVar ((ncoef->ccoef));

	ncoef->varnum = col;
	EGLPNUM_TYPENAME_EGlpNumCopy (ncoef->pcoef, lp->cz[col]);
	EGLPNUM_TYPENAME_EGlpNumCopy (ncoef->ccoef, newcoef);
	EGLPNUM_TYPENAME_EGlpNumCopy (lp->cz[col], newcoef);
	ncoef->next = lp->cchanges;
	lp->cchanges = ncoef;
	EGLPNUM_TYPENAME_EGlpNumAddTo (lp->dz[lp->vindex[col]], ncoef->ccoef);
	EGLPNUM_TYPENAME_EGlpNumSubTo (lp->dz[lp->vindex[col]], ncoef->pcoef);
	lp->ncchange++;

CLEANUP:
	if (rval)
	{
		EGLPNUM_TYPENAME_EGlpNumClearVar ((ncoef->pcoef));
		EGLPNUM_TYPENAME_EGlpNumClearVar ((ncoef->ccoef));
		ILL_IFFREE (ncoef, EGLPNUM_TYPENAME_coefinfo);
	}
	EG_RETURN (rval);
}

void EGLPNUM_TYPENAME_ILLfct_unroll_coef_change (
	EGLPNUM_TYPENAME_lpinfo * lp)
{
	int bascoef = 0;
	EGLPNUM_TYPENAME_coefinfo *cptr = (EGLPNUM_TYPENAME_coefinfo *) lp->cchanges;
	EGLPNUM_TYPENAME_coefinfo *nptr = 0;

	while (lp->ncchange != 0)
	{
		EGLPNUM_TYPENAME_EGlpNumCopy (lp->cz[cptr->varnum], cptr->pcoef);
		if (lp->vstat[cptr->varnum] != STAT_BASIC)
		{
			EGLPNUM_TYPENAME_EGlpNumAddTo (lp->dz[lp->vindex[cptr->varnum]], cptr->pcoef);
			EGLPNUM_TYPENAME_EGlpNumSubTo (lp->dz[lp->vindex[cptr->varnum]], cptr->ccoef);
		}
		else
			bascoef++;

		nptr = cptr->next;
		EGLPNUM_TYPENAME_EGlpNumClearVar ((cptr->pcoef));
		EGLPNUM_TYPENAME_EGlpNumClearVar ((cptr->ccoef));
		ILL_IFFREE (cptr, EGLPNUM_TYPENAME_coefinfo);
		cptr = nptr;
		lp->ncchange--;
	}
	lp->cchanges = cptr;
	if (bascoef)
	{
		EGLPNUM_TYPENAME_ILLfct_compute_piz (lp);
		EGLPNUM_TYPENAME_ILLfct_compute_dz (lp);
	}
}

/* feasibility routines */
void EGLPNUM_TYPENAME_ILLfct_check_pfeasible (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_feas_info * fs,
	const EGLPNUM_TYPE ftol)
{
	int i, col;
	EGLPNUM_TYPE infeas, err1, err2;

	EGLPNUM_TYPENAME_EGlpNumInitVar (infeas);
	EGLPNUM_TYPENAME_EGlpNumInitVar (err1);
	EGLPNUM_TYPENAME_EGlpNumInitVar (err2);
	EGLPNUM_TYPENAME_EGlpNumZero (infeas);
	fs->pstatus = PRIMAL_FEASIBLE;
	EGLPNUM_TYPENAME_EGlpNumZero (fs->totinfeas);
	ILL_IFTRACE ("%s:tol %la\n", __func__, EGLPNUM_TYPENAME_EGlpNumToLf (ftol));

	for (i = 0; i < lp->nrows; i++)
	{
		col = lp->baz[i];
		EGLPNUM_TYPENAME_EGlpNumCopyDiff (err1, lp->xbz[i], lp->uz[col]);
		EGLPNUM_TYPENAME_EGlpNumCopyDiff (err2, lp->lz[col], lp->xbz[i]);
		if (EGLPNUM_TYPENAME_EGlpNumIsLess (ftol, err1)
				&& EGLPNUM_TYPENAME_EGlpNumIsNeq (lp->uz[col], EGLPNUM_TYPENAME_INFTY, EGLPNUM_TYPENAME_oneLpNum))
		{
			EGLPNUM_TYPENAME_EGlpNumAddTo (infeas, err1);
			WARNINGL (QSE_WLVL, EGLPNUM_TYPENAME_EGlpNumIsLess (EGLPNUM_TYPENAME_INFTY, err1),
							 "This is imposible lu = %15lg xbz = %15lg" " EGLPNUM_TYPENAME_INFTY = %15lg",
							 EGLPNUM_TYPENAME_EGlpNumToLf (lp->uz[col]), EGLPNUM_TYPENAME_EGlpNumToLf (lp->xbz[i]),
							 EGLPNUM_TYPENAME_EGlpNumToLf (EGLPNUM_TYPENAME_INFTY));
			lp->bfeas[i] = 1;
		}
		else if (EGLPNUM_TYPENAME_EGlpNumIsLess (ftol, err2)
						 && EGLPNUM_TYPENAME_EGlpNumIsNeq (lp->lz[col], EGLPNUM_TYPENAME_NINFTY, EGLPNUM_TYPENAME_oneLpNum))
		{
			EGLPNUM_TYPENAME_EGlpNumAddTo (infeas, err2);
			WARNINGL (QSE_WLVL, EGLPNUM_TYPENAME_EGlpNumIsLess (EGLPNUM_TYPENAME_INFTY, err2),
							 "This is imposible lz = %15lg xbz = %15lg" " EGLPNUM_TYPENAME_NINFTY = %15lg",
							 EGLPNUM_TYPENAME_EGlpNumToLf (lp->lz[col]), EGLPNUM_TYPENAME_EGlpNumToLf (lp->xbz[i]),
							 EGLPNUM_TYPENAME_EGlpNumToLf (EGLPNUM_TYPENAME_NINFTY));
			lp->bfeas[i] = -1;
		}
		else
			lp->bfeas[i] = 0;
	}
	if (EGLPNUM_TYPENAME_EGlpNumIsNeqqZero (infeas))
	{
		fs->pstatus = PRIMAL_INFEASIBLE;
		EGLPNUM_TYPENAME_EGlpNumCopy (fs->totinfeas, infeas);
		ILL_IFTRACE ("%s:inf %la\n", __func__, EGLPNUM_TYPENAME_EGlpNumToLf (infeas));
		if (EGLPNUM_TYPENAME_EGlpNumIsLessZero (fs->totinfeas))
		{
			QSlog("Negative infeasibility, Imposible! %lf %la",
									EGLPNUM_TYPENAME_EGlpNumToLf (infeas), EGLPNUM_TYPENAME_EGlpNumToLf (infeas));
		}
	}
	EGLPNUM_TYPENAME_EGlpNumCopy (lp->pinfeas, infeas);
	EGLPNUM_TYPENAME_EGlpNumClearVar (infeas);
	EGLPNUM_TYPENAME_EGlpNumClearVar (err1);
	EGLPNUM_TYPENAME_EGlpNumClearVar (err2);
}

/* feasibility routines */
void EGLPNUM_TYPENAME_ILLfct_check_pIpfeasible (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_feas_info * fs,
	EGLPNUM_TYPE ftol)
{
	int i, col;
	int ninf = 0;

	fs->pstatus = PRIMAL_FEASIBLE;
	EGLPNUM_TYPENAME_EGlpNumZero (fs->totinfeas);

	for (i = 0; i < lp->nrows; i++)
	{
		if (!EGLPNUM_TYPENAME_EGlpNumIsNeqZero (lp->xbz[i], ftol))
			continue;
		col = lp->baz[i];
		if (EGLPNUM_TYPENAME_EGlpNumIsGreatZero(lp->xbz[i]) &&
				EGLPNUM_TYPENAME_EGlpNumIsNeqq (lp->uz[col], EGLPNUM_TYPENAME_INFTY))
		{
			ninf++;
		}
		else if (EGLPNUM_TYPENAME_EGlpNumIsLessZero (lp->xbz[i]) &&
						 EGLPNUM_TYPENAME_EGlpNumIsNeqq (lp->lz[col], EGLPNUM_TYPENAME_NINFTY))
		{
			ninf++;
		}
	}
	if (ninf != 0)
		fs->pstatus = PRIMAL_INFEASIBLE;
}

void EGLPNUM_TYPENAME_ILLfct_check_dfeasible (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_feas_info * fs,
	const EGLPNUM_TYPE ftol)
{
	int j, col;
	EGLPNUM_TYPE infeas;

	EGLPNUM_TYPENAME_EGlpNumInitVar (infeas);
	EGLPNUM_TYPENAME_EGlpNumZero (infeas);
	fs->dstatus = DUAL_FEASIBLE;
	EGLPNUM_TYPENAME_EGlpNumZero (fs->totinfeas);

	for (j = 0; j < lp->nnbasic; j++)
	{
		lp->dfeas[j] = 0;
		if (!EGLPNUM_TYPENAME_EGlpNumIsNeqZero (lp->dz[j], ftol))
			continue;
		col = lp->nbaz[j];

		if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED)
			continue;

		if (EGLPNUM_TYPENAME_EGlpNumIsLessZero (lp->dz[j]) &&
				(lp->vstat[col] == STAT_LOWER || lp->vstat[col] == STAT_ZERO))
		{
			EGLPNUM_TYPENAME_EGlpNumSubTo (infeas, lp->dz[j]);
			lp->dfeas[j] = -1;
		}
		else if (EGLPNUM_TYPENAME_EGlpNumIsGreatZero (lp->dz[j]) &&
						 (lp->vstat[col] == STAT_UPPER || lp->vstat[col] == STAT_ZERO))
		{
			EGLPNUM_TYPENAME_EGlpNumAddTo (infeas, lp->dz[j]);
			lp->dfeas[j] = 1;
		}
	}

	if (EGLPNUM_TYPENAME_EGlpNumIsNeqqZero (infeas))
	{
		EGLPNUM_TYPENAME_EGlpNumCopy (fs->totinfeas, infeas);
		fs->dstatus = DUAL_INFEASIBLE;
		ILL_IFTRACE ("%s:inf %la\n", __func__, EGLPNUM_TYPENAME_EGlpNumToLf (infeas));
		if (EGLPNUM_TYPENAME_EGlpNumIsLessZero (fs->totinfeas))
		{
			QSlog("Negative infeasibility, Imposible! %lf %la",
									EGLPNUM_TYPENAME_EGlpNumToLf (infeas), EGLPNUM_TYPENAME_EGlpNumToLf (infeas));
		}
	}
	EGLPNUM_TYPENAME_EGlpNumCopy (lp->dinfeas, infeas);
	EGLPNUM_TYPENAME_EGlpNumClearVar (infeas);
}

void EGLPNUM_TYPENAME_ILLfct_check_pIdfeasible (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_feas_info * fs,
	EGLPNUM_TYPE ftol)
{
	int j, col;
	int ninf = 0;
	EGLPNUM_TYPE *dz = lp->pIdz;

	fs->dstatus = DUAL_FEASIBLE;

	for (j = 0; j < lp->nnbasic; j++)
	{
		if (!EGLPNUM_TYPENAME_EGlpNumIsNeqZero (dz[j], ftol))
			continue;
		col = lp->nbaz[j];
		if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED)
			continue;

		if (EGLPNUM_TYPENAME_EGlpNumIsLessZero (dz[j]) &&
				(lp->vstat[col] == STAT_LOWER || lp->vstat[col] == STAT_ZERO))
			ninf++;
		else if (EGLPNUM_TYPENAME_EGlpNumIsGreatZero (dz[j]) &&
						 (lp->vstat[col] == STAT_UPPER || lp->vstat[col] == STAT_ZERO))
			ninf++;
	}

	if (ninf != 0)
		fs->dstatus = DUAL_INFEASIBLE;
}

void EGLPNUM_TYPENAME_ILLfct_dual_adjust (
	EGLPNUM_TYPENAME_lpinfo * lp,
	const EGLPNUM_TYPE ftol)
{
	int j, col;

	for (j = 0; j < lp->nnbasic; j++)
	{
		if (!EGLPNUM_TYPENAME_EGlpNumIsNeqZero (lp->dz[j], ftol))
			continue;
		col = lp->nbaz[j];
		if (EGLPNUM_TYPENAME_EGlpNumIsLessZero (lp->dz[j]) &&
				EGLPNUM_TYPENAME_EGlpNumIsNeqq (lp->uz[col], EGLPNUM_TYPENAME_INFTY))
			lp->vstat[col] = STAT_UPPER;
		else if (EGLPNUM_TYPENAME_EGlpNumIsGreatZero (lp->dz[j]) &&
						 EGLPNUM_TYPENAME_EGlpNumIsNeqq (lp->lz[col], EGLPNUM_TYPENAME_NINFTY))
			lp->vstat[col] = STAT_LOWER;
	}
}

void EGLPNUM_TYPENAME_ILLfct_dphaseI_simple_update (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPE ftol)
{
	int j, col;

	for (j = 0; j < lp->nnbasic; j++)
	{
		if (!EGLPNUM_TYPENAME_EGlpNumIsNeqZero (lp->dz[j], ftol))
			continue;
		col = lp->nbaz[j];
		if (EGLPNUM_TYPENAME_EGlpNumIsLessZero (lp->dz[j] ) && lp->vtype[col] == VBOUNDED)
			lp->vstat[col] = STAT_UPPER;
		else if (EGLPNUM_TYPENAME_EGlpNumIsGreatZero (lp->dz[j]) && lp->vtype[col] == VBOUNDED)
			lp->vstat[col] = STAT_LOWER;
	}
}

/* set status values */
void EGLPNUM_TYPENAME_ILLfct_set_status_values (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int pstatus,
	int dstatus,
	int ptype,
	int dtype)
{
	if (dstatus == DUAL_FEASIBLE && dtype == PHASEII)
	{
		if (!lp->ncchange)
		{
			lp->probstat.dual_feasible = 1;
			lp->basisstat.dual_feasible = 1;
			lp->basisstat.dual_infeasible = 0;
		}
	}
	if (dstatus == DUAL_INFEASIBLE && dtype == PHASEII)
	{
		if (!lp->ncchange)
		{
			lp->basisstat.dual_feasible = 0;
			lp->basisstat.dual_infeasible = 1;
		}
		if (pstatus == PRIMAL_FEASIBLE && ptype == PHASEI)
			if (!lp->ncchange)
				lp->probstat.dual_infeasible = 1;
	}
	if (pstatus == PRIMAL_FEASIBLE && ptype == PHASEII)
	{
		if (!lp->nbchange)
		{
			lp->probstat.primal_feasible = 1;
			lp->basisstat.primal_feasible = 1;
			lp->basisstat.primal_infeasible = 0;
		}
	}
	if (pstatus == PRIMAL_INFEASIBLE && ptype == PHASEII)
	{
		lp->basisstat.primal_feasible = 0;
		lp->basisstat.primal_infeasible = 1;

		if (dstatus == DUAL_FEASIBLE && dtype == PHASEI)
			lp->probstat.primal_infeasible = 1;
	}
	if (pstatus == PRIMAL_UNBOUNDED)
	{
		if (!lp->nbchange)
		{
			lp->probstat.primal_unbounded = 1;
			lp->basisstat.primal_unbounded = 1;
			lp->probstat.dual_infeasible = 1;
			lp->basisstat.dual_infeasible = 1;
			lp->basisstat.dual_feasible = 0;
		}
	}
	if (dstatus == DUAL_UNBOUNDED)
	{
		if (!lp->ncchange)
		{
			lp->probstat.dual_unbounded = 1;
			lp->basisstat.dual_unbounded = 1;
			lp->probstat.primal_infeasible = 1;
			lp->basisstat.primal_infeasible = 1;
			lp->basisstat.primal_feasible = 0;
		}
	}
	if (lp->probstat.primal_feasible && lp->probstat.dual_feasible)
		lp->probstat.optimal = 1;

	if (lp->basisstat.primal_feasible && lp->basisstat.dual_feasible)
		lp->basisstat.optimal = 1;
	else
		lp->basisstat.optimal = 0;
}

void EGLPNUM_TYPENAME_ILLfct_init_counts (
	EGLPNUM_TYPENAME_lpinfo * lp)
{
	int i;
	EGLPNUM_TYPENAME_count_struct *c = lp->cnts;

#define C_VALUE(a) (1.0+(double)(a)/(PARAM_HEAP_RATIO*ILLutil_our_log2(a)))
	EGLPNUM_TYPENAME_EGlpNumSet (c->y_ravg, C_VALUE (lp->nrows));
	EGLPNUM_TYPENAME_EGlpNumSet (c->za_ravg, C_VALUE (lp->nnbasic));
	ILL_IFTRACE ("%s:%la\n", __func__, EGLPNUM_TYPENAME_EGlpNumToLf (c->za_ravg));
#undef C_VALUE
	c->ynz_cnt = 0;
	c->num_y = 0;
	c->znz_cnt = 0;
	c->num_z = 0;
	c->zanz_cnt = 0;
	c->num_za = 0;
	c->pnorm_cnt = 0;
	c->dnorm_cnt = 0;
	c->pinz_cnt = 0;
	c->num_pi = 0;
	c->pi1nz_cnt = 0;
	c->num_pi1 = 0;
	c->upnz_cnt = 0;
	c->num_up = 0;
	c->pupv_cnt = 0;
	c->dupv_cnt = 0;
	c->pI_iter = 0;
	c->pII_iter = 0;
	c->dI_iter = 0;
	c->dII_iter = 0;
	c->tot_iter = 0;
	for (i = 0; i < 10; i++)
	{
		c->pivpI[i] = 0;
		c->pivpII[i] = 0;
		c->pivdI[i] = 0;
		c->pivdII[i] = 0;
	}
}

static void update_piv_values (
	EGLPNUM_TYPENAME_count_struct * c,
	int phase,
	const EGLPNUM_TYPE piv2)
{
	int i = 0;
	EGLPNUM_TYPE v, piv;

	if (!EGLPNUM_TYPENAME_EGlpNumIsNeqqZero(piv2))
		return;
	EGLPNUM_TYPENAME_EGlpNumInitVar (v);
	EGLPNUM_TYPENAME_EGlpNumInitVar (piv);
	EGLPNUM_TYPENAME_EGlpNumCopyAbs (piv, piv2);
	EGLPNUM_TYPENAME_EGlpNumOne (v);
	while (EGLPNUM_TYPENAME_EGlpNumIsLess (piv, v) && i < 9)
	{
		EGLPNUM_TYPENAME_EGlpNumDivUiTo (v, 10);
		i++;
	}
	switch (phase)
	{
	case PRIMAL_PHASEI:
		c->pivpI[i]++;
		break;
	case PRIMAL_PHASEII:
		c->pivpII[i]++;
		break;
	case DUAL_PHASEI:
		c->pivdI[i]++;
		break;
	case DUAL_PHASEII:
		c->pivdII[i]++;
		break;
	default:
		break;
	}
	EGLPNUM_TYPENAME_EGlpNumClearVar (v);
	EGLPNUM_TYPENAME_EGlpNumClearVar (piv);
}

void EGLPNUM_TYPENAME_ILLfct_update_counts (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int f,
	int upi,
	const EGLPNUM_TYPE upd)
{
	EGLPNUM_TYPENAME_count_struct *c = lp->cnts;

	switch (f)
	{
	case CNT_PPHASE1ITER:
		c->pI_iter++;
		c->tot_iter++;
		break;
	case CNT_PPHASE2ITER:
		c->pII_iter++;
		c->tot_iter++;
		break;
	case CNT_DPHASE1ITER:
		c->dI_iter++;
		c->tot_iter++;
		break;
	case CNT_DPHASE2ITER:
		c->dII_iter++;
		c->tot_iter++;
		break;
	case CNT_YNZ:
		c->ynz_cnt += upi;
		c->num_y++;
		break;
	case CNT_ZANZ:
		c->zanz_cnt += upi;
		c->num_za++;
		break;
	case CNT_PINZ:
		c->pinz_cnt += upi;
		c->num_pi++;
		break;
	case CNT_P1PINZ:
		c->pi1nz_cnt += upi;
		c->num_pi1++;
		break;
	case CNT_UPNZ:
		c->upnz_cnt += upi;
		c->num_up++;
		break;
	case CNT_PIPIV:
		update_piv_values (c, PRIMAL_PHASEI, upd);
		break;
	case CNT_PIIPIV:
		update_piv_values (c, PRIMAL_PHASEII, upd);
		break;
	case CNT_DIPIV:
		update_piv_values (c, DUAL_PHASEI, upd);
		break;
	case CNT_DIIPIV:
		update_piv_values (c, DUAL_PHASEII, upd);
		break;
	case CNT_YRAVG:
		EGLPNUM_TYPENAME_EGlpNumMultUiTo (c->y_ravg, c->tot_iter);
		EGLPNUM_TYPENAME_EGlpNumAddUiTo (c->y_ravg, upi);
		EGLPNUM_TYPENAME_EGlpNumDivUiTo (c->y_ravg, c->tot_iter + 1);
		break;
	case CNT_ZARAVG:
		ILL_IFTRACE ("%s:%d:%d:%d:%la:%la", __func__, f, c->tot_iter, upi,
								 EGLPNUM_TYPENAME_EGlpNumToLf (upd), EGLPNUM_TYPENAME_EGlpNumToLf (c->za_ravg));
		EGLPNUM_TYPENAME_EGlpNumMultUiTo (c->za_ravg, c->tot_iter);
		EGLPNUM_TYPENAME_EGlpNumAddUiTo (c->za_ravg, upi);
		EGLPNUM_TYPENAME_EGlpNumDivUiTo (c->za_ravg, c->tot_iter + 1);
		ILL_IFTRACE (":%la\n", EGLPNUM_TYPENAME_EGlpNumToLf (c->za_ravg));
		break;
	}
}

void EGLPNUM_TYPENAME_ILLfct_print_counts (
	EGLPNUM_TYPENAME_lpinfo * lp)
{
	int i, niter;
	EGLPNUM_TYPENAME_count_struct *c = lp->cnts;

	c->tot_iter = c->pI_iter + c->pII_iter + c->dI_iter + c->dII_iter;
	niter = (c->tot_iter == 0) ? 1 : c->tot_iter;
	QSlog("Counts for problem %s", lp->O->probname);
	if (c->num_y != 0)
		QSlog("avg ynz = %.2f", (double) c->ynz_cnt / c->num_y);
	if (c->num_z != 0)
		QSlog("avg znz = %.2f", (double) c->znz_cnt / c->num_z);
	if (c->num_za != 0)
		QSlog("avg zanz = %.2f", (double) c->zanz_cnt / c->num_za);
	QSlog("avg pnorm = %.2f", (double) c->pnorm_cnt / lp->nnbasic);
	QSlog("avg dnorm = %.2f", (double) c->dnorm_cnt / lp->nrows);
	if (c->num_pi != 0)
		QSlog("avg pinz = %.2f", (double) c->pinz_cnt / c->num_pi);
	if (c->num_pi1 != 0)
		QSlog("avg piInz = %.2f", (double) c->pi1nz_cnt / c->num_pi1);
	if (c->num_up != 0)
		QSlog("avg upnz = %.2f", (double) c->upnz_cnt / c->num_up);

	for (i = 0; i < 10; i++)
		QSlog("piv 1.0e-%d : %d %d %d %d",
								i, c->pivpI[i], c->pivpII[i], c->pivdI[i], c->pivdII[i]);
}


/* c <- a + t*b */
static void add_vectors (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_svector * a,
	EGLPNUM_TYPENAME_svector * b,
	EGLPNUM_TYPENAME_svector * c,
	const EGLPNUM_TYPE t)
{
	int i, r, l;
	EGLPNUM_TYPENAME_svector *w = &(lp->work);

	for (i = 0; i < b->nzcnt; i++)
	{
		r = b->indx[i];
		w->indx[i] = r;
		EGLPNUM_TYPENAME_EGlpNumCopy (w->coef[r], t);
		EGLPNUM_TYPENAME_EGlpNumMultTo (w->coef[r], b->coef[i]);
		lp->iwork[r] = 1;
	}
	l = b->nzcnt;

	for (i = 0; i < a->nzcnt; i++)
	{
		r = a->indx[i];
		if (lp->iwork[r] == 0)
			w->indx[l++] = r;
		EGLPNUM_TYPENAME_EGlpNumAddTo (w->coef[r], a->coef[i]);
	}
	for (i = 0; i < l; i++)
	{
		r = w->indx[i];
		c->indx[i] = r;
		EGLPNUM_TYPENAME_EGlpNumCopy (c->coef[i], w->coef[r]);
		EGLPNUM_TYPENAME_EGlpNumZero (w->coef[r]);
		lp->iwork[r] = 0;
	}
	w->nzcnt = 0;
	c->nzcnt = l;
}

void EGLPNUM_TYPENAME_ILLfct_update_pfeas (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int lindex,
	EGLPNUM_TYPENAME_svector * srhs)
{
	int i, k, r;
	int col, nz = 0;
	int cbnd, f;
	int *perm = lp->upd.perm;
	int *ix = lp->upd.ix;
	int tctr = lp->upd.tctr;
	EGLPNUM_TYPE *t = lp->upd.t;
	EGLPNUM_TYPE tz, *dty, ntmp;
	EGLPNUM_TYPE *l, *x, *u, *pftol = &(lp->tol->ip_tol);

	EGLPNUM_TYPENAME_EGlpNumInitVar (tz);
	EGLPNUM_TYPENAME_EGlpNumInitVar (ntmp);
	dty = &(lp->upd.dty);
	EGLPNUM_TYPENAME_EGlpNumZero (*dty);
	EGLPNUM_TYPENAME_EGlpNumCopyAbs (tz, lp->upd.tz);
	EGLPNUM_TYPENAME_EGlpNumDivUiTo (tz, 100);
	EGLPNUM_TYPENAME_EGlpNumAddTo (tz, lp->upd.tz);
	ILL_IFTRACE ("%s:%d", __func__, tctr);
	for (i = 0; i < tctr && EGLPNUM_TYPENAME_EGlpNumIsLeq (t[perm[i]], tz); i++)
	{
		cbnd = ix[perm[i]] % 10;
		ILL_IFTRACE (":%d", cbnd);
		if (cbnd == BBOUND)
			continue;
		k = ix[perm[i]] / 10;
		r = lp->yjz.indx[k];
		ILL_IFTRACE (":%d:%d:%d", k, r, lp->iwork[r]);

		if (lp->iwork[r] != 1)
		{
			lp->iwork[r] = 1;
			x = &(lp->xbz[r]);
			col = lp->baz[r];
			l = &(lp->lz[col]);
			u = &(lp->uz[col]);

			if (r != lindex)
			{
				f = 0;
				EGLPNUM_TYPENAME_EGlpNumCopyDiff (ntmp, *l, *x);
				if (EGLPNUM_TYPENAME_EGlpNumIsNeqq (*l, EGLPNUM_TYPENAME_NINFTY) && EGLPNUM_TYPENAME_EGlpNumIsLess (*pftol, ntmp))
					f = -1;
				else
				{
					EGLPNUM_TYPENAME_EGlpNumCopyDiff (ntmp, *x, *u);
					if (EGLPNUM_TYPENAME_EGlpNumIsNeqq (*u, EGLPNUM_TYPENAME_INFTY) && EGLPNUM_TYPENAME_EGlpNumIsLess (*pftol, ntmp))
						f = 1;
				}

				ILL_IFTRACE (":%d:%d", f, lp->bfeas[r]);
				if (f != lp->bfeas[r])
				{
					srhs->indx[nz] = r;
					EGLPNUM_TYPENAME_EGlpNumSet (srhs->coef[nz], (double) (f - lp->bfeas[r]));
					EGLPNUM_TYPENAME_EGlpNumAddInnProdTo (*dty, srhs->coef[nz], lp->yjz.coef[k]);
					nz++;
					lp->bfeas[r] = f;
				}
			}
			else
			{
				lp->bfeas[r] = 0;
			}
		}
	}
	while (--i >= 0)
	{
		cbnd = ix[perm[i]] % 10;
		if (cbnd == BBOUND)
			continue;
		k = ix[perm[i]] / 10;
		r = lp->yjz.indx[k];
		lp->iwork[r] = 0;
	}
	srhs->nzcnt = nz;
	ILL_IFTRACE (":%d\n", nz);
	EGLPNUM_TYPENAME_EGlpNumClearVar (tz);
	EGLPNUM_TYPENAME_EGlpNumClearVar (ntmp);
}

void EGLPNUM_TYPENAME_ILLfct_compute_ppIzz (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_svector * srhs,
	EGLPNUM_TYPENAME_svector * ssoln)
{
	if (srhs->nzcnt != 0)
	{
		ILL_IFTRACE ("%s:\n", __func__);
		EGLPNUM_TYPENAME_ILLbasis_row_solve (lp, srhs, ssoln);
	}
}

void EGLPNUM_TYPENAME_ILLfct_update_ppI_prices (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_price_info * pinf,
	EGLPNUM_TYPENAME_svector * srhs,
	EGLPNUM_TYPENAME_svector * ssoln,
	int eindex,
	int lindex,
	const EGLPNUM_TYPE alpha)
{
	EGLPNUM_TYPE ntmp;

	EGLPNUM_TYPENAME_EGlpNumInitVar (ntmp);
	EGLPNUM_TYPENAME_EGlpNumCopy (ntmp, alpha);
	ILL_IFTRACE ("%s:\n", __func__);
	if (lindex == -1)
	{
		if (srhs->nzcnt != 0)
		{
			EGLPNUM_TYPENAME_ILLfct_update_pIpiz (lp, ssoln, EGLPNUM_TYPENAME_oneLpNum);
			if (pinf->p_strategy == COMPLETE_PRICING)
			{
				EGLPNUM_TYPENAME_ILLfct_compute_zA (lp, ssoln, &(lp->zA));
				EGLPNUM_TYPENAME_ILLfct_update_pIdz (lp, &(lp->zA), -1, EGLPNUM_TYPENAME_oneLpNum);
			}
		}
		else
		{
			if (pinf->p_strategy == COMPLETE_PRICING)
				EGLPNUM_TYPENAME_ILLprice_compute_dual_inf (lp, pinf, &eindex, 1, PRIMAL_PHASEI);
			else
				EGLPNUM_TYPENAME_ILLprice_update_mpartial_price (lp, pinf, PRIMAL_PHASEI, COL_PRICING);
			EGLPNUM_TYPENAME_EGlpNumClearVar (ntmp);
			return;
		}
	}
	else
	{
		if (srhs->nzcnt == 0)
		{
			EGLPNUM_TYPENAME_ILLfct_update_pIpiz (lp, &(lp->zz), ntmp);
			if (pinf->p_strategy == COMPLETE_PRICING)
				EGLPNUM_TYPENAME_ILLfct_update_pIdz (lp, &(lp->zA), eindex, ntmp);
		}
		else
		{
			EGLPNUM_TYPENAME_EGlpNumCopyFrac (ntmp, lp->upd.dty, lp->upd.piv);
			EGLPNUM_TYPENAME_EGlpNumSubTo (ntmp, alpha);
			EGLPNUM_TYPENAME_EGlpNumSign (ntmp);
			add_vectors (lp, ssoln, &(lp->zz), &(lp->zz), ntmp);
			EGLPNUM_TYPENAME_ILLfct_update_pIpiz (lp, &(lp->zz), EGLPNUM_TYPENAME_oneLpNum);
			if (pinf->p_strategy == COMPLETE_PRICING)
			{
				EGLPNUM_TYPENAME_ILLfct_compute_zA (lp, &(lp->zz), &(lp->zA));
				EGLPNUM_TYPENAME_ILLfct_update_pIdz (lp, &(lp->zA), eindex, EGLPNUM_TYPENAME_oneLpNum);
			}
		}
		EGLPNUM_TYPENAME_EGlpNumSet (lp->pIdz[eindex], (double) (lp->upd.fs));
		EGLPNUM_TYPENAME_EGlpNumAddTo (lp->pIdz[eindex], ntmp);
		EGLPNUM_TYPENAME_EGlpNumSign (lp->pIdz[eindex]);
	}
	if (pinf->p_strategy == COMPLETE_PRICING)
	{
		EGLPNUM_TYPENAME_ILLprice_compute_dual_inf (lp, pinf, lp->zA.indx, lp->zA.nzcnt,
															 PRIMAL_PHASEI);
		if (eindex > -1)
			EGLPNUM_TYPENAME_ILLprice_compute_dual_inf (lp, pinf, &eindex, 1, PRIMAL_PHASEI);
		EGLPNUM_TYPENAME_ILLfct_update_counts (lp, CNT_ZARAVG, lp->zA.nzcnt, EGLPNUM_TYPENAME_zeroLpNum);
	}
	else
		EGLPNUM_TYPENAME_ILLprice_update_mpartial_price (lp, pinf, PRIMAL_PHASEI, COL_PRICING);
	EGLPNUM_TYPENAME_EGlpNumClearVar (ntmp);
	return;
}

void EGLPNUM_TYPENAME_ILLfct_update_dfeas (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int eindex,
	EGLPNUM_TYPENAME_svector * srhs)
{
	int i, j, k, c;
	int cbnd, col, nz = 0;
	int vs, vt, f;
	int delta;
	int *perm = lp->upd.perm;
	int *ix = lp->upd.ix;
	int tctr = lp->upd.tctr;
	int mcnt, mbeg;
	EGLPNUM_TYPE *t = lp->upd.t;
	EGLPNUM_TYPE *w = lp->work.coef;
	EGLPNUM_TYPE tz;
	EGLPNUM_TYPE *dty = &(lp->upd.dty);
	EGLPNUM_TYPE *dftol = &(lp->tol->id_tol);
	EGLPNUM_TYPE dj;

	EGLPNUM_TYPENAME_EGlpNumInitVar (dj);
	EGLPNUM_TYPENAME_EGlpNumInitVar (tz);
	EGLPNUM_TYPENAME_EGlpNumZero (*dty);
	EGLPNUM_TYPENAME_EGlpNumCopy (tz, lp->upd.tz);
	EGLPNUM_TYPENAME_EGlpNumMultUiTo (tz, 101);
	EGLPNUM_TYPENAME_EGlpNumDivUiTo (tz, 100);

	for (j = 0; j < tctr && EGLPNUM_TYPENAME_EGlpNumIsLeq (t[perm[j]], tz); j++)
	{
		k = ix[perm[j]] / 10;
		c = lp->zA.indx[k];

		if (lp->iwork[c] != 1)
		{
			lp->iwork[c] = 1;
			cbnd = ix[perm[j]] % 10;
			col = lp->nbaz[c];
			EGLPNUM_TYPENAME_EGlpNumCopy (dj, lp->dz[c]);
			vs = lp->vstat[col];
			vt = lp->vtype[col];

			if (cbnd == BSKIP)
			{
				if (!EGLPNUM_TYPENAME_EGlpNumIsNeqZero (dj, *dftol));
				else if (EGLPNUM_TYPENAME_EGlpNumIsLessZero (dj) && vs == STAT_LOWER)
					lp->vstat[col] = STAT_UPPER;
				else if (EGLPNUM_TYPENAME_EGlpNumIsGreatZero (dj) && vs == STAT_UPPER)
					lp->vstat[col] = STAT_LOWER;
			}
			else if (c != eindex)
			{
				if (!EGLPNUM_TYPENAME_EGlpNumIsNeqZero (dj, *dftol))
					f = 0;
				else if (EGLPNUM_TYPENAME_EGlpNumIsLessZero (dj) &&
								 (vs == STAT_LOWER || vs == STAT_ZERO))
					f = -1;
				else if (EGLPNUM_TYPENAME_EGlpNumIsGreatZero (dj) &&
								 (vs == STAT_UPPER || vs == STAT_ZERO))
					f = 1;
				else
					f = 0;

				if (f != lp->dfeas[c])
				{
					delta = f - lp->dfeas[c];
					mcnt = lp->matcnt[col];
					mbeg = lp->matbeg[col];
					EGLPNUM_TYPENAME_EGlpNumSet (dj, (double) (delta));
					for (i = 0; i < mcnt; i++)
						EGLPNUM_TYPENAME_EGlpNumAddInnProdTo (w[lp->matind[mbeg + i]], dj,
																 lp->matval[mbeg + i]);
					EGLPNUM_TYPENAME_EGlpNumAddInnProdTo (*dty, dj, lp->zA.coef[k]);
					nz = 1;
					lp->dfeas[c] = f;
				}
			}
			else
			{
				lp->dfeas[c] = 0;
			}
		}
	}
	while (--j >= 0)
	{
		k = ix[perm[j]] / 10;
		c = lp->zA.indx[k];
		lp->iwork[c] = 0;
	}

	if (nz)
	{
		for (i = 0, nz = 0; i < lp->nrows; i++)
			if (EGLPNUM_TYPENAME_EGlpNumIsNeqqZero (w[i]))
			{
				EGLPNUM_TYPENAME_EGlpNumCopy (srhs->coef[nz], w[i]);
				srhs->indx[nz] = i;
				nz++;
				EGLPNUM_TYPENAME_EGlpNumZero (w[i]);
			}
	}

	srhs->nzcnt = nz;
	EGLPNUM_TYPENAME_EGlpNumClearVar (dj);
	EGLPNUM_TYPENAME_EGlpNumClearVar (tz);
}

void EGLPNUM_TYPENAME_ILLfct_compute_dpIy (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_svector * srhs,
	EGLPNUM_TYPENAME_svector * ssoln)
{
	if (srhs->nzcnt != 0)
	{
		EGLPNUM_TYPENAME_ILLbasis_column_solve (lp, srhs, ssoln);
	}
}

void EGLPNUM_TYPENAME_ILLfct_update_dpI_prices (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_price_info * pinf,
	EGLPNUM_TYPENAME_svector * srhs,
	EGLPNUM_TYPENAME_svector * ssoln,
	int lindex,
	EGLPNUM_TYPE alpha)
{
	int i;
	EGLPNUM_TYPE ntmp;

	EGLPNUM_TYPENAME_EGlpNumInitVar (ntmp);
	EGLPNUM_TYPENAME_EGlpNumZero (ntmp);

	if (srhs->nzcnt == 0)
	{
		EGLPNUM_TYPENAME_ILLfct_update_xz (lp, alpha, -1, -1);
	}
	else
	{
		EGLPNUM_TYPENAME_EGlpNumCopyFrac (ntmp, lp->upd.dty, lp->upd.piv);
		EGLPNUM_TYPENAME_EGlpNumAddTo (ntmp, alpha);
		EGLPNUM_TYPENAME_EGlpNumSign (ntmp);
		add_vectors (lp, ssoln, &(lp->yjz), &(lp->yjz), ntmp);
		EGLPNUM_TYPENAME_EGlpNumSign (ntmp);
		for (i = 0; i < lp->yjz.nzcnt; i++)
			EGLPNUM_TYPENAME_EGlpNumAddTo (lp->xbz[lp->yjz.indx[i]], lp->yjz.coef[i]);
	}
	EGLPNUM_TYPENAME_EGlpNumSet (lp->xbz[lindex], ((double) (-lp->upd.fs)));
	EGLPNUM_TYPENAME_EGlpNumAddTo (lp->xbz[lindex], ntmp);

	if (pinf->d_strategy == COMPLETE_PRICING)
	{
		EGLPNUM_TYPENAME_ILLprice_compute_primal_inf (lp, pinf, lp->yjz.indx, lp->yjz.nzcnt,
																 DUAL_PHASEI);
		EGLPNUM_TYPENAME_ILLprice_compute_primal_inf (lp, pinf, &lindex, 1, DUAL_PHASEI);
		EGLPNUM_TYPENAME_ILLfct_update_counts (lp, CNT_YRAVG, lp->yjz.nzcnt, EGLPNUM_TYPENAME_zeroLpNum);
	}
	else
		EGLPNUM_TYPENAME_ILLprice_update_mpartial_price (lp, pinf, DUAL_PHASEI, ROW_PRICING);
	EGLPNUM_TYPENAME_EGlpNumClearVar (ntmp);
}

void EGLPNUM_TYPENAME_ILLfct_update_dIIfeas (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int eindex,
	EGLPNUM_TYPENAME_svector * srhs)
{
	int j, k;
	int col, indx, vs;
	int *perm = lp->upd.perm;
	int *ix = lp->upd.ix;
	int tctr = lp->upd.tctr;
	EGLPNUM_TYPE *zAj, *l, *u;
	EGLPNUM_TYPE *dty = &(lp->upd.dty);
	EGLPNUM_TYPE *t_max = &(lp->upd.tz);
	EGLPNUM_TYPE *t = lp->upd.t;
	EGLPNUM_TYPE delta;
	EGLPNUM_TYPENAME_svector a;

	EGLPNUM_TYPENAME_EGlpNumInitVar (delta);
	EGLPNUM_TYPENAME_EGlpNumZero (delta);
	EGLPNUM_TYPENAME_EGlpNumZero (*dty);

	srhs->nzcnt = 0;
	for (j = 0; j < tctr && EGLPNUM_TYPENAME_EGlpNumIsLeq (t[perm[j]], *t_max); j++)
	{
		k = ix[perm[j]];
		indx = lp->zA.indx[k];

		if (indx != eindex)
		{
			zAj = &(lp->zA.coef[k]);
			col = lp->nbaz[indx];
			l = &(lp->lz[col]);
			u = &(lp->uz[col]);
			vs = lp->vstat[col];
			if (vs == STAT_UPPER)
				EGLPNUM_TYPENAME_EGlpNumCopyDiff (delta, *l, *u);
			else
				EGLPNUM_TYPENAME_EGlpNumCopyDiff (delta, *u, *l);
			EGLPNUM_TYPENAME_EGlpNumAddInnProdTo (*dty, delta, *zAj);
			lp->vstat[col] = (vs == STAT_UPPER) ? STAT_LOWER : STAT_UPPER;

			a.nzcnt = lp->matcnt[col];
			a.indx = &(lp->matind[lp->matbeg[col]]);
			a.coef = &(lp->matval[lp->matbeg[col]]);
			add_vectors (lp, srhs, &a, srhs, delta);
		}
	}
	EGLPNUM_TYPENAME_EGlpNumClearVar (delta);
}

void EGLPNUM_TYPENAME_ILLfct_compute_dpIIy (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_svector * srhs,
	EGLPNUM_TYPENAME_svector * ssoln)
{
	if (srhs->nzcnt != 0)
	{
		EGLPNUM_TYPENAME_ILLbasis_column_solve (lp, srhs, ssoln);
	}
}

void EGLPNUM_TYPENAME_ILLfct_update_dpII_prices (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_price_info * pinf,
	EGLPNUM_TYPENAME_svector * srhs,
	EGLPNUM_TYPENAME_svector * ssoln,
	/*int eindex,*/
	int lindex,
	EGLPNUM_TYPE eval,
	EGLPNUM_TYPE alpha)
{
	int i;
	EGLPNUM_TYPENAME_svector *u;

	if (srhs->nzcnt == 0)
	{
		EGLPNUM_TYPENAME_ILLfct_update_xz (lp, alpha, -1, -1);
		u = &(lp->yjz);
	}
	else
	{
		if (ssoln->nzcnt != 0)
			for (i = 0; i < ssoln->nzcnt; i++)
				EGLPNUM_TYPENAME_EGlpNumSubTo (lp->xbz[ssoln->indx[i]], ssoln->coef[i]);
		EGLPNUM_TYPENAME_ILLfct_update_xz (lp, alpha, -1, -1);
		add_vectors (lp, ssoln, &(lp->yjz), ssoln, EGLPNUM_TYPENAME_oneLpNum);
		u = ssoln;
	}
	EGLPNUM_TYPENAME_EGlpNumCopySum (lp->xbz[lindex], eval, alpha);

	if (pinf->d_strategy == COMPLETE_PRICING)
	{
		EGLPNUM_TYPENAME_ILLprice_compute_primal_inf (lp, pinf, u->indx, u->nzcnt, DUAL_PHASEII);
		EGLPNUM_TYPENAME_ILLprice_compute_primal_inf (lp, pinf, &lindex, 1, DUAL_PHASEII);
		EGLPNUM_TYPENAME_ILLfct_update_counts (lp, CNT_YRAVG, u->nzcnt, EGLPNUM_TYPENAME_zeroLpNum);
	}
	else
		EGLPNUM_TYPENAME_ILLprice_update_mpartial_price (lp, pinf, DUAL_PHASEII, ROW_PRICING);
}

int EGLPNUM_TYPENAME_ILLfct_test_pivot (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int indx,
	int indxtype,
	EGLPNUM_TYPE piv_val)
{
	int i;
	EGLPNUM_TYPE pval, ntmp;

	EGLPNUM_TYPENAME_EGlpNumInitVar (pval);
	EGLPNUM_TYPENAME_EGlpNumInitVar (ntmp);
	EGLPNUM_TYPENAME_EGlpNumZero (pval);

	if (indxtype == ROW_PIVOT)
	{
		for (i = 0; i < lp->yjz.nzcnt; i++)
			if (lp->yjz.indx[i] == indx)
			{
				EGLPNUM_TYPENAME_EGlpNumCopy (pval, lp->yjz.coef[i]);
				break;
			}
	}
	else
	{
		for (i = 0; i < lp->zA.nzcnt; i++)
			if (lp->zA.indx[i] == indx)
			{
				EGLPNUM_TYPENAME_EGlpNumCopy (pval, lp->zA.coef[i]);
				break;
			}
	}
	EGLPNUM_TYPENAME_EGlpNumCopyDiff (ntmp, pval, piv_val);
	EGLPNUM_TYPENAME_EGlpNumDivTo (ntmp, piv_val);
	if (EGLPNUM_TYPENAME_EGlpNumIsLessZero (ntmp))
		EGLPNUM_TYPENAME_EGlpNumSign (ntmp);
	if (EGLPNUM_TYPENAME_EGlpNumIsLess (EGLPNUM_TYPENAME_ALTPIV_TOLER, ntmp))
	{
#if FCT_DEBUG > 1
		if (indxtype == ROW_PIVOT)
			QSlog("y_i = %.8f, z_j = %.8f %la %la", EGLPNUM_TYPENAME_EGlpNumToLf (pval),
									EGLPNUM_TYPENAME_EGlpNumToLf (piv_val), EGLPNUM_TYPENAME_EGlpNumToLf (EGLPNUM_TYPENAME_ALTPIV_TOLER),
									EGLPNUM_TYPENAME_EGlpNumToLf (ntmp));
		else
			QSlog("z_j = %.8f, y_i = %.8f", EGLPNUM_TYPENAME_EGlpNumToLf (pval),
									EGLPNUM_TYPENAME_EGlpNumToLf (piv_val));
#endif
		EGLPNUM_TYPENAME_EGlpNumClearVar (ntmp);
		EGLPNUM_TYPENAME_EGlpNumClearVar (pval);
		return 1;
	}
	EGLPNUM_TYPENAME_EGlpNumClearVar (pval);
	EGLPNUM_TYPENAME_EGlpNumClearVar (ntmp);
	return 0;
}

#if FCT_DEBUG > 0

void EGLPNUM_TYPENAME_fct_test_workvector (
	EGLPNUM_TYPENAME_lpinfo * lp)
{
	int i, err = 0;

	for (i = 0; i < lp->ncols; i++)
	{
		if (EGLPNUM_TYPENAME_EGlpNumIsNeqqZero (lp->work.coef[i]))
		{
			err++;
			EGLPNUM_TYPENAME_EGlpNumZero (lp->work.coef[i]);
		}
		if (lp->iwork[i] != 0)
		{
			err++;
			lp->iwork[i] = 0;
		}
	}
	if (err)
		QSlog("bad work vector, err=%d", err);
}

void EGLPNUM_TYPENAME_fct_test_pfeasible (
	EGLPNUM_TYPENAME_lpinfo * lp)
{
	int i, col;
	int err = 0;
	EGLPNUM_TYPE *ftol = &(lp->tol->pfeas_tol);

	for (i = 0; i < lp->nrows; i++)
	{
		col = lp->baz[i];

		if (EGLPNUM_TYPENAME_EGlpNumIsNeqq (lp->uz[col], EGLPNUM_TYPENAME_INFTY)
				&& EGLPNUM_TYPENAME_EGlpNumIsSumLess (*ftol, lp->uz[col], lp->xbz[i]))
		{
			if (lp->bfeas[i] != 1)
			{
				err++;
				lp->bfeas[i] = 1;
			}
		}
		else if (EGLPNUM_TYPENAME_EGlpNumIsNeqq (lp->lz[col], EGLPNUM_TYPENAME_NINFTY)
						 && EGLPNUM_TYPENAME_EGlpNumIsSumLess (lp->xbz[i], *ftol, lp->lz[col]))
		{
			if (lp->bfeas[i] != -1)
			{
				err++;
				lp->bfeas[i] = -1;
			}
		}
		/* else if (lp->bfeas[i] != 0) {err++; lp->bfeas[i] = 0;} */
	}
	if (err != 0)
		QSlog("test_pfeas err =%d", err);
}

void EGLPNUM_TYPENAME_fct_test_dfeasible (
	EGLPNUM_TYPENAME_lpinfo * lp)
{
	int j, col;
	int err = 0;
	EGLPNUM_TYPE *ftol = &(lp->tol->dfeas_tol);
	EGLPNUM_TYPE mftol[1];

	EGLPNUM_TYPENAME_EGlpNumInitVar (mftol[0]);
	EGLPNUM_TYPENAME_EGlpNumCopyNeg (mftol[0], *ftol);

	for (j = 0; j < lp->nnbasic; j++)
	{
		col = lp->nbaz[j];

		if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED)
			continue;
		if (EGLPNUM_TYPENAME_EGlpNumIsLess (lp->dz[j], mftol[0]) &&
				(lp->vstat[col] == STAT_LOWER || lp->vstat[col] == STAT_ZERO))
		{
			if (lp->dfeas[j] != -1)
			{
				err++;
				lp->dfeas[j] = -1;
			}
		}
		if (EGLPNUM_TYPENAME_EGlpNumIsLess (*ftol, lp->dz[j]) &&
				(lp->vstat[col] == STAT_UPPER || lp->vstat[col] == STAT_ZERO))
		{
			if (lp->dfeas[j] != 1)
			{
				err++;
				lp->dfeas[j] = 1;
			}
		}
		/* else if (lp->dfeas[j] != 0) {err++; lp->dfeas[j] = 0;} */
	}
	if (err != 0)
		QSlog("test_dfeas err =%d", err);
}

void EGLPNUM_TYPENAME_fct_test_pI_x (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_price_info * p)
{
	int i;
	int ern = 0;
	EGLPNUM_TYPE *x;
	EGLPNUM_TYPE err, diff;

	EGLPNUM_TYPENAME_EGlpNumInitVar (err);
	EGLPNUM_TYPENAME_EGlpNumInitVar (diff);
	EGLPNUM_TYPENAME_EGlpNumZero (err);
	x = EGLPNUM_TYPENAME_EGlpNumAllocArray (lp->nrows);

	for (i = 0; i < lp->nrows; i++)
		EGLPNUM_TYPENAME_EGlpNumCopy (x[i], lp->xbz[i]);
	EGLPNUM_TYPENAME_ILLfct_compute_phaseI_xbz (lp);
	for (i = 0; i < lp->nrows; i++)
	{
		EGLPNUM_TYPENAME_EGlpNumCopyDiff (diff, x[i], lp->xbz[i]);
		if (EGLPNUM_TYPENAME_EGlpNumIsLessZero (diff))
			EGLPNUM_TYPENAME_EGlpNumSign (diff);
		if (EGLPNUM_TYPENAME_EGlpNumIsLess (EGLPNUM_TYPENAME_PFEAS_TOLER, diff))
		{
			EGLPNUM_TYPENAME_EGlpNumAddTo (err, diff);
			ern++;
			QSlog("bad i = %d", i);
		}
	}
	if (EGLPNUM_TYPENAME_EGlpNumIsNeqqZero (err))
		QSlog("dI x err = %.7f, ern = %d", EGLPNUM_TYPENAME_EGlpNumToLf (err), ern);
	EGLPNUM_TYPENAME_ILLprice_compute_primal_inf (lp, p, NULL, 0, DUAL_PHASEI);
	EGLPNUM_TYPENAME_EGlpNumFreeArray (x);
	EGLPNUM_TYPENAME_EGlpNumClearVar (diff);
	EGLPNUM_TYPENAME_EGlpNumClearVar (err);
}

void EGLPNUM_TYPENAME_fct_test_pII_x (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_price_info * p)
{
	int i;
	int ern = 0;
	EGLPNUM_TYPE *x;
	EGLPNUM_TYPE err, diff;

	EGLPNUM_TYPENAME_EGlpNumInitVar (err);
	EGLPNUM_TYPENAME_EGlpNumInitVar (diff);
	EGLPNUM_TYPENAME_EGlpNumZero (err);
	x = EGLPNUM_TYPENAME_EGlpNumAllocArray (lp->nrows);

	for (i = 0; i < lp->nrows; i++)
		EGLPNUM_TYPENAME_EGlpNumCopy (x[i], lp->xbz[i]);
	EGLPNUM_TYPENAME_ILLfct_compute_xbz (lp);
	for (i = 0; i < lp->nrows; i++)
	{
		EGLPNUM_TYPENAME_EGlpNumCopyDiff (diff, x[i], lp->xbz[i]);
		if (EGLPNUM_TYPENAME_EGlpNumIsLessZero (diff ))
			EGLPNUM_TYPENAME_EGlpNumSign (diff);
		if (EGLPNUM_TYPENAME_EGlpNumIsLess (EGLPNUM_TYPENAME_PFEAS_TOLER, diff))
		{
			EGLPNUM_TYPENAME_EGlpNumAddTo (err, diff);
			ern++;
			QSlog("bad i = %d", i);
		}
	}
	if (EGLPNUM_TYPENAME_EGlpNumIsNeqqZero (err))
		QSlog("dII x err = %.7f, ern = %d", EGLPNUM_TYPENAME_EGlpNumToLf (err), ern);
	EGLPNUM_TYPENAME_ILLprice_compute_primal_inf (lp, p, NULL, 0, DUAL_PHASEII);
	EGLPNUM_TYPENAME_EGlpNumFreeArray (x);
	EGLPNUM_TYPENAME_EGlpNumClearVar (diff);
	EGLPNUM_TYPENAME_EGlpNumClearVar (err);
}

void EGLPNUM_TYPENAME_fct_test_pI_pi_dz (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_price_info * p)
{
	int i;
	int ern = 0;
	EGLPNUM_TYPE *pidz;
	EGLPNUM_TYPE err, diff;

	EGLPNUM_TYPENAME_EGlpNumInitVar (err);
	EGLPNUM_TYPENAME_EGlpNumInitVar (diff);
	pidz = EGLPNUM_TYPENAME_EGlpNumAllocArray (lp->ncols);
	EGLPNUM_TYPENAME_EGlpNumZero (err);

	for (i = 0; i < lp->nrows; i++)
		EGLPNUM_TYPENAME_EGlpNumCopy (pidz[i], lp->pIpiz[i]);
	EGLPNUM_TYPENAME_ILLfct_compute_phaseI_piz (lp);
	for (i = 0; i < lp->nrows; i++)
	{
		EGLPNUM_TYPENAME_EGlpNumCopyDiff (diff, pidz[i], lp->pIpiz[i]);
		if (EGLPNUM_TYPENAME_EGlpNumIsLessZero (diff))
			EGLPNUM_TYPENAME_EGlpNumSign (diff);
		if (EGLPNUM_TYPENAME_EGlpNumIsLess (EGLPNUM_TYPENAME_DFEAS_TOLER, diff))
		{
			EGLPNUM_TYPENAME_EGlpNumAddTo (err, diff);
			ern++;
		}
	}
	if (EGLPNUM_TYPENAME_EGlpNumIsNeqqZero (err))
		QSlog("pI pi err = %.7f, ern = %d", EGLPNUM_TYPENAME_EGlpNumToLf (err), ern);

	EGLPNUM_TYPENAME_EGlpNumZero (err);
	ern = 0;
	for (i = 0; i < lp->nnbasic; i++)
		EGLPNUM_TYPENAME_EGlpNumCopy (pidz[i], lp->pIdz[i]);
	EGLPNUM_TYPENAME_ILLfct_compute_phaseI_dz (lp);
	for (i = 0; i < lp->nnbasic; i++)
	{
		EGLPNUM_TYPENAME_EGlpNumCopyDiff (diff, pidz[i], lp->pIdz[i]);
		if (EGLPNUM_TYPENAME_EGlpNumIsLessZero (diff))
			EGLPNUM_TYPENAME_EGlpNumSign (diff);
		if (EGLPNUM_TYPENAME_EGlpNumIsLess (EGLPNUM_TYPENAME_DFEAS_TOLER, diff))
		{
			EGLPNUM_TYPENAME_EGlpNumAddTo (err, diff);
			ern++;
		}
	}
	if (EGLPNUM_TYPENAME_EGlpNumIsNeqqZero (err))
		QSlog("pI dz err = %.7f, ern = %d", EGLPNUM_TYPENAME_EGlpNumToLf (err), ern);
	EGLPNUM_TYPENAME_ILLprice_compute_dual_inf (lp, p, NULL, 0, PRIMAL_PHASEI);
	EGLPNUM_TYPENAME_EGlpNumClearVar (err);
	EGLPNUM_TYPENAME_EGlpNumClearVar (diff);
	EGLPNUM_TYPENAME_EGlpNumFreeArray (pidz);
}

void EGLPNUM_TYPENAME_fct_test_pII_pi_dz (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_price_info * p)
{
	int i;
	int ern = 0;
	EGLPNUM_TYPE *pidz;
	EGLPNUM_TYPE err, diff;

	EGLPNUM_TYPENAME_EGlpNumInitVar (err);
	EGLPNUM_TYPENAME_EGlpNumInitVar (diff);
	EGLPNUM_TYPENAME_EGlpNumZero (err);
	pidz = EGLPNUM_TYPENAME_EGlpNumAllocArray (lp->ncols);

	for (i = 0; i < lp->nrows; i++)
		EGLPNUM_TYPENAME_EGlpNumCopy (pidz[i], lp->piz[i]);
	EGLPNUM_TYPENAME_ILLfct_compute_piz (lp);
	for (i = 0; i < lp->nrows; i++)
	{
		EGLPNUM_TYPENAME_EGlpNumCopyDiff (diff, pidz[i], lp->piz[i]);
		if (EGLPNUM_TYPENAME_EGlpNumIsLessZero (diff))
			EGLPNUM_TYPENAME_EGlpNumSign (diff);
		if (EGLPNUM_TYPENAME_EGlpNumIsLess (EGLPNUM_TYPENAME_DFEAS_TOLER, diff))
		{
			EGLPNUM_TYPENAME_EGlpNumAddTo (err, diff);
			ern++;
		}
	}
	if (EGLPNUM_TYPENAME_EGlpNumIsNeqqZero (err))
		QSlog("pII pi err = %.7f, ern = %d", EGLPNUM_TYPENAME_EGlpNumToLf (err), ern);

	EGLPNUM_TYPENAME_EGlpNumZero (err);
	ern = 0;
	for (i = 0; i < lp->nnbasic; i++)
		EGLPNUM_TYPENAME_EGlpNumCopy (pidz[i], lp->dz[i]);
	EGLPNUM_TYPENAME_ILLfct_compute_dz (lp);
	for (i = 0; i < lp->nnbasic; i++)
	{
		EGLPNUM_TYPENAME_EGlpNumCopyDiff (diff, pidz[i], lp->dz[i]);
		if (EGLPNUM_TYPENAME_EGlpNumIsLessZero (diff))
			EGLPNUM_TYPENAME_EGlpNumSign (diff);
		if (EGLPNUM_TYPENAME_EGlpNumIsLess (EGLPNUM_TYPENAME_DFEAS_TOLER, diff))
		{
			EGLPNUM_TYPENAME_EGlpNumAddTo (err, diff);
			ern++;
		}
	}
	if (EGLPNUM_TYPENAME_EGlpNumIsNeqqZero (err))
		QSlog("pII dz err = %.7f, ern = %d", EGLPNUM_TYPENAME_EGlpNumToLf (err), ern);
	/*
	 * EGLPNUM_TYPENAME_ILLprice_compute_dual_inf (lp, p, NULL, 0, PRIMAL_PHASEII);
	 */
	EGLPNUM_TYPENAME_EGlpNumClearVar (err);
	EGLPNUM_TYPENAME_EGlpNumClearVar (diff);
	EGLPNUM_TYPENAME_EGlpNumFreeArray (pidz);
}

#endif
