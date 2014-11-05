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

#include "qs_config.h"
#include "iqsutil.h"
#include "lpdefs.h"
#include "stddefs.h"
#include "basis.h"
#include "fct.h"
#include "price.h"
#include "ratio.h"
#include "dstruct.h"

bndinfo *ILLfct_new_bndinfo (
	void)
{
	bndinfo *nbnd = (bndinfo *) malloc (sizeof (bndinfo));

	if (!nbnd)
	{
		fprintf (stderr, "not enough memory, in %s\n", __func__);
		exit (1);
	}
	EGlpNumInitVar ((nbnd->pbound));
	EGlpNumInitVar ((nbnd->cbound));
	return nbnd;
}

void ILLfct_free_bndinfo (
	bndinfo * binfo)
{
	EGlpNumClearVar ((binfo->pbound));
	EGlpNumClearVar ((binfo->cbound));
	ILL_IFFREE (binfo, bndinfo);
	return;
}

static int compute_zA1 (
	lpinfo * lp,
	svector * z,
	svector * zA,
	EGlpNum_t ztoler),
/*
  compute_zA2 (lpinfo * lp,
							 svector * z,
							 svector * zA,
							 const EGlpNum_t* ztoler), */
  compute_zA3 (
	lpinfo * lp,
	svector * z,
	svector * zA,
	EGlpNum_t ztoler),
  expand_var_bounds (
	lpinfo * lp,
	EGlpNum_t ftol,
	int *chgb),
  expand_var_coefs (
	lpinfo * lp,
	EGlpNum_t ftol,
	int *chgc);

static void update_piv_values (
	count_struct * c,
	int phase,
	const EGlpNum_t piv),
/*  copy_vectors (svector * a,
								svector * b),*/
  add_vectors (
	lpinfo * lp,
	svector * a,
	svector * b,
	svector * c,
	const EGlpNum_t t);

static double my_rand (
	int bound,
	ILLrandstate * r);


void ILLfct_load_workvector (
	lpinfo * lp,
	svector * s)
{
	int i;

	for (i = 0; i < s->nzcnt; i++)
	{
		lp->work.indx[i] = s->indx[i];
		EGlpNumCopy (lp->work.coef[s->indx[i]], s->coef[i]);
	}
	lp->work.nzcnt = s->nzcnt;
}

void ILLfct_zero_workvector (
	lpinfo * lp)
{
	int i;

	for (i = 0; i < lp->work.nzcnt; i++)
		EGlpNumZero (lp->work.coef[lp->work.indx[i]]);
	lp->work.nzcnt = 0;
}

void ILLfct_set_variable_type (
	lpinfo * lp)
{
	int j;

	for (j = 0; j < lp->ncols; j++)
	{

		if (lp->matcnt[j] == 1 && lp->O->rowmap[lp->matind[lp->matbeg[j]]] == j)
			lp->vclass[j] = CLASS_LOGICAL;
		else
			lp->vclass[j] = CLASS_STRUCT;
		switch ((EGlpNumIsEqqual (lp->uz[j], INFTY) ? 1U : 0U) |
						(EGlpNumIsEqqual (lp->lz[j], NINFTY) ? 2U : 0U))
		{
		case 0:
			if (EGlpNumIsLess (lp->lz[j], lp->uz[j]))
				lp->vtype[j] = VBOUNDED;
			else if (!EGlpNumIsNeqqZero (lp->lz[j]) &&
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

void ILLfct_compute_pobj (
	lpinfo * lp)
{
	int i, j;
	int col;
	EGlpNum_t sum;

	EGlpNumInitVar (sum);
	EGlpNumZero (sum);

	for (i = 0; i < lp->nrows; i++)
		EGlpNumAddInnProdTo (sum, lp->cz[lp->baz[i]], lp->xbz[i]);

	for (j = 0; j < lp->nnbasic; j++)
	{
		col = lp->nbaz[j];
		if (lp->vstat[col] == STAT_UPPER)
			EGlpNumAddInnProdTo (sum, lp->cz[col], lp->uz[col]);
		else if (lp->vstat[col] == STAT_LOWER)
			EGlpNumAddInnProdTo (sum, lp->cz[col], lp->lz[col]);
	}
	EGlpNumCopy (lp->pobjval, sum);
	EGlpNumCopy (lp->objval, sum);
	EGlpNumClearVar (sum);
}

void ILLfct_compute_dobj (
	lpinfo * lp)
{
	int i, j;
	int col;
	EGlpNum_t sum;

	EGlpNumInitVar (sum);
	EGlpNumZero (sum);

	for (i = 0; i < lp->nrows; i++)
		EGlpNumAddInnProdTo (sum, lp->piz[i], lp->bz[i]);

	for (j = 0; j < lp->nnbasic; j++)
	{
		col = lp->nbaz[j];
		if (lp->vstat[col] == STAT_UPPER)
			EGlpNumAddInnProdTo (sum, lp->dz[j], lp->uz[col]);
		else if (lp->vstat[col] == STAT_LOWER)
			EGlpNumAddInnProdTo (sum, lp->dz[j], lp->lz[col]);
	}
	EGlpNumCopy (lp->dobjval, sum);
	EGlpNumCopy (lp->objval, sum);
	EGlpNumClearVar (sum);
}

void ILLfct_compute_xbz (
	lpinfo * lp)
{
	int i, j, r;
	int col, mcnt, mbeg;
	svector *srhs = &(lp->srhs);
	svector *ssoln = &(lp->ssoln);
	EGlpNum_t xval;

	EGlpNumInitVar (xval);

	for (i = 0; i < lp->nrows; i++)
	{
		EGlpNumZero (lp->xbz[i]);
		EGlpNumCopy (srhs->coef[i], lp->bz[i]);
	}
	for (j = 0; j < lp->nnbasic; j++)
	{
		col = lp->nbaz[j];
		EGlpNumZero (xval);
		if (lp->vstat[col] == STAT_UPPER && EGlpNumIsNeqqZero (lp->uz[col]))
			EGlpNumCopy (xval, lp->uz[col]);
		else if (lp->vstat[col] == STAT_LOWER && EGlpNumIsNeqqZero (lp->lz[col]))
			EGlpNumCopy (xval, lp->lz[col]);

		if (EGlpNumIsNeqqZero (xval))
		{
			mcnt = lp->matcnt[col];
			mbeg = lp->matbeg[col];
			for (i = 0; i < mcnt; i++)
				EGlpNumSubInnProdTo (srhs->coef[lp->matind[mbeg + i]], xval,
														 lp->matval[mbeg + i]);
		}
	}
	for (i = 0, r = 0; i < lp->nrows; i++)
		if (EGlpNumIsNeqqZero (srhs->coef[i]))
		{
			EGlpNumCopy (srhs->coef[r], srhs->coef[i]);
			srhs->indx[r] = i;
			r++;
		}
	srhs->nzcnt = r;

	ILLbasis_column_solve (lp, srhs, ssoln);
	for (i = 0; i < ssoln->nzcnt; i++)
		EGlpNumCopy (lp->xbz[ssoln->indx[i]], ssoln->coef[i]);
	EGlpNumClearVar (xval);
}

void ILLfct_compute_piz (
	lpinfo * lp)
{
	int i, r;
	svector *srhs = &(lp->srhs);
	svector *ssoln = &(lp->ssoln);

	for (i = 0, r = 0; i < lp->nrows; i++)
	{
		EGlpNumZero (lp->piz[i]);
		if (EGlpNumIsNeqqZero (lp->cz[lp->baz[i]]))
		{
			srhs->indx[r] = i;
			EGlpNumCopy (srhs->coef[r], lp->cz[lp->baz[i]]);
			r++;
		}
	}
	srhs->nzcnt = r;

	ILLbasis_row_solve (lp, srhs, ssoln);
	for (i = 0; i < ssoln->nzcnt; i++)
		EGlpNumCopy (lp->piz[ssoln->indx[i]], ssoln->coef[i]);
}

void ILLfct_compute_dz (
	lpinfo * lp)
{
	int i, j;
	int col;
	int mcnt, mbeg;
	EGlpNum_t sum;

	EGlpNumInitVar (sum);

	for (j = 0; j < lp->nnbasic; j++)
	{
		EGlpNumZero (sum);
		col = lp->nbaz[j];
		mcnt = lp->matcnt[col];
		mbeg = lp->matbeg[col];
		for (i = 0; i < mcnt; i++)
			EGlpNumAddInnProdTo (sum, lp->piz[lp->matind[mbeg + i]],
													 lp->matval[mbeg + i]);
		EGlpNumCopyDiff (lp->dz[j], lp->cz[col], sum);
	}
	EGlpNumClearVar (sum);
}

void ILLfct_compute_phaseI_xbz (
	lpinfo * lp)
{
	int i, j, r;
	int col, mcnt, mbeg;
	svector *srhs = &(lp->srhs);
	svector *ssoln = &(lp->ssoln);

	for (i = 0; i < lp->nrows; i++)
	{
		EGlpNumZero (lp->xbz[i]);
		EGlpNumZero (srhs->coef[i]);
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
					EGlpNumSubTo (srhs->coef[lp->matind[mbeg + i]], lp->matval[mbeg + i]);
			else
				for (i = 0; i < mcnt; i++)
					EGlpNumAddTo (srhs->coef[lp->matind[mbeg + i]], lp->matval[mbeg + i]);
		}
	}
	for (i = 0, r = 0; i < lp->nrows; i++)
		if (EGlpNumIsNeqqZero (srhs->coef[i]))
		{
			EGlpNumCopy (srhs->coef[r], srhs->coef[i]);
			srhs->indx[r] = i;
			r++;
		}
	srhs->nzcnt = r;

	ILLbasis_column_solve (lp, srhs, ssoln);
	for (i = 0; i < ssoln->nzcnt; i++)
		EGlpNumCopy (lp->xbz[ssoln->indx[i]], ssoln->coef[i]);
}

void ILLfct_compute_phaseI_piz (
	lpinfo * lp)
{
	int i, r;
	svector *srhs = &(lp->srhs);
	svector *ssoln = &(lp->ssoln);

	for (i = 0, r = 0; i < lp->nrows; i++)
	{
		EGlpNumZero (lp->pIpiz[i]);
		if (lp->bfeas[i] != 0)
		{
			srhs->indx[r] = i;
			EGlpNumSet (srhs->coef[r], (double) lp->bfeas[i]);
			r++;
		}
	}
	srhs->nzcnt = r;

	ILLbasis_row_solve (lp, srhs, ssoln);
	for (i = 0; i < ssoln->nzcnt; i++)
		EGlpNumCopy (lp->pIpiz[ssoln->indx[i]], ssoln->coef[i]);
	ILLfct_update_counts (lp, CNT_P1PINZ, ssoln->nzcnt, zeroLpNum);
}

void ILLfct_compute_phaseI_dz (
	lpinfo * lp)
{
	int i, j;
	int col;
	int mcnt, mbeg;
	EGlpNum_t sum;

	EGlpNumInitVar (sum);
	ILL_IFTRACE ("%s\n", __func__);

	for (j = 0; j < lp->nnbasic; j++)
	{
		EGlpNumZero (sum);
		col = lp->nbaz[j];
		mcnt = lp->matcnt[col];
		mbeg = lp->matbeg[col];
		for (i = 0; i < mcnt; i++)
			EGlpNumAddInnProdTo (sum, lp->pIpiz[lp->matind[mbeg + i]],
													 lp->matval[mbeg + i]);
		EGlpNumCopyNeg (lp->pIdz[j], sum);
		ILL_IFTRACE ("%d:%d:%lf:%la\n", j, col, EGlpNumToLf (sum),
								 EGlpNumToLf (sum));
	}
	EGlpNumClearVar (sum);
}

void ILLfct_compute_yz (
	lpinfo * lp,
	svector * yz,
	svector * updz,
	int col)
{
	svector a;

	a.nzcnt = lp->matcnt[col];
	a.indx = &(lp->matind[lp->matbeg[col]]);
	a.coef = &(lp->matval[lp->matbeg[col]]);

	ILLfactor_set_factor_dparam (lp->f, QS_FACTOR_SZERO_TOL, PIVZ_TOLER);
	if (updz)
		ILLbasis_column_solve_update (lp, &a, updz, yz);
	else
		ILLbasis_column_solve (lp, &a, yz);
	ILLfactor_set_factor_dparam (lp->f, QS_FACTOR_SZERO_TOL, SZERO_TOLER);
}

void ILLfct_compute_zz (
	lpinfo * lp,
	svector * zz,
	int row)
{
	ILLfct_compute_binvrow (lp, zz, row, PIVZ_TOLER);
}

void ILLfct_compute_binvrow (
	lpinfo * lp,
	svector * zz,
	int row,
	EGlpNum_t ztoler)
{
	svector a;
	EGlpNum_t e;

	EGlpNumInitVar (e);
	EGlpNumOne (e);

	a.nzcnt = 1;
	a.coef = &e;
	a.indx = &row;

	if (EGlpNumIsGreatZero (ztoler))
		ILLfactor_set_factor_dparam (lp->f, QS_FACTOR_SZERO_TOL, ztoler);
	ILLbasis_row_solve (lp, &a, zz);
	if (EGlpNumIsGreatZero (ztoler))
		ILLfactor_set_factor_dparam (lp->f, QS_FACTOR_SZERO_TOL, SZERO_TOLER);
	EGlpNumClearVar (e);
}

void ILLfct_compute_psteep_upv (
	lpinfo * lp,
	svector * swz)
{
	ILLbasis_row_solve (lp, &(lp->yjz), swz);
}

void ILLfct_compute_dsteep_upv (
	lpinfo * lp,
	svector * swz)
{
	ILLbasis_column_solve (lp, &(lp->zz), swz);
}

static int compute_zA1 (
	lpinfo * lp,
	svector * z,
	svector * zA,
	EGlpNum_t ztoler)
{
	int rval = 0;
	int i, j, nz = 0;
	int col, mcnt, mbeg;
	EGlpNum_t sum;
	EGlpNum_t *v = 0;

	EGlpNumInitVar (sum);
	v = EGlpNumAllocArray (lp->nrows);

	for (i = 0; i < lp->nrows; i++)
		EGlpNumZero (v[i]);
	for (i = 0; i < z->nzcnt; i++)
		EGlpNumCopy (v[z->indx[i]], z->coef[i]);

	for (j = 0; j < lp->nnbasic; j++)
	{
		EGlpNumZero (sum);
		col = lp->nbaz[j];
		mcnt = lp->matcnt[col];
		mbeg = lp->matbeg[col];
		for (i = 0; i < mcnt; i++)
			EGlpNumAddInnProdTo (sum, v[lp->matind[mbeg + i]], lp->matval[mbeg + i]);

		if (EGlpNumIsNeqZero (sum, ztoler))
		{
			EGlpNumCopy (zA->coef[nz], sum);
			zA->indx[nz] = j;
			nz++;
		}
	}
	zA->nzcnt = nz;

	EGlpNumClearVar (sum);
	EGlpNumFreeArray (v);
	EG_RETURN (rval);
}


static int compute_zA3 (
	lpinfo * lp,
	svector * z,
	svector * zA,
	EGlpNum_t ztoler)
{
	int rval = 0;
	int i, j, k, ix;
	int nz = 0;
	int row, col;
	int rcnt, rbeg;
	EGlpNum_t val;

	EGlpNumInitVar (val);
	k = 0;
	for (i = 0; i < z->nzcnt; i++)
	{
		row = z->indx[i];
		EGlpNumCopy (val, z->coef[i]);
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
				EGlpNumAddInnProdTo (lp->work.coef[ix], val, lp->rowval[rbeg + j]);
			}
		}
	}
	for (j = 0; j < k; j++)
	{
		ix = lp->work.indx[j];
		EGlpNumCopy (val, lp->work.coef[ix]);
		EGlpNumZero (lp->work.coef[ix]);
		lp->iwork[ix] = 0;
		if (EGlpNumIsNeqZero (val, ztoler))
		{
			EGlpNumCopy (zA->coef[nz], val);
			zA->indx[nz] = ix;
			nz++;
		}
	}
	zA->nzcnt = nz;
	EGlpNumClearVar (val);
	EG_RETURN (rval);
}

int ILLfct_compute_zA (
	lpinfo * lp,
	svector * z,
	svector * zA)
{
	if (z->nzcnt < lp->nrows / 2)
		return compute_zA3 (lp, z, zA, PIVZ_TOLER);
	else
		return compute_zA1 (lp, z, zA, PIVZ_TOLER);
}

/* compute v^T A */
void ILLfct_compute_vA (
	lpinfo * lp,
	svector * v,
	EGlpNum_t * vA)
{
	int i, j;
	int row, col;
	int rcnt, rbeg;
	EGlpNum_t val;

	EGlpNumInitVar (val);

	for (j = 0; j < lp->ncols; j++)
		EGlpNumZero (vA[j]);

	for (i = 0; i < v->nzcnt; i++)
	{
		row = v->indx[i];
		EGlpNumCopy (val, v->coef[i]);
		rcnt = lp->rowcnt[row];
		rbeg = lp->rowbeg[row];
		for (j = 0; j < rcnt; j++)
		{
			col = lp->rowind[rbeg + j];
			EGlpNumAddInnProdTo (vA[col], val, lp->rowval[rbeg + j]);
		}
	}

	for (j = 0; j < lp->ncols; j++)
		if (!EGlpNumIsNeqZero (vA[j], SZERO_TOLER))
			EGlpNumZero (vA[j]);

	EGlpNumClearVar (val);
	return;
}

/* update information */

/*
1) lvstat - new status of leaving var.
*/
void ILLfct_update_basis_info (
	lpinfo * lp,
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

void ILLfct_update_xz (
	lpinfo * lp,
	EGlpNum_t tz,
	int eindex,
	int lindex)
{
	int i, evar, estat;

	ILL_IFTRACE ("%s:%la:%d:%d:%d\n", __func__, EGlpNumToLf (tz), eindex,
							 lindex, lp->yjz.nzcnt);

	if (EGlpNumIsNeqqZero (tz))
		for (i = 0; i < lp->yjz.nzcnt; i++)
			EGlpNumSubInnProdTo (lp->xbz[lp->yjz.indx[i]], tz, lp->yjz.coef[i]);

	if (lindex >= 0)
	{															/* variable leaves basis */
		evar = lp->nbaz[eindex];
		estat = lp->vstat[evar];
		if (estat == STAT_LOWER)
			EGlpNumCopySum (lp->xbz[lindex], lp->lz[evar], tz);
		else if (estat == STAT_UPPER)
			EGlpNumCopySum (lp->xbz[lindex], lp->uz[evar], tz);
		else if (estat == STAT_ZERO)
			EGlpNumCopy (lp->xbz[lindex], tz);
	}
}

void ILLfct_update_piz (
	lpinfo * lp,
	EGlpNum_t alpha)
{
	int i;

	for (i = 0; i < lp->zz.nzcnt; i++)
		EGlpNumAddInnProdTo (lp->piz[lp->zz.indx[i]], alpha, lp->zz.coef[i]);
}

void ILLfct_update_pIpiz (
	lpinfo * lp,
	svector * z,
	const EGlpNum_t alpha)
{
	int i;

	if (!EGlpNumIsNeqqZero (alpha))
		return;
	if (EGlpNumIsEqqual (alpha, oneLpNum))
	{
		for (i = 0; i < z->nzcnt; i++)
			EGlpNumAddTo (lp->pIpiz[z->indx[i]], z->coef[i]);
	}
	else
	{
		for (i = 0; i < z->nzcnt; i++)
			EGlpNumAddInnProdTo (lp->pIpiz[z->indx[i]], alpha, z->coef[i]);
	}
}

void ILLfct_update_dz (
	lpinfo * lp,
	int eindex,
	EGlpNum_t alpha)
{
	int i;

	for (i = 0; i < lp->zA.nzcnt; i++)
		EGlpNumSubInnProdTo (lp->dz[lp->zA.indx[i]], alpha, lp->zA.coef[i]);
	EGlpNumCopyNeg (lp->dz[eindex], alpha);
}

void ILLfct_update_pIdz (
	lpinfo * lp,
	svector * zA,
	int eindex,
	const EGlpNum_t alpha)
{
	int i;

	if (!EGlpNumIsNeqqZero (alpha))
		return;

	if (EGlpNumIsEqqual (alpha, oneLpNum))
	{
		for (i = 0; i < zA->nzcnt; i++)
			EGlpNumSubTo (lp->pIdz[zA->indx[i]], zA->coef[i]);
	}
	else
	{
		for (i = 0; i < zA->nzcnt; i++)
			EGlpNumSubInnProdTo (lp->pIdz[zA->indx[i]], alpha, zA->coef[i]);
	}
	if (eindex > -1)
		EGlpNumCopyNeg (lp->pIdz[eindex], alpha);
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
	lpinfo * lp,
	EGlpNum_t ftol,
	int *chgb)
{
	int rval = 0;
	int i, col, nchg = 0;
	EGlpNum_t newb, cftol;
	EGlpNum_t *x, *l, *u;
	ILLrandstate r;

	EGlpNumInitVar (newb);
	EGlpNumInitVar (cftol);
	EGlpNumCopyAbs (cftol, ftol);
	EGlpNumDivUiTo (cftol, 10);

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
		EGlpNumCopyDiff (newb, *x, ftol);
		if (EGlpNumIsNeqq (*l, NINFTY) && EGlpNumIsLess (newb, *l))
		{
			EGlpNumSet (newb, -1.0 * (my_rand (50, &(lp->rstate)) + 1.0));
			EGlpNumMultTo (newb, cftol);
			if (EGlpNumIsLess (*x, *l))
				EGlpNumAddTo (newb, *x);
			else
				EGlpNumAddTo (newb, *l);
			rval = ILLfct_bound_shift (lp, col, BOUND_LOWER, newb);
			CHECKRVALG (rval, CLEANUP);
			nchg++;
		}
		EGlpNumCopySum (newb, *x, ftol);
		if (EGlpNumIsNeqq (*u, INFTY) && EGlpNumIsLess (*u, newb))
		{
			EGlpNumSet (newb, my_rand (50, &(lp->rstate)) + 1.0);
			EGlpNumMultTo (newb, cftol);
			if (EGlpNumIsLess (*x, *u))
				EGlpNumAddTo (newb, *u);
			else
				EGlpNumAddTo (newb, *x);
			rval = ILLfct_bound_shift (lp, col, BOUND_UPPER, newb);
			CHECKRVALG (rval, CLEANUP);
			nchg++;
		}
	}
	*chgb = nchg;

CLEANUP:
	EGlpNumClearVar (newb);
	EGlpNumClearVar (cftol);
	EG_RETURN (rval);
}

static int expand_phaseI_bounds (
	lpinfo * lp,
	int *chgb)
{
	int rval = 0;
	int i, col, nchg = 0;
	EGlpNum_t newb, cftol;
	EGlpNum_t *u, *l, *x;
	ILLrandstate r;

	EGlpNumInitVar (newb);
	EGlpNumInitVar (cftol);
	EGlpNumCopyAbs (cftol, lp->tol->ip_tol);
	EGlpNumDivUiTo (cftol, 10);
	ILLutil_sprand (1, &r);

	for (i = 0; i < lp->nrows; i++)
	{
		col = lp->baz[i];
		if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFREE)
			continue;
		x = &(lp->xbz[i]);
		l = &(lp->lz[col]);
		u = &(lp->uz[col]);

		if (EGlpNumIsNeqq (*l, NINFTY) && EGlpNumIsEqual (*x, *l, cftol))
		{
			EGlpNumSet (newb, my_rand (50, &(lp->rstate)) + 1.0);
			EGlpNumMultTo (newb, cftol);
			EGlpNumSign (newb);
			EGlpNumAddTo (newb, *l);
			rval = ILLfct_bound_shift (lp, col, BOUND_LOWER, newb);
			CHECKRVALG (rval, CLEANUP);
			nchg++;
		}
		if (EGlpNumIsNeqq (*u, INFTY) && EGlpNumIsEqual (*x, *u, cftol))
		{
			EGlpNumSet (newb, my_rand (50, &(lp->rstate)) + 1.0);
			EGlpNumMultTo (newb, cftol);
			EGlpNumAddTo (newb, *u);
			rval = ILLfct_bound_shift (lp, col, BOUND_UPPER, newb);
			CHECKRVALG (rval, CLEANUP);
			nchg++;
		}
	}
	*chgb = nchg;

CLEANUP:
	EGlpNumClearVar (newb);
	EGlpNumClearVar (cftol);
	EG_RETURN (rval);
}

int ILLfct_adjust_viol_bounds (
	lpinfo * lp)
{
	int rval = 0;
	int chgb = 0;
	EGlpNum_t tol;

	EGlpNumInitVar (tol);
	EGlpNumCopyNeg (tol, lp->tol->pfeas_tol);
	rval = expand_var_bounds (lp, tol, &chgb);
#if FCT_DEBUG > 0
	if (rval == 0)
		printf ("adjusting %d bounds\n", chgb);
#endif
	EGlpNumClearVar (tol);
	EG_RETURN (rval);
}

int ILLfct_perturb_bounds (
	lpinfo * lp)
{
	int rval = 0;
	int chgb = 0;

	rval = expand_var_bounds (lp, lp->tol->ip_tol, &chgb);
#if FCT_DEBUG > 0
	if (rval == 0)
		printf ("perturbing %d bounds\n", chgb);
#endif
	EG_RETURN (rval);
}

int ILLfct_perturb_phaseI_bounds (
	lpinfo * lp)
{
	int rval = 0;
	int chgb = 0;

	rval = expand_phaseI_bounds (lp, &chgb);
#if FCT_DEBUG > 0
	if (rval == 0)
		printf ("perturbing %d phase I bounds\n", chgb);
#endif
	EG_RETURN (rval);
}

int ILLfct_bound_shift (
	lpinfo * lp,
	int col,
	int bndtype,
	EGlpNum_t newbnd)
{
	int rval = 0;
	bndinfo *nbnd = 0;

	ILL_IFTRACE ("\n%s:%d:%d:%la", __func__, col, bndtype, EGlpNumToLf (newbnd));
	nbnd = ILLfct_new_bndinfo ();

	nbnd->varnum = col;
	nbnd->btype = bndtype;
	if (bndtype == BOUND_LOWER)
	{
		EGlpNumCopy (nbnd->pbound, lp->lz[col]);
		EGlpNumCopy (nbnd->cbound, newbnd);
		EGlpNumCopy (lp->lz[col], newbnd);
	}
	else
	{
		EGlpNumCopy (nbnd->pbound, lp->uz[col]);
		EGlpNumCopy (nbnd->cbound, newbnd);
		EGlpNumCopy (lp->uz[col], newbnd);
	}
	ILL_IFTRACE (":%la", EGlpNumToLf (nbnd->pbound));
	if (lp->vtype[col] == VFIXED || lp->vtype[col] == VARTIFICIAL)
	{
		/* printf ("changing f/a bound\n"); */
		if (EGlpNumIsLess (lp->lz[col], lp->uz[col]))
			lp->vtype[col] = VBOUNDED;
	}

	nbnd->next = lp->bchanges;
	lp->bchanges = nbnd;
	lp->nbchange++;

//CLEANUP:
	if (rval)
		ILLfct_free_bndinfo (nbnd);
	ILL_IFTRACE ("\n");
	EG_RETURN (rval);
}

void ILLfct_unroll_bound_change (
	lpinfo * lp)
{
	int col;
	int changex = 0;
	bndinfo *bptr = lp->bchanges;
	bndinfo *nptr = 0;

	ILL_IFTRACE ("%s:", __func__);

	while (lp->nbchange != 0)
	{
		col = bptr->varnum;
		ILL_IFTRACE (":%d", col);

		if (bptr->btype == BOUND_UPPER)
			EGlpNumCopy (lp->uz[col], bptr->pbound);
		else
			EGlpNumCopy (lp->lz[col], bptr->pbound);

		if (lp->vtype[col] == VBOUNDED)
		{
			if (EGlpNumIsEqqual (lp->lz[col], lp->uz[col]))
				lp->vtype[col] = (!EGlpNumIsNeqqZero (lp->lz[col])) ?
					VARTIFICIAL : VFIXED;
		}

		if (lp->vstat[col] != STAT_BASIC)
		{
			if ((bptr->btype == BOUND_UPPER && lp->vstat[col] == STAT_UPPER) ||
					(bptr->btype == BOUND_LOWER && lp->vstat[col] == STAT_LOWER))
				changex++;
		}
		nptr = bptr->next;
		EGlpNumClearVar ((bptr->cbound));
		EGlpNumClearVar ((bptr->pbound));
		ILL_IFFREE (bptr, bndinfo);
		bptr = nptr;
		lp->nbchange--;
	}
	lp->bchanges = bptr;
	ILL_IFTRACE ("\n");
	if (changex)
		ILLfct_compute_xbz (lp);
}

static int expand_var_coefs (
	lpinfo * lp,
	EGlpNum_t ftol,
	int *chgc)
{
	int rval = 0;
	int i, col, vs, vt;
	int nchg = 0;
	EGlpNum_t newc, cftol, mftol[1];
	EGlpNum_t *c, *dj;
	ILLrandstate r;

	EGlpNumInitVar (newc);
	EGlpNumInitVar (cftol);
	EGlpNumInitVar (mftol[0]);
	EGlpNumCopyAbs (cftol, ftol);
	EGlpNumDivUiTo (cftol, 10);
	EGlpNumCopyNeg (mftol[0], ftol);
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
			EGlpNumCopyDiff (newc, *c, *dj);
			rval = ILLfct_coef_shift (lp, col, newc);
			CHECKRVALG (rval, CLEANUP);
			nchg++;
			break;
		case STAT_LOWER:
			if (EGlpNumIsLess (*dj, ftol))
			{
				EGlpNumSet (newc, my_rand (50, &(lp->rstate)) + 1.0);
				EGlpNumMultTo (newc, cftol);
				EGlpNumAddTo (newc, *c);
				if (EGlpNumIsLessZero (*dj))
					EGlpNumSubTo (newc, *dj);
				rval = ILLfct_coef_shift (lp, col, newc);
				CHECKRVALG (rval, CLEANUP);
				nchg++;
			}
			break;
		case STAT_UPPER:
			if (EGlpNumIsLess (mftol[0], *dj))
			{
				EGlpNumSet (newc, my_rand (50, &(lp->rstate)) + 1.0);
				EGlpNumMultTo (newc, cftol);
				EGlpNumSign (newc);
				EGlpNumAddTo (newc, *c);
				if (EGlpNumIsGreatZero (*dj))
					EGlpNumSubTo (newc, *dj);
				rval = ILLfct_coef_shift (lp, col, newc);
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
	EGlpNumClearVar (mftol[0]);
	EGlpNumClearVar (newc);
	EGlpNumClearVar (cftol);
	EG_RETURN (rval);
}

int ILLfct_adjust_viol_coefs (
	lpinfo * lp)
{
	int rval = 0;
	int chgc = 0;
	EGlpNum_t tol;

	EGlpNumInitVar (tol);
	EGlpNumCopyNeg (tol, lp->tol->dfeas_tol);

	rval = expand_var_coefs (lp, tol, &chgc);
#if FCT_DEBUG > 0
	if (rval == 0)
		printf ("perturbing %d coefs\n", chgc);
#endif
	EGlpNumClearVar (tol);
	EG_RETURN (rval);
}

int ILLfct_perturb_coefs (
	lpinfo * lp)
{
	int rval = 0;
	int chgc = 0;

	rval = expand_var_coefs (lp, lp->tol->id_tol, &chgc);
#if FCT_DEBUG > 0
	if (rval == 0)
		printf ("perturbing %d coefs\n", chgc);
#endif
	EG_RETURN (rval);
}

int ILLfct_coef_shift (
	lpinfo * lp,
	int col,
	EGlpNum_t newcoef)
{
	int rval = 0;
	coefinfo *ncoef = 0;

	ILL_SAFE_MALLOC (ncoef, 1, coefinfo);
	EGlpNumInitVar ((ncoef->pcoef));
	EGlpNumInitVar ((ncoef->ccoef));

	ncoef->varnum = col;
	EGlpNumCopy (ncoef->pcoef, lp->cz[col]);
	EGlpNumCopy (ncoef->ccoef, newcoef);
	EGlpNumCopy (lp->cz[col], newcoef);
	ncoef->next = lp->cchanges;
	lp->cchanges = ncoef;
	EGlpNumAddTo (lp->dz[lp->vindex[col]], ncoef->ccoef);
	EGlpNumSubTo (lp->dz[lp->vindex[col]], ncoef->pcoef);
	lp->ncchange++;

CLEANUP:
	if (rval)
	{
		EGlpNumClearVar ((ncoef->pcoef));
		EGlpNumClearVar ((ncoef->ccoef));
		ILL_IFFREE (ncoef, coefinfo);
	}
	EG_RETURN (rval);
}

void ILLfct_unroll_coef_change (
	lpinfo * lp)
{
	int bascoef = 0;
	coefinfo *cptr = (coefinfo *) lp->cchanges;
	coefinfo *nptr = 0;

	while (lp->ncchange != 0)
	{
		EGlpNumCopy (lp->cz[cptr->varnum], cptr->pcoef);
		if (lp->vstat[cptr->varnum] != STAT_BASIC)
		{
			EGlpNumAddTo (lp->dz[lp->vindex[cptr->varnum]], cptr->pcoef);
			EGlpNumSubTo (lp->dz[lp->vindex[cptr->varnum]], cptr->ccoef);
		}
		else
			bascoef++;

		nptr = cptr->next;
		EGlpNumClearVar ((cptr->pcoef));
		EGlpNumClearVar ((cptr->ccoef));
		ILL_IFFREE (cptr, coefinfo);
		cptr = nptr;
		lp->ncchange--;
	}
	lp->cchanges = cptr;
	if (bascoef)
	{
		ILLfct_compute_piz (lp);
		ILLfct_compute_dz (lp);
	}
}

/* feasibility routines */
void ILLfct_check_pfeasible (
	lpinfo * lp,
	feas_info * fs,
	const EGlpNum_t ftol)
{
	int i, col;
	EGlpNum_t infeas, err1, err2;

	EGlpNumInitVar (infeas);
	EGlpNumInitVar (err1);
	EGlpNumInitVar (err2);
	EGlpNumZero (infeas);
	fs->pstatus = PRIMAL_FEASIBLE;
	EGlpNumZero (fs->totinfeas);
	ILL_IFTRACE ("%s:tol %la\n", __func__, EGlpNumToLf (ftol));

	for (i = 0; i < lp->nrows; i++)
	{
		col = lp->baz[i];
		EGlpNumCopyDiff (err1, lp->xbz[i], lp->uz[col]);
		EGlpNumCopyDiff (err2, lp->lz[col], lp->xbz[i]);
		if (EGlpNumIsLess (ftol, err1)
				&& EGlpNumIsNeq (lp->uz[col], INFTY, oneLpNum))
		{
			EGlpNumAddTo (infeas, err1);
			WARNINGL (QSE_WLVL, EGlpNumIsLess (INFTY, err1),
							 "This is imposible lu = %15lg xbz = %15lg" " INFTY = %15lg",
							 EGlpNumToLf (lp->uz[col]), EGlpNumToLf (lp->xbz[i]),
							 EGlpNumToLf (INFTY));
			lp->bfeas[i] = 1;
		}
		else if (EGlpNumIsLess (ftol, err2)
						 && EGlpNumIsNeq (lp->lz[col], NINFTY, oneLpNum))
		{
			EGlpNumAddTo (infeas, err2);
			WARNINGL (QSE_WLVL, EGlpNumIsLess (INFTY, err2),
							 "This is imposible lz = %15lg xbz = %15lg" " NINFTY = %15lg",
							 EGlpNumToLf (lp->lz[col]), EGlpNumToLf (lp->xbz[i]),
							 EGlpNumToLf (NINFTY));
			lp->bfeas[i] = -1;
		}
		else
			lp->bfeas[i] = 0;
	}
	if (EGlpNumIsNeqqZero (infeas))
	{
		fs->pstatus = PRIMAL_INFEASIBLE;
		EGlpNumCopy (fs->totinfeas, infeas);
		ILL_IFTRACE ("%s:inf %la\n", __func__, EGlpNumToLf (infeas));
		if (EGlpNumIsLessZero (fs->totinfeas))
		{
			printf ("Negative infeasibility, Imposible! %lf %la\n",
							EGlpNumToLf (infeas), EGlpNumToLf (infeas));
		}
	}
	EGlpNumCopy (lp->pinfeas, infeas);
	EGlpNumClearVar (infeas);
	EGlpNumClearVar (err1);
	EGlpNumClearVar (err2);
}

/* feasibility routines */
void ILLfct_check_pIpfeasible (
	lpinfo * lp,
	feas_info * fs,
	EGlpNum_t ftol)
{
	int i, col;
	int ninf = 0;

	fs->pstatus = PRIMAL_FEASIBLE;
	EGlpNumZero (fs->totinfeas);

	for (i = 0; i < lp->nrows; i++)
	{
		if (!EGlpNumIsNeqZero (lp->xbz[i], ftol))
			continue;
		col = lp->baz[i];
		if (EGlpNumIsGreatZero(lp->xbz[i]) &&
				EGlpNumIsNeqq (lp->uz[col], INFTY))
		{
			ninf++;
		}
		else if (EGlpNumIsLessZero (lp->xbz[i]) &&
						 EGlpNumIsNeqq (lp->lz[col], NINFTY))
		{
			ninf++;
		}
	}
	if (ninf != 0)
		fs->pstatus = PRIMAL_INFEASIBLE;
}

void ILLfct_check_dfeasible (
	lpinfo * lp,
	feas_info * fs,
	const EGlpNum_t ftol)
{
	int j, col;
	EGlpNum_t infeas;

	EGlpNumInitVar (infeas);
	EGlpNumZero (infeas);
	fs->dstatus = DUAL_FEASIBLE;
	EGlpNumZero (fs->totinfeas);

	for (j = 0; j < lp->nnbasic; j++)
	{
		lp->dfeas[j] = 0;
		if (!EGlpNumIsNeqZero (lp->dz[j], ftol))
			continue;
		col = lp->nbaz[j];

		if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED)
			continue;

		if (EGlpNumIsLessZero (lp->dz[j]) &&
				(lp->vstat[col] == STAT_LOWER || lp->vstat[col] == STAT_ZERO))
		{
			EGlpNumSubTo (infeas, lp->dz[j]);
			lp->dfeas[j] = -1;
		}
		else if (EGlpNumIsGreatZero (lp->dz[j]) &&
						 (lp->vstat[col] == STAT_UPPER || lp->vstat[col] == STAT_ZERO))
		{
			EGlpNumAddTo (infeas, lp->dz[j]);
			lp->dfeas[j] = 1;
		}
	}

	if (EGlpNumIsNeqqZero (infeas))
	{
		EGlpNumCopy (fs->totinfeas, infeas);
		fs->dstatus = DUAL_INFEASIBLE;
		ILL_IFTRACE ("%s:inf %la\n", __func__, EGlpNumToLf (infeas));
		if (EGlpNumIsLessZero (fs->totinfeas))
		{
			printf ("Negative infeasibility, Imposible! %lf %la\n",
							EGlpNumToLf (infeas), EGlpNumToLf (infeas));
		}
	}
	EGlpNumCopy (lp->dinfeas, infeas);
	EGlpNumClearVar (infeas);
}

void ILLfct_check_pIdfeasible (
	lpinfo * lp,
	feas_info * fs,
	EGlpNum_t ftol)
{
	int j, col;
	int ninf = 0;
	EGlpNum_t *dz = lp->pIdz;

	fs->dstatus = DUAL_FEASIBLE;

	for (j = 0; j < lp->nnbasic; j++)
	{
		if (!EGlpNumIsNeqZero (dz[j], ftol))
			continue;
		col = lp->nbaz[j];
		if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED)
			continue;

		if (EGlpNumIsLessZero (dz[j]) &&
				(lp->vstat[col] == STAT_LOWER || lp->vstat[col] == STAT_ZERO))
			ninf++;
		else if (EGlpNumIsGreatZero (dz[j]) &&
						 (lp->vstat[col] == STAT_UPPER || lp->vstat[col] == STAT_ZERO))
			ninf++;
	}

	if (ninf != 0)
		fs->dstatus = DUAL_INFEASIBLE;
}

void ILLfct_dual_adjust (
	lpinfo * lp,
	const EGlpNum_t ftol)
{
	int j, col;

	for (j = 0; j < lp->nnbasic; j++)
	{
		if (!EGlpNumIsNeqZero (lp->dz[j], ftol))
			continue;
		col = lp->nbaz[j];
		if (EGlpNumIsLessZero (lp->dz[j]) &&
				EGlpNumIsNeqq (lp->uz[col], INFTY))
			lp->vstat[col] = STAT_UPPER;
		else if (EGlpNumIsGreatZero (lp->dz[j]) &&
						 EGlpNumIsNeqq (lp->lz[col], NINFTY))
			lp->vstat[col] = STAT_LOWER;
	}
}

void ILLfct_dphaseI_simple_update (
	lpinfo * lp,
	EGlpNum_t ftol)
{
	int j, col;

	for (j = 0; j < lp->nnbasic; j++)
	{
		if (!EGlpNumIsNeqZero (lp->dz[j], ftol))
			continue;
		col = lp->nbaz[j];
		if (EGlpNumIsLessZero (lp->dz[j] ) && lp->vtype[col] == VBOUNDED)
			lp->vstat[col] = STAT_UPPER;
		else if (EGlpNumIsGreatZero (lp->dz[j]) && lp->vtype[col] == VBOUNDED)
			lp->vstat[col] = STAT_LOWER;
	}
}

/* set status values */
void ILLfct_set_status_values (
	lpinfo * lp,
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

void ILLfct_init_counts (
	lpinfo * lp)
{
	int i;
	count_struct *c = lp->cnts;

#define C_VALUE(a) (1.0+(double)(a)/(PARAM_HEAP_RATIO*ILLutil_our_log2(a)))
	EGlpNumSet (c->y_ravg, C_VALUE (lp->nrows));
	EGlpNumSet (c->za_ravg, C_VALUE (lp->nnbasic));
	ILL_IFTRACE ("%s:%la\n", __func__, EGlpNumToLf (c->za_ravg));
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
	count_struct * c,
	int phase,
	const EGlpNum_t piv2)
{
	int i = 0;
	EGlpNum_t v, piv;

	if (!EGlpNumIsNeqqZero(piv2))
		return;
	EGlpNumInitVar (v);
	EGlpNumInitVar (piv);
	EGlpNumCopyAbs (piv, piv2);
	EGlpNumOne (v);
	while (EGlpNumIsLess (piv, v) && i < 9)
	{
		EGlpNumDivUiTo (v, 10);
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
	EGlpNumClearVar (v);
	EGlpNumClearVar (piv);
}

void ILLfct_update_counts (
	lpinfo * lp,
	int f,
	int upi,
	const EGlpNum_t upd)
{
	count_struct *c = lp->cnts;

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
		EGlpNumMultUiTo (c->y_ravg, c->tot_iter);
		EGlpNumAddUiTo (c->y_ravg, upi);
		EGlpNumDivUiTo (c->y_ravg, c->tot_iter + 1);
		break;
	case CNT_ZARAVG:
		ILL_IFTRACE ("%s:%d:%d:%d:%la:%la", __func__, f, c->tot_iter, upi,
								 EGlpNumToLf (upd), EGlpNumToLf (c->za_ravg));
		EGlpNumMultUiTo (c->za_ravg, c->tot_iter);
		EGlpNumAddUiTo (c->za_ravg, upi);
		EGlpNumDivUiTo (c->za_ravg, c->tot_iter + 1);
		ILL_IFTRACE (":%la\n", EGlpNumToLf (c->za_ravg));
		break;
	}
}

void ILLfct_print_counts (
	lpinfo * lp)
{
	int i, niter;
	count_struct *c = lp->cnts;

	c->tot_iter = c->pI_iter + c->pII_iter + c->dI_iter + c->dII_iter;
	niter = (c->tot_iter == 0) ? 1 : c->tot_iter;
	printf ("Counts for problem %s\n", lp->O->probname);
	if (c->num_y != 0)
		printf ("avg ynz = %.2f\n", (double) c->ynz_cnt / c->num_y);
	if (c->num_z != 0)
		printf ("avg znz = %.2f\n", (double) c->znz_cnt / c->num_z);
	if (c->num_za != 0)
		printf ("avg zanz = %.2f\n", (double) c->zanz_cnt / c->num_za);
	printf ("avg pnorm = %.2f\n", (double) c->pnorm_cnt / lp->nnbasic);
	printf ("avg dnorm = %.2f\n", (double) c->dnorm_cnt / lp->nrows);
	if (c->num_pi != 0)
		printf ("avg pinz = %.2f\n", (double) c->pinz_cnt / c->num_pi);
	if (c->num_pi1 != 0)
		printf ("avg piInz = %.2f\n", (double) c->pi1nz_cnt / c->num_pi1);
	if (c->num_up != 0)
		printf ("avg upnz = %.2f\n", (double) c->upnz_cnt / c->num_up);

	for (i = 0; i < 10; i++)
		printf ("piv 1.0e-%d : %d %d %d %d\n",
						i, c->pivpI[i], c->pivpII[i], c->pivdI[i], c->pivdII[i]);
}


/* c <- a + t*b */
static void add_vectors (
	lpinfo * lp,
	svector * a,
	svector * b,
	svector * c,
	const EGlpNum_t t)
{
	int i, r, l;
	svector *w = &(lp->work);

	for (i = 0; i < b->nzcnt; i++)
	{
		r = b->indx[i];
		w->indx[i] = r;
		EGlpNumCopy (w->coef[r], t);
		EGlpNumMultTo (w->coef[r], b->coef[i]);
		lp->iwork[r] = 1;
	}
	l = b->nzcnt;

	for (i = 0; i < a->nzcnt; i++)
	{
		r = a->indx[i];
		if (lp->iwork[r] == 0)
			w->indx[l++] = r;
		EGlpNumAddTo (w->coef[r], a->coef[i]);
	}
	for (i = 0; i < l; i++)
	{
		r = w->indx[i];
		c->indx[i] = r;
		EGlpNumCopy (c->coef[i], w->coef[r]);
		EGlpNumZero (w->coef[r]);
		lp->iwork[r] = 0;
	}
	w->nzcnt = 0;
	c->nzcnt = l;
}

void ILLfct_update_pfeas (
	lpinfo * lp,
	int lindex,
	svector * srhs)
{
	int i, k, r;
	int col, nz = 0;
	int cbnd, f;
	int *perm = lp->upd.perm;
	int *ix = lp->upd.ix;
	int tctr = lp->upd.tctr;
	EGlpNum_t *t = lp->upd.t;
	EGlpNum_t tz, *dty, ntmp;
	EGlpNum_t *l, *x, *u, *pftol = &(lp->tol->ip_tol);

	EGlpNumInitVar (tz);
	EGlpNumInitVar (ntmp);
	dty = &(lp->upd.dty);
	EGlpNumZero (*dty);
	EGlpNumCopyAbs (tz, lp->upd.tz);
	EGlpNumDivUiTo (tz, 100);
	EGlpNumAddTo (tz, lp->upd.tz);
	ILL_IFTRACE ("%s:%d", __func__, tctr);
	for (i = 0; i < tctr && EGlpNumIsLeq (t[perm[i]], tz); i++)
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
				EGlpNumCopyDiff (ntmp, *l, *x);
				if (EGlpNumIsNeqq (*l, NINFTY) && EGlpNumIsLess (*pftol, ntmp))
					f = -1;
				else
				{
					EGlpNumCopyDiff (ntmp, *x, *u);
					if (EGlpNumIsNeqq (*u, INFTY) && EGlpNumIsLess (*pftol, ntmp))
						f = 1;
				}

				ILL_IFTRACE (":%d:%d", f, lp->bfeas[r]);
				if (f != lp->bfeas[r])
				{
					srhs->indx[nz] = r;
					EGlpNumSet (srhs->coef[nz], (double) (f - lp->bfeas[r]));
					EGlpNumAddInnProdTo (*dty, srhs->coef[nz], lp->yjz.coef[k]);
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
	EGlpNumClearVar (tz);
	EGlpNumClearVar (ntmp);
}

void ILLfct_compute_ppIzz (
	lpinfo * lp,
	svector * srhs,
	svector * ssoln)
{
	if (srhs->nzcnt != 0)
	{
		ILL_IFTRACE ("%s:\n", __func__);
		ILLbasis_row_solve (lp, srhs, ssoln);
	}
}

void ILLfct_update_ppI_prices (
	lpinfo * lp,
	price_info * pinf,
	svector * srhs,
	svector * ssoln,
	int eindex,
	int lindex,
	const EGlpNum_t alpha)
{
	EGlpNum_t ntmp;

	EGlpNumInitVar (ntmp);
	EGlpNumCopy (ntmp, alpha);
	ILL_IFTRACE ("%s:\n", __func__);
	if (lindex == -1)
	{
		if (srhs->nzcnt != 0)
		{
			ILLfct_update_pIpiz (lp, ssoln, oneLpNum);
			if (pinf->p_strategy == COMPLETE_PRICING)
			{
				ILLfct_compute_zA (lp, ssoln, &(lp->zA));
				ILLfct_update_pIdz (lp, &(lp->zA), -1, oneLpNum);
			}
		}
		else
		{
			if (pinf->p_strategy == COMPLETE_PRICING)
				ILLprice_compute_dual_inf (lp, pinf, &eindex, 1, PRIMAL_PHASEI);
			else
				ILLprice_update_mpartial_price (lp, pinf, PRIMAL_PHASEI, COL_PRICING);
			EGlpNumClearVar (ntmp);
			return;
		}
	}
	else
	{
		if (srhs->nzcnt == 0)
		{
			ILLfct_update_pIpiz (lp, &(lp->zz), ntmp);
			if (pinf->p_strategy == COMPLETE_PRICING)
				ILLfct_update_pIdz (lp, &(lp->zA), eindex, ntmp);
		}
		else
		{
			EGlpNumCopyFrac (ntmp, lp->upd.dty, lp->upd.piv);
			EGlpNumSubTo (ntmp, alpha);
			EGlpNumSign (ntmp);
			add_vectors (lp, ssoln, &(lp->zz), &(lp->zz), ntmp);
			ILLfct_update_pIpiz (lp, &(lp->zz), oneLpNum);
			if (pinf->p_strategy == COMPLETE_PRICING)
			{
				ILLfct_compute_zA (lp, &(lp->zz), &(lp->zA));
				ILLfct_update_pIdz (lp, &(lp->zA), eindex, oneLpNum);
			}
		}
		EGlpNumSet (lp->pIdz[eindex], (double) (lp->upd.fs));
		EGlpNumAddTo (lp->pIdz[eindex], ntmp);
		EGlpNumSign (lp->pIdz[eindex]);
	}
	if (pinf->p_strategy == COMPLETE_PRICING)
	{
		ILLprice_compute_dual_inf (lp, pinf, lp->zA.indx, lp->zA.nzcnt,
															 PRIMAL_PHASEI);
		if (eindex > -1)
			ILLprice_compute_dual_inf (lp, pinf, &eindex, 1, PRIMAL_PHASEI);
		ILLfct_update_counts (lp, CNT_ZARAVG, lp->zA.nzcnt, zeroLpNum);
	}
	else
		ILLprice_update_mpartial_price (lp, pinf, PRIMAL_PHASEI, COL_PRICING);
	EGlpNumClearVar (ntmp);
	return;
}

void ILLfct_update_dfeas (
	lpinfo * lp,
	int eindex,
	svector * srhs)
{
	int i, j, k, c;
	int cbnd, col, nz = 0;
	int vs, vt, f;
	int delta;
	int *perm = lp->upd.perm;
	int *ix = lp->upd.ix;
	int tctr = lp->upd.tctr;
	int mcnt, mbeg;
	EGlpNum_t *t = lp->upd.t;
	EGlpNum_t *w = lp->work.coef;
	EGlpNum_t tz;
	EGlpNum_t *dty = &(lp->upd.dty);
	EGlpNum_t *dftol = &(lp->tol->id_tol);
	EGlpNum_t dj;

	EGlpNumInitVar (dj);
	EGlpNumInitVar (tz);
	EGlpNumZero (*dty);
	EGlpNumCopy (tz, lp->upd.tz);
	EGlpNumMultUiTo (tz, 101);
	EGlpNumDivUiTo (tz, 100);

	for (j = 0; j < tctr && EGlpNumIsLeq (t[perm[j]], tz); j++)
	{
		k = ix[perm[j]] / 10;
		c = lp->zA.indx[k];

		if (lp->iwork[c] != 1)
		{
			lp->iwork[c] = 1;
			cbnd = ix[perm[j]] % 10;
			col = lp->nbaz[c];
			EGlpNumCopy (dj, lp->dz[c]);
			vs = lp->vstat[col];
			vt = lp->vtype[col];

			if (cbnd == BSKIP)
			{
				if (!EGlpNumIsNeqZero (dj, *dftol));
				else if (EGlpNumIsLessZero (dj) && vs == STAT_LOWER)
					lp->vstat[col] = STAT_UPPER;
				else if (EGlpNumIsGreatZero (dj) && vs == STAT_UPPER)
					lp->vstat[col] = STAT_LOWER;
			}
			else if (c != eindex)
			{
				if (!EGlpNumIsNeqZero (dj, *dftol))
					f = 0;
				else if (EGlpNumIsLessZero (dj) &&
								 (vs == STAT_LOWER || vs == STAT_ZERO))
					f = -1;
				else if (EGlpNumIsGreatZero (dj) &&
								 (vs == STAT_UPPER || vs == STAT_ZERO))
					f = 1;
				else
					f = 0;

				if (f != lp->dfeas[c])
				{
					delta = f - lp->dfeas[c];
					mcnt = lp->matcnt[col];
					mbeg = lp->matbeg[col];
					EGlpNumSet (dj, (double) (delta));
					for (i = 0; i < mcnt; i++)
						EGlpNumAddInnProdTo (w[lp->matind[mbeg + i]], dj,
																 lp->matval[mbeg + i]);
					EGlpNumAddInnProdTo (*dty, dj, lp->zA.coef[k]);
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
			if (EGlpNumIsNeqqZero (w[i]))
			{
				EGlpNumCopy (srhs->coef[nz], w[i]);
				srhs->indx[nz] = i;
				nz++;
				EGlpNumZero (w[i]);
			}
	}

	srhs->nzcnt = nz;
	EGlpNumClearVar (dj);
	EGlpNumClearVar (tz);
}

void ILLfct_compute_dpIy (
	lpinfo * lp,
	svector * srhs,
	svector * ssoln)
{
	if (srhs->nzcnt != 0)
	{
		ILLbasis_column_solve (lp, srhs, ssoln);
	}
}

void ILLfct_update_dpI_prices (
	lpinfo * lp,
	price_info * pinf,
	svector * srhs,
	svector * ssoln,
	int lindex,
	EGlpNum_t alpha)
{
	int i;
	EGlpNum_t ntmp;

	EGlpNumInitVar (ntmp);
	EGlpNumZero (ntmp);

	if (srhs->nzcnt == 0)
	{
		ILLfct_update_xz (lp, alpha, -1, -1);
	}
	else
	{
		EGlpNumCopyFrac (ntmp, lp->upd.dty, lp->upd.piv);
		EGlpNumAddTo (ntmp, alpha);
		EGlpNumSign (ntmp);
		add_vectors (lp, ssoln, &(lp->yjz), &(lp->yjz), ntmp);
		EGlpNumSign (ntmp);
		for (i = 0; i < lp->yjz.nzcnt; i++)
			EGlpNumAddTo (lp->xbz[lp->yjz.indx[i]], lp->yjz.coef[i]);
	}
	EGlpNumSet (lp->xbz[lindex], ((double) (-lp->upd.fs)));
	EGlpNumAddTo (lp->xbz[lindex], ntmp);

	if (pinf->d_strategy == COMPLETE_PRICING)
	{
		ILLprice_compute_primal_inf (lp, pinf, lp->yjz.indx, lp->yjz.nzcnt,
																 DUAL_PHASEI);
		ILLprice_compute_primal_inf (lp, pinf, &lindex, 1, DUAL_PHASEI);
		ILLfct_update_counts (lp, CNT_YRAVG, lp->yjz.nzcnt, zeroLpNum);
	}
	else
		ILLprice_update_mpartial_price (lp, pinf, DUAL_PHASEI, ROW_PRICING);
	EGlpNumClearVar (ntmp);
}

void ILLfct_update_dIIfeas (
	lpinfo * lp,
	int eindex,
	svector * srhs)
{
	int j, k;
	int col, indx, vs;
	int *perm = lp->upd.perm;
	int *ix = lp->upd.ix;
	int tctr = lp->upd.tctr;
	EGlpNum_t *zAj, *l, *u;
	EGlpNum_t *dty = &(lp->upd.dty);
	EGlpNum_t *t_max = &(lp->upd.tz);
	EGlpNum_t *t = lp->upd.t;
	EGlpNum_t delta;
	svector a;

	EGlpNumInitVar (delta);
	EGlpNumZero (delta);
	EGlpNumZero (*dty);

	srhs->nzcnt = 0;
	for (j = 0; j < tctr && EGlpNumIsLeq (t[perm[j]], *t_max); j++)
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
				EGlpNumCopyDiff (delta, *l, *u);
			else
				EGlpNumCopyDiff (delta, *u, *l);
			EGlpNumAddInnProdTo (*dty, delta, *zAj);
			lp->vstat[col] = (vs == STAT_UPPER) ? STAT_LOWER : STAT_UPPER;

			a.nzcnt = lp->matcnt[col];
			a.indx = &(lp->matind[lp->matbeg[col]]);
			a.coef = &(lp->matval[lp->matbeg[col]]);
			add_vectors (lp, srhs, &a, srhs, delta);
		}
	}
	EGlpNumClearVar (delta);
}

void ILLfct_compute_dpIIy (
	lpinfo * lp,
	svector * srhs,
	svector * ssoln)
{
	if (srhs->nzcnt != 0)
	{
		ILLbasis_column_solve (lp, srhs, ssoln);
	}
}

void ILLfct_update_dpII_prices (
	lpinfo * lp,
	price_info * pinf,
	svector * srhs,
	svector * ssoln,
	/*int eindex,*/
	int lindex,
	EGlpNum_t eval,
	EGlpNum_t alpha)
{
	int i;
	svector *u;

	if (srhs->nzcnt == 0)
	{
		ILLfct_update_xz (lp, alpha, -1, -1);
		u = &(lp->yjz);
	}
	else
	{
		if (ssoln->nzcnt != 0)
			for (i = 0; i < ssoln->nzcnt; i++)
				EGlpNumSubTo (lp->xbz[ssoln->indx[i]], ssoln->coef[i]);
		ILLfct_update_xz (lp, alpha, -1, -1);
		add_vectors (lp, ssoln, &(lp->yjz), ssoln, oneLpNum);
		u = ssoln;
	}
	EGlpNumCopySum (lp->xbz[lindex], eval, alpha);

	if (pinf->d_strategy == COMPLETE_PRICING)
	{
		ILLprice_compute_primal_inf (lp, pinf, u->indx, u->nzcnt, DUAL_PHASEII);
		ILLprice_compute_primal_inf (lp, pinf, &lindex, 1, DUAL_PHASEII);
		ILLfct_update_counts (lp, CNT_YRAVG, u->nzcnt, zeroLpNum);
	}
	else
		ILLprice_update_mpartial_price (lp, pinf, DUAL_PHASEII, ROW_PRICING);
}

int ILLfct_test_pivot (
	lpinfo * lp,
	int indx,
	int indxtype,
	EGlpNum_t piv_val)
{
	int i;
	EGlpNum_t pval, ntmp;

	EGlpNumInitVar (pval);
	EGlpNumInitVar (ntmp);
	EGlpNumZero (pval);

	if (indxtype == ROW_PIVOT)
	{
		for (i = 0; i < lp->yjz.nzcnt; i++)
			if (lp->yjz.indx[i] == indx)
			{
				EGlpNumCopy (pval, lp->yjz.coef[i]);
				break;
			}
	}
	else
	{
		for (i = 0; i < lp->zA.nzcnt; i++)
			if (lp->zA.indx[i] == indx)
			{
				EGlpNumCopy (pval, lp->zA.coef[i]);
				break;
			}
	}
	EGlpNumCopyDiff (ntmp, pval, piv_val);
	EGlpNumDivTo (ntmp, piv_val);
	if (EGlpNumIsLessZero (ntmp))
		EGlpNumSign (ntmp);
	if (EGlpNumIsLess (ALTPIV_TOLER, ntmp))
	{
#if FCT_DEBUG > 1
		if (indxtype == ROW_PIVOT)
			printf ("y_i = %.8f, z_j = %.8f %la %la\n", EGlpNumToLf (pval),
							EGlpNumToLf (piv_val), EGlpNumToLf (ALTPIV_TOLER),
							EGlpNumToLf (ntmp));
		else
			printf ("z_j = %.8f, y_i = %.8f\n", EGlpNumToLf (pval),
							EGlpNumToLf (piv_val));
#endif
		EGlpNumClearVar (ntmp);
		EGlpNumClearVar (pval);
		return 1;
	}
	EGlpNumClearVar (pval);
	EGlpNumClearVar (ntmp);
	return 0;
}

#if FCT_DEBUG > 0

void fct_test_workvector (
	lpinfo * lp)
{
	int i, err = 0;

	for (i = 0; i < lp->ncols; i++)
	{
		if (EGlpNumIsNeqqZero (lp->work.coef[i]))
		{
			err++;
			EGlpNumZero (lp->work.coef[i]);
		}
		if (lp->iwork[i] != 0)
		{
			err++;
			lp->iwork[i] = 0;
		}
	}
	if (err)
		printf ("bad work vector, err=%d\n", err);
}

void fct_test_pfeasible (
	lpinfo * lp)
{
	int i, col;
	int err = 0;
	EGlpNum_t *ftol = &(lp->tol->pfeas_tol);

	for (i = 0; i < lp->nrows; i++)
	{
		col = lp->baz[i];

		if (EGlpNumIsNeqq (lp->uz[col], INFTY)
				&& EGlpNumIsSumLess (*ftol, lp->uz[col], lp->xbz[i]))
		{
			if (lp->bfeas[i] != 1)
			{
				err++;
				lp->bfeas[i] = 1;
			}
		}
		else if (EGlpNumIsNeqq (lp->lz[col], NINFTY)
						 && EGlpNumIsSumLess (lp->xbz[i], *ftol, lp->lz[col]))
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
		printf ("test_pfeas err =%d\n", err);
}

void fct_test_dfeasible (
	lpinfo * lp)
{
	int j, col;
	int err = 0;
	EGlpNum_t *ftol = &(lp->tol->dfeas_tol);
	EGlpNum_t mftol[1];

	EGlpNumInitVar (mftol[0]);
	EGlpNumCopyNeg (mftol[0], *ftol);

	for (j = 0; j < lp->nnbasic; j++)
	{
		col = lp->nbaz[j];

		if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED)
			continue;
		if (EGlpNumIsLess (lp->dz[j], mftol[0]) &&
				(lp->vstat[col] == STAT_LOWER || lp->vstat[col] == STAT_ZERO))
		{
			if (lp->dfeas[j] != -1)
			{
				err++;
				lp->dfeas[j] = -1;
			}
		}
		if (EGlpNumIsLess (*ftol, lp->dz[j]) &&
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
		printf ("test_dfeas err =%d\n", err);
}

void fct_test_pI_x (
	lpinfo * lp,
	price_info * p)
{
	int i;
	int ern = 0;
	EGlpNum_t *x;
	EGlpNum_t err, diff;

	EGlpNumInitVar (err);
	EGlpNumInitVar (diff);
	EGlpNumZero (err);
	x = EGlpNumAllocArray (lp->nrows);

	for (i = 0; i < lp->nrows; i++)
		EGlpNumCopy (x[i], lp->xbz[i]);
	ILLfct_compute_phaseI_xbz (lp);
	for (i = 0; i < lp->nrows; i++)
	{
		EGlpNumCopyDiff (diff, x[i], lp->xbz[i]);
		if (EGlpNumIsLessZero (diff))
			EGlpNumSign (diff);
		if (EGlpNumIsLess (PFEAS_TOLER, diff))
		{
			EGlpNumAddTo (err, diff);
			ern++;
			printf ("bad i = %d\n", i);
		}
	}
	if (EGlpNumIsNeqqZero (err))
		printf ("dI x err = %.7f, ern = %d\n", EGlpNumToLf (err), ern);
	ILLprice_compute_primal_inf (lp, p, NULL, 0, DUAL_PHASEI);
	EGlpNumFreeArray (x);
	EGlpNumClearVar (diff);
	EGlpNumClearVar (err);
}

void fct_test_pII_x (
	lpinfo * lp,
	price_info * p)
{
	int i;
	int ern = 0;
	EGlpNum_t *x;
	EGlpNum_t err, diff;

	EGlpNumInitVar (err);
	EGlpNumInitVar (diff);
	EGlpNumZero (err);
	x = EGlpNumAllocArray (lp->nrows);

	for (i = 0; i < lp->nrows; i++)
		EGlpNumCopy (x[i], lp->xbz[i]);
	ILLfct_compute_xbz (lp);
	for (i = 0; i < lp->nrows; i++)
	{
		EGlpNumCopyDiff (diff, x[i], lp->xbz[i]);
		if (EGlpNumIsLessZero (diff ))
			EGlpNumSign (diff);
		if (EGlpNumIsLess (PFEAS_TOLER, diff))
		{
			EGlpNumAddTo (err, diff);
			ern++;
			printf ("bad i = %d\n", i);
		}
	}
	if (EGlpNumIsNeqqZero (err))
		printf ("dII x err = %.7f, ern = %d\n", EGlpNumToLf (err), ern);
	ILLprice_compute_primal_inf (lp, p, NULL, 0, DUAL_PHASEII);
	EGlpNumFreeArray (x);
	EGlpNumClearVar (diff);
	EGlpNumClearVar (err);
}

void fct_test_pI_pi_dz (
	lpinfo * lp,
	price_info * p)
{
	int i;
	int ern = 0;
	EGlpNum_t *pidz;
	EGlpNum_t err, diff;

	EGlpNumInitVar (err);
	EGlpNumInitVar (diff);
	pidz = EGlpNumAllocArray (lp->ncols);
	EGlpNumZero (err);

	for (i = 0; i < lp->nrows; i++)
		EGlpNumCopy (pidz[i], lp->pIpiz[i]);
	ILLfct_compute_phaseI_piz (lp);
	for (i = 0; i < lp->nrows; i++)
	{
		EGlpNumCopyDiff (diff, pidz[i], lp->pIpiz[i]);
		if (EGlpNumIsLessZero (diff))
			EGlpNumSign (diff);
		if (EGlpNumIsLess (DFEAS_TOLER, diff))
		{
			EGlpNumAddTo (err, diff);
			ern++;
		}
	}
	if (EGlpNumIsNeqqZero (err))
		printf ("pI pi err = %.7f, ern = %d\n", EGlpNumToLf (err), ern);

	EGlpNumZero (err);
	ern = 0;
	for (i = 0; i < lp->nnbasic; i++)
		EGlpNumCopy (pidz[i], lp->pIdz[i]);
	ILLfct_compute_phaseI_dz (lp);
	for (i = 0; i < lp->nnbasic; i++)
	{
		EGlpNumCopyDiff (diff, pidz[i], lp->pIdz[i]);
		if (EGlpNumIsLessZero (diff))
			EGlpNumSign (diff);
		if (EGlpNumIsLess (DFEAS_TOLER, diff))
		{
			EGlpNumAddTo (err, diff);
			ern++;
		}
	}
	if (EGlpNumIsNeqqZero (err))
		printf ("pI dz err = %.7f, ern = %d\n", EGlpNumToLf (err), ern);
	ILLprice_compute_dual_inf (lp, p, NULL, 0, PRIMAL_PHASEI);
	EGlpNumClearVar (err);
	EGlpNumClearVar (diff);
	EGlpNumFreeArray (pidz);
}

void fct_test_pII_pi_dz (
	lpinfo * lp,
	price_info * p)
{
	int i;
	int ern = 0;
	EGlpNum_t *pidz;
	EGlpNum_t err, diff;

	EGlpNumInitVar (err);
	EGlpNumInitVar (diff);
	EGlpNumZero (err);
	pidz = EGlpNumAllocArray (lp->ncols);

	for (i = 0; i < lp->nrows; i++)
		EGlpNumCopy (pidz[i], lp->piz[i]);
	ILLfct_compute_piz (lp);
	for (i = 0; i < lp->nrows; i++)
	{
		EGlpNumCopyDiff (diff, pidz[i], lp->piz[i]);
		if (EGlpNumIsLessZero (diff))
			EGlpNumSign (diff);
		if (EGlpNumIsLess (DFEAS_TOLER, diff))
		{
			EGlpNumAddTo (err, diff);
			ern++;
		}
	}
	if (EGlpNumIsNeqqZero (err))
		printf ("pII pi err = %.7f, ern = %d\n", EGlpNumToLf (err), ern);

	EGlpNumZero (err);
	ern = 0;
	for (i = 0; i < lp->nnbasic; i++)
		EGlpNumCopy (pidz[i], lp->dz[i]);
	ILLfct_compute_dz (lp);
	for (i = 0; i < lp->nnbasic; i++)
	{
		EGlpNumCopyDiff (diff, pidz[i], lp->dz[i]);
		if (EGlpNumIsLessZero (diff))
			EGlpNumSign (diff);
		if (EGlpNumIsLess (DFEAS_TOLER, diff))
		{
			EGlpNumAddTo (err, diff);
			ern++;
		}
	}
	if (EGlpNumIsNeqqZero (err))
		printf ("pII dz err = %.7f, ern = %d\n", EGlpNumToLf (err), ern);
	/*
	 * ILLprice_compute_dual_inf (lp, p, NULL, 0, PRIMAL_PHASEII);
	 */
	EGlpNumClearVar (err);
	EGlpNumClearVar (diff);
	EGlpNumFreeArray (pidz);
}

#endif
