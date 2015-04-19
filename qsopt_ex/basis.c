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

/* RCS_INFO = "$RCSfile: basis.c,v $ $Revision: 1.2 $ $Date: 2003/11/05 16:49:52 $"; */
static int TRACE = 0;

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "qs_config.h"
#include "logging-private.h"

#include "eg_lpnum.h"
#include "eg_io.h"
#include "except.h"
#include "util.h"

#include "sortrus_EGLPNUM_TYPENAME.h"
#include "lpdefs_EGLPNUM_TYPENAME.h"
#include "qstruct_EGLPNUM_TYPENAME.h"
#include "qsopt_EGLPNUM_TYPENAME.h"
#include "basis_EGLPNUM_TYPENAME.h"
#include "fct_EGLPNUM_TYPENAME.h"
#include "lp_EGLPNUM_TYPENAME.h"
#include "lib_EGLPNUM_TYPENAME.h"

//#define DJZERO_TOLER EGLPNUM_TYPENAME_PFEAS_TOLER
#define BASIS_STATS 0
//#define BASIS_DEBUG 10
#define BASIS_DEBUG 0

void EGLPNUM_TYPENAME_ILLbasis_init_vardata (
	EGLPNUM_TYPENAME_var_data * vd)
{
	memset (vd, 0, sizeof (EGLPNUM_TYPENAME_var_data));
	EGLPNUM_TYPENAME_EGlpNumInitVar (vd->cmax);
}

void EGLPNUM_TYPENAME_ILLbasis_clear_vardata (
	EGLPNUM_TYPENAME_var_data * vd)
{
	EGLPNUM_TYPENAME_EGlpNumClearVar (vd->cmax);
	memset (vd, 0, sizeof (EGLPNUM_TYPENAME_var_data));
}

static void get_var_info (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_var_data * v);

static int init_slack_basis (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int *vstat,
	int *irow,
	int *rrow,
	int *unitcol,
	int *icol,
	int *rcol),
  get_initial_basis1 (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int *vstat),
  get_initial_basis2 (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int *vstat),
  set_basis_indices (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int *vstat),
  choose_basis (
	int algorithm,
	EGLPNUM_TYPE pinf1,
	EGLPNUM_TYPE dinf1,
	EGLPNUM_TYPE pinf2,
	EGLPNUM_TYPE dinf2);

void EGLPNUM_TYPENAME_ILLbasis_init_basisinfo (
	EGLPNUM_TYPENAME_lpinfo * lp)
{
	lp->baz = 0;
	lp->nbaz = 0;
	lp->vstat = 0;
	lp->vindex = 0;
	lp->f = 0;
}

void EGLPNUM_TYPENAME_ILLbasis_free_basisinfo (
	EGLPNUM_TYPENAME_lpinfo * lp)
{
	ILL_IFFREE (lp->baz, int);
	ILL_IFFREE (lp->nbaz, int);
	ILL_IFFREE (lp->vstat, int);
	ILL_IFFREE (lp->vindex, int);

	if (lp->f)
	{
		EGLPNUM_TYPENAME_ILLfactor_free_factor_work (lp->f);
		EGLPNUM_TYPENAME_EGlpNumClearVar (lp->f->fzero_tol);
		EGLPNUM_TYPENAME_EGlpNumClearVar (lp->f->szero_tol);
		EGLPNUM_TYPENAME_EGlpNumClearVar (lp->f->partial_tol);
		EGLPNUM_TYPENAME_EGlpNumClearVar (lp->f->maxelem_orig);
		EGLPNUM_TYPENAME_EGlpNumClearVar (lp->f->maxelem_factor);
		EGLPNUM_TYPENAME_EGlpNumClearVar (lp->f->maxelem_cur);
		EGLPNUM_TYPENAME_EGlpNumClearVar (lp->f->partial_cur);
		ILL_IFFREE (lp->f, EGLPNUM_TYPENAME_factor_work);
	}
}

int EGLPNUM_TYPENAME_ILLbasis_build_basisinfo (
	EGLPNUM_TYPENAME_lpinfo * lp)
{
	int rval = 0;

	ILL_SAFE_MALLOC (lp->baz, lp->O->nrows, int);
	ILL_SAFE_MALLOC (lp->nbaz, lp->O->ncols, int);
	ILL_SAFE_MALLOC (lp->vstat, lp->O->ncols, int);
	ILL_SAFE_MALLOC (lp->vindex, lp->O->ncols, int);

	lp->fbasisid = -1;

CLEANUP:
	if (rval)
		EGLPNUM_TYPENAME_ILLbasis_free_basisinfo (lp);
	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLbasis_load (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_ILLlp_basis * B)
{
	int rval = 0;
	char *cstat = B->cstat;
	char *rstat = B->rstat;
	int *structmap = lp->O->structmap;
	int *rowmap = lp->O->rowmap;
	char *sense = lp->O->sense;
	int i, j, ncols = lp->O->ncols, nrows = lp->O->nrows, nstruct = lp->O->nstruct;
	int basic = 0, nonbasic = 0;

	EGLPNUM_TYPENAME_ILLbasis_free_basisinfo (lp);
	EGLPNUM_TYPENAME_ILLbasis_init_basisinfo (lp);
	rval = EGLPNUM_TYPENAME_ILLbasis_build_basisinfo (lp);
	CHECKRVALG (rval, CLEANUP);

	for (i = 0; i < nstruct; i++)
	{
		j = structmap[i];
		if (cstat[i] == QS_COL_BSTAT_BASIC)
		{
			lp->vstat[j] = STAT_BASIC;
			lp->baz[basic] = j;
			lp->vindex[j] = basic;
			basic++;
		}
		else
		{
			lp->nbaz[nonbasic] = j;
			lp->vindex[j] = nonbasic;
			nonbasic++;
			switch (cstat[i])
			{
			case QS_COL_BSTAT_LOWER:
				lp->vstat[j] = STAT_LOWER;
				break;
			case QS_COL_BSTAT_UPPER:
				lp->vstat[j] = STAT_UPPER;
				break;
			case QS_COL_BSTAT_FREE:
				lp->vstat[j] = STAT_ZERO;
				break;
			default:
				QSlog("unknown col basis stat 1: %c", cstat[i]);
				rval = 1;
				goto CLEANUP;
			}
		}
	}

	for (i = 0; i < nrows; i++)
	{
		j = rowmap[i];
		if (sense[i] == 'R')
		{
			if (rstat[i] == QS_ROW_BSTAT_BASIC)
			{
				lp->vstat[j] = STAT_BASIC;
				lp->baz[basic] = j;
				lp->vindex[j] = basic;
				basic++;
			}
			else
			{
				lp->nbaz[nonbasic] = j;
				lp->vindex[j] = nonbasic;
				nonbasic++;
				switch (rstat[i])
				{
				case QS_ROW_BSTAT_LOWER:
					lp->vstat[j] = STAT_LOWER;
					break;
				case QS_ROW_BSTAT_UPPER:
					lp->vstat[j] = STAT_UPPER;
					break;
				default:
					QSlog("unknown range basis stat 2");
					rval = 1;
					goto CLEANUP;
				}
			}
		}
		else
		{
			switch (rstat[i])
			{
			case QS_ROW_BSTAT_BASIC:
				lp->vstat[j] = STAT_BASIC;
				lp->baz[basic] = j;
				lp->vindex[j] = basic;
				basic++;
				break;
			case QS_ROW_BSTAT_LOWER:
				lp->vstat[j] = STAT_LOWER;
				lp->nbaz[nonbasic] = j;
				lp->vindex[j] = nonbasic;
				nonbasic++;
				break;
			default:
				QSlog("unknown row basis stat 3");
				rval = 1;
				goto CLEANUP;
			}
		}
	}

	if (basic + nonbasic != ncols)
	{
		QSlog("error in counts in ILLopt_load_basis");
		rval = 1;
		goto CLEANUP;
	}

	if (lp->fbasisid != 0)
		lp->basisid = 0;
	else
		lp->basisid = 1;

CLEANUP:

	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLbasis_tableau_row (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int row,
	EGLPNUM_TYPE * brow,
	EGLPNUM_TYPE * trow,
	EGLPNUM_TYPE * rhs,
	int strict)
{
	int rval = 0;
	int i;
	int singular = 0;
	int indx;
	EGLPNUM_TYPE coef;
	EGLPNUM_TYPE sum;
	EGLPNUM_TYPENAME_svector z, zA;

	EGLPNUM_TYPENAME_EGlpNumInitVar (coef);
	EGLPNUM_TYPENAME_EGlpNumInitVar (sum);
	EGLPNUM_TYPENAME_EGlpNumZero (sum);

	EGLPNUM_TYPENAME_ILLsvector_init (&z);
	EGLPNUM_TYPENAME_ILLsvector_init (&zA);

	if (lp->basisid == -1)
	{
		QSlog("EGLPNUM_TYPENAME_ILLbasis_tableau_row: no basis");
		rval = E_GENERAL_ERROR;
		ILL_CLEANUP;
	}
	if (lp->fbasisid != lp->basisid)
	{															/* Needs to be changed */
		rval = EGLPNUM_TYPENAME_ILLbasis_factor (lp, &singular);
		CHECKRVALG (rval, CLEANUP);
		if (singular)
		{
			MESSAGE (__QS_SB_VERB, "Singular Basis found!");
			rval = E_BASIS_SINGULAR;
			ILL_CLEANUP;
		}
	}
	if (brow == NULL)
	{
		QSlog("No array for basis inverse row");
		rval = E_GENERAL_ERROR;
		ILL_CLEANUP;
	}

	rval = EGLPNUM_TYPENAME_ILLsvector_alloc (&z, lp->nrows);
	CHECKRVALG (rval, CLEANUP);
	EGLPNUM_TYPENAME_ILLfct_compute_zz (lp, &z, row);

	for (i = 0; i < lp->O->nrows; i++)
		EGLPNUM_TYPENAME_EGlpNumZero (brow[i]);
	for (i = 0; i < z.nzcnt; i++)
	{
		indx = z.indx[i];
		EGLPNUM_TYPENAME_EGlpNumCopy (coef, z.coef[i]);
		EGLPNUM_TYPENAME_EGlpNumCopy (brow[indx], coef);
		EGLPNUM_TYPENAME_EGlpNumAddInnProdTo (sum, coef, lp->bz[indx]);
	}

	if (rhs != NULL)
		EGLPNUM_TYPENAME_EGlpNumCopy (*rhs, sum);
	if (trow != NULL)
	{
		if (!strict)
		{
			rval = EGLPNUM_TYPENAME_ILLsvector_alloc (&zA, lp->ncols);
			if (rval)
				ILL_CLEANUP;
			ILL_IFTRACE ("%s:\n", __func__);
			rval = EGLPNUM_TYPENAME_ILLfct_compute_zA (lp, &z, &zA);
			CHECKRVALG (rval, CLEANUP);

			for (i = 0; i < lp->ncols; i++)
				EGLPNUM_TYPENAME_EGlpNumZero (trow[i]);
			for (i = 0; i < zA.nzcnt; i++)
				EGLPNUM_TYPENAME_EGlpNumCopy (trow[lp->nbaz[zA.indx[i]]], zA.coef[i]);
			EGLPNUM_TYPENAME_EGlpNumOne (trow[lp->baz[row]]);
		}
		else
		{
			EGLPNUM_TYPENAME_ILLfct_compute_vA (lp, &z, trow);
		}
	}

#if BASIS_DEBUG > 0
	if (rhs != NULL && trow != NULL)
	{
		EGLPNUM_TYPE *tr = NULL;

		EGLPNUM_TYPENAME_EGlpNumZero (sum);
		if (strict)
			tr = trow;
		else
		{
			tr = EGLPNUM_TYPENAME_EGlpNumAllocArray (lp->ncols);
			EGLPNUM_TYPENAME_ILLfct_compute_vA (lp, &z, tr);
		}
		for (i = 0; i < lp->nrows; i++)
			if (EGLPNUM_TYPENAME_EGlpNumIsGreatZero (tr[lp->baz[i]]))
				EGLPNUM_TYPENAME_EGlpNumAddTo (sum, tr[lp->baz[i]]);
			else
				EGLPNUM_TYPENAME_EGlpNumSubTo (sum, tr[lp->baz[i]]);
		EGLPNUM_TYPENAME_EGlpNumCopy (coef, EGLPNUM_TYPENAME_oneLpNum);
		EGLPNUM_TYPENAME_EGlpNumSubTo (coef, sum);
		if (EGLPNUM_TYPENAME_EGlpNumIsLessZero (coef))
			EGLPNUM_TYPENAME_EGlpNumSign (coef);
		if (EGLPNUM_TYPENAME_EGlpNumIsLess (EGLPNUM_TYPENAME_PIVZ_TOLER, coef))
			QSlog("tableau: bas computed = %.12f", EGLPNUM_TYPENAME_EGlpNumToLf (sum));
		if (!strict)
			EGLPNUM_TYPENAME_EGlpNumFreeArray (tr);
#if BASIS_DEBUG > 1
		EGLPNUM_TYPENAME_EGlpNumZero (sum);
		for (i = 0; i < lp->ncols; i++)
		{
			if (lp->vstat[i] == STAT_BASIC)
				EGLPNUM_TYPENAME_EGlpNumAddInnProdTo (sum, lp->xbz[lp->vindex[i]], trow[i]);
			else if (lp->vstat[i] == STAT_UPPER)
				EGLPNUM_TYPENAME_EGlpNumAddInnProdTo (sum, lp->uz[i], trow[i]);
			else if (lp->vstat[i] == STAT_LOWER)
				EGLPNUM_TYPENAME_EGlpNumAddInnProdTo (sum, lp->lz[i], trow[i]);
		}
		EGLPNUM_TYPENAME_EGlpNumSet (coef, 1e-10);
		if (EGLPNUM_TYPENAME_EGlpNumIsNeq (sum, *rhs, coef))
			QSlog("tableau rhs = %.9f, computed = %.9f",
									EGLPNUM_TYPENAME_EGlpNumToLf (*rhs), EGLPNUM_TYPENAME_EGlpNumToLf (sum));
#endif
	}
#endif

CLEANUP:
	EGLPNUM_TYPENAME_ILLsvector_free (&z);
	EGLPNUM_TYPENAME_ILLsvector_free (&zA);
	EGLPNUM_TYPENAME_EGlpNumClearVar (coef);
	EGLPNUM_TYPENAME_EGlpNumClearVar (sum);
	return rval;
}

static void get_var_info (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_var_data * v)
{
	int i = 0;

	v->nartif = 0;
	v->nslacks = 0;
	v->nfree = 0;
	v->nbndone = 0;
	v->nbounded = 0;
	v->nfixed = 0;
	EGLPNUM_TYPENAME_EGlpNumCopy (v->cmax, EGLPNUM_TYPENAME_NINFTY);

	for (i = 0; i < lp->ncols; i++)
	{
		switch (lp->vtype[i])
		{
		case VARTIFICIAL:
			v->nartif++;
			break;
		case VFREE:
			v->nfree++;
			break;
		case VLOWER:
		case VUPPER:
			if (lp->vclass[i] == CLASS_LOGICAL)
				v->nslacks++;
			else
				v->nbndone++;
			break;

		case VFIXED:
			v->nfixed++;
		case VBOUNDED:
			if (lp->vclass[i] == CLASS_LOGICAL)
				v->nslacks++;
			else
				v->nbounded++;
			break;
		}
		EGLPNUM_TYPENAME_EGlpNumSetToMaxAbs (v->cmax, lp->cz[i]);
	}

#if BASIS_STATS > 0
	QSlog("cols = %d, acols = %d, total  = %d, nrows = %d, nlog = %d",
							lp->ncols, lp->ncols - lp->nrows,
							v->nartif + v->nfree + v->nslacks + v->nbndone + v->nbounded,
							lp->nrows, v->nartif + v->nslacks);
#endif
}

static int init_slack_basis (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int *vstat,
	int *irow,
	int *rrow,
	int *unitcol,
	int *icol,
	int *rcol)
{
	int j, r, vt;
	int nslacks = 0;

	for (j = 0; j < lp->ncols; j++)
	{
		r = lp->matind[lp->matbeg[j]];
		vt = lp->vtype[j];

		if ((vt == VUPPER || vt == VLOWER || vt == VBOUNDED || vt == VFIXED) &&
				lp->vclass[j] == CLASS_LOGICAL)
		{

			vstat[j] = STAT_BASIC;
			irow[r] = 1;
			rrow[r] = 1;
			unitcol[r] = j;
			if (icol != NULL)
			{
				icol[j] = 1;
				rcol[j] = 1;
			}
			nslacks++;
		}
		else if (vt == VARTIFICIAL)
		{
			unitcol[r] = j;
			vstat[j] = STAT_UPPER;
		}
		else if (vt == VFREE)
			vstat[j] = STAT_ZERO;
		else if (vt == VFIXED || vt == VUPPER)
			vstat[j] = STAT_UPPER;
		else if (vt == VLOWER)
			vstat[j] = STAT_LOWER;
		else if (vt == VBOUNDED)
		{
			if (fabs (EGLPNUM_TYPENAME_EGlpNumToLf (lp->lz[j])) < fabs (EGLPNUM_TYPENAME_EGlpNumToLf (lp->uz[j])))
				vstat[j] = STAT_LOWER;
			else
				vstat[j] = STAT_UPPER;
		}
	}
	return nslacks;
}

static int primal_col_select (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int *vstat,
	int *irow,
	int *rrow,
	int *unitcol,
	EGLPNUM_TYPE * v,
	int *perm,
	int *porder,
	int nbelem,
	int pcols)
{
	int i, j, k, tr, r = 0;
	int mcnt, mbeg;
	int *matbeg = lp->matbeg;
	int *matcnt = lp->matcnt;
	int *matind = lp->matind;
	EGLPNUM_TYPE *matval = lp->matval;
	EGLPNUM_TYPE alpha, val, maxelem;

	EGLPNUM_TYPENAME_EGlpNumInitVar (alpha);
	EGLPNUM_TYPENAME_EGlpNumInitVar (val);
	EGLPNUM_TYPENAME_EGlpNumInitVar (maxelem);

	for (k = 0; k < pcols; k++)
	{
		j = porder[perm[k]];
		mcnt = matcnt[j];
		mbeg = matbeg[j];

		EGLPNUM_TYPENAME_EGlpNumCopy (alpha, EGLPNUM_TYPENAME_NINFTY);
		EGLPNUM_TYPENAME_EGlpNumCopy (maxelem, EGLPNUM_TYPENAME_NINFTY);

		for (i = 0; i < mcnt; i++)
		{
			EGLPNUM_TYPENAME_EGlpNumCopyAbs (val, matval[mbeg + i]);
			if (EGLPNUM_TYPENAME_EGlpNumIsLess (maxelem, val))
				EGLPNUM_TYPENAME_EGlpNumCopy (maxelem, val);
			if (rrow[matind[mbeg + i]] == 0 && EGLPNUM_TYPENAME_EGlpNumIsLess (alpha, val))
			{
				EGLPNUM_TYPENAME_EGlpNumCopy (alpha, val);
				r = matind[mbeg + i];
			}
		}
		EGLPNUM_TYPENAME_EGlpNumCopy (val, maxelem);
		EGLPNUM_TYPENAME_EGlpNumMultTo (val, EGLPNUM_TYPENAME_PARAM_IBASIS_RPIVOT);
		if (EGLPNUM_TYPENAME_EGlpNumIsLess (val, alpha))
		{
			vstat[j] = STAT_BASIC;
			nbelem++;
			irow[r] = 1;
			EGLPNUM_TYPENAME_EGlpNumCopy (v[r], alpha);
			for (i = 0; i < mcnt; i++)
				if (EGLPNUM_TYPENAME_EGlpNumIsNeqqZero (matval[mbeg + i]))
					rrow[matind[mbeg + i]]++;
		}
		else
		{
			EGLPNUM_TYPENAME_EGlpNumCopy (alpha, EGLPNUM_TYPENAME_NINFTY);
			for (i = 0; i < mcnt; i++)
			{
				tr = matind[mbeg + i];
				EGLPNUM_TYPENAME_EGlpNumCopyAbs (val, matval[mbeg + i]);
				EGLPNUM_TYPENAME_EGlpNumDivTo (val, EGLPNUM_TYPENAME_PARAM_IBASIS_RTRIANG);
				if (EGLPNUM_TYPENAME_EGlpNumIsNeqq (v[tr], EGLPNUM_TYPENAME_INFTY) && EGLPNUM_TYPENAME_EGlpNumIsLess (v[tr], val))
				{
					EGLPNUM_TYPENAME_EGlpNumZero (alpha);
					break;
				}
				EGLPNUM_TYPENAME_EGlpNumCopyAbs (val, matval[mbeg + i]);
				if (irow[tr] == 0 && EGLPNUM_TYPENAME_EGlpNumIsLess (alpha, val))
				{
					EGLPNUM_TYPENAME_EGlpNumCopy (alpha, val);
					r = tr;
				}
			}
			if (EGLPNUM_TYPENAME_EGlpNumIsNeqqZero (alpha) && EGLPNUM_TYPENAME_EGlpNumIsNeqq (alpha, EGLPNUM_TYPENAME_NINFTY))
			{
				vstat[j] = STAT_BASIC;
				nbelem++;
				irow[r] = 1;
				EGLPNUM_TYPENAME_EGlpNumCopy (v[r], alpha);
				for (i = 0; i < mcnt; i++)
					if (EGLPNUM_TYPENAME_EGlpNumIsNeqqZero (matval[mbeg + i]))
						rrow[matind[mbeg + i]]++;
			}
		}
	}
#if BASIS_STATS > 0
	QSlog("nartifs = %d", lp->nrows - nbelem);
#endif

	if (nbelem < lp->nrows)
	{
		for (i = 0; i < lp->nrows; i++)
		{
			if (irow[i] == 0)
			{
				if (unitcol[i] != -1)
				{
					vstat[unitcol[i]] = STAT_BASIC;
					nbelem++;
				}
				else
				{
					QSlog("Error: Not enough artificials");
					return -1;
				}
			}
		}
	}
	EGLPNUM_TYPENAME_EGlpNumClearVar (alpha);
	EGLPNUM_TYPENAME_EGlpNumClearVar (val);
	EGLPNUM_TYPENAME_EGlpNumClearVar (maxelem);
	return nbelem;
}

/* This is an implementation of the initial basis procedure
   in: "Implementing the simplex method: the initial basis", by
   Bob Bixby.
   Goals: choose initial variables to go into basis which satisfy:
   1) vars are slacks, 2) vars have freedom to move
   3) initial submatrix is nonsingular, 4) low objective function
   contribution.
*/
static int get_initial_basis1 (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int *vstat)
{
	int rval = 0;
	int i, j, tot1 = 0, tot2 = 0;
	int nbelem = 0, nslacks = 0;
	int tfree = 0, tbndone = 0;
	int tbounded = 0;
	int *irow = NULL, *rrow = NULL;
	int *perm = NULL, *porder = NULL;
	int *unitcol = NULL;
	EGLPNUM_TYPE cmax;
	EGLPNUM_TYPE *v = NULL;
	EGLPNUM_TYPE *qpenalty = NULL;
	EGLPNUM_TYPENAME_var_data vd;

	EGLPNUM_TYPENAME_ILLbasis_init_vardata (&vd);
	EGLPNUM_TYPENAME_EGlpNumInitVar (cmax);

	get_var_info (lp, &vd);
	if (!EGLPNUM_TYPENAME_EGlpNumIsNeqqZero (vd.cmax))
		EGLPNUM_TYPENAME_EGlpNumOne (cmax);
	else
	{
		EGLPNUM_TYPENAME_EGlpNumCopy (cmax, vd.cmax);
		EGLPNUM_TYPENAME_EGlpNumMultUiTo (cmax, 1000);
	}

	ILL_SAFE_MALLOC (irow, lp->nrows, int);
	ILL_SAFE_MALLOC (rrow, lp->nrows, int);

	v = EGLPNUM_TYPENAME_EGlpNumAllocArray (lp->nrows);
	ILL_SAFE_MALLOC (unitcol, lp->nrows, int);

	for (i = 0; i < lp->nrows; i++)
	{
		unitcol[i] = -1;
		EGLPNUM_TYPENAME_EGlpNumCopy (v[i], EGLPNUM_TYPENAME_INFTY);
		irow[i] = 0;
		rrow[i] = 0;
	}

	nslacks = init_slack_basis (lp, vstat, irow, rrow, unitcol, NULL, NULL);
	if (nslacks != vd.nslacks)
	{
		QSlog("complain: incorrect basis info(slacks)");
		rval = E_SIMPLEX_ERROR;
		ILL_CLEANUP;
	}
	if (nslacks == lp->nrows)
		ILL_CLEANUP;
	nbelem = nslacks;
	if (nbelem < lp->nrows)
	{
		for (i = 0; i < lp->nrows; i++)
		{
			if (irow[i] == 0)
			{
				if (unitcol[i] != -1)
				{
					vstat[unitcol[i]] = STAT_BASIC;
					nbelem++;
				}
				else
				{
					QSlog("Error: Not enough artificials");
					return -1;
				}
			}
		}
	}
	ILL_CLEANUP;

	tot1 = vd.nfree + vd.nbndone;
	tot2 = vd.nfree + vd.nbndone + vd.nbounded;
	ILL_SAFE_MALLOC (perm, tot2, int);
	ILL_SAFE_MALLOC (porder, tot2, int);

	qpenalty = EGLPNUM_TYPENAME_EGlpNumAllocArray (tot2);

	for (j = 0; j < lp->ncols; j++)
	{
		if (vstat[j] == STAT_BASIC)
			continue;

		switch (lp->vtype[j])
		{
		case VFREE:
			porder[tfree] = j;
			perm[tfree] = tfree;
			EGLPNUM_TYPENAME_EGlpNumCopyFrac (qpenalty[tfree], lp->cz[j], cmax);
			tfree++;
			break;

		case VLOWER:
		case VUPPER:
			porder[vd.nfree + tbndone] = j;
			perm[vd.nfree + tbndone] = tbndone;
			EGLPNUM_TYPENAME_EGlpNumCopyFrac (qpenalty[vd.nfree + tbndone], lp->cz[j], cmax);
			if (lp->vtype[j] == VLOWER)
				EGLPNUM_TYPENAME_EGlpNumAddTo (qpenalty[vd.nfree + tbndone], lp->lz[j]);
			else
				EGLPNUM_TYPENAME_EGlpNumSubTo (qpenalty[vd.nfree + tbndone], lp->uz[j]);
			tbndone++;
			break;

		case VFIXED:
		case VBOUNDED:
			porder[tot1 + tbounded] = j;
			perm[tot1 + tbounded] = tbounded;
			EGLPNUM_TYPENAME_EGlpNumCopyFrac (qpenalty[tot1 + tbndone], lp->cz[j], cmax);
			EGLPNUM_TYPENAME_EGlpNumAddTo (qpenalty[tot1 + tbndone], lp->lz[j]);
			EGLPNUM_TYPENAME_EGlpNumSubTo (qpenalty[tot1 + tbndone], lp->uz[j]);
			tbounded++;
			break;
		}
	}
	if (tfree != vd.nfree || tbndone != vd.nbndone || tbounded != vd.nbounded)
	{
		QSlog("complain: incorrect basis info");
		rval = E_SIMPLEX_ERROR;
		ILL_CLEANUP;
	}

	EGLPNUM_TYPENAME_ILLutil_EGlpNum_perm_quicksort (perm, qpenalty, vd.nfree);
	EGLPNUM_TYPENAME_ILLutil_EGlpNum_perm_quicksort (perm + vd.nfree, qpenalty + vd.nfree,
																	vd.nbndone);
	EGLPNUM_TYPENAME_ILLutil_EGlpNum_perm_quicksort (perm + tot1, qpenalty + tot1, vd.nbounded);

	for (i = 0; i < vd.nbndone; i++)
		perm[vd.nfree + i] += vd.nfree;
	for (i = 0; i < vd.nbounded; i++)
		perm[tot1 + i] += tot1;

	nbelem =
		primal_col_select (lp, vstat, irow, rrow, unitcol, v, perm, porder, nbelem,
											 tot2);
	if (nbelem != lp->nrows)
	{
		QSlog("complain: incorrect final basis size");
		rval = E_SIMPLEX_ERROR;
		ILL_CLEANUP;
	}

CLEANUP:
	EGLPNUM_TYPENAME_EGlpNumClearVar (cmax);
	if (rval)
		EGLPNUM_TYPENAME_ILLbasis_free_basisinfo (lp);
	ILL_IFFREE (irow, int);
	ILL_IFFREE (rrow, int);

	EGLPNUM_TYPENAME_EGlpNumFreeArray (v);
	ILL_IFFREE (perm, int);
	ILL_IFFREE (porder, int);
	ILL_IFFREE (unitcol, int);

	EGLPNUM_TYPENAME_EGlpNumFreeArray (qpenalty);
	EGLPNUM_TYPENAME_ILLbasis_clear_vardata (&vd);
	EG_RETURN (rval);
}

static int get_initial_basis2 (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int *vstat)
{
	int rval = 0;
	int i, j, k, tot1, tot2;
	int rbeg, rcnt, mcnt;
	int nbelem = 0, nslacks = 0;
	int tfree = 0, tbndone = 0;
	int tbounded = 0;
	int *irow = NULL, *rrow = NULL;
	int *perm = NULL, *porder = NULL;
	int *unitcol = NULL;
	EGLPNUM_TYPE *v = NULL;
	EGLPNUM_TYPE *qpenalty = NULL;
	int col = 0, s_i = 0, selc = 0;
	int *icol = NULL, *rcol = NULL;
	int *plen = NULL;
	EGLPNUM_TYPE *dj = NULL;
	EGLPNUM_TYPENAME_var_data vd;
	EGLPNUM_TYPE seldj;
	EGLPNUM_TYPE selv;
	EGLPNUM_TYPE c_dj;
	EGLPNUM_TYPE cmax;

	EGLPNUM_TYPENAME_EGlpNumInitVar (seldj);
	EGLPNUM_TYPENAME_EGlpNumInitVar (selv);
	EGLPNUM_TYPENAME_EGlpNumInitVar (c_dj);
	EGLPNUM_TYPENAME_EGlpNumInitVar (cmax);
	EGLPNUM_TYPENAME_EGlpNumZero (c_dj);
	EGLPNUM_TYPENAME_EGlpNumZero (selv);
	EGLPNUM_TYPENAME_EGlpNumZero (seldj);
	EGLPNUM_TYPENAME_ILLbasis_init_vardata (&vd);

	get_var_info (lp, &vd);

	ILL_SAFE_MALLOC (irow, lp->nrows, int);
	ILL_SAFE_MALLOC (rrow, lp->nrows, int);

	v = EGLPNUM_TYPENAME_EGlpNumAllocArray (lp->nrows);
	ILL_SAFE_MALLOC (unitcol, lp->nrows, int);
	ILL_SAFE_MALLOC (icol, lp->ncols, int);
	ILL_SAFE_MALLOC (rcol, lp->ncols, int);

	dj = EGLPNUM_TYPENAME_EGlpNumAllocArray (lp->ncols);

	for (i = 0; i < lp->nrows; i++)
	{
		unitcol[i] = -1;
		EGLPNUM_TYPENAME_EGlpNumCopy (v[i], EGLPNUM_TYPENAME_INFTY);
		irow[i] = 0;
		rrow[i] = 0;
	}
	/* assign all d_j */
	for (i = 0; i < lp->ncols; i++)
	{
		icol[i] = 0;
		rcol[i] = 0;
		EGLPNUM_TYPENAME_EGlpNumCopy (dj[i], lp->cz[i]);
	}

	nslacks = init_slack_basis (lp, vstat, irow, rrow, unitcol, icol, rcol);
	if (nslacks != vd.nslacks)
	{
		QSlog("complain: incorrect basis info");
		rval = E_SIMPLEX_ERROR;
		ILL_CLEANUP;
	}
	if (nslacks == lp->nrows)
		ILL_CLEANUP;
	nbelem = nslacks;

	/* allocate maximum required space for perm etc. */
	ILL_SAFE_MALLOC (perm, lp->ncols, int);
	ILL_SAFE_MALLOC (porder, lp->ncols, int);
	ILL_SAFE_MALLOC (plen, lp->nrows, int);

	qpenalty = EGLPNUM_TYPENAME_EGlpNumAllocArray (lp->ncols);

	/* find all unit rows and record lengths */
	for (i = 0; i < lp->nrows; i++)
	{
		if (irow[i] != 1)
		{
			rbeg = lp->rowbeg[i];
			rcnt = lp->rowcnt[i];
			for (j = 0; j < rcnt; j++)
			{
				EGLPNUM_TYPENAME_EGlpNumCopyAbs (cmax, lp->rowval[rbeg + j]);
				if (EGLPNUM_TYPENAME_EGlpNumIsNeqq (cmax, EGLPNUM_TYPENAME_oneLpNum))
					break;
			}
			if (j == rcnt)
			{
				perm[s_i] = s_i;
				porder[s_i] = i;
				plen[s_i] = rcnt;
				s_i++;
			}
		}
	}

	/*sort all unit rows */
	ILLutil_int_perm_quicksort (perm, plen, s_i);

	/* now go through the unit rows */
	for (k = 0; k < s_i; k++)
	{
		i = porder[perm[k]];
		rbeg = lp->rowbeg[i];
		rcnt = lp->rowcnt[i];
		selc = -1;
		EGLPNUM_TYPENAME_EGlpNumCopy (seldj, EGLPNUM_TYPENAME_INFTY);
		EGLPNUM_TYPENAME_EGlpNumZero (selv);

		/* for every row s_i, compute min {d_j : d_j <0 , j is u or l or fr} */
		for (j = 0; j < rcnt; j++)
		{
			col = lp->rowind[rbeg + j];
			if (rcol[col] == 1)
				break;
			if (EGLPNUM_TYPENAME_EGlpNumIsLessZero (dj[col]))
			{
				if (EGLPNUM_TYPENAME_EGlpNumIsLess (dj[col], seldj))
				{
					selc = col;
					EGLPNUM_TYPENAME_EGlpNumCopy (seldj, dj[col]);
					EGLPNUM_TYPENAME_EGlpNumCopy (selv, lp->rowval[rbeg + j]);
				}
			}
		}
		/* select pivot element and update all d_j's */
		if (selc != -1)
		{
			nbelem++;
			irow[i] = 1;
			rrow[i] = 1;
			icol[selc] = 1;
			EGLPNUM_TYPENAME_EGlpNumCopyFrac (c_dj, dj[selc], selv);
			vstat[selc] = STAT_BASIC;
			for (j = 0; j < rcnt; j++)
			{
				col = lp->rowind[rbeg + j];
				EGLPNUM_TYPENAME_EGlpNumSubInnProdTo (dj[col], lp->rowval[rbeg + j], c_dj);
				rcol[col] = 1;
			}
		}
	}
#if BASIS_STATS > 0
	QSlog("unit rows = %d", s_i);
	QSlog("nslacks %d, unit rows selected = %d", nslacks, nbelem - nslacks);
#endif
	/* now go through remaining cols with dj = 0 */
	tot1 = vd.nfree + vd.nbndone;

	if (!EGLPNUM_TYPENAME_EGlpNumIsNeqqZero (vd.cmax))
		EGLPNUM_TYPENAME_EGlpNumOne (cmax);
	else
	{
		EGLPNUM_TYPENAME_EGlpNumCopy (cmax, vd.cmax);
		EGLPNUM_TYPENAME_EGlpNumMultUiTo (cmax, 1000);
	}
	for (j = 0; j < lp->ncols; j++)
	{
		if (vstat[j] == STAT_BASIC)
			continue;
		if (icol[j] == 1 || EGLPNUM_TYPENAME_EGlpNumIsNeqZero (dj[j], EGLPNUM_TYPENAME_BD_TOLER))
			continue;
		mcnt = lp->matcnt[j];

		EGLPNUM_TYPENAME_EGlpNumSet (c_dj, (double) mcnt);
		switch (lp->vtype[j])
		{
		case VFREE:
			porder[tfree] = j;
			perm[tfree] = tfree;
			EGLPNUM_TYPENAME_EGlpNumCopyFrac (qpenalty[tfree], lp->cz[j], cmax);
			EGLPNUM_TYPENAME_EGlpNumAddTo (qpenalty[tfree], c_dj);
			tfree++;
			break;

		case VLOWER:
		case VUPPER:
			porder[vd.nfree + tbndone] = j;
			perm[vd.nfree + tbndone] = tbndone;
			EGLPNUM_TYPENAME_EGlpNumCopyFrac (qpenalty[vd.nfree + tbndone], lp->cz[j], cmax);
			EGLPNUM_TYPENAME_EGlpNumAddTo (qpenalty[vd.nfree + tbndone], c_dj);
			if (lp->vtype[j] == VLOWER)
				EGLPNUM_TYPENAME_EGlpNumAddTo (qpenalty[vd.nfree + tbndone], lp->lz[j]);
			else
				EGLPNUM_TYPENAME_EGlpNumSubTo (qpenalty[vd.nfree + tbndone], lp->uz[j]);
			tbndone++;
			break;

		case VFIXED:
		case VBOUNDED:
			porder[tot1 + tbounded] = j;
			perm[tot1 + tbounded] = tbounded;
			EGLPNUM_TYPENAME_EGlpNumCopyFrac (qpenalty[tot1 + tbounded], lp->cz[j], cmax);
			EGLPNUM_TYPENAME_EGlpNumAddTo (qpenalty[tot1 + tbounded], lp->lz[j]);
			EGLPNUM_TYPENAME_EGlpNumSubTo (qpenalty[tot1 + tbounded], lp->uz[j]);
			EGLPNUM_TYPENAME_EGlpNumAddTo (qpenalty[tot1 + tbounded], c_dj);
			tbounded++;
			break;
		}
	}
#if BASIS_STATS > 0
	QSlog("bfree %d, bone %d, bbnd %d", tfree, tbndone, tbounded);
#endif

	EGLPNUM_TYPENAME_ILLutil_EGlpNum_perm_quicksort (perm, qpenalty, tfree);
	EGLPNUM_TYPENAME_ILLutil_EGlpNum_perm_quicksort (perm + vd.nfree, qpenalty + vd.nfree,
																	tbndone);
	EGLPNUM_TYPENAME_ILLutil_EGlpNum_perm_quicksort (perm + tot1, qpenalty + tot1, tbounded);

	tot2 = tfree + tbndone;
	for (i = 0; i < tbndone; i++)
	{
		perm[tfree + i] = perm[vd.nfree + i] + tfree;
		porder[tfree + i] = porder[vd.nfree + i];
	}
	for (i = 0; i < tbounded; i++)
	{
		perm[tot2 + i] = perm[tot1 + i] + tot2;
		porder[tot2 + i] = porder[tot1 + i];
	}
	tot2 += tbounded;

	nbelem =
		primal_col_select (lp, vstat, irow, rrow, unitcol, v, perm, porder, nbelem,
											 tot2);
	if (nbelem != lp->nrows)
	{
		QSlog("complain: incorrect final basis size");
		rval = E_SIMPLEX_ERROR;
		ILL_CLEANUP;
	}

CLEANUP:
	if (rval)
		EGLPNUM_TYPENAME_ILLbasis_free_basisinfo (lp);

	ILL_IFFREE (irow, int);
	ILL_IFFREE (rrow, int);

	EGLPNUM_TYPENAME_EGlpNumFreeArray (v);
	ILL_IFFREE (unitcol, int);
	ILL_IFFREE (icol, int);
	ILL_IFFREE (rcol, int);

	EGLPNUM_TYPENAME_EGlpNumFreeArray (dj);
	ILL_IFFREE (perm, int);
	ILL_IFFREE (porder, int);
	ILL_IFFREE (plen, int);

	EGLPNUM_TYPENAME_EGlpNumFreeArray (qpenalty);
	EGLPNUM_TYPENAME_EGlpNumClearVar (seldj);
	EGLPNUM_TYPENAME_EGlpNumClearVar (selv);
	EGLPNUM_TYPENAME_EGlpNumClearVar (c_dj);
	EGLPNUM_TYPENAME_EGlpNumClearVar (cmax);
	EGLPNUM_TYPENAME_ILLbasis_clear_vardata (&vd);
	EG_RETURN (rval);
}

static int set_basis_indices (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int *vstat)
{
	int i, b = 0, nb = 0;
	int vs;

	for (i = 0; i < lp->ncols; i++)
	{
		vs = vstat[i];
		lp->vstat[i] = vs;

		if (vs == STAT_BASIC)
		{
			lp->baz[b] = i;
			lp->vindex[i] = b;
			b++;
		}
		else if (vs == STAT_UPPER || vs == STAT_LOWER || vs == STAT_ZERO)
		{
			lp->nbaz[nb] = i;
			lp->vindex[i] = nb;
			nb++;
		}
		else
		{
			QSlog("Error in basis creation");
			return E_SIMPLEX_ERROR;
		}
	}
	if (b != lp->nrows)
	{
		QSlog("Error 2 in basis creation");
		return E_SIMPLEX_ERROR;
	}
	else if (nb != lp->nnbasic)
	{
		QSlog("Error 3 in basis creation");
		return E_SIMPLEX_ERROR;
	}
	return 0;
}

int EGLPNUM_TYPENAME_ILLbasis_get_initial (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int algorithm)
{
	int rval = 0;
	int *vstat = NULL;

	EGLPNUM_TYPENAME_ILLbasis_free_basisinfo (lp);
	EGLPNUM_TYPENAME_ILLbasis_init_basisinfo (lp);
	rval = EGLPNUM_TYPENAME_ILLbasis_build_basisinfo (lp);
	CHECKRVALG (rval, CLEANUP);

	ILL_SAFE_MALLOC (vstat, lp->ncols, int);

	if (algorithm == PRIMAL_SIMPLEX)
		rval = get_initial_basis1 (lp, vstat);
	else
		rval = get_initial_basis2 (lp, vstat);

	if (rval == E_SIMPLEX_ERROR)
	{
#ifdef HAVE_LIBZ
		EGioFile_t *f = EGioOpen ("bad.lp.gz", "w");
#else
#ifdef HAVE_LIBBZ2
		EGioFile_t *f = EGioOpen ("bad.lp.bz2", "w");
#else
		EGioFile_t *f = EGioOpen ("bad.lp", "w");
#endif
#endif
		int tval = EGLPNUM_TYPENAME_ILLwrite_lp_file (lp->O, f, NULL);
		if (tval)
		{
			QSlog("Error writing bad lp");
		}
		if (f != NULL)
			EGioClose (f);
	}
	CHECKRVALG (rval, CLEANUP);

	rval = set_basis_indices (lp, vstat);
	CHECKRVALG (rval, CLEANUP);
	lp->basisid = 0;

CLEANUP:
	ILL_IFFREE (vstat, int);

	EG_RETURN (rval);
}

static int choose_basis (
	int algorithm,
	EGLPNUM_TYPE pinf1,
	EGLPNUM_TYPE dinf1,
	EGLPNUM_TYPE pinf2,
	EGLPNUM_TYPE dinf2)
{
/* We changed the constant definitions outside here, the actual numbers are
 * asigned in lpdata.c. the values are as follows:
 * EGLPNUM_TYPENAME_CB_EPS = 0.001;
 * EGLPNUM_TYPENAME_CB_PRI_RLIMIT = 0.25;
 * EGLPNUM_TYPENAME_CB_INF_RATIO = 10.0; 
 * */
	int choice = 1;
	EGLPNUM_TYPE rp, rd;

	if (algorithm == PRIMAL_SIMPLEX)
	{
		EGLPNUM_TYPENAME_EGlpNumInitVar (rp);
		EGLPNUM_TYPENAME_EGlpNumInitVar (rd);
		EGLPNUM_TYPENAME_EGlpNumCopyDiff (rp, pinf1, pinf2);
		EGLPNUM_TYPENAME_EGlpNumCopyDiff (rd, dinf1, dinf2);
		if (EGLPNUM_TYPENAME_EGlpNumIsLeq (rp, EGLPNUM_TYPENAME_CB_EPS) && EGLPNUM_TYPENAME_EGlpNumIsLeq (rd, EGLPNUM_TYPENAME_CB_EPS))
			choice = 1;
		else
		{
			EGLPNUM_TYPENAME_EGlpNumSign (rp);
			EGLPNUM_TYPENAME_EGlpNumSign (rd);
			if (EGLPNUM_TYPENAME_EGlpNumIsLeq (rp, EGLPNUM_TYPENAME_CB_EPS) && EGLPNUM_TYPENAME_EGlpNumIsLeq (rd, EGLPNUM_TYPENAME_CB_EPS))
				choice = 2;
			else if (EGLPNUM_TYPENAME_EGlpNumIsLess (pinf1, pinf2) && EGLPNUM_TYPENAME_EGlpNumIsLess (dinf2, dinf1))
			{
				choice = 1;
				EGLPNUM_TYPENAME_EGlpNumCopyFrac (rp, pinf1, pinf2);
				EGLPNUM_TYPENAME_EGlpNumCopyFrac (rd, dinf2, dinf1);
				EGLPNUM_TYPENAME_EGlpNumMultTo (rd, EGLPNUM_TYPENAME_CB_INF_RATIO);
				if (EGLPNUM_TYPENAME_EGlpNumIsLess (EGLPNUM_TYPENAME_CB_PRI_RLIMIT, rp) && (EGLPNUM_TYPENAME_EGlpNumIsLess (rd, rp)))
					choice = 2;
			}
			else if (EGLPNUM_TYPENAME_EGlpNumIsLess (pinf2, pinf1) && EGLPNUM_TYPENAME_EGlpNumIsLess (dinf1, dinf2))
			{
				choice = 2;
				EGLPNUM_TYPENAME_EGlpNumCopyFrac (rp, pinf2, pinf1);
				EGLPNUM_TYPENAME_EGlpNumCopyFrac (rd, dinf1, dinf2);
				EGLPNUM_TYPENAME_EGlpNumMultTo (rd, EGLPNUM_TYPENAME_CB_INF_RATIO);
				if (EGLPNUM_TYPENAME_EGlpNumIsLess (EGLPNUM_TYPENAME_CB_PRI_RLIMIT, rp) && EGLPNUM_TYPENAME_EGlpNumIsLess (rd, rp))
					choice = 1;
			}
			else
				choice = 1;
		}
		EGLPNUM_TYPENAME_EGlpNumClearVar (rp);
		EGLPNUM_TYPENAME_EGlpNumClearVar (rd);
	}
	ILL_IFTRACE ("%s:%d\n", __func__, choice);
	return choice;
}

int EGLPNUM_TYPENAME_ILLbasis_get_cinitial (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int algorithm)
{
	int rval = 0;
	int *vstat1 = NULL;
	int *vstat2 = NULL;
	int singular;
	int choice = 0;

#if BASIS_STATS > 0
	int i, nz1 = 0, nz2 = 0;
#endif
	EGLPNUM_TYPE pinf1, pinf2, dinf1, dinf2;
	EGLPNUM_TYPENAME_feas_info fi;

	EGLPNUM_TYPENAME_EGlpNumInitVar (pinf1);
	EGLPNUM_TYPENAME_EGlpNumInitVar (pinf2);
	EGLPNUM_TYPENAME_EGlpNumInitVar (dinf1);
	EGLPNUM_TYPENAME_EGlpNumInitVar (dinf2);
	EGLPNUM_TYPENAME_EGlpNumInitVar (fi.totinfeas);

	EGLPNUM_TYPENAME_ILLbasis_free_basisinfo (lp);
	EGLPNUM_TYPENAME_ILLbasis_init_basisinfo (lp);
	rval = EGLPNUM_TYPENAME_ILLbasis_build_basisinfo (lp);
	CHECKRVALG (rval, CLEANUP);

	ILL_SAFE_MALLOC (vstat1, lp->ncols, int);
	ILL_SAFE_MALLOC (vstat2, lp->ncols, int);

	if (algorithm != PRIMAL_SIMPLEX)
	{
		rval = get_initial_basis2 (lp, vstat2);
		CHECKRVALG (rval, CLEANUP);
		rval = set_basis_indices (lp, vstat2);
		lp->basisid = 0;
		ILL_CLEANUP;
	}

	rval = get_initial_basis1 (lp, vstat1);
	CHECKRVALG (rval, CLEANUP);
	rval = get_initial_basis2 (lp, vstat2);
	CHECKRVALG (rval, CLEANUP);
	lp->basisid = 0;

	/* handle first basis */
	rval = set_basis_indices (lp, vstat1);
	CHECKRVALG (rval, CLEANUP);
#if BASIS_STATS > 0
	for (i = 0; i < lp->nrows; i++)
		nz1 += lp->matcnt[lp->baz[i]];
#endif
	rval = EGLPNUM_TYPENAME_ILLbasis_factor (lp, &singular);
	if (singular)
		MESSAGE (__QS_SB_VERB, "Singular Basis found!");
	CHECKRVALG (rval, CLEANUP);

	EGLPNUM_TYPENAME_ILLfct_compute_piz (lp);
	EGLPNUM_TYPENAME_ILLfct_compute_dz (lp);
	EGLPNUM_TYPENAME_ILLfct_dual_adjust (lp, EGLPNUM_TYPENAME_zeroLpNum);
	EGLPNUM_TYPENAME_ILLfct_compute_xbz (lp);

	EGLPNUM_TYPENAME_ILLfct_check_pfeasible (lp, &fi, lp->tol->pfeas_tol);
	EGLPNUM_TYPENAME_ILLfct_check_dfeasible (lp, &fi, lp->tol->dfeas_tol);
	EGLPNUM_TYPENAME_EGlpNumCopy (pinf1, lp->pinfeas);
	EGLPNUM_TYPENAME_EGlpNumCopy (dinf1, lp->dinfeas);
	/*
	 * EGLPNUM_TYPENAME_ILLfct_compute_pobj (lp);  obj1p = lp->objval;
	 * EGLPNUM_TYPENAME_ILLfct_compute_dobj (lp);  obj1d = lp->objval;
	 */

	/* handle second basis */
	rval = set_basis_indices (lp, vstat2);
	CHECKRVALG (rval, CLEANUP);
#if BASIS_STATS > 0
	for (i = 0; i < lp->nrows; i++)
		nz2 += lp->matcnt[lp->baz[i]];
#endif
	rval = EGLPNUM_TYPENAME_ILLbasis_factor (lp, &singular);
	if (singular)
		MESSAGE (__QS_SB_VERB, "Singular Basis found!");
	CHECKRVALG (rval, CLEANUP);

	EGLPNUM_TYPENAME_ILLfct_compute_piz (lp);
	EGLPNUM_TYPENAME_ILLfct_compute_dz (lp);
	EGLPNUM_TYPENAME_ILLfct_dual_adjust (lp, EGLPNUM_TYPENAME_zeroLpNum);
	EGLPNUM_TYPENAME_ILLfct_compute_xbz (lp);

	EGLPNUM_TYPENAME_ILLfct_check_pfeasible (lp, &fi, lp->tol->pfeas_tol);
	EGLPNUM_TYPENAME_ILLfct_check_dfeasible (lp, &fi, lp->tol->dfeas_tol);
	EGLPNUM_TYPENAME_EGlpNumCopy (pinf2, lp->pinfeas);
	EGLPNUM_TYPENAME_EGlpNumCopy (dinf2, lp->dinfeas);

#if BASIS_STATS > 0
	QSlog("b1: nz %d pinf %.2f dinf %.2f", nz1, EGLPNUM_TYPENAME_EGlpNumToLf (pinf1),
							EGLPNUM_TYPENAME_EGlpNumToLf (dinf1));
	QSlog("b2: nz %d pinf %.2f dinf %.2f", nz2, EGLPNUM_TYPENAME_EGlpNumToLf (pinf2),
					EGLPNUM_TYPENAME_EGlpNumToLf (dinf2));
#endif
	choice = choose_basis (algorithm, pinf1, dinf1, pinf2, dinf2);
	if (choice == 1)
	{
		lp->fbasisid = -1;
		rval = set_basis_indices (lp, vstat1);
		CHECKRVALG (rval, CLEANUP);
	}

CLEANUP:
	if (rval == E_SIMPLEX_ERROR)
	{
#ifdef HAVE_LIBZ
		EGioFile_t *fil = EGioOpen ("bad.lp.gz", "w");
#else
#ifdef HAVE_LIBBZ2
		EGioFile_t *fil = EGioOpen ("bad.lp.bz2", "w");
#else
		EGioFile_t *fil = EGioOpen ("bad.lp", "w");
#endif
#endif
		int tval = EGLPNUM_TYPENAME_ILLwrite_lp_file (lp->O, fil, NULL);

		if (tval)
		{
			QSlog("Error writing bad lp");
		}
		if (fil != NULL)
			EGioClose (fil);
	}
	ILL_IFFREE (vstat1, int);
	ILL_IFFREE (vstat2, int);

	EGLPNUM_TYPENAME_EGlpNumClearVar (pinf1);
	EGLPNUM_TYPENAME_EGlpNumClearVar (pinf2);
	EGLPNUM_TYPENAME_EGlpNumClearVar (dinf1);
	EGLPNUM_TYPENAME_EGlpNumClearVar (dinf2);
	EGLPNUM_TYPENAME_EGlpNumClearVar (fi.totinfeas);
	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLbasis_factor (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int *singular)
{
	int rval = 0;
	int i;
	int eindex;
	int lindex;
	int ltype;
	int lvstat;
	int nsing = 0;
	int *singr = 0;
	int *singc = 0;

	*singular = 0;
	do
	{
		if (lp->f)
		{
			EGLPNUM_TYPENAME_ILLfactor_free_factor_work (lp->f);
		}
		else
		{
			ILL_SAFE_MALLOC (lp->f, 1, EGLPNUM_TYPENAME_factor_work);
			EGLPNUM_TYPENAME_EGlpNumInitVar (lp->f->fzero_tol);
			EGLPNUM_TYPENAME_EGlpNumInitVar (lp->f->szero_tol);
			EGLPNUM_TYPENAME_EGlpNumInitVar (lp->f->partial_tol);
			EGLPNUM_TYPENAME_EGlpNumInitVar (lp->f->maxelem_orig);
			EGLPNUM_TYPENAME_EGlpNumInitVar (lp->f->maxelem_factor);
			EGLPNUM_TYPENAME_EGlpNumInitVar (lp->f->maxelem_cur);
			EGLPNUM_TYPENAME_EGlpNumInitVar (lp->f->partial_cur);
			EGLPNUM_TYPENAME_ILLfactor_init_factor_work (lp->f);
		}
		rval = EGLPNUM_TYPENAME_ILLfactor_create_factor_work (lp->f, lp->O->nrows);
		CHECKRVALG (rval, CLEANUP);

		rval = EGLPNUM_TYPENAME_ILLfactor (lp->f, lp->baz, lp->matbeg, lp->matcnt,
											lp->matind, lp->matval, &nsing, &singr, &singc);
		CHECKRVALG (rval, CLEANUP);

		if (nsing != 0)
		{
			*singular = 1;
			MESSAGE (__QS_SB_VERB, "Found singular basis!");
			for (i = 0; i < nsing; i++)
			{
				eindex = lp->vindex[lp->O->rowmap[singr[i]]];
				lindex = singc[i];
				ltype = lp->vtype[lp->baz[lindex]];

				if (ltype == VBOUNDED || ltype == VLOWER || ltype == VARTIFICIAL)
					lvstat = STAT_LOWER;
				else if (ltype == VUPPER)
					lvstat = STAT_UPPER;
				else
					lvstat = STAT_ZERO;

				EGLPNUM_TYPENAME_ILLfct_update_basis_info (lp, eindex, lindex, lvstat);
				lp->basisid++;
			}
			ILL_IFFREE (singr, int);
			ILL_IFFREE (singc, int);
		}

	} while (nsing != 0);

	lp->fbasisid = lp->basisid;

CLEANUP:
	ILL_IFFREE (singr, int);
	ILL_IFFREE (singc, int);

	if (rval)
		QSlog("Error: unknown in %s", __func__);
	EG_RETURN (rval);
}

int EGLPNUM_TYPENAME_ILLbasis_refactor (
	EGLPNUM_TYPENAME_lpinfo * lp)
{
	int sing = 0;
	int rval = 0;

	rval = EGLPNUM_TYPENAME_ILLbasis_factor (lp, &sing);
	if (sing)
	{
		MESSAGE (__QS_SB_VERB, "Singular Basis found!");
		rval = QS_LP_CHANGE_PREC;
		return rval;
	}
	EG_RETURN (rval);
}

void EGLPNUM_TYPENAME_ILLbasis_column_solve (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_svector * rhs,
	EGLPNUM_TYPENAME_svector * soln)
{
	EGLPNUM_TYPENAME_ILLfactor_ftran (lp->f, rhs, soln);
}

void EGLPNUM_TYPENAME_ILLbasis_column_solve_update (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_svector * rhs,
	EGLPNUM_TYPENAME_svector * upd,
	EGLPNUM_TYPENAME_svector * soln)
{
	EGLPNUM_TYPENAME_ILLfactor_ftran_update (lp->f, rhs, upd, soln);
}

void EGLPNUM_TYPENAME_ILLbasis_row_solve (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_svector * rhs,
	EGLPNUM_TYPENAME_svector * soln)
{
	EGLPNUM_TYPENAME_ILLfactor_btran (lp->f, rhs, soln);
}

int EGLPNUM_TYPENAME_ILLbasis_update (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_svector * y,
	int lindex,
	int *refactor,
	int *singular)
{
#if 0														/* To always refactor, change 0 to 1 */
	*refactor = 1;
	return EGLPNUM_TYPENAME_ILLbasis_factor (lp, singular);
#else

	int rval = 0;

	*refactor = 0;
	rval = EGLPNUM_TYPENAME_ILLfactor_update (lp->f, y, lindex, refactor);
	if (rval == E_FACTOR_BLOWUP || rval == E_UPDATE_SINGULAR_ROW
			|| rval == E_UPDATE_SINGULAR_COL)
	{
/* Bico - comment out for dist
       QSlog("Warning: numerically bad basis in EGLPNUM_TYPENAME_ILLfactor_update");
*/
		*refactor = 1;
		rval = 0;
	}
	if (rval == E_UPDATE_NOSPACE)
	{
		*refactor = 1;
		rval = 0;
	}

	if (*refactor)
	{
		rval = EGLPNUM_TYPENAME_ILLbasis_factor (lp, singular);
		if (*singular)
			MESSAGE (__QS_SB_VERB, "Singular Basis found!");
	}
	if (rval)
	{
		EGioFile_t *eout = 0;
		int tval;

		QSlog("write bad lp to factor.lp");
#ifdef HAVE_LIBZ
		eout = EGioOpen ("factor.lp.gz", "w");
#else
#ifdef HAVE_LIBBZ2
		eout = EGioOpen ("factor.lp.bz2", "w");
#else
		eout = EGioOpen ("factor.lp", "w");
#endif
#endif
		if (!eout)
		{
			QSlog("could not open file to write bad factor lp");
		}
		else
		{
			tval = EGLPNUM_TYPENAME_ILLwrite_lp_file (lp->O, eout, NULL);
			if (tval)
			{
				QSlog("error while writing bad factor lp");
			}
			EGioClose (eout);
		}

		QSlog("write bad basis to factor.bas");
		tval = EGLPNUM_TYPENAME_ILLlib_writebasis (lp, 0, "factor.bas");
		if (tval)
		{
			QSlog("error while writing factor basis");
		}
	}

	EG_RETURN (rval);
#endif
}
