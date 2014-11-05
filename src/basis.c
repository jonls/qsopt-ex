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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "qs_config.h"
#include "config.h"
#include "sortrus.h"
#include "iqsutil.h"
#include "lpdefs.h"
#include "qstruct.h"
#include "qsopt.h"
#include "basis.h"
#include "fct.h"
#include "lp.h"
#include "lib.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

//#define DJZERO_TOLER PFEAS_TOLER
#define BASIS_STATS 0
//#define BASIS_DEBUG 10
#define BASIS_DEBUG 0

void ILLbasis_init_vardata (
	var_data * vd)
{
	memset (vd, 0, sizeof (var_data));
	EGlpNumInitVar (vd->cmax);
}

void ILLbasis_clear_vardata (
	var_data * vd)
{
	EGlpNumClearVar (vd->cmax);
	memset (vd, 0, sizeof (var_data));
}

static void get_var_info (
	lpinfo * lp,
	var_data * v);

static int init_slack_basis (
	lpinfo * lp,
	int *vstat,
	int *irow,
	int *rrow,
	int *unitcol,
	int *icol,
	int *rcol),
  get_initial_basis1 (
	lpinfo * lp,
	int *vstat),
  get_initial_basis2 (
	lpinfo * lp,
	int *vstat),
  set_basis_indices (
	lpinfo * lp,
	int *vstat),
  choose_basis (
	int algorithm,
	EGlpNum_t pinf1,
	EGlpNum_t dinf1,
	EGlpNum_t pinf2,
	EGlpNum_t dinf2);

void ILLbasis_init_basisinfo (
	lpinfo * lp)
{
	lp->baz = 0;
	lp->nbaz = 0;
	lp->vstat = 0;
	lp->vindex = 0;
	lp->f = 0;
}

void ILLbasis_free_basisinfo (
	lpinfo * lp)
{
	ILL_IFFREE (lp->baz, int);
	ILL_IFFREE (lp->nbaz, int);
	ILL_IFFREE (lp->vstat, int);
	ILL_IFFREE (lp->vindex, int);

	if (lp->f)
	{
		ILLfactor_free_factor_work (lp->f);
		EGlpNumClearVar (lp->f->fzero_tol);
		EGlpNumClearVar (lp->f->szero_tol);
		EGlpNumClearVar (lp->f->partial_tol);
		EGlpNumClearVar (lp->f->maxelem_orig);
		EGlpNumClearVar (lp->f->maxelem_factor);
		EGlpNumClearVar (lp->f->maxelem_cur);
		EGlpNumClearVar (lp->f->partial_cur);
		ILL_IFFREE (lp->f, factor_work);
	}
}

int ILLbasis_build_basisinfo (
	lpinfo * lp)
{
	int rval = 0;

	ILL_SAFE_MALLOC (lp->baz, lp->O->nrows, int);
	ILL_SAFE_MALLOC (lp->nbaz, lp->O->ncols, int);
	ILL_SAFE_MALLOC (lp->vstat, lp->O->ncols, int);
	ILL_SAFE_MALLOC (lp->vindex, lp->O->ncols, int);

	lp->fbasisid = -1;

CLEANUP:
	if (rval)
		ILLbasis_free_basisinfo (lp);
	EG_RETURN (rval);
}

int ILLbasis_load (
	lpinfo * lp,
	ILLlp_basis * B)
{
	int rval = 0;
	char *cstat = B->cstat;
	char *rstat = B->rstat;
	int *structmap = lp->O->structmap;
	int *rowmap = lp->O->rowmap;
	char *sense = lp->O->sense;
	int i, j, ncols = lp->O->ncols, nrows = lp->O->nrows, nstruct = lp->O->nstruct;
	int basic = 0, nonbasic = 0;

	ILLbasis_free_basisinfo (lp);
	ILLbasis_init_basisinfo (lp);
	rval = ILLbasis_build_basisinfo (lp);
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
				fprintf (stderr, "unknown col basis stat 1: %c\n", cstat[i]);
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
					fprintf (stderr, "unknown range basis stat 2\n");
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
				fprintf (stderr, "unknown row basis stat 3\n");
				rval = 1;
				goto CLEANUP;
			}
		}
	}

	if (basic + nonbasic != ncols)
	{
		fprintf (stderr, "error in counts in ILLopt_load_basis\n");
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

int ILLbasis_tableau_row (
	lpinfo * lp,
	int row,
	EGlpNum_t * brow,
	EGlpNum_t * trow,
	EGlpNum_t * rhs,
	int strict)
{
	int rval = 0;
	int i;
	int singular = 0;
	int indx;
	EGlpNum_t coef;
	EGlpNum_t sum;
	svector z, zA;

	EGlpNumInitVar (coef);
	EGlpNumInitVar (sum);
	EGlpNumZero (sum);

	ILLsvector_init (&z);
	ILLsvector_init (&zA);

	if (lp->basisid == -1)
	{
		fprintf (stderr, "ILLbasis_tableau_row: no basis\n");
		rval = E_GENERAL_ERROR;
		ILL_CLEANUP;
	}
	if (lp->fbasisid != lp->basisid)
	{															/* Needs to be changed */
		rval = ILLbasis_factor (lp, &singular);
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
		fprintf (stderr, "No array for basis inverse row\n");
		rval = E_GENERAL_ERROR;
		ILL_CLEANUP;
	}

	rval = ILLsvector_alloc (&z, lp->nrows);
	CHECKRVALG (rval, CLEANUP);
	ILLfct_compute_zz (lp, &z, row);

	for (i = 0; i < lp->O->nrows; i++)
		EGlpNumZero (brow[i]);
	for (i = 0; i < z.nzcnt; i++)
	{
		indx = z.indx[i];
		EGlpNumCopy (coef, z.coef[i]);
		EGlpNumCopy (brow[indx], coef);
		EGlpNumAddInnProdTo (sum, coef, lp->bz[indx]);
	}

	if (rhs != NULL)
		EGlpNumCopy (*rhs, sum);
	if (trow != NULL)
	{
		if (!strict)
		{
			rval = ILLsvector_alloc (&zA, lp->ncols);
			if (rval)
				ILL_CLEANUP;
			ILL_IFTRACE ("%s:\n", __func__);
			rval = ILLfct_compute_zA (lp, &z, &zA);
			CHECKRVALG (rval, CLEANUP);

			for (i = 0; i < lp->ncols; i++)
				EGlpNumZero (trow[i]);
			for (i = 0; i < zA.nzcnt; i++)
				EGlpNumCopy (trow[lp->nbaz[zA.indx[i]]], zA.coef[i]);
			EGlpNumOne (trow[lp->baz[row]]);
		}
		else
		{
			ILLfct_compute_vA (lp, &z, trow);
		}
	}

#if BASIS_DEBUG > 0
	if (rhs != NULL && trow != NULL)
	{
		EGlpNum_t *tr = NULL;

		EGlpNumZero (sum);
		if (strict)
			tr = trow;
		else
		{
			tr = EGlpNumAllocArray (lp->ncols);
			ILLfct_compute_vA (lp, &z, tr);
		}
		for (i = 0; i < lp->nrows; i++)
			if (EGlpNumIsGreatZero (tr[lp->baz[i]]))
				EGlpNumAddTo (sum, tr[lp->baz[i]]);
			else
				EGlpNumSubTo (sum, tr[lp->baz[i]]);
		EGlpNumCopy (coef, oneLpNum);
		EGlpNumSubTo (coef, sum);
		if (EGlpNumIsLessZero (coef))
			EGlpNumSign (coef);
		if (EGlpNumIsLess (PIVZ_TOLER, coef))
			fprintf (stderr, "tableau: bas computed = %.12f\n", EGlpNumToLf (sum));
		if (!strict)
			EGlpNumFreeArray (tr);
#if BASIS_DEBUG > 1
		EGlpNumZero (sum);
		for (i = 0; i < lp->ncols; i++)
		{
			if (lp->vstat[i] == STAT_BASIC)
				EGlpNumAddInnProdTo (sum, lp->xbz[lp->vindex[i]], trow[i]);
			else if (lp->vstat[i] == STAT_UPPER)
				EGlpNumAddInnProdTo (sum, lp->uz[i], trow[i]);
			else if (lp->vstat[i] == STAT_LOWER)
				EGlpNumAddInnProdTo (sum, lp->lz[i], trow[i]);
		}
		EGlpNumSet (coef, 1e-10);
		if (EGlpNumIsNeq (sum, *rhs, coef))
			fprintf (stderr, "tableau rhs = %.9f, computed = %.9f\n",
							 EGlpNumToLf (*rhs), EGlpNumToLf (sum));
#endif
	}
#endif

CLEANUP:
	ILLsvector_free (&z);
	ILLsvector_free (&zA);
	EGlpNumClearVar (coef);
	EGlpNumClearVar (sum);
	return rval;
}

static void get_var_info (
	lpinfo * lp,
	var_data * v)
{
	int i = 0;

	v->nartif = 0;
	v->nslacks = 0;
	v->nfree = 0;
	v->nbndone = 0;
	v->nbounded = 0;
	v->nfixed = 0;
	EGlpNumCopy (v->cmax, NINFTY);

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
		EGlpNumSetToMaxAbs (v->cmax, lp->cz[i]);
	}

#if BASIS_STATS > 0
	printf ("cols = %d, acols = %d, total  = %d, nrows = %d, nlog = %d\n",
					lp->ncols, lp->ncols - lp->nrows,
					v->nartif + v->nfree + v->nslacks + v->nbndone + v->nbounded,
					lp->nrows, v->nartif + v->nslacks);
#endif
}

static int init_slack_basis (
	lpinfo * lp,
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
			if (fabs (EGlpNumToLf (lp->lz[j])) < fabs (EGlpNumToLf (lp->uz[j])))
				vstat[j] = STAT_LOWER;
			else
				vstat[j] = STAT_UPPER;
		}
	}
	return nslacks;
}

static int primal_col_select (
	lpinfo * lp,
	int *vstat,
	int *irow,
	int *rrow,
	int *unitcol,
	EGlpNum_t * v,
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
	EGlpNum_t *matval = lp->matval;
	EGlpNum_t alpha, val, maxelem;

	EGlpNumInitVar (alpha);
	EGlpNumInitVar (val);
	EGlpNumInitVar (maxelem);

	for (k = 0; k < pcols; k++)
	{
		j = porder[perm[k]];
		mcnt = matcnt[j];
		mbeg = matbeg[j];

		EGlpNumCopy (alpha, NINFTY);
		EGlpNumCopy (maxelem, NINFTY);

		for (i = 0; i < mcnt; i++)
		{
			EGlpNumCopyAbs (val, matval[mbeg + i]);
			if (EGlpNumIsLess (maxelem, val))
				EGlpNumCopy (maxelem, val);
			if (rrow[matind[mbeg + i]] == 0 && EGlpNumIsLess (alpha, val))
			{
				EGlpNumCopy (alpha, val);
				r = matind[mbeg + i];
			}
		}
		EGlpNumCopy (val, maxelem);
		EGlpNumMultTo (val, PARAM_IBASIS_RPIVOT);
		if (EGlpNumIsLess (val, alpha))
		{
			vstat[j] = STAT_BASIC;
			nbelem++;
			irow[r] = 1;
			EGlpNumCopy (v[r], alpha);
			for (i = 0; i < mcnt; i++)
				if (EGlpNumIsNeqqZero (matval[mbeg + i]))
					rrow[matind[mbeg + i]]++;
		}
		else
		{
			EGlpNumCopy (alpha, NINFTY);
			for (i = 0; i < mcnt; i++)
			{
				tr = matind[mbeg + i];
				EGlpNumCopyAbs (val, matval[mbeg + i]);
				EGlpNumDivTo (val, PARAM_IBASIS_RTRIANG);
				if (EGlpNumIsNeqq (v[tr], INFTY) && EGlpNumIsLess (v[tr], val))
				{
					EGlpNumZero (alpha);
					break;
				}
				EGlpNumCopyAbs (val, matval[mbeg + i]);
				if (irow[tr] == 0 && EGlpNumIsLess (alpha, val))
				{
					EGlpNumCopy (alpha, val);
					r = tr;
				}
			}
			if (EGlpNumIsNeqqZero (alpha) && EGlpNumIsNeqq (alpha, NINFTY))
			{
				vstat[j] = STAT_BASIC;
				nbelem++;
				irow[r] = 1;
				EGlpNumCopy (v[r], alpha);
				for (i = 0; i < mcnt; i++)
					if (EGlpNumIsNeqqZero (matval[mbeg + i]))
						rrow[matind[mbeg + i]]++;
			}
		}
	}
#if BASIS_STATS > 0
	printf ("nartifs = %d\n", lp->nrows - nbelem);
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
					fprintf (stderr, "Error: Not enough artificials\n");
					return -1;
				}
			}
		}
	}
	EGlpNumClearVar (alpha);
	EGlpNumClearVar (val);
	EGlpNumClearVar (maxelem);
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
	lpinfo * lp,
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
	EGlpNum_t cmax;
	EGlpNum_t *v = NULL;
	EGlpNum_t *qpenalty = NULL;
	var_data vd;

	ILLbasis_init_vardata (&vd);
	EGlpNumInitVar (cmax);

	get_var_info (lp, &vd);
	if (!EGlpNumIsNeqqZero (vd.cmax))
		EGlpNumOne (cmax);
	else
	{
		EGlpNumCopy (cmax, vd.cmax);
		EGlpNumMultUiTo (cmax, 1000);
	}

	ILL_SAFE_MALLOC (irow, lp->nrows, int);
	ILL_SAFE_MALLOC (rrow, lp->nrows, int);

	v = EGlpNumAllocArray (lp->nrows);
	ILL_SAFE_MALLOC (unitcol, lp->nrows, int);

	for (i = 0; i < lp->nrows; i++)
	{
		unitcol[i] = -1;
		EGlpNumCopy (v[i], INFTY);
		irow[i] = 0;
		rrow[i] = 0;
	}

	nslacks = init_slack_basis (lp, vstat, irow, rrow, unitcol, NULL, NULL);
	if (nslacks != vd.nslacks)
	{
		printf ("complain: incorrect basis info(slacks)\n");
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
					fprintf (stderr, "Error: Not enough artificials\n");
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

	qpenalty = EGlpNumAllocArray (tot2);

	for (j = 0; j < lp->ncols; j++)
	{
		if (vstat[j] == STAT_BASIC)
			continue;

		switch (lp->vtype[j])
		{
		case VFREE:
			porder[tfree] = j;
			perm[tfree] = tfree;
			EGlpNumCopyFrac (qpenalty[tfree], lp->cz[j], cmax);
			tfree++;
			break;

		case VLOWER:
		case VUPPER:
			porder[vd.nfree + tbndone] = j;
			perm[vd.nfree + tbndone] = tbndone;
			EGlpNumCopyFrac (qpenalty[vd.nfree + tbndone], lp->cz[j], cmax);
			if (lp->vtype[j] == VLOWER)
				EGlpNumAddTo (qpenalty[vd.nfree + tbndone], lp->lz[j]);
			else
				EGlpNumSubTo (qpenalty[vd.nfree + tbndone], lp->uz[j]);
			tbndone++;
			break;

		case VFIXED:
		case VBOUNDED:
			porder[tot1 + tbounded] = j;
			perm[tot1 + tbounded] = tbounded;
			EGlpNumCopyFrac (qpenalty[tot1 + tbndone], lp->cz[j], cmax);
			EGlpNumAddTo (qpenalty[tot1 + tbndone], lp->lz[j]);
			EGlpNumSubTo (qpenalty[tot1 + tbndone], lp->uz[j]);
			tbounded++;
			break;
		}
	}
	if (tfree != vd.nfree || tbndone != vd.nbndone || tbounded != vd.nbounded)
	{
		printf ("complain: incorrect basis info \n");
		rval = E_SIMPLEX_ERROR;
		ILL_CLEANUP;
	}

	ILLutil_EGlpNum_perm_quicksort (perm, qpenalty, vd.nfree);
	ILLutil_EGlpNum_perm_quicksort (perm + vd.nfree, qpenalty + vd.nfree,
																	vd.nbndone);
	ILLutil_EGlpNum_perm_quicksort (perm + tot1, qpenalty + tot1, vd.nbounded);

	for (i = 0; i < vd.nbndone; i++)
		perm[vd.nfree + i] += vd.nfree;
	for (i = 0; i < vd.nbounded; i++)
		perm[tot1 + i] += tot1;

	nbelem =
		primal_col_select (lp, vstat, irow, rrow, unitcol, v, perm, porder, nbelem,
											 tot2);
	if (nbelem != lp->nrows)
	{
		printf ("complain: incorrect final basis size\n");
		rval = E_SIMPLEX_ERROR;
		ILL_CLEANUP;
	}

CLEANUP:
	EGlpNumClearVar (cmax);
	if (rval)
		ILLbasis_free_basisinfo (lp);
	ILL_IFFREE (irow, int);
	ILL_IFFREE (rrow, int);

	EGlpNumFreeArray (v);
	ILL_IFFREE (perm, int);
	ILL_IFFREE (porder, int);
	ILL_IFFREE (unitcol, int);

	EGlpNumFreeArray (qpenalty);
	ILLbasis_clear_vardata (&vd);
	EG_RETURN (rval);
}

static int get_initial_basis2 (
	lpinfo * lp,
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
	EGlpNum_t *v = NULL;
	EGlpNum_t *qpenalty = NULL;
	int col = 0, s_i = 0, selc = 0;
	int *icol = NULL, *rcol = NULL;
	int *plen = NULL;
	EGlpNum_t *dj = NULL;
	var_data vd;
	EGlpNum_t seldj;
	EGlpNum_t selv;
	EGlpNum_t c_dj;
	EGlpNum_t cmax;

	EGlpNumInitVar (seldj);
	EGlpNumInitVar (selv);
	EGlpNumInitVar (c_dj);
	EGlpNumInitVar (cmax);
	EGlpNumZero (c_dj);
	EGlpNumZero (selv);
	EGlpNumZero (seldj);
	ILLbasis_init_vardata (&vd);

	get_var_info (lp, &vd);

	ILL_SAFE_MALLOC (irow, lp->nrows, int);
	ILL_SAFE_MALLOC (rrow, lp->nrows, int);

	v = EGlpNumAllocArray (lp->nrows);
	ILL_SAFE_MALLOC (unitcol, lp->nrows, int);
	ILL_SAFE_MALLOC (icol, lp->ncols, int);
	ILL_SAFE_MALLOC (rcol, lp->ncols, int);

	dj = EGlpNumAllocArray (lp->ncols);

	for (i = 0; i < lp->nrows; i++)
	{
		unitcol[i] = -1;
		EGlpNumCopy (v[i], INFTY);
		irow[i] = 0;
		rrow[i] = 0;
	}
	/* assign all d_j */
	for (i = 0; i < lp->ncols; i++)
	{
		icol[i] = 0;
		rcol[i] = 0;
		EGlpNumCopy (dj[i], lp->cz[i]);
	}

	nslacks = init_slack_basis (lp, vstat, irow, rrow, unitcol, icol, rcol);
	if (nslacks != vd.nslacks)
	{
		printf ("complain: incorrect basis info\n");
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

	qpenalty = EGlpNumAllocArray (lp->ncols);

	/* find all unit rows and record lengths */
	for (i = 0; i < lp->nrows; i++)
	{
		if (irow[i] != 1)
		{
			rbeg = lp->rowbeg[i];
			rcnt = lp->rowcnt[i];
			for (j = 0; j < rcnt; j++)
			{
				EGlpNumCopyAbs (cmax, lp->rowval[rbeg + j]);
				if (EGlpNumIsNeqq (cmax, oneLpNum))
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
		EGlpNumCopy (seldj, INFTY);
		EGlpNumZero (selv);

		/* for every row s_i, compute min {d_j : d_j <0 , j is u or l or fr} */
		for (j = 0; j < rcnt; j++)
		{
			col = lp->rowind[rbeg + j];
			if (rcol[col] == 1)
				break;
			if (EGlpNumIsLessZero (dj[col]))
			{
				if (EGlpNumIsLess (dj[col], seldj))
				{
					selc = col;
					EGlpNumCopy (seldj, dj[col]);
					EGlpNumCopy (selv, lp->rowval[rbeg + j]);
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
			EGlpNumCopyFrac (c_dj, dj[selc], selv);
			vstat[selc] = STAT_BASIC;
			for (j = 0; j < rcnt; j++)
			{
				col = lp->rowind[rbeg + j];
				EGlpNumSubInnProdTo (dj[col], lp->rowval[rbeg + j], c_dj);
				rcol[col] = 1;
			}
		}
	}
#if BASIS_STATS > 0
	printf ("unit rows = %d\n", s_i);
	printf ("nslacks %d, unit rows selected = %d\n", nslacks, nbelem - nslacks);
#endif
	/* now go through remaining cols with dj = 0 */
	tot1 = vd.nfree + vd.nbndone;

	if (!EGlpNumIsNeqqZero (vd.cmax))
		EGlpNumOne (cmax);
	else
	{
		EGlpNumCopy (cmax, vd.cmax);
		EGlpNumMultUiTo (cmax, 1000);
	}
	for (j = 0; j < lp->ncols; j++)
	{
		if (vstat[j] == STAT_BASIC)
			continue;
		if (icol[j] == 1 || EGlpNumIsNeqZero (dj[j], BD_TOLER))
			continue;
		mcnt = lp->matcnt[j];

		EGlpNumSet (c_dj, (double) mcnt);
		switch (lp->vtype[j])
		{
		case VFREE:
			porder[tfree] = j;
			perm[tfree] = tfree;
			EGlpNumCopyFrac (qpenalty[tfree], lp->cz[j], cmax);
			EGlpNumAddTo (qpenalty[tfree], c_dj);
			tfree++;
			break;

		case VLOWER:
		case VUPPER:
			porder[vd.nfree + tbndone] = j;
			perm[vd.nfree + tbndone] = tbndone;
			EGlpNumCopyFrac (qpenalty[vd.nfree + tbndone], lp->cz[j], cmax);
			EGlpNumAddTo (qpenalty[vd.nfree + tbndone], c_dj);
			if (lp->vtype[j] == VLOWER)
				EGlpNumAddTo (qpenalty[vd.nfree + tbndone], lp->lz[j]);
			else
				EGlpNumSubTo (qpenalty[vd.nfree + tbndone], lp->uz[j]);
			tbndone++;
			break;

		case VFIXED:
		case VBOUNDED:
			porder[tot1 + tbounded] = j;
			perm[tot1 + tbounded] = tbounded;
			EGlpNumCopyFrac (qpenalty[tot1 + tbounded], lp->cz[j], cmax);
			EGlpNumAddTo (qpenalty[tot1 + tbounded], lp->lz[j]);
			EGlpNumSubTo (qpenalty[tot1 + tbounded], lp->uz[j]);
			EGlpNumAddTo (qpenalty[tot1 + tbounded], c_dj);
			tbounded++;
			break;
		}
	}
#if BASIS_STATS > 0
	printf ("bfree %d, bone %d, bbnd %d\n", tfree, tbndone, tbounded);
#endif

	ILLutil_EGlpNum_perm_quicksort (perm, qpenalty, tfree);
	ILLutil_EGlpNum_perm_quicksort (perm + vd.nfree, qpenalty + vd.nfree,
																	tbndone);
	ILLutil_EGlpNum_perm_quicksort (perm + tot1, qpenalty + tot1, tbounded);

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
		printf ("complain: incorrect final basis size\n");
		rval = E_SIMPLEX_ERROR;
		ILL_CLEANUP;
	}

CLEANUP:
	if (rval)
		ILLbasis_free_basisinfo (lp);

	ILL_IFFREE (irow, int);
	ILL_IFFREE (rrow, int);

	EGlpNumFreeArray (v);
	ILL_IFFREE (unitcol, int);
	ILL_IFFREE (icol, int);
	ILL_IFFREE (rcol, int);

	EGlpNumFreeArray (dj);
	ILL_IFFREE (perm, int);
	ILL_IFFREE (porder, int);
	ILL_IFFREE (plen, int);

	EGlpNumFreeArray (qpenalty);
	EGlpNumClearVar (seldj);
	EGlpNumClearVar (selv);
	EGlpNumClearVar (c_dj);
	EGlpNumClearVar (cmax);
	ILLbasis_clear_vardata (&vd);
	EG_RETURN (rval);
}

static int set_basis_indices (
	lpinfo * lp,
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
			fprintf (stderr, "Error in basis creation\n");
			return E_SIMPLEX_ERROR;
		}
	}
	if (b != lp->nrows)
	{
		fprintf (stderr, "Error 2 in basis creation\n");
		return E_SIMPLEX_ERROR;
	}
	else if (nb != lp->nnbasic)
	{
		fprintf (stderr, "Error 3 in basis creation\n");
		return E_SIMPLEX_ERROR;
	}
	return 0;
}

int ILLbasis_get_initial (
	lpinfo * lp,
	int algorithm)
{
	int rval = 0;
	int *vstat = NULL;

	ILLbasis_free_basisinfo (lp);
	ILLbasis_init_basisinfo (lp);
	rval = ILLbasis_build_basisinfo (lp);
	CHECKRVALG (rval, CLEANUP);

	ILL_SAFE_MALLOC (vstat, lp->ncols, int);

	if (algorithm == PRIMAL_SIMPLEX)
		rval = get_initial_basis1 (lp, vstat);
	else
		rval = get_initial_basis2 (lp, vstat);

	if (rval == E_SIMPLEX_ERROR)
	{
		#ifdef HAVE_ZLIB_H
		EGioFile_t *f = EGioOpen ("bad.lp.gz", "w");
		#else
		#ifdef HAVE_BZLIB_H
		EGioFile_t *f = EGioOpen ("bad.lp.bz2", "w");
		#else
		EGioFile_t *f = EGioOpen ("bad.lp", "w");
		#endif
		#endif
		int tval = ILLwrite_lp_file (lp->O, f, NULL);
		if (tval)
		{
			fprintf (stderr, "Error writing bad lp\n");
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
	EGlpNum_t pinf1,
	EGlpNum_t dinf1,
	EGlpNum_t pinf2,
	EGlpNum_t dinf2)
{
/* We changed the constant definitions outside here, the actual numbers are
 * asigned in lpdata.c. the values are as follows:
 * CB_EPS = 0.001;
 * CB_PRI_RLIMIT = 0.25;
 * CB_INF_RATIO = 10.0; 
 * */
	int choice = 1;
	EGlpNum_t rp, rd;

	if (algorithm == PRIMAL_SIMPLEX)
	{
		EGlpNumInitVar (rp);
		EGlpNumInitVar (rd);
		EGlpNumCopyDiff (rp, pinf1, pinf2);
		EGlpNumCopyDiff (rd, dinf1, dinf2);
		if (EGlpNumIsLeq (rp, CB_EPS) && EGlpNumIsLeq (rd, CB_EPS))
			choice = 1;
		else
		{
			EGlpNumSign (rp);
			EGlpNumSign (rd);
			if (EGlpNumIsLeq (rp, CB_EPS) && EGlpNumIsLeq (rd, CB_EPS))
				choice = 2;
			else if (EGlpNumIsLess (pinf1, pinf2) && EGlpNumIsLess (dinf2, dinf1))
			{
				choice = 1;
				EGlpNumCopyFrac (rp, pinf1, pinf2);
				EGlpNumCopyFrac (rd, dinf2, dinf1);
				EGlpNumMultTo (rd, CB_INF_RATIO);
				if (EGlpNumIsLess (CB_PRI_RLIMIT, rp) && (EGlpNumIsLess (rd, rp)))
					choice = 2;
			}
			else if (EGlpNumIsLess (pinf2, pinf1) && EGlpNumIsLess (dinf1, dinf2))
			{
				choice = 2;
				EGlpNumCopyFrac (rp, pinf2, pinf1);
				EGlpNumCopyFrac (rd, dinf1, dinf2);
				EGlpNumMultTo (rd, CB_INF_RATIO);
				if (EGlpNumIsLess (CB_PRI_RLIMIT, rp) && EGlpNumIsLess (rd, rp))
					choice = 1;
			}
			else
				choice = 1;
		}
		EGlpNumClearVar (rp);
		EGlpNumClearVar (rd);
	}
	ILL_IFTRACE ("%s:%d\n", __func__, choice);
	return choice;
}

int ILLbasis_get_cinitial (
	lpinfo * lp,
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
	EGlpNum_t pinf1, pinf2, dinf1, dinf2;
	feas_info fi;

	EGlpNumInitVar (pinf1);
	EGlpNumInitVar (pinf2);
	EGlpNumInitVar (dinf1);
	EGlpNumInitVar (dinf2);
	EGlpNumInitVar (fi.totinfeas);

	ILLbasis_free_basisinfo (lp);
	ILLbasis_init_basisinfo (lp);
	rval = ILLbasis_build_basisinfo (lp);
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
	rval = ILLbasis_factor (lp, &singular);
	if (singular)
		MESSAGE (__QS_SB_VERB, "Singular Basis found!");
	CHECKRVALG (rval, CLEANUP);

	ILLfct_compute_piz (lp);
	ILLfct_compute_dz (lp);
	ILLfct_dual_adjust (lp, zeroLpNum);
	ILLfct_compute_xbz (lp);

	ILLfct_check_pfeasible (lp, &fi, lp->tol->pfeas_tol);
	ILLfct_check_dfeasible (lp, &fi, lp->tol->dfeas_tol);
	EGlpNumCopy (pinf1, lp->pinfeas);
	EGlpNumCopy (dinf1, lp->dinfeas);
	/*
	 * ILLfct_compute_pobj (lp);  obj1p = lp->objval;
	 * ILLfct_compute_dobj (lp);  obj1d = lp->objval;
	 */

	/* handle second basis */
	rval = set_basis_indices (lp, vstat2);
	CHECKRVALG (rval, CLEANUP);
#if BASIS_STATS > 0
	for (i = 0; i < lp->nrows; i++)
		nz2 += lp->matcnt[lp->baz[i]];
#endif
	rval = ILLbasis_factor (lp, &singular);
	if (singular)
		MESSAGE (__QS_SB_VERB, "Singular Basis found!");
	CHECKRVALG (rval, CLEANUP);

	ILLfct_compute_piz (lp);
	ILLfct_compute_dz (lp);
	ILLfct_dual_adjust (lp, zeroLpNum);
	ILLfct_compute_xbz (lp);

	ILLfct_check_pfeasible (lp, &fi, lp->tol->pfeas_tol);
	ILLfct_check_dfeasible (lp, &fi, lp->tol->dfeas_tol);
	EGlpNumCopy (pinf2, lp->pinfeas);
	EGlpNumCopy (dinf2, lp->dinfeas);

#if BASIS_STATS > 0
	printf ("b1: nz %d pinf %.2f dinf %.2f\n", nz1, EGlpNumToLf (pinf1),
					EGlpNumToLf (dinf1));
	printf ("b2: nz %d pinf %.2f dinf %.2f\n", nz2, EGlpNumToLf (pinf2),
					EGlpNumToLf (dinf2));
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
		#ifdef HAVE_ZLIB_H
		EGioFile_t *fil = EGioOpen ("bad.lp.gz", "w");
		#else
		#ifdef HAVE_BZLIB_H
		EGioFile_t *fil = EGioOpen ("bad.lp.bz2", "w");
		#else
		EGioFile_t *fil = EGioOpen ("bad.lp", "w");
		#endif
		#endif
		int tval = ILLwrite_lp_file (lp->O, fil, NULL);

		if (tval)
		{
			fprintf (stderr, "Error writing bad lp\n");
		}
		if (fil != NULL)
			EGioClose (fil);
	}
	ILL_IFFREE (vstat1, int);
	ILL_IFFREE (vstat2, int);

	EGlpNumClearVar (pinf1);
	EGlpNumClearVar (pinf2);
	EGlpNumClearVar (dinf1);
	EGlpNumClearVar (dinf2);
	EGlpNumClearVar (fi.totinfeas);
	EG_RETURN (rval);
}

int ILLbasis_factor (
	lpinfo * lp,
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
			ILLfactor_free_factor_work (lp->f);
		}
		else
		{
			ILL_SAFE_MALLOC (lp->f, 1, factor_work);
			EGlpNumInitVar (lp->f->fzero_tol);
			EGlpNumInitVar (lp->f->szero_tol);
			EGlpNumInitVar (lp->f->partial_tol);
			EGlpNumInitVar (lp->f->maxelem_orig);
			EGlpNumInitVar (lp->f->maxelem_factor);
			EGlpNumInitVar (lp->f->maxelem_cur);
			EGlpNumInitVar (lp->f->partial_cur);
			ILLfactor_init_factor_work (lp->f);
		}
		rval = ILLfactor_create_factor_work (lp->f, lp->O->nrows);
		CHECKRVALG (rval, CLEANUP);

		rval = ILLfactor (lp->f, lp->baz, lp->matbeg, lp->matcnt,
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

				ILLfct_update_basis_info (lp, eindex, lindex, lvstat);
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
		fprintf (stderr, "Error: unknown in %s\n", __func__);
	EG_RETURN (rval);
}

int ILLbasis_refactor (
	lpinfo * lp)
{
	int sing = 0;
	int rval = 0;

	rval = ILLbasis_factor (lp, &sing);
	if (sing)
	{
		MESSAGE (__QS_SB_VERB, "Singular Basis found!");
		rval = QS_LP_CHANGE_PREC;
		return rval;
	}
	EG_RETURN (rval);
}

void ILLbasis_column_solve (
	lpinfo * lp,
	svector * rhs,
	svector * soln)
{
	ILLfactor_ftran (lp->f, rhs, soln);
}

void ILLbasis_column_solve_update (
	lpinfo * lp,
	svector * rhs,
	svector * upd,
	svector * soln)
{
	ILLfactor_ftran_update (lp->f, rhs, upd, soln);
}

void ILLbasis_row_solve (
	lpinfo * lp,
	svector * rhs,
	svector * soln)
{
	ILLfactor_btran (lp->f, rhs, soln);
}

int ILLbasis_update (
	lpinfo * lp,
	svector * y,
	int lindex,
	int *refactor,
	int *singular)
{
#if 0														/* To always refactor, change 0 to 1 */
	*refactor = 1;
	return ILLbasis_factor (lp, singular);
#else

	int rval = 0;

	*refactor = 0;
	rval = ILLfactor_update (lp->f, y, lindex, refactor);
	if (rval == E_FACTOR_BLOWUP || rval == E_UPDATE_SINGULAR_ROW
			|| rval == E_UPDATE_SINGULAR_COL)
	{
/* Bico - comment out for dist
       fprintf(stderr, "Warning: numerically bad basis in ILLfactor_update\n");
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
		rval = ILLbasis_factor (lp, singular);
		if (*singular)
			MESSAGE (__QS_SB_VERB, "Singular Basis found!");
	}
	if (rval)
	{
		EGioFile_t *eout = 0;
		int tval;

		printf ("write bad lp to factor.lp\n");
		fflush (stdout);
		#ifdef HAVE_ZLIB_H
		eout = EGioOpen ("factor.lp.gz", "w");
		#else
		#ifdef HAVE_BZLIB_H
		eout = EGioOpen ("factor.lp.bz2", "w");
		#else
		eout = EGioOpen ("factor.lp", "w");
		#endif
		#endif
		if (!eout)
		{
			fprintf (stderr, "could not open file to write bad factor lp\n");
		}
		else
		{
			tval = ILLwrite_lp_file (lp->O, eout, NULL);
			if (tval)
			{
				fprintf (stderr, "error while writing bad factor lp\n");
			}
			EGioClose (eout);
		}

		printf ("write bad basis to factor.bas\n");
		fflush (stdout);
		tval = ILLlib_writebasis (lp, 0, "factor.bas");
		if (tval)
		{
			fprintf (stderr, "error while writing factor basis\n");
		}
	}

	EG_RETURN (rval);
#endif
}
