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

/* "$RCSfile: ratio.c,v $ $Revision: 1.2 $ $Date: 2003/11/05 16:49:52 $"; */
static int TRACE = 0;

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "eg_lpnum.h"
#include "eg_io.h"

#include "sortrus.h"
#include "stddefs.h"
#include "iqsutil.h"
#include "lpdefs.h"
#include "ratio.h"
#include "fct.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

void ILLratio_pI_test (
	lpinfo * lp,
	int eindex,
	int dir,
	ratio_res * rs)
{
	int i = 0, k = 0;
	int col, ecol;
	int cbnd, indx = 0;
	int tctr = 0;
	int *perm = lp->upd.perm;
	int *ix = lp->upd.ix;
	EGlpNum_t *pivtol = &(lp->tol->pivot_tol);
	EGlpNum_t *dftol = &(lp->tol->id_tol);

	 /*HHH*/ EGlpNum_t * t = lp->upd.t;
	EGlpNum_t t_i, delta, y_ij, rcost, nrcost, ntmp;
	EGlpNum_t *x, *l, *u;

	 /*HHH*/ EGlpNumInitVar (t_i);
	EGlpNumInitVar (delta);
	EGlpNumInitVar (y_ij);
	EGlpNumInitVar (rcost);
	EGlpNumInitVar (nrcost);
	EGlpNumInitVar (ntmp);
	EGlpNumZero (t_i);
	EGlpNumZero (y_ij);
	EGlpNumZero (delta);
	rs->lindex = -1;
	EGlpNumZero (rs->tz);
	EGlpNumZero (rs->pivotval);
	rs->ratio_stat = RATIO_FAILED;
	rs->lvstat = -1;
	ecol = lp->nbaz[eindex];
	ILL_IFTRACE2 ("%s:%d:%d:%d:%d", __func__, eindex, dir, ecol,
								(VBOUNDED == lp->vtype[ecol]));
	if (lp->vtype[ecol] == VBOUNDED)
	{
		EGlpNumCopyDiff (t[0], lp->uz[ecol], lp->lz[ecol]);
		ix[0] = BBOUND;
		ILL_IFTRACE2 (":%d[%d](%la,%la,%la)\n", ix[tctr], tctr,
									EGlpNumToLf (t[tctr]), EGlpNumToLf (lp->uz[ecol]),
									EGlpNumToLf (lp->lz[ecol]));
		tctr++;
	}
	ILL_IFTRACE2 (":%d", lp->yjz.nzcnt);
	for (k = 0; k < lp->yjz.nzcnt; k++)
	{
		EGlpNumCopy (y_ij, lp->yjz.coef[k]);
		if (!EGlpNumIsNeqZero (y_ij, *pivtol))
			continue;

		i = lp->yjz.indx[k];
		x = &(lp->xbz[i]);
		col = lp->baz[i];
		l = &(lp->lz[col]);
		u = &(lp->uz[col]);

		if ((dir == VINCREASE && EGlpNumIsGreatZero (y_ij)) ||
				(dir == VDECREASE && EGlpNumIsLessZero (y_ij)))
		{
			if (EGlpNumIsLessZero (y_ij))
				EGlpNumSign (y_ij);
			ILL_IFTRACE2 (":%d", lp->bfeas[i]);
			if (lp->bfeas[i] > 0)
			{
				EGlpNumCopyDiffRatio (t[tctr], *x, *u, y_ij);
				ix[tctr] = 10 * k + BATOUPPER;
				ILL_IFTRACE2 (":%d[%d](%la)\n", ix[tctr], tctr, EGlpNumToLf (t[tctr]));
				tctr++;
				if (EGlpNumIsNeqq (*l, NINFTY))
				{
					EGlpNumCopyDiffRatio (t[tctr], *x, *l, y_ij);
					ix[tctr] = 10 * k + BATOLOWER;
					ILL_IFTRACE2 (":%d[%d](%la)\n", ix[tctr], tctr,
												EGlpNumToLf (t[tctr]));
					tctr++;
				}
			}
			else if (lp->bfeas[i] == 0)
			{
				if (EGlpNumIsNeqq (*l, NINFTY))
				{
					EGlpNumCopyDiffRatio (t[tctr], *x, *l, y_ij);
					ix[tctr] = 10 * k + BATOLOWER;
					ILL_IFTRACE2 (":%d[%d](%la)\n", ix[tctr], tctr,
												EGlpNumToLf (t[tctr]));
					tctr++;
				}
			}
		}
		else if ((dir == VINCREASE && EGlpNumIsLessZero (y_ij)) ||
						 (dir == VDECREASE && EGlpNumIsGreatZero (y_ij)))
		{
			if (EGlpNumIsLessZero (y_ij))
				EGlpNumSign (y_ij);
			ILL_IFTRACE2 (":%d", lp->bfeas[i]);
			if (lp->bfeas[i] < 0)
			{
				EGlpNumCopyDiffRatio (t[tctr], *l, *x, y_ij);
				ix[tctr] = 10 * k + BBTOLOWER;
				ILL_IFTRACE2 (":%d[%d](%la)\n", ix[tctr], tctr, EGlpNumToLf (t[tctr]));
				tctr++;
				if (EGlpNumIsNeqq (*u, INFTY))
				{
					EGlpNumCopyDiffRatio (t[tctr], *u, *x, y_ij);
					ix[tctr] = 10 * k + BBTOUPPER;
					ILL_IFTRACE2 (":%d[%d](%la)\n", ix[tctr], tctr,
												EGlpNumToLf (t[tctr]));
					tctr++;
				}
			}
			else if (lp->bfeas[i] == 0)
			{
				if (EGlpNumIsNeqq (*u, INFTY))
				{
					EGlpNumCopyDiffRatio (t[tctr], *u, *x, y_ij);
					ix[tctr] = 10 * k + BBTOUPPER;
					ILL_IFTRACE2 (":%d[%d](%la)\n", ix[tctr], tctr,
												EGlpNumToLf (t[tctr]));
					tctr++;
				}
			}
		}
	}
	if (tctr == 0)
	{
		rs->ratio_stat = RATIO_FAILED;
		ILL_CLEANUP;
	}

	for (i = 0; i < tctr; i++)
		perm[i] = i;
	ILLutil_EGlpNum_perm_quicksort (perm, t, tctr);

	EGlpNumZero (lp->upd.c_obj);
	EGlpNumCopy (rcost, lp->pIdz[eindex]);
	ILL_IFTRACE2 ("\n%s:%d:%lf", __func__, tctr, EGlpNumToLf (rcost));
	for (i = 0; i < tctr; i++)
	{
		EGlpNumCopy (t_i, t[perm[i]]);
		EGlpNumCopy (ntmp, t_i);
		EGlpNumSubTo (ntmp, delta);
		EGlpNumAddInnProdTo (lp->upd.c_obj, ntmp, rcost);
		EGlpNumCopy (delta, t_i);
		ILL_IFTRACE2 (":%d:%lf", perm[i], EGlpNumToLf (delta));
		 /*HHH*/ cbnd = ix[perm[i]] % 10;
		if (cbnd != BBOUND)
		{
			k = ix[perm[i]] / 10;
			EGlpNumCopy (y_ij, lp->yjz.coef[k]);
			indx = lp->yjz.indx[k];
			ILL_IFTRACE2 (":%d", indx);
		}

		switch (cbnd)
		{
		case BBOUND:
			rs->ratio_stat = RATIO_NOBCHANGE;
			EGlpNumCopy (rs->tz, t_i);
			if (dir != VINCREASE)
				EGlpNumSign (rs->tz);
			ILL_CLEANUP;

		case BATOLOWER:
		case BATOUPPER:
			EGlpNumAddTo (rcost, y_ij);
			break;
		case BBTOLOWER:
		case BBTOUPPER:
			EGlpNumSubTo (rcost, y_ij);
			break;
		}
		EGlpNumCopyNeg (nrcost, rcost);
		if ((dir == VINCREASE && EGlpNumIsLeq (nrcost, *dftol)) ||
				(dir == VDECREASE && EGlpNumIsLeq (rcost, *dftol)))
		{
			/* change 5 to -1 if t_i > 0 is required below */
			if (EGlpNumIsLessZero (t_i) && i > 5)
			{
				/* printf ("pIhell %.5f %d\n", t_i, i); */
				EGlpNumDivUiTo (t_i, 2);
				rs->ratio_stat = RATIO_NEGATIVE;
				EGlpNumZero (rs->tz);
				ILL_CLEANUP;
			}
			rs->lindex = indx;
			rs->ratio_stat = RATIO_BCHANGE;
			if (cbnd == BATOLOWER || cbnd == BBTOLOWER)
				rs->lvstat = STAT_LOWER;
			else
				rs->lvstat = STAT_UPPER;

			EGlpNumCopy (rs->pivotval, y_ij);
			EGlpNumCopy (rs->tz, t_i);
			if (dir != VINCREASE)
				EGlpNumSign (rs->tz);
			ILL_CLEANUP;
		}
	}

CLEANUP:
	ILLfct_update_counts (lp, CNT_PIPIV, 0, rs->pivotval);
	ILL_IFTRACE2 (":tctr %d:%d\n", tctr, rs->ratio_stat);
	lp->upd.tctr = tctr;
	lp->upd.i = i;
	EGlpNumCopy (lp->upd.tz, t_i);
	EGlpNumCopy (lp->upd.piv, rs->pivotval);
	if (dir == VDECREASE)
		EGlpNumSign (lp->upd.c_obj);
	if (rs->lindex != -1)
		lp->upd.fs = lp->bfeas[rs->lindex];
	EGlpNumClearVar (t_i);
	EGlpNumClearVar (delta);
	EGlpNumClearVar (y_ij);
	EGlpNumClearVar (rcost);
	EGlpNumClearVar (nrcost);
	EGlpNumClearVar (ntmp);
}

void ILLratio_pII_test (
	lpinfo * lp,
	int eindex,
	int dir,
	ratio_res * rs)
{
	int i, k, indx, col, ecol;
	EGlpNum_t *x, *l, *u, t_max, ayi_max, yi_max, ay_ij, y_ij, t_i, t_z;
	EGlpNum_t *pivtol = &(lp->tol->pivot_tol);
	EGlpNum_t *pftol = &(lp->tol->pfeas_tol);

	EGlpNumInitVar (y_ij);
	EGlpNumInitVar (ay_ij);
	EGlpNumInitVar (t_i);
	EGlpNumInitVar (t_z);
	EGlpNumInitVar (t_max);
	EGlpNumInitVar (yi_max);
	EGlpNumInitVar (ayi_max);
	 /*HHH*/ rs->boundch = 0;
	rs->lindex = -1;
	EGlpNumZero (rs->tz);
	rs->ratio_stat = RATIO_FAILED;
	rs->lvstat = -1;
	EGlpNumZero (rs->pivotval);
	EGlpNumZero (rs->lbound);
	ecol = lp->nbaz[eindex];

	for (k = 0, EGlpNumCopy (t_max, INFTY); k < lp->yjz.nzcnt; k++)
	{
		EGlpNumCopy (y_ij, lp->yjz.coef[k]);
		EGlpNumCopyAbs (ay_ij, y_ij);
		if (!EGlpNumIsNeqZero (y_ij, *pivtol))
			continue;

		EGlpNumCopy (t_i, INFTY);
		i = lp->yjz.indx[k];
		x = &(lp->xbz[i]);
		col = lp->baz[i];
		l = &(lp->lz[col]);
		u = &(lp->uz[col]);

		if ((dir == VINCREASE && EGlpNumIsGreatZero (y_ij)) ||
				(dir == VDECREASE && EGlpNumIsLessZero (y_ij)))
		{
			if (EGlpNumIsNeqq (*l, NINFTY))
			{
				EGlpNumCopyDiff (t_i, *x, *l);
				EGlpNumAddTo (t_i, *pftol);
				EGlpNumDivTo (t_i, ay_ij);
			}
		}
		else if ((dir == VINCREASE && EGlpNumIsLessZero (y_ij)) ||
						 (dir == VDECREASE && EGlpNumIsGreatZero (y_ij)))
		{
			if (EGlpNumIsNeqq (*u, INFTY))
			{
				EGlpNumCopySum (t_i, *u, *pftol);
				EGlpNumSubTo (t_i, *x);
				EGlpNumDivTo (t_i, ay_ij);
			}
		}
		if (EGlpNumIsEqqual (t_i, INFTY))
			continue;

		if (EGlpNumIsLess (t_i, t_max))
		{
			/*HHH tind = i; yval = fabs (y_ij); tval = t_i - pftol/fabs(y_ij); */
			EGlpNumCopy (t_max, t_i);
		}
	}
	/* we use yi_max as temporal variable here */
	EGlpNumCopyDiff (yi_max, lp->uz[ecol], lp->lz[ecol]);
	if (lp->vtype[ecol] == VBOUNDED && EGlpNumIsLeq (yi_max, t_max))
	{

		EGlpNumCopy (t_max, yi_max);
		rs->ratio_stat = RATIO_NOBCHANGE;
		EGlpNumCopy (rs->tz, t_max);
		if (dir != VINCREASE)
			EGlpNumSign (rs->tz);
		ILL_CLEANUP;
	}

	if (EGlpNumIsLeq (INFTY, t_max))
	{
		rs->ratio_stat = RATIO_UNBOUNDED;
		ILL_CLEANUP;
	}
	/*if (EGlpNumIsLess (t_max, zeroLpNum))
	 * printf ("pIIhell\n");
	 */
	indx = -1;
	EGlpNumZero (t_z);
	EGlpNumZero (yi_max);
	EGlpNumZero (ayi_max);
	ILL_IFTRACE2 (":%d", lp->yjz.nzcnt);
	for (k = 0; k < lp->yjz.nzcnt; k++)
	{
		EGlpNumCopy (y_ij, lp->yjz.coef[k]);
		EGlpNumCopyAbs (ay_ij, y_ij);
		if (!EGlpNumIsNeqZero (y_ij, *pivtol))
			continue;

		EGlpNumCopy (t_i, INFTY);
		i = lp->yjz.indx[k];
		x = &(lp->xbz[i]);
		col = lp->baz[i];
		l = &(lp->lz[col]);
		u = &(lp->uz[col]);

		if ((dir == VINCREASE && EGlpNumIsGreatZero (y_ij)) ||
				(dir == VDECREASE && EGlpNumIsLessZero (y_ij)))
		{
			if (EGlpNumIsNeqq (*l, NINFTY))
				EGlpNumCopyDiffRatio (t_i, *x, *l, ay_ij);
		}
		else if ((dir == VINCREASE && EGlpNumIsLessZero (y_ij)) ||
						 (dir == VDECREASE && EGlpNumIsGreatZero (y_ij)))
		{
			if (EGlpNumIsNeqq (*u, INFTY))
				EGlpNumCopyDiffRatio (t_i, *u, *x, ay_ij);
		}

		if (EGlpNumIsLeq (t_i, t_max))
		{
			if (EGlpNumIsLess (ayi_max, ay_ij))
			{
				EGlpNumCopy (yi_max, y_ij);
				EGlpNumCopy (ayi_max, ay_ij);
				indx = i;
				EGlpNumCopy (t_z, t_i);
				ILL_IFTRACE2 (":%d:%lf:%lf:%lf:%lf", indx, EGlpNumToLf (t_i),
											EGlpNumToLf (t_max), EGlpNumToLf (ayi_max),
											EGlpNumToLf (ay_ij));
			}
		}
	}

	if (indx < 0)
	{
		rs->ratio_stat = RATIO_FAILED;
	}
	else
	{
		/*
		 * if (tind != rs->lindex){
		 * HHHprintf ("tmax %e tval = %e yval = %e tind = %d\n", t_max, tval, yval, tind);
		 * HHHprintf ("h tval = %e yval = %e tind = %d\n",rs->tz, yi_max, rs->lindex);
		 * }
		 */
		ILL_IFTRACE2 (":%d", indx);
		rs->lindex = indx;
		EGlpNumCopy (rs->tz, t_z);
		EGlpNumCopy (rs->pivotval, yi_max);
		rs->ratio_stat = RATIO_BCHANGE;

		if (dir == VINCREASE)
			rs->lvstat =
				(EGlpNumIsGreatZero (yi_max)) ? STAT_LOWER : STAT_UPPER;
		else
			rs->lvstat =
				(EGlpNumIsGreatZero (yi_max)) ? STAT_UPPER : STAT_LOWER;

		if (EGlpNumIsLessZero (rs->tz))
		{
			ILL_IFTRACE2 ("need to change bound, tz=%la\n", EGlpNumToLf (rs->tz));
			EGlpNumCopyAbs (rs->tz, t_max);
			EGlpNumDivUiTo (rs->tz, 10);
			rs->boundch = 1;
			EGlpNumCopy (rs->lbound, lp->xbz[rs->lindex]);
			if (rs->lvstat == STAT_LOWER)
				EGlpNumSubInnProdTo (rs->lbound, rs->tz, ayi_max);
			else
				EGlpNumAddInnProdTo (rs->lbound, rs->tz, ayi_max);
		}
		if (dir == VDECREASE)
			EGlpNumSign (rs->tz);
	}
CLEANUP:
	ILLfct_update_counts (lp, CNT_PIIPIV, 0, rs->pivotval);
	EGlpNumClearVar (y_ij);
	EGlpNumClearVar (ay_ij);
	EGlpNumClearVar (t_i);
	EGlpNumClearVar (t_z);
	EGlpNumClearVar (t_max);
	EGlpNumClearVar (yi_max);
	EGlpNumClearVar (ayi_max);
}

#define GET_XY_DRATIOTEST \
      if (lp->vstat[col] == STAT_UPPER){ \
				EGlpNumCopyNeg(x,lp->dz[j]);\
        EGlpNumCopy(y, *zAj);\
      } \
      else{ \
         EGlpNumCopy(x, lp->dz[j]); \
         EGlpNumCopyNeg(y, *zAj);\
      } \
      if (lvstat == STAT_UPPER) \
         EGlpNumSign(y);


void ILLratio_dI_test (
	lpinfo * lp,
	int lindex,
	int lvstat,
	ratio_res * rs)
{
	int j = 0, k;
	int col;
	int cbnd, indx;
	int tctr = 0;
	int *perm = lp->upd.perm;
	int *ix = lp->upd.ix;
	EGlpNum_t *t = lp->upd.t;
	EGlpNum_t *zAj, x, y, t_j, theta, rcost, delta;
	EGlpNum_t *pftol = &(lp->tol->ip_tol);
	EGlpNum_t *pivtol = &(lp->tol->pivot_tol);

	EGlpNumInitVar (x);
	EGlpNumInitVar (y);
	EGlpNumInitVar (t_j);
	EGlpNumInitVar (theta);
	EGlpNumInitVar (rcost);
	EGlpNumInitVar (delta);
	EGlpNumZero (delta);
	EGlpNumZero (t_j);
	EGlpNumZero (rs->tz);
	 /*HHH*/ rs->eindex = -1;
	rs->ratio_stat = RATIO_FAILED;
	EGlpNumZero (rs->pivotval);

	for (k = 0; k < lp->zA.nzcnt; k++)
	{
		zAj = &(lp->zA.coef[k]);
		if (!EGlpNumIsNeqZero (*zAj, *pivtol))
			continue;

		EGlpNumCopy (t_j, INFTY);
		j = lp->zA.indx[k];
		col = lp->nbaz[j];

		if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED)
			continue;

		GET_XY_DRATIOTEST;

		if (EGlpNumIsLessZero (y))
		{
			if (lp->dfeas[j] != 0 && lp->vstat[col] != STAT_ZERO)
			{
				EGlpNumCopyFrac (t[tctr], x, y);
				ix[tctr] = 10 * k + BBTOLOWER;
				tctr++;
			}
			else if (lp->vstat[col] == STAT_ZERO)
			{
				if (lp->dfeas[j] < 0)
				{
					EGlpNumCopyFrac (t[tctr], x, y);
					ix[tctr] = 10 * k + BBTOLOWER;
					tctr++;
				}
				if (lp->dfeas[j] <= 0)
				{
					EGlpNumCopyFrac (t[tctr], x, y);
					ix[tctr] = 10 * k + BBTOUPPER;
					tctr++;
				}
			}
		}
		else
		{
			if (lp->dfeas[j] > 0)
			{
				if (lp->vstat[col] == STAT_ZERO)
				{
					EGlpNumCopyFrac (t[tctr], x, y);
					ix[tctr] = 10 * k + BATOUPPER;
					tctr++;
					EGlpNumCopyFrac (t[tctr], x, y);
					ix[tctr] = 10 * k + BATOLOWER;
					tctr++;
				}
			}
			else if (lp->dfeas[j] == 0)
			{
				EGlpNumCopyFrac (t[tctr], x, y);
				if (lp->vtype[col] == VBOUNDED)
					ix[tctr] = 10 * k + BSKIP;
				else
					ix[tctr] = 10 * k + BATOLOWER;
				tctr++;
			}
		}
	}

	if (tctr == 0)
	{
		rs->ratio_stat = RATIO_FAILED;
		ILL_CLEANUP;
	}

	for (j = 0; j < tctr; j++)
		perm[j] = j;
	ILLutil_EGlpNum_perm_quicksort (perm, t, tctr);

	EGlpNumZero (lp->upd.c_obj);
	EGlpNumCopy (rcost, lp->xbz[lindex]);
	if (lvstat == STAT_LOWER)
		EGlpNumSign (rcost);
	for (j = 0; j < tctr; j++)
	{
		cbnd = ix[perm[j]] % 10;
		if (cbnd == BSKIP)
			continue;

		EGlpNumCopy (t_j, t[perm[j]]);
		EGlpNumCopy (x, t_j);
		EGlpNumSubTo (x, delta);
		EGlpNumAddInnProdTo (lp->upd.c_obj, x, rcost);
		EGlpNumCopy (delta, t_j);
		k = ix[perm[j]] / 10;
		zAj = &(lp->zA.coef[k]);
		indx = lp->zA.indx[k];

		if (lp->vstat[lp->nbaz[indx]] == STAT_LOWER
				|| lp->vstat[lp->nbaz[indx]] == STAT_ZERO)
			EGlpNumCopyNeg (theta, *zAj);
		else
			EGlpNumCopy (theta, *zAj);

		if (lvstat == STAT_UPPER)
			EGlpNumSign (theta);

		switch (cbnd)
		{
		case BATOLOWER:
		case BATOUPPER:
			EGlpNumSubTo (rcost, theta);
			break;
		case BBTOLOWER:
		case BBTOUPPER:
			EGlpNumAddTo (rcost, theta);
			break;
		}
		if (EGlpNumIsLeq (rcost, *pftol))
		{
			/* if (t_j < 0.0) printf ("dIhell\n"); */
			rs->eindex = indx;
			EGlpNumCopy (rs->tz, t_j);
			EGlpNumCopy (rs->pivotval, *zAj);
			rs->ratio_stat = RATIO_BCHANGE;
			ILL_CLEANUP;
		}
	}

CLEANUP:
	ILLfct_update_counts (lp, CNT_DIPIV, 0, rs->pivotval);
	ILL_IFTRACE2 ("%s:tctr %d\n", __func__, tctr);
	lp->upd.tctr = tctr;
	lp->upd.i = j;
	EGlpNumCopyAbs (lp->upd.tz, t_j);
	EGlpNumCopy (lp->upd.piv, rs->pivotval);
	if (rs->eindex != -1)
		lp->upd.fs = lp->dfeas[rs->eindex];
	EGlpNumClearVar (x);
	EGlpNumClearVar (y);
	EGlpNumClearVar (t_j);
	EGlpNumClearVar (theta);
	EGlpNumClearVar (rcost);
	EGlpNumClearVar (delta);
}

void ILLratio_dII_test (
	lpinfo * lp,
	/*int lindex,*/
	int lvstat,
	ratio_res * rs)
{
	int j, k, indx;
	int col, ecol;
	EGlpNum_t *zAj, azAj, az_max, x, y, t_j, z_max, t_max, t_z;
	EGlpNum_t *dftol = &(lp->tol->dfeas_tol);
	EGlpNum_t *pivtol = &(lp->tol->pivot_tol);

	EGlpNumInitVar (x);
	EGlpNumInitVar (y);
	EGlpNumInitVar (t_j);
	EGlpNumInitVar (z_max);
	EGlpNumInitVar (t_max);
	EGlpNumInitVar (az_max);
	EGlpNumInitVar (azAj);
	EGlpNumInitVar (t_z);
	EGlpNumZero (t_j);
	rs->coeffch = 0;
	EGlpNumZero (rs->ecoeff);
	rs->eindex = -1;
	rs->ratio_stat = RATIO_FAILED;
	ILL_IFTRACE2 ("%s:tctr %d\n", __func__, 0);
	lp->upd.tctr = 0;
	EGlpNumZero (lp->upd.dty);
	for (k = 0, EGlpNumCopy (t_max, INFTY); k < lp->zA.nzcnt; k++)
	{
		zAj = &(lp->zA.coef[k]);
		if (!EGlpNumIsNeqZero (*zAj, *pivtol))
			continue;

		EGlpNumCopy (t_j, INFTY);
		j = lp->zA.indx[k];
		col = lp->nbaz[j];

		if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED)
			continue;

		GET_XY_DRATIOTEST;

//#warning adding/substracting tolerances to used value, is it rigght?
		if (EGlpNumIsGreatZero (y))
		{
			//t_j = (x + dftol) / y;
			EGlpNumCopySum (t_j, x, *dftol);
			EGlpNumDivTo (t_j, y);
		}
		else
		{
//#warning adding/substracting tolerances to used value, is it rigght?
			if (lp->vstat[col] == STAT_ZERO)
				EGlpNumCopyDiffRatio (t_j, x, *dftol, y);
		}
		//if (t_j == INFTY)
		if (EGlpNumIsEqqual (t_j, INFTY))
			continue;

		if (EGlpNumIsLess (t_j, t_max))
			EGlpNumCopy (t_max, t_j);
	}

	if (EGlpNumIsLeq (INFTY, t_max))
	{
		rs->ratio_stat = RATIO_UNBOUNDED;
		ILL_CLEANUP;
	}
	/* if (t_max < 0.0) printf ("dIIhell\n"); */

	indx = -1;
	EGlpNumZero (t_z);
	EGlpNumZero (z_max);
	EGlpNumZero (az_max);

	for (k = 0; k < lp->zA.nzcnt; k++)
	{
		zAj = &(lp->zA.coef[k]);
		EGlpNumCopyAbs (azAj, *zAj);
		if (!EGlpNumIsNeqZero (*zAj, *pivtol))
			continue;

		EGlpNumCopy (t_j, INFTY);
		j = lp->zA.indx[k];
		col = lp->nbaz[j];

		if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED)
			continue;

		GET_XY_DRATIOTEST;

		if (EGlpNumIsGreatZero (y) || lp->vstat[col] == STAT_ZERO)
			EGlpNumCopyFrac (t_j, x, y);

		if (EGlpNumIsLeq (t_j, t_max) && (EGlpNumIsLess (az_max, azAj)))
		{
			EGlpNumCopy (z_max, *zAj);
			EGlpNumCopy (az_max, azAj);
			indx = j;
			EGlpNumCopy (t_z, t_j);
		}
	}


	if (indx < 0)
	{
		rs->ratio_stat = RATIO_FAILED;
	}
	else
	{
		rs->eindex = indx;
		EGlpNumCopy (rs->tz, t_z);
		EGlpNumCopy (rs->pivotval, z_max);
		rs->ratio_stat = RATIO_BCHANGE;

		if (EGlpNumIsLessZero (rs->tz))
		{
			EGlpNumCopyAbs (rs->tz, t_max);
			EGlpNumDivUiTo (rs->tz, 20);
			rs->coeffch = 1;
			ecol = lp->nbaz[indx];
			EGlpNumCopyDiff (rs->ecoeff, lp->cz[ecol], lp->dz[indx]);
			switch (lp->vstat[ecol])
			{
			case STAT_LOWER:
				EGlpNumAddInnProdTo (rs->ecoeff, rs->tz, az_max);
				break;
			case STAT_UPPER:
				EGlpNumSubInnProdTo (rs->ecoeff, rs->tz, az_max);
				break;
			default:
				EGlpNumZero (rs->tz);
				break;
			}
		}
	}

CLEANUP:
	ILLfct_update_counts (lp, CNT_DIIPIV, 0, rs->pivotval);
	EGlpNumCopy (lp->upd.piv, rs->pivotval);
	EGlpNumClearVar (x);
	EGlpNumClearVar (y);
	EGlpNumClearVar (t_j);
	EGlpNumClearVar (z_max);
	EGlpNumClearVar (t_max);
	EGlpNumClearVar (t_z);
	EGlpNumClearVar (az_max);
	EGlpNumClearVar (azAj);
}

void ILLratio_longdII_test (
	lpinfo * lp,
	int lindex,
	int lvstat,
	ratio_res * rs)
{
	int j, k, indx = 0, tctr = 0;
	int col, ecol;
	int vs, bnd_exist = 0;
	int *perm = lp->upd.perm;
	int *ix = lp->upd.ix;
	int b_indx = -1;
	EGlpNum_t *t = lp->upd.t;
	EGlpNum_t *l,
		*u,
		*xb,
		*zAj = 0,
		x,
		y,
		t_j,
		z_max,
		t_max, t_z, theta, rcost, delta, zb_val, tb_val, az_max, azb_val, azAj;
	EGlpNum_t *pftol = &(lp->tol->pfeas_tol);
	EGlpNum_t *dftol = &(lp->tol->dfeas_tol);
	EGlpNum_t *pivtol = &(lp->tol->pivot_tol);

	EGlpNumInitVar (x);
	EGlpNumInitVar (azAj);
	EGlpNumInitVar (y);
	EGlpNumInitVar (t_j);
	EGlpNumInitVar (z_max);
	EGlpNumInitVar (az_max);
	EGlpNumInitVar (t_max);
	EGlpNumInitVar (t_z);
	EGlpNumInitVar (theta);
	EGlpNumInitVar (rcost);
	EGlpNumInitVar (delta);
	EGlpNumInitVar (zb_val);
	EGlpNumInitVar (azb_val);
	EGlpNumInitVar (tb_val);
	EGlpNumZero (t_j);
	EGlpNumZero (delta);
	EGlpNumZero (zb_val);
	EGlpNumZero (azb_val);
	EGlpNumCopy (tb_val, NINFTY);
//#warning not sure about THIS line
	EGlpNumZero (rs->pivotval);

	rs->coeffch = 0;
	rs->eindex = -1;
	rs->ratio_stat = RATIO_FAILED;

	ILL_IFTRACE2 ("%s:tctr %d\n", __func__, 0);
	lp->upd.tctr = 0;
	lp->upd.i = 0;
	EGlpNumZero (lp->upd.tz);
	EGlpNumZero (lp->upd.piv);
	EGlpNumZero (lp->upd.c_obj);
	EGlpNumZero (lp->upd.dty);

	xb = &(lp->xbz[lindex]);
	col = lp->baz[lindex];
	l = &(lp->lz[col]);
	u = &(lp->uz[col]);
	//rcost = (lvstat == STAT_LOWER) ? l - xb : xb - u;
	if (lvstat == STAT_LOWER)
		EGlpNumCopyDiff (rcost, *l, *xb);
	else
		EGlpNumCopyDiff (rcost, *xb, *u);

	for (k = 0, EGlpNumCopy (t_max, INFTY); k < lp->zA.nzcnt; k++)
	{
		zAj = &(lp->zA.coef[k]);
		if (!EGlpNumIsNeqZero (*zAj, *pivtol))
			continue;

		EGlpNumCopy (t_j, INFTY);
		j = lp->zA.indx[k];
		col = lp->nbaz[j];

		if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED)
			continue;
		if (lp->vtype[col] == VBOUNDED)
		{
			bnd_exist++;
			continue;
		}

		GET_XY_DRATIOTEST;

		if (EGlpNumIsGreatZero (y))
		{
			//t_j = (x + dftol) / y;
//#warning Using tolerances to add to result, is it right?
			EGlpNumCopySum (t_j, x, *dftol);
			EGlpNumDivTo (t_j, y);
		}
		else
		{
			if (lp->vstat[col] == STAT_ZERO)
				EGlpNumCopyDiffRatio (t_j, x, *dftol, y);
		}
		if (EGlpNumIsEqqual (t_j, INFTY))
			continue;

		if (EGlpNumIsLess (t_j, t_max))
			EGlpNumCopy (t_max, t_j);
	}
	if (EGlpNumIsLessZero (t_max))
	{
		/*printf ("dIIhell, %.4f\n", t_max); */
		rs->ratio_stat = RATIO_NEGATIVE;
		ILL_CLEANUP;
	}

	if (bnd_exist == 0 && EGlpNumIsLeq (INFTY, t_max))
	{
		rs->ratio_stat = RATIO_UNBOUNDED;
		/*
		 * printf ("x = %.8f, b = %.2f \n", lp->xbz[lindex], (lvstat == STAT_LOWER ) ? lp->lz[lp->baz[lindex]] : lp->uz[lp->baz[lindex]]);
		 */
		ILL_CLEANUP;
	}

	if (bnd_exist != 0)
	{
		for (k = 0; k < lp->zA.nzcnt; k++)
		{
			zAj = &(lp->zA.coef[k]);
			if (!EGlpNumIsNeqZero (*zAj, *pivtol))
				continue;

			EGlpNumCopy (t_j, INFTY);
			j = lp->zA.indx[k];
			col = lp->nbaz[j];

			if (lp->vtype[col] != VBOUNDED)
				continue;

			GET_XY_DRATIOTEST;

			if (EGlpNumIsGreatZero (y))
			{
				EGlpNumCopyFrac (t_j, x, y);
				if (EGlpNumIsLeq (t_j, t_max))
				{
					EGlpNumCopy (t[tctr], t_j);
					ix[tctr] = k;
					tctr++;
				}
			}
		}
	}

	if (tctr != 0)
	{
		for (j = 0; j < tctr; j++)
			perm[j] = j;
		ILLutil_EGlpNum_perm_quicksort (perm, t, tctr);

		for (j = 0; j < tctr; j++)
		{

			EGlpNumCopy (t_j, t[perm[j]]);
			/* we use x as temporal storage */
			//lp->upd.c_obj += (t_j - delta) * rcost;
			EGlpNumCopy (x, t_j);
			EGlpNumSubTo (x, delta);
			EGlpNumAddInnProdTo (lp->upd.c_obj, x, rcost);
			EGlpNumCopy (delta, t_j);
			 /*HHH*/ k = ix[perm[j]];
			zAj = &(lp->zA.coef[k]);
			indx = lp->zA.indx[k];
			col = lp->nbaz[indx];
			l = &(lp->lz[col]);
			u = &(lp->uz[col]);
			vs = lp->vstat[col];
			//theta = (vs == STAT_UPPER) ? (l - u) * zAj : (u - l) * zAj;
			EGlpNumCopyDiff (theta, *l, *u);
			EGlpNumMultTo (theta, *zAj);
			if (vs != STAT_UPPER)
				EGlpNumSign (theta);
			if (lvstat == STAT_LOWER)
				EGlpNumAddTo (rcost, theta);
			else
				EGlpNumSubTo (rcost, theta);

			if (EGlpNumIsLeq (rcost, *pftol))
			{
				rs->eindex = indx;
				EGlpNumCopy (rs->tz, t_j);
				EGlpNumCopy (rs->pivotval, *zAj);
				rs->ratio_stat = RATIO_BCHANGE;

				if (EGlpNumIsLessZero (rs->tz))
				{
					EGlpNumZero (rs->tz);
					rs->coeffch = 1;
					//rs->ecoeff = lp->cz[col] - lp->dz[indx];
					EGlpNumCopyDiff (rs->ecoeff, lp->cz[col], lp->dz[indx]);
					//lp->upd.c_obj += (rs->tz - delta) * rcost; note ts->tz == 0;
					EGlpNumSubInnProdTo (lp->upd.c_obj, delta, rcost);
				}
				ILL_IFTRACE2 ("%s:tctr %d\n", __func__, tctr);
				lp->upd.tctr = tctr;
				lp->upd.i = j;
				EGlpNumCopy (lp->upd.tz, rs->tz);
				ILL_CLEANUP;
			}
		}
		ILL_IFTRACE2 ("%s:tctr %d\n", __func__, tctr);
		lp->upd.tctr = tctr;
		lp->upd.i = tctr;
		EGlpNumCopy (lp->upd.tz, t_j);
		EGlpNumCopy (zb_val, *zAj);
		EGlpNumCopyAbs (azb_val, zb_val);
		EGlpNumCopy (tb_val, t_j);
		b_indx = indx;
	}

	if (bnd_exist != 0 && EGlpNumIsLeq (INFTY, t_max))
	{
		rs->ratio_stat = RATIO_UNBOUNDED;
		/* printf ("rcost: %.8f\n", rcost); */
		ILL_CLEANUP;
	}

	EGlpNumZero (z_max);
	EGlpNumZero (az_max);
	indx = -1;
	EGlpNumZero (t_z);
	for (k = 0; k < lp->zA.nzcnt; k++)
	{
		zAj = &(lp->zA.coef[k]);
		EGlpNumCopyAbs (azAj, *zAj);
		if (!EGlpNumIsNeqZero (*zAj, *pivtol))
			continue;

		EGlpNumCopy (t_j, INFTY);
		j = lp->zA.indx[k];
		col = lp->nbaz[j];

		if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED ||
				lp->vtype[col] == VBOUNDED)
			continue;

		GET_XY_DRATIOTEST;

		if (EGlpNumIsGreatZero (y) || lp->vstat[col] == STAT_ZERO)
			EGlpNumCopyFrac (t_j, x, y);

		if (EGlpNumIsLeq (t_j, t_max))
		{
			if (EGlpNumIsLess (az_max, azAj))
			{
				EGlpNumCopy (z_max, *zAj);
				EGlpNumCopy (az_max, azAj);
				indx = j;
				EGlpNumCopy (t_z, t_j);
			}
		}
	}

	if (indx < 0)
	{
		rs->ratio_stat = RATIO_FAILED;
		ILL_CLEANUP;
	}
	if ((tctr == 0) || (EGlpNumIsLessZero (tb_val)) ||
			(tctr != 0 && EGlpNumIsLeq (tb_val, t_z) &&
			 EGlpNumIsLeq (azb_val, az_max)))
	{
		/* we use x as temporal vvariable */
		/* lp->upd.c_obj += (t_z - delta) * rcost; */
		EGlpNumCopyDiff (x, t_z, delta);
		EGlpNumAddInnProdTo (lp->upd.c_obj, x, rcost);
		EGlpNumCopy (delta, t_z);
		rs->eindex = indx;
		EGlpNumCopy (rs->tz, t_z);
		EGlpNumCopy (rs->pivotval, z_max);
		rs->ratio_stat = RATIO_BCHANGE;
	}
	/* For now */
	else if (tctr != 0)
	{
		rs->eindex = b_indx;
		EGlpNumCopy (rs->tz, tb_val);
		EGlpNumCopy (rs->pivotval, zb_val);
		rs->ratio_stat = RATIO_BCHANGE;
		lp->upd.i -= 1;
	}

	if (EGlpNumIsLessZero (rs->tz))
	{
		/* if (tctr != 0) printf ("despite long step\n"); */
		/* rs->tz = fabs (t_max / 20.0); */
		EGlpNumCopyAbs (rs->tz, t_max);
		EGlpNumDivUiTo (rs->tz, 20);
		rs->coeffch = 1;

		ecol = lp->nbaz[indx];
		if (lp->vstat[ecol] == STAT_LOWER)
		{
			/*rs->ecoeff = lp->cz[ecol] - lp->dz[indx] + rs->tz * fabs (z_max); */
			EGlpNumCopy (rs->ecoeff, az_max);
			EGlpNumMultTo (rs->ecoeff, rs->tz);
			EGlpNumAddTo (rs->ecoeff, lp->cz[ecol]);
			EGlpNumSubTo (rs->ecoeff, lp->dz[indx]);
		}
		else if (lp->vstat[ecol] == STAT_UPPER)
		{
			/*rs->ecoeff = lp->cz[ecol] - lp->dz[indx] - rs->tz * fabs (z_max); */
			EGlpNumCopy (rs->ecoeff, az_max);
			EGlpNumMultTo (rs->ecoeff, rs->tz);
			EGlpNumSign (rs->ecoeff);
			EGlpNumAddTo (rs->ecoeff, lp->cz[ecol]);
			EGlpNumSubTo (rs->ecoeff, lp->dz[indx]);
		}
		else
		{
			/*rs->ecoeff = lp->cz[ecol] - lp->dz[indx]; */
			EGlpNumCopyDiff (rs->ecoeff, lp->cz[ecol], lp->dz[indx]);
			EGlpNumZero (rs->tz);
		}
		/* we use x as temporal storage */
		/*lp->upd.c_obj += (rs->tz - delta) * rcost; */
		EGlpNumCopy (x, rs->tz);
		EGlpNumSubTo (x, delta);
		EGlpNumAddInnProdTo (lp->upd.c_obj, x, rcost);
	}

CLEANUP:
	ILLfct_update_counts (lp, CNT_DIIPIV, 0, rs->pivotval);
	EGlpNumCopy (lp->upd.piv, rs->pivotval);
	EGlpNumClearVar (x);
	EGlpNumClearVar (y);
	EGlpNumClearVar (t_j);
	EGlpNumClearVar (z_max);
	EGlpNumClearVar (az_max);
	EGlpNumClearVar (t_max);
	EGlpNumClearVar (t_z);
	EGlpNumClearVar (theta);
	EGlpNumClearVar (rcost);
	EGlpNumClearVar (delta);
	EGlpNumClearVar (zb_val);
	EGlpNumClearVar (azb_val);
	EGlpNumClearVar (tb_val);
	EGlpNumClearVar (azAj);
}

void ILLratio_pivotin_test (
	lpinfo * lp,
	int *rlist,
	int rcnt,
	ratio_res * rs)
{
	int i, k, col;
	EGlpNum_t *x, *l, *u;
	EGlpNum_t ay_ij,
		at_i, at_l, at_u, ayi_max, y_ij, t_i, t_l, t_u, t_max, yi_max;
	EGlpNum_t *pivtol = &(lp->tol->pivot_tol);

	if (rcnt <= 0 || rs == NULL)
		return;
	EGlpNumInitVar (ay_ij);
	EGlpNumInitVar (at_i);
	EGlpNumInitVar (at_l);
	EGlpNumInitVar (at_u);
	EGlpNumInitVar (ayi_max);
	EGlpNumInitVar (t_max);
	EGlpNumInitVar (y_ij);
	EGlpNumInitVar (t_i);
	EGlpNumInitVar (t_l);
	EGlpNumInitVar (t_u);
	EGlpNumInitVar (yi_max);
	rs->boundch = 0;
	rs->lindex = -1;
	EGlpNumZero (rs->tz);
	rs->ratio_stat = RATIO_FAILED;
	rs->lvstat = -1;
	EGlpNumZero (rs->pivotval);
	EGlpNumZero (rs->lbound);

	for (i = 0; i < rcnt; i++)
		lp->iwork[rlist[i]] = 1;

	for (k = 0, EGlpNumCopy (t_max, INFTY); k < lp->yjz.nzcnt; k++)
	{
		EGlpNumCopy (y_ij, lp->yjz.coef[k]);
		if (!EGlpNumIsNeqZero (y_ij, *pivtol))
			continue;

		i = lp->yjz.indx[k];
		if (lp->iwork[lp->baz[i]] == 1)
			continue;
		x = &(lp->xbz[i]);
		col = lp->baz[i];
		l = &(lp->lz[col]);
		u = &(lp->uz[col]);
		EGlpNumCopy (t_u, INFTY);
		EGlpNumCopy (at_u, INFTY);
		EGlpNumCopy (t_l, NINFTY);
		EGlpNumCopy (at_l, INFTY);

		if (EGlpNumIsNeqq (*l, NINFTY))
		{
			EGlpNumCopyDiffRatio (t_l, *x, *l, y_ij);
			EGlpNumCopyAbs (at_l, t_l);
			if (EGlpNumIsLess (at_l, t_max))
				EGlpNumCopy (t_max, at_l);
		}
		if (EGlpNumIsNeqq (*u, INFTY))
		{
			EGlpNumCopyDiffRatio (t_u, *x, *u, y_ij);
			EGlpNumCopyAbs (at_u, t_u);
			if (EGlpNumIsLess (at_u, t_max))
				EGlpNumCopy (t_max, at_u);
		}
	}

	if (EGlpNumIsLeq (INFTY, t_max))
	{
		rs->ratio_stat = RATIO_UNBOUNDED;
		ILL_CLEANUP;
	}

	EGlpNumZero (yi_max);
	EGlpNumZero (ayi_max);
	EGlpNumMultUiTo (t_max, 101);
	EGlpNumDivUiTo (t_max, 100);
	for (k = 0; k < lp->yjz.nzcnt; k++)
	{
		EGlpNumCopy (y_ij, lp->yjz.coef[k]);
		EGlpNumCopyAbs (ay_ij, y_ij);
		if (!EGlpNumIsNeqZero (y_ij, *pivtol))
			continue;

		i = lp->yjz.indx[k];
		if (lp->iwork[lp->baz[i]] == 1)
			continue;
		x = &(lp->xbz[i]);
		col = lp->baz[i];
		l = &(lp->lz[col]);
		u = &(lp->uz[col]);

		EGlpNumCopy (t_u, INFTY);
		EGlpNumCopy (at_u, t_u);
		EGlpNumCopy (t_l, NINFTY);
		EGlpNumCopy (at_l, t_u);
		if (EGlpNumIsNeqq (*l, NINFTY))
		{
			EGlpNumCopyDiffRatio (t_l, *x, *l, y_ij);
			EGlpNumCopyAbs (at_l, t_l);
		}
		if (EGlpNumIsNeqq (*u, INFTY))
		{
			EGlpNumCopyDiffRatio (t_u, *x, *u, y_ij);
			EGlpNumCopyAbs (at_u, t_u);
		}
		//t_i = (fabs (t_l) < fabs (t_u)) ? t_l : t_u;
		if (EGlpNumIsLess (at_l, at_u))
		{
			EGlpNumCopy (t_i, t_l);
			EGlpNumCopy (at_i, at_l);
		}
		else
		{
			EGlpNumCopy (t_i, t_u);
			EGlpNumCopy (at_i, at_u);
		}
		/*if (fabs (t_i) <= t_max + t_max * (1.0e-2)) */
		if (EGlpNumIsLeq (at_i, t_max))
		{
			if (EGlpNumIsLess (ayi_max, ay_ij))
			{
				EGlpNumCopy (yi_max, y_ij);
				EGlpNumCopy (ayi_max, ay_ij);
				rs->lindex = i;
				EGlpNumCopy (rs->tz, t_i);
				rs->lvstat = (EGlpNumIsLess (at_l, at_u)) ? STAT_LOWER : STAT_UPPER;
			}
		}
	}

	if (rs->lindex < 0)
	{
		rs->ratio_stat = RATIO_FAILED;
	}
	else
	{
		rs->ratio_stat = RATIO_BCHANGE;
		EGlpNumCopy (rs->pivotval, yi_max);
	}
CLEANUP:
	for (i = 0; i < rcnt; i++)
		lp->iwork[rlist[i]] = 0;
	EGlpNumClearVar (t_max);
	EGlpNumClearVar (ay_ij);
	EGlpNumClearVar (at_i);
	EGlpNumClearVar (at_l);
	EGlpNumClearVar (at_u);
	EGlpNumClearVar (ayi_max);
	EGlpNumClearVar (y_ij);
	EGlpNumClearVar (t_i);
	EGlpNumClearVar (t_l);
	EGlpNumClearVar (t_u);
	EGlpNumClearVar (yi_max);
	return;
}
