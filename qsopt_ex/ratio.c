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

#include "logging-private.h"

#include "eg_lpnum.h"
#include "eg_io.h"
#include "except.h"
#include "trace.h"

#include "sortrus_EGLPNUM_TYPENAME.h"
#include "stddefs.h"
#include "lpdefs_EGLPNUM_TYPENAME.h"
#include "ratio_EGLPNUM_TYPENAME.h"
#include "fct_EGLPNUM_TYPENAME.h"


void EGLPNUM_TYPENAME_ILLratio_pI_test (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int eindex,
	int dir,
	EGLPNUM_TYPENAME_ratio_res * rs)
{
	int i = 0, k = 0;
	int col, ecol;
	int cbnd, indx = 0;
	int tctr = 0;
	int *perm = lp->upd.perm;
	int *ix = lp->upd.ix;
	EGLPNUM_TYPE *pivtol = &(lp->tol->pivot_tol);
	EGLPNUM_TYPE *dftol = &(lp->tol->id_tol);

	 /*HHH*/ EGLPNUM_TYPE * t = lp->upd.t;
	EGLPNUM_TYPE t_i, delta, y_ij, rcost, nrcost, ntmp;
	EGLPNUM_TYPE *x, *l, *u;

	 /*HHH*/ EGLPNUM_TYPENAME_EGlpNumInitVar (t_i);
	EGLPNUM_TYPENAME_EGlpNumInitVar (delta);
	EGLPNUM_TYPENAME_EGlpNumInitVar (y_ij);
	EGLPNUM_TYPENAME_EGlpNumInitVar (rcost);
	EGLPNUM_TYPENAME_EGlpNumInitVar (nrcost);
	EGLPNUM_TYPENAME_EGlpNumInitVar (ntmp);
	EGLPNUM_TYPENAME_EGlpNumZero (t_i);
	EGLPNUM_TYPENAME_EGlpNumZero (y_ij);
	EGLPNUM_TYPENAME_EGlpNumZero (delta);
	rs->lindex = -1;
	EGLPNUM_TYPENAME_EGlpNumZero (rs->tz);
	EGLPNUM_TYPENAME_EGlpNumZero (rs->pivotval);
	rs->ratio_stat = RATIO_FAILED;
	rs->lvstat = -1;
	ecol = lp->nbaz[eindex];
	ILL_IFTRACE2 ("%s:%d:%d:%d:%d", __func__, eindex, dir, ecol,
								(VBOUNDED == lp->vtype[ecol]));
	if (lp->vtype[ecol] == VBOUNDED)
	{
		EGLPNUM_TYPENAME_EGlpNumCopyDiff (t[0], lp->uz[ecol], lp->lz[ecol]);
		ix[0] = BBOUND;
		ILL_IFTRACE2 (":%d[%d](%la,%la,%la)\n", ix[tctr], tctr,
									EGLPNUM_TYPENAME_EGlpNumToLf (t[tctr]), EGLPNUM_TYPENAME_EGlpNumToLf (lp->uz[ecol]),
									EGLPNUM_TYPENAME_EGlpNumToLf (lp->lz[ecol]));
		tctr++;
	}
	ILL_IFTRACE2 (":%d", lp->yjz.nzcnt);
	for (k = 0; k < lp->yjz.nzcnt; k++)
	{
		EGLPNUM_TYPENAME_EGlpNumCopy (y_ij, lp->yjz.coef[k]);
		if (!EGLPNUM_TYPENAME_EGlpNumIsNeqZero (y_ij, *pivtol))
			continue;

		i = lp->yjz.indx[k];
		x = &(lp->xbz[i]);
		col = lp->baz[i];
		l = &(lp->lz[col]);
		u = &(lp->uz[col]);

		if ((dir == VINCREASE && EGLPNUM_TYPENAME_EGlpNumIsGreatZero (y_ij)) ||
				(dir == VDECREASE && EGLPNUM_TYPENAME_EGlpNumIsLessZero (y_ij)))
		{
			if (EGLPNUM_TYPENAME_EGlpNumIsLessZero (y_ij))
				EGLPNUM_TYPENAME_EGlpNumSign (y_ij);
			ILL_IFTRACE2 (":%d", lp->bfeas[i]);
			if (lp->bfeas[i] > 0)
			{
				EGLPNUM_TYPENAME_EGlpNumCopyDiffRatio (t[tctr], *x, *u, y_ij);
				ix[tctr] = 10 * k + BATOUPPER;
				ILL_IFTRACE2 (":%d[%d](%la)\n", ix[tctr], tctr, EGLPNUM_TYPENAME_EGlpNumToLf (t[tctr]));
				tctr++;
				if (EGLPNUM_TYPENAME_EGlpNumIsNeqq (*l, EGLPNUM_TYPENAME_NINFTY))
				{
					EGLPNUM_TYPENAME_EGlpNumCopyDiffRatio (t[tctr], *x, *l, y_ij);
					ix[tctr] = 10 * k + BATOLOWER;
					ILL_IFTRACE2 (":%d[%d](%la)\n", ix[tctr], tctr,
												EGLPNUM_TYPENAME_EGlpNumToLf (t[tctr]));
					tctr++;
				}
			}
			else if (lp->bfeas[i] == 0)
			{
				if (EGLPNUM_TYPENAME_EGlpNumIsNeqq (*l, EGLPNUM_TYPENAME_NINFTY))
				{
					EGLPNUM_TYPENAME_EGlpNumCopyDiffRatio (t[tctr], *x, *l, y_ij);
					ix[tctr] = 10 * k + BATOLOWER;
					ILL_IFTRACE2 (":%d[%d](%la)\n", ix[tctr], tctr,
												EGLPNUM_TYPENAME_EGlpNumToLf (t[tctr]));
					tctr++;
				}
			}
		}
		else if ((dir == VINCREASE && EGLPNUM_TYPENAME_EGlpNumIsLessZero (y_ij)) ||
						 (dir == VDECREASE && EGLPNUM_TYPENAME_EGlpNumIsGreatZero (y_ij)))
		{
			if (EGLPNUM_TYPENAME_EGlpNumIsLessZero (y_ij))
				EGLPNUM_TYPENAME_EGlpNumSign (y_ij);
			ILL_IFTRACE2 (":%d", lp->bfeas[i]);
			if (lp->bfeas[i] < 0)
			{
				EGLPNUM_TYPENAME_EGlpNumCopyDiffRatio (t[tctr], *l, *x, y_ij);
				ix[tctr] = 10 * k + BBTOLOWER;
				ILL_IFTRACE2 (":%d[%d](%la)\n", ix[tctr], tctr, EGLPNUM_TYPENAME_EGlpNumToLf (t[tctr]));
				tctr++;
				if (EGLPNUM_TYPENAME_EGlpNumIsNeqq (*u, EGLPNUM_TYPENAME_INFTY))
				{
					EGLPNUM_TYPENAME_EGlpNumCopyDiffRatio (t[tctr], *u, *x, y_ij);
					ix[tctr] = 10 * k + BBTOUPPER;
					ILL_IFTRACE2 (":%d[%d](%la)\n", ix[tctr], tctr,
												EGLPNUM_TYPENAME_EGlpNumToLf (t[tctr]));
					tctr++;
				}
			}
			else if (lp->bfeas[i] == 0)
			{
				if (EGLPNUM_TYPENAME_EGlpNumIsNeqq (*u, EGLPNUM_TYPENAME_INFTY))
				{
					EGLPNUM_TYPENAME_EGlpNumCopyDiffRatio (t[tctr], *u, *x, y_ij);
					ix[tctr] = 10 * k + BBTOUPPER;
					ILL_IFTRACE2 (":%d[%d](%la)\n", ix[tctr], tctr,
												EGLPNUM_TYPENAME_EGlpNumToLf (t[tctr]));
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
	EGLPNUM_TYPENAME_ILLutil_EGlpNum_perm_quicksort (perm, t, tctr);

	EGLPNUM_TYPENAME_EGlpNumZero (lp->upd.c_obj);
	EGLPNUM_TYPENAME_EGlpNumCopy (rcost, lp->pIdz[eindex]);
	ILL_IFTRACE2 ("\n%s:%d:%lf", __func__, tctr, EGLPNUM_TYPENAME_EGlpNumToLf (rcost));
	for (i = 0; i < tctr; i++)
	{
		EGLPNUM_TYPENAME_EGlpNumCopy (t_i, t[perm[i]]);
		EGLPNUM_TYPENAME_EGlpNumCopy (ntmp, t_i);
		EGLPNUM_TYPENAME_EGlpNumSubTo (ntmp, delta);
		EGLPNUM_TYPENAME_EGlpNumAddInnProdTo (lp->upd.c_obj, ntmp, rcost);
		EGLPNUM_TYPENAME_EGlpNumCopy (delta, t_i);
		ILL_IFTRACE2 (":%d:%lf", perm[i], EGLPNUM_TYPENAME_EGlpNumToLf (delta));
		 /*HHH*/ cbnd = ix[perm[i]] % 10;
		if (cbnd != BBOUND)
		{
			k = ix[perm[i]] / 10;
			EGLPNUM_TYPENAME_EGlpNumCopy (y_ij, lp->yjz.coef[k]);
			indx = lp->yjz.indx[k];
			ILL_IFTRACE2 (":%d", indx);
		}

		switch (cbnd)
		{
		case BBOUND:
			rs->ratio_stat = RATIO_NOBCHANGE;
			EGLPNUM_TYPENAME_EGlpNumCopy (rs->tz, t_i);
			if (dir != VINCREASE)
				EGLPNUM_TYPENAME_EGlpNumSign (rs->tz);
			ILL_CLEANUP;

		case BATOLOWER:
		case BATOUPPER:
			EGLPNUM_TYPENAME_EGlpNumAddTo (rcost, y_ij);
			break;
		case BBTOLOWER:
		case BBTOUPPER:
			EGLPNUM_TYPENAME_EGlpNumSubTo (rcost, y_ij);
			break;
		}
		EGLPNUM_TYPENAME_EGlpNumCopyNeg (nrcost, rcost);
		if ((dir == VINCREASE && EGLPNUM_TYPENAME_EGlpNumIsLeq (nrcost, *dftol)) ||
				(dir == VDECREASE && EGLPNUM_TYPENAME_EGlpNumIsLeq (rcost, *dftol)))
		{
			/* change 5 to -1 if t_i > 0 is required below */
			if (EGLPNUM_TYPENAME_EGlpNumIsLessZero (t_i) && i > 5)
			{
				/* QSlog("pIhell %.5f %d", t_i, i); */
				EGLPNUM_TYPENAME_EGlpNumDivUiTo (t_i, 2);
				rs->ratio_stat = RATIO_NEGATIVE;
				EGLPNUM_TYPENAME_EGlpNumZero (rs->tz);
				ILL_CLEANUP;
			}
			rs->lindex = indx;
			rs->ratio_stat = RATIO_BCHANGE;
			if (cbnd == BATOLOWER || cbnd == BBTOLOWER)
				rs->lvstat = STAT_LOWER;
			else
				rs->lvstat = STAT_UPPER;

			EGLPNUM_TYPENAME_EGlpNumCopy (rs->pivotval, y_ij);
			EGLPNUM_TYPENAME_EGlpNumCopy (rs->tz, t_i);
			if (dir != VINCREASE)
				EGLPNUM_TYPENAME_EGlpNumSign (rs->tz);
			ILL_CLEANUP;
		}
	}

CLEANUP:
	EGLPNUM_TYPENAME_ILLfct_update_counts (lp, CNT_PIPIV, 0, rs->pivotval);
	ILL_IFTRACE2 (":tctr %d:%d\n", tctr, rs->ratio_stat);
	lp->upd.tctr = tctr;
	lp->upd.i = i;
	EGLPNUM_TYPENAME_EGlpNumCopy (lp->upd.tz, t_i);
	EGLPNUM_TYPENAME_EGlpNumCopy (lp->upd.piv, rs->pivotval);
	if (dir == VDECREASE)
		EGLPNUM_TYPENAME_EGlpNumSign (lp->upd.c_obj);
	if (rs->lindex != -1)
		lp->upd.fs = lp->bfeas[rs->lindex];
	EGLPNUM_TYPENAME_EGlpNumClearVar (t_i);
	EGLPNUM_TYPENAME_EGlpNumClearVar (delta);
	EGLPNUM_TYPENAME_EGlpNumClearVar (y_ij);
	EGLPNUM_TYPENAME_EGlpNumClearVar (rcost);
	EGLPNUM_TYPENAME_EGlpNumClearVar (nrcost);
	EGLPNUM_TYPENAME_EGlpNumClearVar (ntmp);
}

void EGLPNUM_TYPENAME_ILLratio_pII_test (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int eindex,
	int dir,
	EGLPNUM_TYPENAME_ratio_res * rs)
{
	int i, k, indx, col, ecol;
	EGLPNUM_TYPE *x, *l, *u, t_max, ayi_max, yi_max, ay_ij, y_ij, t_i, t_z;
	EGLPNUM_TYPE *pivtol = &(lp->tol->pivot_tol);
	EGLPNUM_TYPE *pftol = &(lp->tol->pfeas_tol);

	EGLPNUM_TYPENAME_EGlpNumInitVar (y_ij);
	EGLPNUM_TYPENAME_EGlpNumInitVar (ay_ij);
	EGLPNUM_TYPENAME_EGlpNumInitVar (t_i);
	EGLPNUM_TYPENAME_EGlpNumInitVar (t_z);
	EGLPNUM_TYPENAME_EGlpNumInitVar (t_max);
	EGLPNUM_TYPENAME_EGlpNumInitVar (yi_max);
	EGLPNUM_TYPENAME_EGlpNumInitVar (ayi_max);
	 /*HHH*/ rs->boundch = 0;
	rs->lindex = -1;
	EGLPNUM_TYPENAME_EGlpNumZero (rs->tz);
	rs->ratio_stat = RATIO_FAILED;
	rs->lvstat = -1;
	EGLPNUM_TYPENAME_EGlpNumZero (rs->pivotval);
	EGLPNUM_TYPENAME_EGlpNumZero (rs->lbound);
	ecol = lp->nbaz[eindex];

	for (k = 0, EGLPNUM_TYPENAME_EGlpNumCopy (t_max, EGLPNUM_TYPENAME_INFTY); k < lp->yjz.nzcnt; k++)
	{
		EGLPNUM_TYPENAME_EGlpNumCopy (y_ij, lp->yjz.coef[k]);
		EGLPNUM_TYPENAME_EGlpNumCopyAbs (ay_ij, y_ij);
		if (!EGLPNUM_TYPENAME_EGlpNumIsNeqZero (y_ij, *pivtol))
			continue;

		EGLPNUM_TYPENAME_EGlpNumCopy (t_i, EGLPNUM_TYPENAME_INFTY);
		i = lp->yjz.indx[k];
		x = &(lp->xbz[i]);
		col = lp->baz[i];
		l = &(lp->lz[col]);
		u = &(lp->uz[col]);

		if ((dir == VINCREASE && EGLPNUM_TYPENAME_EGlpNumIsGreatZero (y_ij)) ||
				(dir == VDECREASE && EGLPNUM_TYPENAME_EGlpNumIsLessZero (y_ij)))
		{
			if (EGLPNUM_TYPENAME_EGlpNumIsNeqq (*l, EGLPNUM_TYPENAME_NINFTY))
			{
				EGLPNUM_TYPENAME_EGlpNumCopyDiff (t_i, *x, *l);
				EGLPNUM_TYPENAME_EGlpNumAddTo (t_i, *pftol);
				EGLPNUM_TYPENAME_EGlpNumDivTo (t_i, ay_ij);
			}
		}
		else if ((dir == VINCREASE && EGLPNUM_TYPENAME_EGlpNumIsLessZero (y_ij)) ||
						 (dir == VDECREASE && EGLPNUM_TYPENAME_EGlpNumIsGreatZero (y_ij)))
		{
			if (EGLPNUM_TYPENAME_EGlpNumIsNeqq (*u, EGLPNUM_TYPENAME_INFTY))
			{
				EGLPNUM_TYPENAME_EGlpNumCopySum (t_i, *u, *pftol);
				EGLPNUM_TYPENAME_EGlpNumSubTo (t_i, *x);
				EGLPNUM_TYPENAME_EGlpNumDivTo (t_i, ay_ij);
			}
		}
		if (EGLPNUM_TYPENAME_EGlpNumIsEqqual (t_i, EGLPNUM_TYPENAME_INFTY))
			continue;

		if (EGLPNUM_TYPENAME_EGlpNumIsLess (t_i, t_max))
		{
			/*HHH tind = i; yval = fabs (y_ij); tval = t_i - pftol/fabs(y_ij); */
			EGLPNUM_TYPENAME_EGlpNumCopy (t_max, t_i);
		}
	}
	/* we use yi_max as temporal variable here */
	EGLPNUM_TYPENAME_EGlpNumCopyDiff (yi_max, lp->uz[ecol], lp->lz[ecol]);
	if (lp->vtype[ecol] == VBOUNDED && EGLPNUM_TYPENAME_EGlpNumIsLeq (yi_max, t_max))
	{

		EGLPNUM_TYPENAME_EGlpNumCopy (t_max, yi_max);
		rs->ratio_stat = RATIO_NOBCHANGE;
		EGLPNUM_TYPENAME_EGlpNumCopy (rs->tz, t_max);
		if (dir != VINCREASE)
			EGLPNUM_TYPENAME_EGlpNumSign (rs->tz);
		ILL_CLEANUP;
	}

	if (EGLPNUM_TYPENAME_EGlpNumIsLeq (EGLPNUM_TYPENAME_INFTY, t_max))
	{
		rs->ratio_stat = RATIO_UNBOUNDED;
		ILL_CLEANUP;
	}
	/*if (EGLPNUM_TYPENAME_EGlpNumIsLess (t_max, EGLPNUM_TYPENAME_zeroLpNum))
	 * QSlog("pIIhell");
	 */
	indx = -1;
	EGLPNUM_TYPENAME_EGlpNumZero (t_z);
	EGLPNUM_TYPENAME_EGlpNumZero (yi_max);
	EGLPNUM_TYPENAME_EGlpNumZero (ayi_max);
	ILL_IFTRACE2 (":%d", lp->yjz.nzcnt);
	for (k = 0; k < lp->yjz.nzcnt; k++)
	{
		EGLPNUM_TYPENAME_EGlpNumCopy (y_ij, lp->yjz.coef[k]);
		EGLPNUM_TYPENAME_EGlpNumCopyAbs (ay_ij, y_ij);
		if (!EGLPNUM_TYPENAME_EGlpNumIsNeqZero (y_ij, *pivtol))
			continue;

		EGLPNUM_TYPENAME_EGlpNumCopy (t_i, EGLPNUM_TYPENAME_INFTY);
		i = lp->yjz.indx[k];
		x = &(lp->xbz[i]);
		col = lp->baz[i];
		l = &(lp->lz[col]);
		u = &(lp->uz[col]);

		if ((dir == VINCREASE && EGLPNUM_TYPENAME_EGlpNumIsGreatZero (y_ij)) ||
				(dir == VDECREASE && EGLPNUM_TYPENAME_EGlpNumIsLessZero (y_ij)))
		{
			if (EGLPNUM_TYPENAME_EGlpNumIsNeqq (*l, EGLPNUM_TYPENAME_NINFTY))
				EGLPNUM_TYPENAME_EGlpNumCopyDiffRatio (t_i, *x, *l, ay_ij);
		}
		else if ((dir == VINCREASE && EGLPNUM_TYPENAME_EGlpNumIsLessZero (y_ij)) ||
						 (dir == VDECREASE && EGLPNUM_TYPENAME_EGlpNumIsGreatZero (y_ij)))
		{
			if (EGLPNUM_TYPENAME_EGlpNumIsNeqq (*u, EGLPNUM_TYPENAME_INFTY))
				EGLPNUM_TYPENAME_EGlpNumCopyDiffRatio (t_i, *u, *x, ay_ij);
		}

		if (EGLPNUM_TYPENAME_EGlpNumIsLeq (t_i, t_max))
		{
			if (EGLPNUM_TYPENAME_EGlpNumIsLess (ayi_max, ay_ij))
			{
				EGLPNUM_TYPENAME_EGlpNumCopy (yi_max, y_ij);
				EGLPNUM_TYPENAME_EGlpNumCopy (ayi_max, ay_ij);
				indx = i;
				EGLPNUM_TYPENAME_EGlpNumCopy (t_z, t_i);
				ILL_IFTRACE2 (":%d:%lf:%lf:%lf:%lf", indx, EGLPNUM_TYPENAME_EGlpNumToLf (t_i),
											EGLPNUM_TYPENAME_EGlpNumToLf (t_max), EGLPNUM_TYPENAME_EGlpNumToLf (ayi_max),
											EGLPNUM_TYPENAME_EGlpNumToLf (ay_ij));
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
		EGLPNUM_TYPENAME_EGlpNumCopy (rs->tz, t_z);
		EGLPNUM_TYPENAME_EGlpNumCopy (rs->pivotval, yi_max);
		rs->ratio_stat = RATIO_BCHANGE;

		if (dir == VINCREASE)
			rs->lvstat =
				(EGLPNUM_TYPENAME_EGlpNumIsGreatZero (yi_max)) ? STAT_LOWER : STAT_UPPER;
		else
			rs->lvstat =
				(EGLPNUM_TYPENAME_EGlpNumIsGreatZero (yi_max)) ? STAT_UPPER : STAT_LOWER;

		if (EGLPNUM_TYPENAME_EGlpNumIsLessZero (rs->tz))
		{
			ILL_IFTRACE2 ("need to change bound, tz=%la\n", EGLPNUM_TYPENAME_EGlpNumToLf (rs->tz));
			EGLPNUM_TYPENAME_EGlpNumCopyAbs (rs->tz, t_max);
			EGLPNUM_TYPENAME_EGlpNumDivUiTo (rs->tz, 10);
			rs->boundch = 1;
			EGLPNUM_TYPENAME_EGlpNumCopy (rs->lbound, lp->xbz[rs->lindex]);
			if (rs->lvstat == STAT_LOWER)
				EGLPNUM_TYPENAME_EGlpNumSubInnProdTo (rs->lbound, rs->tz, ayi_max);
			else
				EGLPNUM_TYPENAME_EGlpNumAddInnProdTo (rs->lbound, rs->tz, ayi_max);
		}
		if (dir == VDECREASE)
			EGLPNUM_TYPENAME_EGlpNumSign (rs->tz);
	}
CLEANUP:
	EGLPNUM_TYPENAME_ILLfct_update_counts (lp, CNT_PIIPIV, 0, rs->pivotval);
	EGLPNUM_TYPENAME_EGlpNumClearVar (y_ij);
	EGLPNUM_TYPENAME_EGlpNumClearVar (ay_ij);
	EGLPNUM_TYPENAME_EGlpNumClearVar (t_i);
	EGLPNUM_TYPENAME_EGlpNumClearVar (t_z);
	EGLPNUM_TYPENAME_EGlpNumClearVar (t_max);
	EGLPNUM_TYPENAME_EGlpNumClearVar (yi_max);
	EGLPNUM_TYPENAME_EGlpNumClearVar (ayi_max);
}

#define GET_XY_DRATIOTEST \
      if (lp->vstat[col] == STAT_UPPER){ \
				EGLPNUM_TYPENAME_EGlpNumCopyNeg(x,lp->dz[j]);\
        EGLPNUM_TYPENAME_EGlpNumCopy(y, *zAj);\
      } \
      else{ \
         EGLPNUM_TYPENAME_EGlpNumCopy(x, lp->dz[j]); \
         EGLPNUM_TYPENAME_EGlpNumCopyNeg(y, *zAj);\
      } \
      if (lvstat == STAT_UPPER) \
         EGLPNUM_TYPENAME_EGlpNumSign(y);


void EGLPNUM_TYPENAME_ILLratio_dI_test (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int lindex,
	int lvstat,
	EGLPNUM_TYPENAME_ratio_res * rs)
{
	int j = 0, k;
	int col;
	int cbnd, indx;
	int tctr = 0;
	int *perm = lp->upd.perm;
	int *ix = lp->upd.ix;
	EGLPNUM_TYPE *t = lp->upd.t;
	EGLPNUM_TYPE *zAj, x, y, t_j, theta, rcost, delta;
	EGLPNUM_TYPE *pftol = &(lp->tol->ip_tol);
	EGLPNUM_TYPE *pivtol = &(lp->tol->pivot_tol);

	EGLPNUM_TYPENAME_EGlpNumInitVar (x);
	EGLPNUM_TYPENAME_EGlpNumInitVar (y);
	EGLPNUM_TYPENAME_EGlpNumInitVar (t_j);
	EGLPNUM_TYPENAME_EGlpNumInitVar (theta);
	EGLPNUM_TYPENAME_EGlpNumInitVar (rcost);
	EGLPNUM_TYPENAME_EGlpNumInitVar (delta);
	EGLPNUM_TYPENAME_EGlpNumZero (delta);
	EGLPNUM_TYPENAME_EGlpNumZero (t_j);
	EGLPNUM_TYPENAME_EGlpNumZero (rs->tz);
	 /*HHH*/ rs->eindex = -1;
	rs->ratio_stat = RATIO_FAILED;
	EGLPNUM_TYPENAME_EGlpNumZero (rs->pivotval);

	for (k = 0; k < lp->zA.nzcnt; k++)
	{
		zAj = &(lp->zA.coef[k]);
		if (!EGLPNUM_TYPENAME_EGlpNumIsNeqZero (*zAj, *pivtol))
			continue;

		EGLPNUM_TYPENAME_EGlpNumCopy (t_j, EGLPNUM_TYPENAME_INFTY);
		j = lp->zA.indx[k];
		col = lp->nbaz[j];

		if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED)
			continue;

		GET_XY_DRATIOTEST;

		if (EGLPNUM_TYPENAME_EGlpNumIsLessZero (y))
		{
			if (lp->dfeas[j] != 0 && lp->vstat[col] != STAT_ZERO)
			{
				EGLPNUM_TYPENAME_EGlpNumCopyFrac (t[tctr], x, y);
				ix[tctr] = 10 * k + BBTOLOWER;
				tctr++;
			}
			else if (lp->vstat[col] == STAT_ZERO)
			{
				if (lp->dfeas[j] < 0)
				{
					EGLPNUM_TYPENAME_EGlpNumCopyFrac (t[tctr], x, y);
					ix[tctr] = 10 * k + BBTOLOWER;
					tctr++;
				}
				if (lp->dfeas[j] <= 0)
				{
					EGLPNUM_TYPENAME_EGlpNumCopyFrac (t[tctr], x, y);
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
					EGLPNUM_TYPENAME_EGlpNumCopyFrac (t[tctr], x, y);
					ix[tctr] = 10 * k + BATOUPPER;
					tctr++;
					EGLPNUM_TYPENAME_EGlpNumCopyFrac (t[tctr], x, y);
					ix[tctr] = 10 * k + BATOLOWER;
					tctr++;
				}
			}
			else if (lp->dfeas[j] == 0)
			{
				EGLPNUM_TYPENAME_EGlpNumCopyFrac (t[tctr], x, y);
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
	EGLPNUM_TYPENAME_ILLutil_EGlpNum_perm_quicksort (perm, t, tctr);

	EGLPNUM_TYPENAME_EGlpNumZero (lp->upd.c_obj);
	EGLPNUM_TYPENAME_EGlpNumCopy (rcost, lp->xbz[lindex]);
	if (lvstat == STAT_LOWER)
		EGLPNUM_TYPENAME_EGlpNumSign (rcost);
	for (j = 0; j < tctr; j++)
	{
		cbnd = ix[perm[j]] % 10;
		if (cbnd == BSKIP)
			continue;

		EGLPNUM_TYPENAME_EGlpNumCopy (t_j, t[perm[j]]);
		EGLPNUM_TYPENAME_EGlpNumCopy (x, t_j);
		EGLPNUM_TYPENAME_EGlpNumSubTo (x, delta);
		EGLPNUM_TYPENAME_EGlpNumAddInnProdTo (lp->upd.c_obj, x, rcost);
		EGLPNUM_TYPENAME_EGlpNumCopy (delta, t_j);
		k = ix[perm[j]] / 10;
		zAj = &(lp->zA.coef[k]);
		indx = lp->zA.indx[k];

		if (lp->vstat[lp->nbaz[indx]] == STAT_LOWER
				|| lp->vstat[lp->nbaz[indx]] == STAT_ZERO)
			EGLPNUM_TYPENAME_EGlpNumCopyNeg (theta, *zAj);
		else
			EGLPNUM_TYPENAME_EGlpNumCopy (theta, *zAj);

		if (lvstat == STAT_UPPER)
			EGLPNUM_TYPENAME_EGlpNumSign (theta);

		switch (cbnd)
		{
		case BATOLOWER:
		case BATOUPPER:
			EGLPNUM_TYPENAME_EGlpNumSubTo (rcost, theta);
			break;
		case BBTOLOWER:
		case BBTOUPPER:
			EGLPNUM_TYPENAME_EGlpNumAddTo (rcost, theta);
			break;
		}
		if (EGLPNUM_TYPENAME_EGlpNumIsLeq (rcost, *pftol))
		{
			/* if (t_j < 0.0) QSlog("dIhell"); */
			rs->eindex = indx;
			EGLPNUM_TYPENAME_EGlpNumCopy (rs->tz, t_j);
			EGLPNUM_TYPENAME_EGlpNumCopy (rs->pivotval, *zAj);
			rs->ratio_stat = RATIO_BCHANGE;
			ILL_CLEANUP;
		}
	}

CLEANUP:
	EGLPNUM_TYPENAME_ILLfct_update_counts (lp, CNT_DIPIV, 0, rs->pivotval);
	ILL_IFTRACE2 ("%s:tctr %d\n", __func__, tctr);
	lp->upd.tctr = tctr;
	lp->upd.i = j;
	EGLPNUM_TYPENAME_EGlpNumCopyAbs (lp->upd.tz, t_j);
	EGLPNUM_TYPENAME_EGlpNumCopy (lp->upd.piv, rs->pivotval);
	if (rs->eindex != -1)
		lp->upd.fs = lp->dfeas[rs->eindex];
	EGLPNUM_TYPENAME_EGlpNumClearVar (x);
	EGLPNUM_TYPENAME_EGlpNumClearVar (y);
	EGLPNUM_TYPENAME_EGlpNumClearVar (t_j);
	EGLPNUM_TYPENAME_EGlpNumClearVar (theta);
	EGLPNUM_TYPENAME_EGlpNumClearVar (rcost);
	EGLPNUM_TYPENAME_EGlpNumClearVar (delta);
}

void EGLPNUM_TYPENAME_ILLratio_dII_test (
	EGLPNUM_TYPENAME_lpinfo * lp,
	/*int lindex,*/
	int lvstat,
	EGLPNUM_TYPENAME_ratio_res * rs)
{
	int j, k, indx;
	int col, ecol;
	EGLPNUM_TYPE *zAj, azAj, az_max, x, y, t_j, z_max, t_max, t_z;
	EGLPNUM_TYPE *dftol = &(lp->tol->dfeas_tol);
	EGLPNUM_TYPE *pivtol = &(lp->tol->pivot_tol);

	EGLPNUM_TYPENAME_EGlpNumInitVar (x);
	EGLPNUM_TYPENAME_EGlpNumInitVar (y);
	EGLPNUM_TYPENAME_EGlpNumInitVar (t_j);
	EGLPNUM_TYPENAME_EGlpNumInitVar (z_max);
	EGLPNUM_TYPENAME_EGlpNumInitVar (t_max);
	EGLPNUM_TYPENAME_EGlpNumInitVar (az_max);
	EGLPNUM_TYPENAME_EGlpNumInitVar (azAj);
	EGLPNUM_TYPENAME_EGlpNumInitVar (t_z);
	EGLPNUM_TYPENAME_EGlpNumZero (t_j);
	rs->coeffch = 0;
	EGLPNUM_TYPENAME_EGlpNumZero (rs->ecoeff);
	rs->eindex = -1;
	rs->ratio_stat = RATIO_FAILED;
	ILL_IFTRACE2 ("%s:tctr %d\n", __func__, 0);
	lp->upd.tctr = 0;
	EGLPNUM_TYPENAME_EGlpNumZero (lp->upd.dty);
	for (k = 0, EGLPNUM_TYPENAME_EGlpNumCopy (t_max, EGLPNUM_TYPENAME_INFTY); k < lp->zA.nzcnt; k++)
	{
		zAj = &(lp->zA.coef[k]);
		if (!EGLPNUM_TYPENAME_EGlpNumIsNeqZero (*zAj, *pivtol))
			continue;

		EGLPNUM_TYPENAME_EGlpNumCopy (t_j, EGLPNUM_TYPENAME_INFTY);
		j = lp->zA.indx[k];
		col = lp->nbaz[j];

		if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED)
			continue;

		GET_XY_DRATIOTEST;

//#warning adding/substracting tolerances to used value, is it rigght?
		if (EGLPNUM_TYPENAME_EGlpNumIsGreatZero (y))
		{
			//t_j = (x + dftol) / y;
			EGLPNUM_TYPENAME_EGlpNumCopySum (t_j, x, *dftol);
			EGLPNUM_TYPENAME_EGlpNumDivTo (t_j, y);
		}
		else
		{
//#warning adding/substracting tolerances to used value, is it rigght?
			if (lp->vstat[col] == STAT_ZERO)
				EGLPNUM_TYPENAME_EGlpNumCopyDiffRatio (t_j, x, *dftol, y);
		}
		//if (t_j == EGLPNUM_TYPENAME_INFTY)
		if (EGLPNUM_TYPENAME_EGlpNumIsEqqual (t_j, EGLPNUM_TYPENAME_INFTY))
			continue;

		if (EGLPNUM_TYPENAME_EGlpNumIsLess (t_j, t_max))
			EGLPNUM_TYPENAME_EGlpNumCopy (t_max, t_j);
	}

	if (EGLPNUM_TYPENAME_EGlpNumIsLeq (EGLPNUM_TYPENAME_INFTY, t_max))
	{
		rs->ratio_stat = RATIO_UNBOUNDED;
		ILL_CLEANUP;
	}
	/* if (t_max < 0.0) QSlog("dIIhell"); */

	indx = -1;
	EGLPNUM_TYPENAME_EGlpNumZero (t_z);
	EGLPNUM_TYPENAME_EGlpNumZero (z_max);
	EGLPNUM_TYPENAME_EGlpNumZero (az_max);

	for (k = 0; k < lp->zA.nzcnt; k++)
	{
		zAj = &(lp->zA.coef[k]);
		EGLPNUM_TYPENAME_EGlpNumCopyAbs (azAj, *zAj);
		if (!EGLPNUM_TYPENAME_EGlpNumIsNeqZero (*zAj, *pivtol))
			continue;

		EGLPNUM_TYPENAME_EGlpNumCopy (t_j, EGLPNUM_TYPENAME_INFTY);
		j = lp->zA.indx[k];
		col = lp->nbaz[j];

		if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED)
			continue;

		GET_XY_DRATIOTEST;

		if (EGLPNUM_TYPENAME_EGlpNumIsGreatZero (y) || lp->vstat[col] == STAT_ZERO)
			EGLPNUM_TYPENAME_EGlpNumCopyFrac (t_j, x, y);

		if (EGLPNUM_TYPENAME_EGlpNumIsLeq (t_j, t_max) && (EGLPNUM_TYPENAME_EGlpNumIsLess (az_max, azAj)))
		{
			EGLPNUM_TYPENAME_EGlpNumCopy (z_max, *zAj);
			EGLPNUM_TYPENAME_EGlpNumCopy (az_max, azAj);
			indx = j;
			EGLPNUM_TYPENAME_EGlpNumCopy (t_z, t_j);
		}
	}


	if (indx < 0)
	{
		rs->ratio_stat = RATIO_FAILED;
	}
	else
	{
		rs->eindex = indx;
		EGLPNUM_TYPENAME_EGlpNumCopy (rs->tz, t_z);
		EGLPNUM_TYPENAME_EGlpNumCopy (rs->pivotval, z_max);
		rs->ratio_stat = RATIO_BCHANGE;

		if (EGLPNUM_TYPENAME_EGlpNumIsLessZero (rs->tz))
		{
			EGLPNUM_TYPENAME_EGlpNumCopyAbs (rs->tz, t_max);
			EGLPNUM_TYPENAME_EGlpNumDivUiTo (rs->tz, 20);
			rs->coeffch = 1;
			ecol = lp->nbaz[indx];
			EGLPNUM_TYPENAME_EGlpNumCopyDiff (rs->ecoeff, lp->cz[ecol], lp->dz[indx]);
			switch (lp->vstat[ecol])
			{
			case STAT_LOWER:
				EGLPNUM_TYPENAME_EGlpNumAddInnProdTo (rs->ecoeff, rs->tz, az_max);
				break;
			case STAT_UPPER:
				EGLPNUM_TYPENAME_EGlpNumSubInnProdTo (rs->ecoeff, rs->tz, az_max);
				break;
			default:
				EGLPNUM_TYPENAME_EGlpNumZero (rs->tz);
				break;
			}
		}
	}

CLEANUP:
	EGLPNUM_TYPENAME_ILLfct_update_counts (lp, CNT_DIIPIV, 0, rs->pivotval);
	EGLPNUM_TYPENAME_EGlpNumCopy (lp->upd.piv, rs->pivotval);
	EGLPNUM_TYPENAME_EGlpNumClearVar (x);
	EGLPNUM_TYPENAME_EGlpNumClearVar (y);
	EGLPNUM_TYPENAME_EGlpNumClearVar (t_j);
	EGLPNUM_TYPENAME_EGlpNumClearVar (z_max);
	EGLPNUM_TYPENAME_EGlpNumClearVar (t_max);
	EGLPNUM_TYPENAME_EGlpNumClearVar (t_z);
	EGLPNUM_TYPENAME_EGlpNumClearVar (az_max);
	EGLPNUM_TYPENAME_EGlpNumClearVar (azAj);
}

void EGLPNUM_TYPENAME_ILLratio_longdII_test (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int lindex,
	int lvstat,
	EGLPNUM_TYPENAME_ratio_res * rs)
{
	int j, k, indx = 0, tctr = 0;
	int col, ecol;
	int vs, bnd_exist = 0;
	int *perm = lp->upd.perm;
	int *ix = lp->upd.ix;
	int b_indx = -1;
	EGLPNUM_TYPE *t = lp->upd.t;
	EGLPNUM_TYPE *l,
		*u,
		*xb,
		*zAj = 0,
		x,
		y,
		t_j,
		z_max,
		t_max, t_z, theta, rcost, delta, zb_val, tb_val, az_max, azb_val, azAj;
	EGLPNUM_TYPE *pftol = &(lp->tol->pfeas_tol);
	EGLPNUM_TYPE *dftol = &(lp->tol->dfeas_tol);
	EGLPNUM_TYPE *pivtol = &(lp->tol->pivot_tol);

	EGLPNUM_TYPENAME_EGlpNumInitVar (x);
	EGLPNUM_TYPENAME_EGlpNumInitVar (azAj);
	EGLPNUM_TYPENAME_EGlpNumInitVar (y);
	EGLPNUM_TYPENAME_EGlpNumInitVar (t_j);
	EGLPNUM_TYPENAME_EGlpNumInitVar (z_max);
	EGLPNUM_TYPENAME_EGlpNumInitVar (az_max);
	EGLPNUM_TYPENAME_EGlpNumInitVar (t_max);
	EGLPNUM_TYPENAME_EGlpNumInitVar (t_z);
	EGLPNUM_TYPENAME_EGlpNumInitVar (theta);
	EGLPNUM_TYPENAME_EGlpNumInitVar (rcost);
	EGLPNUM_TYPENAME_EGlpNumInitVar (delta);
	EGLPNUM_TYPENAME_EGlpNumInitVar (zb_val);
	EGLPNUM_TYPENAME_EGlpNumInitVar (azb_val);
	EGLPNUM_TYPENAME_EGlpNumInitVar (tb_val);
	EGLPNUM_TYPENAME_EGlpNumZero (t_j);
	EGLPNUM_TYPENAME_EGlpNumZero (delta);
	EGLPNUM_TYPENAME_EGlpNumZero (zb_val);
	EGLPNUM_TYPENAME_EGlpNumZero (azb_val);
	EGLPNUM_TYPENAME_EGlpNumCopy (tb_val, EGLPNUM_TYPENAME_NINFTY);
//#warning not sure about THIS line
	EGLPNUM_TYPENAME_EGlpNumZero (rs->pivotval);

	rs->coeffch = 0;
	rs->eindex = -1;
	rs->ratio_stat = RATIO_FAILED;

	ILL_IFTRACE2 ("%s:tctr %d\n", __func__, 0);
	lp->upd.tctr = 0;
	lp->upd.i = 0;
	EGLPNUM_TYPENAME_EGlpNumZero (lp->upd.tz);
	EGLPNUM_TYPENAME_EGlpNumZero (lp->upd.piv);
	EGLPNUM_TYPENAME_EGlpNumZero (lp->upd.c_obj);
	EGLPNUM_TYPENAME_EGlpNumZero (lp->upd.dty);

	xb = &(lp->xbz[lindex]);
	col = lp->baz[lindex];
	l = &(lp->lz[col]);
	u = &(lp->uz[col]);
	//rcost = (lvstat == STAT_LOWER) ? l - xb : xb - u;
	if (lvstat == STAT_LOWER)
		EGLPNUM_TYPENAME_EGlpNumCopyDiff (rcost, *l, *xb);
	else
		EGLPNUM_TYPENAME_EGlpNumCopyDiff (rcost, *xb, *u);

	for (k = 0, EGLPNUM_TYPENAME_EGlpNumCopy (t_max, EGLPNUM_TYPENAME_INFTY); k < lp->zA.nzcnt; k++)
	{
		zAj = &(lp->zA.coef[k]);
		if (!EGLPNUM_TYPENAME_EGlpNumIsNeqZero (*zAj, *pivtol))
			continue;

		EGLPNUM_TYPENAME_EGlpNumCopy (t_j, EGLPNUM_TYPENAME_INFTY);
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

		if (EGLPNUM_TYPENAME_EGlpNumIsGreatZero (y))
		{
			//t_j = (x + dftol) / y;
//#warning Using tolerances to add to result, is it right?
			EGLPNUM_TYPENAME_EGlpNumCopySum (t_j, x, *dftol);
			EGLPNUM_TYPENAME_EGlpNumDivTo (t_j, y);
		}
		else
		{
			if (lp->vstat[col] == STAT_ZERO)
				EGLPNUM_TYPENAME_EGlpNumCopyDiffRatio (t_j, x, *dftol, y);
		}
		if (EGLPNUM_TYPENAME_EGlpNumIsEqqual (t_j, EGLPNUM_TYPENAME_INFTY))
			continue;

		if (EGLPNUM_TYPENAME_EGlpNumIsLess (t_j, t_max))
			EGLPNUM_TYPENAME_EGlpNumCopy (t_max, t_j);
	}
	if (EGLPNUM_TYPENAME_EGlpNumIsLessZero (t_max))
	{
		/*QSlog("dIIhell, %.4f", t_max); */
		rs->ratio_stat = RATIO_NEGATIVE;
		ILL_CLEANUP;
	}

	if (bnd_exist == 0 && EGLPNUM_TYPENAME_EGlpNumIsLeq (EGLPNUM_TYPENAME_INFTY, t_max))
	{
		rs->ratio_stat = RATIO_UNBOUNDED;
		/*
		 * QSlog("x = %.8f, b = %.2f", lp->xbz[lindex], (lvstat == STAT_LOWER ) ? lp->lz[lp->baz[lindex]] : lp->uz[lp->baz[lindex]]);
		 */
		ILL_CLEANUP;
	}

	if (bnd_exist != 0)
	{
		for (k = 0; k < lp->zA.nzcnt; k++)
		{
			zAj = &(lp->zA.coef[k]);
			if (!EGLPNUM_TYPENAME_EGlpNumIsNeqZero (*zAj, *pivtol))
				continue;

			EGLPNUM_TYPENAME_EGlpNumCopy (t_j, EGLPNUM_TYPENAME_INFTY);
			j = lp->zA.indx[k];
			col = lp->nbaz[j];

			if (lp->vtype[col] != VBOUNDED)
				continue;

			GET_XY_DRATIOTEST;

			if (EGLPNUM_TYPENAME_EGlpNumIsGreatZero (y))
			{
				EGLPNUM_TYPENAME_EGlpNumCopyFrac (t_j, x, y);
				if (EGLPNUM_TYPENAME_EGlpNumIsLeq (t_j, t_max))
				{
					EGLPNUM_TYPENAME_EGlpNumCopy (t[tctr], t_j);
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
		EGLPNUM_TYPENAME_ILLutil_EGlpNum_perm_quicksort (perm, t, tctr);

		for (j = 0; j < tctr; j++)
		{

			EGLPNUM_TYPENAME_EGlpNumCopy (t_j, t[perm[j]]);
			/* we use x as temporal storage */
			//lp->upd.c_obj += (t_j - delta) * rcost;
			EGLPNUM_TYPENAME_EGlpNumCopy (x, t_j);
			EGLPNUM_TYPENAME_EGlpNumSubTo (x, delta);
			EGLPNUM_TYPENAME_EGlpNumAddInnProdTo (lp->upd.c_obj, x, rcost);
			EGLPNUM_TYPENAME_EGlpNumCopy (delta, t_j);
			 /*HHH*/ k = ix[perm[j]];
			zAj = &(lp->zA.coef[k]);
			indx = lp->zA.indx[k];
			col = lp->nbaz[indx];
			l = &(lp->lz[col]);
			u = &(lp->uz[col]);
			vs = lp->vstat[col];
			//theta = (vs == STAT_UPPER) ? (l - u) * zAj : (u - l) * zAj;
			EGLPNUM_TYPENAME_EGlpNumCopyDiff (theta, *l, *u);
			EGLPNUM_TYPENAME_EGlpNumMultTo (theta, *zAj);
			if (vs != STAT_UPPER)
				EGLPNUM_TYPENAME_EGlpNumSign (theta);
			if (lvstat == STAT_LOWER)
				EGLPNUM_TYPENAME_EGlpNumAddTo (rcost, theta);
			else
				EGLPNUM_TYPENAME_EGlpNumSubTo (rcost, theta);

			if (EGLPNUM_TYPENAME_EGlpNumIsLeq (rcost, *pftol))
			{
				rs->eindex = indx;
				EGLPNUM_TYPENAME_EGlpNumCopy (rs->tz, t_j);
				EGLPNUM_TYPENAME_EGlpNumCopy (rs->pivotval, *zAj);
				rs->ratio_stat = RATIO_BCHANGE;

				if (EGLPNUM_TYPENAME_EGlpNumIsLessZero (rs->tz))
				{
					EGLPNUM_TYPENAME_EGlpNumZero (rs->tz);
					rs->coeffch = 1;
					//rs->ecoeff = lp->cz[col] - lp->dz[indx];
					EGLPNUM_TYPENAME_EGlpNumCopyDiff (rs->ecoeff, lp->cz[col], lp->dz[indx]);
					//lp->upd.c_obj += (rs->tz - delta) * rcost; note ts->tz == 0;
					EGLPNUM_TYPENAME_EGlpNumSubInnProdTo (lp->upd.c_obj, delta, rcost);
				}
				ILL_IFTRACE2 ("%s:tctr %d\n", __func__, tctr);
				lp->upd.tctr = tctr;
				lp->upd.i = j;
				EGLPNUM_TYPENAME_EGlpNumCopy (lp->upd.tz, rs->tz);
				ILL_CLEANUP;
			}
		}
		ILL_IFTRACE2 ("%s:tctr %d\n", __func__, tctr);
		lp->upd.tctr = tctr;
		lp->upd.i = tctr;
		EGLPNUM_TYPENAME_EGlpNumCopy (lp->upd.tz, t_j);
		EGLPNUM_TYPENAME_EGlpNumCopy (zb_val, *zAj);
		EGLPNUM_TYPENAME_EGlpNumCopyAbs (azb_val, zb_val);
		EGLPNUM_TYPENAME_EGlpNumCopy (tb_val, t_j);
		b_indx = indx;
	}

	if (bnd_exist != 0 && EGLPNUM_TYPENAME_EGlpNumIsLeq (EGLPNUM_TYPENAME_INFTY, t_max))
	{
		rs->ratio_stat = RATIO_UNBOUNDED;
		/* QSlog("rcost: %.8f", rcost); */
		ILL_CLEANUP;
	}

	EGLPNUM_TYPENAME_EGlpNumZero (z_max);
	EGLPNUM_TYPENAME_EGlpNumZero (az_max);
	indx = -1;
	EGLPNUM_TYPENAME_EGlpNumZero (t_z);
	for (k = 0; k < lp->zA.nzcnt; k++)
	{
		zAj = &(lp->zA.coef[k]);
		EGLPNUM_TYPENAME_EGlpNumCopyAbs (azAj, *zAj);
		if (!EGLPNUM_TYPENAME_EGlpNumIsNeqZero (*zAj, *pivtol))
			continue;

		EGLPNUM_TYPENAME_EGlpNumCopy (t_j, EGLPNUM_TYPENAME_INFTY);
		j = lp->zA.indx[k];
		col = lp->nbaz[j];

		if (lp->vtype[col] == VARTIFICIAL || lp->vtype[col] == VFIXED ||
				lp->vtype[col] == VBOUNDED)
			continue;

		GET_XY_DRATIOTEST;

		if (EGLPNUM_TYPENAME_EGlpNumIsGreatZero (y) || lp->vstat[col] == STAT_ZERO)
			EGLPNUM_TYPENAME_EGlpNumCopyFrac (t_j, x, y);

		if (EGLPNUM_TYPENAME_EGlpNumIsLeq (t_j, t_max))
		{
			if (EGLPNUM_TYPENAME_EGlpNumIsLess (az_max, azAj))
			{
				EGLPNUM_TYPENAME_EGlpNumCopy (z_max, *zAj);
				EGLPNUM_TYPENAME_EGlpNumCopy (az_max, azAj);
				indx = j;
				EGLPNUM_TYPENAME_EGlpNumCopy (t_z, t_j);
			}
		}
	}

	if (indx < 0)
	{
		rs->ratio_stat = RATIO_FAILED;
		ILL_CLEANUP;
	}
	if ((tctr == 0) || (EGLPNUM_TYPENAME_EGlpNumIsLessZero (tb_val)) ||
			(tctr != 0 && EGLPNUM_TYPENAME_EGlpNumIsLeq (tb_val, t_z) &&
			 EGLPNUM_TYPENAME_EGlpNumIsLeq (azb_val, az_max)))
	{
		/* we use x as temporal vvariable */
		/* lp->upd.c_obj += (t_z - delta) * rcost; */
		EGLPNUM_TYPENAME_EGlpNumCopyDiff (x, t_z, delta);
		EGLPNUM_TYPENAME_EGlpNumAddInnProdTo (lp->upd.c_obj, x, rcost);
		EGLPNUM_TYPENAME_EGlpNumCopy (delta, t_z);
		rs->eindex = indx;
		EGLPNUM_TYPENAME_EGlpNumCopy (rs->tz, t_z);
		EGLPNUM_TYPENAME_EGlpNumCopy (rs->pivotval, z_max);
		rs->ratio_stat = RATIO_BCHANGE;
	}
	/* For now */
	else if (tctr != 0)
	{
		rs->eindex = b_indx;
		EGLPNUM_TYPENAME_EGlpNumCopy (rs->tz, tb_val);
		EGLPNUM_TYPENAME_EGlpNumCopy (rs->pivotval, zb_val);
		rs->ratio_stat = RATIO_BCHANGE;
		lp->upd.i -= 1;
	}

	if (EGLPNUM_TYPENAME_EGlpNumIsLessZero (rs->tz))
	{
		/* if (tctr != 0) QSlog("despite long step"); */
		/* rs->tz = fabs (t_max / 20.0); */
		EGLPNUM_TYPENAME_EGlpNumCopyAbs (rs->tz, t_max);
		EGLPNUM_TYPENAME_EGlpNumDivUiTo (rs->tz, 20);
		rs->coeffch = 1;

		ecol = lp->nbaz[indx];
		if (lp->vstat[ecol] == STAT_LOWER)
		{
			/*rs->ecoeff = lp->cz[ecol] - lp->dz[indx] + rs->tz * fabs (z_max); */
			EGLPNUM_TYPENAME_EGlpNumCopy (rs->ecoeff, az_max);
			EGLPNUM_TYPENAME_EGlpNumMultTo (rs->ecoeff, rs->tz);
			EGLPNUM_TYPENAME_EGlpNumAddTo (rs->ecoeff, lp->cz[ecol]);
			EGLPNUM_TYPENAME_EGlpNumSubTo (rs->ecoeff, lp->dz[indx]);
		}
		else if (lp->vstat[ecol] == STAT_UPPER)
		{
			/*rs->ecoeff = lp->cz[ecol] - lp->dz[indx] - rs->tz * fabs (z_max); */
			EGLPNUM_TYPENAME_EGlpNumCopy (rs->ecoeff, az_max);
			EGLPNUM_TYPENAME_EGlpNumMultTo (rs->ecoeff, rs->tz);
			EGLPNUM_TYPENAME_EGlpNumSign (rs->ecoeff);
			EGLPNUM_TYPENAME_EGlpNumAddTo (rs->ecoeff, lp->cz[ecol]);
			EGLPNUM_TYPENAME_EGlpNumSubTo (rs->ecoeff, lp->dz[indx]);
		}
		else
		{
			/*rs->ecoeff = lp->cz[ecol] - lp->dz[indx]; */
			EGLPNUM_TYPENAME_EGlpNumCopyDiff (rs->ecoeff, lp->cz[ecol], lp->dz[indx]);
			EGLPNUM_TYPENAME_EGlpNumZero (rs->tz);
		}
		/* we use x as temporal storage */
		/*lp->upd.c_obj += (rs->tz - delta) * rcost; */
		EGLPNUM_TYPENAME_EGlpNumCopy (x, rs->tz);
		EGLPNUM_TYPENAME_EGlpNumSubTo (x, delta);
		EGLPNUM_TYPENAME_EGlpNumAddInnProdTo (lp->upd.c_obj, x, rcost);
	}

CLEANUP:
	EGLPNUM_TYPENAME_ILLfct_update_counts (lp, CNT_DIIPIV, 0, rs->pivotval);
	EGLPNUM_TYPENAME_EGlpNumCopy (lp->upd.piv, rs->pivotval);
	EGLPNUM_TYPENAME_EGlpNumClearVar (x);
	EGLPNUM_TYPENAME_EGlpNumClearVar (y);
	EGLPNUM_TYPENAME_EGlpNumClearVar (t_j);
	EGLPNUM_TYPENAME_EGlpNumClearVar (z_max);
	EGLPNUM_TYPENAME_EGlpNumClearVar (az_max);
	EGLPNUM_TYPENAME_EGlpNumClearVar (t_max);
	EGLPNUM_TYPENAME_EGlpNumClearVar (t_z);
	EGLPNUM_TYPENAME_EGlpNumClearVar (theta);
	EGLPNUM_TYPENAME_EGlpNumClearVar (rcost);
	EGLPNUM_TYPENAME_EGlpNumClearVar (delta);
	EGLPNUM_TYPENAME_EGlpNumClearVar (zb_val);
	EGLPNUM_TYPENAME_EGlpNumClearVar (azb_val);
	EGLPNUM_TYPENAME_EGlpNumClearVar (tb_val);
	EGLPNUM_TYPENAME_EGlpNumClearVar (azAj);
}

void EGLPNUM_TYPENAME_ILLratio_pivotin_test (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int *rlist,
	int rcnt,
	EGLPNUM_TYPENAME_ratio_res * rs)
{
	int i, k, col;
	EGLPNUM_TYPE *x, *l, *u;
	EGLPNUM_TYPE ay_ij,
		at_i, at_l, at_u, ayi_max, y_ij, t_i, t_l, t_u, t_max, yi_max;
	EGLPNUM_TYPE *pivtol = &(lp->tol->pivot_tol);

	if (rcnt <= 0 || rs == NULL)
		return;
	EGLPNUM_TYPENAME_EGlpNumInitVar (ay_ij);
	EGLPNUM_TYPENAME_EGlpNumInitVar (at_i);
	EGLPNUM_TYPENAME_EGlpNumInitVar (at_l);
	EGLPNUM_TYPENAME_EGlpNumInitVar (at_u);
	EGLPNUM_TYPENAME_EGlpNumInitVar (ayi_max);
	EGLPNUM_TYPENAME_EGlpNumInitVar (t_max);
	EGLPNUM_TYPENAME_EGlpNumInitVar (y_ij);
	EGLPNUM_TYPENAME_EGlpNumInitVar (t_i);
	EGLPNUM_TYPENAME_EGlpNumInitVar (t_l);
	EGLPNUM_TYPENAME_EGlpNumInitVar (t_u);
	EGLPNUM_TYPENAME_EGlpNumInitVar (yi_max);
	rs->boundch = 0;
	rs->lindex = -1;
	EGLPNUM_TYPENAME_EGlpNumZero (rs->tz);
	rs->ratio_stat = RATIO_FAILED;
	rs->lvstat = -1;
	EGLPNUM_TYPENAME_EGlpNumZero (rs->pivotval);
	EGLPNUM_TYPENAME_EGlpNumZero (rs->lbound);

	for (i = 0; i < rcnt; i++)
		lp->iwork[rlist[i]] = 1;

	for (k = 0, EGLPNUM_TYPENAME_EGlpNumCopy (t_max, EGLPNUM_TYPENAME_INFTY); k < lp->yjz.nzcnt; k++)
	{
		EGLPNUM_TYPENAME_EGlpNumCopy (y_ij, lp->yjz.coef[k]);
		if (!EGLPNUM_TYPENAME_EGlpNumIsNeqZero (y_ij, *pivtol))
			continue;

		i = lp->yjz.indx[k];
		if (lp->iwork[lp->baz[i]] == 1)
			continue;
		x = &(lp->xbz[i]);
		col = lp->baz[i];
		l = &(lp->lz[col]);
		u = &(lp->uz[col]);
		EGLPNUM_TYPENAME_EGlpNumCopy (t_u, EGLPNUM_TYPENAME_INFTY);
		EGLPNUM_TYPENAME_EGlpNumCopy (at_u, EGLPNUM_TYPENAME_INFTY);
		EGLPNUM_TYPENAME_EGlpNumCopy (t_l, EGLPNUM_TYPENAME_NINFTY);
		EGLPNUM_TYPENAME_EGlpNumCopy (at_l, EGLPNUM_TYPENAME_INFTY);

		if (EGLPNUM_TYPENAME_EGlpNumIsNeqq (*l, EGLPNUM_TYPENAME_NINFTY))
		{
			EGLPNUM_TYPENAME_EGlpNumCopyDiffRatio (t_l, *x, *l, y_ij);
			EGLPNUM_TYPENAME_EGlpNumCopyAbs (at_l, t_l);
			if (EGLPNUM_TYPENAME_EGlpNumIsLess (at_l, t_max))
				EGLPNUM_TYPENAME_EGlpNumCopy (t_max, at_l);
		}
		if (EGLPNUM_TYPENAME_EGlpNumIsNeqq (*u, EGLPNUM_TYPENAME_INFTY))
		{
			EGLPNUM_TYPENAME_EGlpNumCopyDiffRatio (t_u, *x, *u, y_ij);
			EGLPNUM_TYPENAME_EGlpNumCopyAbs (at_u, t_u);
			if (EGLPNUM_TYPENAME_EGlpNumIsLess (at_u, t_max))
				EGLPNUM_TYPENAME_EGlpNumCopy (t_max, at_u);
		}
	}

	if (EGLPNUM_TYPENAME_EGlpNumIsLeq (EGLPNUM_TYPENAME_INFTY, t_max))
	{
		rs->ratio_stat = RATIO_UNBOUNDED;
		ILL_CLEANUP;
	}

	EGLPNUM_TYPENAME_EGlpNumZero (yi_max);
	EGLPNUM_TYPENAME_EGlpNumZero (ayi_max);
	EGLPNUM_TYPENAME_EGlpNumMultUiTo (t_max, 101);
	EGLPNUM_TYPENAME_EGlpNumDivUiTo (t_max, 100);
	for (k = 0; k < lp->yjz.nzcnt; k++)
	{
		EGLPNUM_TYPENAME_EGlpNumCopy (y_ij, lp->yjz.coef[k]);
		EGLPNUM_TYPENAME_EGlpNumCopyAbs (ay_ij, y_ij);
		if (!EGLPNUM_TYPENAME_EGlpNumIsNeqZero (y_ij, *pivtol))
			continue;

		i = lp->yjz.indx[k];
		if (lp->iwork[lp->baz[i]] == 1)
			continue;
		x = &(lp->xbz[i]);
		col = lp->baz[i];
		l = &(lp->lz[col]);
		u = &(lp->uz[col]);

		EGLPNUM_TYPENAME_EGlpNumCopy (t_u, EGLPNUM_TYPENAME_INFTY);
		EGLPNUM_TYPENAME_EGlpNumCopy (at_u, t_u);
		EGLPNUM_TYPENAME_EGlpNumCopy (t_l, EGLPNUM_TYPENAME_NINFTY);
		EGLPNUM_TYPENAME_EGlpNumCopy (at_l, t_u);
		if (EGLPNUM_TYPENAME_EGlpNumIsNeqq (*l, EGLPNUM_TYPENAME_NINFTY))
		{
			EGLPNUM_TYPENAME_EGlpNumCopyDiffRatio (t_l, *x, *l, y_ij);
			EGLPNUM_TYPENAME_EGlpNumCopyAbs (at_l, t_l);
		}
		if (EGLPNUM_TYPENAME_EGlpNumIsNeqq (*u, EGLPNUM_TYPENAME_INFTY))
		{
			EGLPNUM_TYPENAME_EGlpNumCopyDiffRatio (t_u, *x, *u, y_ij);
			EGLPNUM_TYPENAME_EGlpNumCopyAbs (at_u, t_u);
		}
		//t_i = (fabs (t_l) < fabs (t_u)) ? t_l : t_u;
		if (EGLPNUM_TYPENAME_EGlpNumIsLess (at_l, at_u))
		{
			EGLPNUM_TYPENAME_EGlpNumCopy (t_i, t_l);
			EGLPNUM_TYPENAME_EGlpNumCopy (at_i, at_l);
		}
		else
		{
			EGLPNUM_TYPENAME_EGlpNumCopy (t_i, t_u);
			EGLPNUM_TYPENAME_EGlpNumCopy (at_i, at_u);
		}
		/*if (fabs (t_i) <= t_max + t_max * (1.0e-2)) */
		if (EGLPNUM_TYPENAME_EGlpNumIsLeq (at_i, t_max))
		{
			if (EGLPNUM_TYPENAME_EGlpNumIsLess (ayi_max, ay_ij))
			{
				EGLPNUM_TYPENAME_EGlpNumCopy (yi_max, y_ij);
				EGLPNUM_TYPENAME_EGlpNumCopy (ayi_max, ay_ij);
				rs->lindex = i;
				EGLPNUM_TYPENAME_EGlpNumCopy (rs->tz, t_i);
				rs->lvstat = (EGLPNUM_TYPENAME_EGlpNumIsLess (at_l, at_u)) ? STAT_LOWER : STAT_UPPER;
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
		EGLPNUM_TYPENAME_EGlpNumCopy (rs->pivotval, yi_max);
	}
CLEANUP:
	for (i = 0; i < rcnt; i++)
		lp->iwork[rlist[i]] = 0;
	EGLPNUM_TYPENAME_EGlpNumClearVar (t_max);
	EGLPNUM_TYPENAME_EGlpNumClearVar (ay_ij);
	EGLPNUM_TYPENAME_EGlpNumClearVar (at_i);
	EGLPNUM_TYPENAME_EGlpNumClearVar (at_l);
	EGLPNUM_TYPENAME_EGlpNumClearVar (at_u);
	EGLPNUM_TYPENAME_EGlpNumClearVar (ayi_max);
	EGLPNUM_TYPENAME_EGlpNumClearVar (y_ij);
	EGLPNUM_TYPENAME_EGlpNumClearVar (t_i);
	EGLPNUM_TYPENAME_EGlpNumClearVar (t_l);
	EGLPNUM_TYPENAME_EGlpNumClearVar (t_u);
	EGLPNUM_TYPENAME_EGlpNumClearVar (yi_max);
	return;
}
