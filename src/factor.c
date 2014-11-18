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

/* RCS_INFO = "$RCSfile: factor.c,v $ $Revision: 1.2 $ $Date: 2003/11/05 16:49:52 $"; */
//static int TRACE = 0;

/* implement a = max(a,abs(b)) and execute the extra code if the update is
 * needed */
#define EGlpNumSetToMaxAbsAndDo(a,b,c) \
	if(EGlpNumIsGreatZero(b))\
	{\
		if(EGlpNumIsLess(a,b)){\
			EGlpNumCopy(a,b);\
			c;\
			}\
	}\
	else\
	{\
		EGlpNumSign(a);\
		if(EGlpNumIsLess(b,a)){\
			EGlpNumCopy(a,b);\
			c;\
			}\
		EGlpNumSign(a);\
	}


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "config.h"

#include "eg_lpnum.h"
#include "eg_io.h"

#include "iqsutil.h"
#include "lpdefs.h"
#include "factor.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif


#undef  RECORD
#undef  DEBUG_FACTOR
#undef  SOLVE_DEBUG

#undef  FACTOR_DEBUG
#undef  UPDATE_DEBUG

#undef TRACK_FACTOR
#undef NOTICE_BLOWUP

#undef  FACTOR_STATS
#undef  UPDATE_STATS
#undef  GROWTH_STATS

#undef  UPDATE_STUDY

#undef  SORT_RESULTS

#ifdef UPDATE_STUDY
int nupdate = 0;
long int colspiketot = 0.0;
long int rowspiketot = 0.0;
long int permshifttot = 0.0;
long int leftetatot = 0.0;
#endif

void ILLfactor_init_factor_work (
	factor_work * f)
{
	f->max_k = 1000;							/* must be less than 46340 (2^15.5) */
	EGlpNumCopy (f->fzero_tol, SZERO_TOLER);	/* 2^-50 */
	EGlpNumCopy (f->szero_tol, SZERO_TOLER);	/* 2^-50 */
	EGlpNumCopy (f->partial_tol, OBJBND_TOLER);	/* 2^-7 */
	f->ur_space_mul = 2.0;
	f->uc_space_mul = 1.1;
	f->lc_space_mul = 1.1;
	f->er_space_mul = 1000.0;
	f->grow_mul = 1.5;
	f->p = 4;
	f->etamax = 100;
	f->minmult = 1e3;
	f->maxmult = 1e5;
	f->updmaxmult = 1e7;
	f->dense_fract = 0.25;
	f->dense_min = 25;
	EGlpNumCopy (f->partial_cur, f->partial_tol);
	f->work_coef = 0;
	f->work_indx = 0;
	f->uc_inf = 0;
	f->ur_inf = 0;
	f->lc_inf = 0;
	f->lr_inf = 0;
	f->er_inf = 0;
	f->ucindx = 0;
	f->ucrind = 0;
	f->uccoef = 0;
	f->urindx = 0;
	f->urcind = 0;
	f->urcoef = 0;
	f->lcindx = 0;
	f->lccoef = 0;
	f->lrindx = 0;
	f->lrcoef = 0;
	f->erindx = 0;
	f->ercoef = 0;
	f->rperm = 0;
	f->rrank = 0;
	f->cperm = 0;
	f->crank = 0;
	f->dmat = 0;
	ILLsvector_init (&f->xtmp);
}

void ILLfactor_free_factor_work (
	factor_work * f)
{
#ifdef UPDATE_STUDY
	if (nupdate)
	{
		MESSAGE(0, "UPDATE STUDY: avg %d upd: %.2f col, %.2f row, %.2f lefteta, "
						"%.2f perm", nupdate, ((double) colspiketot) / nupdate, 
						((double) rowspiketot) / nupdate, ((double) leftetatot) / nupdate,
						((double) permshifttot) / nupdate);
	}
#endif
	EGlpNumFreeArray (f->work_coef);
	ILL_IFFREE (f->work_indx, int);

	ILL_IFFREE (f->uc_inf, uc_info);
	if (f->dim + f->max_k > 0 && f->ur_inf)
	{
		unsigned int i = f->dim + f->max_k + 1;

		while (i--)
			EGlpNumClearVar (f->ur_inf[i].max);
	}
	ILL_IFFREE (f->ur_inf, ur_info);
	ILL_IFFREE (f->lc_inf, lc_info);
	ILL_IFFREE (f->lr_inf, lr_info);
	ILL_IFFREE (f->er_inf, er_info);
	ILL_IFFREE (f->ucindx, int);
	ILL_IFFREE (f->ucrind, int);

	EGlpNumFreeArray (f->uccoef);
	ILL_IFFREE (f->urindx, int);
	ILL_IFFREE (f->urcind, int);

	EGlpNumFreeArray (f->urcoef);
	ILL_IFFREE (f->lcindx, int);

	EGlpNumFreeArray (f->lccoef);
	ILL_IFFREE (f->lrindx, int);

	EGlpNumFreeArray (f->lrcoef);
	ILL_IFFREE (f->erindx, int);

	EGlpNumFreeArray (f->ercoef);
	ILL_IFFREE (f->rperm, int);
	ILL_IFFREE (f->rrank, int);
	ILL_IFFREE (f->cperm, int);
	ILL_IFFREE (f->crank, int);

	EGlpNumFreeArray (f->dmat);
	ILLsvector_free (&f->xtmp);
}

int ILLfactor_set_factor_iparam (
	factor_work * f,
	int param,
	int val)
{
	switch (param)
	{
	case QS_FACTOR_MAX_K:
		f->max_k = val;
		break;
	case QS_FACTOR_P:
		f->p = val;
		break;
	case QS_FACTOR_ETAMAX:
		f->etamax = val;
		break;
	case QS_FACTOR_DENSE_MIN:
		f->dense_min = val;
		break;
	default:
		fprintf (stderr, "Invalid param %d in ILLfactor_set_factor_iparam\n",
						 param);
		return 1;
	}
	return 0;
}

int ILLfactor_set_factor_dparam (
	factor_work * f,
	int param,
	EGlpNum_t val)
{
	switch (param)
	{
	case QS_FACTOR_FZERO_TOL:
		EGlpNumCopy (f->fzero_tol, val);
		break;
	case QS_FACTOR_SZERO_TOL:
		EGlpNumCopy (f->szero_tol, val);
		break;
	case QS_FACTOR_UR_SPACE_MUL:
		f->ur_space_mul = EGlpNumToLf (val);
		break;
	case QS_FACTOR_UC_SPACE_MUL:
		f->uc_space_mul = EGlpNumToLf (val);
		break;
	case QS_FACTOR_LC_SPACE_MUL:
		f->lc_space_mul = EGlpNumToLf (val);
		break;
	case QS_FACTOR_LR_SPACE_MUL:
		f->lr_space_mul = EGlpNumToLf (val);
		break;
	case QS_FACTOR_ER_SPACE_MUL:
		f->er_space_mul = EGlpNumToLf (val);
		break;
	case QS_FACTOR_GROW_MUL:
		f->grow_mul = EGlpNumToLf (val);
		break;
	case QS_FACTOR_MAXMULT:
		f->maxmult = EGlpNumToLf (val);
		break;
	case QS_FACTOR_UPDMAXMULT:
		f->updmaxmult = EGlpNumToLf (val);
		break;
	case QS_FACTOR_DENSE_FRACT:
		f->dense_fract = EGlpNumToLf (val);
		break;
	case QS_FACTOR_PARTIAL_TOL:
		EGlpNumCopy (f->partial_tol, val);
		EGlpNumCopy (f->partial_cur, val);
		break;
	default:
		fprintf (stderr, "Invalid param %d in ILLfactor_set_factor_dparam\n",
						 param);
		return 1;
	}
	return 0;
}

int ILLfactor_create_factor_work (
	factor_work * f,
	int dim)
{
	int i;
	int rval;

	f->dim = dim;
	f->etacnt = 0;
	f->work_coef = EGlpNumAllocArray (dim);
	ILL_SAFE_MALLOC (f->work_indx, dim, int);

	ILL_SAFE_MALLOC (f->uc_inf, dim + (f->max_k + 1), uc_info);
	ILL_SAFE_MALLOC (f->ur_inf, dim + (f->max_k + 1), ur_info);
	ILL_SAFE_MALLOC (f->lc_inf, dim, lc_info);
	ILL_SAFE_MALLOC (f->lr_inf, dim, lr_info);
	ILL_SAFE_MALLOC (f->rperm, dim, int);
	ILL_SAFE_MALLOC (f->rrank, dim, int);
	ILL_SAFE_MALLOC (f->cperm, dim, int);
	ILL_SAFE_MALLOC (f->crank, dim, int);

	for (i = dim + f->max_k + 1; i--;)
		EGlpNumInitVar (f->ur_inf[i].max);

	for (i = 0; i < dim; i++)
	{
		EGlpNumZero (f->work_coef[i]);
		f->work_indx[i] = 0;
		f->uc_inf[i].nzcnt = 0;
		f->ur_inf[i].nzcnt = 0;
		f->lc_inf[i].nzcnt = 0;
		f->lr_inf[i].nzcnt = 0;
		f->rperm[i] = i;
		f->rrank[i] = i;
		f->cperm[i] = i;
		f->crank[i] = i;
	}
	for (i = 0; i <= f->max_k; i++)
	{
		f->uc_inf[dim + i].nzcnt = i;
		f->uc_inf[dim + i].next = dim + i;
		f->uc_inf[dim + i].prev = dim + i;
		f->ur_inf[dim + i].nzcnt = i;
		f->ur_inf[dim + i].next = dim + i;
		f->ur_inf[dim + i].prev = dim + i;
	}

	rval = ILLsvector_alloc (&f->xtmp, dim);
	CHECKRVALG (rval, CLEANUP);

	rval = 0;

CLEANUP:
	if (rval)
	{
		ILLfactor_free_factor_work (f);
	}
	EG_RETURN (rval);
}
#ifdef FACTOR_DEBUG
static void dump_matrix (
	factor_work * f,
	int remaining)
{
	int dim = f->dim;
	ur_info *ur_inf = f->ur_inf;
	uc_info *uc_inf = f->uc_inf;
	lc_info *lc_inf = f->lc_inf;
	lr_info *lr_inf = f->lr_inf;
	er_info *er_inf = f->er_inf;
	int nzcnt;
	int beg;

	int i;
	int j;

	for (i = 0; i < dim; i++)
	{
		if (!remaining || ur_inf[i].next >= 0)
		{
			printf ("Row %d %d (max %.3f):", i, f->rrank[i],
							EGlpNumToLf (ur_inf[i].max));
			nzcnt = ur_inf[i].nzcnt;
			beg = ur_inf[i].rbeg;
			for (j = 0; j < nzcnt; j++)
			{
				if (j == ur_inf[i].pivcnt)
				{
					printf (" |");
				}
				printf (" %.3f*%d", EGlpNumToLf (f->urcoef[beg + j]),
								f->urindx[beg + j]);
				if (f->urcind)
					printf ("@%d", f->urcind[beg + j]);
			}
			printf ("\n");
		}
	}
	if (f->dmat)
	{
		int start = 0;

		if (remaining)
			start = f->stage - f->dense_base;
		printf ("Dcols at %d %d - %d    :", f->stage - f->dense_base,
						f->dense_base + start, f->nstages);
		for (j = start; j < f->dcols; j++)
		{
			printf (" %5d", f->cperm[j + f->dense_base]);
		}
		printf ("\n");
		for (i = start; i < f->drows; i++)
		{
			printf ("DRow %d %d (max %.3f):", i,
							f->rperm[i + f->dense_base],
							EGlpNumToLf (ur_inf[f->rperm[i + f->dense_base]].max));
			for (j = start; j < f->dcols; j++)
			{
				if (j == f->drows)
				{
					printf (" |");
				}
				printf (" %.3f", EGlpNumToLf (f->dmat[i * f->dcols + j]));
			}
			printf ("\n");
		}
	}

	if (!remaining)
	{
		for (i = 0; i < f->stage; i++)
		{
			printf ("L col %d:", lc_inf[i].c);
			nzcnt = lc_inf[i].nzcnt;
			beg = lc_inf[i].cbeg;
			for (j = 0; j < nzcnt; j++)
			{
				printf (" %.3f*%d", EGlpNumToLf (f->lccoef[beg + j]),
								f->lcindx[beg + j]);
			}
			printf ("\n");
		}
		for (i = f->nstages; i < f->dim; i++)
		{
			printf ("L col %d:", lc_inf[i].c);
			nzcnt = lc_inf[i].nzcnt;
			beg = lc_inf[i].cbeg;
			for (j = 0; j < nzcnt; j++)
			{
				printf (" %.3f*%d", EGlpNumToLf (f->lccoef[beg + j]),
								f->lcindx[beg + j]);
			}
			printf ("\n");
		}
		for (i = 0; i < f->dim; i++)
		{
			if (!lr_inf[i].nzcnt)
				continue;
			printf ("L row %d:", lr_inf[i].r);
			nzcnt = lr_inf[i].nzcnt;
			beg = lr_inf[i].rbeg;
			for (j = 0; j < nzcnt; j++)
			{
				printf (" %.3f*%d", EGlpNumToLf (f->lrcoef[beg + j]),
								f->lrindx[beg + j]);
			}
			printf ("\n");
		}
	}

	if (!remaining)
	{
		for (i = 0; i < f->etacnt; i++)
		{
			printf ("Eta row %d:", f->er_inf[i].r);
			nzcnt = er_inf[i].nzcnt;
			beg = er_inf[i].rbeg;
			for (j = 0; j < nzcnt; j++)
			{
				printf (" %.3f*%d", EGlpNumToLf (f->ercoef[beg + j]),
								f->erindx[beg + j]);
			}
			printf ("\n");
		}
	}

	for (i = 0; i < dim; i++)
	{
		if (!remaining || uc_inf[i].next >= 0)
		{
			printf ("Col %d %d:", i, f->crank[i]);
			nzcnt = uc_inf[i].nzcnt;
			beg = uc_inf[i].cbeg;
			for (j = 0; j < nzcnt; j++)
			{
				if (f->uccoef != 0)
				{
					printf (" %.3f*%d", EGlpNumToLf (f->uccoef[beg + j]),
									f->ucindx[beg + j]);
					if (f->ucrind)
						printf ("@%d", f->ucrind[beg + j]);
				}
				else
				{
					printf (" %d", f->ucindx[beg + j]);
				}
			}
			printf ("\n");
		}
	}

	if (!remaining)
	{
		printf ("rperm:");
		for (i = 0; i < dim; i++)
		{
			if (i == f->nstages)
				printf ("|");
			if (i == f->stage)
				printf ("|");
			printf (" %d", f->rperm[i]);
		}
		printf ("\n");

		printf ("cperm:");
		for (i = 0; i < dim; i++)
		{
			if (i == f->nstages)
				printf ("|");
			if (i == f->stage)
				printf ("|");
			printf (" %d", f->cperm[i]);
		}
		printf ("\n");
	}

	printf ("Rows by nzcnt:\n");
	for (i = 0; i <= f->max_k; i++)
	{
		if (ur_inf[dim + i].next != dim + i)
		{
			printf ("%d:", i);
			for (j = ur_inf[dim + i].next; j != dim + i; j = ur_inf[j].next)
			{
				printf (" %d", j);
			}
			printf ("\n");
		}
	}

	printf ("Cols by nzcnt:\n");
	for (i = 0; i <= f->max_k; i++)
	{
		if (uc_inf[dim + i].next != dim + i)
		{
			printf ("%d:", i);
			for (j = uc_inf[dim + i].next; j != dim + i; j = uc_inf[j].next)
			{
				printf (" %d", j);
			}
			printf ("\n");
		}
	}

	printf ("\n");
	fflush (stdout);
}
#endif

#ifdef SORT_RESULTS
static void sort_vector2 (
	int nzcnt,
	int *indx,
	EGlpNum_t * coef)
{
	int i;
	int j;
	int itmp;
	EGlpNum_t ctmp;

	EGlpNumInitVar (ctmp);

	for (i = 1; i < nzcnt; i++)
	{
		itmp = indx[i];
		EGlpNumCopy (ctmp, coef[i]);
		for (j = i; j >= 1 && indx[j - 1] > itmp; j--)
		{
			indx[j] = indx[j - 1];
			EGlpNumCopy (coef[j], coef[j - 1]);
		}
		indx[j] = itmp;
		EGlpNumCopy (coef[j], ctmp);
	}
	EGlpNumClearVar (ctmp);
}

static void sort_vector (
	svector * x)
{
	sort_vector2 (x->nzcnt, x->indx, x->coef);
}
#endif

#ifdef DEBUG_FACTOR
static int check_matrix (
	factor_work * f)
{
	ur_info *ur_inf = f->ur_inf;
	uc_info *uc_inf = f->uc_inf;
	int rbeg;
	int nzcnt;
	int cbeg;
	int c;
	int r;
	int j;
	int nerr = 0;

	for (r = 0; r < f->dim; r++)
	{
		nzcnt = ur_inf[r].nzcnt;
		rbeg = ur_inf[r].rbeg;
		for (j = 0; j < nzcnt; j++)
		{
			c = f->urindx[rbeg + j];
			cbeg = uc_inf[c].cbeg;
			if (f->ucindx[cbeg + f->urcind[rbeg + j]] != r)
			{
				MESSAGE(0,"index mismatch, row %d column %d", r, c);
				nerr++;
			}
			if (fabs(EGlpNumToLf(f->uccoef[cbeg + f->urcind[rbeg + j]]) - EGlpNumToLf(f->urcoef[rbeg + j]))>1000*EGlpNumToLf(epsLpNum))
			{
				MESSAGE(0,"coef mismatch, row %d column %d", r, c);
				nerr++;
			}
		}
	}
	if (f->urindx[f->ur_space] != 0)
	{
		MESSAGE(0,"last urindx entry %d != 0", f->urindx[f->ur_space]);
		nerr++;
	}

	for (c = 0; c < f->dim; c++)
	{
		nzcnt = uc_inf[c].nzcnt;
		cbeg = uc_inf[c].cbeg;
		for (j = 0; j < nzcnt; j++)
		{
			r = f->ucindx[cbeg + j];
			rbeg = ur_inf[r].rbeg;
			if (f->urindx[rbeg + f->ucrind[cbeg + j]] != c)
			{
				MESSAGE(0,"index mismatch, column %d row %d", c, r);
				nerr++;
			}
			if (f->urcoef[rbeg + f->ucrind[cbeg + j]] != f->uccoef[cbeg + j])
			{
				MESSAGE(0,"coef mismatch, column %d row %d", c, r);
				nerr++;
			}
		}
	}
	if (f->ucindx[f->uc_space] != 0)
	{
		MESSAGE(0,"last ucindx entry %d != 0", f->ucindx[f->uc_space]);
		nerr++;
	}
	if (nerr)
	{
		dump_matrix (f, 0);
		return E_CHECK_FAILED;
	}
	return 0;
}
#endif

#ifdef FACTOR_STATS
static void dump_factor_stats (
	factor_work * f)
{
	int dim = f->dim;
	int ecnt = f->etacnt;
	ur_info *ur_inf = f->ur_inf;
	lc_info *lc_inf = f->lc_inf;
	er_info *er_inf = f->er_inf;
	EGlpNum_t *urcoef = f->urcoef;
	EGlpNum_t *lccoef = f->lccoef;
	EGlpNum_t *ercoef = f->ercoef;
	int lnzcnt = 0;
	int unzcnt = 0;
	int enzcnt = 0;
	int nzcnt;
	int beg;
	EGlpNum_t umax;
	EGlpNum_t lmax;
	EGlpNum_t emax;
	int i;
	int j;

	EGlpNumInitVar (umax);
	EGlpNumInitVar (lmax);
	EGlpNumInitVar (emax);
	EGlpNumZero (umax);
	for (i = 0; i < dim; i++)
	{
		nzcnt = ur_inf[i].nzcnt;
		beg = ur_inf[i].rbeg;
		unzcnt += nzcnt;
		for (j = 0; j < nzcnt; j++)
		{
			EGlpNumSetToMaxAbs (umax, urcoef[beg + j]);
		}
	}
	EGlpNumZero (lmax);
	for (i = 0; i < dim; i++)
	{
		nzcnt = lc_inf[i].nzcnt;
		beg = lc_inf[i].cbeg;
		lnzcnt += nzcnt;
		for (j = 0; j < nzcnt; j++)
		{
			EGlpNumSetToMaxAbs (lmax, lccoef[beg + j]);
		}
	}
	EGlpNumZero (emax);
	for (i = 0; i < ecnt; i++)
	{
		nzcnt = er_inf[i].nzcnt;
		beg = er_inf[i].rbeg;
		enzcnt += nzcnt;
		for (j = 0; j < nzcnt; j++)
		{
			EGlpNumSetToMaxAbs (emax, ercoef[beg + j]);
		}
	}
	MESSAGE(0, "factor U %d nzs %.3e max L %d nzs %.3e max E %d nzs %.3e max",
					unzcnt, EGlpNumToLf (umax), lnzcnt, EGlpNumToLf (lmax), enzcnt,
					EGlpNumToLf (emax));
	fflush (stdout);
	EGlpNumClearVar (umax);
	EGlpNumClearVar (lmax);
	EGlpNumClearVar (emax);
}
#endif

static void clear_work (
	factor_work * f)
{
	int i;
	int dim = f->dim;
	EGlpNum_t *work_coef = f->work_coef;

	for (i = 0; i < dim; i++)
	{
		EGlpNumZero (work_coef[i]);
	}
}

static void load_row (
	factor_work * f,
	int r)
{
	EGlpNum_t *prow_urcoef = f->urcoef + f->ur_inf[r].rbeg;
	int *prow_urindx = f->urindx + f->ur_inf[r].rbeg;
	int prow_nzcnt = f->ur_inf[r].nzcnt;
	EGlpNum_t *work_coef = f->work_coef;
	int *work_indx = f->work_indx;
	int i;
	int j;

	for (i = 0; i < prow_nzcnt; i++)
	{
		j = prow_urindx[i];
		EGlpNumCopy (work_coef[j], prow_urcoef[i]);
		work_indx[j] = 1;
	}
}

static void clear_row (
	factor_work * f,
	int r)
{
	int *prow_urindx = f->urindx + f->ur_inf[r].rbeg;
	int prow_nzcnt = f->ur_inf[r].nzcnt;
	EGlpNum_t *work_coef = f->work_coef;
	int *work_indx = f->work_indx;
	int i;
	int j;

	for (i = 0; i < prow_nzcnt; i++)
	{
		j = prow_urindx[i];
		EGlpNumZero (work_coef[j]);
		work_indx[j] = 0;
	}
}

static int make_ur_space (
	factor_work * f,
	int space)
{
	EGlpNum_t *new_urcoef = 0;
	int *new_urindx = 0;
	int *new_urcind = 0;
	EGlpNum_t *urcoef = f->urcoef;
	int *urindx = f->urindx;
	int *urcind = f->urcind;
	int minspace;
	ur_info *ur_inf = f->ur_inf;
	int dim = f->dim;
	int new_nzcnt = 0, old_nzcnt;
	int rbeg;
	int nzcnt;
	int i;
	int j;
	int rval;

	minspace = f->ur_space;
	nzcnt = space;
	for (i = 0; i < dim; i++)
		nzcnt += ur_inf[i].nzcnt;
	old_nzcnt = nzcnt;
	while (nzcnt * 2 >= minspace)
	{
		minspace = 1 + minspace * f->grow_mul;
	}

#ifdef GROWTH_STATS
	printf ("make_ur_space growing from %d to %d...", f->ur_space, minspace);
	fflush (stdout);
#endif
	new_urcoef = EGlpNumAllocArray (minspace);
	ILL_SAFE_MALLOC (new_urindx, minspace + 1, int);

	if (urcind)
	{
		ILL_SAFE_MALLOC (new_urcind, minspace, int);
	}

	if (urcind)
	{
		for (j = 0; j < dim; j++)
		{
			rbeg = ur_inf[j].rbeg;
			nzcnt = ur_inf[j].nzcnt;
			ur_inf[j].rbeg = new_nzcnt;
			for (i = 0; i < nzcnt; i++)
			{
				new_urindx[new_nzcnt] = urindx[rbeg + i];
				EGlpNumCopy (new_urcoef[new_nzcnt], urcoef[rbeg + i]);
				new_urcind[new_nzcnt] = urcind[rbeg + i];
				new_nzcnt++;
			}
		}
	}
	else
	{
		for (j = 0; j < dim; j++)
		{
			rbeg = ur_inf[j].rbeg;
			nzcnt = ur_inf[j].nzcnt;
			ur_inf[j].rbeg = new_nzcnt;
			for (i = 0; i < nzcnt; i++)
			{
				new_urindx[new_nzcnt] = urindx[rbeg + i];
				EGlpNumCopy (new_urcoef[new_nzcnt], urcoef[rbeg + i]);
				new_nzcnt++;
			}
		}
	}

	for (i = new_nzcnt; i < minspace; i++)
	{
		new_urindx[i] = -1;
	}
	new_urindx[minspace] = 0;
	EGlpNumFreeArray (f->urcoef);
	f->urcoef = new_urcoef;
	new_urcoef = 0;

	ILL_IFFREE (f->urindx, int);

	f->urindx = new_urindx;
	new_urindx = 0;

	ILL_IFFREE (f->urcind, int);

	f->urcind = new_urcind;
	new_urcind = 0;

	f->ur_freebeg = new_nzcnt;
	f->ur_space = minspace;

#ifdef GROWTH_STATS
	MESSAGE (0,"%d/%d nonzeros", new_nzcnt, old_nzcnt);
	dump_factor_stats (f);
#endif

	rval = 0;

CLEANUP:
	ILL_IFFREE (new_urcoef, EGlpNum_t);
	ILL_IFFREE (new_urindx, int);
	ILL_IFFREE (new_urcind, int);

	EG_RETURN (rval);
}

static int make_uc_space (
	factor_work * f,
	int space)
{
	EGlpNum_t *new_uccoef = 0;
	int *new_ucindx = 0;
	int *new_ucrind = 0;
	int uc_freebeg = f->uc_freebeg;
	EGlpNum_t *uccoef = f->uccoef;
	int *ucindx = f->ucindx;
	int *ucrind = f->ucrind;
	int minspace = uc_freebeg + space;
	uc_info *uc_inf = f->uc_inf;
	int dim = f->dim;
	int new_nzcnt = 0;
	int cbeg;
	int nzcnt;
	int i;
	int j;
	int rval;

	minspace = f->uc_space;
	nzcnt = space;
	for( i = 0 ; i < dim ; i++) nzcnt += uc_inf[i].nzcnt;
	while(nzcnt*2 >= minspace)
	{
		minspace = 10 + (f->grow_mul * minspace);
	}

#ifdef GROWTH_STATS
	MESSAGE (0,"make_uc_space growing from %d to %d...", f->uc_space, minspace);
#endif

	ILL_SAFE_MALLOC (new_ucindx, minspace + 1, int);

	if (ucrind)
	{
		new_uccoef = EGlpNumAllocArray (minspace);
		ILL_SAFE_MALLOC (new_ucrind, minspace, int);
	}

	if (ucrind)
	{
		for (j = 0; j < dim; j++)
		{
			cbeg = uc_inf[j].cbeg;
			nzcnt = uc_inf[j].nzcnt;
			uc_inf[j].cbeg = new_nzcnt;
			for (i = 0; i < nzcnt; i++)
			{
				new_ucindx[new_nzcnt] = ucindx[cbeg + i];
				EGlpNumCopy (new_uccoef[new_nzcnt], uccoef[cbeg + i]);
				new_ucrind[new_nzcnt] = ucrind[cbeg + i];
				new_nzcnt++;
			}
		}
	}
	else
	{
		for (j = 0; j < dim; j++)
		{
			cbeg = uc_inf[j].cbeg;
			nzcnt = uc_inf[j].nzcnt;
			uc_inf[j].cbeg = new_nzcnt;
			for (i = 0; i < nzcnt; i++)
			{
				new_ucindx[new_nzcnt] = ucindx[cbeg + i];
				new_nzcnt++;
			}
		}
	}

	for (i = new_nzcnt; i < minspace; i++)
	{
		new_ucindx[i] = -1;
	}
	new_ucindx[minspace] = 0;

	EGlpNumFreeArray (f->uccoef);
	f->uccoef = new_uccoef;
	new_uccoef = 0;

	ILL_IFFREE (f->ucindx, int);

	f->ucindx = new_ucindx;
	new_ucindx = 0;

	ILL_IFFREE (f->ucrind, int);

	f->ucrind = new_ucrind;
	new_ucrind = 0;

	f->uc_freebeg = new_nzcnt;
	f->uc_space = minspace;

#ifdef GROWTH_STATS
	MESSAGE (0,"%d nonzeros", new_nzcnt);
	dump_factor_stats (f);
#endif

	rval = 0;

CLEANUP:
	ILL_IFFREE (new_uccoef, EGlpNum_t);
	ILL_IFFREE (new_ucindx, int);
	ILL_IFFREE (new_ucrind, int);

	EG_RETURN (rval);
}

static int make_lc_space (
	factor_work * f,
	int space)
{
	EGlpNum_t *new_lccoef = 0;
	int *new_lcindx = 0;
	int lc_freebeg = f->lc_freebeg;
	EGlpNum_t *lccoef = f->lccoef;
	int *lcindx = f->lcindx;
	int minspace = lc_freebeg + space;
	int i;
	int rval;

	if (f->lc_space * f->grow_mul > minspace)
	{
		minspace = f->lc_space * f->grow_mul;
	}

#ifdef GROWTH_STATS
	MESSAGE (0,"make_lc_space growing from %d to %d...", f->lc_space, minspace);
#endif

	new_lccoef = EGlpNumAllocArray (minspace);
	ILL_SAFE_MALLOC (new_lcindx, minspace, int);

	for (i = 0; i < lc_freebeg; i++)
	{
		EGlpNumCopy (new_lccoef[i], lccoef[i]);
		new_lcindx[i] = lcindx[i];
	}

	EGlpNumFreeArray (lccoef);
	f->lccoef = new_lccoef;
	new_lccoef = 0;

	ILL_IFFREE (lcindx, int);

	f->lcindx = new_lcindx;
	new_lcindx = 0;

	f->lc_space = minspace;

#ifdef GROWTH_STATS
	dump_factor_stats (f);
#endif

	rval = 0;

CLEANUP:
	ILL_IFFREE (new_lccoef, EGlpNum_t);
	ILL_IFFREE (new_lcindx, int);

	EG_RETURN (rval);
}

static void set_col_nz (
	factor_work * f,
	int c)
{
	uc_info *uc_inf = f->uc_inf;
	int nzcnt = uc_inf[c].nzcnt;
	int max_k = f->max_k;
	int dim = f->dim;

	if (uc_inf[c].next >= 0)
	{
		uc_inf[uc_inf[c].next].prev = uc_inf[c].prev;
		uc_inf[uc_inf[c].prev].next = uc_inf[c].next;

		if (nzcnt >= max_k)
			nzcnt = max_k;
		uc_inf[c].next = uc_inf[dim + nzcnt].next;
		uc_inf[c].prev = dim + nzcnt;
		uc_inf[dim + nzcnt].next = c;
		uc_inf[uc_inf[c].next].prev = c;
	}
}

static void set_row_nz (
	factor_work * f,
	int r)
{
	ur_info *ur_inf = f->ur_inf;
	int nzcnt = ur_inf[r].pivcnt;
	int max_k = f->max_k;
	int dim = f->dim;

	if (ur_inf[r].next >= 0)
	{
		ur_inf[ur_inf[r].next].prev = ur_inf[r].prev;
		ur_inf[ur_inf[r].prev].next = ur_inf[r].next;

		if (nzcnt >= max_k)
			nzcnt = max_k;
		ur_inf[r].next = ur_inf[dim + nzcnt].next;
		ur_inf[r].prev = dim + nzcnt;
		ur_inf[dim + nzcnt].next = r;
		ur_inf[ur_inf[r].next].prev = r;
	}
}

static void remove_col_nz (
	factor_work * f,
	int r,
	int c)
{
	uc_info *uc_inf = f->uc_inf;
	int *ucindx = f->ucindx + uc_inf[c].cbeg;
	int nzcnt = uc_inf[c].nzcnt;
	int i;

	for (i = 0; i < nzcnt; i++)
	{
		if (ucindx[i] == r)
		{
			--nzcnt;
			ucindx[i] = ucindx[nzcnt];
			ucindx[nzcnt] = -1;
			break;
		}
	}
	uc_inf[c].nzcnt = nzcnt;

	set_col_nz (f, c);
}

static void remove_row_nz (
	factor_work * f,
	int r,
	int c)
{
	ur_info *ur_inf = f->ur_inf;
	int *urindx = f->urindx + ur_inf[r].rbeg;
	EGlpNum_t *urcoef = f->urcoef + ur_inf[r].rbeg;
	int pivcnt = ur_inf[r].pivcnt;
	EGlpNum_t max;
	int tind;
	EGlpNum_t tcoef;
	int i;

	EGlpNumInitVar (tcoef);
	EGlpNumInitVar (max);
	EGlpNumZero (max);

	for (i = 0; i < pivcnt; i++)
	{
		if (urindx[i] == c)
		{
			--pivcnt;
			ILL_SWAP (urindx[i], urindx[pivcnt], tind);
			EGLPNUM_SWAP (urcoef[i], urcoef[pivcnt], tcoef);
			--i;
		}
		else
		{
			EGlpNumSetToMaxAbs (max, urcoef[i]);
		}
	}
	ur_inf[r].pivcnt = pivcnt;
	EGlpNumCopy (ur_inf[r].max, max);
	set_row_nz (f, r);
	EGlpNumClearVar (max);
	EGlpNumClearVar (tcoef);
}

static int add_col_nz (
	factor_work * f,
	int r,
	int c)
{
	uc_info *uc_inf = f->uc_inf;
	int cbeg = uc_inf[c].cbeg;
	int nzcnt = uc_inf[c].nzcnt;
	int uc_freebeg = f->uc_freebeg;
	int *ucindx = f->ucindx;
	int i;
	int rval = 0;

	if (uc_inf[c].next == -1)
	{
		return 0;
	}

	if (ucindx[cbeg + nzcnt] == -1)
	{
		ucindx[cbeg + nzcnt] = r;
		uc_inf[c].nzcnt++;
		if (nzcnt + cbeg == uc_freebeg)
		{
			f->uc_freebeg = uc_freebeg + 1;
		}
	}
	else
	{
		if (uc_freebeg + nzcnt + 1 >= f->uc_space)
		{
			rval = make_uc_space (f, nzcnt + 1);
			CHECKRVALG (rval, CLEANUP);
			uc_freebeg = f->uc_freebeg;
			cbeg = uc_inf[c].cbeg;
			ucindx = f->ucindx;
		}
		for (i = 0; i < nzcnt; i++)
		{
			ucindx[uc_freebeg + i] = ucindx[cbeg + i];
			ucindx[cbeg + i] = -1;
		}
		ucindx[uc_freebeg + nzcnt] = r;
		uc_inf[c].cbeg = uc_freebeg;
		uc_inf[c].nzcnt++;
		f->uc_freebeg = uc_freebeg + nzcnt + 1;
	}

	set_col_nz (f, c);
CLEANUP:
	EG_RETURN (rval);
}

static void disable_col (
	factor_work * f,
	int c)
{
	uc_info *uc_inf = f->uc_inf;

	if (uc_inf[c].next >= 0)
	{
		uc_inf[uc_inf[c].next].prev = uc_inf[c].prev;
		uc_inf[uc_inf[c].prev].next = uc_inf[c].next;

		uc_inf[c].next = -2;
		uc_inf[c].prev = -2;
	}
}

static void remove_col (
	factor_work * f,
	int c)
{
	uc_info *uc_inf = f->uc_inf;
	int cbeg = uc_inf[c].cbeg;
	int nzcnt = uc_inf[c].nzcnt;
	int *ucindx = f->ucindx;
	int i;

	for (i = 0; i < nzcnt; i++)
	{
		ucindx[cbeg + i] = -1;
	}
	uc_inf[c].cbeg = 0;
	uc_inf[c].nzcnt = 0;

	if (uc_inf[c].next >= 0)
	{
		uc_inf[uc_inf[c].next].prev = uc_inf[c].prev;
		uc_inf[uc_inf[c].prev].next = uc_inf[c].next;

		uc_inf[c].next = -1;
		uc_inf[c].prev = -1;
	}
}

static void remove_row (
	factor_work * f,
	int r)
{
	ur_info *ur_inf = f->ur_inf;

	if (ur_inf[r].next >= 0)
	{
		ur_inf[ur_inf[r].next].prev = ur_inf[r].prev;
		ur_inf[ur_inf[r].prev].next = ur_inf[r].next;

		ur_inf[r].next = -1;
		ur_inf[r].prev = -1;
	}
}

static void find_coef (
	factor_work * f,
	int r,
	int c,
	EGlpNum_t * coef)
{
	EGlpNum_t *prow_urcoef = f->urcoef + f->ur_inf[r].rbeg;
	int *prow_urindx = f->urindx + f->ur_inf[r].rbeg;
	int i;
	int prow_nzcnt = f->ur_inf[r].nzcnt;

	EGlpNumZero (*coef);
	for (i = 0; i < prow_nzcnt; i++)
	{
		if (prow_urindx[i] == c)
		{
			EGlpNumCopy (*coef, prow_urcoef[i]);
			return;
		}
	}
	fprintf (stderr, "Coefficient not found\n");
	return;
}

static int elim_row (
	factor_work * f,
	int elim_r,
	int r,
	int c,
	EGlpNum_t * p_pivot_coef)
{
	ur_info *ur_inf = f->ur_inf;
	EGlpNum_t *work_coef = f->work_coef;
	int *work_indx = f->work_indx;
	EGlpNum_t *urcoef = f->urcoef;
	int *urindx = f->urindx;
	int prow_beg = ur_inf[r].rbeg;
	int prow_nzcnt = ur_inf[r].nzcnt;
	int prow_pivcnt = ur_inf[r].pivcnt;
	int fill = ur_inf[elim_r].nzcnt;
	int cancel = 0;
	EGlpNum_t max;
	int erow_beg;
	int erow_nzcnt;
	int erow_pivcnt;
	EGlpNum_t x;
	int i;
	int j;
	int rval = 0;
	EGlpNum_t elim_coef;

	EGlpNumInitVar (max);
	EGlpNumInitVar (x);
	EGlpNumInitVar (elim_coef);
	EGlpNumZero (max);
	find_coef (f, r, c, &elim_coef);
	EGlpNumDivTo (elim_coef, work_coef[c]);
	EGlpNumCopy (*p_pivot_coef, elim_coef);

	for (i = 0; i < prow_nzcnt; i++)
	{
		j = urindx[prow_beg + i];
		if (work_indx[j] == 1)
		{
			EGlpNumCopy (x, urcoef[prow_beg + i]);
			EGlpNumSubInnProdTo (x, elim_coef, work_coef[j]);
			if ((!(EGlpNumIsNeqZero (x, f->fzero_tol))) || j == c)
			{
				cancel++;
				if (j != c)
				{
					remove_col_nz (f, r, j);
				}
				if (i < prow_pivcnt)
				{
					prow_pivcnt--;
					prow_nzcnt--;
					urindx[prow_beg + i] = urindx[prow_beg + prow_pivcnt];
					EGlpNumCopy (urcoef[prow_beg + i], urcoef[prow_beg + prow_pivcnt]);
					if (prow_pivcnt != prow_nzcnt)
					{
						urindx[prow_beg + prow_pivcnt] = urindx[prow_beg + prow_nzcnt];
						EGlpNumCopy (urcoef[prow_beg + prow_pivcnt],
												 urcoef[prow_beg + prow_nzcnt]);
					}
				}
				else
				{
					prow_nzcnt--;
					urindx[prow_beg + i] = urindx[prow_beg + prow_nzcnt];
					EGlpNumCopy (urcoef[prow_beg + i], urcoef[prow_beg + prow_nzcnt]);
				}
				urindx[prow_beg + prow_nzcnt] = -1;
				i--;
			}
			else
			{
				EGlpNumCopy (urcoef[prow_beg + i], x);
				if (i < prow_pivcnt)
				{
					EGlpNumSetToMaxAbs (max, x);
				}
			}
			work_indx[j] = 0;
			fill--;
		}
		else
		{
			if (i < prow_pivcnt)
			{
				EGlpNumSetToMaxAbs (max, urcoef[prow_beg + i]);
			}
		}
	}

	if (fill > 0)
	{
		ur_inf[r].nzcnt = prow_nzcnt;
		ur_inf[r].pivcnt = prow_pivcnt;
		if (fill > cancel)
		{
			int ur_freebeg = f->ur_freebeg;

			if (ur_freebeg + prow_nzcnt + fill >= f->ur_space)
			{
				rval = make_ur_space (f, prow_nzcnt + fill);
				CHECKRVALG (rval, CLEANUP);
				urcoef = f->urcoef;
				urindx = f->urindx;
				ur_freebeg = f->ur_freebeg;
				prow_beg = f->ur_inf[r].rbeg;
			}
			for (i = 0; i < prow_nzcnt; i++)
			{
				urindx[ur_freebeg + i] = urindx[prow_beg + i];
				EGlpNumCopy (urcoef[ur_freebeg + i], urcoef[prow_beg + i]);
				urindx[prow_beg + i] = -1;
			}
			ur_inf[r].rbeg = ur_freebeg;
			f->ur_freebeg = ur_freebeg + prow_nzcnt + fill;
			prow_beg = ur_freebeg;
		}

		erow_beg = ur_inf[elim_r].rbeg;
		erow_nzcnt = ur_inf[elim_r].nzcnt;
		erow_pivcnt = ur_inf[elim_r].pivcnt;

		for (i = 0; i < erow_pivcnt; i++)
		{
			j = urindx[erow_beg + i];
			if (work_indx[j] == 1)
			{
				EGlpNumCopyNeg (x, elim_coef);
				EGlpNumMultTo (x, urcoef[erow_beg + i]);
				if (EGlpNumIsNeqZero (x, f->fzero_tol))
				{
					rval = add_col_nz (f, r, j);
					CHECKRVALG (rval, CLEANUP);
					if (prow_pivcnt != prow_nzcnt)
					{
						urindx[prow_beg + prow_nzcnt] = urindx[prow_beg + prow_pivcnt];
						EGlpNumCopy (urcoef[prow_beg + prow_nzcnt],
												 urcoef[prow_beg + prow_pivcnt]);
					}
					urindx[prow_beg + prow_pivcnt] = j;
					EGlpNumCopy (urcoef[prow_beg + prow_pivcnt], x);
					EGlpNumSetToMaxAbs (max, x);
					prow_pivcnt++;
					prow_nzcnt++;
				}
			}
			else
			{
				work_indx[j] = 1;
			}
		}
		for (i = erow_pivcnt; i < erow_nzcnt; i++)
		{
			j = urindx[erow_beg + i];
			if (work_indx[j] == 1)
			{
				EGlpNumCopyNeg (x, elim_coef);
				EGlpNumMultTo (x, urcoef[erow_beg + i]);
				if (EGlpNumIsNeqZero (x, f->fzero_tol))
				{
					rval = add_col_nz (f, r, j);
					CHECKRVALG (rval, CLEANUP);
					urindx[prow_beg + prow_nzcnt] = j;
					EGlpNumCopy (urcoef[prow_beg + prow_nzcnt], x);
					prow_nzcnt++;
				}
			}
			else
			{
				work_indx[j] = 1;
			}
		}
	}
	else
	{
		erow_nzcnt = ur_inf[elim_r].nzcnt;
		erow_beg = ur_inf[elim_r].rbeg;
		for (i = 0; i < erow_nzcnt; i++)
		{
			j = urindx[erow_beg + i];
			work_indx[j] = 1;
		}
	}

	ur_inf[r].nzcnt = prow_nzcnt;
	ur_inf[r].pivcnt = prow_pivcnt;
	EGlpNumCopy (ur_inf[r].max, max);

	set_row_nz (f, r);
CLEANUP:
	EGlpNumClearVar (elim_coef);
	EGlpNumClearVar (x);
	EGlpNumClearVar (max);
	EG_RETURN (rval);
}

#define SETPERM(f,s,r,c) {                    \
        f->rperm[f->rrank[r]] = f->rperm[s];  \
        f->rrank[f->rperm[s]] = f->rrank[r];  \
        f->rperm[s] = r;                      \
        f->rrank[r] = s;                      \
                                              \
        f->cperm[f->crank[c]] = f->cperm[s];  \
        f->crank[f->cperm[s]] = f->crank[c];  \
        f->cperm[s] = c;                      \
        f->crank[c] = s;                      \
}

static int elim (
	factor_work * f,
	int r,
	int c)
{
	uc_info *uc_inf = f->uc_inf;
	ur_info *ur_inf = f->ur_inf;
	lc_info *lc_inf = f->lc_inf;
	int *urindx;
	int *ucindx;
	int *lcindx;
	EGlpNum_t *urcoef;
	EGlpNum_t *lccoef;
	EGlpNum_t pivot_coef;
	int nzcnt;
	int lc_freebeg;
	int s = f->stage;
	int i;
	int j;
	int rval = 0;

	EGlpNumInitVar (pivot_coef);

	if (uc_inf[c].nzcnt == 1)
	{
		/* col singleton */
		SETPERM (f, s, r, c);

		lc_inf[s].cbeg = -1;
		lc_inf[s].c = r;
		lc_inf[s].nzcnt = 0;
		f->stage++;

		urindx = f->urindx + ur_inf[r].rbeg;
		urcoef = f->urcoef + ur_inf[r].rbeg;
		nzcnt = ur_inf[r].nzcnt;
		for (i = 0; i < nzcnt; i++)
		{
			j = urindx[i];
			remove_col_nz (f, r, j);
			if (j == c)
			{
				urindx[i] = urindx[0];
				urindx[0] = c;
				EGLPNUM_SWAP (urcoef[0], urcoef[i], pivot_coef);
			}
		}
		remove_row (f, r);
		remove_col (f, c);
	}
	else if (ur_inf[r].nzcnt == 1)
	{
		/* row singleton */
		--(f->nstages);
		SETPERM (f, f->nstages, r, c);

		lc_inf[f->nstages].cbeg = -1;
		lc_inf[f->nstages].c = r;
		lc_inf[f->nstages].nzcnt = 0;

		ucindx = f->ucindx + uc_inf[c].cbeg;
		nzcnt = uc_inf[c].nzcnt;
		for (i = 0; i < nzcnt; i++)
		{
			j = ucindx[i];
			remove_row_nz (f, j, c);
		}
		remove_row (f, r);
		remove_col (f, c);
	}
	else
	{
		SETPERM (f, s, r, c);
		f->stage++;

		nzcnt = uc_inf[c].nzcnt;
		if (f->lc_freebeg + nzcnt >= f->lc_space)
		{
			rval = make_lc_space (f, nzcnt);
			CHECKRVALG (rval, CLEANUP);
		}
		lc_freebeg = f->lc_freebeg;
		lc_inf[s].cbeg = lc_freebeg;
		lc_inf[s].c = r;
		lcindx = f->lcindx;
		lccoef = f->lccoef;
		load_row (f, r);
		ucindx = f->ucindx + uc_inf[c].cbeg;
		for (i = 0; i < nzcnt; i++)
		{
			j = f->ucindx[uc_inf[c].cbeg + i];
			if (j != r)
			{
				rval = elim_row (f, r, j, c, &pivot_coef);
				CHECKRVALG (rval, CLEANUP);
				lcindx[lc_freebeg] = j;
				EGlpNumCopy (lccoef[lc_freebeg], pivot_coef);
				lc_freebeg++;
#ifdef TRACK_FACTOR
				EGlpNumSetToMaxAbs (f->maxelem_factor, pivot_coef);
				if (EGlpNumIsLess (f->maxelem_factor, ur_inf[r].max))
					EGlpNumCopy (f->maxelem_factor, ur_inf[r].max);
#endif /* TRACK_FACTOR */
			}
		}
		lc_inf[s].nzcnt = lc_freebeg - lc_inf[s].cbeg;
		f->lc_freebeg = lc_freebeg;

		clear_row (f, r);

		urindx = f->urindx + ur_inf[r].rbeg;
		urcoef = f->urcoef + ur_inf[r].rbeg;
		nzcnt = ur_inf[r].nzcnt;
		for (i = 0; i < nzcnt; i++)
		{
			j = urindx[i];
			remove_col_nz (f, r, j);
			if (j == c)
			{
				urindx[i] = urindx[0];
				urindx[0] = c;
				EGLPNUM_SWAP (urcoef[0], urcoef[i], pivot_coef);
			}
		}
		remove_row (f, r);
		remove_col (f, c);
	}
CLEANUP:
	EGlpNumClearVar (pivot_coef);
	EG_RETURN (rval);
}

static void find_pivot_column (
	factor_work * f,
	int c,
	int *p_r)
{
	uc_info *uc_inf = f->uc_inf;
	ur_info *ur_inf = f->ur_inf;
	int *ucindx = f->ucindx;
	int nzcnt = uc_inf[c].nzcnt;
	int cbeg = uc_inf[c].cbeg;
	EGlpNum_t num_tmp[2];
	int bestnz = -1;
	int i;
	int r;

	EGlpNumInitVar (num_tmp[0]);
	EGlpNumInitVar (num_tmp[1]);

	*p_r = -1;
	for (i = 0; i < nzcnt; i++)
	{
		r = ucindx[cbeg + i];
		if((bestnz == -1 || ur_inf[r].pivcnt < bestnz))
		{
			find_coef (f, r, c, num_tmp);
			if(EGlpNumIsLessZero(num_tmp[0]))
				EGlpNumSign (num_tmp[0]);
			EGlpNumCopy (num_tmp[1], f->partial_cur);
			EGlpNumMultTo (num_tmp[1], ur_inf[r].max);
			if(EGlpNumIsLeq (num_tmp[1], num_tmp[0]))
			{
				bestnz = ur_inf[r].pivcnt;
				*p_r = r;
			}
		}
	}
	EGlpNumClearVar (num_tmp[0]);
	EGlpNumClearVar (num_tmp[1]);
}

static void find_pivot_row (
	factor_work * f,
	int r,
	int *p_c)
{
	uc_info *uc_inf = f->uc_inf;
	ur_info *ur_inf = f->ur_inf;
	int *urindx = f->urindx;
	EGlpNum_t *urcoef = f->urcoef;
	int pivcnt = ur_inf[r].pivcnt;
	int rbeg = ur_inf[r].rbeg;
	EGlpNum_t thresh[2];
	int bestnz = -1;
	int i;
	int c;

	EGlpNumInitVar (thresh[0]);
	EGlpNumInitVar (thresh[1]);
	EGlpNumCopy (thresh[0], f->partial_cur);
	EGlpNumMultTo (thresh[0], ur_inf[r].max);
	*p_c = -1;
	for (i = 0; i < pivcnt; i++)
	{
		c = urindx[rbeg + i];
		if ((bestnz == -1 || uc_inf[c].nzcnt < bestnz))
		{
			EGlpNumCopyAbs (thresh[1], urcoef[rbeg + i]);
			if(EGlpNumIsLeq (thresh[0], thresh[1]))
			{
				bestnz = uc_inf[c].nzcnt;
				*p_c = c;
			}
		}
	}
	EGlpNumClearVar (thresh[0]);
	EGlpNumClearVar (thresh[1]);
}

static int find_pivot (
	factor_work * f,
	int *p_r,
	int *p_c)
{
	uc_info *uc_inf = f->uc_inf;
	ur_info *ur_inf = f->ur_inf;
	int dim = f->dim;
	int max_k = f->max_k;
	int p = f->p;
	int c;
	int r;
	int mm = 0;
	int n = 0;
	int m;
	int k = 2;

	if (uc_inf[dim + 1].next != dim + 1)
	{
		c = uc_inf[dim + 1].next;
		r = f->ucindx[uc_inf[c].cbeg];
		*p_c = c;
		*p_r = r;
		return 0;
	}
	else if (ur_inf[dim + 1].next != dim + 1)
	{
		r = ur_inf[dim + 1].next;
		c = f->urindx[ur_inf[r].rbeg];
		*p_c = c;
		*p_r = r;
		return 0;
	}
	*p_r = -1;
	*p_c = -1;
	for (; k <= max_k && (mm == 0 || mm > (k - 1) * (k - 1)); k++)
	{
		if (uc_inf[dim + k].next != dim + k)
		{
			for (c = uc_inf[dim + k].next; c != dim + k; c = uc_inf[c].next)
			{
				find_pivot_column (f, c, &r);
				if (r >= 0)
				{
					m = (uc_inf[c].nzcnt - 1) * (ur_inf[r].pivcnt - 1);
					if (mm == 0 || m < mm)
					{
						mm = m;
						*p_c = c;
						*p_r = r;
						if (mm <= (k - 1) * (k - 1))
						{
							return 0;
						}
					}
				}
				else
				{
					c = uc_inf[c].prev;
					disable_col (f, uc_inf[c].next);
				}
				n++;
				if (n >= p && mm != 0)
				{
					return 0;
				}
			}
		}

		if (ur_inf[dim + k].next != dim + k)
		{
			for (r = ur_inf[dim + k].next; r != dim + k; r = ur_inf[r].next)
			{
				find_pivot_row (f, r, &c);
				if (c >= 0)
				{
					m = (uc_inf[c].nzcnt - 1) * (ur_inf[r].pivcnt - 1);
					if (mm == 0 || m < mm)
					{
						mm = m;
						*p_c = c;
						*p_r = r;
						if (mm <= k * (k - 1))
						{
							return 0;
						}
					}
				}
				n++;
				if (n >= p && mm != 0)
				{
					return 0;
				}
			}
		}
	}
	if (mm != 0)
	{
		return 0;
	}
	else
	{
		//fprintf (stderr, "No acceptable pivot found\n");
		return E_NO_PIVOT;
	}
}

static int create_factor_space (
	factor_work * f)
{
	uc_info *uc_inf = f->uc_inf;
	ur_info *ur_inf = f->ur_inf;
	int dim = f->dim;
	int nzcnt;
	int i;
	int rval;

	nzcnt = 0;
	for (i = 0; i < dim; i++)
	{
		nzcnt += ur_inf[i].nzcnt;
	}

	if (f->ucindx == 0)
	{
		f->uc_space = nzcnt * f->uc_space_mul;
		ILL_SAFE_MALLOC (f->ucindx, f->uc_space + 1, int);
	}

	if (f->urindx == 0 || f->urcoef == 0)
	{
		ILL_IFFREE (f->urindx, int);

		EGlpNumFreeArray (f->urcoef);
		f->ur_space = nzcnt * f->ur_space_mul;
		ILL_SAFE_MALLOC (f->urindx, f->ur_space + 1, int);

		f->urcoef = EGlpNumAllocArray (f->ur_space);
	}

	if (f->lcindx == 0 || f->lccoef == 0)
	{
		ILL_IFFREE (f->lcindx, int);

		EGlpNumFreeArray (f->lccoef);
		f->lc_space = nzcnt * f->lc_space_mul;
		ILL_SAFE_MALLOC (f->lcindx, f->lc_space, int);

		f->lccoef = EGlpNumAllocArray (f->lc_space);
	}

	nzcnt = 0;
	for (i = 0; i < dim; i++)
	{
		ur_inf[i].rbeg = nzcnt;
		nzcnt += ur_inf[i].nzcnt;
		ur_inf[i].nzcnt = ur_inf[i].rbeg;
	}
	f->ur_freebeg = nzcnt;

	nzcnt = 0;
	for (i = 0; i < dim; i++)
	{
		uc_inf[i].cbeg = nzcnt;
		nzcnt += uc_inf[i].nzcnt;
		uc_inf[i].nzcnt = uc_inf[i].cbeg;
	}
	f->uc_freebeg = nzcnt;

	f->lc_freebeg = 0;

	rval = 0;
CLEANUP:
	EG_RETURN (rval);
}

static int init_matrix (
	factor_work * f,
	int *basis,
	int *cbeg,
	int *clen,
	int *in_ucindx,
	EGlpNum_t * in_uccoef)
{
	uc_info *uc_inf = f->uc_inf;
	ur_info *ur_inf = f->ur_inf;
	int dim = f->dim;
	int max_k = f->max_k;
	int *ucindx;
	int *urindx;
	EGlpNum_t *urcoef;
	int nzcnt;
	int beg;
	int i;
	int j;
	int r;
	int rval = 0;
	EGlpNum_t v;
	EGlpNum_t max;

	EGlpNumInitVar (v);
	EGlpNumInitVar (max);

	for (i = 0; i < dim; i++)
	{
		ur_inf[i].nzcnt = 0;
	}
	for (i = 0; i < dim; i++)
	{
		nzcnt = clen[basis[i]];
		beg = cbeg[basis[i]];
		uc_inf[i].nzcnt = nzcnt;
		for (j = 0; j < nzcnt; j++)
		{
			r = in_ucindx[beg + j];
			ur_inf[r].nzcnt++;
		}
	}

	rval = create_factor_space (f);
	CHECKRVALG (rval, CLEANUP);

	urindx = f->urindx;
	ucindx = f->ucindx;
	urcoef = f->urcoef;

	for (i = 0; i < dim; i++)
	{
		nzcnt = clen[basis[i]];
		beg = cbeg[basis[i]];
		for (j = 0; j < nzcnt; j++)
		{
			EGlpNumCopy (v, in_uccoef[beg + j]);
			if (!(EGlpNumIsNeqZero (v, f->fzero_tol)))
				continue;
			r = in_ucindx[beg + j];
			ucindx[uc_inf[i].nzcnt++] = r;
			urindx[ur_inf[r].nzcnt] = i;
			EGlpNumCopy (urcoef[ur_inf[r].nzcnt], v);
			ur_inf[r].nzcnt++;
		}
	}

	for (i = 0; i < dim; i++)
	{
		uc_inf[i].nzcnt -= uc_inf[i].cbeg;
		ur_inf[i].nzcnt -= ur_inf[i].rbeg;
	}

	j = f->uc_space;
	for (i = f->uc_freebeg; i < j; i++)
	{
		ucindx[i] = -1;
	}
	ucindx[j] = 0;

	j = f->ur_space;
	for (i = f->ur_freebeg; i < j; i++)
	{
		urindx[i] = -1;
	}
	urindx[j] = 0;

	for (i = 0; i < dim; i++)
	{
		nzcnt = ur_inf[i].nzcnt;
		ur_inf[i].pivcnt = nzcnt;
		beg = ur_inf[i].rbeg;
		EGlpNumZero (max);
		for (j = 0; j < nzcnt; j++)
		{
			EGlpNumSetToMaxAbs (max, urcoef[beg + j]);
		}
		EGlpNumCopy (ur_inf[i].max, max);
	}

	for (i = 0; i <= max_k; i++)
	{
		ur_inf[dim + i].next = dim + i;
		ur_inf[dim + i].prev = dim + i;
		uc_inf[dim + i].next = dim + i;
		uc_inf[dim + i].prev = dim + i;
	}

	for (i = 0; i < dim; i++)
	{
		nzcnt = uc_inf[i].nzcnt;
		if (nzcnt >= max_k)
			nzcnt = max_k;
		uc_inf[i].next = uc_inf[dim + nzcnt].next;
		uc_inf[i].prev = dim + nzcnt;
		uc_inf[dim + nzcnt].next = i;
		uc_inf[uc_inf[i].next].prev = i;

		nzcnt = ur_inf[i].pivcnt;
		if (nzcnt >= max_k)
			nzcnt = max_k;
		ur_inf[i].next = ur_inf[dim + nzcnt].next;
		ur_inf[i].prev = dim + nzcnt;
		ur_inf[dim + nzcnt].next = i;
		ur_inf[ur_inf[i].next].prev = i;
	}

#ifdef TRACK_FACTOR
	EGlpNumZero (max);
	nzcnt = 0;
	for (i = 0; i < dim; i++)
	{
		if (EGlpNumIsLess (max, ur_inf[i].max))
			EGlpNumCopy (max, ur_inf[i].max);
		nzcnt += ur_inf[i].nzcnt;
	}

	EGlpNumCopy (f->maxelem_orig, max);
	f->nzcnt_orig = nzcnt;
	EGlpNumCopy (f->maxelem_factor, f->maxelem_orig);
	f->nzcnt_factor = f->nzcnt_orig;
#endif /* TRACK_FACTOR */

	/* sentinal for column space */
	ucindx[f->uc_space] = 0;

	clear_work (f);

CLEANUP:
	EGlpNumClearVar (max);
	EGlpNumClearVar (v);
	EG_RETURN (rval);
}

static int build_iteration_u_data (
	factor_work * f)
{
	int dim = f->dim;
	ur_info *ur_inf = f->ur_inf;
	uc_info *uc_inf = f->uc_inf;
	EGlpNum_t *uccoef = 0;
	int *ucindx = 0;
	int *urindx = f->urindx;
	EGlpNum_t *urcoef = f->urcoef;
	int *ucrind = 0;
	int *urcind = 0;
	int nzcnt;
	int beg;
	int cbeg;
	int cnzcnt;
	int uc_space = f->uc_space;
	int er_space;
	int i;
	int j;
	int k;
	int rval;

	nzcnt = 0;
	for (i = 0; i < dim; i++)
	{
		nzcnt += ur_inf[i].nzcnt;
	}

#ifdef TRACK_FACTOR
	f->nzcnt_factor = nzcnt;
#endif /* TRACK_FACTOR */

	EGlpNumFreeArray (f->uccoef);
	uccoef = EGlpNumAllocArray (nzcnt);
	f->uccoef = uccoef;

	ILL_IFFREE (f->ucrind, int);
	ILL_SAFE_MALLOC (ucrind, nzcnt, int);

	f->ucrind = ucrind;

	ILL_IFFREE (f->urcind, int);
	ILL_SAFE_MALLOC (urcind, f->ur_space, int);

	f->urcind = urcind;

	if (uc_space < nzcnt)
	{
		ILL_IFFREE (f->ucindx, int);
		ILL_SAFE_MALLOC (f->ucindx, nzcnt + 1, int);
	}
	f->uc_space = nzcnt;
	uc_space = nzcnt;
	ucindx = f->ucindx;

	for (i = 0; i < dim; i++)
	{
		uc_inf[i].nzcnt = 0;
	}

	for (i = 0; i < dim; i++)
	{
		nzcnt = ur_inf[i].nzcnt;
		beg = ur_inf[i].rbeg;
		for (j = 0; j < nzcnt; j++)
		{
			uc_inf[urindx[beg + j]].nzcnt++;
		}
		ur_inf[i].delay = 0;
	}

	nzcnt = 0;
	for (i = 0; i < dim; i++)
	{
		uc_inf[i].cbeg = nzcnt;
		nzcnt += uc_inf[i].nzcnt;
		uc_inf[i].nzcnt = 0;
		uc_inf[i].delay = 0;
	}

	f->uc_freebeg = nzcnt;
	for (i = nzcnt; i < uc_space; i++)
	{
		ucindx[i] = -1;
	}
	ucindx[uc_space] = 0;

	for (i = 0; i < dim; i++)
	{
		nzcnt = ur_inf[i].nzcnt;
		beg = ur_inf[i].rbeg;
		k = urindx[beg];
		cbeg = uc_inf[k].cbeg;
		cnzcnt = uc_inf[k].nzcnt;
		if (cnzcnt != 0)
		{
			ucindx[cbeg + cnzcnt] = ucindx[cbeg];
			EGlpNumCopy (uccoef[cbeg + cnzcnt], uccoef[cbeg]);
			ucrind[cbeg + cnzcnt] = ucrind[cbeg];
			urcind[ur_inf[ucindx[cbeg]].rbeg + ucrind[cbeg]] = cnzcnt;
		}
		ucindx[cbeg] = i;
		EGlpNumCopy (uccoef[cbeg], urcoef[beg]);
		ucrind[cbeg] = 0;
		urcind[beg] = 0;
		uc_inf[k].nzcnt = cnzcnt + 1;
		for (j = 1; j < nzcnt; j++)
		{
			k = urindx[beg + j];
			cbeg = uc_inf[k].cbeg;
			cnzcnt = uc_inf[k].nzcnt;
			ucindx[cbeg + cnzcnt] = i;
			EGlpNumCopy (uccoef[cbeg + cnzcnt], urcoef[beg + j]);
			ucrind[cbeg + cnzcnt] = j;
			urcind[beg + j] = cnzcnt;
			uc_inf[k].nzcnt++;
		}
	}

	for (i = 0; i < dim; i++)
	{
		f->rrank[f->rperm[i]] = i;
	}

	nzcnt = f->ur_space;

	for (i = f->ur_freebeg; i < nzcnt; i++)
	{
		urindx[i] = -1;
	}
	urindx[nzcnt] = 0;

	clear_work (f);

	er_space = f->er_space_mul * f->etamax;
	ILL_SAFE_MALLOC (f->er_inf, f->etamax, er_info);
	ILL_SAFE_MALLOC (f->erindx, er_space, int);

	f->ercoef = EGlpNumAllocArray (er_space);
	f->etacnt = 0;
	f->er_freebeg = 0;
	f->er_space = er_space;

	rval = 0;

CLEANUP:
	EG_RETURN (rval);
}

static int build_iteration_l_data (
	factor_work * f)
{
	int dim = f->dim;
	lc_info *lc_inf = f->lc_inf;
	lr_info *lr_inf = f->lr_inf;
	EGlpNum_t *lrcoef = 0;
	int *lrindx = 0;
	EGlpNum_t *lccoef = f->lccoef;
	int *lcindx = f->lcindx;
	int nzcnt;
	int beg;
	int rnzcnt;
	int rbeg;
	int i;
	int j;
	int k;
	int c;
	int rval;

	nzcnt = 0;
	for (i = 0; i < dim; i++)
	{
		nzcnt += lc_inf[i].nzcnt;
		lr_inf[i].nzcnt = 0;
		lr_inf[i].delay = 0;
		lc_inf[lc_inf[i].c].crank = i;
	}

	EGlpNumFreeArray (f->lrcoef);
	if (nzcnt)
	{
		lrcoef = EGlpNumAllocArray (nzcnt);
		f->lrcoef = lrcoef;
	}

	ILL_IFFREE (f->lrindx, int);
	ILL_SAFE_MALLOC (lrindx, nzcnt + 1, int);

	f->lrindx = lrindx;

	for (i = 0; i < dim; i++)
	{
		nzcnt = lc_inf[i].nzcnt;
		beg = lc_inf[i].cbeg;
		lc_inf[i].delay = 0;
		for (j = 0; j < nzcnt; j++)
		{
			lr_inf[lc_inf[lcindx[beg + j]].crank].nzcnt++;
		}
	}

	nzcnt = 0;
	for (i = 0; i < dim; i++)
	{
		lr_inf[i].rbeg = nzcnt;
		nzcnt += lr_inf[i].nzcnt;
		lr_inf[i].nzcnt = 0;
		lr_inf[i].r = lc_inf[i].c;
		lr_inf[lr_inf[i].r].rrank = i;
	}

	for (i = 0; i < dim; i++)
	{
		nzcnt = lc_inf[i].nzcnt;
		beg = lc_inf[i].cbeg;
		c = lc_inf[i].c;
		for (j = 0; j < nzcnt; j++)
		{
			k = lc_inf[lcindx[beg + j]].crank;
			rbeg = lr_inf[k].rbeg;
			rnzcnt = lr_inf[k].nzcnt;
			lrindx[rbeg + rnzcnt] = c;
			EGlpNumCopy (lrcoef[rbeg + rnzcnt], lccoef[beg + j]);
			lr_inf[k].nzcnt++;
		}
	}

#ifdef TRACK_FACTOR
	nzcnt = f->nzcnt_factor;
	for (i = 0; i < dim; i++)
	{
		nzcnt += lc_inf[i].nzcnt;
	}
	f->nzcnt_factor = nzcnt;

	EGlpNumCopy (f->maxelem_cur, f->maxelem_factor);
	f->nzcnt_cur = f->nzcnt_factor;

/*
    dump_factor_stats (f);
    printf ("orig max  %e nzcnt %d\n", f->maxelem_orig, f->nzcnt_orig);
    printf ("f maxelem %e nzcnt %d\n", f->maxelem_cur, f->nzcnt_cur);
*/
#endif /* TRACK_FACTOR */

	rval = 0;

CLEANUP:
	EG_RETURN (rval);
}

static int handle_singularity (
	factor_work * f)
{
	int rval = 0;
	int nsing;
	int *singr = 0;
	int *singc = 0;
	int i;

	if (f->p_nsing == 0 || f->p_singr == 0 || f->p_singc == 0)
	{
		fprintf (stderr, "singular basis, but no place for singularity data\n");
		return E_SING_NO_DATA;
	}

	nsing = f->nstages - f->stage;
	ILL_SAFE_MALLOC (singr, nsing, int);
	ILL_SAFE_MALLOC (singc, nsing, int);

	for (i = f->stage; i < f->nstages; i++)
	{
		singr[i - f->stage] = f->rperm[i];
		singc[i - f->stage] = f->cperm[i];
	}
	*f->p_nsing = nsing;
	*f->p_singr = singr;
	*f->p_singc = singc;
	singr = 0;
	singc = 0;

CLEANUP:
	ILL_IFFREE (singr, int);
	ILL_IFFREE (singc, int);

	EG_RETURN (rval);
}

static int dense_build_matrix (
	factor_work * f)
{
	EGlpNum_t *dmat = 0;
	int stage = f->stage;
	int drows = f->nstages - stage;
	int dcols = f->dim - stage;
	int dsize = drows * dcols;
	int *crank = f->crank;
	EGlpNum_t *urcoef = f->urcoef;
	int *urindx = f->urindx;
	int nzcnt;
	int beg;
	int i;
	int r;
	int j;
	int rval = 0;

	dmat = EGlpNumAllocArray (dsize);

	for (i = 0; i < dsize; i++)
		EGlpNumZero (dmat[i]);

	for (i = 0; i < drows; i++)
	{
		r = f->rperm[i + stage];
		nzcnt = f->ur_inf[r].nzcnt;
		beg = f->ur_inf[r].rbeg;
		for (j = 0; j < nzcnt; j++)
		{
			EGlpNumCopy (dmat[i * dcols - stage + crank[urindx[beg + j]]],
									 urcoef[beg + j]);
		}
	}

	f->drows = drows;
	f->dcols = dcols;
	f->dense_base = f->stage;
	f->dmat = dmat;
	dmat = 0;

//CLEANUP:
	EGlpNumFreeArray (dmat);
	EG_RETURN (rval);
}

static int dense_find_pivot (
	factor_work * f,
	int *p_r,
	int *p_c)
{
	int dcols = f->dcols;
	int drows = f->drows;
	EGlpNum_t *dmat = f->dmat;
	int dense_base = f->dense_base;
	int s = f->stage - dense_base;
	ur_info *ur_inf = f->ur_inf;
	int *rperm = f->rperm;
	EGlpNum_t maxval;
	int max_r;
	int max_c;
	int i;

	EGlpNumInitVar (maxval);
	EGlpNumZero (maxval);
	max_r = -1;
	for (i = s; i < drows; i++)
	{
		if (EGlpNumIsLess (maxval, ur_inf[rperm[dense_base + i]].max))
		{
			EGlpNumCopy (maxval, ur_inf[rperm[dense_base + i]].max);
			max_r = i;
		}
	}
	if (max_r == -1)
	{
		return E_NO_PIVOT;
	}

	EGlpNumZero (maxval);
	max_c = -1;
	for (i = s; i < drows; i++)
	{
		EGlpNumSetToMaxAbsAndDo (maxval, dmat[max_r * dcols + i], max_c = i);
	}
	if (max_c == -1)
	{
		return E_NO_PIVOT;
	}
	*p_r = max_r;
	*p_c = max_c;

	EGlpNumClearVar (maxval);
	return 0;
}

static void dense_swap (
	factor_work * f,
	int r,
	int c)
{
	int dcols = f->dcols;
	int drows = f->drows;
	EGlpNum_t *dmat = f->dmat;
	int dense_base = f->dense_base;
	int s = f->stage - dense_base;
	int i;
	EGlpNum_t v;

	EGlpNumInitVar (v);

	if (r != s)
	{
		ILL_SWAP (f->rperm[dense_base + s], f->rperm[dense_base + r], i);
		f->rrank[f->rperm[dense_base + s]] = dense_base + s;
		f->rrank[f->rperm[dense_base + r]] = dense_base + r;
		for (i = 0; i < dcols; i++)
		{
			EGLPNUM_SWAP (dmat[s * dcols + i], dmat[r * dcols + i], v);
		}
	}
	if (c != s)
	{
		ILL_SWAP (f->cperm[dense_base + s], f->cperm[dense_base + c], i);
		f->crank[f->cperm[dense_base + s]] = dense_base + s;
		f->crank[f->cperm[dense_base + c]] = dense_base + c;
		for (i = 0; i < drows; i++)
		{
			EGLPNUM_SWAP (dmat[i * dcols + s], dmat[i * dcols + c], v);
		}
	}
	EGlpNumClearVar (v);
}

static void dense_elim (
	factor_work * f,
	int r,
	int c)
{
	int dcols = f->dcols;
	int drows = f->drows;
	EGlpNum_t *dmat = f->dmat;
	int dense_base = f->dense_base;
	int s = f->stage - dense_base;
	ur_info *ur_inf = f->ur_inf;
	int *rperm = f->rperm;
	int i;
	int j;
	EGlpNum_t pivval;
	EGlpNum_t max;
	EGlpNum_t v;
	EGlpNum_t w;

#ifdef TRACK_FACTOR
	EGlpNum_t maxelem_factor;

	EGlpNumInitVar (maxelem_factor);
	EGlpNumCopy (maxelem_factor, f->maxelem_factor);
#endif
	EGlpNumInitVar (pivval);
	EGlpNumInitVar (max);
	EGlpNumInitVar (v);
	EGlpNumInitVar (w);

	dense_swap (f, r, c);
	f->stage++;
	EGlpNumCopyFrac (pivval, oneLpNum, dmat[s * dcols + s]);
	for (i = s + 1; i < drows; i++)
	{
		EGlpNumCopy (v, dmat[i * dcols + s]);
		if (EGlpNumIsNeqqZero (v))
		{
			EGlpNumMultTo (v, pivval);
			if (EGlpNumIsNeqZero (v, f->fzero_tol))
			{
				EGlpNumCopy (dmat[i * dcols + s], v);
#ifdef TRACK_FACTOR
				EGlpNumSetToMaxAbs (maxelem_factor, v);
#endif
				EGlpNumZero (max);
				for (j = s + 1; j < drows; j++)
				{
					EGlpNumCopy (w, dmat[i * dcols + j]);
					EGlpNumSubInnProdTo (w, v, dmat[s * dcols + j]);
					EGlpNumCopy (dmat[i * dcols + j], w);
					EGlpNumSetToMaxAbs (max, w);
				}
				for (j = drows; j < dcols; j++)
				{
					EGlpNumCopy (w, dmat[i * dcols + j]);
					EGlpNumSubInnProdTo (w, v, dmat[s * dcols + j]);
					EGlpNumCopy (dmat[i * dcols + j], w);
				}
				EGlpNumCopy (ur_inf[rperm[dense_base + i]].max, max);
#ifdef TRACK_FACTOR
				if (EGlpNumIsLess (maxelem_factor, max))
					EGlpNumCopy (maxelem_factor, max);
#endif
			}
			else
			{
				EGlpNumZero (dmat[i * dcols + s]);
			}
		}
	}
#ifdef TRACK_FACTOR
	EGlpNumCopy (f->maxelem_factor, maxelem_factor);
	EGlpNumClearVar (maxelem_factor);
#endif
	EGlpNumClearVar (pivval);
	EGlpNumClearVar (max);
	EGlpNumClearVar (v);
	EGlpNumClearVar (w);
}

static int dense_replace_row (
	factor_work * f,
	int i)
{
	int dcols = f->dcols;
	int dense_base = f->dense_base;
	EGlpNum_t *dmat = f->dmat + i * dcols;
	EGlpNum_t *urcoef;
	ur_info *ur_inf = f->ur_inf;
	int *cperm = f->cperm;
	int r = f->rperm[dense_base + i];
	int *urindx;
	int nzcnt;
	int beg;
	int j;
	int rval = 0;

	nzcnt = 0;
	for (j = i; j < dcols; j++)
	{
		if (EGlpNumIsNeqZero (dmat[j], f->fzero_tol))
		{
			nzcnt++;
		}
	}
	if (nzcnt > ur_inf[r].nzcnt)
	{
		if (ur_inf[r].rbeg + ur_inf[r].nzcnt == f->ur_freebeg)
		{
			f->ur_freebeg = ur_inf[r].rbeg;
		}
		ur_inf[r].nzcnt = 0;
		if (f->ur_freebeg + nzcnt > f->ur_space)
		{
			rval = make_ur_space (f, nzcnt);
			CHECKRVALG (rval, CLEANUP);
		}
		ur_inf[r].rbeg = f->ur_freebeg;
		f->ur_freebeg += nzcnt;
	}
	beg = ur_inf[r].rbeg;
	urcoef = f->urcoef;
	urindx = f->urindx;
	for (j = i; j < dcols; j++)
	{
		if (EGlpNumIsNeqZero (dmat[j], f->fzero_tol))
		{
			EGlpNumCopy (urcoef[beg], dmat[j]);
			urindx[beg] = cperm[dense_base + j];
			beg++;
		}
	}
	ur_inf[r].nzcnt = beg - ur_inf[r].rbeg;
CLEANUP:
	EG_RETURN (rval);
}

static int dense_create_col (
	factor_work * f,
	int i)
{
	int dcols = f->dcols;
	int drows = f->drows;
	int dense_base = f->dense_base;
	EGlpNum_t *dmat = f->dmat;
	EGlpNum_t *lccoef;
	lc_info *lc_inf = f->lc_inf;
	int *rperm = f->rperm;
	int *lcindx;
	int nzcnt;
	int beg;
	int j;
	int rval = 0;

	nzcnt = 0;
	for (j = i + 1; j < drows; j++)
	{
		if (EGlpNumIsNeqZero (dmat[j * dcols + i], f->fzero_tol))
		{
			nzcnt++;
		}
	}

	if (f->lc_freebeg + nzcnt >= f->lc_space)
	{
		rval = make_lc_space (f, nzcnt);
		CHECKRVALG (rval, CLEANUP);
	}
	beg = f->lc_freebeg;
	lc_inf[dense_base + i].cbeg = beg;
	lc_inf[dense_base + i].c = rperm[dense_base + i];
	lcindx = f->lcindx;
	lccoef = f->lccoef;

	for (j = i + 1; j < drows; j++)
	{
		if (EGlpNumIsNeqZero (dmat[j * dcols + i], f->fzero_tol))
		{
			EGlpNumCopy (lccoef[beg], dmat[j * dcols + i]);
			lcindx[beg] = rperm[dense_base + j];
			beg++;
		}
	}
	lc_inf[dense_base + i].nzcnt = beg - lc_inf[dense_base + i].cbeg;
	f->lc_freebeg = beg;
CLEANUP:
	EG_RETURN (rval);
}

static int dense_replace (
	factor_work * f)
{
	int drows = f->drows;
	int rval = 0;
	int i;

	for (i = 0; i < drows; i++)
	{
		rval = dense_replace_row (f, i);
		CHECKRVALG (rval, CLEANUP);
		rval = dense_create_col (f, i);
		CHECKRVALG (rval, CLEANUP);
	}
	EGlpNumFreeArray (f->dmat);
	f->drows = 0;
	f->dcols = 0;
CLEANUP:
	EG_RETURN (rval);
}

static int dense_factor (
	factor_work * f)
{
	int r;
	int c;
	int rval = 0;

#ifdef TRACK_FACTOR
#ifdef NOTICE_BLOWUP
	double tmpsize;
#endif
#endif

/*
    printf ("dense kernel, %d rows, %d  cols...\n", f->nstages - f->stage,
            f->dim - f->stage);
    fflush (stdout);
*/

	rval = dense_build_matrix (f);
	CHECKRVALG (rval, CLEANUP);

#ifdef FACTOR_DEBUG
#if (FACTOR_DEBUG+0>1)
	MESSAGE (0,"before Dense ILLfactor");
	dump_matrix (f, 1);
#endif
#endif

	while (f->stage < f->nstages)
	{
		r = f->stage - f->dense_base;
		rval = dense_find_pivot (f, &r, &c);
		if (rval == E_NO_PIVOT)
		{
			rval = handle_singularity (f);
			CHECKRVALG (rval, CLEANUP);
			return E_SINGULAR_INTERNAL;
		}
		else
		{
			CHECKRVALG (rval, CLEANUP);
		}
#ifdef FACTOR_DEBUG
#if (FACTOR_DEBUG+0>2)
		MESSAGE (0,"dense pivot elem: %d %d", r, c);
#endif
#endif /* FACTOR_DEBUG */
		dense_elim (f, r, c);

#ifdef TRACK_FACTOR
#ifdef NOTICE_BLOWUP
		tmpsize = f->maxmult * EGlpNumToLf (f->maxelem_orig);
		if (tmpsize < EGlpNumToLf (f->maxelem_factor) &&
				EGlpNumIsLess (f->partial_cur, oneLpNum))
		{
			return E_FACTOR_BLOWUP;
		}
#endif /* NOTICE_BLOWUP */
#endif /* TRACK_FACTOR */

#ifdef FACTOR_DEBUG
#if (FACTOR_DEBUG+0>1)
		MESSAGE (0,"After dense pivot stage %d (%d) of %d (%d)",
						f->stage - f->dense_base, f->stage,
						f->nstages - f->dense_base, f->nstages);
#endif
#if (FACTOR_DEBUG+0>2)
		dump_matrix (f, 1);
#endif
#endif /* FACTOR_DEBUG */
	}

#ifdef FACTOR_DEBUG
	MESSAGE (0,"After dense ILLfactor:\n");
	dump_matrix (f, 0);
#endif /* FACTOR_DEBUG */

	rval = dense_replace (f);
	CHECKRVALG (rval, CLEANUP);

#ifdef FACTOR_DEBUG
	MESSAGE (0,"After replacement:\n");
	dump_matrix (f, 0);
#endif /* FACTOR_DEBUG */

CLEANUP:
	EG_RETURN (rval);
}

#ifdef RECORD
EGioFile_t *fsave = 0;
int fsavecnt = 0;
#endif /* RECORD */

static int ILLfactor_try (
	factor_work * f,
	int *basis,
	int *cbeg,
	int *clen,
	int *cindx,
	EGlpNum_t * ccoef)
{
	int rval = 0;
	int r;
	int c;

#ifdef TRACK_FACTOR
#ifdef NOTICE_BLOWUP
	EGlpNum_t tmpsize;

	EGlpNumInitVar (tmpsize);
#endif
#endif

#ifdef RECORD
	{
		int ncol = 0;
		int nzcnt = 0;
		int dim = f->dim;
		int i;
		int j;
		char fnambuf[40];

		for (i = 0; i < dim; i++)
		{
			if (basis[i] > ncol)
				ncol = basis[i];
		}
		ncol++;
		for (i = 0; i < ncol; i++)
		{
			nzcnt += clen[i];
		}
		if (fsave)
			EGioClose (fsave);
		#if HAVE_ZLIB_H
		sprintf (fnambuf, "prob.mat.%d.gz", fsavecnt);
		#elif HAVE_BZLIB_H
		sprintf (fnambuf, "prob.mat.%d.bz2", fsavecnt);
		#else
		sprintf (fnambuf, "prob.mat.%d", fsavecnt);
		#endif
		fsavecnt++;
		fsave = EGioOpen (fnambuf, "w");
		EGioPrintf (fsave, "%d %d %d\n", f->dim, ncol, nzcnt);
		for (i = 0; i < dim; i++)
		{
			EGioPrintf (fsave, "%d ", basis[i]);
		}
		EGioPrintf (fsave, "\n");
		for (i = 0; i < ncol; i++)
		{
			EGioPrintf (fsave, "%d", clen[i]);
			for (j = 0; j < clen[i]; j++)
			{
				EGioPrintf (fsave, " %d %.16lg", cindx[cbeg[i] + j],
								 EGlpNumToLf (ccoef[cbeg[i] + j]));
			}
			EGioPrintf (fsave, "\n");
		}
		EGioPrintf (fsave, "\n");
		EGioFlush (fsave);
	}
#endif /* RECORD */

	rval = init_matrix (f, basis, cbeg, clen, cindx, ccoef);
	CHECKRVALG (rval, CLEANUP);

	f->stage = 0;
	f->nstages = f->dim;

#ifdef FACTOR_DEBUG
	MESSAGE (0,"Initial matrix:");
#if (FACTOR_DEBUG+0>1)
	dump_matrix (f, 0);
#endif
#endif /* FACTOR_DEBUG */
#ifdef FACTOR_STATS
	printf ("Initial matrix: ");
	dump_factor_stats (f);
#endif /* FACTOR_STATS */

	while (f->stage < f->nstages)
	{
		rval = find_pivot (f, &r, &c);
		if (rval == E_NO_PIVOT)
		{
			rval = handle_singularity (f);
			CHECKRVALG (rval, CLEANUP);
			return 0;
		}
		else
		{
			CHECKRVALG (rval, CLEANUP);
		}
		if (f->ur_inf[r].pivcnt > f->dense_fract * (f->nstages - f->stage) &&
				f->uc_inf[c].nzcnt > f->dense_fract * (f->nstages - f->stage) &&
				f->nstages - f->stage > f->dense_min)
		{
			rval = dense_factor (f);
			if (rval == E_SINGULAR_INTERNAL)
				return 0;
			if (rval)
				return rval;
			break;
		}
#ifdef FACTOR_DEBUG
		MESSAGE (0,"pivot elem: %d %d", r, c);
#endif /* FACTOR_DEBUG */
		rval = elim (f, r, c);
		CHECKRVALG (rval, CLEANUP);

#ifdef TRACK_FACTOR
#ifdef NOTICE_BLOWUP
		EGlpNumSet (tmpsize, f->maxmult);
		EGlpNumMultTo (tmpsize, f->maxelem_orig);
		if (EGlpNumIsLess (tmpsize, f->maxelem_factor) &&
				EGlpNumIsLess (f->partial_cur, oneLpNum))
		{
			return E_FACTOR_BLOWUP;
		}
#endif /* NOTICE_BLOWUP */
#endif /* TRACK_FACTOR */

#ifdef FACTOR_DEBUG
#if (FACTOR_DEBUG+0>3)
		MESSAGE (0,"After pivot stage %d of %d", f->stage, f->nstages);
		dump_matrix (f, 0);
#endif
#endif /* FACTOR_DEBUG */
	}

	rval = build_iteration_u_data (f);
	CHECKRVALG (rval, CLEANUP);

	rval = build_iteration_l_data (f);
	CHECKRVALG (rval, CLEANUP);

#ifdef TRACK_FACTOR
#ifdef NOTICE_BLOWUP
	EGlpNumSet (tmpsize, f->minmult);
	EGlpNumMultTo (tmpsize, f->maxelem_orig);
	if (EGlpNumIsLess (f->maxelem_factor, tmpsize) &&
			EGlpNumIsLess (f->partial_tol, f->partial_cur))
	{
		if (EGlpNumIsGreaDbl (f->partial_cur, 0.5))
		{
			EGlpNumSet (f->partial_cur, 0.5);
		}
		else if (EGlpNumIsGreaDbl (f->partial_cur, 0.25))
		{
			EGlpNumSet (f->partial_cur, 0.25);
		}
		else if (EGlpNumIsGreaDbl (f->partial_cur, 0.1))
		{
			EGlpNumSet (f->partial_cur, 0.1);
		}
		else
		{
			EGlpNumDivUiTo (f->partial_cur, 10);
		}
		if (EGlpNumIsLess (f->partial_cur, f->partial_tol))
		{
			EGlpNumCopy (f->partial_cur, f->partial_tol);
		}
/*  Bico - comment out for dist 
        fprintf (stderr, "factor good, lowering partial tolerance to %.2f\n",
                 f->partial_cur);
*/
	}
#endif /* NOTICE_BLOWUP */
#endif /* TRACK_FACTOR */

#ifdef FACTOR_DEBUG
	MESSAGE(0,"Factored matrix:");
#if (FACTOR_DEBUG+0>1)
	dump_matrix (f, 0);
#endif
#endif /* FACTOR_DEBUG */

#ifdef FACTOR_STATS
	printf ("Factored matrix: ");
	dump_factor_stats (f);
#endif /* FACTOR_STATS */
CLEANUP:
#ifdef TRACK_FACTOR
#ifdef NOTICE_BLOWUP
	EGlpNumClearVar (tmpsize);
#endif
#endif
	EG_RETURN (rval);
}

int ILLfactor (
	factor_work * f,
	int *basis,
	int *cbeg,
	int *clen,
	int *cindx,
	EGlpNum_t * ccoef,
	int *p_nsing,
	int **p_singr,
	int **p_singc)
{
	int rval;

	f->p_nsing = p_nsing;
	f->p_singr = p_singr;
	f->p_singc = p_singc;
	*p_nsing = 0;

AGAIN:
	rval = ILLfactor_try (f, basis, cbeg, clen, cindx, ccoef);
	if (rval == E_FACTOR_BLOWUP)
	{
		if (EGlpNumIsLessDbl (f->partial_cur, 0.1))
		{
			EGlpNumMultUiTo (f->partial_cur, 10);
		}
		else if (EGlpNumIsLessDbl (f->partial_cur, 0.25))
		{
			EGlpNumSet (f->partial_cur, 0.25);
		}
		else if (EGlpNumIsLessDbl (f->partial_cur, 0.5))
		{
			EGlpNumSet (f->partial_cur, 0.5);
		}
		else if (EGlpNumIsLess (f->partial_cur, oneLpNum))
		{
			EGlpNumOne (f->partial_cur);
		}
		else
		{
			EG_RETURN (rval);
		}
/* Bico - comment out for dist
        fprintf (stderr, "factor blowup, changing partial tolerance to %.2f\n",
                 f->partial_cur);
*/
		goto AGAIN;
	}
	EG_RETURN (rval);
}

static void ILLfactor_ftranl (
	factor_work * f,
	EGlpNum_t * a)
{
	int *lcindx = f->lcindx;
	lc_info *lc_inf = f->lc_inf;
	EGlpNum_t *lccoef = f->lccoef;
	int dim = f->dim;
	int beg;
	int nzcnt;
	int i;
	int j;
	EGlpNum_t v;

	EGlpNumInitVar (v);

	for (i = 0; i < dim; i++)
	{
		EGlpNumCopy (v, a[lc_inf[i].c]);
		if (EGlpNumIsNeqqZero (v))
		{
			nzcnt = lc_inf[i].nzcnt;
			beg = lc_inf[i].cbeg;
			for (j = 0; j < nzcnt; j++)
			{
				EGlpNumSubInnProdTo (a[lcindx[beg + j]], v, lccoef[beg + j]);
			}
		}
#ifdef SOLVE_DEBUG
#if (SOLVE_DEBUG+0 > 1)
		fpritf (stderr,"ILLfactor_ftran a after l %d:", i);
		for (j = 0; j < f->dim; j++)
		{
			fprintf (stderr," %.3f", EGlpNumToLf (a[j]));
		}
		MESSAGE(0," ");
#endif
#endif /* SOLVE_DEBUG */
	}
#ifdef SOLVE_DEBUG
#if (SOLVE_DEBUG+0 <= 1)
	printf ("ILLfactor_ftran a after l:");
	for (j = 0; j < f->dim; j++)
	{
		printf (" %.3f", EGlpNumToLf (a[j]));
	}
	printf ("\n");
	fflush (stdout);
#endif
#endif /* SOLVE_DEBUG */
	EGlpNumClearVar (v);
}

#if 0
static void ftranl3_delay (
	factor_work * f,
	int c)
{
	lc_info *lc_inf = f->lc_inf;
	int nzcnt;
	int *indx;
	int i;

	c = lc_inf[c].crank;
	nzcnt = lc_inf[c].nzcnt;
	indx = f->lcindx + lc_inf[c].cbeg;
	for (i = 0; i < nzcnt; i++)
	{
		c = indx[i];
		if (lc_inf[c].delay++ == 0)
		{
			ftranl3_delay (f, c);
		}
	}
}
#endif

static void ftranl3_delay2 (
	factor_work * f,
	int c)
{
	lc_info *lc_inf = f->lc_inf;
	int nzcnt;
	int *indx;
	int i;
	int last;

	do
	{
		c = lc_inf[c].crank;
		nzcnt = lc_inf[c].nzcnt;
		indx = f->lcindx + lc_inf[c].cbeg;
		last = -1;
		for (i = 0; i < nzcnt; i++)
		{
			c = indx[i];
			if (lc_inf[c].delay++ == 0)
			{
				if (last >= 0)
				{
					ftranl3_delay2 (f, last);
				}
				last = c;
			}
		}
		c = last;
	} while (c >= 0);
}

#if 0
static void ftranl3_process (
	factor_work * f,
	int c,
	svector * x)
{
	lc_info *lc_inf = f->lc_inf;
	EGlpNum_t *work = f->work_coef;
	int nzcnt;
	int *indx;
	int i;
	EGlpNum_t *coef;
	EGlpNum_t v;

	EGlpNumInitVar (v);

	EGlpNumCopy (v, work[c]);
	EGlpNumZero (work[c]);
	if (EGlpNumIsNeqqZero (v))
	{
		x->indx[x->nzcnt] = c;
		EGlpNumCopy (x->coef[x->nzcnt], v);
		x->nzcnt++;
	}
	c = lc_inf[c].crank;
	nzcnt = lc_inf[c].nzcnt;
	indx = f->lcindx + lc_inf[c].cbeg;
	coef = f->lccoef + lc_inf[c].cbeg;
	for (i = 0; i < nzcnt; i++)
	{
		c = indx[i];
		EGlpNumSubInnProdTo (work[c], v, coef[i]);
		if (--lc_inf[c].delay == 0)
		{
			ftranl3_process (f, c, x);
		}
	}
	EGlpNumClearVar (v);
}
#endif

static void ftranl3_process2 (
	factor_work * f,
	int c,
	svector * x)
{
	lc_info *lc_inf = f->lc_inf;
	EGlpNum_t *work = f->work_coef;
	int nzcnt;
	int *indx;
	EGlpNum_t *coef;
	int i;
	int last;
	EGlpNum_t v;

	EGlpNumInitVar (v);

	do
	{
		EGlpNumCopy (v, work[c]);
		EGlpNumZero (work[c]);
		if (EGlpNumIsNeqqZero (v))
		{
			x->indx[x->nzcnt] = c;
			EGlpNumCopy (x->coef[x->nzcnt], v);
			x->nzcnt++;
		}
		c = lc_inf[c].crank;
		nzcnt = lc_inf[c].nzcnt;
		indx = f->lcindx + lc_inf[c].cbeg;
		coef = f->lccoef + lc_inf[c].cbeg;
		last = -1;
		for (i = 0; i < nzcnt; i++)
		{
			c = indx[i];
			EGlpNumSubInnProdTo (work[c], v, coef[i]);
			if (--lc_inf[c].delay == 0)
			{
				if (last >= 0)
				{
					ftranl3_process2 (f, last, x);
				}
				last = c;
			}
		}
		c = last;
	} while (c >= 0);
	EGlpNumClearVar (v);
}

static void ILLfactor_ftranl3 (
	factor_work * f,
	svector * a,
	svector * x)
{
	EGlpNum_t *work = f->work_coef;
	int anzcnt = a->nzcnt;
	int *aindx = a->indx;
	EGlpNum_t *acoef = a->coef;
	lc_info *lc_inf = f->lc_inf;
	int i;

	for (i = 0; i < anzcnt; i++)
	{
		if (lc_inf[aindx[i]].delay++ == 0)
		{
			ftranl3_delay2 (f, aindx[i]);
		}
		EGlpNumCopy (work[aindx[i]], acoef[i]);
	}
	x->nzcnt = 0;
	for (i = 0; i < anzcnt; i++)
	{
		if (--lc_inf[aindx[i]].delay == 0)
		{
			ftranl3_process2 (f, aindx[i], x);
		}
	}
#ifdef SOLVE_DEBUG
	printf ("ILLfactor_ftran x after l3:");
	for (i = 0; i < x->nzcnt; i++)
	{
		printf (" %.3f*%d", EGlpNumToLf (x->coef[i]), x->indx[i]);
	}
	printf ("\n");
	fflush (stdout);
#endif /* SOLVE_DEBUG */
}

static void ILLfactor_ftrane (
	factor_work * f,
	EGlpNum_t * a)
{
	int *erindx = f->erindx;
	EGlpNum_t *ercoef = f->ercoef;
	er_info *er_inf = f->er_inf;
	int etacnt = f->etacnt;
	int beg;
	int nzcnt;
	int i;
	int j;
	EGlpNum_t v;

	EGlpNumInitVar (v);

	for (i = 0; i < etacnt; i++)
	{
		EGlpNumCopy (v, a[er_inf[i].r]);
		nzcnt = er_inf[i].nzcnt;
		beg = er_inf[i].rbeg;
		for (j = 0; j < nzcnt; j++)
		{
			EGlpNumSubInnProdTo (v, ercoef[beg + j], a[erindx[beg + j]]);
		}
		EGlpNumCopy (a[er_inf[i].r], v);
#ifdef SOLVE_DEBUG
#if (SOLVE_DEBUG+0 > 1)
		printf ("ILLfactor_ftran a after eta %d:", i);
		for (j = 0; j < f->dim; j++)
		{
			printf (" %.3f", EGlpNumToLf (a[j]));
		}
		printf ("\n");
		fflush (stdout);
#endif
#endif /* SOLVE_DEBUG */
	}
#ifdef SOLVE_DEBUG
#if (SOLVE_DEBUG+0 <= 1)
	printf ("ILLfactor_ftran a after eta:");
	for (j = 0; j < f->dim; j++)
	{
		printf (" %.3f", EGlpNumToLf (a[j]));
	}
	printf ("\n");
	fflush (stdout);
#endif
#endif /* SOLVE_DEBUG */
	EGlpNumClearVar (v);
}

static void ILLfactor_ftrane2 (
	factor_work * f,
	svector * a)
{
	int *erindx = f->erindx;
	EGlpNum_t *ercoef = f->ercoef;
	er_info *er_inf = f->er_inf;
	int etacnt = f->etacnt;
	int beg;
	int nzcnt;
	int anzcnt = a->nzcnt;
	int *aindx = a->indx;
	EGlpNum_t *acoef = a->coef;
	EGlpNum_t *work_coef = f->work_coef;
	int *work_indx = f->work_indx;
	int i;
	int j;
	int r;
	EGlpNum_t v;

	EGlpNumInitVar (v);

	for (i = 0; i < anzcnt; i++)
	{
		EGlpNumCopy (work_coef[aindx[i]], acoef[i]);
		work_indx[aindx[i]] = i + 1;
	}
	for (i = 0; i < etacnt; i++)
	{
		r = er_inf[i].r;
		EGlpNumCopy (v, work_coef[r]);
		nzcnt = er_inf[i].nzcnt;
		beg = er_inf[i].rbeg;
		for (j = 0; j < nzcnt; j++)
		{
			EGlpNumSubInnProdTo (v, ercoef[beg + j], work_coef[erindx[beg + j]]);
		}
		if (EGlpNumIsNeqqZero (v))
		{
			EGlpNumCopy (work_coef[r], v);
			if (work_indx[r] == 0)
			{
				EGlpNumCopy (acoef[anzcnt], v);
				aindx[anzcnt] = r;
				work_indx[r] = anzcnt + 1;
				anzcnt++;
			}
			else
			{
				EGlpNumCopy (acoef[work_indx[r] - 1], v);
			}
		}
		else
		{
			EGlpNumZero (work_coef[r]);
			if (work_indx[r])
			{
				EGlpNumZero (acoef[work_indx[r] - 1]);
			}
		}
#ifdef SOLVE_DEBUG
#if (SOLVE_DEBUG+0 > 1)
		printf ("ILLfactor_ftran a after eta2 %d:", i);
		for (j = 0; j < anzcnt; j++)
		{
			printf (" %.3f*%d", EGlpNumToLf (acoef[j]), aindx[j]);
		}
		printf ("\n");
		fflush (stdout);
#endif
#endif /* SOLVE_DEBUG */
	}
	i = 0;
	while (i < anzcnt)
	{
		EGlpNumZero (work_coef[aindx[i]]);
		work_indx[aindx[i]] = 0;
		if (EGlpNumIsNeqZero (acoef[i], f->fzero_tol))
		{
			/*if (acoef[i] > fzero_tol || acoef[i] < -fzero_tol) */
			i++;
		}
		else
		{
			--anzcnt;
			EGlpNumCopy (acoef[i], acoef[anzcnt]);
			aindx[i] = aindx[anzcnt];
		}
	}
	a->nzcnt = anzcnt;

#ifdef SOLVE_DEBUG
#if (SOLVE_DEBUG+0 <= 1)
	printf ("ILLfactor_ftran a after eta2:");
	for (j = 0; j < anzcnt; j++)
	{
		printf (" %.3f*%d", EGlpNumToLf (acoef[j]), aindx[j]);
	}
	printf ("\n");
	fflush (stdout);
#endif
#endif /* SOLVE_DEBUG */
	EGlpNumClearVar (v);
}

static void ILLfactor_ftranu (
	factor_work * f,
	EGlpNum_t * a,
	svector * x)
{
	int *ucindx = f->ucindx;
	EGlpNum_t *uccoef = f->uccoef;
	uc_info *uc_inf = f->uc_inf;
	int *cperm = f->cperm;
	int *rperm = f->rperm;
	int dim = f->dim;
	int xnzcnt = 0;
	int *xindx = x->indx;
	EGlpNum_t *xcoef = x->coef;
	int nzcnt;
	int beg;
	int i;
	int j;
	EGlpNum_t v;

	EGlpNumInitVar (v);

	for (i = dim - 1; i >= 0; i--)
	{
		EGlpNumCopy (v, a[rperm[i]]);
		if (EGlpNumIsNeqqZero (v))	/*((v = a[rperm[i]]) != 0.0) */
		{
			j = cperm[i];
			beg = uc_inf[j].cbeg;
			EGlpNumDivTo (v, uccoef[beg]);
			if (EGlpNumIsNeqZero (v, f->szero_tol))
			{
				/*if (v > szero_tol || v < -szero_tol) */
				xindx[xnzcnt] = j;
				EGlpNumCopy (xcoef[xnzcnt], v);
				xnzcnt++;
			}
			nzcnt = uc_inf[j].nzcnt;
			for (j = 1; j < nzcnt; j++)
			{
				EGlpNumSubInnProdTo (a[ucindx[beg + j]], v, uccoef[beg + j]);
			}
			EGlpNumZero (a[rperm[i]]);
		}
	}
	x->nzcnt = xnzcnt;
#ifdef SOLVE_DEBUG
	printf ("ILLfactor_ftran x after u:");
	for (j = 0; j < x->nzcnt; j++)
	{
		printf (" %.3f*%d", EGlpNumToLf (x->coef[j]), x->indx[j]);
	}
	printf ("\n");
	fflush (stdout);
#endif /* SOLVE_DEBUG */
	EGlpNumClearVar (v);
}


#if 0
static void ftranu3_delay (
	factor_work * f,
	int c)
{
	uc_info *uc_inf = f->uc_inf;
	int nzcnt;
	int *indx;
	int i;

	c = f->cperm[f->rrank[c]];
	nzcnt = uc_inf[c].nzcnt;
	indx = f->ucindx + uc_inf[c].cbeg;
	for (i = 1; i < nzcnt; i++)
	{
		c = indx[i];
		if (uc_inf[c].delay++ == 0)
		{
			ftranu3_delay (f, c);
		}
	}
}
#endif

static void ftranu3_delay2 (
	factor_work * f,
	int c)
{
	uc_info *uc_inf = f->uc_inf;
	int nzcnt;
	int *indx;
	int i;
	int last;

	do
	{
		c = f->cperm[f->rrank[c]];
		nzcnt = uc_inf[c].nzcnt;
		indx = f->ucindx + uc_inf[c].cbeg;
		last = -1;
		for (i = 1; i < nzcnt; i++)
		{
			c = indx[i];
			if (uc_inf[c].delay++ == 0)
			{
				if (last >= 0)
				{
					ftranu3_delay2 (f, last);
				}
				last = c;
			}
		}
		c = last;
	} while (c >= 0);
}

#if 0
static void ftranu3_process (
	factor_work * f,
	int c,
	svector * x)
{
	uc_info *uc_inf = f->uc_inf;
	EGlpNum_t *work = f->work_coef;
	int nzcnt;
	int *indx;
	EGlpNum_t *coef;
	int i;
	EGlpNum_t v;

	EGlpNumInitVar (v);

	EGlpNumCopy (v, work[c]);
	EGlpNumZero (work[c]);
	c = f->cperm[f->rrank[c]];
	nzcnt = uc_inf[c].nzcnt;
	indx = f->ucindx + uc_inf[c].cbeg;
	coef = f->uccoef + uc_inf[c].cbeg;
	EGlpNumDivTo (v, coef[0]);
	if (EGlpNumIsNeqZero (v, f->szero_tol))
		/*if (v > szero_tol || v < -szero_tol) */
	{
		x->indx[x->nzcnt] = c;
		EGlpNumCopy (x->coef[x->nzcnt], v);
		x->nzcnt++;
	}
	for (i = 1; i < nzcnt; i++)
	{
		c = indx[i];
		EGlpNumSubInnProdTo (work[c], v, coef[i]);
		if (--uc_inf[c].delay == 0)
		{
			ftranu3_process (f, c, x);
		}
	}
	EGlpNumClearVar (v);
}
#endif

static void ftranu3_process2 (
	factor_work * f,
	int c,
	svector * x)
{
	uc_info *uc_inf = f->uc_inf;
	EGlpNum_t *work = f->work_coef;
	int nzcnt;
	int *indx;
	EGlpNum_t *coef;
	int i;
	int last;
	EGlpNum_t v;

	EGlpNumInitVar (v);

	do
	{
		EGlpNumCopy (v, work[c]);
		EGlpNumZero (work[c]);
		c = f->cperm[f->rrank[c]];
		nzcnt = uc_inf[c].nzcnt;
		indx = f->ucindx + uc_inf[c].cbeg;
		coef = f->uccoef + uc_inf[c].cbeg;
		EGlpNumDivTo (v, coef[0]);
		if (EGlpNumIsNeqZero (v, f->szero_tol))
			/*if (v > szero_tol || v < -szero_tol) */
		{
			x->indx[x->nzcnt] = c;
			EGlpNumCopy (x->coef[x->nzcnt], v);
			x->nzcnt++;
		}
		last = -1;
		for (i = 1; i < nzcnt; i++)
		{
			c = indx[i];
			EGlpNumSubInnProdTo (work[c], v, coef[i]);
			if (--uc_inf[c].delay == 0)
			{
				if (last >= 0)
				{
					ftranu3_process2 (f, last, x);
				}
				last = c;
			}
		}
		c = last;
	} while (c >= 0);
	EGlpNumClearVar (v);
}

static void ILLfactor_ftranu3 (
	factor_work * f,
	svector * a,
	svector * x)
{
	EGlpNum_t *work = f->work_coef;
	int anzcnt = a->nzcnt;
	int *aindx = a->indx;
	EGlpNum_t *acoef = a->coef;
	uc_info *uc_inf = f->uc_inf;
	int i;

	for (i = 0; i < anzcnt; i++)
	{
		if (uc_inf[aindx[i]].delay++ == 0)
		{
			ftranu3_delay2 (f, aindx[i]);
		}
		EGlpNumCopy (work[aindx[i]], acoef[i]);
	}
	x->nzcnt = 0;
	for (i = 0; i < anzcnt; i++)
	{
		if (--uc_inf[aindx[i]].delay == 0)
		{
			ftranu3_process2 (f, aindx[i], x);
		}
	}
#ifdef SOLVE_DEBUG
	printf ("ILLfactor_ftran x after u3:");
	for (i = 0; i < x->nzcnt; i++)
	{
		printf (" %.3f*%d", EGlpNumToLf (x->coef[i]), x->indx[i]);
	}
	printf ("\n");
	fflush (stdout);
#endif /* SOLVE_DEBUG */
}

/* ILLfactor_ftran solves Bx=a for x */
void ILLfactor_ftran (
	factor_work * f,
	svector * a,
	svector * x)
{
	int i;
	int nzcnt;
	int sparse;
	int *aindx;
	EGlpNum_t *acoef;
	EGlpNum_t *work_coef = f->work_coef;

#ifdef RECORD
	{
		EGioPrintf (fsave, "f %d", a->nzcnt);
		for (i = 0; i < a->nzcnt; i++)
		{
			EGioPrintf (fsave, " %d %.16e", a->indx[i], EGlpNumToLf (a->coef[i]));
		}
		EGioPrintf (fsave, "\n");
		EGioFlush (fsave);
	}
#endif /* RECORD */
#ifdef DEBUG_FACTOR
	{
		printf ("ILLfactor_ftran a:");
		for (i = 0; i < a->nzcnt; i++)
		{
			printf (" %d %la", a->indx[i], EGlpNumToLf (a->coef[i]));
		}
		printf ("\n");
		fflush (stdout);
	}
#endif /* DEBUG_FACTOR */

	if (a->nzcnt >= SPARSE_FACTOR * f->dim)
	{
		nzcnt = a->nzcnt;
		aindx = a->indx;
		acoef = a->coef;
		for (i = 0; i < nzcnt; i++)
		{
			EGlpNumCopy (work_coef[aindx[i]], acoef[i]);
		}
		sparse = 0;
	}
	else
	{
		sparse = 1;
	}

	if (sparse)
	{
		ILLfactor_ftranl3 (f, a, &f->xtmp);
		if (f->xtmp.nzcnt >= SPARSE_FACTOR * f->dim)
		{
			nzcnt = f->xtmp.nzcnt;
			aindx = f->xtmp.indx;
			acoef = f->xtmp.coef;

			for (i = 0; i < nzcnt; i++)
			{
				EGlpNumCopy (work_coef[aindx[i]], acoef[i]);
			}
			sparse = 0;
		}
	}
	else
	{
		ILLfactor_ftranl (f, work_coef);
	}

	if (sparse)
	{
		ILLfactor_ftrane2 (f, &f->xtmp);
		if (f->xtmp.nzcnt >= SPARSE_FACTOR * f->dim)
		{
			nzcnt = f->xtmp.nzcnt;
			aindx = f->xtmp.indx;
			acoef = f->xtmp.coef;

			for (i = 0; i < nzcnt; i++)
			{
				EGlpNumCopy (work_coef[aindx[i]], acoef[i]);
			}
			sparse = 0;
		}
	}
	else
	{
		ILLfactor_ftrane (f, work_coef);
	}

	if (sparse)
	{
		ILLfactor_ftranu3 (f, &f->xtmp, x);
	}
	else
	{
		ILLfactor_ftranu (f, work_coef, x);
	}

#ifdef SORT_RESULTS
	sort_vector (x);
#endif

#ifdef DEBUG_FACTOR
	{
		printf ("ILLfactor_ftran x:");
		for (i = 0; i < x->nzcnt; i++)
		{
			printf (" %d %la", x->indx[i], EGlpNumToLf (x->coef[i]));
		}
		printf ("\n");
		fflush (stdout);
	}
#endif /* DEBUG_FACTOR */
	return;
}

/* ILLfactor_ftran_update solves Bx=a for x, and also returns upd, where Ux=upd */
void ILLfactor_ftran_update (
	factor_work * f,
	svector * a,
	svector * upd,
	svector * x)
{
	int i;
	int nzcnt;
	int dim;
	int sparse;
	int *aindx;
	EGlpNum_t *acoef;
	EGlpNum_t *work_coef = f->work_coef;

#ifdef RECORD
	{
		EGioPrintf (fsave, "F %d", a->nzcnt);
		for (i = 0; i < a->nzcnt; i++)
		{
			EGioPrintf (fsave, " %d %.16e", a->indx[i], EGlpNumToLf (a->coef[i]));
		}
		EGioPrintf (fsave, "\n");
		EGioFlush (fsave);
	}
#endif /* RECORD */
#ifdef DEBUG_FACTOR
	{
		printf ("ILLfactor_ftran_update a:");
		for (i = 0; i < a->nzcnt; i++)
		{
			printf (" %d %.3f", a->indx[i], EGlpNumToLf (a->coef[i]));
		}
		printf ("\n");
		fflush (stdout);
	}
#endif /* DEBUG_FACTOR */

	if (a->nzcnt >= SPARSE_FACTOR * f->dim)
	{
		aindx = a->indx;
		acoef = a->coef;
		nzcnt = a->nzcnt;

		for (i = 0; i < nzcnt; i++)
		{
			EGlpNumCopy (work_coef[aindx[i]], acoef[i]);
		}
		sparse = 0;
	}
	else
	{
		sparse = 1;
	}

	if (sparse)
	{
		ILLfactor_ftranl3 (f, a, upd);
		if (upd->nzcnt >= SPARSE_FACTOR * f->dim)
		{
			nzcnt = upd->nzcnt;
			aindx = upd->indx;
			acoef = upd->coef;

			for (i = 0; i < nzcnt; i++)
			{
				EGlpNumCopy (work_coef[aindx[i]], acoef[i]);
			}
			sparse = 0;
		}
	}
	else
	{
		ILLfactor_ftranl (f, work_coef);
	}

	if (sparse)
	{
		ILLfactor_ftrane2 (f, upd);
		if (upd->nzcnt >= SPARSE_FACTOR * f->dim)
		{
			nzcnt = upd->nzcnt;
			aindx = upd->indx;
			acoef = upd->coef;

			for (i = 0; i < nzcnt; i++)
			{
				EGlpNumCopy (work_coef[aindx[i]], acoef[i]);
			}
			sparse = 0;
		}
	}
	else
	{
		ILLfactor_ftrane (f, work_coef);
		nzcnt = 0;
		dim = f->dim;
		aindx = upd->indx;
		acoef = upd->coef;
		for (i = 0; i < dim; i++)
		{
			if (EGlpNumIsNeqqZero (work_coef[i]))
			{
				if (EGlpNumIsNeqZero (work_coef[i], f->szero_tol))
					/*if(work_coef[i] > szero_tol || work_coef[i] < -szero_tol) */
				{
					aindx[nzcnt] = i;
					EGlpNumCopy (acoef[nzcnt], work_coef[i]);
					nzcnt++;
				}
			}
		}
		upd->nzcnt = nzcnt;
	}

	if (sparse)
	{
		ILLfactor_ftranu3 (f, upd, x);
	}
	else
	{
		ILLfactor_ftranu (f, work_coef, x);
	}

#ifdef SORT_RESULTS
	sort_vector (upd);
	sort_vector (x);
#endif

#ifdef DEBUG_FACTOR
	{
		printf ("ILLfactor_ftran update x:");
		for (i = 0; i < x->nzcnt; i++)
		{
			printf (" %d %.3f", x->indx[i], EGlpNumToLf (x->coef[i]));
		}
		printf ("\n");
		fflush (stdout);
	}
#endif /* DEBUG_FACTOR */
}


static void ILLfactor_btranl2 (
	factor_work * f,
	EGlpNum_t * x)
{
	int *lrindx = f->lrindx;
	EGlpNum_t *lrcoef = f->lrcoef;
	lr_info *lr_inf = f->lr_inf;
	int dim = f->dim;
	int nzcnt;
	int beg;
	int i;
	int j;
	EGlpNum_t v;

	EGlpNumInitVar (v);

	for (i = dim - 1; i >= 0; i--)
	{
#ifdef SOLVE_DEBUG
#if (SOLVE_DEBUG+0 > 1)
		printf ("ILLfactor_btran x before l2 %d:", i);
		for (j = 0; j < f->dim; j++)
		{
			printf (" %.3f", EGlpNumToLf (x[j]));
		}
		printf ("\n");
		fflush (stdout);
#endif
#endif /* SOLVE_DEBUG */
		EGlpNumCopy (v, x[lr_inf[i].r]);
		if (EGlpNumIsNeqqZero (v))
		{
			nzcnt = lr_inf[i].nzcnt;
			beg = lr_inf[i].rbeg;
			for (j = 0; j < nzcnt; j++)
			{
				EGlpNumSubInnProdTo (x[lrindx[beg + j]], v, lrcoef[beg + j]);
			}
		}
	}
#ifdef SOLVE_DEBUG
#if (SOLVE_DEBUG+0 <= 1)
	printf ("ILLfactor_btran x after l2:");
	for (j = 0; j < f->dim; j++)
	{
		printf (" %.3f", EGlpNumToLf (x[j]));
	}
	printf ("\n");
	fflush (stdout);
#endif
#endif /* SOLVE_DEBUG */
	EGlpNumClearVar (v);
}

#if 0
static void btranl3_delay (
	factor_work * f,
	int r)
{
	lr_info *lr_inf = f->lr_inf;
	int nzcnt;
	int *indx;
	int i;

	r = lr_inf[r].rrank;
	nzcnt = lr_inf[r].nzcnt;
	indx = f->lrindx + lr_inf[r].rbeg;
	for (i = 0; i < nzcnt; i++)
	{
		r = indx[i];
		if (lr_inf[r].delay++ == 0)
		{
			btranl3_delay (f, r);
		}
	}
}
#endif

static void btranl3_delay2 (
	factor_work * f,
	int r)
{
	lr_info *lr_inf = f->lr_inf;
	int nzcnt;
	int *indx;
	int i;
	int last;

	do
	{
		r = lr_inf[r].rrank;
		nzcnt = lr_inf[r].nzcnt;
		indx = f->lrindx + lr_inf[r].rbeg;
		last = -1;
		for (i = 0; i < nzcnt; i++)
		{
			r = indx[i];
			if (lr_inf[r].delay++ == 0)
			{
				if (last >= 0)
				{
					btranl3_delay2 (f, last);
				}
				last = r;
			}
		}
		r = last;
	} while (r >= 0);
}

#if 0
static void btranl3_process (
	factor_work * f,
	int r,
	svector * x)
{
	lr_info *lr_inf = f->lr_inf;
	EGlpNum_t *work = f->work_coef;
	int nzcnt;
	int *indx;
	EGlpNum_t *coef;
	int i;
	EGlpNum_t v;

	EGlpNumInitVar (v);

	EGlpNumCopy (v, work[r]);
	EGlpNumZero (work[r]);
	if (EGlpNumIsNeqZero (v, f->szero_tol))
		/*if (v > szero_tol || v < -szero_tol) */
	{
		x->indx[x->nzcnt] = r;
		EGlpNumCopy (x->coef[x->nzcnt], v);
		x->nzcnt++;
	}
	r = lr_inf[r].rrank;
	nzcnt = lr_inf[r].nzcnt;
	indx = f->lrindx + lr_inf[r].rbeg;
	coef = f->lrcoef + lr_inf[r].rbeg;
	for (i = 0; i < nzcnt; i++)
	{
		r = indx[i];
		EGlpNumSubInnProdTo (work[r], v, coef[i]);
		if (--lr_inf[r].delay == 0)
		{
			btranl3_process (f, r, x);
		}
	}
	EGlpNumClearVar (v);
}
#endif

static void btranl3_process2 (
	factor_work * f,
	int r,
	svector * x)
{
	lr_info *lr_inf = f->lr_inf;
	EGlpNum_t *work = f->work_coef;
	int nzcnt;
	int *indx;
	EGlpNum_t *coef;
	int i;
	int last;
	EGlpNum_t v;

	EGlpNumInitVar (v);

	do
	{
		EGlpNumCopy (v, work[r]);
		EGlpNumZero (work[r]);
		if (EGlpNumIsNeqZero (v, f->szero_tol))
			/*if (v > szero_tol || v < -szero_tol) */
		{
			x->indx[x->nzcnt] = r;
			EGlpNumCopy (x->coef[x->nzcnt], v);
			x->nzcnt++;
		}
		r = lr_inf[r].rrank;
		nzcnt = lr_inf[r].nzcnt;
		indx = f->lrindx + lr_inf[r].rbeg;
		coef = f->lrcoef + lr_inf[r].rbeg;
		last = -1;
		for (i = 0; i < nzcnt; i++)
		{
			r = indx[i];
			EGlpNumSubInnProdTo (work[r], v, coef[i]);
			if (--lr_inf[r].delay == 0)
			{
				if (last >= 0)
				{
					btranl3_process2 (f, last, x);
				}
				last = r;
			}
		}
		r = last;
	} while (r >= 0);
	EGlpNumClearVar (v);
}

static void ILLfactor_btranl3 (
	factor_work * f,
	svector * a,
	svector * x)
{
	EGlpNum_t *work = f->work_coef;
	int anzcnt = a->nzcnt;
	int *aindx = a->indx;
	EGlpNum_t *acoef = a->coef;
	lr_info *lr_inf = f->lr_inf;
	int i;

	for (i = 0; i < anzcnt; i++)
	{
		if (lr_inf[aindx[i]].delay++ == 0)
		{
			btranl3_delay2 (f, aindx[i]);
		}
		EGlpNumCopy (work[aindx[i]], acoef[i]);
	}
	x->nzcnt = 0;
	for (i = 0; i < anzcnt; i++)
	{
		if (--lr_inf[aindx[i]].delay == 0)
		{
			btranl3_process2 (f, aindx[i], x);
		}
	}
#ifdef SOLVE_DEBUG
	printf ("ILLfactor_btran x after l3:");
	for (i = 0; i < x->nzcnt; i++)
	{
		printf (" %.3f*%d", EGlpNumToLf (x->coef[i]), x->indx[i]);
	}
	printf ("\n");
	fflush (stdout);
#endif /* SOLVE_DEBUG */
}

static void ILLfactor_btrane (
	factor_work * f,
	EGlpNum_t * x)
{
	int *erindx = f->erindx;
	EGlpNum_t *ercoef = f->ercoef;
	er_info *er_inf = f->er_inf;
	int etacnt = f->etacnt;
	int beg;
	int nzcnt;
	int i;
	int j;
	EGlpNum_t v;

	EGlpNumInitVar (v);

	for (i = etacnt - 1; i >= 0; i--)
	{
#ifdef SOLVE_DEBUG
#if (SOLVE_DEBUG+0 > 1)
		printf ("ILLfactor_btran x before eta %d:", i);
		for (j = 0; j < f->dim; j++)
		{
			printf (" %.3f", EGlpNumToLf (x[j]));
		}
		printf ("\n");
		fflush (stdout);
#endif
#endif /* SOLVE_DEBUG */
		EGlpNumCopy (v, x[er_inf[i].r]);
		if (EGlpNumIsNeqqZero (v))
		{
			nzcnt = er_inf[i].nzcnt;
			beg = er_inf[i].rbeg;
			for (j = 0; j < nzcnt; j++)
			{
				EGlpNumSubInnProdTo (x[erindx[beg + j]], v, ercoef[beg + j]);
			}
		}
	}
#ifdef SOLVE_DEBUG
#if (SOLVE_DEBUG+0 <= 1)
	printf ("ILLfactor_btran x after eta:");
	for (j = 0; j < f->dim; j++)
	{
		printf (" %.3f", EGlpNumToLf (x[j]));
	}
	printf ("\n");
	fflush (stdout);
#endif
#endif /* SOLVE_DEBUG */
	EGlpNumClearVar (v);
}

static void ILLfactor_btrane2 (
	factor_work * f,
	svector * x)
{
	int *erindx = f->erindx;
	EGlpNum_t *ercoef = f->ercoef;
	er_info *er_inf = f->er_inf;
	int etacnt = f->etacnt;
	int beg;
	int nzcnt;
	int xnzcnt = x->nzcnt;
	int *xindx = x->indx;
	EGlpNum_t *xcoef = x->coef;
	EGlpNum_t *work_coef = f->work_coef;
	int *work_indx = f->work_indx;
	int i;
	int j;
	EGlpNum_t v;

	EGlpNumInitVar (v);

	for (i = 0; i < xnzcnt; i++)
	{
		EGlpNumCopy (work_coef[xindx[i]], xcoef[i]);
		work_indx[xindx[i]] = i + 1;
	}
	for (i = etacnt - 1; i >= 0; i--)
	{
#ifdef SOLVE_DEBUG
#if (SOLVE_DEBUG+0 > 1)
		printf ("ILLfactor_btran x before eta2 %d:", i);
		for (j = 0; j < xnzcnt; j++)
		{
			printf (" %.3f*%d", EGlpNumToLf (work_coef[xindx[j]]), xindx[j]);
		}
		printf ("\n");
		fflush (stdout);
#endif
#endif /* SOLVE_DEBUG */
		EGlpNumCopy (v, work_coef[er_inf[i].r]);
		if (EGlpNumIsNeqqZero (v))
		{
			nzcnt = er_inf[i].nzcnt;
			beg = er_inf[i].rbeg;
			for (j = 0; j < nzcnt; j++)
			{
				if (work_indx[erindx[beg + j]] == 0)
				{
					work_indx[erindx[beg + j]] = xnzcnt;
					xindx[xnzcnt++] = erindx[beg + j];
				}
				EGlpNumSubInnProdTo (work_coef[erindx[beg + j]], v, ercoef[beg + j]);
			}
		}
	}

	j = 0;
	while (j < xnzcnt)
	{
		EGlpNumCopy (xcoef[j], work_coef[xindx[j]]);
		EGlpNumZero (work_coef[xindx[j]]);
		work_indx[xindx[j]] = 0;
		if (!EGlpNumIsNeqqZero (xcoef[j]))
		{
			--xnzcnt;
			xindx[j] = xindx[xnzcnt];
		}
		else
		{
			j++;
		}
	}
	x->nzcnt = xnzcnt;

#ifdef SOLVE_DEBUG
#if (SOLVE_DEBUG+0 <= 1)
	printf ("ILLfactor_btran x after eta2:");
	for (j = 0; j < xnzcnt; j++)
	{
		printf (" %.3f*%d", EGlpNumToLf (xcoef[j]), xindx[j]);
	}
	printf ("\n");
	fflush (stdout);
#endif
#endif /* SOLVE_DEBUG */
	EGlpNumClearVar (v);
}

static void ILLfactor_btranu (
	factor_work * f,
	EGlpNum_t * a,
	svector * x)
{
	int *urindx = f->urindx;
	EGlpNum_t *urcoef = f->urcoef;
	ur_info *ur_inf = f->ur_inf;
	int *rperm = f->rperm;
	int *cperm = f->cperm;
	int dim = f->dim;
	int xnzcnt = 0;
	int *xindx = x->indx;
	EGlpNum_t *xcoef = x->coef;
	int nzcnt;
	int beg;
	int i;
	int j;
	EGlpNum_t v;

	EGlpNumInitVar (v);

	for (i = 0; i < dim; i++)
	{
		EGlpNumCopy (v, a[cperm[i]]);
		if (EGlpNumIsNeqqZero (v))
		{
			j = rperm[i];
			beg = ur_inf[j].rbeg;
			EGlpNumDivTo (v, urcoef[beg]);
			if (EGlpNumIsNeqZero (v, f->szero_tol))	/*
																							 * if (v > szero_tol || v < -szero_tol) */
			{
				xindx[xnzcnt] = j;
				EGlpNumCopy (xcoef[xnzcnt], v);
				xnzcnt++;
			}
			nzcnt = ur_inf[j].nzcnt;
			for (j = 1; j < nzcnt; j++)
			{
				EGlpNumSubInnProdTo (a[urindx[beg + j]], v, urcoef[beg + j]);
			}
			EGlpNumZero (a[cperm[i]]);
		}
	}
	x->nzcnt = xnzcnt;
#ifdef SOLVE_DEBUG
	printf ("ILLfactor_btran x after u:");
	for (i = 0; i < x->nzcnt; i++)
	{
		printf (" %.3f*%d", EGlpNumToLf (x->coef[i]), x->indx[i]);
	}
	printf ("\n");
	fflush (stdout);
#endif /* SOLVE_DEBUG */
	EGlpNumClearVar (v);
}


#if 0
static void btranu3_delay (
	factor_work * f,
	int r)
{
	ur_info *ur_inf = f->ur_inf;
	int nzcnt;
	int *indx;
	int i;

	r = f->rperm[f->crank[r]];
	nzcnt = ur_inf[r].nzcnt;
	indx = f->urindx + ur_inf[r].rbeg;
	for (i = 1; i < nzcnt; i++)
	{
		r = indx[i];
		if (ur_inf[r].delay++ == 0)
		{
			btranu3_delay (f, r);
		}
	}
}
#endif

static void btranu3_delay2 (
	factor_work * f,
	int r)
{
	ur_info *ur_inf = f->ur_inf;
	int nzcnt;
	int *indx;
	int i;
	int last;

	do
	{
		r = f->rperm[f->crank[r]];
		nzcnt = ur_inf[r].nzcnt;
		indx = f->urindx + ur_inf[r].rbeg;
		last = -1;
		for (i = 1; i < nzcnt; i++)
		{
			r = indx[i];
			if (ur_inf[r].delay++ == 0)
			{
				if (last >= 0)
				{
					btranu3_delay2 (f, last);
				}
				last = r;
			}
		}
		r = last;
	} while (r >= 0);
}

#if 0
static void btranu3_process (
	factor_work * f,
	int r,
	svector * x)
{
	ur_info *ur_inf = f->ur_inf;
	EGlpNum_t *work = f->work_coef;
	int nzcnt;
	int *indx;
	EGlpNum_t *coef;
	int i;
	EGlpNum_t v;

	EGlpNumInitVar (v);

	EGlpNumCopy (v, work[r]);
	EGlpNumZero (work[r]);
	r = f->rperm[f->crank[r]];
	nzcnt = ur_inf[r].nzcnt;
	indx = f->urindx + ur_inf[r].rbeg;
	coef = f->urcoef + ur_inf[r].rbeg;
	EGlpNumDivTo (v, coef[0]);
	if (EGlpNumIsNeqqZero (v))
	{
		x->indx[x->nzcnt] = r;
		EGlpNumCopy (x->coef[x->nzcnt], v);
		x->nzcnt++;
	}
	for (i = 1; i < nzcnt; i++)
	{
		r = indx[i];
		EGlpNumSubInnProdTo (work[r], v, coef[i]);
		if (--ur_inf[r].delay == 0)
		{
			btranu3_process (f, r, x);
		}
	}
	EGlpNumClearVar (v);
}
#endif

static void btranu3_process2 (
	factor_work * f,
	int r,
	svector * x)
{
	ur_info *ur_inf = f->ur_inf;
	EGlpNum_t *work = f->work_coef;
	int nzcnt;
	int *indx;
	EGlpNum_t *coef;
	int i;
	int last;
	EGlpNum_t v;

	EGlpNumInitVar (v);

	do
	{
		EGlpNumCopy (v, work[r]);
		EGlpNumZero (work[r]);
		r = f->rperm[f->crank[r]];
		nzcnt = ur_inf[r].nzcnt;
		indx = f->urindx + ur_inf[r].rbeg;
		coef = f->urcoef + ur_inf[r].rbeg;
		EGlpNumDivTo (v, coef[0]);
		if (EGlpNumIsNeqqZero (v))
		{
			x->indx[x->nzcnt] = r;
			EGlpNumCopy (x->coef[x->nzcnt], v);
			x->nzcnt++;
		}
		last = -1;
		for (i = 1; i < nzcnt; i++)
		{
			r = indx[i];
			EGlpNumSubInnProdTo (work[r], v, coef[i]);
			if (--ur_inf[r].delay == 0)
			{
				if (last >= 0)
				{
					btranu3_process2 (f, last, x);
				}
				last = r;
			}
		}
		r = last;
	} while (r >= 0);
	EGlpNumClearVar (v);
}

static void ILLfactor_btranu3 (
	factor_work * f,
	svector * a,
	svector * x)
{
	EGlpNum_t *work = f->work_coef;
	int anzcnt = a->nzcnt;
	int *aindx = a->indx;
	EGlpNum_t *acoef = a->coef;
	ur_info *ur_inf = f->ur_inf;
	int i;

	for (i = 0; i < anzcnt; i++)
	{
		if (ur_inf[aindx[i]].delay++ == 0)
		{
			btranu3_delay2 (f, aindx[i]);
		}
		EGlpNumCopy (work[aindx[i]], acoef[i]);
	}
	x->nzcnt = 0;
	for (i = 0; i < anzcnt; i++)
	{
		if (--ur_inf[aindx[i]].delay == 0)
		{
			btranu3_process2 (f, aindx[i], x);
		}
	}
#ifdef SOLVE_DEBUG
	printf ("ILLfactor_btran x after u3:");
	for (i = 0; i < x->nzcnt; i++)
	{
		printf (" %.3f*%d", EGlpNumToLf (x->coef[i]), x->indx[i]);
	}
	printf ("\n");
	fflush (stdout);
#endif /* SOLVE_DEBUG */
}

/* ILLfactor_btran solves x^tB=a^t (or, B^t x = a) for x */
void ILLfactor_btran (
	factor_work * f,
	svector * a,
	svector * x)
{
	int i;
	int nzcnt;
	int sparse;
	int *aindx = a->indx;
	EGlpNum_t *acoef = a->coef;
	EGlpNum_t *work_coef = f->work_coef;
	int dim = f->dim;

#ifdef RECORD
	{
		EGioPrintf (fsave, "b %d", a->nzcnt);
		for (i = 0; i < a->nzcnt; i++)
		{
			EGioPrintf (fsave, " %d %.16e", a->indx[i], EGlpNumToLf (a->coef[i]));
		}
		EGioPrintf (fsave, "\n");
		EGioFlush (fsave);
	}
#endif /* RECORD */
#ifdef DEBUG_FACTOR
	{
		printf ("ILLfactor_btran a:");
		for (i = 0; i < a->nzcnt; i++)
		{
			printf (" %d %.3f", a->indx[i], EGlpNumToLf (a->coef[i]));
		}
		printf ("\n");
		fflush (stdout);
	}
#endif /* DEBUG_FACTOR */

	if (a->nzcnt >= SPARSE_FACTOR * f->dim)
	{
		aindx = a->indx;
		acoef = a->coef;
		work_coef = f->work_coef;
		nzcnt = a->nzcnt;
		for (i = 0; i < nzcnt; i++)
		{
			EGlpNumCopy (work_coef[aindx[i]], acoef[i]);
		}
		sparse = 0;
	}
	else
	{
		sparse = 1;
	}

	if (sparse)
	{
		ILLfactor_btranu3 (f, a, &f->xtmp);
	}
	else
	{
		ILLfactor_btranu (f, work_coef, &f->xtmp);
	}

	if (f->xtmp.nzcnt >= SPARSE_FACTOR * f->dim)
	{
		aindx = f->xtmp.indx;
		acoef = f->xtmp.coef;
		work_coef = f->work_coef;
		nzcnt = f->xtmp.nzcnt;
		for (i = 0; i < nzcnt; i++)
		{
			EGlpNumCopy (work_coef[aindx[i]], acoef[i]);
		}
		sparse = 0;
	}
	else
	{
		sparse = 1;
	}

	if (sparse)
	{
		ILLfactor_btrane2 (f, &f->xtmp);
		if (f->xtmp.nzcnt >= SPARSE_FACTOR * f->dim)
		{
			aindx = f->xtmp.indx;
			acoef = f->xtmp.coef;
			work_coef = f->work_coef;
			nzcnt = f->xtmp.nzcnt;
			for (i = 0; i < nzcnt; i++)
			{
				EGlpNumCopy (work_coef[aindx[i]], acoef[i]);
			}
			sparse = 0;
		}
	}
	else
	{
		ILLfactor_btrane (f, work_coef);
	}

	if (sparse)
	{
		ILLfactor_btranl3 (f, &f->xtmp, x);
	}
	else
	{
		ILLfactor_btranl2 (f, work_coef);
		dim = f->dim;
		nzcnt = 0;
		aindx = x->indx;
		acoef = x->coef;
		for (i = 0; i < dim; i++)
		{
			if (EGlpNumIsNeqqZero (work_coef[i]))
			{
				if (EGlpNumIsNeqZero (work_coef[i], f->szero_tol))
					/*if (work_coef[i] > szero_tol || work_coef[i] < -szero_tol) */
				{
					aindx[nzcnt] = i;
					EGlpNumCopy (acoef[nzcnt], work_coef[i]);
					nzcnt++;
				}
				EGlpNumZero (work_coef[i]);
			}
		}
		x->nzcnt = nzcnt;
	}

#ifdef SORT_RESULTS
	sort_vector (x);
#endif

#ifdef DEBUG_FACTOR
	{
		printf ("ILLfactor_btran x:");
		for (i = 0; i < x->nzcnt; i++)
		{
			printf (" %d %.3f", x->indx[i], EGlpNumToLf (x->coef[i]));
		}
		printf ("\n");
		fflush (stdout);
	}
#endif /* DEBUG_FACTOR */
	return;
}

static int expand_col (
	factor_work * f,
	int col)
{
	uc_info *uc_inf = f->uc_inf + col;
	int uc_freebeg = f->uc_freebeg;
	int nzcnt = uc_inf->nzcnt;
	int cbeg;
	EGlpNum_t *uccoef;
	int *ucindx;
	int *ucrind;
	int i;
	int rval = 0;

	if (uc_freebeg + nzcnt + 1 >= f->uc_space)
	{
		rval = make_uc_space (f, nzcnt + 1);
		CHECKRVALG (rval, CLEANUP);
		uc_freebeg = f->uc_freebeg;
	}
	cbeg = uc_inf->cbeg;
	uccoef = f->uccoef;
	ucindx = f->ucindx;
	ucrind = f->ucrind;

	for (i = 0; i < nzcnt; i++)
	{
		EGlpNumCopy (uccoef[uc_freebeg + i], uccoef[cbeg + i]);
		ucindx[uc_freebeg + i] = ucindx[cbeg + i];
		ucrind[uc_freebeg + i] = ucrind[cbeg + i];
		ucindx[cbeg + i] = -1;
	}

	uc_inf->cbeg = uc_freebeg;
	f->uc_freebeg = uc_freebeg + nzcnt;
CLEANUP:
	EG_RETURN (rval);
}

static int expand_row (
	factor_work * f,
	int row)
{
	ur_info *ur_inf = f->ur_inf + row;
	int ur_freebeg = f->ur_freebeg;
	int nzcnt = ur_inf->nzcnt;
	int rbeg;
	EGlpNum_t *urcoef;
	int *urindx;
	int *urcind;
	int i;
	int rval = 0;

	if (ur_freebeg + nzcnt + 1 >= f->ur_space)
	{
		rval = make_ur_space (f, nzcnt + 1);
		CHECKRVALG (rval, CLEANUP);
		ur_freebeg = f->ur_freebeg;
	}
	rbeg = ur_inf->rbeg;
	urcoef = f->urcoef;
	urindx = f->urindx;
	urcind = f->urcind;

	for (i = 0; i < nzcnt; i++)
	{
		EGlpNumCopy (urcoef[ur_freebeg + i], urcoef[rbeg + i]);
		urindx[ur_freebeg + i] = urindx[rbeg + i];
		urcind[ur_freebeg + i] = urcind[rbeg + i];
		urindx[rbeg + i] = -1;
	}

	ur_inf->rbeg = ur_freebeg;
	f->ur_freebeg = ur_freebeg + nzcnt;
CLEANUP:
	EG_RETURN (rval);
}

static int add_nonzero (
	factor_work * f,
	int row,
	int col,
	EGlpNum_t val)
{
	ur_info *ur_inf = f->ur_inf + row;
	uc_info *uc_inf = f->uc_inf + col;
	int cnzcnt = uc_inf->nzcnt;
	int rnzcnt = ur_inf->nzcnt;
	int cloc = uc_inf->cbeg + cnzcnt;
	int rloc = ur_inf->rbeg + rnzcnt;
	int rval = 0;

	if (f->ucindx[cloc] != -1)
	{
		rval = expand_col (f, col);
		CHECKRVALG (rval, CLEANUP);
		cloc = uc_inf->cbeg + cnzcnt;
	}
	TESTG ((rval = (rloc < 0 || rloc > f->ur_space)), CLEANUP,
				 "rloc %d outside boundaries [0:%d]", rloc, f->ur_space);
	if (f->urindx[rloc] != -1)
	{
		rval = expand_row (f, row);
		CHECKRVALG (rval, CLEANUP);
		rloc = ur_inf->rbeg + rnzcnt;
	}
	f->ucindx[cloc] = row;
	EGlpNumCopy (f->uccoef[cloc], val);
	f->ucrind[cloc] = rnzcnt;
	f->urindx[rloc] = col;
	EGlpNumCopy (f->urcoef[rloc], val);
	f->urcind[rloc] = cnzcnt;

	if (cloc == f->uc_freebeg)
		f->uc_freebeg++;
	if (rloc == f->ur_freebeg)
		f->ur_freebeg++;

	uc_inf->nzcnt = cnzcnt + 1;
	ur_inf->nzcnt = rnzcnt + 1;
CLEANUP:
	EG_RETURN (rval);
}

static int delete_nonzero_row (
	factor_work * f,
	int row,
	int ind)
{
	ur_info *ur_inf = f->ur_inf;
	EGlpNum_t *urcoef = f->urcoef;
	int *urindx = f->urindx;
	int *urcind = f->urcind;
	int *ucrind = f->ucrind;
	int rbeg = ur_inf[row].rbeg;
	int nzcnt = ur_inf[row].nzcnt - 1;
	int cbeg, rval = 0;
	#ifdef DEBUG_FACTOR
	TESTG((rval=(nzcnt<0)),CLEANUP,"Deleting empty row %d ind %d!",row, ind);
	#endif

	if (ind != nzcnt)
	{
		EGlpNumCopy (urcoef[rbeg + ind], urcoef[rbeg + nzcnt]);
		urindx[rbeg + ind] = urindx[rbeg + nzcnt];
		urcind[rbeg + ind] = urcind[rbeg + nzcnt];
		cbeg = f->uc_inf[urindx[rbeg + nzcnt]].cbeg;
		ucrind[cbeg + urcind[rbeg + nzcnt]] = ind;
		urindx[rbeg + nzcnt] = -1;
	}
	ur_inf[row].nzcnt = nzcnt;
	#ifdef DEBUG_FACTOR
	CLEANUP:
	#endif
	return rval;
}

static void delete_nonzero_col (
	factor_work * f,
	int col,
	int ind)
{
	uc_info *uc_inf = f->uc_inf;
	EGlpNum_t *uccoef = f->uccoef;
	int *ucindx = f->ucindx;
	int *ucrind = f->ucrind;
	int *urcind = f->urcind;
	int cbeg = uc_inf[col].cbeg;
	int nzcnt = uc_inf[col].nzcnt - 1;
	int rbeg;

	if (ind != nzcnt)
	{
		EGlpNumCopy (uccoef[cbeg + ind], uccoef[cbeg + nzcnt]);
		ucindx[cbeg + ind] = ucindx[cbeg + nzcnt];
		ucrind[cbeg + ind] = ucrind[cbeg + nzcnt];
		rbeg = f->ur_inf[ucindx[cbeg + nzcnt]].rbeg;
		urcind[rbeg + ucrind[cbeg + nzcnt]] = ind;
		ucindx[cbeg + nzcnt] = -1;
	}
	uc_inf[col].nzcnt = nzcnt;
}

static int delete_column (
	factor_work * f,
	int col)
{
	uc_info *uc_inf = f->uc_inf;
	int beg = uc_inf[col].cbeg;
	int nzcnt = uc_inf[col].nzcnt;
	int *ucindx = f->ucindx + beg;
	int *ucrind = f->ucrind + beg;
	int i, rval = 0;
	#ifdef DEBUG_FACTOR
	rval = check_matrix(f);
	TESTG(rval,CLEANUP,"Corrupted Matrix");
	#endif

	for (i = 0; i < nzcnt; i++)
	{
		rval = delete_nonzero_row (f, ucindx[i], ucrind[i]);
		CHECKRVALG(rval,CLEANUP);
		ucindx[i] = -1;
	}
	uc_inf[col].nzcnt = 0;

#ifdef TRACK_FACTOR
	f->nzcnt_cur -= nzcnt;
#endif

	#ifdef DEBUG_FACTOR
	rval = check_matrix(f);
	TESTG(rval,CLEANUP,"Corrupted Matrix");
	#endif /* DEBUG_FACTOR */
	CLEANUP:
	EG_RETURN(rval);
}

static int delete_row (
	factor_work * f,
	int row,
	svector * x)
{
	ur_info *ur_inf = f->ur_inf;
	int beg = ur_inf[row].rbeg;
	int nzcnt = ur_inf[row].nzcnt;
	int *urindx = f->urindx + beg;
	EGlpNum_t *urcoef = f->urcoef + beg;
	int *urcind = f->urcind + beg;
	int i,rval=0;

	#ifdef DEBUG_FACTOR
	rval = check_matrix(f);
	TESTG(rval,CLEANUP,"Corrupted Matrix");
	#endif /* DEBUG_FACTOR */

	for (i = 0; i < nzcnt; i++)
	{
		x->indx[i] = urindx[i];
		EGlpNumCopy (x->coef[i], urcoef[i]);
		delete_nonzero_col (f, urindx[i], urcind[i]);
		urindx[i] = -1;
	}
	x->nzcnt = nzcnt;
	ur_inf[row].nzcnt = 0;

#ifdef TRACK_FACTOR
	f->nzcnt_cur -= nzcnt;
#endif

	#ifdef DEBUG_FACTOR
	rval = check_matrix(f);
	TESTG(rval,CLEANUP,"Corrupted Matrix");
	CLEANUP:
	#endif /* DEBUG_FACTOR */
	return rval;
}

static int create_column (
	factor_work * f,
	svector * a,
	int col,
	int *p_last_rank)
{
	int *rrank = f->rrank;
	int nzcnt = a->nzcnt;
	int *aindx = a->indx;
	EGlpNum_t *acoef = a->coef;
	int i;
	int j;
	int rval = 0;
	int last_rank = -1;

#ifdef TRACK_FACTOR
	EGlpNum_t max;

	EGlpNumInitVar (max);
	EGlpNumCopy (max, f->maxelem_cur);
#endif /* TRACK_FACTOR */

	last_rank = 0;

	for (i = 0; i < nzcnt; i++)
	{
		rval = add_nonzero (f, aindx[i], col, acoef[i]);
		CHECKRVALG (rval, CLEANUP);
#ifdef TRACK_FACTOR
		EGlpNumSetToMaxAbs (max, acoef[i]);
#endif /* TRACK_FACTOR */
		j = rrank[aindx[i]];
		if (j > last_rank)
			last_rank = j;
	}
	*p_last_rank = last_rank;

#ifdef TRACK_FACTOR
	f->nzcnt_cur += nzcnt;
	EGlpNumCopy (f->maxelem_cur, max);
	EGlpNumClearVar (max);
#endif /* TRACK_FACTOR */

	#ifdef DEBUG_FACTOR
	rval = check_matrix(f);
	TESTG(rval,CLEANUP,"Corrupted Matrix");
	#endif /* DEBUG_FACTOR */

	CLEANUP:
	EG_RETURN (rval);
}

#ifdef UPDATE_STUDY
static int column_rank (
	factor_work * f,
	int col)
{
	int *cperm = f->cperm;
	int dim = f->dim;
	int i;

	for (i = 0; i < dim; i++)
	{
		if (cperm[i] == col)
		{
			return i;
		}
	}
	return 0;
}
#endif

static void shift_permutations (
	factor_work * f,
	int rank_p,
	int rank_r)
{
	int *cperm = f->cperm;
	int *crank = f->crank;
	int *rperm = f->rperm;
	int *rrank = f->rrank;
	int col_p = cperm[rank_p];
	int row_p = rperm[rank_p];
	int i;

	for (i = rank_p; i < rank_r; i++)
	{
		cperm[i] = cperm[i + 1];
		crank[cperm[i]] = i;
		rperm[i] = rperm[i + 1];
		rrank[rperm[i]] = i;
	}
	cperm[rank_r] = col_p;
	crank[col_p] = rank_r;
	rperm[rank_r] = row_p;
	rrank[row_p] = rank_r;
}

static int eliminate_row (
	factor_work * f,
	int rank_p,
	int rank_r)
{
	ur_info *ur_inf = f->ur_inf;
	int *rperm = f->rperm;
	int *cperm = f->cperm;
	int *urindx = f->urindx;
	EGlpNum_t *urcoef = f->urcoef;
	int *erindx = f->erindx;
	EGlpNum_t *ercoef = f->ercoef;
	EGlpNum_t *work_coef = f->work_coef;
	int er_freebeg = f->er_freebeg;
	int er_space = f->er_space;
	int beg;
	int nzcnt;
	int i;
	int j;
	int c;
	int r;
	EGlpNum_t pivot_mul;

#ifdef TRACK_FACTOR
	EGlpNum_t max;

	EGlpNumInitVar (max);
	EGlpNumCopy (max, f->maxelem_cur);
#endif /* TRACK_FACTOR */
	EGlpNumInitVar (pivot_mul);

	for (i = rank_p; i < rank_r; i++)
	{
		c = cperm[i];
		if (EGlpNumIsNeqZero (work_coef[c], f->fzero_tol))	/*
																												 * if (work_coef[c] > fzero_tol || work_coef[c] < -fzero_tol) */
		{
			r = rperm[i];
			beg = ur_inf[r].rbeg;
			nzcnt = ur_inf[r].nzcnt;
			EGlpNumCopyFrac (pivot_mul, work_coef[c], urcoef[beg]);
			EGlpNumZero (work_coef[c]);
			for (j = 1; j < nzcnt; j++)
			{
				EGlpNumSubInnProdTo (work_coef[urindx[beg + j]], pivot_mul, urcoef[beg + j]);	/* 0.85 */
			}
			if (er_freebeg >= er_space)
			{
				/* fprintf (stderr, "no space in eliminate_row\n"); */
#ifdef TRACK_FACTOR
				EGlpNumClearVar (max);
#endif
				EGlpNumClearVar (pivot_mul);
				return E_UPDATE_NOSPACE;
			}
			erindx[er_freebeg] = r;
			EGlpNumCopy (ercoef[er_freebeg], pivot_mul);
#ifdef TRACK_FACTOR
			EGlpNumSetToMaxAbs (max, pivot_mul);
#endif /* TRACK_FACTOR */
			er_freebeg++;
		}
		else
		{
			EGlpNumZero (work_coef[c]);
		}
	}
	f->er_freebeg = er_freebeg;
#ifdef TRACK_FACTOR
	EGlpNumCopy (f->maxelem_cur, max);
	EGlpNumClearVar (max);
#endif /* TRACK_FACTOR */
	EGlpNumClearVar (pivot_mul);
	return 0;
}

static int create_row (
	factor_work * f,
	EGlpNum_t * a,
	int row,
	int minrank)
{
	int *cperm = f->cperm;
	int dim = f->dim;
	int i;
	int j;
	int rval = 0;

#ifdef TRACK_FACTOR
	EGlpNum_t max;

	EGlpNumInitVar (max);
	EGlpNumCopy (max, f->maxelem_cur);
#endif /* TRACK_FACTOR */

	for (i = minrank; i < dim; i++)
	{
		if (EGlpNumIsNeqqZero (a[cperm[i]]))
		{
			j = cperm[i];
			if (EGlpNumIsNeqZero (a[j], f->fzero_tol))	/*
																									 * if (a[j] > fzero_tol || a[j] < -fzero_tol) */
			{
				rval = add_nonzero (f, row, j, a[j]);
				CHECKRVALG (rval, CLEANUP);
#ifdef TRACK_FACTOR
				EGlpNumSetToMaxAbs (max, a[j]);
#endif /* TRACK_FACTOR */
			}
			EGlpNumZero (a[j]);
		}
	}

#ifdef TRACK_FACTOR
	f->nzcnt_cur += f->ur_inf[row].nzcnt;
	EGlpNumCopy (f->maxelem_cur, max);
	EGlpNumClearVar (max);
#endif /* TRACK_FACTOR */

	#ifdef DEBUG_FACTOR
	rval = check_matrix(f);
	TESTG(rval,CLEANUP,"Corrupted Matrix");
	#endif /* DEBUG_FACTOR */
	CLEANUP:
	EG_RETURN (rval);
}

static void serow_delay (
	factor_work * f,
	int r,
	int rank_r)
{
	ur_info *ur_inf = f->ur_inf;
	int *crank = f->crank;
	int nzcnt;
	int *indx;
	int i;
	int last;

	do
	{
		r = f->rperm[crank[r]];
		nzcnt = ur_inf[r].nzcnt;
		indx = f->urindx + ur_inf[r].rbeg;
		last = -1;
		for (i = 1; i < nzcnt; i++)
		{
			r = indx[i];
			if (ur_inf[r].delay++ == 0 && crank[r] < rank_r)
			{
				if (last >= 0)
				{
					serow_delay (f, last, rank_r);
				}
				last = r;
			}
		}
		r = last;
	} while (r >= 0);
}

static int serow_process (
	factor_work * f,
	int r,
	svector * newr,
	int rank_r)
{
	ur_info *ur_inf = f->ur_inf;
	EGlpNum_t *work = f->work_coef;
	int nzcnt;
	int *indx;
	EGlpNum_t *coef;
	int i;
	EGlpNum_t v;
	int last;
	int rval;

	EGlpNumInitVar (v);

	do
	{
		EGlpNumCopy (v, work[r]);
		EGlpNumZero (work[r]);
		if (f->crank[r] >= rank_r)
		{
			if (EGlpNumIsNeqZero (v, f->fzero_tol))	/*
																							 * if (v > fzero_tol || v < -fzero_tol) */
			{
				/* stash this nonzero in the resulting row */
#ifdef TRACK_FACTOR
				EGlpNumSetToMaxAbs (f->maxelem_cur, v);
#endif /* TRACK_FACTOR */
				newr->indx[newr->nzcnt] = r;
				EGlpNumCopy (newr->coef[newr->nzcnt], v);
				newr->nzcnt++;
				EGlpNumClearVar (v);
				return 0;
			}
			else
			{
				EGlpNumClearVar (v);
				return 0;
			}
		}
		r = f->rperm[f->crank[r]];
		nzcnt = ur_inf[r].nzcnt;
		indx = f->urindx + ur_inf[r].rbeg;
		coef = f->urcoef + ur_inf[r].rbeg;
		EGlpNumDivTo (v, coef[0]);
		if (EGlpNumIsNeqZero (v, f->fzero_tol))	/*
																						 * if (v > fzero_tol || v < -fzero_tol) */
		{
			/* stash v in eta */
			if (f->er_freebeg >= f->er_space)
			{
				/* fprintf (stderr, "no space in eliminate_row\n"); */
				EGlpNumClearVar (v);
				return E_UPDATE_NOSPACE;
			}
			f->erindx[f->er_freebeg] = r;
			EGlpNumCopy (f->ercoef[f->er_freebeg], v);
#ifdef TRACK_FACTOR
			EGlpNumSetToMaxAbs (f->maxelem_cur, v);
#endif /* TRACK_FACTOR */
			f->er_freebeg++;
		}
		last = -1;
		for (i = 1; i < nzcnt; i++)
		{
			r = indx[i];
			EGlpNumSubInnProdTo (work[r], v, coef[i]);
			if (--ur_inf[r].delay == 0)
			{
				if (last >= 0)
				{
					rval = serow_process (f, last, newr, rank_r);
					if (rval)
					{
						EGlpNumClearVar (v);
						return rval;
					}
				}
				last = r;
			}
		}
		r = last;
	} while (r >= 0);
	EGlpNumClearVar (v);
	return 0;
}

static int sparse_eliminate_row (
	factor_work * f,
	svector * x,
	int row_p,
	int rank_r)
{
	EGlpNum_t *work = f->work_coef;
	int xnzcnt = x->nzcnt;
	int *xindx = x->indx;
	EGlpNum_t *xcoef = x->coef;
	ur_info *ur_inf = f->ur_inf;
	int *crank = f->crank;
	int i;
	int j;
	int rval = 0;
	svector newr;

	newr.indx = 0;
	newr.coef = 0;

	for (i = 0; i < xnzcnt; i++)
	{
		j = xindx[i];
		if (ur_inf[j].delay++ == 0 && crank[j] < rank_r)
		{
			serow_delay (f, j, rank_r);
		}
		EGlpNumCopy (work[j], xcoef[i]);
	}

	newr.nzcnt = 0;
	ILL_SAFE_MALLOC (newr.indx, f->dim, int);

	newr.coef = EGlpNumAllocArray (f->dim);

	for (i = 0; i < xnzcnt; i++)
	{
		j = xindx[i];
		if (--ur_inf[j].delay == 0)
		{
			rval = serow_process (f, j, &newr, rank_r);
			CHECKRVALG (rval, CLEANUP);
		}
	}

	for (i = 0; i < newr.nzcnt; i++)
	{
		rval = add_nonzero (f, row_p, newr.indx[i], newr.coef[i]);
		CHECKRVALG (rval, CLEANUP);
	}

#ifdef TRACK_FACTOR
	f->nzcnt_cur += newr.nzcnt;
#endif /* TRACK_FACTOR */

CLEANUP:
	EGlpNumFreeArray (newr.coef);
	ILL_IFFREE (newr.indx, int);

	/* Bico 031210 - chg from ILL_RETURN */
	EG_RETURN (rval);
}

static int move_pivot_row (
	factor_work * f,
	int r,
	int c)
{
	ur_info *ur_inf = f->ur_inf + r;
	uc_info *uc_inf = f->uc_inf;
	int beg = ur_inf->rbeg;
	int nzcnt = ur_inf->nzcnt;
	int *urindx = f->urindx;
	int *urcind = f->urcind;
	int *ucrind = f->ucrind;
	EGlpNum_t *urcoef = f->urcoef;
	EGlpNum_t dt;
	int it;
	int i;

	if (urindx[beg] == c)
		return 0;
	EGlpNumInitVar (dt);

	for (i = 1; i < nzcnt; i++)
	{
		if (urindx[beg + i] == c)
		{
			EGLPNUM_SWAP (urcoef[beg], urcoef[beg + i], dt);
			ILL_SWAP (urcind[beg], urcind[beg + i], it);
			urindx[beg + i] = urindx[beg];
			urindx[beg] = c;
			ucrind[uc_inf[c].cbeg + urcind[beg]] = 0;
			ucrind[uc_inf[urindx[beg + i]].cbeg + urcind[beg + i]] = i;
			EGlpNumClearVar (dt);
			return 0;
		}
	}
	MESSAGE (__QS_SB_VERB, "pivot row nonzero not found");
	EGlpNumClearVar (dt);
	return E_UPDATE_SINGULAR_ROW;
}

static int move_pivot_col (
	factor_work * f,
	int c,
	int r)
{
	uc_info *uc_inf = f->uc_inf + c;
	ur_info *ur_inf = f->ur_inf;
	int beg = uc_inf->cbeg;
	int nzcnt = uc_inf->nzcnt;
	int *ucindx = f->ucindx;
	int *ucrind = f->ucrind;
	int *urcind = f->urcind;
	EGlpNum_t *uccoef = f->uccoef;
	EGlpNum_t dt;
	int i, it;

	if (ucindx[beg] == r)
		return 0;
	EGlpNumInitVar (dt);

	for (i = 1; i < nzcnt; i++)
	{
		if (ucindx[beg + i] == r)
		{
			EGLPNUM_SWAP (uccoef[beg], uccoef[beg + i], dt);
			ILL_SWAP (ucrind[beg], ucrind[beg + i], it);
			ucindx[beg + i] = ucindx[beg];
			ucindx[beg] = r;
			urcind[ur_inf[r].rbeg + ucrind[beg]] = 0;
			urcind[ur_inf[ucindx[beg + i]].rbeg + ucrind[beg + i]] = i;
			EGlpNumClearVar (dt);
			return 0;
		}
	}
	MESSAGE(__QS_SB_VERB, "pivot col nonzero not found");
	EGlpNumClearVar (dt);
	return E_UPDATE_SINGULAR_COL;
}

static int move_pivot (
	factor_work * f,
	int rank_r)
{
	int r = f->rperm[rank_r];
	int c = f->cperm[rank_r];
	int rval = 0;

	rval = move_pivot_row (f, r, c);
	CHECKRVALG (rval, CLEANUP);

	#ifdef DEBUG_FACTOR
	rval = check_matrix(f);
	TESTG(rval,CLEANUP,"Corrupted Matrix");
	#endif /* DEBUG_FACTOR */

	rval = move_pivot_col (f, c, r);
	if(rval != E_UPDATE_SINGULAR_COL) CHECKRVALG (rval, CLEANUP);
	else goto CLEANUP;

	#ifdef DEBUG_FACTOR
	rval = check_matrix(f);
	TESTG(rval,CLEANUP,"Corrupted Matrix");
	#endif /* DEBUG_FACTOR */

	CLEANUP:
	if(rval != E_UPDATE_SINGULAR_COL) EG_RETURN (rval);							/* Bico 031209 - chg from RETURN */
	return rval;
}

int ILLfactor_update (
	factor_work * f,
	svector * a,
	int col_p,
	int *p_refact)
{
	int row_p;
	int rank_r = 0;
	int rank_p = 0;
	int rval = 0;
	int nzcnt;
	int *aindx;
	EGlpNum_t *acoef;
	EGlpNum_t *work_coef = f->work_coef;

#ifdef TRACK_FACTOR
#ifdef NOTICE_BLOWUP
	EGlpNum_t tmpsize;
#endif
#endif
	int i;

#ifdef RECORD
	{
		EGioPrintf (fsave, "u %d %d", col_p, a->nzcnt);
		for (i = 0; i < a->nzcnt; i++)
		{
			EGioPrintf (fsave, " %d %.16e", a->indx[i], EGlpNumToLf (a->coef[i]));
		}
		EGioPrintf (fsave, "\n");
		EGioFlush (fsave);
	}
#endif /* RECORD */

#ifdef DEBUG_FACTOR
	{
		printf ("ILLfactor_update col %d:", col_p);
		for (i = 0; i < a->nzcnt; i++)
		{
			printf (" %.3f*%d", EGlpNumToLf (a->coef[i]), a->indx[i]);
		}
		printf ("\n");
		fflush (stdout);
	}
#endif /* DEBUG_FACTOR */

	#ifdef DEBUG_FACTOR
	rval = check_matrix(f);
	TESTG(rval,CLEANUP,"Corrupted Matrix");
	#endif /* DEBUG_FACTOR */

	if (f->etacnt >= f->etamax)
	{
		*p_refact = 1;
		return 0;
	}

#ifdef UPDATE_STUDY
	nupdate++;
#endif

	row_p = f->ucindx[f->uc_inf[col_p].cbeg];

	rval = delete_column (f, col_p);
	CHECKRVALG (rval, CLEANUP);

	rval = create_column (f, a, col_p, &rank_r);
	/* if (rval) fprintf (stderr, "create_column failed\n"); */
	CHECKRVALG (rval, CLEANUP);

	rank_p = f->crank[col_p];
#ifdef UPDATE_STUDY
	if (rank_p != f->rrank[row_p] || rank_p != column_rank (f, col_p))
	{
		printf ("rank_p %d rrank[row_p] %d column_rank(f,col_p) %d\n",
						rank_p, f->rrank[row_p], column_rank (f, col_p));
	}
	if (rank_r > rank_p)
	{
		permshifttot += rank_r - rank_p;
	}
	for (i = 0; i < a->nzcnt; i++)
	{
		if (f->rrank[a->indx[i]] > rank_p)
			colspiketot++;
	}
	for (i = 0; i < f->ur_inf[row_p].nzcnt; i++)
	{
		if (f->crank[f->urindx[f->ur_inf[row_p].rbeg + i]] <= rank_r &&
				f->crank[f->urindx[f->ur_inf[row_p].rbeg + i]] != rank_p)
		{
			rowspiketot++;
		}
	}
#endif

	shift_permutations (f, rank_p, rank_r);

	rval = delete_row (f, row_p, &f->xtmp);
	CHECKRVALG(rval,CLEANUP);

	f->er_inf[f->etacnt].rbeg = f->er_freebeg;
	f->er_inf[f->etacnt].r = row_p;

	if (f->xtmp.nzcnt >= SPARSE_FACTOR * f->dim)
	{
		nzcnt = f->xtmp.nzcnt;
		aindx = f->xtmp.indx;
		acoef = f->xtmp.coef;

		for (i = 0; i < nzcnt; i++)
		{
			EGlpNumCopy (work_coef[aindx[i]], acoef[i]);
		}

		rval = eliminate_row (f, rank_p, rank_r);
		/* if (rval) fprintf (stderr, "eliminate_row failed\n"); */
		CHECKRVALG (rval, CLEANUP);

		rval = create_row (f, f->work_coef, row_p, rank_r);
		/* if (rval) fprintf (stderr, "create_row failed\n"); */
		CHECKRVALG (rval, CLEANUP);
	}
	else
	{
		rval = sparse_eliminate_row (f, &f->xtmp, row_p, rank_r);
		/* if (rval) fprintf (stderr, "sparse_eliminate_row failed\n"); */
		CHECKRVALG (rval, CLEANUP);
	}

	if (f->er_freebeg - f->er_inf[f->etacnt].rbeg > 0)
	{
		f->er_inf[f->etacnt].nzcnt = f->er_freebeg - f->er_inf[f->etacnt].rbeg;
#ifdef TRACK_FACTOR
		f->nzcnt_cur += f->er_inf[f->etacnt].nzcnt;
#endif /* TRACK_FACTOR */
#ifdef UPDATE_STUDY
		leftetatot += f->er_inf[f->etacnt].nzcnt;
#endif

#ifdef SORT_RESULTS
		sort_vector2 (f->er_inf[f->etacnt].nzcnt,
									f->erindx + f->er_inf[f->etacnt].rbeg,
									f->ercoef + f->er_inf[f->etacnt].rbeg);
#endif

		f->etacnt++;
	}

	rval = move_pivot (f, rank_r);
	/* if (rval) fprintf (stderr, "move_pivot failed\n"); */
	if(rval != E_UPDATE_SINGULAR_COL) CHECKRVALG (rval, CLEANUP);
	else goto CLEANUP;

#ifdef UPDATE_DEBUG
	printf ("Updated factorization:\n");
#if (UPDATE_DEBUG+0>1)
	dump_matrix (f, 0);
#endif
	fflush (stdout);
#endif /* UPDATE_DEBUG */

#ifdef TRACK_FACTOR
#ifdef NOTICE_BLOWUP
	EGlpNumInitVar (tmpsize);
	EGlpNumSet (tmpsize, f->updmaxmult);
	EGlpNumMultTo (tmpsize, f->maxelem_orig);
	if (EGlpNumIsLess (tmpsize, f->maxelem_cur))
	{
/* Bico - comment out for dist 
        fprintf (stderr, "factor_update blowup max cur %e max orig %e\n",
                 f->maxelem_cur, f->maxelem_orig);
*/
		EGlpNumClearVar (tmpsize);
		return E_FACTOR_BLOWUP;
	}
	EGlpNumClearVar (tmpsize);
#endif /* NOTICE_BLOWUP */
#endif
#ifdef UPDATE_STATS
	dump_factor_stats (f);
#endif
CLEANUP:
	if(rval != E_UPDATE_SINGULAR_COL) EG_RETURN (rval);							/* Bico 031209 - chg from RETURN */
	return rval;
}
