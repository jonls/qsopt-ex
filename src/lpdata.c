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

/* RCS_INFO = "$RCSfile: lpdata.c,v $ $Revision: 1.2 $ $Date: 2003/11/05 16:49:52 $"; */
/****************************************************************************/
/*                                                                          */
/*               Routines for Manipulating and Writing LPdata               */
/*                                                                          */
/*  EXPORTED FUNCTIONS                                                      */
/*                                                                          */
/*    int ILLlpdata_buildrows (ILLlpdata *lp, int **rowbeg, int **rowcnt,     */
/*            int **rowind, double **rowval, int include_logicals)          */
/*      - include_logicals:  if nonzero, then logical variables will be     */
/*          included in the row data                                        */
/*                                                                          */
/*                                                                          */
/*    All _init routines initialize fields of allocated structure to        */
/*    appropiate default values                                             */
/*    The _free routines free structures contained in pareameter structure  */
/*    but not the parameter itself.                                         */
/*    The _alloc routines check whether given parameter is NULL; they either*/
/*    print an error message or fill structure with default values or the   */
/*    given paremeter values.                                               */
/*                                                                          */
/*    void ILLlpdata_init (ILLlpdata *lp)                                   */
/*    void ILLlpdata_free (ILLlpdata *lp)                                   */
/*                                                                          */
/*    void ILLlp_basis_init (ILLlp_basis *B)                                */
/*    void ILLlp_basis_free (ILLlp_basis *B)                                */
/*    int ILLlp_basis_alloc (ILLlp_basis *B, int nstruct, int nrows)        */
/*                                                                          */
/*    void ILLlp_cache_init (ILLlp_cache *C)                                */
/*    void ILLlp_cache_free (ILLlp_cache *C)                                */
/*    int ILLlp_cache_alloc (ILLlp_cache *C, int nstruct, int nrows)        */
/*                                                                          */
/*    void ILLlp_sinfo_init (ILLlp_sinfo *sinfo)                            */
/*    void ILLlp_sinfo_free (ILLlp_sinfo *sinfo)                            */
/*                                                                          */
/*    int ILLlp_rows_init(ILLlp_rows *lprows, ILLlpdata *lp,                */
/*                                           int include_logicals)          */
/*                                                                          */
/****************************************************************************/

#include <stdlib.h>
#include <stdarg.h>

#include "qs_config.h"

#include "eg_lpnum.h"
#include "eg_io.h"

#include "iqsutil.h"
#include "lpdata.h"
#include "qstruct.h"
#include "qsopt.h"
#include "lp.h"
#include "mps.h"
#include "rawlp.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

//static int TRACE = 0;

EGlpNum_t PARAM_IBASIS_RPIVOT;
EGlpNum_t PARAM_IBASIS_RTRIANG;
EGlpNum_t PARAM_MIN_DNORM;
EGlpNum_t PFEAS_TOLER;
EGlpNum_t BD_TOLER;
EGlpNum_t DFEAS_TOLER;
EGlpNum_t PIVOT_TOLER;
EGlpNum_t SZERO_TOLER;
EGlpNum_t PIVZ_TOLER;
EGlpNum_t OBJBND_TOLER;
EGlpNum_t DBNDPIV_TOLER;
EGlpNum_t DBNDPIV_RATIO;
EGlpNum_t ALTPIV_TOLER;
//EGlpNum_t DJZERO_TOLER;
EGlpNum_t PROGRESS_ZERO;				/*   1e-7 */
EGlpNum_t PROGRESS_THRESH;			/*   1e-5 */
EGlpNum_t CB_EPS;
EGlpNum_t CB_INF_RATIO;
EGlpNum_t CB_PRI_RLIMIT;
EGlpNum_t ILL_MAXDOUBLE;
EGlpNum_t ILL_MINDOUBLE;

/* ========================================================================= */
int __QSEX_SETUP = 0;
/* ========================================================================= */
void ILLstart ( void)
{
	if (__QSEX_SETUP)
		return;
	EGlpNumInitVar (PARAM_IBASIS_RPIVOT);
	EGlpNumInitVar (PARAM_IBASIS_RTRIANG);
	EGlpNumInitVar (PARAM_MIN_DNORM);
	EGlpNumInitVar (PFEAS_TOLER);
	EGlpNumInitVar (BD_TOLER);
	EGlpNumInitVar (DFEAS_TOLER);
	EGlpNumInitVar (PIVOT_TOLER);
	EGlpNumInitVar (SZERO_TOLER);
	EGlpNumInitVar (PIVZ_TOLER);
	EGlpNumInitVar (OBJBND_TOLER);
	EGlpNumInitVar (DBNDPIV_TOLER);
	EGlpNumInitVar (DBNDPIV_RATIO);
	EGlpNumInitVar (ALTPIV_TOLER);
	//EGlpNumInitVar (DJZERO_TOLER);
	EGlpNumInitVar (PROGRESS_ZERO);	/*            1e-7 */
	EGlpNumInitVar (PROGRESS_THRESH);	/*          1e-5 */
	EGlpNumInitVar (CB_PRI_RLIMIT);
	EGlpNumInitVar (CB_INF_RATIO);
	EGlpNumInitVar (CB_EPS);
	EGlpNumInitVar (ILL_MAXDOUBLE);
	EGlpNumInitVar (ILL_MINDOUBLE);
	/* parameters that do depend on the tolerance to zero */
	EGlpNumSet (PARAM_MIN_DNORM, 4.5036e-9);
	EGlpNumMultTo (PARAM_MIN_DNORM, epsLpNum);
	EGlpNumSet (PFEAS_TOLER, 4.5036e9);
	EGlpNumMultTo (PFEAS_TOLER, epsLpNum);
	EGlpNumSet (BD_TOLER, 4.5036e8);
	EGlpNumMultTo (BD_TOLER, epsLpNum);
	EGlpNumSet (DFEAS_TOLER, 4.5036e9);
	EGlpNumMultTo (DFEAS_TOLER, epsLpNum);
	EGlpNumSet (PIVOT_TOLER, 4.5036e5);
	EGlpNumMultTo (PIVOT_TOLER, epsLpNum);
	EGlpNumSet (SZERO_TOLER, 4.5036);
	EGlpNumMultTo (SZERO_TOLER, epsLpNum);
	EGlpNumSet (PIVZ_TOLER, 4.5036e3);
	EGlpNumMultTo (PIVZ_TOLER, epsLpNum);
	EGlpNumSet (OBJBND_TOLER, 4.5036e13);
	EGlpNumMultTo (OBJBND_TOLER, epsLpNum);
	EGlpNumSet (ALTPIV_TOLER, 4.5036e7);
	EGlpNumMultTo (ALTPIV_TOLER, epsLpNum);
	EGlpNumSet (PROGRESS_ZERO, 4.5036e8);
	EGlpNumMultTo (PROGRESS_ZERO, epsLpNum);
	EGlpNumSet (PROGRESS_THRESH, 4.5036e10);
	EGlpNumMultTo (PROGRESS_THRESH, epsLpNum);
#if VERBOSE_LEVEL <= DEBUG
	MESSAGE (VERBOSE_LEVEL, "Setting PARAM_MIN_DNORM to %lg", EGlpNumToLf (PARAM_MIN_DNORM));
	MESSAGE (VERBOSE_LEVEL, "Setting PFEAS_TOLER to %lg", EGlpNumToLf (PFEAS_TOLER));
	MESSAGE (VERBOSE_LEVEL, "Setting BD_TOLER to %lg", EGlpNumToLf (BD_TOLER));
	MESSAGE (VERBOSE_LEVEL, "Setting DFEAS_TOLER to %lg", EGlpNumToLf (DFEAS_TOLER));
	MESSAGE (VERBOSE_LEVEL, "Setting PIVOT_TOLER to %lg", EGlpNumToLf (PIVOT_TOLER));
	MESSAGE (VERBOSE_LEVEL, "Setting SZERO_TOLER to %lg", EGlpNumToLf (SZERO_TOLER));
	MESSAGE (VERBOSE_LEVEL, "Setting PIVZ_TOLER to %lg", EGlpNumToLf (PIVZ_TOLER));
	MESSAGE (VERBOSE_LEVEL, "Setting OBJBND_TOLER to %lg", EGlpNumToLf (OBJBND_TOLER));
	MESSAGE (VERBOSE_LEVEL, "Setting ALTPIV_TOLER to %lg", EGlpNumToLf (ALTPIV_TOLER));
	MESSAGE (VERBOSE_LEVEL, "Setting PROGRESS_ZERO to %lg", EGlpNumToLf (PROGRESS_ZERO));
	MESSAGE (VERBOSE_LEVEL, "Setting PROGRESS_THRESH to %lg", EGlpNumToLf (PROGRESS_THRESH));
#endif
	/* parameters that do not depend on the tolerance to zero */
	EGlpNumSet (ILL_MAXDOUBLE, 1e150);
	EGlpNumSet (ILL_MINDOUBLE, -1e150);
	EGlpNumSet (PARAM_IBASIS_RPIVOT, 0.98);
	EGlpNumSet (PARAM_IBASIS_RTRIANG, 0.01);
	EGlpNumSet (DBNDPIV_TOLER, 1e-3);
	EGlpNumSet (DBNDPIV_RATIO, 1e-2);
	//EGlpNumSet (DJZERO_TOLER, 1e-8);
	EGlpNumSet (CB_EPS, 0.001);
	EGlpNumSet (CB_INF_RATIO, 10.0);
	EGlpNumSet (CB_PRI_RLIMIT, 0.25);
	__QSEX_SETUP = 1;
}

/* ========================================================================= */
void ILLchange_precision (
	void)
{
	EGlpNumClearVar (PFEAS_TOLER);
	EGlpNumClearVar (BD_TOLER);
	EGlpNumClearVar (DFEAS_TOLER);
	EGlpNumClearVar (PIVOT_TOLER);
	EGlpNumClearVar (SZERO_TOLER);
	EGlpNumClearVar (PIVZ_TOLER);
	EGlpNumClearVar (OBJBND_TOLER);
	EGlpNumClearVar (ALTPIV_TOLER);
	EGlpNumClearVar (PARAM_MIN_DNORM);
	EGlpNumClearVar (PROGRESS_ZERO);
	EGlpNumClearVar (PROGRESS_THRESH);
	EGlpNumInitVar (PROGRESS_ZERO);
	EGlpNumInitVar (PROGRESS_THRESH);
	EGlpNumInitVar (PFEAS_TOLER);
	EGlpNumInitVar (BD_TOLER);
	EGlpNumInitVar (DFEAS_TOLER);
	EGlpNumInitVar (PIVOT_TOLER);
	EGlpNumInitVar (SZERO_TOLER);
	EGlpNumInitVar (PIVZ_TOLER);
	EGlpNumInitVar (OBJBND_TOLER);
	EGlpNumInitVar (ALTPIV_TOLER);
	EGlpNumInitVar (PARAM_MIN_DNORM);
	/* parameters that do depend on the tolerance to zero */
	EGlpNumSet (PARAM_MIN_DNORM, 4.5036e-9);
	EGlpNumMultTo (PARAM_MIN_DNORM, epsLpNum);
	EGlpNumSet (PFEAS_TOLER, 4.5036e9);
	EGlpNumMultTo (PFEAS_TOLER, epsLpNum);
	EGlpNumSet (BD_TOLER, 4.5036e8);
	EGlpNumMultTo (BD_TOLER, epsLpNum);
	EGlpNumSet (DFEAS_TOLER, 4.5036e9);
	EGlpNumMultTo (DFEAS_TOLER, epsLpNum);
	EGlpNumSet (PIVOT_TOLER, 4.5036e5);
	EGlpNumMultTo (PIVOT_TOLER, epsLpNum);
	EGlpNumSet (SZERO_TOLER, 4.5036);
	EGlpNumMultTo (SZERO_TOLER, epsLpNum);
	EGlpNumSet (PIVZ_TOLER, 4.5036e3);
	EGlpNumMultTo (PIVZ_TOLER, epsLpNum);
	EGlpNumSet (OBJBND_TOLER, 4.5036e13);
	EGlpNumMultTo (OBJBND_TOLER, epsLpNum);
	EGlpNumSet (ALTPIV_TOLER, 4.5036e7);
	EGlpNumMultTo (ALTPIV_TOLER, epsLpNum);
	EGlpNumSet (PROGRESS_ZERO, 4.5036e8);
	EGlpNumMultTo (PROGRESS_ZERO, epsLpNum);
	EGlpNumSet (PROGRESS_THRESH, 4.5036e10);
	EGlpNumMultTo (PROGRESS_THRESH, epsLpNum);
}

/* ========================================================================= */
void ILLend ( void)
{
	if (!__QSEX_SETUP)
		return;
	EGlpNumClearVar (PARAM_IBASIS_RPIVOT);
	EGlpNumClearVar (PARAM_IBASIS_RTRIANG);
	EGlpNumClearVar (PARAM_MIN_DNORM);
	EGlpNumClearVar (PFEAS_TOLER);
	EGlpNumClearVar (BD_TOLER);
	EGlpNumClearVar (DFEAS_TOLER);
	EGlpNumClearVar (PIVOT_TOLER);
	EGlpNumClearVar (SZERO_TOLER);
	EGlpNumClearVar (PIVZ_TOLER);
	EGlpNumClearVar (OBJBND_TOLER);
	EGlpNumClearVar (DBNDPIV_TOLER);
	EGlpNumClearVar (DBNDPIV_RATIO);
	EGlpNumClearVar (ALTPIV_TOLER);
	//EGlpNumClearVar (DJZERO_TOLER);
	EGlpNumClearVar (PROGRESS_ZERO);	/*            1e-7 */
	EGlpNumClearVar (PROGRESS_THRESH);	/*          1e-5 */
	EGlpNumClearVar (CB_EPS);
	EGlpNumClearVar (CB_INF_RATIO);
	EGlpNumClearVar (CB_PRI_RLIMIT);
	EGlpNumClearVar (ILL_MAXDOUBLE);
	EGlpNumClearVar (ILL_MINDOUBLE);
	__QSEX_SETUP = 0;
}

QSdata *ILLread (
	qsline_reader * file,
	const char *fname,
	int isMps)
{
	int rval = 0;
	QSdata *p = 0;
	ILLlpdata *lp;
	rawlpdata rawlp;

	ILL_FAILfalse (file != NULL, NULL);
	ILL_FAILfalse (fname != NULL, NULL);

	p = QScreate_prob (fname, QS_MIN);
	ILL_CHECKnull (p, NULL);
	ILL_IFFREE (p->qslp->probname, char);

	lp = p->qslp;

	ILLinit_rawlpdata (&rawlp, file->error_collector);
	ILLlpdata_init (lp);

	if (isMps != 0)
	{
		rval = ILLread_mps (file, fname, &rawlp);
	}
	else
	{
		rval = ILLread_lp (file, fname, &rawlp);
	}
	CHECKRVALG (rval, CLEANUP);

	rval = ILLrawlpdata_to_lpdata (&rawlp, lp);
	CHECKRVALG (rval, CLEANUP);

CLEANUP:
	ILLfree_rawlpdata (&rawlp);
	if (rval != 0)
	{
		QSfree_prob (p);
		p = 0;
	}
	return p;
}

void ILLlpdata_init (
	ILLlpdata * lp)
{
	if (lp)
	{
		lp->nrows = 0;
		lp->ncols = 0;
		lp->nstruct = 0;
		lp->nzcount = 0;
		lp->rowsize = 0;
		lp->colsize = 0;
		lp->structsize = 0;
		lp->objsense = ILL_MIN;
		lp->sense = 0;
		lp->obj = 0;
		lp->rhs = 0;
		lp->rangeval = 0;
		lp->lower = 0;
		lp->upper = 0;

		ILLmatrix_init (&lp->A);
		ILLmatrix_init (&lp->sos);
		lp->rA = 0;
		lp->is_sos_mem = NULL;
		lp->refrowname = NULL;
		lp->refind = -1;

		lp->colnames = 0;
		ILLsymboltab_init (&lp->coltab);
		lp->rownames = 0;
		ILLsymboltab_init (&lp->rowtab);
		lp->objname = 0;

		lp->probname = 0;
		lp->intmarker = 0;
		lp->structmap = 0;
		lp->rowmap = 0;
		lp->basis = 0;
		/*lp->presolve   = 0; */
		lp->sinfo = 0;

		ILLstring_reporter_init (&lp->reporter, ILL_fprintf, stdout);
	}
}

void ILLlpdata_free (
	ILLlpdata * lp)
{
	int i;

	if (lp)
	{
		ILL_IFFREE (lp->sense, char);

		EGlpNumFreeArray (lp->obj);
		EGlpNumFreeArray (lp->rhs);
		EGlpNumFreeArray (lp->rangeval);
		EGlpNumFreeArray (lp->lower);
		EGlpNumFreeArray (lp->upper);
		ILLmatrix_free (&lp->A);
		if (lp->rA)
		{
			ILLlp_rows_clear (lp->rA);
			ILL_IFFREE (lp->rA, ILLlp_rows);
		}
		ILL_IFFREE (lp->is_sos_mem, int);
		ILL_IFFREE (lp->refrowname, char);

		ILLmatrix_free (&lp->sos);
		if (lp->colnames)
		{
			for (i = 0; i < lp->nstruct; i++)
			{
				ILL_IFFREE (lp->colnames[i], char);
			}
			ILL_IFFREE (lp->colnames, char *);
		}
		ILLsymboltab_free (&lp->coltab);
		if (lp->rownames)
		{
			for (i = 0; i < lp->nrows; i++)
			{
				ILL_IFFREE (lp->rownames[i], char);
			}
			ILL_IFFREE (lp->rownames, char *);
		}
		ILLsymboltab_free (&lp->rowtab);
		ILL_IFFREE (lp->objname, char);
		ILL_IFFREE (lp->probname, char);
		ILL_IFFREE (lp->intmarker, char);
		ILL_IFFREE (lp->structmap, int);
		ILL_IFFREE (lp->rowmap, int);

		if (lp->sinfo)
		{
			ILLlp_sinfo_free (lp->sinfo);
			ILL_IFFREE (lp->sinfo, ILLlp_sinfo);
		}
		ILLlpdata_init (lp);
	}
}

void ILLlp_basis_init (
	ILLlp_basis * B)
{
	if (B)
	{
		B->cstat = 0;
		B->rstat = 0;
		B->rownorms = 0;
		B->colnorms = 0;
		B->nstruct = 0;
		B->nrows = 0;
	}
}

void ILLlp_basis_free (
	ILLlp_basis * B)
{
	if (B)
	{
		ILL_IFFREE (B->cstat, char);
		ILL_IFFREE (B->rstat, char);

		EGlpNumFreeArray (B->rownorms);
		EGlpNumFreeArray (B->colnorms);
		B->nstruct = 0;
		B->nrows = 0;
	}
}

int ILLlp_basis_alloc (
	ILLlp_basis * B,
	int nstruct,
	int nrows)
{
	int rval = 0;

	ILL_FAILtrue (B == NULL, "ILLlp_basis_alloc called without a basis");

	B->nstruct = nstruct;
	B->nrows = nrows;

	if (nstruct > 0)
	{
		ILL_SAFE_MALLOC (B->cstat, nstruct, char);
	}

	if (nrows > 0)
	{
		ILL_SAFE_MALLOC (B->rstat, nrows, char);
	}

CLEANUP:

	if (rval)
	{
		ILLlp_basis_free (B);
	}

	EG_RETURN (rval);
}

void ILLlp_cache_init (
	ILLlp_cache * C)
{
	if (C)
	{
		C->x = 0;
		C->rc = 0;
		C->pi = 0;
		C->slack = 0;
		C->nstruct = 0;
		C->nrows = 0;
		C->status = 0;
		EGlpNumZero (C->val);
	}
}

void ILLlp_cache_free (
	ILLlp_cache * C)
{
	if (C)
	{
		EGlpNumFreeArray (C->x);
		EGlpNumFreeArray (C->rc);
		EGlpNumFreeArray (C->pi);
		EGlpNumFreeArray (C->slack);
		C->nstruct = 0;
		C->nrows = 0;
		C->status = 0;
	}
}

int ILLlp_cache_alloc (
	ILLlp_cache * C,
	int nstruct,
	int nrows)
{
	int rval = 0;

	ILL_FAILtrue (C == NULL, "ILLlp_cache_alloc called without a cache");

	C->nstruct = nstruct;
	C->nrows = nrows;

	if (nstruct > 0)
	{
		C->x = EGlpNumAllocArray (nstruct);
		C->rc = EGlpNumAllocArray (nstruct);
	}

	if (nrows > 0)
	{
		C->pi = EGlpNumAllocArray (nrows);
		C->slack = EGlpNumAllocArray (nrows);
	}

CLEANUP:

	if (rval)
	{
		ILLlp_cache_free (C);
	}

	EG_RETURN (rval);
}


int ILLlp_rows_init (
	ILLlp_rows * lprows,
	ILLlpdata * lp,
	int include_logicals)
{
	int rval = 0;
	int i, k, st;
	int *beg, *cnt, *ind;
	EGlpNum_t *val;
	ILLmatrix *A;
	char *hit = 0;
	int *inv_structmap = 0;

	/* If logicals are not included, then the columns are ordered as in */
	/* lp->structmap.  Otherwise, the columns are ordered as in the     */
	/* matrix structure.                                                */

	if (lprows != NULL)
	{
		lprows->rowbeg = 0;
		lprows->rowcnt = 0;
		lprows->rowind = 0;
		lprows->rowval = 0;
	}

	ILL_FAILfalse ((lp != NULL) && (lprows != NULL),
								 "called with a NULL pointer");

	A = &lp->A;

	if (lp->nrows > 0)
	{
		if (include_logicals == 0)
		{
			ILL_FAILtrue (lp->rowmap == NULL, "Programming error.");
			ILL_SAFE_MALLOC (hit, lp->ncols, char);

			for (i = 0; i < lp->ncols; i++)
			{
				hit[i] = 0;
			}
			for (i = 0; i < lp->nrows; i++)
			{
				hit[lp->rowmap[i]] = 1;
			}

			ILL_SAFE_MALLOC (inv_structmap, lp->ncols, int);

			for (i = 0; i < lp->nstruct; i++)
			{
				inv_structmap[lp->structmap[i]] = i;
			}
		}

		ILL_SAFE_MALLOC (lprows->rowbeg, lp->nrows, int);
		ILL_SAFE_MALLOC (lprows->rowcnt, lp->nrows, int);

		if (((include_logicals != 0) && lp->nzcount > 0) ||
				((include_logicals == 0) && lp->nzcount > lp->nrows))
		{
			if (include_logicals != 0)
			{
				ILL_SAFE_MALLOC (lprows->rowind, lp->nzcount, int);

				lprows->rowval = EGlpNumAllocArray (lp->nzcount);
			}
			else
			{
				ILL_SAFE_MALLOC (lprows->rowind, lp->nzcount - lp->nrows, int);

				lprows->rowval = EGlpNumAllocArray (lp->nzcount - lp->nrows);
			}
		}

		beg = lprows->rowbeg;
		cnt = lprows->rowcnt;
		ind = lprows->rowind;
		val = lprows->rowval;

		for (i = 0; i < lp->nrows; i++)
		{
			cnt[i] = 0;
		}

		for (i = 0; i < lp->ncols; i++)
		{
			if ((include_logicals != 0) || hit[i] == 0)
			{
				k = A->matbeg[i];
				st = k + A->matcnt[i];
				for (; k < st; k++)
				{
					cnt[A->matind[k]]++;
				}
			}
		}

		for (i = 0, k = 0; i < lp->nrows; i++)
		{
			beg[i] = k;
			k += cnt[i];
		}

		for (i = 0; i < lp->ncols; i++)
		{
			if ((include_logicals != 0) || hit[i] == 0)
			{
				k = A->matbeg[i];
				st = k + A->matcnt[i];
				for (; k < st; k++)
				{
					if (include_logicals != 0)
					{
						ind[beg[A->matind[k]]] = i;
					}
					else
					{
						ind[beg[A->matind[k]]] = inv_structmap[i];
					}
					EGlpNumCopy (val[beg[A->matind[k]]], A->matval[k]);
					beg[A->matind[k]]++;
				}
			}
		}

		for (i = 0, k = 0; i < lp->nrows; i++)
		{
			beg[i] = k;
			k += cnt[i];
		}
	}
CLEANUP:

	if (rval)
	{
		ILLlp_rows_clear (lprows);
	}
	ILL_IFFREE (hit, char);
	ILL_IFFREE (inv_structmap, int);

	EG_RETURN (rval);
}

void ILLlp_rows_clear (
	ILLlp_rows * lprows)
{
	if (lprows != NULL)
	{
		ILL_IFFREE (lprows->rowbeg, int);
		ILL_IFFREE (lprows->rowcnt, int);
		ILL_IFFREE (lprows->rowind, int);

		EGlpNumFreeArray (lprows->rowval);
	}
}

static int wr_line (
	ILLlpdata * lp,
	const char *format,
	va_list argptr)
{
	char buffer[ILL_namebufsize];
	int rval = 0;

	rval = vsprintf (buffer, format, argptr);
	if (rval > 0)
	{
                /* Bico -- OPTERON DEBUGGING 051005  */
                /* Replaced ILLstring_report by the explicit call to */
                /* fprintf.                                          */
                /*rval = fprintf (lp->reporter.dest, buffer);
                if (rval < 0) rval = 1;
                else          rval = 0;
								*/
		/* daespino -- BACK to ILLstring_report to support compresed files 090909
		 * */
		rval = ILLstring_report(buffer, &lp->reporter);
	}
	return rval;
}

int ILLprint_report (
	ILLlpdata * lp,
	const char *format,
	...)
{
	va_list marker;
	int rval = 0;

	va_start (marker, format);		/* ANSI style */
	rval = wr_line (lp, format, marker);
	va_end (marker);							/* Reset variable arguments.      */
	return rval;
}
