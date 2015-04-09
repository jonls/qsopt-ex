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
/*    int ILLlpdata_buildrows (EGLPNUM_TYPENAME_ILLlpdata *lp, int **rowbeg, int **rowcnt,     */
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
/*    void EGLPNUM_TYPENAME_ILLlpdata_init (EGLPNUM_TYPENAME_ILLlpdata *lp)                                   */
/*    void EGLPNUM_TYPENAME_ILLlpdata_free (EGLPNUM_TYPENAME_ILLlpdata *lp)                                   */
/*                                                                          */
/*    void EGLPNUM_TYPENAME_ILLlp_basis_init (EGLPNUM_TYPENAME_ILLlp_basis *B)                                */
/*    void EGLPNUM_TYPENAME_ILLlp_basis_free (EGLPNUM_TYPENAME_ILLlp_basis *B)                                */
/*    int EGLPNUM_TYPENAME_ILLlp_basis_alloc (EGLPNUM_TYPENAME_ILLlp_basis *B, int nstruct, int nrows)        */
/*                                                                          */
/*    void EGLPNUM_TYPENAME_ILLlp_cache_init (EGLPNUM_TYPENAME_ILLlp_cache *C)                                */
/*    void EGLPNUM_TYPENAME_ILLlp_cache_free (EGLPNUM_TYPENAME_ILLlp_cache *C)                                */
/*    int EGLPNUM_TYPENAME_ILLlp_cache_alloc (EGLPNUM_TYPENAME_ILLlp_cache *C, int nstruct, int nrows)        */
/*                                                                          */
/*    void EGLPNUM_TYPENAME_ILLlp_sinfo_init (EGLPNUM_TYPENAME_ILLlp_sinfo *sinfo)                            */
/*    void EGLPNUM_TYPENAME_ILLlp_sinfo_free (EGLPNUM_TYPENAME_ILLlp_sinfo *sinfo)                            */
/*                                                                          */
/*    int EGLPNUM_TYPENAME_ILLlp_rows_init(EGLPNUM_TYPENAME_ILLlp_rows *lprows, EGLPNUM_TYPENAME_ILLlpdata *lp,                */
/*                                           int include_logicals)          */
/*                                                                          */
/****************************************************************************/

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <stdlib.h>
#include <stdarg.h>

#include "qs_config.h"
#include "logging-private.h"

#include "allocrus.h"
#include "eg_lpnum.h"
#include "eg_io.h"
#include "except.h"

#include "lpdata_EGLPNUM_TYPENAME.h"
#include "qstruct_EGLPNUM_TYPENAME.h"
#include "qsopt_EGLPNUM_TYPENAME.h"
#include "lp_EGLPNUM_TYPENAME.h"
#include "mps_EGLPNUM_TYPENAME.h"
#include "rawlp_EGLPNUM_TYPENAME.h"

//static int TRACE = 0;

EGLPNUM_TYPE EGLPNUM_TYPENAME_PARAM_IBASIS_RPIVOT;
EGLPNUM_TYPE EGLPNUM_TYPENAME_PARAM_IBASIS_RTRIANG;
EGLPNUM_TYPE EGLPNUM_TYPENAME_PARAM_MIN_DNORM;
EGLPNUM_TYPE EGLPNUM_TYPENAME_PFEAS_TOLER;
EGLPNUM_TYPE EGLPNUM_TYPENAME_BD_TOLER;
EGLPNUM_TYPE EGLPNUM_TYPENAME_DFEAS_TOLER;
EGLPNUM_TYPE EGLPNUM_TYPENAME_PIVOT_TOLER;
EGLPNUM_TYPE EGLPNUM_TYPENAME_SZERO_TOLER;
EGLPNUM_TYPE EGLPNUM_TYPENAME_PIVZ_TOLER;
EGLPNUM_TYPE EGLPNUM_TYPENAME_OBJBND_TOLER;
EGLPNUM_TYPE EGLPNUM_TYPENAME_DBNDPIV_TOLER;
EGLPNUM_TYPE EGLPNUM_TYPENAME_DBNDPIV_RATIO;
EGLPNUM_TYPE EGLPNUM_TYPENAME_ALTPIV_TOLER;
//EGLPNUM_TYPE DJZERO_TOLER;
EGLPNUM_TYPE EGLPNUM_TYPENAME_PROGRESS_ZERO;				/*   1e-7 */
EGLPNUM_TYPE EGLPNUM_TYPENAME_PROGRESS_THRESH;			/*   1e-5 */
EGLPNUM_TYPE EGLPNUM_TYPENAME_CB_EPS;
EGLPNUM_TYPE EGLPNUM_TYPENAME_CB_INF_RATIO;
EGLPNUM_TYPE EGLPNUM_TYPENAME_CB_PRI_RLIMIT;
EGLPNUM_TYPE EGLPNUM_TYPENAME_ILL_MAXDOUBLE;
EGLPNUM_TYPE EGLPNUM_TYPENAME_ILL_MINDOUBLE;

/* ========================================================================= */
int EGLPNUM_TYPENAME___QSEX_SETUP = 0;
/* ========================================================================= */
void EGLPNUM_TYPENAME_ILLstart ( void)
{
	if (EGLPNUM_TYPENAME___QSEX_SETUP)
		return;
	EGLPNUM_TYPENAME_EGlpNumInitVar (EGLPNUM_TYPENAME_PARAM_IBASIS_RPIVOT);
	EGLPNUM_TYPENAME_EGlpNumInitVar (EGLPNUM_TYPENAME_PARAM_IBASIS_RTRIANG);
	EGLPNUM_TYPENAME_EGlpNumInitVar (EGLPNUM_TYPENAME_PARAM_MIN_DNORM);
	EGLPNUM_TYPENAME_EGlpNumInitVar (EGLPNUM_TYPENAME_PFEAS_TOLER);
	EGLPNUM_TYPENAME_EGlpNumInitVar (EGLPNUM_TYPENAME_BD_TOLER);
	EGLPNUM_TYPENAME_EGlpNumInitVar (EGLPNUM_TYPENAME_DFEAS_TOLER);
	EGLPNUM_TYPENAME_EGlpNumInitVar (EGLPNUM_TYPENAME_PIVOT_TOLER);
	EGLPNUM_TYPENAME_EGlpNumInitVar (EGLPNUM_TYPENAME_SZERO_TOLER);
	EGLPNUM_TYPENAME_EGlpNumInitVar (EGLPNUM_TYPENAME_PIVZ_TOLER);
	EGLPNUM_TYPENAME_EGlpNumInitVar (EGLPNUM_TYPENAME_OBJBND_TOLER);
	EGLPNUM_TYPENAME_EGlpNumInitVar (EGLPNUM_TYPENAME_DBNDPIV_TOLER);
	EGLPNUM_TYPENAME_EGlpNumInitVar (EGLPNUM_TYPENAME_DBNDPIV_RATIO);
	EGLPNUM_TYPENAME_EGlpNumInitVar (EGLPNUM_TYPENAME_ALTPIV_TOLER);
	//EGLPNUM_TYPENAME_EGlpNumInitVar (DJZERO_TOLER);
	EGLPNUM_TYPENAME_EGlpNumInitVar (EGLPNUM_TYPENAME_PROGRESS_ZERO);	/*            1e-7 */
	EGLPNUM_TYPENAME_EGlpNumInitVar (EGLPNUM_TYPENAME_PROGRESS_THRESH);	/*          1e-5 */
	EGLPNUM_TYPENAME_EGlpNumInitVar (EGLPNUM_TYPENAME_CB_PRI_RLIMIT);
	EGLPNUM_TYPENAME_EGlpNumInitVar (EGLPNUM_TYPENAME_CB_INF_RATIO);
	EGLPNUM_TYPENAME_EGlpNumInitVar (EGLPNUM_TYPENAME_CB_EPS);
	EGLPNUM_TYPENAME_EGlpNumInitVar (EGLPNUM_TYPENAME_ILL_MAXDOUBLE);
	EGLPNUM_TYPENAME_EGlpNumInitVar (EGLPNUM_TYPENAME_ILL_MINDOUBLE);
	/* parameters that do depend on the tolerance to zero */
	EGLPNUM_TYPENAME_EGlpNumSet (EGLPNUM_TYPENAME_PARAM_MIN_DNORM, 4.5036e-9);
	EGLPNUM_TYPENAME_EGlpNumMultTo (EGLPNUM_TYPENAME_PARAM_MIN_DNORM, EGLPNUM_TYPENAME_epsLpNum);
	EGLPNUM_TYPENAME_EGlpNumSet (EGLPNUM_TYPENAME_PFEAS_TOLER, 4.5036e9);
	EGLPNUM_TYPENAME_EGlpNumMultTo (EGLPNUM_TYPENAME_PFEAS_TOLER, EGLPNUM_TYPENAME_epsLpNum);
	EGLPNUM_TYPENAME_EGlpNumSet (EGLPNUM_TYPENAME_BD_TOLER, 4.5036e8);
	EGLPNUM_TYPENAME_EGlpNumMultTo (EGLPNUM_TYPENAME_BD_TOLER, EGLPNUM_TYPENAME_epsLpNum);
	EGLPNUM_TYPENAME_EGlpNumSet (EGLPNUM_TYPENAME_DFEAS_TOLER, 4.5036e9);
	EGLPNUM_TYPENAME_EGlpNumMultTo (EGLPNUM_TYPENAME_DFEAS_TOLER, EGLPNUM_TYPENAME_epsLpNum);
	EGLPNUM_TYPENAME_EGlpNumSet (EGLPNUM_TYPENAME_PIVOT_TOLER, 4.5036e5);
	EGLPNUM_TYPENAME_EGlpNumMultTo (EGLPNUM_TYPENAME_PIVOT_TOLER, EGLPNUM_TYPENAME_epsLpNum);
	EGLPNUM_TYPENAME_EGlpNumSet (EGLPNUM_TYPENAME_SZERO_TOLER, 4.5036);
	EGLPNUM_TYPENAME_EGlpNumMultTo (EGLPNUM_TYPENAME_SZERO_TOLER, EGLPNUM_TYPENAME_epsLpNum);
	EGLPNUM_TYPENAME_EGlpNumSet (EGLPNUM_TYPENAME_PIVZ_TOLER, 4.5036e3);
	EGLPNUM_TYPENAME_EGlpNumMultTo (EGLPNUM_TYPENAME_PIVZ_TOLER, EGLPNUM_TYPENAME_epsLpNum);
	EGLPNUM_TYPENAME_EGlpNumSet (EGLPNUM_TYPENAME_OBJBND_TOLER, 4.5036e13);
	EGLPNUM_TYPENAME_EGlpNumMultTo (EGLPNUM_TYPENAME_OBJBND_TOLER, EGLPNUM_TYPENAME_epsLpNum);
	EGLPNUM_TYPENAME_EGlpNumSet (EGLPNUM_TYPENAME_ALTPIV_TOLER, 4.5036e7);
	EGLPNUM_TYPENAME_EGlpNumMultTo (EGLPNUM_TYPENAME_ALTPIV_TOLER, EGLPNUM_TYPENAME_epsLpNum);
	EGLPNUM_TYPENAME_EGlpNumSet (EGLPNUM_TYPENAME_PROGRESS_ZERO, 4.5036e8);
	EGLPNUM_TYPENAME_EGlpNumMultTo (EGLPNUM_TYPENAME_PROGRESS_ZERO, EGLPNUM_TYPENAME_epsLpNum);
	EGLPNUM_TYPENAME_EGlpNumSet (EGLPNUM_TYPENAME_PROGRESS_THRESH, 4.5036e10);
	EGLPNUM_TYPENAME_EGlpNumMultTo (EGLPNUM_TYPENAME_PROGRESS_THRESH, EGLPNUM_TYPENAME_epsLpNum);
#if VERBOSE_LEVEL <= DEBUG
	MESSAGE (VERBOSE_LEVEL, "Setting EGLPNUM_TYPENAME_PARAM_MIN_DNORM to %lg", EGLPNUM_TYPENAME_EGlpNumToLf (EGLPNUM_TYPENAME_PARAM_MIN_DNORM));
	MESSAGE (VERBOSE_LEVEL, "Setting EGLPNUM_TYPENAME_PFEAS_TOLER to %lg", EGLPNUM_TYPENAME_EGlpNumToLf (EGLPNUM_TYPENAME_PFEAS_TOLER));
	MESSAGE (VERBOSE_LEVEL, "Setting EGLPNUM_TYPENAME_BD_TOLER to %lg", EGLPNUM_TYPENAME_EGlpNumToLf (EGLPNUM_TYPENAME_BD_TOLER));
	MESSAGE (VERBOSE_LEVEL, "Setting EGLPNUM_TYPENAME_DFEAS_TOLER to %lg", EGLPNUM_TYPENAME_EGlpNumToLf (EGLPNUM_TYPENAME_DFEAS_TOLER));
	MESSAGE (VERBOSE_LEVEL, "Setting EGLPNUM_TYPENAME_PIVOT_TOLER to %lg", EGLPNUM_TYPENAME_EGlpNumToLf (EGLPNUM_TYPENAME_PIVOT_TOLER));
	MESSAGE (VERBOSE_LEVEL, "Setting EGLPNUM_TYPENAME_SZERO_TOLER to %lg", EGLPNUM_TYPENAME_EGlpNumToLf (EGLPNUM_TYPENAME_SZERO_TOLER));
	MESSAGE (VERBOSE_LEVEL, "Setting EGLPNUM_TYPENAME_PIVZ_TOLER to %lg", EGLPNUM_TYPENAME_EGlpNumToLf (EGLPNUM_TYPENAME_PIVZ_TOLER));
	MESSAGE (VERBOSE_LEVEL, "Setting EGLPNUM_TYPENAME_OBJBND_TOLER to %lg", EGLPNUM_TYPENAME_EGlpNumToLf (EGLPNUM_TYPENAME_OBJBND_TOLER));
	MESSAGE (VERBOSE_LEVEL, "Setting EGLPNUM_TYPENAME_ALTPIV_TOLER to %lg", EGLPNUM_TYPENAME_EGlpNumToLf (EGLPNUM_TYPENAME_ALTPIV_TOLER));
	MESSAGE (VERBOSE_LEVEL, "Setting EGLPNUM_TYPENAME_PROGRESS_ZERO to %lg", EGLPNUM_TYPENAME_EGlpNumToLf (EGLPNUM_TYPENAME_PROGRESS_ZERO));
	MESSAGE (VERBOSE_LEVEL, "Setting EGLPNUM_TYPENAME_PROGRESS_THRESH to %lg", EGLPNUM_TYPENAME_EGlpNumToLf (EGLPNUM_TYPENAME_PROGRESS_THRESH));
#endif
	/* parameters that do not depend on the tolerance to zero */
	EGLPNUM_TYPENAME_EGlpNumSet (EGLPNUM_TYPENAME_ILL_MAXDOUBLE, 1e150);
	EGLPNUM_TYPENAME_EGlpNumSet (EGLPNUM_TYPENAME_ILL_MINDOUBLE, -1e150);
	EGLPNUM_TYPENAME_EGlpNumSet (EGLPNUM_TYPENAME_PARAM_IBASIS_RPIVOT, 0.98);
	EGLPNUM_TYPENAME_EGlpNumSet (EGLPNUM_TYPENAME_PARAM_IBASIS_RTRIANG, 0.01);
	EGLPNUM_TYPENAME_EGlpNumSet (EGLPNUM_TYPENAME_DBNDPIV_TOLER, 1e-3);
	EGLPNUM_TYPENAME_EGlpNumSet (EGLPNUM_TYPENAME_DBNDPIV_RATIO, 1e-2);
	//EGLPNUM_TYPENAME_EGlpNumSet (DJZERO_TOLER, 1e-8);
	EGLPNUM_TYPENAME_EGlpNumSet (EGLPNUM_TYPENAME_CB_EPS, 0.001);
	EGLPNUM_TYPENAME_EGlpNumSet (EGLPNUM_TYPENAME_CB_INF_RATIO, 10.0);
	EGLPNUM_TYPENAME_EGlpNumSet (EGLPNUM_TYPENAME_CB_PRI_RLIMIT, 0.25);
	EGLPNUM_TYPENAME___QSEX_SETUP = 1;
}

/* ========================================================================= */
void EGLPNUM_TYPENAME_ILLchange_precision (
	void)
{
	EGLPNUM_TYPENAME_EGlpNumClearVar (EGLPNUM_TYPENAME_PFEAS_TOLER);
	EGLPNUM_TYPENAME_EGlpNumClearVar (EGLPNUM_TYPENAME_BD_TOLER);
	EGLPNUM_TYPENAME_EGlpNumClearVar (EGLPNUM_TYPENAME_DFEAS_TOLER);
	EGLPNUM_TYPENAME_EGlpNumClearVar (EGLPNUM_TYPENAME_PIVOT_TOLER);
	EGLPNUM_TYPENAME_EGlpNumClearVar (EGLPNUM_TYPENAME_SZERO_TOLER);
	EGLPNUM_TYPENAME_EGlpNumClearVar (EGLPNUM_TYPENAME_PIVZ_TOLER);
	EGLPNUM_TYPENAME_EGlpNumClearVar (EGLPNUM_TYPENAME_OBJBND_TOLER);
	EGLPNUM_TYPENAME_EGlpNumClearVar (EGLPNUM_TYPENAME_ALTPIV_TOLER);
	EGLPNUM_TYPENAME_EGlpNumClearVar (EGLPNUM_TYPENAME_PARAM_MIN_DNORM);
	EGLPNUM_TYPENAME_EGlpNumClearVar (EGLPNUM_TYPENAME_PROGRESS_ZERO);
	EGLPNUM_TYPENAME_EGlpNumClearVar (EGLPNUM_TYPENAME_PROGRESS_THRESH);
	EGLPNUM_TYPENAME_EGlpNumInitVar (EGLPNUM_TYPENAME_PROGRESS_ZERO);
	EGLPNUM_TYPENAME_EGlpNumInitVar (EGLPNUM_TYPENAME_PROGRESS_THRESH);
	EGLPNUM_TYPENAME_EGlpNumInitVar (EGLPNUM_TYPENAME_PFEAS_TOLER);
	EGLPNUM_TYPENAME_EGlpNumInitVar (EGLPNUM_TYPENAME_BD_TOLER);
	EGLPNUM_TYPENAME_EGlpNumInitVar (EGLPNUM_TYPENAME_DFEAS_TOLER);
	EGLPNUM_TYPENAME_EGlpNumInitVar (EGLPNUM_TYPENAME_PIVOT_TOLER);
	EGLPNUM_TYPENAME_EGlpNumInitVar (EGLPNUM_TYPENAME_SZERO_TOLER);
	EGLPNUM_TYPENAME_EGlpNumInitVar (EGLPNUM_TYPENAME_PIVZ_TOLER);
	EGLPNUM_TYPENAME_EGlpNumInitVar (EGLPNUM_TYPENAME_OBJBND_TOLER);
	EGLPNUM_TYPENAME_EGlpNumInitVar (EGLPNUM_TYPENAME_ALTPIV_TOLER);
	EGLPNUM_TYPENAME_EGlpNumInitVar (EGLPNUM_TYPENAME_PARAM_MIN_DNORM);
	/* parameters that do depend on the tolerance to zero */
	EGLPNUM_TYPENAME_EGlpNumSet (EGLPNUM_TYPENAME_PARAM_MIN_DNORM, 4.5036e-9);
	EGLPNUM_TYPENAME_EGlpNumMultTo (EGLPNUM_TYPENAME_PARAM_MIN_DNORM, EGLPNUM_TYPENAME_epsLpNum);
	EGLPNUM_TYPENAME_EGlpNumSet (EGLPNUM_TYPENAME_PFEAS_TOLER, 4.5036e9);
	EGLPNUM_TYPENAME_EGlpNumMultTo (EGLPNUM_TYPENAME_PFEAS_TOLER, EGLPNUM_TYPENAME_epsLpNum);
	EGLPNUM_TYPENAME_EGlpNumSet (EGLPNUM_TYPENAME_BD_TOLER, 4.5036e8);
	EGLPNUM_TYPENAME_EGlpNumMultTo (EGLPNUM_TYPENAME_BD_TOLER, EGLPNUM_TYPENAME_epsLpNum);
	EGLPNUM_TYPENAME_EGlpNumSet (EGLPNUM_TYPENAME_DFEAS_TOLER, 4.5036e9);
	EGLPNUM_TYPENAME_EGlpNumMultTo (EGLPNUM_TYPENAME_DFEAS_TOLER, EGLPNUM_TYPENAME_epsLpNum);
	EGLPNUM_TYPENAME_EGlpNumSet (EGLPNUM_TYPENAME_PIVOT_TOLER, 4.5036e5);
	EGLPNUM_TYPENAME_EGlpNumMultTo (EGLPNUM_TYPENAME_PIVOT_TOLER, EGLPNUM_TYPENAME_epsLpNum);
	EGLPNUM_TYPENAME_EGlpNumSet (EGLPNUM_TYPENAME_SZERO_TOLER, 4.5036);
	EGLPNUM_TYPENAME_EGlpNumMultTo (EGLPNUM_TYPENAME_SZERO_TOLER, EGLPNUM_TYPENAME_epsLpNum);
	EGLPNUM_TYPENAME_EGlpNumSet (EGLPNUM_TYPENAME_PIVZ_TOLER, 4.5036e3);
	EGLPNUM_TYPENAME_EGlpNumMultTo (EGLPNUM_TYPENAME_PIVZ_TOLER, EGLPNUM_TYPENAME_epsLpNum);
	EGLPNUM_TYPENAME_EGlpNumSet (EGLPNUM_TYPENAME_OBJBND_TOLER, 4.5036e13);
	EGLPNUM_TYPENAME_EGlpNumMultTo (EGLPNUM_TYPENAME_OBJBND_TOLER, EGLPNUM_TYPENAME_epsLpNum);
	EGLPNUM_TYPENAME_EGlpNumSet (EGLPNUM_TYPENAME_ALTPIV_TOLER, 4.5036e7);
	EGLPNUM_TYPENAME_EGlpNumMultTo (EGLPNUM_TYPENAME_ALTPIV_TOLER, EGLPNUM_TYPENAME_epsLpNum);
	EGLPNUM_TYPENAME_EGlpNumSet (EGLPNUM_TYPENAME_PROGRESS_ZERO, 4.5036e8);
	EGLPNUM_TYPENAME_EGlpNumMultTo (EGLPNUM_TYPENAME_PROGRESS_ZERO, EGLPNUM_TYPENAME_epsLpNum);
	EGLPNUM_TYPENAME_EGlpNumSet (EGLPNUM_TYPENAME_PROGRESS_THRESH, 4.5036e10);
	EGLPNUM_TYPENAME_EGlpNumMultTo (EGLPNUM_TYPENAME_PROGRESS_THRESH, EGLPNUM_TYPENAME_epsLpNum);
}

/* ========================================================================= */
void EGLPNUM_TYPENAME_ILLend ( void)
{
	if (!EGLPNUM_TYPENAME___QSEX_SETUP)
		return;
	EGLPNUM_TYPENAME_EGlpNumClearVar (EGLPNUM_TYPENAME_PARAM_IBASIS_RPIVOT);
	EGLPNUM_TYPENAME_EGlpNumClearVar (EGLPNUM_TYPENAME_PARAM_IBASIS_RTRIANG);
	EGLPNUM_TYPENAME_EGlpNumClearVar (EGLPNUM_TYPENAME_PARAM_MIN_DNORM);
	EGLPNUM_TYPENAME_EGlpNumClearVar (EGLPNUM_TYPENAME_PFEAS_TOLER);
	EGLPNUM_TYPENAME_EGlpNumClearVar (EGLPNUM_TYPENAME_BD_TOLER);
	EGLPNUM_TYPENAME_EGlpNumClearVar (EGLPNUM_TYPENAME_DFEAS_TOLER);
	EGLPNUM_TYPENAME_EGlpNumClearVar (EGLPNUM_TYPENAME_PIVOT_TOLER);
	EGLPNUM_TYPENAME_EGlpNumClearVar (EGLPNUM_TYPENAME_SZERO_TOLER);
	EGLPNUM_TYPENAME_EGlpNumClearVar (EGLPNUM_TYPENAME_PIVZ_TOLER);
	EGLPNUM_TYPENAME_EGlpNumClearVar (EGLPNUM_TYPENAME_OBJBND_TOLER);
	EGLPNUM_TYPENAME_EGlpNumClearVar (EGLPNUM_TYPENAME_DBNDPIV_TOLER);
	EGLPNUM_TYPENAME_EGlpNumClearVar (EGLPNUM_TYPENAME_DBNDPIV_RATIO);
	EGLPNUM_TYPENAME_EGlpNumClearVar (EGLPNUM_TYPENAME_ALTPIV_TOLER);
	//EGLPNUM_TYPENAME_EGlpNumClearVar (DJZERO_TOLER);
	EGLPNUM_TYPENAME_EGlpNumClearVar (EGLPNUM_TYPENAME_PROGRESS_ZERO);	/*            1e-7 */
	EGLPNUM_TYPENAME_EGlpNumClearVar (EGLPNUM_TYPENAME_PROGRESS_THRESH);	/*          1e-5 */
	EGLPNUM_TYPENAME_EGlpNumClearVar (EGLPNUM_TYPENAME_CB_EPS);
	EGLPNUM_TYPENAME_EGlpNumClearVar (EGLPNUM_TYPENAME_CB_INF_RATIO);
	EGLPNUM_TYPENAME_EGlpNumClearVar (EGLPNUM_TYPENAME_CB_PRI_RLIMIT);
	EGLPNUM_TYPENAME_EGlpNumClearVar (EGLPNUM_TYPENAME_ILL_MAXDOUBLE);
	EGLPNUM_TYPENAME_EGlpNumClearVar (EGLPNUM_TYPENAME_ILL_MINDOUBLE);
	EGLPNUM_TYPENAME___QSEX_SETUP = 0;
}

EGLPNUM_TYPENAME_QSdata *EGLPNUM_TYPENAME_ILLread (
	EGLPNUM_TYPENAME_qsline_reader * file,
	const char *fname,
	int isMps)
{
	int rval = 0;
	EGLPNUM_TYPENAME_QSdata *p = 0;
	EGLPNUM_TYPENAME_ILLlpdata *lp;
	EGLPNUM_TYPENAME_rawlpdata rawlp;

	ILL_FAILfalse (file != NULL, NULL);
	ILL_FAILfalse (fname != NULL, NULL);

	p = EGLPNUM_TYPENAME_QScreate_prob (fname, QS_MIN);
	ILL_CHECKnull (p, NULL);
	ILL_IFFREE (p->qslp->probname, char);

	lp = p->qslp;

	EGLPNUM_TYPENAME_ILLinit_rawlpdata (&rawlp, file->error_collector);
	EGLPNUM_TYPENAME_ILLlpdata_init (lp);

	if (isMps != 0)
	{
		rval = EGLPNUM_TYPENAME_ILLread_mps (file, fname, &rawlp);
	}
	else
	{
		rval = EGLPNUM_TYPENAME_ILLread_lp (file, fname, &rawlp);
	}
	CHECKRVALG (rval, CLEANUP);

	rval = EGLPNUM_TYPENAME_ILLrawlpdata_to_lpdata (&rawlp, lp);
	CHECKRVALG (rval, CLEANUP);

CLEANUP:
	EGLPNUM_TYPENAME_ILLfree_rawlpdata (&rawlp);
	if (rval != 0)
	{
		EGLPNUM_TYPENAME_QSfree_prob (p);
		p = 0;
	}
	return p;
}

static int EGLPNUM_TYPENAME_ILLlpdata_log(
	void *dest, const char *s)
{
	if (s != NULL) {
		QSlog("%s", s);
	}
	return 0;
}

void EGLPNUM_TYPENAME_ILLlpdata_init (
	EGLPNUM_TYPENAME_ILLlpdata * lp)
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
		lp->objsense = EGLPNUM_TYPENAME_ILL_MIN;
		lp->sense = 0;
		lp->obj = 0;
		lp->rhs = 0;
		lp->rangeval = 0;
		lp->lower = 0;
		lp->upper = 0;

		EGLPNUM_TYPENAME_ILLmatrix_init (&lp->A);
		EGLPNUM_TYPENAME_ILLmatrix_init (&lp->sos);
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

		ILLstring_reporter_init(
			&lp->reporter, EGLPNUM_TYPENAME_ILLlpdata_log, NULL);
	}
}

void EGLPNUM_TYPENAME_ILLlpdata_free (
	EGLPNUM_TYPENAME_ILLlpdata * lp)
{
	int i;

	if (lp)
	{
		ILL_IFFREE (lp->sense, char);

		EGLPNUM_TYPENAME_EGlpNumFreeArray (lp->obj);
		EGLPNUM_TYPENAME_EGlpNumFreeArray (lp->rhs);
		EGLPNUM_TYPENAME_EGlpNumFreeArray (lp->rangeval);
		EGLPNUM_TYPENAME_EGlpNumFreeArray (lp->lower);
		EGLPNUM_TYPENAME_EGlpNumFreeArray (lp->upper);
		EGLPNUM_TYPENAME_ILLmatrix_free (&lp->A);
		if (lp->rA)
		{
			EGLPNUM_TYPENAME_ILLlp_rows_clear (lp->rA);
			ILL_IFFREE (lp->rA, EGLPNUM_TYPENAME_ILLlp_rows);
		}
		ILL_IFFREE (lp->is_sos_mem, int);
		ILL_IFFREE (lp->refrowname, char);

		EGLPNUM_TYPENAME_ILLmatrix_free (&lp->sos);
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
			EGLPNUM_TYPENAME_ILLlp_sinfo_free (lp->sinfo);
			ILL_IFFREE (lp->sinfo, EGLPNUM_TYPENAME_ILLlp_sinfo);
		}
		EGLPNUM_TYPENAME_ILLlpdata_init (lp);
	}
}

void EGLPNUM_TYPENAME_ILLlp_basis_init (
	EGLPNUM_TYPENAME_ILLlp_basis * B)
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

void EGLPNUM_TYPENAME_ILLlp_basis_free (
	EGLPNUM_TYPENAME_ILLlp_basis * B)
{
	if (B)
	{
		ILL_IFFREE (B->cstat, char);
		ILL_IFFREE (B->rstat, char);

		EGLPNUM_TYPENAME_EGlpNumFreeArray (B->rownorms);
		EGLPNUM_TYPENAME_EGlpNumFreeArray (B->colnorms);
		B->nstruct = 0;
		B->nrows = 0;
	}
}

int EGLPNUM_TYPENAME_ILLlp_basis_alloc (
	EGLPNUM_TYPENAME_ILLlp_basis * B,
	int nstruct,
	int nrows)
{
	int rval = 0;

	ILL_FAILtrue (B == NULL, "EGLPNUM_TYPENAME_ILLlp_basis_alloc called without a basis");

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
		EGLPNUM_TYPENAME_ILLlp_basis_free (B);
	}

	EG_RETURN (rval);
}

void EGLPNUM_TYPENAME_ILLlp_cache_init (
	EGLPNUM_TYPENAME_ILLlp_cache * C)
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
		EGLPNUM_TYPENAME_EGlpNumZero (C->val);
	}
}

void EGLPNUM_TYPENAME_ILLlp_cache_free (
	EGLPNUM_TYPENAME_ILLlp_cache * C)
{
	if (C)
	{
		EGLPNUM_TYPENAME_EGlpNumFreeArray (C->x);
		EGLPNUM_TYPENAME_EGlpNumFreeArray (C->rc);
		EGLPNUM_TYPENAME_EGlpNumFreeArray (C->pi);
		EGLPNUM_TYPENAME_EGlpNumFreeArray (C->slack);
		C->nstruct = 0;
		C->nrows = 0;
		C->status = 0;
	}
}

int EGLPNUM_TYPENAME_ILLlp_cache_alloc (
	EGLPNUM_TYPENAME_ILLlp_cache * C,
	int nstruct,
	int nrows)
{
	int rval = 0;

	ILL_FAILtrue (C == NULL, "EGLPNUM_TYPENAME_ILLlp_cache_alloc called without a cache");

	C->nstruct = nstruct;
	C->nrows = nrows;

	if (nstruct > 0)
	{
		C->x = EGLPNUM_TYPENAME_EGlpNumAllocArray (nstruct);
		C->rc = EGLPNUM_TYPENAME_EGlpNumAllocArray (nstruct);
	}

	if (nrows > 0)
	{
		C->pi = EGLPNUM_TYPENAME_EGlpNumAllocArray (nrows);
		C->slack = EGLPNUM_TYPENAME_EGlpNumAllocArray (nrows);
	}

CLEANUP:

	if (rval)
	{
		EGLPNUM_TYPENAME_ILLlp_cache_free (C);
	}

	EG_RETURN (rval);
}


int EGLPNUM_TYPENAME_ILLlp_rows_init (
	EGLPNUM_TYPENAME_ILLlp_rows * lprows,
	EGLPNUM_TYPENAME_ILLlpdata * lp,
	int include_logicals)
{
	int rval = 0;
	int i, k, st;
	int *beg, *cnt, *ind;
	EGLPNUM_TYPE *val;
	EGLPNUM_TYPENAME_ILLmatrix *A;
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

				lprows->rowval = EGLPNUM_TYPENAME_EGlpNumAllocArray (lp->nzcount);
			}
			else
			{
				ILL_SAFE_MALLOC (lprows->rowind, lp->nzcount - lp->nrows, int);

				lprows->rowval = EGLPNUM_TYPENAME_EGlpNumAllocArray (lp->nzcount - lp->nrows);
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
					EGLPNUM_TYPENAME_EGlpNumCopy (val[beg[A->matind[k]]], A->matval[k]);
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
		EGLPNUM_TYPENAME_ILLlp_rows_clear (lprows);
	}
	ILL_IFFREE (hit, char);
	ILL_IFFREE (inv_structmap, int);

	EG_RETURN (rval);
}

void EGLPNUM_TYPENAME_ILLlp_rows_clear (
	EGLPNUM_TYPENAME_ILLlp_rows * lprows)
{
	if (lprows != NULL)
	{
		ILL_IFFREE (lprows->rowbeg, int);
		ILL_IFFREE (lprows->rowcnt, int);
		ILL_IFFREE (lprows->rowind, int);

		EGLPNUM_TYPENAME_EGlpNumFreeArray (lprows->rowval);
	}
}

static int wr_line (
	EGLPNUM_TYPENAME_ILLlpdata * lp,
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

int EGLPNUM_TYPENAME_ILLprint_report (
	EGLPNUM_TYPENAME_ILLlpdata * lp,
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
