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

/* $RCSfile: simplex_EGLPNUM_TYPENAME.h,v $ $Revision: 1.3 $ $Date: 2003/11/05 16:57:39 $ */

#ifndef EGLPNUM_TYPENAME___SIMPLEX_H
#define EGLPNUM_TYPENAME___SIMPLEX_H

struct itcnt_t;

#include "lpdata_EGLPNUM_TYPENAME.h"
#include "basicdefs.h"

typedef struct EGLPNUM_TYPENAME_param_info
{
	int origalgo;
	int pphaseI;
	int pphaseII;
	int dphaseI;
	int dphaseII;
	int p_strategy;
	int d_strategy;
}
EGLPNUM_TYPENAME_param_info;

typedef struct EGLPNUM_TYPENAME_iter_info
{
	int newphase;
	int nextphase;
	int nextstep;
	int sdisplay;
	int itercnt;
	int solstatus;
	int curtime;
	int rounds;
	int chkobj;
	int nosolve;
	int noprog;
	int inner;
	int algorithm;
	int resumeid;
	int pricetype;
	int n_restart;
	int n_pivot_fail;
	EGLPNUM_TYPE prevobj;
	EGLPNUM_TYPE objtol;
	EGLPNUM_TYPENAME_param_info oldinfo;
}
EGLPNUM_TYPENAME_iter_info;

void EGLPNUM_TYPENAME_ILLsimplex_init_lpinfo ( EGLPNUM_TYPENAME_lpinfo * lp),
		EGLPNUM_TYPENAME_ILLsimplex_free_lpinfo ( EGLPNUM_TYPENAME_lpinfo * lp),
		EGLPNUM_TYPENAME_ILLsimplex_load_lpinfo ( EGLPNUM_TYPENAME_ILLlpdata * qslp, EGLPNUM_TYPENAME_lpinfo * lp),
		EGLPNUM_TYPENAME_ILLsimplex_set_bound ( EGLPNUM_TYPENAME_lpinfo * lp, const EGLPNUM_TYPE * objbound, int sense);
void EGLPNUM_TYPENAME_free_internal_lpinfo ( EGLPNUM_TYPENAME_lpinfo * lp);
void EGLPNUM_TYPENAME_init_internal_lpinfo ( EGLPNUM_TYPENAME_lpinfo * lp);
int EGLPNUM_TYPENAME_build_internal_lpinfo ( EGLPNUM_TYPENAME_lpinfo * lp);
int EGLPNUM_TYPENAME_ILLsimplex_retest_psolution ( EGLPNUM_TYPENAME_lpinfo * lp, EGLPNUM_TYPENAME_price_info * p, int phase,
			EGLPNUM_TYPENAME_feas_info * fs),
		EGLPNUM_TYPENAME_ILLsimplex_retest_dsolution ( EGLPNUM_TYPENAME_lpinfo * lp, EGLPNUM_TYPENAME_price_info * p, int phase,
			EGLPNUM_TYPENAME_feas_info * fs),
		EGLPNUM_TYPENAME_ILLsimplex_solution ( EGLPNUM_TYPENAME_lpinfo * lp, EGLPNUM_TYPE * xz, EGLPNUM_TYPE * piz,
			EGLPNUM_TYPE * dz, EGLPNUM_TYPE * objval),
		EGLPNUM_TYPENAME_ILLsimplex_infcertificate ( EGLPNUM_TYPENAME_lpinfo * lp, EGLPNUM_TYPE * pi),
		EGLPNUM_TYPENAME_ILLsimplex ( EGLPNUM_TYPENAME_lpinfo * lp, int algorithm, EGLPNUM_TYPENAME_ILLlp_basis * B, 
			EGLPNUM_TYPENAME_price_info * pinf, int *sol_status, int sdisplay, itcnt_t* itcnt),
		EGLPNUM_TYPENAME_ILLsimplex_pivotin ( EGLPNUM_TYPENAME_lpinfo * lp, EGLPNUM_TYPENAME_price_info * pinf, int rcnt, 
			int *rlist, int pivot_opt, int *basis_mod);

#endif /* EGLPNUM_TYPENAME___SIMPLEX_H */
