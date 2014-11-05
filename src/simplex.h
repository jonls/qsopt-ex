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

/* $RCSfile: simplex.h,v $ $Revision: 1.3 $ $Date: 2003/11/05 16:57:39 $ */

#ifndef __SIMPLEX_H
#define __SIMPLEX_H

struct itcnt_t;
#include "config.h"
#include "lpdata.h"
#include "basicdefs.h"
typedef struct param_info
{
	int origalgo;
	int pphaseI;
	int pphaseII;
	int dphaseI;
	int dphaseII;
	int p_strategy;
	int d_strategy;
}
param_info;

typedef struct iter_info
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
	EGlpNum_t prevobj;
	EGlpNum_t objtol;
	param_info oldinfo;
}
iter_info;

void ILLsimplex_init_lpinfo ( lpinfo * lp),
		ILLsimplex_free_lpinfo ( lpinfo * lp),
		ILLsimplex_load_lpinfo ( ILLlpdata * qslp, lpinfo * lp),
		ILLsimplex_set_bound ( lpinfo * lp, const EGlpNum_t * objbound, int sense);
void free_internal_lpinfo ( lpinfo * lp);
void init_internal_lpinfo ( lpinfo * lp);
int build_internal_lpinfo ( lpinfo * lp);
int ILLsimplex_retest_psolution ( lpinfo * lp, price_info * p, int phase,
			feas_info * fs),
		ILLsimplex_retest_dsolution ( lpinfo * lp, price_info * p, int phase,
			feas_info * fs),
		ILLsimplex_solution ( lpinfo * lp, EGlpNum_t * xz, EGlpNum_t * piz,
			EGlpNum_t * dz, EGlpNum_t * objval),
		ILLsimplex_infcertificate ( lpinfo * lp, EGlpNum_t * pi),
		ILLsimplex ( lpinfo * lp, int algorithm, ILLlp_basis * B, 
			price_info * pinf, int *sol_status, int sdisplay, itcnt_t* itcnt),
		ILLsimplex_pivotin ( lpinfo * lp, price_info * pinf, int rcnt, 
			int *rlist, int pivot_opt, int *basis_mod);

#endif /* __SIMPLEX_H */
