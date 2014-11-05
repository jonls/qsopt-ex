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

/* RCSINFO $Id: factor.h,v 1.3 2003/11/05 16:57:39 meven Exp $ */
#ifndef __QS_FACTOR_H_
#define __QS_FACTOR_H_
#include "basicdefs.h"
#include "qs_config.h"
#include "dstruct.h"

typedef char QSbool;

typedef struct uc_info
{
	int cbeg;
	int nzcnt;
	int next;
	int prev;
	int delay;
}
uc_info;

typedef struct ur_info
{
	EGlpNum_t max;
	int rbeg;
	int nzcnt;
	int pivcnt;
	int next;
	int prev;
	int delay;
}
ur_info;

typedef struct lc_info
{
	int cbeg;
	int nzcnt;
	int c;
	int crank;
	int delay;
}
lc_info;

typedef struct lr_info
{
	int rbeg;
	int nzcnt;
	int r;
	int rrank;
	int delay;
}
lr_info;

typedef struct er_info
{
	int rbeg;
	int nzcnt;
	int r;
}
er_info;

typedef struct factor_work
{
	int max_k;
	EGlpNum_t fzero_tol;
	EGlpNum_t szero_tol;
	EGlpNum_t partial_tol;
	double ur_space_mul;
	double uc_space_mul;
	double lc_space_mul;
	double lr_space_mul;
	double er_space_mul;
	double grow_mul;
	int p;
	int etamax;
	double minmult;
	double maxmult;
	double updmaxmult;
	double dense_fract;
	int dense_min;

	EGlpNum_t maxelem_orig;
	int nzcnt_orig;
	EGlpNum_t maxelem_factor;
	int nzcnt_factor;
	EGlpNum_t maxelem_cur;
	int nzcnt_cur;

	EGlpNum_t partial_cur;

	int dim;
	int stage;
	int nstages;
	int etacnt;
	EGlpNum_t *work_coef;
	int *work_indx;
	uc_info *uc_inf;
	ur_info *ur_inf;
	lc_info *lc_inf;
	lr_info *lr_inf;
	er_info *er_inf;
	int *ucindx;									/* row index for column data */
	int *ucrind;									/* index of column in row data */
	EGlpNum_t *uccoef;						/* coefficient for column data */
	int *urindx;									/* col index for row data */
	int *urcind;									/* index of row in column data */
	EGlpNum_t *urcoef;						/* coefficient for row data */
	int *lcindx;									/* row index for L data */
	EGlpNum_t *lccoef;						/* coefficient for L row data */
	int *lrindx;									/* col index for L data */
	EGlpNum_t *lrcoef;						/* coefficient for L col data */
	int *erindx;									/* col index for eta data */
	EGlpNum_t *ercoef;						/* coefficient for eta data */
	int *rperm;
	int *rrank;
	int *cperm;
	int *crank;
	svector xtmp;
	int ur_freebeg;
	int ur_space;
	int uc_freebeg;
	int uc_space;
	int lc_freebeg;
	int lc_space;
	int lr_freebeg;
	int lr_space;
	int er_freebeg;
	int er_space;

	int *p_nsing;
	int **p_singr;
	int **p_singc;

	EGlpNum_t *dmat;
	int drows;
	int dcols;
	int dense_base;
}
factor_work;

void ILLfactor_init_factor_work (
	factor_work * f),
  ILLfactor_free_factor_work (
	factor_work * f),
  ILLfactor_ftran (
	factor_work * f,
	svector * a,
	svector * x),
  ILLfactor_ftran_update (
	factor_work * f,
	svector * a,
	svector * upd,
	svector * x),
  ILLfactor_btran (
	factor_work * f,
	svector * a,
	svector * x);

int ILLfactor_create_factor_work (
	factor_work * f,
	int dim),
  ILLfactor_set_factor_iparam (
	factor_work * f,
	int param,
	int val),
  ILLfactor_set_factor_dparam (
	factor_work * f,
	int param,
	EGlpNum_t val),
  ILLfactor (
	factor_work * f,
	int *basis,
	int *cbeg,
	int *clen,
	int *cindx,
	EGlpNum_t * ccoef,
	int *p_nsing,
	int **p_singr,
	int **p_singc),
  ILLfactor_update (
	factor_work * f,
	svector * a,
	int col,
	int *p_refact);

#endif /* __QS_FACTOR_H_ */
