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

/* RCSINFO $Id: factor_EGLPNUM_TYPENAME.h,v 1.3 2003/11/05 16:57:39 meven Exp $ */
#ifndef EGLPNUM_TYPENAME___QS_FACTOR_H_
#define EGLPNUM_TYPENAME___QS_FACTOR_H_

#include "basicdefs.h"
#include "dstruct_EGLPNUM_TYPENAME.h"

typedef char EGLPNUM_TYPENAME_QSbool;

typedef struct EGLPNUM_TYPENAME_uc_info
{
	int cbeg;
	int nzcnt;
	int next;
	int prev;
	int delay;
}
EGLPNUM_TYPENAME_uc_info;

typedef struct EGLPNUM_TYPENAME_ur_info
{
	EGLPNUM_TYPE max;
	int rbeg;
	int nzcnt;
	int pivcnt;
	int next;
	int prev;
	int delay;
}
EGLPNUM_TYPENAME_ur_info;

typedef struct EGLPNUM_TYPENAME_lc_info
{
	int cbeg;
	int nzcnt;
	int c;
	int crank;
	int delay;
}
EGLPNUM_TYPENAME_lc_info;

typedef struct EGLPNUM_TYPENAME_lr_info
{
	int rbeg;
	int nzcnt;
	int r;
	int rrank;
	int delay;
}
EGLPNUM_TYPENAME_lr_info;

typedef struct EGLPNUM_TYPENAME_er_info
{
	int rbeg;
	int nzcnt;
	int r;
}
EGLPNUM_TYPENAME_er_info;

typedef struct EGLPNUM_TYPENAME_factor_work
{
	int max_k;
	EGLPNUM_TYPE fzero_tol;
	EGLPNUM_TYPE szero_tol;
	EGLPNUM_TYPE partial_tol;
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

	EGLPNUM_TYPE maxelem_orig;
	int nzcnt_orig;
	EGLPNUM_TYPE maxelem_factor;
	int nzcnt_factor;
	EGLPNUM_TYPE maxelem_cur;
	int nzcnt_cur;

	EGLPNUM_TYPE partial_cur;

	int dim;
	int stage;
	int nstages;
	int etacnt;
	EGLPNUM_TYPE *work_coef;
	int *work_indx;
	EGLPNUM_TYPENAME_uc_info *uc_inf;
	EGLPNUM_TYPENAME_ur_info *ur_inf;
	EGLPNUM_TYPENAME_lc_info *lc_inf;
	EGLPNUM_TYPENAME_lr_info *lr_inf;
	EGLPNUM_TYPENAME_er_info *er_inf;
	int *ucindx;									/* row index for column data */
	int *ucrind;									/* index of column in row data */
	EGLPNUM_TYPE *uccoef;						/* coefficient for column data */
	int *urindx;									/* col index for row data */
	int *urcind;									/* index of row in column data */
	EGLPNUM_TYPE *urcoef;						/* coefficient for row data */
	int *lcindx;									/* row index for L data */
	EGLPNUM_TYPE *lccoef;						/* coefficient for L row data */
	int *lrindx;									/* col index for L data */
	EGLPNUM_TYPE *lrcoef;						/* coefficient for L col data */
	int *erindx;									/* col index for eta data */
	EGLPNUM_TYPE *ercoef;						/* coefficient for eta data */
	int *rperm;
	int *rrank;
	int *cperm;
	int *crank;
	EGLPNUM_TYPENAME_svector xtmp;
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

	EGLPNUM_TYPE *dmat;
	int drows;
	int dcols;
	int dense_base;
}
EGLPNUM_TYPENAME_factor_work;

void EGLPNUM_TYPENAME_ILLfactor_init_factor_work (
	EGLPNUM_TYPENAME_factor_work * f),
  EGLPNUM_TYPENAME_ILLfactor_free_factor_work (
	EGLPNUM_TYPENAME_factor_work * f),
  EGLPNUM_TYPENAME_ILLfactor_ftran (
	EGLPNUM_TYPENAME_factor_work * f,
	EGLPNUM_TYPENAME_svector * a,
	EGLPNUM_TYPENAME_svector * x),
  EGLPNUM_TYPENAME_ILLfactor_ftran_update (
	EGLPNUM_TYPENAME_factor_work * f,
	EGLPNUM_TYPENAME_svector * a,
	EGLPNUM_TYPENAME_svector * upd,
	EGLPNUM_TYPENAME_svector * x),
  EGLPNUM_TYPENAME_ILLfactor_btran (
	EGLPNUM_TYPENAME_factor_work * f,
	EGLPNUM_TYPENAME_svector * a,
	EGLPNUM_TYPENAME_svector * x);

int EGLPNUM_TYPENAME_ILLfactor_create_factor_work (
	EGLPNUM_TYPENAME_factor_work * f,
	int dim),
  EGLPNUM_TYPENAME_ILLfactor_set_factor_iparam (
	EGLPNUM_TYPENAME_factor_work * f,
	int param,
	int val),
  EGLPNUM_TYPENAME_ILLfactor_set_factor_dparam (
	EGLPNUM_TYPENAME_factor_work * f,
	int param,
	EGLPNUM_TYPE val),
  EGLPNUM_TYPENAME_ILLfactor (
	EGLPNUM_TYPENAME_factor_work * f,
	int *basis,
	int *cbeg,
	int *clen,
	int *cindx,
	EGLPNUM_TYPE * ccoef,
	int *p_nsing,
	int **p_singr,
	int **p_singc),
  EGLPNUM_TYPENAME_ILLfactor_update (
	EGLPNUM_TYPENAME_factor_work * f,
	EGLPNUM_TYPENAME_svector * a,
	int col,
	int *p_refact);

#endif /* EGLPNUM_TYPENAME___QS_FACTOR_H_ */
