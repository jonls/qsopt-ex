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

/* RCSINFO $Id: basis.h,v 1.3 2003/11/05 16:57:39 meven Exp $ */
#ifndef __BASIS_H
#define __BASIS_H

#include "config.h"
#include "dstruct.h"
#include "lpdefs.h"
#include "lpdata.h"

#if 0
#if EGLPNUM_TYPE != DBL_TYPE && EGLPNUM_TYPE != LDBL_TYPE
extern EGlpNum_t CB_PRI_RLIMIT;	/* = 0.25 */
extern EGlpNum_t CB_INF_RATIO;	/* = 10.0 */
extern EGlpNum_t CB_EPS;				/* = 0.001 */
#endif
#endif

typedef struct var_data
{
	int nartif;
	int nslacks;
	int nfree;
	int nbndone;
	int nbounded;
	int nfixed;
	EGlpNum_t cmax;
}
var_data;

void ILLbasis_init_vardata (
	var_data * vd);
void ILLbasis_clear_vardata (
	var_data * vd);

int ILLbasis_build_basisinfo (
	lpinfo * lp),
  ILLbasis_get_initial (
	lpinfo * lp,
	int algorithm),
  ILLbasis_get_cinitial (
	lpinfo * lp,
	int algorithm),
  ILLbasis_load (
	lpinfo * lp,
	ILLlp_basis * B),
  ILLbasis_tableau_row (
	lpinfo * lp,
	int row,
	EGlpNum_t * brow,
	EGlpNum_t * trow,
	EGlpNum_t * rhs,
	int strict),
  ILLbasis_factor (
	lpinfo * lp,
	int *singular),
  ILLbasis_refactor (
	lpinfo * lp),
  ILLbasis_update (
	lpinfo * lp,
	svector * y,
	int lindex,
	int *refactor,
	int *singular);

void ILLbasis_column_solve (
	lpinfo * lp,
	svector * rhs,
	svector * soln),
  ILLbasis_column_solve_update (
	lpinfo * lp,
	svector * rhs,
	svector * upd,
	svector * soln),
  ILLbasis_row_solve (
	lpinfo * lp,
	svector * rhs,
	svector * soln),
  ILLbasis_free_basisinfo (
	lpinfo * lp),
  ILLbasis_free_fbasisinfo (
	lpinfo * lp),
  ILLbasis_init_basisinfo (
	lpinfo * lp);

#endif /* __BASIS_H */
