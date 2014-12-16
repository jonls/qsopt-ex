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

/* RCSINFO $Id: basis_EGLPNUM_TYPENAME.h,v 1.3 2003/11/05 16:57:39 meven Exp $ */
#ifndef EGLPNUM_TYPENAME___BASIS_H
#define EGLPNUM_TYPENAME___BASIS_H

#include "dstruct_EGLPNUM_TYPENAME.h"
#include "lpdefs_EGLPNUM_TYPENAME.h"
#include "lpdata_EGLPNUM_TYPENAME.h"


typedef struct EGLPNUM_TYPENAME_var_data
{
	int nartif;
	int nslacks;
	int nfree;
	int nbndone;
	int nbounded;
	int nfixed;
	EGLPNUM_TYPE cmax;
}
EGLPNUM_TYPENAME_var_data;

void EGLPNUM_TYPENAME_ILLbasis_init_vardata (
	EGLPNUM_TYPENAME_var_data * vd);
void EGLPNUM_TYPENAME_ILLbasis_clear_vardata (
	EGLPNUM_TYPENAME_var_data * vd);

int EGLPNUM_TYPENAME_ILLbasis_build_basisinfo (
	EGLPNUM_TYPENAME_lpinfo * lp),
  EGLPNUM_TYPENAME_ILLbasis_get_initial (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int algorithm),
  EGLPNUM_TYPENAME_ILLbasis_get_cinitial (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int algorithm),
  EGLPNUM_TYPENAME_ILLbasis_load (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_ILLlp_basis * B),
  EGLPNUM_TYPENAME_ILLbasis_tableau_row (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int row,
	EGLPNUM_TYPE * brow,
	EGLPNUM_TYPE * trow,
	EGLPNUM_TYPE * rhs,
	int strict),
  EGLPNUM_TYPENAME_ILLbasis_factor (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int *singular),
  EGLPNUM_TYPENAME_ILLbasis_refactor (
	EGLPNUM_TYPENAME_lpinfo * lp),
  EGLPNUM_TYPENAME_ILLbasis_update (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_svector * y,
	int lindex,
	int *refactor,
	int *singular);

void EGLPNUM_TYPENAME_ILLbasis_column_solve (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_svector * rhs,
	EGLPNUM_TYPENAME_svector * soln),
  EGLPNUM_TYPENAME_ILLbasis_column_solve_update (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_svector * rhs,
	EGLPNUM_TYPENAME_svector * upd,
	EGLPNUM_TYPENAME_svector * soln),
  EGLPNUM_TYPENAME_ILLbasis_row_solve (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_svector * rhs,
	EGLPNUM_TYPENAME_svector * soln),
  EGLPNUM_TYPENAME_ILLbasis_free_basisinfo (
	EGLPNUM_TYPENAME_lpinfo * lp),
  EGLPNUM_TYPENAME_ILLbasis_free_fbasisinfo (
	EGLPNUM_TYPENAME_lpinfo * lp),
  EGLPNUM_TYPENAME_ILLbasis_init_basisinfo (
	EGLPNUM_TYPENAME_lpinfo * lp);

#endif /* EGLPNUM_TYPENAME___BASIS_H */
