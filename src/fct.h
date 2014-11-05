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

/* RCSINFO $Id: fct.h,v 1.3 2003/11/05 16:57:39 meven Exp $ */
#ifndef __FUNCTIONS_H
#define __FUNCTIONS_H
#include "basicdefs.h"
int ILLfct_compute_zA (
	lpinfo * lp,
	svector * z,
	svector * zA),
  ILLfct_compute_wz (
	lpinfo * lp,
	EGlpNum_t * wz),
  ILLfct_adjust_viol_bounds (
	lpinfo * lp),
  ILLfct_perturb_bounds (
	lpinfo * lp),
  ILLfct_perturb_phaseI_bounds (
	lpinfo * lp),
  ILLfct_bound_shift (
	lpinfo * lp,
	int col,
	int bndtype,
	EGlpNum_t newbnd),
  ILLfct_adjust_viol_coefs (
	lpinfo * lp),
  ILLfct_perturb_coefs (
	lpinfo * lp),
  ILLfct_coef_shift (
	lpinfo * lp,
	int col,
	EGlpNum_t newcoef),
  ILLfct_test_pivot (
	lpinfo * lp,
	int indx,
	int indxtype,
	EGlpNum_t piv_val);

void ILLfct_load_workvector (
	lpinfo * lp,
	svector * s),
  ILLfct_zero_workvector (
	lpinfo * lp),
  ILLfct_set_variable_type (
	lpinfo * lp),
  ILLfct_compute_pobj (
	lpinfo * lp),
  ILLfct_compute_dobj (
	lpinfo * lp),
  ILLfct_compute_xbz (
	lpinfo * lp),
  ILLfct_compute_piz (
	lpinfo * lp),
  ILLfct_compute_phaseI_xbz (
	lpinfo * lp),
  ILLfct_compute_phaseI_piz (
	lpinfo * lp),
  ILLfct_compute_dz (
	lpinfo * lp),
  ILLfct_compute_phaseI_dz (
	lpinfo * lp),
  ILLfct_compute_yz (
	lpinfo * lp,
	svector * yz,
	svector * updz,
	int ecol),
  ILLfct_compute_zz (
	lpinfo * lp,
	svector * zz,
	int lindex),
  ILLfct_compute_binvrow (
	lpinfo * lp,
	svector * zz,
	int row,
	EGlpNum_t ztoler),
  ILLfct_compute_vA (
	lpinfo * lp,
	svector * v,
	EGlpNum_t * vA),
  ILLfct_compute_psteep_upv (
	lpinfo * lp,
	svector * swz),
  ILLfct_compute_dsteep_upv (
	lpinfo * lp,
	svector * swz),
  ILLfct_update_basis_info (
	lpinfo * lp,
	int eindex,
	int lindex,
	int lvstat),
  ILLfct_update_xz (
	lpinfo * lp,
	EGlpNum_t tz,
	int eindex,
	int lindex),
  ILLfct_update_piz (
	lpinfo * lp,
	EGlpNum_t alpha),
  ILLfct_update_pIpiz (
	lpinfo * lp,
	svector * z,
	const EGlpNum_t alpha),
  ILLfct_update_dz (
	lpinfo * lp,
	int eindex,
	EGlpNum_t alpha),
  ILLfct_update_pIdz (
	lpinfo * lp,
	svector * zA,
	int eindex,
	const EGlpNum_t alpha),
  ILLfct_unroll_bound_change (
	lpinfo * lp),
  ILLfct_unroll_coef_change (
	lpinfo * lp),
  ILLfct_check_pfeasible (
	lpinfo * lp,
	feas_info * fs,
	const EGlpNum_t ftol),
  ILLfct_check_pIpfeasible (
	lpinfo * lp,
	feas_info * fs,
	EGlpNum_t ftol),
  ILLfct_check_dfeasible (
	lpinfo * lp,
	feas_info * fs,
	const EGlpNum_t ftol),
  ILLfct_dual_adjust (
	lpinfo * lp,
	const EGlpNum_t ftol),
  ILLfct_dphaseI_simple_update (
	lpinfo * lp,
	EGlpNum_t ftol),
  ILLfct_set_status_values (
	lpinfo * lp,
	int pstatus,
	int dstatus,
	int ptype,
	int dtype),
  ILLfct_init_counts (
	lpinfo * lp),
  ILLfct_update_counts (
	lpinfo * lp,
	int f,
	int upi,
	const EGlpNum_t upd),
  ILLfct_print_counts (
	lpinfo * lp),
  ILLfct_check_pIdfeasible (
	lpinfo * lp,
	feas_info * fs,
	EGlpNum_t ftol),
  ILLfct_update_pfeas (
	lpinfo * lp,
	int lindex,
	svector * srhs),
  ILLfct_compute_ppIzz (
	lpinfo * lp,
	svector * srhs,
	svector * ssoln),
  ILLfct_update_ppI_prices (
	lpinfo * lp,
	price_info * pinf,
	svector * srhs,
	svector * ssoln,
	int eindex,
	int lindex,
	const EGlpNum_t alpha),
  ILLfct_update_dfeas (
	lpinfo * lp,
	int eindex,
	svector * srhs),
  ILLfct_compute_dpIy (
	lpinfo * lp,
	svector * srhs,
	svector * ssoln),
  ILLfct_update_dpI_prices (
	lpinfo * lp,
	price_info * pinf,
	svector * srhs,
	svector * ssoln,
	int lindex,
	EGlpNum_t alpha),
  ILLfct_update_dIIfeas (
	lpinfo * lp,
	int eindex,
	svector * srhs),
  ILLfct_compute_dpIIy (
	lpinfo * lp,
	svector * srhs,
	svector * ssoln),
  ILLfct_update_dpII_prices (
	lpinfo * lp,
	price_info * pinf,
	svector * srhs,
	svector * ssoln,
	/*int eindex,*/
	int lindex,
	EGlpNum_t eval,
	EGlpNum_t alpha);

void fct_test_workvector (
	lpinfo * lp),
  fct_test_pfeasible (
	lpinfo * lp),
  fct_test_dfeasible (
	lpinfo * lp),
  fct_test_pI_x (
	lpinfo * lp,
	price_info * p),
  fct_test_pII_x (
	lpinfo * lp,
	price_info * p),
  fct_test_pI_pi_dz (
	lpinfo * lp,
	price_info * p),
  fct_test_pII_pi_dz (
	lpinfo * lp,
	price_info * p);

bndinfo *ILLfct_new_bndinfo (
	void);
void ILLfct_free_bndinfo (
	bndinfo * binfo);

#endif /* __FUNCTIONS_H */
