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

/* RCSINFO $Id: fct_EGLPNUM_TYPENAME.h,v 1.3 2003/11/05 16:57:39 meven Exp $ */
#ifndef EGLPNUM_TYPENAME___FUNCTIONS_H
#define EGLPNUM_TYPENAME___FUNCTIONS_H

#include "basicdefs.h"
#include "lpdefs_EGLPNUM_TYPENAME.h"


int EGLPNUM_TYPENAME_ILLfct_compute_zA (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_svector * z,
	EGLPNUM_TYPENAME_svector * zA),
  EGLPNUM_TYPENAME_ILLfct_compute_wz (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPE * wz),
  EGLPNUM_TYPENAME_ILLfct_adjust_viol_bounds (
	EGLPNUM_TYPENAME_lpinfo * lp),
  EGLPNUM_TYPENAME_ILLfct_perturb_bounds (
	EGLPNUM_TYPENAME_lpinfo * lp),
  EGLPNUM_TYPENAME_ILLfct_perturb_phaseI_bounds (
	EGLPNUM_TYPENAME_lpinfo * lp),
  EGLPNUM_TYPENAME_ILLfct_bound_shift (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int col,
	int bndtype,
	EGLPNUM_TYPE newbnd),
  EGLPNUM_TYPENAME_ILLfct_adjust_viol_coefs (
	EGLPNUM_TYPENAME_lpinfo * lp),
  EGLPNUM_TYPENAME_ILLfct_perturb_coefs (
	EGLPNUM_TYPENAME_lpinfo * lp),
  EGLPNUM_TYPENAME_ILLfct_coef_shift (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int col,
	EGLPNUM_TYPE newcoef),
  EGLPNUM_TYPENAME_ILLfct_test_pivot (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int indx,
	int indxtype,
	EGLPNUM_TYPE piv_val);

void EGLPNUM_TYPENAME_ILLfct_load_workvector (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_svector * s),
  EGLPNUM_TYPENAME_ILLfct_zero_workvector (
	EGLPNUM_TYPENAME_lpinfo * lp),
  EGLPNUM_TYPENAME_ILLfct_set_variable_type (
	EGLPNUM_TYPENAME_lpinfo * lp),
  EGLPNUM_TYPENAME_ILLfct_compute_pobj (
	EGLPNUM_TYPENAME_lpinfo * lp),
  EGLPNUM_TYPENAME_ILLfct_compute_dobj (
	EGLPNUM_TYPENAME_lpinfo * lp),
  EGLPNUM_TYPENAME_ILLfct_compute_xbz (
	EGLPNUM_TYPENAME_lpinfo * lp),
  EGLPNUM_TYPENAME_ILLfct_compute_piz (
	EGLPNUM_TYPENAME_lpinfo * lp),
  EGLPNUM_TYPENAME_ILLfct_compute_phaseI_xbz (
	EGLPNUM_TYPENAME_lpinfo * lp),
  EGLPNUM_TYPENAME_ILLfct_compute_phaseI_piz (
	EGLPNUM_TYPENAME_lpinfo * lp),
  EGLPNUM_TYPENAME_ILLfct_compute_dz (
	EGLPNUM_TYPENAME_lpinfo * lp),
  EGLPNUM_TYPENAME_ILLfct_compute_phaseI_dz (
	EGLPNUM_TYPENAME_lpinfo * lp),
  EGLPNUM_TYPENAME_ILLfct_compute_yz (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_svector * yz,
	EGLPNUM_TYPENAME_svector * updz,
	int ecol),
  EGLPNUM_TYPENAME_ILLfct_compute_zz (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_svector * zz,
	int lindex),
  EGLPNUM_TYPENAME_ILLfct_compute_binvrow (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_svector * zz,
	int row,
	EGLPNUM_TYPE ztoler),
  EGLPNUM_TYPENAME_ILLfct_compute_vA (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_svector * v,
	EGLPNUM_TYPE * vA),
  EGLPNUM_TYPENAME_ILLfct_compute_psteep_upv (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_svector * swz),
  EGLPNUM_TYPENAME_ILLfct_compute_dsteep_upv (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_svector * swz),
  EGLPNUM_TYPENAME_ILLfct_update_basis_info (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int eindex,
	int lindex,
	int lvstat),
  EGLPNUM_TYPENAME_ILLfct_update_xz (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPE tz,
	int eindex,
	int lindex),
  EGLPNUM_TYPENAME_ILLfct_update_piz (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPE alpha),
  EGLPNUM_TYPENAME_ILLfct_update_pIpiz (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_svector * z,
	const EGLPNUM_TYPE alpha),
  EGLPNUM_TYPENAME_ILLfct_update_dz (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int eindex,
	EGLPNUM_TYPE alpha),
  EGLPNUM_TYPENAME_ILLfct_update_pIdz (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_svector * zA,
	int eindex,
	const EGLPNUM_TYPE alpha),
  EGLPNUM_TYPENAME_ILLfct_unroll_bound_change (
	EGLPNUM_TYPENAME_lpinfo * lp),
  EGLPNUM_TYPENAME_ILLfct_unroll_coef_change (
	EGLPNUM_TYPENAME_lpinfo * lp),
  EGLPNUM_TYPENAME_ILLfct_check_pfeasible (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_feas_info * fs,
	const EGLPNUM_TYPE ftol),
  EGLPNUM_TYPENAME_ILLfct_check_pIpfeasible (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_feas_info * fs,
	EGLPNUM_TYPE ftol),
  EGLPNUM_TYPENAME_ILLfct_check_dfeasible (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_feas_info * fs,
	const EGLPNUM_TYPE ftol),
  EGLPNUM_TYPENAME_ILLfct_dual_adjust (
	EGLPNUM_TYPENAME_lpinfo * lp,
	const EGLPNUM_TYPE ftol),
  EGLPNUM_TYPENAME_ILLfct_dphaseI_simple_update (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPE ftol),
  EGLPNUM_TYPENAME_ILLfct_set_status_values (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int pstatus,
	int dstatus,
	int ptype,
	int dtype),
  EGLPNUM_TYPENAME_ILLfct_init_counts (
	EGLPNUM_TYPENAME_lpinfo * lp),
  EGLPNUM_TYPENAME_ILLfct_update_counts (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int f,
	int upi,
	const EGLPNUM_TYPE upd),
  EGLPNUM_TYPENAME_ILLfct_print_counts (
	EGLPNUM_TYPENAME_lpinfo * lp),
  EGLPNUM_TYPENAME_ILLfct_check_pIdfeasible (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_feas_info * fs,
	EGLPNUM_TYPE ftol),
  EGLPNUM_TYPENAME_ILLfct_update_pfeas (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int lindex,
	EGLPNUM_TYPENAME_svector * srhs),
  EGLPNUM_TYPENAME_ILLfct_compute_ppIzz (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_svector * srhs,
	EGLPNUM_TYPENAME_svector * ssoln),
  EGLPNUM_TYPENAME_ILLfct_update_ppI_prices (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_price_info * pinf,
	EGLPNUM_TYPENAME_svector * srhs,
	EGLPNUM_TYPENAME_svector * ssoln,
	int eindex,
	int lindex,
	const EGLPNUM_TYPE alpha),
  EGLPNUM_TYPENAME_ILLfct_update_dfeas (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int eindex,
	EGLPNUM_TYPENAME_svector * srhs),
  EGLPNUM_TYPENAME_ILLfct_compute_dpIy (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_svector * srhs,
	EGLPNUM_TYPENAME_svector * ssoln),
  EGLPNUM_TYPENAME_ILLfct_update_dpI_prices (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_price_info * pinf,
	EGLPNUM_TYPENAME_svector * srhs,
	EGLPNUM_TYPENAME_svector * ssoln,
	int lindex,
	EGLPNUM_TYPE alpha),
  EGLPNUM_TYPENAME_ILLfct_update_dIIfeas (
	EGLPNUM_TYPENAME_lpinfo * lp,
	int eindex,
	EGLPNUM_TYPENAME_svector * srhs),
  EGLPNUM_TYPENAME_ILLfct_compute_dpIIy (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_svector * srhs,
	EGLPNUM_TYPENAME_svector * ssoln),
  EGLPNUM_TYPENAME_ILLfct_update_dpII_prices (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_price_info * pinf,
	EGLPNUM_TYPENAME_svector * srhs,
	EGLPNUM_TYPENAME_svector * ssoln,
	/*int eindex,*/
	int lindex,
	EGLPNUM_TYPE eval,
	EGLPNUM_TYPE alpha);

void EGLPNUM_TYPENAME_fct_test_workvector (
	EGLPNUM_TYPENAME_lpinfo * lp),
  EGLPNUM_TYPENAME_fct_test_pfeasible (
	EGLPNUM_TYPENAME_lpinfo * lp),
  EGLPNUM_TYPENAME_fct_test_dfeasible (
	EGLPNUM_TYPENAME_lpinfo * lp),
  EGLPNUM_TYPENAME_fct_test_pI_x (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_price_info * p),
  EGLPNUM_TYPENAME_fct_test_pII_x (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_price_info * p),
  EGLPNUM_TYPENAME_fct_test_pI_pi_dz (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_price_info * p),
  EGLPNUM_TYPENAME_fct_test_pII_pi_dz (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_price_info * p);

EGLPNUM_TYPENAME_bndinfo *EGLPNUM_TYPENAME_ILLfct_new_bndinfo (
	void);
void EGLPNUM_TYPENAME_ILLfct_free_bndinfo (
	EGLPNUM_TYPENAME_bndinfo * binfo);

#endif /* EGLPNUM_TYPENAME___FUNCTIONS_H */
