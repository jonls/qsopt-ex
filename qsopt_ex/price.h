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

/*  $RCSfile: price_EGLPNUM_TYPENAME.h,v $ $Revision: 1.3 $ $Date: 2003/11/05 16:57:39 $" */
#ifndef EGLPNUM_TYPENAME___PRICE_H
#define EGLPNUM_TYPENAME___PRICE_H

#include "dstruct_EGLPNUM_TYPENAME.h"
#include "basicdefs.h"

typedef struct EGLPNUM_TYPENAME_price_res
{
	int eindex;
	int dir;
	int lindex;
	int lvstat;
	int price_stat;
	EGLPNUM_TYPE dinfeas;
	EGLPNUM_TYPE pinfeas;
}
EGLPNUM_TYPENAME_price_res;

int EGLPNUM_TYPENAME_ILLprice_test_for_heap (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	EGLPNUM_TYPENAME_price_info * const pinf,
	int const nkeys,
	EGLPNUM_TYPE * keylist,
	int const algo,
	int const upd),
  EGLPNUM_TYPENAME_ILLprice_build_heap (
	EGLPNUM_TYPENAME_price_info * const pinf,
	int const nkeys,
	EGLPNUM_TYPE * keylist),
  EGLPNUM_TYPENAME_ILLprice_build_pricing_info (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	EGLPNUM_TYPENAME_price_info * const pinf,
	int const phase),
  EGLPNUM_TYPENAME_ILLprice_update_pricing_info (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	EGLPNUM_TYPENAME_price_info * const pinf,
	int const phase,
	EGLPNUM_TYPENAME_svector * const wz,
	int const eindex,
	int const lindex,
	EGLPNUM_TYPE y),
  EGLPNUM_TYPENAME_ILLprice_get_price (
	EGLPNUM_TYPENAME_price_info * const p,
	int const phase),
  EGLPNUM_TYPENAME_ILLprice_build_mpartial_info (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	EGLPNUM_TYPENAME_price_info * const pinf,
	int const pricetype),
  EGLPNUM_TYPENAME_ILLprice_build_pdevex_norms (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	EGLPNUM_TYPENAME_p_devex_info * const pdinfo,
	int const reinit),
  EGLPNUM_TYPENAME_ILLprice_update_pdevex_norms (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	EGLPNUM_TYPENAME_p_devex_info * const pdinfo,
	int const eindex,
	EGLPNUM_TYPE yl),
  EGLPNUM_TYPENAME_ILLprice_build_psteep_norms (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	EGLPNUM_TYPENAME_p_steep_info * const psinfo),
  EGLPNUM_TYPENAME_ILLprice_build_ddevex_norms (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	EGLPNUM_TYPENAME_d_devex_info * const ddinfo,
	int const reinit),
  EGLPNUM_TYPENAME_ILLprice_update_ddevex_norms (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	EGLPNUM_TYPENAME_d_devex_info * const ddinfo,
	int const eindex,
	EGLPNUM_TYPE yl),
  EGLPNUM_TYPENAME_ILLprice_build_dsteep_norms (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	EGLPNUM_TYPENAME_d_steep_info * const dsinfo),
  EGLPNUM_TYPENAME_ILLprice_get_dsteep_norms (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	int const count,
	int *constrowind,
	EGLPNUM_TYPE * const norms),
  EGLPNUM_TYPENAME_ILLprice_get_rownorms (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	EGLPNUM_TYPENAME_price_info * const pinf,
	EGLPNUM_TYPE * const rnorms),
  EGLPNUM_TYPENAME_ILLprice_get_colnorms (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	EGLPNUM_TYPENAME_price_info * const pinf,
	EGLPNUM_TYPE * const cnorms),
  EGLPNUM_TYPENAME_ILLprice_get_newnorms (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	int const nelems,
	EGLPNUM_TYPE * const norms,
	int *const matcnt,
	int *const matbeg,
	int *const matind,
	EGLPNUM_TYPE * const matval,
	int const option),
  EGLPNUM_TYPENAME_ILLprice_get_new_rownorms (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	int const newrows,
	EGLPNUM_TYPE * const rnorms,
	int *const rmatcnt,
	int *const rmatbeg,
	int *const rmatind,
	EGLPNUM_TYPE * const rmatval),
  EGLPNUM_TYPENAME_ILLprice_get_new_colnorms (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	int const newrows,
	EGLPNUM_TYPE * const rnorms,
	int *const matcnt,
	int *const matbeg,
	int *const matind,
	EGLPNUM_TYPE * const matval),
  EGLPNUM_TYPENAME_ILLprice_load_rownorms (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	EGLPNUM_TYPE * const rnorms,
	EGLPNUM_TYPENAME_price_info * const pinf),
  EGLPNUM_TYPENAME_ILLprice_load_colnorms (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	EGLPNUM_TYPE * const cnorms,
	EGLPNUM_TYPENAME_price_info * const pinf);


void EGLPNUM_TYPENAME_ILLprice_free_heap (
	EGLPNUM_TYPENAME_price_info * const pinf),
  EGLPNUM_TYPENAME_ILLprice_init_pricing_info (
	EGLPNUM_TYPENAME_price_info * const pinf),
  EGLPNUM_TYPENAME_ILLprice_free_pricing_info (
	EGLPNUM_TYPENAME_price_info * const pinf),
  EGLPNUM_TYPENAME_ILLprice_free_mpartial_info (
	EGLPNUM_TYPENAME_mpart_info * p),
  EGLPNUM_TYPENAME_ILLprice_init_mpartial_price (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	EGLPNUM_TYPENAME_price_info * const pinf,
	int const phase,
	int const pricetype),
  EGLPNUM_TYPENAME_ILLprice_update_mpartial_price (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	EGLPNUM_TYPENAME_price_info * const pinf,
	int const phase,
	int const pricetype),
  EGLPNUM_TYPENAME_ILLprice_delete_onempart_price (
	/*EGLPNUM_TYPENAME_lpinfo * const lp,*/
	EGLPNUM_TYPENAME_price_info * const pinf,
	int const indx,
	int const pricetype),
  EGLPNUM_TYPENAME_ILLprice_mpartial_group (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	EGLPNUM_TYPENAME_mpart_info * const p,
	int const phase,
	int const g,
	int const pricetype),
  EGLPNUM_TYPENAME_ILLprice_column (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	int const ix,
	int const phase,
	EGLPNUM_TYPENAME_price_res * const pr),
  EGLPNUM_TYPENAME_ILLprice_row (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	int const ix,
	int const phase,
	EGLPNUM_TYPENAME_price_res * const pr),
  EGLPNUM_TYPENAME_ILLprice_update_psteep_norms (
	EGLPNUM_TYPENAME_lpinfo * lp,
	EGLPNUM_TYPENAME_p_steep_info * psinfo,
	EGLPNUM_TYPENAME_svector * wz,
	int eindex,
	EGLPNUM_TYPE yl),
  EGLPNUM_TYPENAME_ILLprice_update_dsteep_norms (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	EGLPNUM_TYPENAME_d_steep_info * const dsinfo,
	EGLPNUM_TYPENAME_svector * const wz,
	int const lindex,
	EGLPNUM_TYPE yl),
  EGLPNUM_TYPENAME_ILLprice_compute_dual_inf (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	EGLPNUM_TYPENAME_price_info * const p,
	int *const ix,
	int const icnt,
	int const phase),
  EGLPNUM_TYPENAME_ILLprice_primal (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	EGLPNUM_TYPENAME_price_info * const pinf,
	EGLPNUM_TYPENAME_price_res * const pr,
	int const phase),
  EGLPNUM_TYPENAME_ILLprice_compute_primal_inf (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	EGLPNUM_TYPENAME_price_info * const p,
	int *const ix,
	int const icnt,
	int const phase),
  EGLPNUM_TYPENAME_ILLprice_dual (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	EGLPNUM_TYPENAME_price_info * const pinf,
	int const phase,
	EGLPNUM_TYPENAME_price_res * const pr);

void EGLPNUM_TYPENAME_test_dsteep_norms (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	EGLPNUM_TYPENAME_price_info * const p);

#endif /* EGLPNUM_TYPENAME___PRICE_H */
