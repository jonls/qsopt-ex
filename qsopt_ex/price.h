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

/*  $RCSfile: price.h,v $ $Revision: 1.3 $ $Date: 2003/11/05 16:57:39 $" */
#ifndef __PRICE_H
#define __PRICE_H

#include "dstruct.h"
#include "basicdefs.h"

typedef struct price_res
{
	int eindex;
	int dir;
	int lindex;
	int lvstat;
	int price_stat;
	EGlpNum_t dinfeas;
	EGlpNum_t pinfeas;
}
price_res;

int ILLprice_test_for_heap (
	lpinfo * const lp,
	price_info * const pinf,
	int const nkeys,
	EGlpNum_t * keylist,
	int const algo,
	int const upd),
  ILLprice_build_heap (
	price_info * const pinf,
	int const nkeys,
	EGlpNum_t * keylist),
  ILLprice_build_pricing_info (
	lpinfo * const lp,
	price_info * const pinf,
	int const phase),
  ILLprice_update_pricing_info (
	lpinfo * const lp,
	price_info * const pinf,
	int const phase,
	svector * const wz,
	int const eindex,
	int const lindex,
	EGlpNum_t y),
  ILLprice_get_price (
	price_info * const p,
	int const phase),
  ILLprice_build_mpartial_info (
	lpinfo * const lp,
	price_info * const pinf,
	int const pricetype),
  ILLprice_build_pdevex_norms (
	lpinfo * const lp,
	p_devex_info * const pdinfo,
	int const reinit),
  ILLprice_update_pdevex_norms (
	lpinfo * const lp,
	p_devex_info * const pdinfo,
	int const eindex,
	EGlpNum_t yl),
  ILLprice_build_psteep_norms (
	lpinfo * const lp,
	p_steep_info * const psinfo),
  ILLprice_build_ddevex_norms (
	lpinfo * const lp,
	d_devex_info * const ddinfo,
	int const reinit),
  ILLprice_update_ddevex_norms (
	lpinfo * const lp,
	d_devex_info * const ddinfo,
	int const eindex,
	EGlpNum_t yl),
  ILLprice_build_dsteep_norms (
	lpinfo * const lp,
	d_steep_info * const dsinfo),
  ILLprice_get_dsteep_norms (
	lpinfo * const lp,
	int const count,
	int *constrowind,
	EGlpNum_t * const norms),
  ILLprice_get_rownorms (
	lpinfo * const lp,
	price_info * const pinf,
	EGlpNum_t * const rnorms),
  ILLprice_get_colnorms (
	lpinfo * const lp,
	price_info * const pinf,
	EGlpNum_t * const cnorms),
  ILLprice_get_newnorms (
	lpinfo * const lp,
	int const nelems,
	EGlpNum_t * const norms,
	int *const matcnt,
	int *const matbeg,
	int *const matind,
	EGlpNum_t * const matval,
	int const option),
  ILLprice_get_new_rownorms (
	lpinfo * const lp,
	int const newrows,
	EGlpNum_t * const rnorms,
	int *const rmatcnt,
	int *const rmatbeg,
	int *const rmatind,
	EGlpNum_t * const rmatval),
  ILLprice_get_new_colnorms (
	lpinfo * const lp,
	int const newrows,
	EGlpNum_t * const rnorms,
	int *const matcnt,
	int *const matbeg,
	int *const matind,
	EGlpNum_t * const matval),
  ILLprice_load_rownorms (
	lpinfo * const lp,
	EGlpNum_t * const rnorms,
	price_info * const pinf),
  ILLprice_load_colnorms (
	lpinfo * const lp,
	EGlpNum_t * const cnorms,
	price_info * const pinf);


void ILLprice_free_heap (
	price_info * const pinf),
  ILLprice_init_pricing_info (
	price_info * const pinf),
  ILLprice_free_pricing_info (
	price_info * const pinf),
  ILLprice_free_mpartial_info (
	mpart_info * p),
  ILLprice_init_mpartial_price (
	lpinfo * const lp,
	price_info * const pinf,
	int const phase,
	int const pricetype),
  ILLprice_update_mpartial_price (
	lpinfo * const lp,
	price_info * const pinf,
	int const phase,
	int const pricetype),
  ILLprice_delete_onempart_price (
	/*lpinfo * const lp,*/
	price_info * const pinf,
	int const indx,
	int const pricetype),
  ILLprice_mpartial_group (
	lpinfo * const lp,
	mpart_info * const p,
	int const phase,
	int const g,
	int const pricetype),
  ILLprice_column (
	lpinfo * const lp,
	int const ix,
	int const phase,
	price_res * const pr),
  ILLprice_row (
	lpinfo * const lp,
	int const ix,
	int const phase,
	price_res * const pr),
  ILLprice_update_psteep_norms (
	lpinfo * lp,
	p_steep_info * psinfo,
	svector * wz,
	int eindex,
	EGlpNum_t yl),
  ILLprice_update_dsteep_norms (
	lpinfo * const lp,
	d_steep_info * const dsinfo,
	svector * const wz,
	int const lindex,
	EGlpNum_t yl),
  ILLprice_compute_dual_inf (
	lpinfo * const lp,
	price_info * const p,
	int *const ix,
	int const icnt,
	int const phase),
  ILLprice_primal (
	lpinfo * const lp,
	price_info * const pinf,
	price_res * const pr,
	int const phase),
  ILLprice_compute_primal_inf (
	lpinfo * const lp,
	price_info * const p,
	int *const ix,
	int const icnt,
	int const phase),
  ILLprice_dual (
	lpinfo * const lp,
	price_info * const pinf,
	int const phase,
	price_res * const pr);

void test_dsteep_norms (
	lpinfo * const lp,
	price_info * const p);

#endif /* __PRICE_H */
