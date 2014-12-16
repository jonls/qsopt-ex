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

/*  $RCSfile: ratio_EGLPNUM_TYPENAME.h,v $ $Revision: 1.3 $ $Date: 2003/11/05 16:57:39 $" */
#ifndef EGLPNUM_TYPENAME___RATIO_H
#define EGLPNUM_TYPENAME___RATIO_H
#include "basicdefs.h"
typedef struct EGLPNUM_TYPENAME_ratio_res
{
	EGLPNUM_TYPE tz;
	int eindex;
	int lindex;
	int lvstat;
	int ratio_stat;
	int boundch;
	int coeffch;
	EGLPNUM_TYPE lbound;
	EGLPNUM_TYPE ecoeff;
	EGLPNUM_TYPE pivotval;
}
EGLPNUM_TYPENAME_ratio_res;

void EGLPNUM_TYPENAME_ILLratio_pI_test (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	int const eindex,
	int const dir,
	EGLPNUM_TYPENAME_ratio_res * const rs),
  EGLPNUM_TYPENAME_ILLratio_pII_test (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	int const eindex,
	int const dir,
	EGLPNUM_TYPENAME_ratio_res * const rs),
  EGLPNUM_TYPENAME_ILLratio_dI_test (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	int const lindex,
	int const lvstat,
	EGLPNUM_TYPENAME_ratio_res * const rs),
  EGLPNUM_TYPENAME_ILLratio_dII_test (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	/*int const lindex,*/
	int const lvstat,
	EGLPNUM_TYPENAME_ratio_res * const rs),
  EGLPNUM_TYPENAME_ILLratio_longdII_test (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	int const lindex,
	int const lvstat,
	EGLPNUM_TYPENAME_ratio_res * const rs),
  EGLPNUM_TYPENAME_ILLratio_pivotin_test (
	EGLPNUM_TYPENAME_lpinfo * const lp,
	int *const rlist,
	int const rcnt,
	EGLPNUM_TYPENAME_ratio_res * const rs);

#endif /* EGLPNUM_TYPENAME___RATIO_H */
