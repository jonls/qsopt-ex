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

/*  $RCSfile: ratio.h,v $ $Revision: 1.3 $ $Date: 2003/11/05 16:57:39 $" */
#ifndef __RATIO_H
#define __RATIO_H
#include "basicdefs.h"
typedef struct ratio_res
{
	EGlpNum_t tz;
	int eindex;
	int lindex;
	int lvstat;
	int ratio_stat;
	int boundch;
	int coeffch;
	EGlpNum_t lbound;
	EGlpNum_t ecoeff;
	EGlpNum_t pivotval;
}
ratio_res;

void ILLratio_pI_test (
	lpinfo * const lp,
	int const eindex,
	int const dir,
	ratio_res * const rs),
  ILLratio_pII_test (
	lpinfo * const lp,
	int const eindex,
	int const dir,
	ratio_res * const rs),
  ILLratio_dI_test (
	lpinfo * const lp,
	int const lindex,
	int const lvstat,
	ratio_res * const rs),
  ILLratio_dII_test (
	lpinfo * const lp,
	/*int const lindex,*/
	int const lvstat,
	ratio_res * const rs),
  ILLratio_longdII_test (
	lpinfo * const lp,
	int const lindex,
	int const lvstat,
	ratio_res * const rs),
  ILLratio_pivotin_test (
	lpinfo * const lp,
	int *const rlist,
	int const rcnt,
	ratio_res * const rs);

#endif /* __RATIO_H */
