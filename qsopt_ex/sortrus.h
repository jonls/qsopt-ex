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

#ifndef EGLPNUM_TYPENAME___SORTRUS_H__
#define EGLPNUM_TYPENAME___SORTRUS_H__
/****************************************************************************/
/*                                                                          */
/*                             sortrus.c                                    */
/*                                                                          */
/****************************************************************************/

#include "urandom.h"

void EGLPNUM_TYPENAME_ILLutil_EGlpNum_perm_quicksort (
	int *perm,
	EGLPNUM_TYPE * len,
	int n),
  EGLPNUM_TYPENAME_ILLutil_EGlpNum_rselect (
	int *arr,
	int l,
	int r,
	int m,
	EGLPNUM_TYPE * coord,
	ILLrandstate * rstate);
#endif
