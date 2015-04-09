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

#ifndef EGLPNUM_TYPENAME___DHEAPS_I_H__
#define EGLPNUM_TYPENAME___DHEAPS_I_H__

#include "eg_lpnum.h"

/****************************************************************************/
/*                                                                          */
/*                             dheaps_i.c                                   */
/*                                                                          */
/****************************************************************************/

typedef struct EGLPNUM_TYPENAME_ILLdheap
{
	EGLPNUM_TYPE *key;
	int *entry;
	int *loc;
	int total_space;
	int size;
}
EGLPNUM_TYPENAME_ILLdheap;

void EGLPNUM_TYPENAME_ILLutil_dheap_free (
	EGLPNUM_TYPENAME_ILLdheap * h),
  EGLPNUM_TYPENAME_ILLutil_dheap_delete (
	EGLPNUM_TYPENAME_ILLdheap * h,
	int i),
  EGLPNUM_TYPENAME_ILLutil_dheap_changekey (
	EGLPNUM_TYPENAME_ILLdheap * h,
	int i,
	EGLPNUM_TYPE * newkey),
  EGLPNUM_TYPENAME_ILLutil_dheap_findmin (
	EGLPNUM_TYPENAME_ILLdheap * h,
	int *i),
  EGLPNUM_TYPENAME_ILLutil_dheap_deletemin (
	EGLPNUM_TYPENAME_ILLdheap * h,
	int *i);

int EGLPNUM_TYPENAME_ILLutil_dheap_init (
	EGLPNUM_TYPENAME_ILLdheap * h,
	int k),
  EGLPNUM_TYPENAME_ILLutil_dheap_resize (
	EGLPNUM_TYPENAME_ILLdheap * h,
	int newsize),
  EGLPNUM_TYPENAME_ILLutil_dheap_insert (
	EGLPNUM_TYPENAME_ILLdheap * h,
	int i);



#endif
