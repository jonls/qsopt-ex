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

#ifndef __DHEAPS_I_H__
#define __DHEAPS_I_H__
/****************************************************************************/
/*                                                                          */
/*                             dheaps_i.c                                   */
/*                                                                          */
/****************************************************************************/

typedef struct ILLdheap
{
	EGlpNum_t *key;
	int *entry;
	int *loc;
	int total_space;
	int size;
}
ILLdheap;

void ILLutil_dheap_free (
	ILLdheap * h),
  ILLutil_dheap_delete (
	ILLdheap * h,
	int i),
  ILLutil_dheap_changekey (
	ILLdheap * h,
	int i,
	EGlpNum_t * newkey),
  ILLutil_dheap_findmin (
	ILLdheap * h,
	int *i),
  ILLutil_dheap_deletemin (
	ILLdheap * h,
	int *i);

int ILLutil_dheap_init (
	ILLdheap * h,
	int k),
  ILLutil_dheap_resize (
	ILLdheap * h,
	int newsize),
  ILLutil_dheap_insert (
	ILLdheap * h,
	int i);



#endif
