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

/* RCSINFO $Id: dstruct.h,v 1.3 2003/11/05 16:57:39 meven Exp $ */
/****************************************************************************/
/*                                                                          */
/*                           svector.h                                      */
/*                                                                          */
/****************************************************************************/

#ifndef __SVECTOR_H
#define __SVECTOR_H

#include "qs_config.h"
#include "eg_io.h"

typedef struct svector
{
	int nzcnt;
	int *indx;
	int size;
	EGlpNum_t *coef;
}
svector;

void ILLsvector_init (
	svector * s),
  ILLsvector_free (
	svector * s);

int ILLsvector_alloc (
	svector * s,
	int nzcnt),
  ILLsvector_copy (
	const svector * s_in,
	svector * s_out);

#endif /* __SVECTOR_H */

/****************************************************************************/
/*                                                                          */
/*                           heap.h                                         */
/*                                                                          */
/****************************************************************************/

#ifndef __HEAP_H
#define __HEAP_H

typedef struct
{
	int *entry;
	int *loc;
	EGlpNum_t *key;
	int hexist;
	int maxsize;
	int size;
}
heap;

void ILLheap_insert (
	heap * const h,
	int const ix),
  ILLheap_modify (
	heap * const h,
	int const ix),
  ILLheap_delete (
	heap * const h,
	int const ix),
  ILLheap_init (
	heap * const h),
  ILLheap_free (
	heap * const h);

int ILLheap_findmin (
	heap * const h),
  ILLheap_build (
	heap * const h,
	int const nelems,
	EGlpNum_t * key);

#endif /* __HEAP_H */

/****************************************************************************/
/*                                                                          */
/*                         matrix.h                                         */
/*                                                                          */
/****************************************************************************/

#ifndef __MATRIX_H
#define __MATRIX_H

typedef struct ILLmatrix
{
	EGlpNum_t *matval;						/* The coefficients.                       */
	int *matcnt;									/* Number of coefs in each col.            */
	int *matind;									/* The row indices of the coefs.           */
	int *matbeg;									/* The start of each col.                  */
	int matcols;									/* Number of columns.                      */
	int matrows;
	int matcolsize;								/* Length of matbeg and matcnt.            */
	int matsize;									/* Length of matind and matval.            */
	int matfree;									/* Free space at end of matind.            */
	/* Note: free elements marked by -1 in     */
	/* matind; we keep at least 1 free at end. */
}
ILLmatrix;

void ILLmatrix_init (
	ILLmatrix * A);
void ILLmatrix_free (
	ILLmatrix * A);
void ILLmatrix_prt (
	EGioFile_t * fd,
	ILLmatrix * A);

#endif /* __MATRIX_H */
