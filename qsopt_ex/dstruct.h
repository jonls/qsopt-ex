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

/* RCSINFO $Id: dstruct_EGLPNUM_TYPENAME.h,v 1.3 2003/11/05 16:57:39 meven Exp $ */
/****************************************************************************/
/*                                                                          */
/*                           EGLPNUM_TYPENAME_svector.h                                      */
/*                                                                          */
/****************************************************************************/

#ifndef EGLPNUM_TYPENAME___SVECTOR_H
#define EGLPNUM_TYPENAME___SVECTOR_H

#include "eg_io.h"

typedef struct EGLPNUM_TYPENAME_svector
{
	int nzcnt;
	int *indx;
	int size;
	EGLPNUM_TYPE *coef;
}
EGLPNUM_TYPENAME_svector;

void EGLPNUM_TYPENAME_ILLsvector_init (
	EGLPNUM_TYPENAME_svector * s),
  EGLPNUM_TYPENAME_ILLsvector_free (
	EGLPNUM_TYPENAME_svector * s);

int EGLPNUM_TYPENAME_ILLsvector_alloc (
	EGLPNUM_TYPENAME_svector * s,
	int nzcnt),
  EGLPNUM_TYPENAME_ILLsvector_copy (
	const EGLPNUM_TYPENAME_svector * s_in,
	EGLPNUM_TYPENAME_svector * s_out);

#endif /* EGLPNUM_TYPENAME___SVECTOR_H */

/****************************************************************************/
/*                                                                          */
/*                           EGLPNUM_TYPENAME_heap.h                                         */
/*                                                                          */
/****************************************************************************/

#ifndef EGLPNUM_TYPENAME___HEAP_H
#define EGLPNUM_TYPENAME___HEAP_H

typedef struct
{
	int *entry;
	int *loc;
	EGLPNUM_TYPE *key;
	int hexist;
	int maxsize;
	int size;
}
EGLPNUM_TYPENAME_heap;

void EGLPNUM_TYPENAME_ILLheap_insert (
	EGLPNUM_TYPENAME_heap * const h,
	int const ix),
  EGLPNUM_TYPENAME_ILLheap_modify (
	EGLPNUM_TYPENAME_heap * const h,
	int const ix),
  EGLPNUM_TYPENAME_ILLheap_delete (
	EGLPNUM_TYPENAME_heap * const h,
	int const ix),
  EGLPNUM_TYPENAME_ILLheap_init (
	EGLPNUM_TYPENAME_heap * const h),
  EGLPNUM_TYPENAME_ILLheap_free (
	EGLPNUM_TYPENAME_heap * const h);

int EGLPNUM_TYPENAME_ILLheap_findmin (
	EGLPNUM_TYPENAME_heap * const h),
  EGLPNUM_TYPENAME_ILLheap_build (
	EGLPNUM_TYPENAME_heap * const h,
	int const nelems,
	EGLPNUM_TYPE * key);

#endif /* EGLPNUM_TYPENAME___HEAP_H */

/****************************************************************************/
/*                                                                          */
/*                         matrix.h                                         */
/*                                                                          */
/****************************************************************************/

#ifndef EGLPNUM_TYPENAME___MATRIX_H
#define EGLPNUM_TYPENAME___MATRIX_H

typedef struct EGLPNUM_TYPENAME_ILLmatrix
{
	EGLPNUM_TYPE *matval;						/* The coefficients.                       */
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
EGLPNUM_TYPENAME_ILLmatrix;

void EGLPNUM_TYPENAME_ILLmatrix_init (
	EGLPNUM_TYPENAME_ILLmatrix * A);
void EGLPNUM_TYPENAME_ILLmatrix_free (
	EGLPNUM_TYPENAME_ILLmatrix * A);
void EGLPNUM_TYPENAME_ILLmatrix_prt (
	EGioFile_t * fd,
	EGLPNUM_TYPENAME_ILLmatrix * A);

#endif /* EGLPNUM_TYPENAME___MATRIX_H */
