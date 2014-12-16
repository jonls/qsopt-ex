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

#ifndef EGLPNUM_TYPENAME___PRIORITY_H__
#define EGLPNUM_TYPENAME___PRIORITY_H__
#include "dheaps_i_EGLPNUM_TYPENAME.h"
/****************************************************************************/
/*                                                                          */
/*                             priority.c                                   */
/*                                                                          */
/****************************************************************************/

typedef struct EGLPNUM_TYPENAME_ILLpriority
{
	EGLPNUM_TYPENAME_ILLdheap EGLPNUM_TYPENAME_heap;
	union EGLPNUM_TYPENAME_ILLpri_data
	{
		void *data;
		int next;
	}
	 *pri_info;
	int space;
	int freelist;
}
EGLPNUM_TYPENAME_ILLpriority;

void EGLPNUM_TYPENAME_ILLutil_priority_free (
	EGLPNUM_TYPENAME_ILLpriority * pri),
  EGLPNUM_TYPENAME_ILLutil_priority_delete (
	EGLPNUM_TYPENAME_ILLpriority * pri,
	int handle),
  EGLPNUM_TYPENAME_ILLutil_priority_changekey (
	EGLPNUM_TYPENAME_ILLpriority * pri,
	int handle,
	EGLPNUM_TYPE * newkey),
  EGLPNUM_TYPENAME_ILLutil_priority_findmin (
	EGLPNUM_TYPENAME_ILLpriority * pri,
	EGLPNUM_TYPE * keyval,
	void **en),
  EGLPNUM_TYPENAME_ILLutil_priority_deletemin (
	EGLPNUM_TYPENAME_ILLpriority * pri,
	EGLPNUM_TYPE * keyval,
	void **en);

int EGLPNUM_TYPENAME_ILLutil_priority_init (
	EGLPNUM_TYPENAME_ILLpriority * pri,
	int k),
  EGLPNUM_TYPENAME_ILLutil_priority_insert (
	EGLPNUM_TYPENAME_ILLpriority * pri,
	void *data,
	EGLPNUM_TYPE * keyval,
	int *handle);



#endif
