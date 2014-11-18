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

/* RCSINFO $Id: priority.c,v 1.2 2003/11/05 16:47:22 meven Exp $ */
/****************************************************************************/
/*                                                                          */
/*  This file is part of CONCORDE                                           */
/*                                                                          */
/*  (c) Copyright 1995--1999 by David Applegate, Robert Bixby,              */
/*  Vasek Chvatal, and William Cook                                         */
/*                                                                          */
/*  Permission is granted for academic research use.  For other uses,       */
/*  contact the authors for licensing options.                              */
/*                                                                          */
/*  Use at your own risk.  We make no guarantees about the                  */
/*  correctness or usefulness of this code.                                 */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/*                      PRIORITY QUEUE ROUTINES                             */
/*                                                                          */
/*                                                                          */
/*                              TSP CODE                                    */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: March 3, 1997                                                     */
/*        March 13, 2002 - Cook Modified for QS)                            */
/*  Reference: R.E. Tarjan, Data Structures and Network Algorithms          */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int ILLutil_priority_init (ILLpriority *pri, int k)                     */
/*    -h should point to a ILLpriority struct.                              */
/*    -k an initial allocation for the priority queue.                      */
/*                                                                          */
/*  void ILLutil_priority_free (ILLpriority *pri)                           */
/*    -frees the spaces allocated for the ILLpriority queue.                */
/*                                                                          */
/*  void ILLutil_priority_findmin (ILLpriority *pri, double *keyval         */
/*      void **en)                                                          */
/*    -en the entry with least key value (NULL if no entries in heap).      */
/*    -if (keyval != NULL), *keyval will be the minimum key value.          */
/*                                                                          */
/*  int ILLutil_priority_insert (ILLpriority *pri, void *data,              */
/*      double keyval, int *handle)                                         */
/*    -adds (data, keyval) to h.                                            */
/*    -handle returns a handle (>= 0) to use when deleting or changing the  */
/*     entry                                                                */
/*                                                                          */
/*  void ILLutil_priority_delete (ILLpriority *pri, int handle)             */
/*    -deletes an entry from the queue.  handle is the value returned by    */
/*     ILLutil_priority_insert.                                             */
/*                                                                          */
/*  void ILLutil_priority_deletemin (ILLpriority *pri, double *keyval,      */
/*       void **en)                                                         */
/*    -like ILLutil_priority_findmin, but also deletes the entry.           */
/*                                                                          */
/*  void ILLutil_priority_changekey (ILLpriority *pri, int handle,          */
/*      double newkey)                                                      */
/*    -changes the key of an entry in the queue.  handle is the value       */
/*     returned by ILLutil_priority_insert.                                 */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/*  NOTES:                                                                  */
/*      These priority queue routines use the ILLdheap routines to maintain */
/*  the priority queue.                                                     */
/*                                                                          */
/****************************************************************************/

#include "qs_config.h"

#include "eg_lpnum.h"

#include "machdefs.h"
#include "priority.h"
#include "allocrus.h"
#include "except.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif


int ILLutil_priority_init (
	ILLpriority * pri,
	int k)
{
	int i;
	int list;
	int rval = 0;

	pri->space = k;
	ILL_SAFE_MALLOC (pri->pri_info, k, union ILLpri_data);

	rval = ILLutil_dheap_init (&pri->heap, k);
	ILL_CLEANUP_IF (rval);

	list = -1;
	for (i = k - 1; i >= 0; i--)
	{
		pri->pri_info[i].next = list;
		list = i;
	}
	pri->freelist = list;

CLEANUP:

	if (rval)
	{
		ILL_IFFREE (pri->pri_info, union ILLpri_data);
	}
	return rval;
}

void ILLutil_priority_free (
	ILLpriority * pri)
{
	ILLutil_dheap_free (&pri->heap);
	ILL_IFFREE (pri->pri_info, union ILLpri_data);

	pri->space = 0;
}

void ILLutil_priority_findmin (
	ILLpriority * pri,
	EGlpNum_t * keyval,
	void **en)
{
	int handle;

	ILLutil_dheap_findmin (&pri->heap, &handle);

	if (handle < 0)
	{
		*en = (void *) NULL;
	}
	else
	{
		if (keyval)
			EGlpNumCopy (*keyval, pri->heap.key[handle]);
		*en = pri->pri_info[handle].data;
	}
}

int ILLutil_priority_insert (
	ILLpriority * pri,
	void *data,
	EGlpNum_t * keyval,
	int *handle)
{
	int newsize;
	int i;
	int list;
	int rval = 0;

	if (pri->freelist == -1)
	{
		/* Change from 1.3 * pri->space to avoid a warning */
		newsize = pri->space + (pri->space / 3);
		if (newsize < pri->space + 1000)
			newsize = pri->space + 1000;
		rval = ILLutil_dheap_resize (&pri->heap, newsize);
		ILL_CLEANUP_IF (rval);

		pri->pri_info =
			EGrealloc (pri->pri_info, sizeof (union ILLpri_data) * newsize);
		//rval = ILLutil_reallocrus_count ((void **) &pri->pri_info, newsize,
		//                                 sizeof (union ILLpri_data));
		//ILL_CLEANUP_IF (rval);

		list = -1;
		for (i = newsize - 1; i >= pri->space; i--)
		{
			pri->pri_info[i].next = list;
			list = i;
		}
		pri->space = newsize;
		pri->freelist = list;
	}

	i = pri->freelist;
	pri->freelist = pri->pri_info[i].next;
	pri->pri_info[i].data = data;
	EGlpNumCopy (pri->heap.key[i], *keyval);
	rval = ILLutil_dheap_insert (&pri->heap, i);
	ILL_CLEANUP_IF (rval);

	if (handle)
		*handle = i;

CLEANUP:

	return rval;
}

void ILLutil_priority_delete (
	ILLpriority * pri,
	int handle)
{
	ILLutil_dheap_delete (&pri->heap, handle);
	pri->pri_info[handle].next = pri->freelist;
	pri->freelist = handle;
}

void ILLutil_priority_deletemin (
	ILLpriority * pri,
	EGlpNum_t * keyval,
	void **en)
{
	int handle;
	void *data;

	ILLutil_dheap_deletemin (&pri->heap, &handle);

	if (handle < 0)
	{
		*en = (void *) NULL;
	}
	else
	{
		if (keyval)
			EGlpNumCopy (*keyval, pri->heap.key[handle]);
		data = pri->pri_info[handle].data;
		pri->pri_info[handle].next = pri->freelist;
		pri->freelist = handle;
		*en = data;
	}
}

void ILLutil_priority_changekey (
	ILLpriority * pri,
	int handle,
	EGlpNum_t * newkey)
{
	ILLutil_dheap_changekey (&pri->heap, handle, newkey);
}
