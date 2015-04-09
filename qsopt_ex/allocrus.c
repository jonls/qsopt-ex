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

/* RCSINFO $Id: allocrus.c,v 1.2 2003/11/05 16:47:22 meven Exp $ */
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
/*                   MEMORY ALLOCATION MACROS                               */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: February 2, 1995 (cofeb16)                                        */
/*                                                                          */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  void *ILLutil_allocrus (size_t size)                                    */
/*    RETURNS a pointer to an allocated block of "size" memory.             */
/*                                                                          */
/*  void ILLutil_freerus (void *ptr)                                        */
/*    FREES ptr.                                                            */
/*                                                                          */
/*  void *ILLutil_reallocrus (void *ptr, size_t size)                       */
/*    REALLOCS ptr to size bytes.                                           */
/*                                                                          */
/*  int ILLutil_reallocrus_scale (void **pptr, int *pnnum, int count,       */
/*      double scale, size_t size)                                          */
/*    void **pptr (a reference to the pointer to the allocated space)       */
/*    int *pnnum (a reference to the number of objects in the               */
/*                allocated space)                                          */
/*    int count (a minimum value for the new nnum)                          */
/*    double scale (a scale factor to apply to nnum)                        */
/*    int size (the size of objects to be realloced)                        */
/*    RETURNS 0 if *pptr was successfully changed to point to at            */
/*            least max(*pnnum*scale, *pnnum+1000, count) objects.          */
/*            *pnnum is changed to the new object count.                    */
/*            Otherwise, prints an error message, leaves *pptr and          */
/*            *pnnum alone, and returns nonzero.                            */
/*                                                                          */
/*  int ILLutil_reallocrus_count (void **pptr, int count,                   */
/*      size_t size)                                                        */
/*    void **pptr (a reference to the pointer to the allocated space)       */
/*    int count (number of objects to be realloced)                         */
/*    int size (the size of the objects to be realloced)                    */
/*    RETURNS 0 is successful, and 1 if the realloc failed.                 */
/*                                                                          */
/*  ILLbigchunkptr *ILLutil_bigchunkalloc (void)                            */
/*         RETURNS a ILLbigchunkptr with the "this_one" field loaded with a */
/*                 a pointer to a bigchunk of memory.                       */
/*    NOTES:                                                                */
/*       The idea is to use bigchunks (the size of a bigchunk is defined    */
/*       by ILL_BIGCHUNK in util.h) to supply local routines with memory    */
/*       for ptrs, so the memory can be shared with other                   */
/*       local routines.                                                    */
/*                                                                          */
/*  ILLutil_bigchunkfree (ILLbigchunkptr *bp)                               */
/*    ACTION: Frees a ILLbigchunkptr.                                       */
/*                                                                          */
/*  void ILLptrworld_init (ILLptrworld *world)                              */
/*     initialize a ILLptrworld with 1 reference                            */
/*                                                                          */
/*  void ILLptrworld_add (ILLptrworld *world)                               */
/*     add a reference to a ILLptrworld                                     */
/*                                                                          */
/*  void ILLptrworld_delete (ILLptrworld *world)                            */
/*     delete a reference to a ptrworld, and free if no more references     */
/*                                                                          */
/****************************************************************************/

#include <stdlib.h>
#include <stdio.h>

#include "logging-private.h"

#include "except.h"
#include "util.h"

int ILLTRACE_MALLOC = 0;

typedef struct ILLbigchunk
{
	char space[ILL_BIGCHUNK];
	ILLbigchunkptr ptr;
}
ILLbigchunk;

void *ILLutil_allocrus (
	size_t size)
{
	void *mem = (void *) NULL;

	if (size == 0)
	{
		//QSlog("Warning: 0 bytes allocated");
	}

	mem = (void *) malloc (size);
	if (mem == (void *) NULL)
	{
		QSlog("Out of memory. Asked for %d bytes", (int) size);
	}
	return mem;
}

void ILLutil_freerus (
	void *p)
{
	if (!p)
	{
		//QSlog("Warning: null pointer freed");
		return;
	}

	free (p);
}

void *ILLutil_reallocrus (
	void *ptr,
	size_t size)
{
	void *newptr;

	if (!ptr)
	{
		return ILLutil_allocrus (size);
	}
	else
	{
		newptr = (void *) realloc (ptr, size);
		if (!newptr)
		{
			QSlog("Out of memory.  Tried to grow to %d bytes", (int) size);
		}
		return newptr;
	}
}

int ILLutil_reallocrus_scale (
	void **pptr,
	int *pnnum,
	int count,
	double scale,
	size_t size)
{
	int rval = 0;
	int newsize = (int) (((double) *pnnum) * scale);
	void *p;

	if (newsize < *pnnum + 1000)
		newsize = *pnnum + 1000;
	if (newsize < count)
		newsize = count;
	p = ILLutil_reallocrus (*pptr, newsize * size);
	if (!p)
	{
		rval = ILL_GENERAL_ERROR;
		ILL_REPRT ("ILLutil_reallocrus_scale failed\n");
		ILL_CLEANUP;
	}
	else
	{
		*pptr = p;
		*pnnum = newsize;
	}
CLEANUP:
	return rval;
}

int ILLutil_reallocrus_count (
	void **pptr,
	int count,
	size_t size)
{
	int rval = 0;
	void *p = ILLutil_reallocrus (*pptr, count * size);

	if (!p)
	{
		rval = ILL_GENERAL_ERROR;
		ILL_REPRT ("ILLutil_reallocrus_count failed\n");
		ILL_CLEANUP;
	}
	else
	{
		*pptr = p;
	}
CLEANUP:
	return rval;
}


ILLbigchunkptr *ILLutil_bigchunkalloc (
	void)
{
	ILLbigchunk *p;

	ILL_NEW_no_rval (p, ILLbigchunk);

	p->ptr.this_chunk = p;
	p->ptr.this_one = (void *) p->space;
CLEANUP:
	if (p == (ILLbigchunk *) NULL)
	{
		return (ILLbigchunkptr *) NULL;
	}
	return &(p->ptr);
}

void ILLutil_bigchunkfree (
	ILLbigchunkptr * bp)
{
	/* This copy is necessary since ILL_FREE zeros its first argument */
	ILLbigchunk *p = bp->this_chunk;

	ILL_IFFREE (p, ILLbigchunk);
}

void ILLptrworld_init (
	ILLptrworld * world)
{
	world->refcount = 1;
	world->freelist = (void *) NULL;
	world->chunklist = (ILLbigchunkptr *) NULL;
}

void ILLptrworld_add (
	ILLptrworld * world)
{
	world->refcount++;
}

void ILLptrworld_delete (
	ILLptrworld * world)
{
	world->refcount--;
	if (world->refcount <= 0)
	{
		ILLbigchunkptr *bp, *bpnext;

		for (bp = world->chunklist; bp; bp = bpnext)
		{
			bpnext = bp->next;
			ILLutil_bigchunkfree (bp);
		}
		world->chunklist = (ILLbigchunkptr *) NULL;
		world->freelist = (void *) NULL;
		world->refcount = 0;
	}
}
