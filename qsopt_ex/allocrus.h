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

#ifndef __ALLOCRUS_H__
#define __ALLOCRUS_H__
/****************************************************************************/
/*                                                                          */
/*                             allocrus.c                                   */
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
/*  Date: February 24, 1995 (cofeb24)                                       */
/*  see allocrus.c                                                          */
/*                                                                          */
/****************************************************************************/

#include <stdlib.h>
#include <stdio.h>

#include "logging-private.h"

extern int ILLTRACE_MALLOC;

#define ILL_UTIL_SAFE_MALLOC(nnum,type,varname)                             \
    (((ILLTRACE_MALLOC) ? QSlog("%s.%d: %s: ILL_UTIL_SAFE_MALLOC: %s = %d * %s\n", __FILE__, __LINE__, __DEV_FUNCTION__, #varname, nnum, #type) : 0), \
     (type *) ILLutil_allocrus (((size_t) (nnum)) * sizeof (type)))

#define ILL_IFFREE(object,type) {                                           \
    if ((object)) {                                                        \
       ILLutil_freerus ((void *) (object));                                 \
       object = (type *) NULL;                                             \
    }}

#define ILL_PTRWORLD_ALLOC_ROUTINE(type, ptr_alloc_r, ptr_bulkalloc_r)        \
                                                                             \
static int ptr_bulkalloc_r (ILLptrworld *world, int nalloc)                   \
{                                                                            \
    ILLbigchunkptr *bp;                                                       \
    int i;                                                                   \
    int count = ILL_BIGCHUNK / sizeof ( type );                               \
    type *p;                                                                 \
                                                                             \
    while (nalloc > 0) {                                                     \
        bp = ILLutil_bigchunkalloc ();                                        \
        if (bp == (ILLbigchunkptr *) NULL) {                                  \
            QSlog("ptr alloc failed\n");                          \
            return 1;                                                        \
        }                                                                    \
        bp->next = world->chunklist ;                                        \
        world->chunklist = bp;                                               \
                                                                             \
        p = ( type * ) bp->this_one;                                         \
        for (i=count-2; i>=0; i--) {                                         \
            p[i].next = &p[i+1];                                             \
        }                                                                    \
        p[count - 1].next = (type *) world->freelist;                        \
        world->freelist = (void *) p;                                        \
        nalloc -= count;                                                     \
    }                                                                        \
    return 0;                                                                \
}                                                                            \
                                                                             \
static type *ptr_alloc_r (ILLptrworld *world)                                 \
{                                                                            \
    type *p;                                                                 \
                                                                             \
    if (world->freelist == (void *) NULL) {                                  \
        if (ptr_bulkalloc_r (world, 1)) {                                    \
            QSlog("ptr alloc failed\n");                          \
            return ( type * ) NULL;                                          \
        }                                                                    \
    }                                                                        \
    p = (type *) world->freelist ;                                           \
    world->freelist = (void *) p->next;                                      \
                                                                             \
    return p;                                                                \
}

#define ILL_PTRWORLD_FREE_ROUTINE(type, ptr_free_r)                           \
                                                                             \
static void ptr_free_r (ILLptrworld *world, type *p)                          \
{                                                                            \
    p->next = (type *) world->freelist ;                                     \
    world->freelist = (void *) p;                                            \
}

#define ILL_PTRWORLD_LISTADD_ROUTINE(type, entrytype, ptr_listadd_r, ptr_alloc_r) \
                                                                             \
static int ptr_listadd_r (type **list, entrytype x, ILLptrworld *world)       \
{                                                                            \
    if (list != (type **) NULL) {                                            \
        type *p = ptr_alloc_r (world);                                       \
                                                                             \
        if (p == (type *) NULL) {                                            \
            QSlog("ptr list add failed\n");                       \
            return 1;                                                        \
        }                                                                    \
        p->this = x;                                                         \
        p->next = *list;                                                     \
        *list = p;                                                           \
    }                                                                        \
    return 0;                                                                \
}

#define ILL_PTRWORLD_LISTFREE_ROUTINE(type, ptr_listfree_r, ptr_free_r)       \
                                                                             \
static void ptr_listfree_r (ILLptrworld *world, type *p)                      \
{                                                                            \
    type *next;                                                              \
                                                                             \
    while (p != (type *) NULL) {                                             \
        next = p->next;                                                      \
        ptr_free_r (world, p);                                               \
        p = next;                                                            \
    }                                                                        \
}

#define ILL_PTRWORLD_LEAKS_ROUTINE(type, ptr_leaks_r, field, fieldtype)       \
                                                                             \
static int ptr_leaks_r (ILLptrworld *world, int *total, int *onlist)          \
{                                                                            \
    int count = ILL_BIGCHUNK / sizeof ( type );                               \
    int duplicates = 0;                                                      \
    type * p;                                                                \
    ILLbigchunkptr *bp;                                                       \
                                                                             \
    *total = 0;                                                              \
    *onlist = 0;                                                             \
                                                                             \
    for (bp = world->chunklist ; bp; bp = bp->next)                          \
        (*total) += count;                                                   \
                                                                             \
    for (p = (type *) world->freelist ; p; p = p->next) {                    \
        (*onlist)++;                                                         \
        p-> field = ( fieldtype ) 0;                                         \
    }                                                                        \
    for (p = (type *) world->freelist ; p; p = p->next) {                    \
        if ((unsigned long) p-> field == (unsigned long) (size_t) 1)                           \
            duplicates++;                                                    \
        else                                                                 \
            p-> field = ( fieldtype ) (size_t) 1;                            \
    }                                                                        \
    if (duplicates) {                                                        \
        QSlog("WARNING: %d duplicates on ptr free list \n",       \
							duplicates);																							\
    }                                                                        \
    return *total - *onlist;                                                 \
}

#define ILL_PTRWORLD_ROUTINES(type, ptr_alloc_r, ptr_bulkalloc_r, ptr_free_r) \
ILL_PTRWORLD_ALLOC_ROUTINE (type, ptr_alloc_r, ptr_bulkalloc_r)               \
ILL_PTRWORLD_FREE_ROUTINE (type, ptr_free_r)

#define ILL_PTRWORLD_LIST_ROUTINES(type, entrytype, ptr_alloc_r, ptr_bulkalloc_r, ptr_free_r, ptr_listadd_r, ptr_listfree_r) \
ILL_PTRWORLD_ROUTINES (type, ptr_alloc_r, ptr_bulkalloc_r, ptr_free_r)        \
ILL_PTRWORLD_LISTADD_ROUTINE (type, entrytype, ptr_listadd_r, ptr_alloc_r)    \
ILL_PTRWORLD_LISTFREE_ROUTINE (type, ptr_listfree_r, ptr_free_r)

#define ILL_BIGCHUNK ((int) ((1<<16) - sizeof (ILLbigchunkptr) - 16))

struct ILLbigchunk;

typedef struct ILLbigchunkptr
{
	void *this_one;
	struct ILLbigchunk *this_chunk;
	struct ILLbigchunkptr *next;
}
ILLbigchunkptr;


typedef struct ILLptrworld
{
	int refcount;
	void *freelist;
	ILLbigchunkptr *chunklist;
}
ILLptrworld;



void *ILLutil_allocrus (
	size_t size),
 *ILLutil_reallocrus (
	void *ptr,
	size_t size),
  ILLutil_freerus (
	void *p),
  ILLutil_bigchunkfree (
	ILLbigchunkptr * bp),
  ILLptrworld_init (
	ILLptrworld * world),
  ILLptrworld_add (
	ILLptrworld * world),
  ILLptrworld_delete (
	ILLptrworld * world);

int ILLutil_reallocrus_scale (
	void **pptr,
	int *pnnum,
	int count,
	double scale,
	size_t size),
  ILLutil_reallocrus_count (
	void **pptr,
	int count,
	size_t size);

ILLbigchunkptr *ILLutil_bigchunkalloc (
	void);



#endif
