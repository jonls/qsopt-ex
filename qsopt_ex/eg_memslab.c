/* EGlib "Efficient General Library" provides some basic structures and
 * algorithms commons in many optimization algorithms.
 *
 * Copyright (C) 2008 Daniel Espinoza.
 * 
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation; either version 2.1 of the License, or (at your
 * option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public 
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA 
 * */
/* ========================================================================= */
/** @file
 * @ingroup EGmemSlab */
/** @addtogroup EGmemSlab */
/** @{ */
/* ========================================================================= */
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h> /* For printf format support */

#include "eg_memslab.h"

#include "logging-private.h"

/* ========================================================================= */
int EGmemSlabPoolSetParam(EGmemSlabPool_t*const pool,
													const int param,
													const int val)
{
	int rval = 0;
	switch(param)
	{
		case EG_MSLBP_FREEFREE:
			__EGmspLock(pool);
			if(val) pool->freefree = 1;
			else pool->freefree = 0;
			__EGmspUnlock(pool);
			break;
		default:
			rval = 1;
			MESSAGE(0,"Unknown parameter %d",param);
			break;
	}
	EG_RETURN(rval);
}
/* ========================================================================= */
void EGmemSlabDisplay(const EGmemSlab_t*const slab)
{
	const size_t n_elem = slab->control.pool ? slab->control.pool->n_elem : (size_t)0;
	const size_t n_elem2 = (n_elem/8)*8;
	register size_t i;
	QSlog("Slab %p:", (const void*const)slab);
	QSlog("\t->base     : %8p", (void*)(slab->control.base));
	QSlog("\t->elem_sz  : %8zd", slab->control.elem_sz);
	QSlog("\t->n_elem   : %8zd", slab->control.n_elem);
	QSlog("\t->slab_cn  : [%8p,%8p]", 
							(void*)(slab->control.slab_cn.prev), 
							(void*)(slab->control.slab_cn.next));
	QSlog("\t->pool     : %8p", (void*)(slab->control.pool));
	QSlog("\t->next     : %8zd", slab->control.next);
	QSlog("\t->next_list:");
	for( i = 0 ; i < n_elem2 ; i+= 8)
	{
		QSlog("\t[%3zu]=%3u [%3zu]=%3u [%3zu]=%3u [%3zu]=%3u "
								"[%3zu]=%3u [%3zu]=%3u [%3zu]=%3u [%3zu]=%3u",
								i, ((unsigned)(slab->next_list[i])),
								i+1, ((unsigned)(slab->next_list[i+1])),
								i+2, ((unsigned)(slab->next_list[i+2])), 
								i+3, ((unsigned)(slab->next_list[i+3])),
								i+4, ((unsigned)(slab->next_list[i+4])),
								i+5, ((unsigned)(slab->next_list[i+5])),
								i+6, ((unsigned)(slab->next_list[i+6])),
								i+7, ((unsigned)(slab->next_list[i+7])));
	}
	QSlog("\t");
	for( ; i < n_elem ; i++)
	{
		QSlog("[%3zu]=%3u ",i, ((unsigned)(slab->next_list[i])));
	}
}
/* ========================================================================= */
void __EGmemSlabInit( EGmemSlab_t*const slab,
										EGmemSlabPool_t*const Pool)
{
	const EGconstructor_f _EGconstr = Pool->constr;
	const size_t elem_sz = Pool->elem_sz;
	const size_t n_elem = Pool->n_elem;
	register size_t i;
	char*base = (char*)(EG_MEM_ALIGN(sizeof(EGmsbControl_t)+Pool->n_elem) +
											(char*)(slab) + (Pool->c_color));
	slab->control.base = base;
	slab->control.elem_sz = elem_sz;
	slab->control.n_elem = 0;
	EGeListAddAfter(&(slab->control.slab_cn),&(Pool->empty));
	slab->control.pool = Pool;
	slab->control.next = 0;
	/* now we initialize all elements and list of next elements */
	for( i = 0 ; i < n_elem; i++)
	{
		slab->next_list[i] = (uint8_t)(i+1);
		if(_EGconstr)
		{
			_EGconstr(base);
			base += elem_sz;
		}
	}
	slab->next_list[n_elem-1] = EG_SLAB_ENDMARK;
	/* if we are profiling, update here */
	#if EG_SLAB_PROFILE <= DEBUG
	slab->control.pool->n_allocs++;
	slab->control.pool->n_slabs++;
	if(slab->control.pool->max_slabs < slab->control.pool->n_slabs)
		slab->control.pool->max_slabs = slab->control.pool->n_slabs;
	#endif
	/* change current color in main pool */
	Pool->c_color += (uint8_t)(EG_MEM_ALIGNMENT);
	if(Pool->c_color > Pool->max_color) Pool->c_color = 0;
	if(EG_SLAB_VERBOSE <= DEBUG)
	{
		QSlog("Initializing slab as:");
		EGmemSlabDisplay(slab);
	}
}
/* ========================================================================= */
void EGmemSlabClear( EGmemSlab_t*slab)
{
	const EGdestructor_f _EGdest = slab->control.pool->dest;
	const size_t elem_sz = slab->control.elem_sz;
	const size_t n_elem = slab->control.pool->n_elem;
	char*base = slab->control.base;
	register size_t i = n_elem;
	if(EG_SLAB_VERBOSE <= DEBUG)
	{
		QSlog("slab before clearing:");
		EGmemSlabDisplay(slab);
	}
	WARNINGL( EG_SLAB_DEBUG, slab->control.n_elem, 
						"Clearing slab at %p with %zd elements of size %zd",
						(void*)slab, slab->control.n_elem, elem_sz);
	/* clear each element */
	if(_EGdest)
	{
		while(i--)
		{
			_EGdest(base);
			base += elem_sz;
		}
	}
	/* remove the slab from it's containing list */
	EGeListDel(&(slab->control.slab_cn));
	/* if we are profiling, update here */
	#if EG_SLAB_PROFILE <= DEBUG
	slab->control.pool->n_slabs--;
	#endif
	/* if we are debugging, poison all data */
	if(EG_SLAB_DEBUG <= DEBUG)
	{
		slab->control.base = (char*)EG_SLAB_POISON;
		slab->control.slab_cn = (EGeList_t){(void*)EG_SLAB_POISON,(void*)EG_SLAB_POISON};
		slab->control.next = EG_SLAB_ENDMARK;
		for( i = n_elem ; i-- ; ) slab->next_list[i] = EG_SLAB_ENDMARK;
	}
	/* and display if needed */
	if(EG_SLAB_VERBOSE <= DEBUG)
	{
		QSlog("slab after clearing:");
		EGmemSlabDisplay(slab);
	}
	if(EG_SLAB_DEBUG <= DEBUG)
		slab->control.pool = (EGmemSlabPool_t*)EG_SLAB_POISON;
}
/* ========================================================================= */
/** @brief Initialize the profiling data of a slab pool */
#if EG_SLAB_PROFILE <= DEBUG
#define EGmemSlabPoolInitProfile(_EGmPl, __sz, __file, __func, __line) ({\
	_EGmPl->file = __file;\
	_EGmPl->func = __func;\
	_EGmPl->line = __line;\
	_EGmPl->real_sz = __sz;\
	_EGmPl->n_slabs = _EGmPl->n_tot = _EGmPl->max_tot = _EGmPl->max_slabs = \
	_EGmPl->ncals = _EGmPl->n_allocs = 0;})
#else
#define EGmemSlabPoolInitProfile(_EGmPl, __sz, __file, __func, __line)
#endif

/* ========================================================================= */
void __EGmemSlabPoolInit( EGmemSlabPool_t*const pool,
													const size_t sz,
													EGconstructor_f constr_fn,
													EGdestructor_f dest_fn,
													const char*const file,
													const char*const func,
													const int line)
{
	const size_t elem_sz = sz < EG_SLAB_LLIMIT ? EG_SLAB_LLIMIT: EG_MEM_ALIGN(sz);
	const size_t n_elem = (EG_SLAB_SIZE - EG_MEM_ALIGN(sizeof(EGmsbControl_t)))/
												(elem_sz +1);
	/* check that the real element size is within bounds */
	if(elem_sz > EG_SLAB_ULIMIT)
	{
		QSlog("ERROR: Trying to initializate slab pool with element size"
								" %zd > %zd (hard upper limit)", elem_sz, EG_SLAB_ULIMIT);
		exit(EXIT_FAILURE);
	}
	/* initialize the structure */
	#if HAVE_EG_THREAD
	pthread_mutex_init(&(pool->mt),0);
	#endif
	__EGmspLock(pool);
	EGeListInit(&(pool->half));
	EGeListInit(&(pool->empty));
	EGeListInit(&(pool->full));
	pool->constr = constr_fn;
	pool->dest = dest_fn;
	pool->elem_sz = (uint16_t)(elem_sz);
	pool->n_elem = (uint8_t)(n_elem);
	pool->c_color = 0;
	pool->max_color = ((uint8_t)(EG_SLAB_SIZE - EG_MEM_ALIGN(sizeof(EGmsbControl_t)+ n_elem) - 
											(elem_sz*n_elem)));
	pool->freefree = 1;
	EGmemSlabPoolInitProfile(pool,sz, file, func, line);
	/* verbose output */
	if(EG_SLAB_VERBOSE <= DEBUG )
		EGmemSlabPoolDisplay(pool);
	__EGmspUnlock(pool);
}
/* ========================================================================= */
void EGmemSlabPoolClear(EGmemSlabPool_t*const Pool)
{
	void* _EGptr;
	__EGmspLock(Pool);
	while(!EGeListIsEmpty(&(Pool->half)))
	{
		_EGptr = EGmemSlabGetSlab(Pool->half.next);
		EGmemSlabClear(_EGptr);
		free(_EGptr);
	}
	while(!EGeListIsEmpty(&(Pool->empty)))
	{
		_EGptr = EGmemSlabGetSlab(Pool->empty.next);
		EGmemSlabClear(_EGptr);
		free(_EGptr);
	}
	while(!EGeListIsEmpty(&(Pool->full)))
	{
		_EGptr = EGmemSlabGetSlab(Pool->full.next);
		EGmemSlabClear(_EGptr);
		free((void*)EG_SLAB_PAGE(_EGptr));
	}
	/* verbose output */
	if(EG_SLAB_VERBOSE <= DEBUG || EG_SLAB_PROFILE <= DEBUG)
	{
		QSlog("After clearing slab pool:");
		EGmemSlabPoolDisplay(Pool);
	}
	__EGmspUnlock(Pool);
}
/* ========================================================================= */
void EGmemSlabPoolShrink(EGmemSlabPool_t*const Pool)
{
	void* _EGptr;
	__EGmspLock(Pool);
	while(!EGeListIsEmpty(&(Pool->empty)))
	{
		_EGptr = EGmemSlabGetSlab(Pool->empty.next);
		EGmemSlabClear(_EGptr);
		free((void*)EG_SLAB_PAGE(_EGptr));
	}
	__EGmspUnlock(Pool);
}
/* ========================================================================= */
void EGmemSlabPoolDisplay(const EGmemSlabPool_t*const pool)
{
	QSlog("Pool %p:", (const void*const)pool);
	QSlog("\t->half      : [%8p,%8p]",
							(void*)(pool->half.prev), 
							(void*)(pool->half.next));
	QSlog("\t->empty     : [%8p,%8p]",
							(void*)(pool->empty.prev), 
							(void*)(pool->empty.next));
	QSlog("\t->full      : [%8p,%8p]",
							(void*)(pool->full.prev), 
							(void*)(pool->full.next));
	QSlog("\t->constr    : %8p", (void*)(pool->constr));
	QSlog("\t->dest      : %8p", (void*)(pool->dest));
	QSlog("\t->elem_sz   : %8"PRIu16, pool->elem_sz);
	QSlog("\t->n_elem    : %8"PRIu8, pool->n_elem);
	QSlog("\t->c_color   : %8"PRIu8, pool->c_color);
	QSlog("\t->max_color : %8"PRIu8, pool->max_color);
	QSlog("\t->freefree  : %8"PRIu8, pool->freefree);
	#if EG_SLAB_PROFILE <= DEBUG
	if(pool->ncals)
	{
		QSlog("\t->file      : %s", pool->file);
		QSlog("\t->func      : %s", pool->func);
		QSlog("\t->line      : %8d", pool->line);
		QSlog("\t->real_sz   : %8"PRIu64, pool->real_sz);
		QSlog("\t->n_slabs   : %8"PRIu64, pool->n_slabs);
		QSlog("\t->n_tot     : %8"PRIu64, pool->n_tot);
		QSlog("\t->max_tot   : %8"PRIu64, pool->max_tot);
		QSlog("\t->max_slabs : %8"PRIu64, pool->max_slabs);
		QSlog("\t->ncals     : %8"PRIu64, pool->ncals);
		QSlog("\t->n_allocs  : %8"PRIu64, pool->n_allocs);
		QSlog("\tEficiency   :");
		QSlog("\t\talloc ratio  : %8lg (%"PRIu64"/%"PRIu64")", 
								((double)pool->ncals)/pool->n_allocs, pool->ncals, pool->n_allocs);
		QSlog("\t\tmemory waste : %8.3lf %%",
								100*(1-((double)pool->max_tot*pool->real_sz)/
										 ((double)pool->max_slabs*EG_SLAB_SIZE)));
	}
	#endif
}
/* ========================================================================= */
