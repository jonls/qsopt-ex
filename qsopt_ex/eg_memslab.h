/* EGlib "Efficient General Library" provides some basic structures and
 * algorithms commons in many optimization algorithms.
 *
 * Copyright (C) 2005,2006,2007,2008,2009,2010,2011 Daniel Espinoza.
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
/** @defgroup EGmemSlab EGmemSlab
 *
 * This is a basic interface for slab pool managment. The idea comes from Slabs
 * as defined in both Linux and Solaris (see "The Slab Allocator: An
 * Object-Caching Kernel Memory Allocator", by Jeff Bonwick, Sun Microsystems),
 * the basic idea is
 * to provide pool for a specific type of object, and to store them in an
 * initialized state, so that initialization and destruction only is done while
 * growing/freeing the memory slabs, thus this approach should provide greater
 * advantages for complitated to initialize structures. and in theory
 * (althought not yet implemented) this structure can be managed so as to
 * provide a shrinkable memory managment on the fly.
 *
 * In this implementation we only allow small caches (i.e. objects must be
 * smaller than EG_SLAB_ULIMIT and internally we don't allocate objects smaller
 * than EG_SLAB_LLIMIT), 
 * within a unique memory page. We could allow in the future for more flexible
 * slabs. This implementation also uses colored slabs (see the paper for
 * further details).
 *
 * Here we can see a schematic drawing of the slab allocator structure and
 * functions:
 *
 * @version 0.9.0
 * @par History:
 * - 2011-03-01
 * 						- Make EGms_t thread-safe, i.e. a single memory pool can serve
 * 							several threads. This makes the default memory manager for GMP
 * 							thread safe as well
 * - 2010-12-27
 * 						- Add short-hand names for functions
 * - 2008-10-06
 * 						- Second implementation
 * - 2005-07-30
 * 						- First Implpementation.
 * */
/** @file
 * @ingroup EGmemSlab */
/** @addtogroup EGmemSlab */
/**  @{  */
/** @example eg_memslab.ex.c */
/* ========================================================================= */
#ifndef __EG_MEM_SLAB_H__
#define __EG_MEM_SLAB_H__

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <unistd.h>

#include "eg_mem.h"
#include "eg_elist.h"

/* ========================================================================= */
/** @name EGmemSlab Parameters
 * @brief parameters for controling slab pool allocation */
/* @{ */
/** @brief control the handle of empty slabs, if set to non-zero, free unused
 * slabs as they become unused, if set to zero, keep all unused slabs until
 * #EGmemSlabPoolShrink is called */
#define EG_MSLBP_FREEFREE 1
/* @} */
/* ========================================================================= */
/* declare the slab pool structure */
struct EGmemSlabPool_t;

/* ========================================================================= */
/** @brief maximum size of the objects that can be allocated via slab. */
#define EG_SLAB_ULIMIT ((size_t)1023)
#define EG_SLAB_LLIMIT ((size_t)16)

/* ========================================================================= */
/** @brief size of the memory slabs (in bytes ) */
#define EG_SLAB_SIZE ((size_t)0x1000)

/* ========================================================================= */
/** @brief mask to detect the position of a piece of memory within a slab */
#define EG_SLAB_MASK (~(EG_SLAB_SIZE-1))

/* ========================================================================= */
/** @brief address used to check consistency if enabled */
#define EG_SLAB_POISON ((size_t)0xdeadbeef)

/* ========================================================================= */
/** @brief if set to one, enable profiling for the slab allocator */
#define EG_SLAB_PROFILE 1000

/* ========================================================================= */
/** @brief local verbose level for the slab allocator, the lower the level, the
 * more information will be printed on screen. */
#define EG_SLAB_VERBOSE 1000

/* ========================================================================= */
/** @brief local debug level for the slab allocator, the lower the level, the
 * more testing will be done. */
#define EG_SLAB_DEBUG 1000

/* ========================================================================= */
/** @brief end of list marker, note that this can not be more than 255 */
#define EG_SLAB_ENDMARK ((uint8_t)255U)

#ifndef EG_SLAB_REDUCE_TO_MALLOC
/* ========================================================================= */
/** @brief if set to one, reduce the slab pool allocator to a simple malloc
 * call */
#define EG_SLAB_REDUCE_TO_MALLOC 0
#endif


/* ========================================================================= */
/** @brief Given a pointer, return a pointer to the beginning of the containing
 * page. */
#define EG_SLAB_PAGE(__ptr) (((size_t)__ptr)&EG_SLAB_MASK)

/* ========================================================================= */
/** @brief structure that holds the information relevant to each slab
 * */
typedef struct
{
	char*base;								/**< @brief Where the data-payload start */
	size_t elem_sz;						/**< @brief byte-size of each element */
	size_t n_elem;						/**< @brief number of used elements */
	EGeList_t slab_cn;				/**< @brief Connector into the list of slabs*/
	struct EGmemSlabPool_t *pool;	/**< Pointer to the slab pool structure */
	size_t next;							/**< @brief next free element */
} EGmsbControl_t;
typedef struct
{
	EGmsbControl_t control;	/**< @brief comon base structure for slabs */
	uint8_t next_list[];		/**< @brief list of free elements, the next 
																element is next_list[0], whenever we reach a
																value of 255, it is the end of the free-list.
																*/
} EGmemSlab_t;

/* ========================================================================= */
/** @brief structure used to store a slab memory pool */
typedef struct EGmemSlabPool_t
{
	EGeList_t half;					/**< Head of the list for half-full slabs */
	EGeList_t empty;				/**< Head of the list for non used slabs */
	EGeList_t full;					/**< Head of the list for fully used slabs*/
	EGconstructor_f constr;	/**< Constructor for the local elements */
	EGdestructor_f dest;		/**< Destructor for the local elements */
	uint16_t elem_sz;				/**< Size of the elements in the slab, including 
															 extra space for pointer to next. */
	uint8_t n_elem;					/**< Total number of elements in each slab */
	uint8_t c_color;				/**< Last used color while creating slabs. */
	uint8_t max_color;			/**< Maximum valid value for colors in this pool */
	uint8_t freefree:1;			/**< if non-zero, free non-used slabs */
	uint8_t pad1:7;					/**< padding */
	uint16_t pad2;					/**< padding */
	char const *file;				/**< File where the structure was initialized */
	char const *func;				/**< Function where the structure was initialized */
	int line;								/**< Line where the structure was initialized */
	uint64_t real_sz;				/**< Actual size of the elements asked by the user */
	uint64_t n_slabs;				/**< Number of slabs */
	uint64_t n_tot;					/**< Total number of elements in use by the user */
	uint64_t max_tot;				/**< Maximum number of elements allocated */
	uint64_t max_slabs;			/**< Maximum number of slabs used */
	uint64_t ncals;					/**< number to alloc calls */
	uint64_t n_allocs;			/**< number slab alloc calls */
	#if HAVE_EG_THREAD
	pthread_mutex_t mt;			/**< mutex for memory slab manager */
	#endif
}
EGmemSlabPool_t;

/* ========================================================================= */
/** @brief slab lock/unlock macros */
#if HAVE_EG_THREAD
#define __EGmspLock(__slab) pthread_mutex_lock(&(__slab->mt))
#define __EGmspUnlock(__slab) pthread_mutex_unlock(&(__slab->mt))
#else
#define __EGmspLock(__slab)
#define __EGmspUnlock(__slab)
#endif
/* ========================================================================= */
/** @brief display given slab structure 
 * @param slab what to display
 * */
void EGmemSlabDisplay(const EGmemSlab_t*const slab);

/* ========================================================================= */
/** @brief display given pool structure
 * @param pool what to display
 * */
void EGmemSlabPoolDisplay(const EGmemSlabPool_t*const pool);

/* ========================================================================= */
/** @brief given a piece of memory that should be within a slab, return the
 * pointer to the related slab structure, remember that the slab structure is
 * at the beggining of the page. */
#define EGmemSlabGetSlab(__ptr) \
		((EGmemSlab_t*)(EG_SLAB_PAGE(__ptr)))

/* ========================================================================= */
/** @brief initialize a slab structure. This include calling the constructor
 * for all elements in the slab. Note that this function asumes that all memory
 * has been previously set to NULL. It will also place this slab in the list 
 * of empty slabs in the given pool.
 * @param slab pointer within the memory range of the slab to be initialized.
 * @param __Pool Slab __Pool where this slab will bellong from. The slab pool
 * should be initialized (i.e. should have a constructor and destructor, and an
 * element size set.
 * */
void __EGmemSlabInit( EGmemSlab_t*const slab,
										EGmemSlabPool_t*const __Pool);
/* ========================================================================= */
/** @brief given an initialized slab, clear all internally allocated memory,
 * and leave the slab ready to be freed by 'free', this include calling the
 * destructor for all elements in the slab.
 * @param slab pointer to an area of the slab memory to be clear.
 * @note If debugging is enabled, then all fields will be poisoned so that
 * subsequent use of this structure will fail (but for the free call). Also, if
 * debugging is enabled, we will check that the slab has no element in use.
 * */
void EGmemSlabClear( EGmemSlab_t*const slab);

/* ========================================================================= */
/** @brief Given a non-full slab, extract a pointer to the next unused element
 * in the slab, and update all internal data. and if it becomes full, then move
 * it to the full list within the pool. Also, if debugging, poison the pointer
 * to the enext element in the returned element. If the slab is not full, and
 * the number of active elements is one, then move the slab to the half-full
 * slab list in the pool.
 * @param slab pointer within a slab memory.
 * @return pointer to a void* of initialized memory by the given contructor in
 * the slab pool.
 * */
#define EGmemSlabPopElement(slab) ({\
	EGmemSlab_t*const _EGmSlb = EGmemSlabGetSlab(slab);\
	EGmemSlabPool_t*const _EGPlRf = _EGmSlb->control.pool;\
	const size_t _EGmSlb_esz = _EGmSlb->control.elem_sz;\
	const size_t _EGmSlb_ne = _EGmSlb->control.next;\
	const size_t _EGmSlb_nn = _EGmSlb->next_list[_EGmSlb_ne];\
	void*const _EGelem = (void*)(_EGmSlb_ne*_EGmSlb_esz + _EGmSlb->control.base);\
	/* now update the slab */\
	_EGmSlb->control.n_elem++;\
	_EGmSlb->control.next = _EGmSlb_nn;\
	EXITL(EG_SLAB_DEBUG,_EGmSlb_ne == EG_SLAB_ENDMARK, "Allocating from full slab");\
	/* if the slab is full, move it to full list */\
	if(_EGmSlb_nn == EG_SLAB_ENDMARK){\
		EGeListMoveAfter(&(_EGmSlb->control.slab_cn),&(_EGPlRf->full));}\
	/* if the slab is first-time used, move to half list */\
	else if(_EGmSlb->control.n_elem == 1U){\
		EGeListMoveAfter(&(_EGmSlb->control.slab_cn),&(_EGPlRf->half));}\
	/* return the element */\
	_EGelem;})

/* ========================================================================= */
/** @brief Given an used object within a slab, give it back to the slab for
 * future use.
 * @param __ptr pointer to the element to be given back to its containing slab.
 * */
#define EGmemSlabPushElement(__ptr) do{\
	char*const _EGmsbPtr = (char*)(__ptr);\
	EGmemSlab_t*const _EGmSlb = EGmemSlabGetSlab(_EGmsbPtr);\
	EGmemSlabPool_t*const _EGPlRf = _EGmSlb->control.pool;\
	__EGmspLock(_EGPlRf);{\
	const size_t _EGmSlb_esz = _EGmSlb->control.elem_sz;\
	const size_t _EGmSlb_ne = _EGmSlb->control.next;\
	const size_t _EGmSlb_nn = ((_EGmsbPtr) - _EGmSlb->control.base)/_EGmSlb_esz;\
	/* if debugging, check for poison in the pointer to the next element in the \
	 * given element */\
	EXITL(EG_SLAB_DEBUG, !_EGmSlb->control.n_elem, "freeing from an empty slab");\
	/* now actually put the element into the slab */\
	_EGmSlb->control.n_elem--;\
	_EGmSlb->control.next = _EGmSlb_nn;\
	_EGmSlb->next_list[_EGmSlb_nn] = _EGmSlb_ne;\
	__EGmsbUPD2(_EGPlRf);\
	/* if the slab is now not being used, update accordingly */\
	if(_EGmSlb_ne == EG_SLAB_ENDMARK){\
		EGeListMoveAfter(&(_EGmSlb->control.slab_cn),&(_EGPlRf->half));}\
	else if(!_EGmSlb->control.n_elem){\
		if(_EGPlRf->freefree){\
			EGmemSlabClear(_EGmSlb);\
			free((void*)_EGmSlb);}\
		else{\
			EGeListMoveAfter(&(_EGmSlb->control.slab_cn),&(_EGPlRf->empty));}}}\
	__EGmspUnlock(_EGPlRf);\
	}while(0)

/* ========================================================================= */
/** @brief initialize a slab pool as an empty pool for elements of the given
 * size, and with te given constructor and destructors. 
 * @param constr_fn constructor fnctioin for the elements to be stored in the
 * pool.
 * @param dest_fn destructor function for the elements to be stored in the
 * pool.
 * @param pool pointer to the slab pool to initialize.
 * @param sz (real) size (in bytes) of the elements to be hold. in the pool.
 * This means that sz is the result of sizeof(TYPE), where TYPE is the
 * structure to be pooled. */
#define EGmemSlabPoolInit(pool, sz, constr_fn, dest_fn) \
				__EGmemSlabPoolInit(pool, sz, constr_fn, dest_fn, \
														__FILE__, __func__, __LINE__)
void __EGmemSlabPoolInit( EGmemSlabPool_t*const pool,
													const size_t sz,
													EGconstructor_f constr_fn,
													EGdestructor_f dest_fn,
													const char*const file,
													const char*const func,
													const int line);
/* ========================================================================= */
/** @brief clear a slab pool and all internal sub-structures and data, no
 * further calls to this structure are posible after this (but for freeing the
 * memory containing this data, or to re-initialize it).
 * @param __Pool slab pool to be cleared.
 * */
void EGmemSlabPoolClear(EGmemSlabPool_t*const __Pool);
/* ========================================================================= */
/** @brief add one to the given pointer, if profiling is enabled, otherwise, do
 * nothing */
#if EG_SLAB_PROFILE <= DEBUG
#define __EGmsbUPD1(__pool) do{\
	__pool->ncals++;\
	__pool->n_tot++;\
	if(__pool->max_tot < __pool->n_tot) __pool->max_tot = __pool->n_tot;}while(0)
#define __EGmsbUPD2(__pool) __pool->n_tot--
#else
#define __EGmsbUPD1(__pool)
#define __EGmsbUPD2(__pool)
#endif

/* ========================================================================= */
/** @brief Given a slab pool, return an element from the pool. 
 * @param __Pool slab pool from where we will get the memory. 
 * @return pointer to an initialize element. */
#if EG_SLAB_REDUCE_TO_MALLOC
#define EGmemSlabPoolAlloc(__Pool) ({\
	EGmemSlabPool_t*const _EGmPl = (__Pool);\
	void*_EGmb = EGmalloc(_EGmPl->elem_sz+sizeof(void*));\
	void**_EGpt = (void**)_EGmb;\
	__EGmsbUPD1(_EGmPl);\
	(*_EGpt) = _EGmPl;\
	_EGmb=((void*)(_EGpt+1));\
	if(_EGmPl->constr) _EGmPl->constr(_EGmb);\
	_EGmb;})
#else
#define EGmemSlabPoolAlloc(__Pool) ({\
	EGmemSlabPool_t*const _EGmPl = (__Pool);\
	void* _EGSmbRf = 0;\
	void* _EGrptr = 0;\
	int __EGmPlerr=0;\
	__EGmspLock(_EGmPl);\
	__EGmsbUPD1(_EGmPl);\
	if(!EGeListIsEmpty(&(_EGmPl->half))){ _EGSmbRf = _EGmPl->half.next;}\
	else if(!EGeListIsEmpty(&(_EGmPl->empty))){ _EGSmbRf = _EGmPl->empty.next;}\
	else{\
		if((__EGmPlerr=posix_memalign(&_EGSmbRf,EG_SLAB_SIZE,EG_SLAB_SIZE))){\
			EXIT(1,"posix_memalign falied with code %d, error %s",__EGmPlerr,\
					strerror(__EGmPlerr));}\
		__EGmemSlabInit(_EGSmbRf,_EGmPl);}\
	_EGrptr = EGmemSlabPopElement(_EGSmbRf);\
	MESSAGE(EG_SLAB_VERBOSE,"Returning %p",_EGrptr);\
	__EGmspUnlock(_EGmPl);\
	_EGrptr;})
#endif

/* ========================================================================= */
/** @brief Given a pointer to an element allocated through a slab pool, give it
 * back to the pool.
 * @param __ptr pointer to be returned to the pool.
 * */
#if EG_SLAB_REDUCE_TO_MALLOC
#define EGmemSlabPoolFree(__ptr) ({\
	void**_EGptr = ((void**)(__ptr))-1;\
	EGmemSlabPool_t*const _EGmPl = (EGmemSlabPool_t*)(*_EGptr);\
	__EGmsbUPD2(_EGmPl);\
	if(_EGmPl->dest) _EGmPl->dest(((void*)(_EGptr+1)));\
	EGfree(_EGptr);})
#else
#define EGmemSlabPoolFree(__ptr) EGmemSlabPushElement(__ptr)
#endif

/* ========================================================================= */
/** @brief Given a slab pool, free all unused slabs 
 * @param pool slab pool to be shrinked. */
void EGmemSlabPoolShrink(EGmemSlabPool_t*const pool);

/* ========================================================================= */
/** @brief set parameters for a slab pool */
int EGmemSlabPoolSetParam(EGmemSlabPool_t*const pool,
													const int param,
													const int val);
/* ========================================================================= */
/** @brief short names for common functions */
#define EGmsAlloc  EGmemSlabPoolAlloc
#define EGmsFree   EGmemSlabPoolFree
#define EGmsInit   EGmemSlabPoolInit
#define EGmsClear  EGmemSlabPoolClear
#define EGms_t     EGmemSlabPool_t
#define EGmsShrink EGmemSlabPoolShrink
/* ========================================================================= */
/* end of eg_memslab.h */
/**  @}  */
#endif
