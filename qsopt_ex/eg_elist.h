/* EGlib "Efficient General Library" provides some basic structures and
 * algorithms commons in many optimization algorithms.
 *
 * Copyright (C) 2005 Daniel Espinoza and Marcos Goycoolea.
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
/** @defgroup EGeList EGeList
 *
 * Here we define the basic interface for a circular linked list where the list
 * is embeded in some other structure. The ideas come from the Linux Kernel
 * implementation of lists. This implementation is based on the philosophy of 
 * embeded structures.
 *
 * 
 * @version 0.0.1
 * @par History:
 * - 2005-08-19
 * 						- Add debugging control
 * - 2005-05-23
 * 						- First Implementation.
 *
 * @note In general, the functions described bellow don't perform consistency
 * checks. It is asumed that the user does know what is he doing.
 *
 * @note If you want to have some debugging control try changing the debug level
 * at compile time, and lowering the debug level asociated to the list function
 * as defined in eg_configure.h.
 *
 * */
/** @file 
 * @ingroup EGeList */
/** @addtogroup EGeList */
/** @{ */
/** @example eg_elist.ex.c 
 * This is a working (althought useless) example on @ref EGeList.
 * */
/* ========================================================================= */
#ifndef __EG_ELIST_H__
#define __EG_ELIST_H__

#include "eg_macros.h"

/* ========================================================================= */
/** @brief debug level for lists */
#define __EL_DEBUG_ 100
/* ========================================================================= */
/** @brief Null-initialized embeded list */
#define EGeListNull ((EGeList_t){0,0})

/* ========================================================================= */
/** 
 * @brief List Node Structure.
 * @par Description:
 * This structure is to store a general node of the list. It is composed by
 * two members, that point to the next and previous structures in the list.  */
typedef struct EGeList_t
{
	struct EGeList_t *next;/**< Pointer to the next structure in the list */
	struct EGeList_t *prev;/**< Pointer to the previous structure in the list */
}
EGeList_t;

/* ========================================================================= */
/** @brief Initialize a given structure to point to itself (in circular
 * fashion).
 * @param __lpt pointer to the list to initialize. 
 * @return the pointer to the list. */
#define EGeListInit(__lpt) ({\
	EGeList_t*const __EGeL_init =(__lpt);\
	__EGeL_init->next = __EGeL_init->prev = __EGeL_init;})

/* ========================================================================= */
/** @brief Insert a __newpt __entry between two known consecutive entries.
 * @par Description:
 * This is only for internal list manipulation, where we know the prev/next
 * entries already.
 * @param __newpt pointer to the list node to insert.
 * @param __prevpt pointer to the node to preceed the __newpt node.
 * @param __nextpt pointer to the node to follow the __newpt node. 
 * @return the address of __newpt.
 * */
#define __EGeListAdd(__newpt,__prevpt,__nextpt) ({\
	EGeList_t*const __EGeL_add_new = (__newpt);\
	EGeList_t*const __EGeL_add_prev = (__prevpt);\
	EGeList_t*const __EGeL_add_next = (__nextpt);\
	__EGeL_add_next->prev = __EGeL_add_new;\
	__EGeL_add_prev->next = __EGeL_add_new;\
	__EGeL_add_new->next = __EGeL_add_next;\
	__EGeL_add_new->prev = __EGeL_add_prev;\
	__EGeL_add_new;})

/* ========================================================================= */
/** @brief Insert a __newpt __entry after the given pointer.
 * @param __newpt pointer to the __newpt list node to insert.
 * @param __head pointer from where the __newpt __entry will follow.
 * @return the pointer to the __newpt __entry in the list. 
 * */
#define EGeListAddAfter(__newpt,__head) __EGeListAdd(__newpt,__head,(__head)->next)

/* ========================================================================= */
/** @brief Insert a __newpt __entry before the given pointer.
 * @param __newpt pointer to the __newpt list node to insert.
 * @param __tailpt pointer that will follow the __newpt __entry in the list.
 * @return the pointer to the __newpt __entry in the list. 
 * */
#define EGeListAddBefore(__newpt,__tailpt) __EGeListAdd(__newpt,(__tailpt)->prev,__tailpt)

/* ========================================================================= */
/** @brief Given two nodes, link them as if they would follow one another in the
 * list (used to delete points from a list).
 * @param __prevpt pointer to the guy to be in first in the list.
 * @param __nextpt pointer to the guy to follow in the list.
 * @par Description:
 * This function is intended to be used only internally, where we know what is
 * what, if you use it is because you also know what is going on.
 * */
#define __EGeListLink(__prevpt,__nextpt) ({\
	EGeList_t* __EGeL_lnk_prev = (__prevpt);\
	EGeList_t* __EGeL_lnk_next = (__nextpt);\
	__EGeL_lnk_prev->next = __EGeL_lnk_next;\
	__EGeL_lnk_next->prev = __EGeL_lnk_prev;\
	0;})

/* ========================================================================= */
/** @brief Given a node, eliminate it from the list it bellongs. but don't
 * change the internal data in the eliminated list (be carefull, if you will 
 * use it afterwards, then you MUST initialize it). If debugging is enabled,
 * then whenever you delete, the connector is reseted to 0xffffffff. What you
 * can count on is that the connector won't be NULL after deleting it from the
 * list, but it's values may be lost if we are debugging.
 * @param __entry __entry to eliminate from the list.
 * @return pointer to the deleted __entry from the list.*/
#define EGeListDel(__entry) ({\
	EGeList_t *const __EGeL_del_entr = (__entry);\
	__EGeListLink(__EGeL_del_entr->prev,__EGeL_del_entr->next);\
	if(__EL_DEBUG_ <= DEBUG) \
		(*__EGeL_del_entr) = (EGeList_t){(EGeList_t*)0xffffffffU,\
																		 (EGeList_t*)0xffffffffU};\
	__EGeL_del_entr;})

/* ========================================================================= */
/** @brief Replace one __entry with another in a list.
 * @param __oldpt __entry to be replaced, note that the pointers stored in next/prev
 * won't be changed, this may possible lead to errors if the __entry is used
 * afterwards without initialization.
 * @param __newpt __newpt __entry in the list.
 * @return pointer to the old replaced member.
 * */
#define EGeListReplace(__oldpt,__newpt) ({\
	EGeList_t* __EGeL_rep_old = (__oldpt);\
	EGeList_t* __EGeL_rep_new = (__newpt);\
	__EGeL_rep_new->next = __EGeL_rep_old->next;\
	__EGeL_rep_new->prev = __EGeL_rep_old->prev;\
	__EGeL_rep_new->next->prev = __EGeL_rep_new;\
	__EGeL_rep_new->prev->next = __EGeL_rep_new;\
	__EGeL_rep_old;})

/* ========================================================================= */
/** @brief Move an element from one list to another (deleting it from the
 * original one).
 * @param __entry element to be removed from it's current list to a position 
 * after the given __head.
 * @param __head element to be before the moved element.
 * */
#define EGeListMoveAfter(__entry,__head) ({\
	__EGeListLink((__entry)->prev,(__entry)->next);\
	EGeListAddAfter(__entry,__head);})

/* ========================================================================= */
/** @brief Move an element from one list to another (deleting it from the
 * original one).
 * @param __entry element to be removed from it's current list to a position
 * before the given __tailpt.
 * @param __tailpt element to be after the moved element.
 * */
#define EGeListMoveBefore(__entry,__tailpt) ({\
	__EGeListLink((__entry)->prev,(__entry)->next);\
	EGeListAddBefore(__entry,__tailpt);})

/* ========================================================================= */
/** @brief test whether a list is empty (i.e. he is its own next pointer) */
#define EGeListIsEmpty(__head) ({\
	EGeList_t* __EGeL_emp_head = (__head);\
	(__EGeL_emp_head == __EGeL_emp_head->next);})

/* ========================================================================= */
/** @brief move all elements in one list to the given location in a second list.
 * Note that this function assumes that the list is represented by a pointer to
 * an EGeList_t structure that act as a marker but that don't bellong to the
 * list, and thus is not included in the joinded list.
 * @param __list marker to the list to be joined with the second. Note that the
 * fields in list won't be reinitialized, so be carefull with that, because the
 * fields are pointing to inconsistent data as it is, if you want to reutilize
 * the list you must call #EGeListInit before.
 * @param __head position from where the list will be spliced in.
 * @note Note that the original list is left in an undefined status, so before
 * use it, it should be re-initialized.
 * */
#define __EGeListSplice(__list,__head) ({\
	EGeList_t* __EGeL_spl_list = (__list);\
	EGeList_t* __EGeL_spl_first = __EGeL_spl_list->next;\
	EGeList_t* __EGeL_spl_last = __EGeL_spl_list->prev;\
	EGeList_t* __EGeL_spl_head = (__head);\
	EGeList_t* __EGeL_spl_at = __EGeL_spl_head->next;\
	__EGeL_spl_first->prev = __EGeL_spl_head;\
	__EGeL_spl_head->next = __EGeL_spl_first;\
	__EGeL_spl_last->next = __EGeL_spl_at;\
	__EGeL_spl_at->prev = __EGeL_spl_last;\
	0;})
#define EGeListSplice(__list,__head) ({if(!EGeListIsEmpty(__list)) __EGeListSplice(__list,__head);})

/* ========================================================================= */
/** @}*/
/* end of eg_elist.h */
#endif
