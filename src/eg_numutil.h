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
#ifndef __EG_NUMUTIL_H__
#define __EG_NUMUTIL_H__
#include "qs_config.h"
#include "eg_macros.h"
#include "eg_nummacros.h"
#include "eg_lpnum.h"

/* ========================================================================= */
/** @defgroup EGlpNumUtil General Number Utilities
 * Here we put some utilities common for different number types but thaat we
 * want to implement as templates, like permutation sorting, inner product of
 * vectors, and so-on..
 * 
 * @par History:
 * Revision 0.0.2
 *  - 2007-10-08
 *  					- Separate template file and independet file into eg_nummacros.h
 *  					- Move EGabs, EGswap, EGmin and EGmax to this file
 *  - 2005-10-31
 *  					- First implementation.
 * */
/** @{*/
/** @file
 * @brief This file provide the user interface and function definitions for
 * general number utilities.
 * */
/* ========================================================================= */
/** @brief compute the inner product of two arrays.
 * @param arr1 first array.
 * @param arr2 second array.
 * @param length number of entries to consider in both arrays, from zero to
 * length - 1.
 * @param rop where to store the result.
 * */
#define EGlpNumInnProd(__rop,__arr1,__arr2,__length) __EGlpNumInnProd((&(__rop)),__arr1,__arr2,__length)
/* ========================================================================= */
/** @brief internal version, this is done to avoid using stdc99 and rely on
 * more basic stdc89 */
void __EGlpNumInnProd(EGlpNum_t*rop,EGlpNum_t*const arr1,EGlpNum_t*const arr2, const size_t length);
/* ========================================================================= */
/** @brief Sort (in increasing order) a sub-set of entries in an array using 
 * quicksort, by permutating the order of the elements in the subset rather 
 * than in the whole original array.
 * @param sz length of the permutation array.
 * @param perm array of indices of elements that we want to sort.
 * @param elem array (of length at least max(perm[k]:k=0,...,sz-1)) containing
 * the elements to be sorted.
 * @note The array of elements is not changed by this function.
 * @note This code is based in concorde's implementation of
 * permutation-quick-sort.
 * */
void EGutilPermSort (const size_t sz,
										 int *const perm,
										 const EGlpNum_t * const elem);

/* ========================================================================= */
/** @brief Sort (in decreasing order) a sub-set of entries in an array using 
 * quicksort, by permutating the order of the elements in the subset rather 
 * than in the whole original array.
 * @param sz length of the permutation array.
 * @param perm array of indices of elements that we want to sort.
 * @param elem array (of length at least max(perm[k]:k=0,...,sz-1)) containing
 * the elements to be sorted.
 * @note The array of elements is not changed by this function.
 * @note This code is based in concorde's implementation of
 * permutation-quick-sort.
 * */
void EGutilPermSort2 (const size_t sz,
										 int*const perm,
										 const EGlpNum_t*const elem);

/* ========================================================================= */
/** @}*/
#endif
