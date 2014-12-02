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

#ifndef __SORTRUS_H__
#define __SORTRUS_H__
/****************************************************************************/
/*                                                                          */
/*                             sortrus.c                                    */
/*                                                                          */
/****************************************************************************/

#include "urandom.h"

void ILLutil_int_array_quicksort (
	int *len,
	int n),
  ILLutil_int_perm_quicksort (
	int *perm,
	int *len,
	int n),
  ILLutil_double_perm_quicksort (
	int *perm,
	double *len,
	int n),
	ILLutil_str_perm_quicksort (
	int *perm,
	char **len,
	int n),
  ILLutil_rselect (
	int *arr,
	int l,
	int r,
	int m,
	double *coord,
	ILLrandstate * rstate);

char *ILLutil_linked_radixsort (
	char *data,
	char *datanext,
	char *dataval,
	int valsize);


#endif
