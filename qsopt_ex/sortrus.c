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

/* RCSINFO $Id: sortrus.c,v 1.2 2003/11/05 16:47:22 meven Exp $ */
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
/*                         SORTING ROUTINES                                 */
/*                                                                          */
/*                             TSP CODE                                     */
/*                                                                          */
/*   Written by:  Applegate, Bixby, Chvatal, and Cook                       */
/*   DATE:  February 24, 1994                                               */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  char *ILLutil_linked_radixsort (char *data, char *datanext,              */
/*      char *dataval, int valsize)                                         */
/*    USAGE:                                                                */
/*      head = (bar *) ILLutil_linked_radixsort ((char *) head,              */
/*         (char *) &(head->next), (char *) &(head->val), sizeof (int));    */
/*    Then head is the start of the linked list in increasing order of      */
/*    val, with next as the field that links the bars.                      */
/*    WARNING: DOES NOT HANDLE NEGATIVE NUMBERS PROPERLY.                   */
/*                                                                          */
/*  void ILLutil_int_array_quicksort (int *len, int n)                       */
/*    len - the array to be sorted                                          */
/*    n - the number of elements in len                                     */
/*    Uses quicksort to put len in increasing order.                        */
/*                                                                          */
/*  void ILLutil_int_perm_quicksort (int *perm, int *len, int n)             */
/*    perm - must be allocated and initialized by the calling routine,      */
/*           it will be arranged in increasing order of len.                */
/*    n - the number of elements in perm and len.                           */
/*                                                                          */
/*  void ILLutil_double_perm_quicksort (int *perm, double *len, int n)       */
/*    perm - must be allocated and initialized by the calling routine,      */
/*           it will be arranged in increasing order of len.                */
/*    n - the number of elements in perm and len.                           */
/*                                                                          */
/*  void ILLutil_rselect (int *arr, int l, int r, int m,                     */
/*      double *coord, ILLrandstate *rstate)                                 */
/*    arr - permutation that will be rearranged                             */
/*    l,r - specify the range of arr that we are interested in              */
/*    m - is the index into l,r that is the break point for the perm        */
/*    coord - gives the keys that determine the ordering                    */
/*                                                                          */
/****************************************************************************/

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "qs_config.h"
#include "logging-private.h"

#include "eg_lpnum.h"

#include "util.h"
#include "except.h"


#define BITS_PER_PASS (8)

#define NBINS (1<<BITS_PER_PASS)


static void select_EGlpNum_split (
	int *arr,
	int n,
	EGLPNUM_TYPE * v,
	int *start,
	int *end,
	EGLPNUM_TYPE * coord),
  select_EGlpNum_sort (
	int *arr,
	int n,
	EGLPNUM_TYPE * coord),
  select_EGlpNum_sort_dsample (
	EGLPNUM_TYPE * samp,
	int n);

void EGLPNUM_TYPENAME_ILLutil_EGlpNum_perm_quicksort (
	int *perm,
	EGLPNUM_TYPE * len,
	int n)
{
	int i, j, temp;
	EGLPNUM_TYPE t;

	if (n <= 1)
		return;

	EGLPNUM_TYPENAME_EGlpNumInitVar (t);
	ILL_SWAP (perm[0], perm[(n - 1) / 2], temp);

	i = 0;
	j = n;
	EGLPNUM_TYPENAME_EGlpNumCopy (t, len[perm[0]]);

	for (;;)
	{
		do
			i++;
		while (i < n && EGLPNUM_TYPENAME_EGlpNumIsLess (len[perm[i]], t));
		do
			j--;
		while (EGLPNUM_TYPENAME_EGlpNumIsLess (t, len[perm[j]]));
		if (j < i)
			break;
		ILL_SWAP (perm[i], perm[j], temp);
	}
	ILL_SWAP (perm[0], perm[j], temp);

	EGLPNUM_TYPENAME_EGlpNumClearVar (t);
	EGLPNUM_TYPENAME_ILLutil_EGlpNum_perm_quicksort (perm, len, j);
	EGLPNUM_TYPENAME_ILLutil_EGlpNum_perm_quicksort (perm + i, len, n - i);
}

/**********  Median - Select Routines **********/

/* NSAMPLES should be odd */
#define NSAMPLES 3
#define SORTSIZE 20


void EGLPNUM_TYPENAME_ILLutil_EGlpNum_rselect (
	int *arr,
	int l,
	int r,
	int m,
	EGLPNUM_TYPE * coord,
	ILLrandstate * rstate)
{
	EGLPNUM_TYPE *samplevals = EGLPNUM_TYPENAME_EGlpNumAllocArray (NSAMPLES);
	int i;
	int st, en;
	int n;

	arr += l;
	n = r - l + 1;
	m -= l;

	while (n > SORTSIZE)
	{
		for (i = 0; i < NSAMPLES; i++)
		{
			EGLPNUM_TYPENAME_EGlpNumCopy (samplevals[i], coord[arr[ILLutil_lprand (rstate) % n]]);
		}
		select_EGlpNum_sort_dsample (samplevals, NSAMPLES);
		select_EGlpNum_split (arr, n, &(samplevals[(NSAMPLES - 1) / 2]),
													&st, &en, coord);
		if (st > m)
		{
			n = st;
		}
		else if (en <= m)
		{
			arr += en;
			n -= en;
			m -= en;
		}
		else
		{
			return;
		}
	}

	select_EGlpNum_sort (arr, n, coord);
	EGLPNUM_TYPENAME_EGlpNumFreeArray (samplevals);
	return;
}

static void select_EGlpNum_split (
	int *arr,
	int n,
	EGLPNUM_TYPE * v,
	int *start,
	int *end,
	EGLPNUM_TYPE * coord)
{
	int i, j, k;
	int t;

	i = 0;
	j = k = n;

	while (i < j)
	{
		if (EGLPNUM_TYPENAME_EGlpNumIsLess (coord[arr[i]], *v))
		{
			i++;
		}
		else if (EGLPNUM_TYPENAME_EGlpNumIsEqqual (coord[arr[i]], *v))
		{
			j--;
			ILL_SWAP (arr[i], arr[j], t);
		}
		else
		{
			j--;
			k--;
			t = arr[i];
			arr[i] = arr[j];
			arr[j] = arr[k];
			arr[k] = t;
		}
	}
	*start = j;
	*end = k;
	return;
}

static void select_EGlpNum_sort (
	int *arr,
	int n,
	EGLPNUM_TYPE * coord)
{
	int i, j;
	int t;

	for (i = 1; i < n; i++)
	{
		t = arr[i];
		for (j = i; j > 0 && EGLPNUM_TYPENAME_EGlpNumIsLess (coord[t], coord[arr[j - 1]]); j--)
		{
			arr[j] = arr[j - 1];
		}
		arr[j] = t;
	}
}

static void select_EGlpNum_sort_dsample (
	EGLPNUM_TYPE * samp,
	int n)
{
	int i, j;
	EGLPNUM_TYPE t;

	EGLPNUM_TYPENAME_EGlpNumInitVar (t);

	for (i = 1; i < n; i++)
	{
		EGLPNUM_TYPENAME_EGlpNumCopy (t, samp[i]);
		for (j = i; j > 0 && EGLPNUM_TYPENAME_EGlpNumIsLess (t, samp[j - 1]); j--)
		{
			EGLPNUM_TYPENAME_EGlpNumCopy (samp[j], samp[j - 1]);
		}
		EGLPNUM_TYPENAME_EGlpNumCopy (samp[j], t);
	}
	EGLPNUM_TYPENAME_EGlpNumClearVar (t);
}




