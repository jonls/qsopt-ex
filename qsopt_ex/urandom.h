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

#ifndef __URANDOM_H__
#define __URANDOM_H__
/****************************************************************************/
/*                                                                          */
/*                             urandom.c                                    */
/*                                                                          */
/****************************************************************************/

/* since urandom's generator does everything modulo ILL_PRANDMAX, if two
 * seeds are congruent mod x and x|ILL_PRANDMAX, then the resulting numbers
 * will be congruent mod x.  One example was if ILL_PRANDMAX = 1000000000 and
 * urandom is used to generate a point set from a 1000x1000 grid, seeds
 * congruent mod 1000 generate the same point set.
 *
 * For this reason, we use 1000000007 (a prime)
 */
#define ILL_PRANDMAX 1000000007

typedef struct ILLrandstate
{
	int a;
	int b;
	int arr[55];
}
ILLrandstate;


void ILLutil_sprand (
	int seed,
	ILLrandstate * r);

int ILLutil_lprand (
	ILLrandstate * r);



#endif
