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

/* RCSINFO $Id: urandom.c,v 1.2 2003/11/05 16:47:22 meven Exp $ */
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
/*              MACHINE INDEPENDENT RANDOM NUMBER GENERATOR                 */
/*                                                                          */
/*                             TSP CODE                                     */
/*                                                                          */
/*                                                                          */
/*  Written by:  DIMACS  (modified for TSP)                                 */
/*  Date: February 7, 1995  (cofeb16)                                       */
/*                                                                          */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  void ILLutil_sprand (int seed, ILLrandstate *r)                         */
/*    - Call once to initialize the generator.                              */
/*                                                                          */
/*  int ILLutil_lprand (QSrandstate *r)                                     */
/*    - Returns an integer in the range 0 to ILL_PRANDMAX - 1.              */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/*    NOTES (from DIMACS):                                                  */
/*        This file contains a set of c-language functions for generating   */
/*    uniform integers.   This is a COMPLETELY PORTABLE generator. It will  */
/*    give IDENTICAL sequences of random numbers for any architecture with  */
/*    at least 30-bit integers, regardless of the integer representation,   */
/*    INT_MAX value, or roundoff/truncation method, etc.                    */
/*        This Truly Remarkable RNG is described more fully in              */
/*    J. Bentley's column, ``The Software Exploratorium ''. It is based on  */
/*    one in Knuth, Vol 2, Section 3.2.2 (Algorithm A).                     */
/*                                                                          */
/****************************************************************************/


#include "util.h"


void ILLutil_sprand (
	int seed,
	ILLrandstate * r)
{
	int i, ii;
	int last, next;
	int *arr = r->arr;

	arr[0] = last = seed;
	next = 1;
	for (i = 1; i < 55; i++)
	{
		ii = (21 * i) % 55;
		arr[ii] = next;
		next = last - next;
		if (next < 0)
			next += ILL_PRANDMAX;
		last = arr[ii];
	}
	r->a = 0;
	r->b = 24;
	for (i = 0; i < 165; i++)
		last = ILLutil_lprand (r);
}


int ILLutil_lprand (
	ILLrandstate * r)
{
	int t;

	if (r->a-- == 0)
		r->a = 54;
	if (r->b-- == 0)
		r->b = 54;

	t = r->arr[r->a] - r->arr[r->b];

	if (t < 0)
		t += ILL_PRANDMAX;

	r->arr[r->a] = t;

	return t;
}


#ifdef      TRY_CODE

/*-----------------------------------------------*/
/* This is a little driver program so you can    */
/* test the code.                                */
/* Typing: a.out 0 3 1                           */
/* should produce                                */
/*     921674862                                 */
/*     250065336                                 */
/*     377506581                                 */
/*  Typing: a.out 1000000 1 2                    */
/*  should produce                               */
/*     57265995                                  */
/*-----------------------------------------------*/

int main (
	int ac,
	char **av)
{
	int i;
	int j;
	int n;
	int m;
	int seed;
	ILLrandstate rstate;

	if (ac < 4)
	{
		fprintf (stderr, "Usage: #discard #print #seed\n");
		return 0;
	}
	m = atoi (av[1]);							/* Number to discard initially */
	n = atoi (av[2]);							/* Number to print */
	seed = atoi (av[3]);					/* Seed */

	ILLutil_sprand (seed, &rstate);

	for (i = 0; i < m; i++)
		j = ILLutil_lprand (&rstate);
	for (i = 0; i < n; i++)
		printf ("%ld\n", ILLutil_lprand (&rstate));
	return 0;
}

#endif /* TRY_CODE */
