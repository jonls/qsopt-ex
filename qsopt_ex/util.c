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

/* RCSINFO $Id: util.c,v 1.2 2003/11/05 16:47:22 meven Exp $ */
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
/*               MISCELLANEOUS UTILITY ROUTINES                             */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: October 12, 1995                                                  */
/*  Date: September 28, 1997                                                */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  unsigned int ILLutil_nextprime (unsigned int x)                         */
/*    FINDS the smallest positive prime >= x                                */
/*                                                                          */
/*  int ILLutil_our_gcd (int a, int b)                                      */
/*    COMPUTES gcd(a,b)                                                     */
/*    -gcd(a,b) is always >= 0                                              */
/*    -a and b can be negative, positive, or zero                           */
/*                                                                          */
/*  int ILLutil_our_lcm (int a, int b)                                      */
/*    COMPUTES lcm(a,b)                                                     */
/*    -lcm(a,b) is always >= 0                                              */
/*    -a and b can be negative, positive, or zero                           */
/*                                                                          */
/*  double ILLutil_our_floor (double x)                                     */
/*    REURNS the greatest integer no larger than x.                         */
/*                                                                          */
/*  double ILLutil_our_ceil (double x)                                      */
/*    REURNS the least integer no smaller than x.                           */
/*                                                                          */
/*  double ILLutil_our_frac (double x)                                      */
/*    REURNS the fractional part of x.                                      */
/*                                                                          */
/*  char *ILLutil_strchr (const char *s, int c)                             */
/*    RETURNS a pointer to the first occurrence of c in s, or NULL if c     */
/*    does not occur in s                                                   */
/*                                                                          */
/* int ILLutil_strcasecmp(const char *s1, const char *s2)                   */
/*    RETURNS the string comparison, iqnoring the case of the letters.      */
/*                                                                          */
/* int ILLutil_strncasecmp(const char *s1, const char *s2, size_t n)        */
/*    RETURNS the string comparision, iqnoring case, looks at max n bytes.  */
/*                                                                          */
/*  char *ILLutil_str(const char *s)                                        */
/*    allocates and returns a copy of s                                     */
/*                                                                          */
/****************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "except.h"
#include "util.h"


static int isprime (
	unsigned int x);


unsigned int ILLutil_nextprime (
	unsigned int x)
{
	if (x < 3)
		return 3;
	x |= 1;
	while (!isprime (x))
		x += 2;
	return x;
}

static int isprime (
	unsigned int p)
{
	unsigned int i;

	if ((p & 1) == 0)
		return 0;
	for (i = 3; i * i <= p; i += 2)
	{
		if (p % i == 0)
			return 0;
	}
	return 1;
}

int ILLutil_our_gcd (
	int a,
	int b)
{
	int c;

	if (a < 0)
		a = -a;
	if (b < 0)
		b = -b;
	if (a > b)
		ILL_SWAP (a, b, c);

	while (a)
	{
		c = b % a;
		b = a;
		a = c;
	}
	return b;
}

int ILLutil_our_lcm (
	int a,
	int b)
{
	int c;

	if (a < 0)
		a = -a;
	if (b < 0)
		b = -b;

	c = ILLutil_our_gcd (a, b);

	return (a / c) * b;
}

double ILLutil_our_floor (
	double x)
{
	return floor (x);
}

double ILLutil_our_ceil (
	double x)
{
	return ceil (x);
}

double ILLutil_our_frac (
	double x)
{
	return x - floor (x);
}

double ILLutil_norm_sqr (
	double *v,
	int len)
{
	int i;
	double sum = 0.0;

	for (i = 0; i < len; i++)
		sum += v[i] * v[i];
	return sum;
}

int ILLutil_our_log2 (
	int a)
{
	int i = 0, j = 1;

	while (j < a)
	{
		j = j << 1;
		i++;
	}
	return i ? i : 1;
}

const char *ILLutil_strchr (
	const char *s,
	int c)
{
	while (*s)
	{
		if (*s == c)
			return s;
		s++;
	}
	return (char *) NULL;
}

int ILLutil_strcasecmp (
	const char *s1,
	const char *s2)
{
	return strcasecmp (s1, s2);
}

int ILLutil_strncasecmp (
	const char *s1,
	const char *s2,
	size_t n)
{
	return strncasecmp (s1, s2, n);
}

int ILLutil_index (
	const char *list[],
	const char *name)
{
	int i;

	for (i = 0; list[i] != NULL; i++)
	{
		if (!strcmp (name, list[i]))
		{
			return i;
		}
	}
	return -1;
}

int ILLutil_array_index (
	char *list[],
	int n,
	const char *name)
{
	int i;

	for (i = 0; i < n; i++)
	{
		if ((list[i] != NULL) && !strcmp (name, list[i]))
		{
			return i;
		}
	}
	return -1;
}

char *ILLutil_str (
	const char *str)
{
	int len;
	char *cpy = NULL;

	if (str != NULL)
	{
		len = strlen (str) + 1;
		ILL_SAFE_MALLOC_no_rval (cpy, len, char);

		strcpy (cpy, str);
	}
CLEANUP:
	return cpy;
}
