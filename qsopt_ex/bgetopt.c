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

/* RCSINFO $Id: bgetopt.c,v 1.2 2003/11/05 16:47:22 meven Exp $ */
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
/*                        PORTABLE GETOPT                                   */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: 1993 (?) (fmfeb02)                                                */
/*  Modified: 15 February 1995 (bico)  - added warning                      */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int ILLutil_bix_getopt (int argc, char **argv, const char *def,          */
/*      int *p_optind, char **p_optarg)                                     */
/*     parse an argument list                                               */
/*                                                                          */
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "bgetopt.h"
#include "logging-private.h"


int ILLutil_bix_getopt (
	int ac,
	char **av,
	const char *def,
	int *p_optind,
	char **p_optarg)
{
	int c;
	char *sp = av[*p_optind];
	char bwarn[2];

	if (*p_optind < 1 || *p_optind >= ac)
	{
		*p_optind = ac;
		return (EOF);
	}
	if ((int) *sp != (int) '-')
		return (EOF);
	if ((int) *(sp + 1) == (int) '-')
	{
		(*p_optind)++;
		return (EOF);
	}
	(av[*p_optind])++;
	sp++;
	while ((int) *sp != (int) *def && (int) *def != (int) '\0')
		def++;
	if ((int) *def == (int) '\0')
	{
		*p_optind = ac;
		bwarn[0] = *sp;							/* Bico: February 8, 1995 */
		bwarn[1] = '\0';
		QSlog("Illegal option: -%s", bwarn);
		return ILL_BIX_GETOPT_UNKNOWN;
	}
	if ((int) *(def + 1) != (int) ':')
	{
		c = *sp;
		if ((int) *(sp + 1) != (int) '\0')
			*sp = '-';
		else
			(*p_optind)++;
		return (c);
	}
	else
	{
		if ((int) *(sp + 1) != (int) '\0')
		{
			*p_optarg = sp + 1;
			c = *sp;
			(*p_optind)++;
			return (c);
		}
		else if (*p_optind >= ac - 1)
		{
			*p_optind = ac;
			return (EOF);
		}
		else
		{
			*p_optarg = av[*p_optind + 1];
			c = *sp;
			*p_optind += 2;
			return (c);
		}
	}
}
