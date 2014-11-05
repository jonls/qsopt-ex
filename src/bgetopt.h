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

#ifndef __BGETOPT_H__
#define __BGETOPT_H__

/****************************************************************************/
/*                                                                          */
/*                             bgetopt.c                                    */
/*                                                                          */
/****************************************************************************/
int ILLutil_bix_getopt (
	int argc,
	char **argv,
	const char *def,
	int *p_optind,
	char **p_optarg);


#define ILL_BIX_GETOPT_UNKNOWN -3038



#endif
