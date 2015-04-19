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

/* RCSINFO $Id: exception.c,v 1.2 2003/11/05 16:47:22 meven Exp $ */
#include <stdio.h>
#include <string.h>

#include "except.h"

#include "logging-private.h"


void ILL_report (
	const char *msg,
	const char *fct,
	const char *file,
	unsigned int line,
	int with_src_info)
{
	if (msg != NULL)
	{
		QSlog("FAILURE: %s", msg);
		if (with_src_info == 1)
		{
			if (fct != NULL)
			{
				QSlog("\tin function %s", fct);
			}
			QSlog("\tin file %s line %d", file, line);
		}
	}
}

/* by default we turn off verbose messages of singular basis */
int __QS_SB_VERB = 1000;
