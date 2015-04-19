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

/* RCSINFO $Id: names.c,v 1.2 2003/11/05 16:47:22 meven Exp $ */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "names.h"

#include "logging-private.h"

#include "util.h"
#include "except.h"
#include "symtab.h"


void ILLfree_names (
	char **names,
	int count)
{
	int i;

	if (names)
	{
		for (i = 0; i < count; i++)
		{
			ILL_IFFREE (names[i], char);
		}
		ILL_IFFREE (names, char *);
	}
}

int ILLgenerate_names (
	char prefix,
	int nnames,
	char ***names)
{
	int rval = 0;
	int i, len;
	char *buf = (char *) NULL;

	*names = (char **) NULL;
	if (nnames == 0)
		goto CLEANUP;

	ILL_SAFE_MALLOC (buf, ILL_namebufsize, char);
	ILL_SAFE_MALLOC (*names, nnames, char *);

	for (i = 0; i < nnames; i++)
	{
		(*names)[i] = (char *) NULL;
	}

	for (i = 0; i < nnames; i++)
	{
		sprintf (buf, "%c%d", prefix, i);
		len = strlen (buf) + 1;
		ILL_SAFE_MALLOC ((*names)[i], len, char);

		strcpy ((*names)[i], buf);
	}

CLEANUP:

	if (rval)
	{
		if (*names)
		{
			ILLfree_names (*names, nnames);
			*names = (char **) NULL;
		}
	}

	ILL_IFFREE (buf, char);

	if (rval)
	{
		QSlog("ILLsymboltab_generate_names failed");
	}
	return rval;
}
