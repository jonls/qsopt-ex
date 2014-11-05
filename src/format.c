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

/* RCS_INFO = "$RCSfile: format_error.c,v $ $Revision: 1.2 $ $Date: 2003/11/05 16:49:52 $"; */

#include "qs_config.h"

#include "eg_fp.h"
#include "eg_io.h"

#include "format.h"
#include "qsopt.h"
#include "iqsutil.h"

int ILLformat_error_create (
	qsformat_error * error,
	int mode,
	const char *desc,
	int lineNum,
	const char *theLine,
	int atPos)
{
	int len;
	int rval = 0;

	error->theLine = NULL;
	error->desc = NULL;
	error->next = NULL;

	ILL_FAILtrue (desc == NULL, "non empty error desc please");
	ILL_FAILtrue (mode >= QS_INPUT_NERROR
								|| mode < 0, "0<= mode <=QS_INPUT_NERROR");
	error->type = mode;
	len = strlen (desc);
	ILL_SAFE_MALLOC (error->desc, len + 1, char);

	strcpy (error->desc, desc);
	error->lineNumber = lineNum;
	if (theLine != NULL)
	{
		len = strlen (theLine);
		ILL_SAFE_MALLOC (error->theLine, len + 2, char);

		strcpy (error->theLine, theLine);
		if (error->theLine[len - 1] != '\n')
		{
			error->theLine[len] = '\n';
			error->theLine[len + 1] = '\0';
		}
	}
	error->at = atPos;
CLEANUP:
	if (rval)
	{
		ILLformat_error_delete (error);
	}
	return rval;
}

void ILLformat_error_delete (
	qsformat_error * error)
{
	ILL_IFFREE (error->desc, char);
	ILL_IFFREE (error->theLine, char);
}

void ILLformat_error_print (
	EGioFile_t * out,
	qsformat_error * error)
{
	int at = error->at;
	int tp = error->type;
	const char *type = "Error";
	const char *line = NULL;
	int i;

	type = QSformat_error_type_string (tp);

	EGioPrintf (out, "%s  line %d pos %d\n",
					 type, QSerror_get_line_number (error), at);
	line = QSerror_get_line (error);
	if (line != NULL)
	{
		EGioPrintf (out, "LINE %s", line);
		if (at >= 0)
		{
			EGioPrintf (out, ".....");
			for (i = 0; i <= (at - 1); i++)
			{
				if (line[i] == '\t')
				{
					EGioPrintf(out, "\t");
				}
				else
				{
					EGioPrintf(out, ".");
				}
			}
			EGioPrintf (out, "^\n");
		}
	}
	else
	{
		EGioPrintf (out, "NO LINE\n");
	}
	EGioPrintf (out, "MSG: %s\n", QSerror_get_desc (error));
}


qserror_collector *ILLerror_collector_new (
	qsadd_error_fct fct,
	void *dest)
{
	int rval = 0;
	qserror_collector *c = NULL;

	ILL_SAFE_MALLOC (c, 1, qserror_collector);
	c->add_error = fct;
	c->dest = dest;

CLEANUP:
	if (rval)
	{
		ILL_IFFREE (c, qserror_collector);
	}
	return c;
}

qserror_collector *ILLerror_memory_collector_new (
	qserror_memory * dest)
{
	return ILLerror_collector_new (ILLadd_error_to_memory, dest);
}

void ILLerror_collector_free (
	qserror_collector * c)
{
	ILL_IFFREE (c, qserror_collector);
}

qserror_memory *ILLerror_memory_create (
	int takeErrorLines)
{
	int rval = 0, i;
	qserror_memory *mem = NULL;

	ILL_SAFE_MALLOC (mem, 1, qserror_memory);
	for (i = 0; i < QS_INPUT_NERROR; i++)
	{
		mem->has_error[i] = 0;
	}
	mem->error_list = NULL;
	mem->nerror = 0;
	mem->hasErrorLines = takeErrorLines;
CLEANUP:
	return mem;
}

void ILLerror_memory_free (
	qserror_memory * mem)
{
	qsformat_error *ths, *nxt;

	if (mem != NULL)
	{
		ths = mem->error_list;
		while (ths != NULL)
		{
			nxt = ths->next;
			ILL_IFFREE (ths, qsformat_error);
			ths = nxt;
		}
		ILL_IFFREE (mem, qserror_memory);
	}
}

int ILLadd_error_to_memory (
	void *dest,
	const qsformat_error * error)
{
	int rval = 0;
	qserror_memory *mem = (qserror_memory *) dest;
	qsformat_error *e = 0;

	ILL_CHECKnull (mem, "must give non NULL qserror_memory");

	ILL_SAFE_MALLOC (e, 1, qsformat_error);
	rval = ILLformat_error_create (e, error->type, error->desc,
																 error->lineNumber,
																 (mem->hasErrorLines) ? error->theLine : NULL,
																 error->at);
	ILL_CLEANUP_IF (rval);
	e->next = mem->error_list;
	mem->error_list = e;
	mem->nerror++;
	mem->has_error[error->type]++;

CLEANUP:
	if (rval)
	{
		ILLformat_error_delete (e);
		ILL_IFFREE (e, qsformat_error);
	}
	return rval;
}
