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

/* RCS_INFO = "$RCSfile: read_mps_state.c,v $ $Revision: 1.2 $ $Date: 2003/11/05 16:49:52 $"; */

/****************************************************************************/
/*                                                                          */
/*               Routines to Support  Reading MPS Files                     */
/*                                                                          */
/****************************************************************************/

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <stdlib.h>
#include <stdarg.h>

#include "qs_config.h"
#include "logging-private.h"

#include "eg_lpnum.h"
#include "eg_io.h"
#include "except.h"

#include "mps_EGLPNUM_TYPENAME.h"
#include "rawlp_EGLPNUM_TYPENAME.h"
#include "read_mps_EGLPNUM_TYPENAME.h"
#include "read_lp_EGLPNUM_TYPENAME.h"						/* for EGLPNUM_TYPENAME_ILLget_value */
#include "format_EGLPNUM_TYPENAME.h"


static int TRACE = 0;

#define END_LINE(p)  (((*(p)) == '$' || (*(p)) == '\n' || (*(p)) == '\0') ? 1 : 0)

static int mps_skip_comment (
	EGLPNUM_TYPENAME_ILLread_mps_state * state);
static char ILLmps_next_field_is_number (
	EGLPNUM_TYPENAME_ILLread_mps_state * state);

int EGLPNUM_TYPENAME_ILLmps_state_init (
	EGLPNUM_TYPENAME_ILLread_mps_state * state,
	EGLPNUM_TYPENAME_qsline_reader * file,
	const char *fname)
{
	int i, rval = 0;

	ILL_FAILtrue (file == 0, "need file");
	state->p = 0;
	state->file_name = fname;
	state->file = file;

	for (i = 0; i < ILL_MPS_N_SECTIONS; i++)
	{
		state->section[i] = 0;
	}
	state->active = ILL_MPS_NONE;
	state->intvar = 0;
	state->sosvar = 0;
	state->line_num = 0;
	state->p = 0;

	state->obj = 0;
	state->line[0] = '\0';
	state->key[0] = '\0';
	state->field[0] = '\0';

CLEANUP:
	ILL_RESULT (rval, "EGLPNUM_TYPENAME_ILLmps_state_init");
}

int EGLPNUM_TYPENAME_ILLmps_next_line (
	EGLPNUM_TYPENAME_ILLread_mps_state * state)
{
	int rval = 0;

	/* if field 3 or 5 start with $ rest of line is interpreted as comment */
	state->line[0] = '\0';
	state->p = 0;
	while (EGLPNUM_TYPENAME_ILLline_reader_get (state->line, ILL_namebufsize - 2, state->file)
				 != 0)
	{
		state->line_num++;
		state->key[0] = '\0';
		state->field[0] = '\0';
		state->field_num = 1;
		state->p = state->line;
		if (!EGLPNUM_TYPENAME_ILL_ISBLANK ((state->line)))
		{
			if (state->line[0] == '*' || state->line[0] == '\n')
			{
				continue;								/* comment or blank line */
			}
			else
			{
				if (sscanf (state->p, "%s", state->key) == 1)
				{
					state->p += strlen (state->key);
					while (EGLPNUM_TYPENAME_ILL_ISBLANK (state->p))
					{
						state->p++;
					}
					if (sscanf (state->p, "%s", state->field) == 1)
					{
						state->p += strlen (state->field);
					}
					else
					{
						ILL_FAILfalse (state->field[0] == '\0', "sscanf problem?");
					}
				}
				else
				{
					ILL_FAILfalse (0, "should almost never happen");
				}
			}
		}
		else
		{
			while (EGLPNUM_TYPENAME_ILL_ISBLANK (state->p))
			{
				state->p++;
			}
			if (sscanf (state->p, "%s", state->field) < 1)
			{
				continue;								/* nothing more on line */
			}
			else
			{
				if (state->field[0] == '\0')
				{
					continue;							/* found empty string */
				}
				state->p += strlen (state->field);
			}
		}
		return 0;
	}
CLEANUP:
	return 1;											/* end of file */
}

/* fields 3,5,7,... may start with '$' signifying comment 
 * ==> if we find a '$' as next non blank and 
 *     we read 2,4,6,... fields successfully so far
 *     we have a comment 
 */
static int mps_skip_comment (
	EGLPNUM_TYPENAME_ILLread_mps_state * state)
{
	int rval;

	while (EGLPNUM_TYPENAME_ILL_ISBLANK (state->p))
	{
		state->p++;
	}
	rval = ((*state->p == '$') && (state->field_num >= 2) &&
					(state->field_num % 2 == 0));
	return rval;
}

int EGLPNUM_TYPENAME_ILLmps_next_field (
	EGLPNUM_TYPENAME_ILLread_mps_state * state)
{
	state->field[0] = '\0';
	if (!mps_skip_comment (state))
	{
		if (sscanf (state->p, "%s", state->field) == 1)
		{
			state->p += strlen (state->field) + 1;
			state->field_num++;
			return 0;
		}
	}
	return 1;											/* no more fields */
}

static char get_double (
	EGLPNUM_TYPENAME_ILLread_mps_state * state,
	int peek,
	EGLPNUM_TYPE * coef)
{
	char ok = 0;
	int len, rval = 0;

	ILL_FAILfalse (state != 0, "must have state");
	if (mps_skip_comment (state))
		return 0;
	len = EGLPNUM_TYPENAME_ILLget_value (state->p, coef);
	if (len > 0)
	{
		if (!peek)
		{
			state->p += len;
			state->field_num++;
		}
		ok = 1;
	}
CLEANUP:
	ILL_RESULT (ok, "get_double");
}

int EGLPNUM_TYPENAME_ILLmps_next_coef (
	EGLPNUM_TYPENAME_ILLread_mps_state * state,
	EGLPNUM_TYPE * coef)
{
	int len = 0;

	if (!mps_skip_comment (state))
	{
		len = get_double (state, 0, coef);
	}
	ILL_RESULT (!(len > 0), "EGLPNUM_TYPENAME_ILLmps_next_coef");
}

int EGLPNUM_TYPENAME_ILLmps_next_bound (
	EGLPNUM_TYPENAME_ILLread_mps_state * state,
	EGLPNUM_TYPE * coef)
{
	int len = 0, sign = 1;
	char c, *p;

	if (!mps_skip_comment (state))
	{
		c = *state->p;
		if (c == '-')
		{
			sign = -1;
			len = 1;
		}
		else
		{
			if (c == '+')
			{
				len = 1;
			}
		}
		if (!strncasecmp (state->p + len, "INFINITY", (size_t) 8))
		{
			len += 8;
		}
		else
		{
			if (!strncasecmp (state->p + len, "INF", (size_t) 3))
			{
				len += 3;
			}
		}
		if (len > 1)
		{
			state->p += len;
			p = state->p;
			mps_skip_comment (state);
			if (!END_LINE (state->p) && p == state->p)
			{
				/* found no blanks so this INF/INFINITY is the prefix 
				 * of something else */
				state->p -= len;
				return 1;								/* no coef found */
			}
			else
			{
				if (sign == 1)
					EGLPNUM_TYPENAME_EGlpNumCopy (*coef, EGLPNUM_TYPENAME_ILL_MAXDOUBLE);
				else
					EGLPNUM_TYPENAME_EGlpNumCopy (*coef, EGLPNUM_TYPENAME_ILL_MINDOUBLE);
				state->field_num++;
				ILL_RESULT (0, "EGLPNUM_TYPENAME_ILLmps_next_bound");
			}
		}
		if (get_double (state, 0, coef))
		{
			ILL_RESULT (0, "EGLPNUM_TYPENAME_ILLmps_next_bound");
		}
		else
		{
			ILL_RESULT (1, "EGLPNUM_TYPENAME_ILLmps_next_bound");	/* no coef found */
		}
	}
	ILL_RETURN (1, "EGLPNUM_TYPENAME_ILLmps_next_bound");
}

static char ILLmps_next_field_is_number (
	EGLPNUM_TYPENAME_ILLread_mps_state * state)
{
	EGLPNUM_TYPE d;
	int len = 0;

	if (!mps_skip_comment (state))
	{
		EGLPNUM_TYPENAME_EGlpNumInitVar (d);
		len = get_double (state, 1, &d);
		EGLPNUM_TYPENAME_EGlpNumClearVar (d);
	}
	return (len > 0);
}

void EGLPNUM_TYPENAME_ILLmps_check_end_of_line (
	EGLPNUM_TYPENAME_ILLread_mps_state * state)
{
	if (!mps_skip_comment (state))
	{
		if (!END_LINE (state->p))
		{
			EGLPNUM_TYPENAME_ILLmps_warn (state, "Extra fields on line.");
		}
	}
}

void EGLPNUM_TYPENAME_ILLmps_set_end_of_line (
	EGLPNUM_TYPENAME_ILLread_mps_state * state)
{
	*state->p = '\n';
}

int EGLPNUM_TYPENAME_ILLmps_set_section (
	EGLPNUM_TYPENAME_ILLread_mps_state * state,
	const ILLmps_section sec)
{
	int rval = 0;

	ILL_FAILfalse (sec != ILL_MPS_NONE, "must be in a proper section");
	if (state->section[sec])
	{
		rval = EGLPNUM_TYPENAME_ILLmps_error (state, "Two %s sections.\n", EGLPNUM_TYPENAME_ILLmps_section_name[sec]);
	}
	state->section[sec]++;
	state->active = sec;
CLEANUP:
	ILL_RESULT (rval, "EGLPNUM_TYPENAME_ILLmps_set_section");
}

int EGLPNUM_TYPENAME_ILLmps_int_sos_mode (
	EGLPNUM_TYPENAME_ILLread_mps_state * state)
{
	if (!strcmp (state->field, "'INTORG'"))
	{
		if (state->intvar)
		{
			return !EGLPNUM_TYPENAME_ILLmps_error (state, "'INTEND' expected.\n");
		}
		else
		{
			state->intvar = 1;
			ILL_RESULT (0, "EGLPNUM_TYPENAME_ILLmps_int_sos_mode");
		}
	}
	if (!strcmp (state->field, "'INTEND'"))
	{
		if (state->intvar)
		{
			state->intvar = 0;
			ILL_RESULT (0, "EGLPNUM_TYPENAME_ILLmps_int_sos_mode");
		}
		else
		{
			return !EGLPNUM_TYPENAME_ILLmps_error (state, "'INTORG' expected.\n");
		}
	}
	if (!strcmp (state->field, "'SOSORG'"))
	{
		if (state->sosvar)
		{
			return !EGLPNUM_TYPENAME_ILLmps_error (state, "'SOSEND' expected.\n");
		}
		else
		{
			state->sosvar = 1;
			ILL_RESULT (0, "EGLPNUM_TYPENAME_ILLmps_int_sos_mode");
		}
	}
	if (!strcmp (state->field, "'SOSEND'"))
	{
		if (state->sosvar)
		{
			state->sosvar = 0;
			ILL_RESULT (0, "EGLPNUM_TYPENAME_ILLmps_int_sos_mode");
		}
		else
		{
			return !EGLPNUM_TYPENAME_ILLmps_error (state, "'SOSORG' expected.\n");
		}
	}
	return EGLPNUM_TYPENAME_ILLmps_error (state, "%s is not a MARKER field.\n", state->field);
}

const char *EGLPNUM_TYPENAME_ILLmps_possibly_blank_name (
	const char *field,
	EGLPNUM_TYPENAME_ILLread_mps_state * state,
	ILLsymboltab * tab)
{
	int ind;

	if (ILLsymboltab_lookup (tab, field, &ind) == 0)
	{
		/* Possibly a blank identifier on the line */
		if (ILLmps_next_field_is_number (state))
		{
			/* assume a blank bound ident */
			return " ";
		}
		else
		{
			return field;
		}
	}
	else
	{
		return field;
	}
}

int EGLPNUM_TYPENAME_ILLmps_empty_key (
	EGLPNUM_TYPENAME_ILLread_mps_state * state)
{
	return state->key[0] == '\0';
}

int EGLPNUM_TYPENAME_ILLmps_empty_field (
	EGLPNUM_TYPENAME_ILLread_mps_state * state)
{
	return state->field[0] == '\0';
}

static void mps_err (
	EGLPNUM_TYPENAME_ILLread_mps_state * state,
	int isError,
	const char *format,
	va_list args)
{
	int rval = 0;
	const char *type = (isError) ? "MPS Error" : "MPS Warning";
	int errtype, slen, at;
	EGLPNUM_TYPENAME_qsformat_error error;
	char error_desc[256];

	ILL_FAILfalse_no_rval (format != 0, "format != 0");
	ILL_FAILfalse_no_rval (format[0] != '\0', "format[0] != '0'");
	ILL_FAILfalse_no_rval (state != 0, "state != 0");
	ILL_FAILfalse_no_rval (state->file != 0, "state->file != 0");

	if (state->p == 0)
	{
		at = -1;
	}
	else
	{
		ILL_FAILfalse (state->p >= state->line, "state->p >= state->line");
		at = state->p - state->line;
	}
	vsprintf (error_desc, format, args);
	slen = strlen (error_desc);
	if ((slen > 0) && error_desc[slen - 1] != '\n')
	{
		error_desc[slen] = '\n';
		error_desc[slen + 1] = '\0';
	}

	if (state->file->error_collector != 0)
	{
		errtype = (isError) ? QS_MPS_FORMAT_ERROR : QS_MPS_FORMAT_WARN;
		EGLPNUM_TYPENAME_ILLformat_error_create (&error, errtype, error_desc,
														(int) (state->line_num), state->line, at);
		EGLPNUM_TYPENAME_ILLformat_error (state->file->error_collector, &error);
		EGLPNUM_TYPENAME_ILLformat_error_delete (&error);
	}
	else
	{
		QSlog("%s %d: %s\t", state->file_name, state->line_num,
								state->line);
		QSlog("%s: %s", type, error_desc);
	}
CLEANUP:
	;
}

int EGLPNUM_TYPENAME_ILLmps_error (
	EGLPNUM_TYPENAME_ILLread_mps_state * state,
	const char *format,
	...)
{
	va_list args;

	va_start (args, format);
	mps_err (state, TRUE, format, args);
	/* ILL_RESULT(1, "EGLPNUM_TYPENAME_ILLmps_error"); */
	return 1;
}

void EGLPNUM_TYPENAME_ILLmps_warn (
	EGLPNUM_TYPENAME_ILLread_mps_state * state,
	const char *format,
	...)
{
	va_list args;

	va_start (args, format);
	if (format != 0)
	{
		mps_err (state, FALSE, format, args);
	}
}
