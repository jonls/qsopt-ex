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

/* RCS_INFO = "$RCSfile: lp.c,v $ $Revision: 1.2 $ $Date: 2003/11/05 16:49:52 $"; */

/****************************************************************************/
/*                                                                          */
/*               Routines for Reading and Writing LP Files                  */
/*                                                                          */
/*  EXPORTED FUNCTIONS                                                      */
/*                                                                          */
/*        int ILLwrite_lp (EGioFile_t *out, ILLlpdata *lp)                  */
/*        int ILLread_lp (qsline_reader f, const char *fname, rawlpdata *lp)*/
/*        int ILLis_lp_name_char (char c, int pos)                          */
/*        int ILLread_constraint_expr(ILLread_lp_state *state,              */
/*          rawlpdata *lp, int rowind, int allowNew)                        */
/*        int ILLread_constraint_name (ILLread_lp_state *state,             */
/*          char **rowname)                                                 */
/*        int ILLread_one_constraint (ILLread_lp_state *state,              */
/*          const char *rowname, rawlpdata *lp, int allowNewCols)           */
/*                                                                          */
/****************************************************************************/

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "qs_config.h"

#include "eg_lpnum.h"
#include "eg_io.h"

/* from the cplex manual: 
 
    not exceed 16 characters, all of which must be alphanumeric
    (a-z, A-Z, 0-9) or one of these symbols: ! " # $ % & ( ) / ,
    . ; ? @ _ ` ' { } | ~. Longer names will be truncated to 16
    characters. A variable name can not begin with a number or a
    period. 

    The letter E or e, alone or followed by other valid symbols,
    or followed by another E or e, should be avoided as this
    notation is reserved for exponential entries. Thus, variables
    can not be named e9, E-24, E8cats, or other names that could
    be interpreted as an exponent. Even variable names such as
    eels or example can cause a read error, depending on their
    placement in an input line. 

    Also, the following characters are not valid in variable
    names (in order to allow for quadratic objective
    information): ^, *, [ and ]. 
*/
/* OUR DEFINTION:
 * -) variables consist of a-z A-Z 0-9 !"#$%(),;.?@_`'{}|~ 
 *    don't start with a digit or '.'
 * return  0 for variable 
 * return -1 for keyword 
 * return  1 otherwise 
 */


int ILLis_lp_name_char (
	int c,
	int pos)
{
	return ((('a' <= c) && (c <= 'z')) ||
					(('A' <= c) && (c <= 'Z')) ||
					((pos > 0) && ('0' <= c) && (c <= '9')) ||
					((pos > 0) && (c == '.')) ||
					(strchr ("!\"#$%&()/,;?@_`'{}|~", c) != NULL));
}


/* 
 * -) anything after '\' is comment 
 * -) Problem is optional and comes first
 * -) Minimize, Maximize, ... comes after Problem
 * -) variables consist of a-z A-Z 0-9!"#$%(),;.?@_`'{}|~ 
 */

/*  NOTES                                                                   */
/*   
 *    don't start with a digit or '.'
 */

static const int LINE_LEN = 256;

#include "iqsutil.h"
#include "lp.h"
#include "rawlp.h"
#include "read_lp.h"
#include "write_lp.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif
//extern EGlpNum_t SZERO_TOLER;
static int TRACE = 0;

static int read_problem_name (
	ILLread_lp_state * state,
	rawlpdata * lp);
static int read_minmax (
	ILLread_lp_state * state,
	rawlpdata * lp);
static int read_objective (
	ILLread_lp_state * state,
	rawlpdata * lp);
static int read_objective (
	ILLread_lp_state * state,
	rawlpdata * lp);
static int read_constraints (
	ILLread_lp_state * state,
	rawlpdata * lp,
	int allowNewCols);
static int read_colname (
	ILLread_lp_state * state,
	ILLsymboltab * coltab,
	int mustHave);
static int read_integer (
	ILLread_lp_state * state,
	rawlpdata * lp);
static int read_bounds (
	ILLread_lp_state * state,
	rawlpdata * lp);
static int add_var (
	rawlpdata * lp,
	ILLread_lp_state * state,
	EGlpNum_t coef,
	int row,
	int allowNew);

/*------------------------------------------------------------------------ 
 * ILLwrite_lp and support routines
 */
static int fix_names (
	qserror_collector * collector,
	char **names,
	int nnames,
	const char *extra,
	int prefix,
	char ***newnames);
static void write_objective (
	ILLlpdata * lp,
	const char *objname,
	char **colnames);
static int write_row (
	ILLlpdata * lp,
	ILLlp_rows * lprows,
	int i,
	char **rownames,
	char **colnames,
	int *colInRow,
	EGlpNum_t * colCoef);
static int write_bounds (
	ILLlpdata * lp,
	char **colnames);
static void write_intvars (
	ILLlpdata * lp,
	char **colnames);

int ILLwrite_lp (
	ILLlpdata * lp,
	qserror_collector * collector)
{
	int rval = 0;
	int i;
	ILLlp_rows lp_rows, *lprows = NULL;
	char **colnames = (char **) NULL;
	char **rownames = (char **) NULL;
	EGlpNum_t *colCoef = NULL;
	int *colInRow = NULL;
	const char *objname;

	ILL_FAILfalse (lp, "called without data\n");
	if (lp->nstruct == 0 || lp->nrows == 0)
	{
		EG_RETURN (rval);
	}
	ILL_FAILfalse (lp->colnames != NULL, "lp->colnames != NULL");
	ILL_FAILfalse (lp->rownames != NULL, "lp->rownames != NULL");
	ILL_FAILfalse (lp->nstruct == lp->coltab.tablesize,
								 "lp coltab has nstruct entries");
	if (lp->objname == (char *) NULL)
	{
		ILL_FAILfalse (lp->nrows == lp->rowtab.tablesize,
									 "lp rowtab should have nrows entries");
	}
	else
	{
		ILL_FAILfalse (lp->nrows + 1 == lp->rowtab.tablesize,
									 "lp rowtab should have nrows+1 entries");
		ILL_FAILfalse (ILLsymboltab_contains (&lp->rowtab, lp->objname),
									 "rowtab must contain objname");
	}

	rval = fix_names (collector, lp->colnames, lp->nstruct, NULL, 'x', &colnames);
	CHECKRVALG (rval, CLEANUP);

	rval = fix_names (collector, lp->rownames, lp->nrows,
										(lp->objname) ? lp->objname : "obj", 'c', &rownames);
	CHECKRVALG (rval, CLEANUP);
	objname = rownames[lp->nrows];

	ILL_FAILtrue (objname == NULL, "OOps, that should never happen");
	CHECKRVALG (rval, CLEANUP);

	if (lp->sos.matcols > 0)
	{
		rval +=
			ILLdata_error (collector, "Can't express SOS information in LP format.");
	}

	write_objective (lp, objname, colnames);

	/* Note, ILLlp_rows_init returns cols ordered by structmap, so we may use 
	 * colnames[i] when pulling i from the matrix data.  */

	lprows = &lp_rows;
	if (ILLlp_rows_init (lprows, lp, 0) != 0)
	{
		rval += 1;
		ILL_FAILtrue (rval, "ILLlp_rows_init failed\n");
	}

	colCoef = EGlpNumAllocArray (lp->nstruct);
	ILL_SAFE_MALLOC (colInRow, lp->nstruct, int);

	for (i = 0; i < lp->nstruct; i++)
	{
		colInRow[i] = -1;
	}

	ILLprint_report (lp, "Subject To\n");
	for (i = 0; i < lp->nrows; i++)
	{
		if (lprows->rowcnt[i] == 0)
		{
			/*
			 * ILLdata_warn (collector, "Not printing  empty row \"%s\".", rownames[i]);
			 */
			continue;
		}
		rval += write_row (lp, lprows, i, rownames, colnames, colInRow, colCoef);
	}

	rval += write_bounds (lp, colnames);

	if (lp->intmarker != NULL)
	{
		write_intvars (lp, colnames);
	}

	ILLprint_report (lp, "End\n");
CLEANUP:
	if (lprows != NULL)
	{
		ILLlp_rows_clear (lprows);
	}
	ILLfree_names (colnames, lp->nstruct);
	ILLfree_names (rownames, lp->nrows + 1);
	EGlpNumFreeArray (colCoef);
	ILL_IFFREE (colInRow, int);

	EG_RETURN (rval);
}

static void write_objective (
	ILLlpdata * lp,
	const char *objname,
	char **colnames)
{
	int ri, i, k, var;
	ILLwrite_lp_state ln, *line = &ln;

	if (lp->probname != NULL)
	{
		ILLprint_report (lp, "Problem\n %s\n", lp->probname);
	}
	if (lp->objsense == ILL_MIN)
	{
		ILLprint_report (lp, "Minimize\n");
	}
	else
	{
		ILLprint_report (lp, "Maximize\n");
	}
	ILLwrite_lp_state_init (line, NULL);
	ILLwrite_lp_state_append (line, " ");
	ILLwrite_lp_state_append (line, objname);
	ILLwrite_lp_state_append (line, ": ");
	ILLwrite_lp_state_save_start (line);

	for (ri = 0, var = 0; ri < lp->nstruct; ri++)
	{
		i = lp->structmap[ri];
		if (EGlpNumIsNeqqZero (lp->obj[i]))
		{
			ILLwrite_lp_state_append_coef (line, lp->obj[i], var);
			ILLwrite_lp_state_append (line, " ");
			ILLwrite_lp_state_append (line, colnames[ri]);
			var++;

			/* we put a least 4 terms on a line 
			 * and then we stop after LINE_LEN or more characters 
			 */
			if ((line->total >= LINE_LEN) && (var >= 4))
			{
				/* see whether there is another term 
				 * if so append a '+' and print line */
				k = ri + 1;
				while (k < lp->nstruct)
				{
					if (EGlpNumIsLessZero (lp->obj[lp->structmap[k]]))
					{
						break;
					}
					else
					{
						if (EGlpNumIsGreatZero (lp->obj[lp->structmap[k]]))
						{
							ILLwrite_lp_state_append (line, " +");
							break;
						}
					}
					k++;
				}
				var = 0;								/* next line does not need to prefix coef with '+' */
				ILLprint_report (lp, "%s\n", line->buf);
				ILLwrite_lp_state_start (line);
			}
		}
	}
	if (var > 0)
	{
		ILLprint_report (lp, "%s\n", line->buf);
	}
}

static void write_the_expr (
	ILLlpdata * lp,
	ILLwrite_lp_state * line,
	char *rowname,
	ILLlp_rows * lprows,
	int row,
	char **colnames,
	int *colInRow,
	EGlpNum_t * colCoef,
	int ncols)
{
	int var, firstVar, k, i;
	EGlpNum_t *coef;

	ILLwrite_lp_state_init (line, NULL);
	if (rowname != NULL)
	{
		ILLwrite_lp_state_append (line, " ");
		ILLwrite_lp_state_append (line, rowname);
		ILLwrite_lp_state_append (line, ": ");
	}
	else
	{
		ILLwrite_lp_state_append (line, "   ");
	}
	ILLwrite_lp_state_save_start (line);

	for (k = lprows->rowbeg[row];
			 k < lprows->rowbeg[row] + lprows->rowcnt[row]; k++)
	{
		i = lprows->rowind[k];
		colInRow[i] = row;
		EGlpNumCopy (colCoef[i], lprows->rowval[k]);
	}
	var = 0;
	firstVar = 1;
	for (i = 0; i < ncols; i++)
	{
		if (colInRow[i] == row)
		{
			if (EGlpNumIsNeqqZero (colCoef[i]))
			{
				coef = &(colCoef[i]);
				if (line->total >= LINE_LEN)
				{
					ILLprint_report (lp, "%s\n", line->buf);
					ILLwrite_lp_state_start (line);
					if ((!firstVar) && !EGlpNumIsLessZero (*coef))
					{
						ILLwrite_lp_state_append (line, " +");
					}
					var = 0;
				}

				ILLwrite_lp_state_append_coef (line, *coef, var);
				ILLwrite_lp_state_append (line, " ");
				ILLwrite_lp_state_append (line, colnames[i]);
				var++;
				firstVar = 0;
			}
		}
	}
}

static int write_row (
	ILLlpdata * lp,
	ILLlp_rows * lprows,
	int i,
	char **rownames,
	char **colnames,
	int *colInRow,
	EGlpNum_t * colCoef)
{
	ILLwrite_lp_state ln, *line = &ln;
	int rval = 0;
	EGlpNum_t ntmp;

	write_the_expr (lp, line, rownames[i], lprows, i, colnames,
									colInRow, colCoef, lp->nstruct);

	switch (lp->sense[i])
	{
	case 'G':
		ILLwrite_lp_state_append (line, " >= ");
		ILLwrite_lp_state_append_number (line, lp->rhs[i]);
		break;
	case 'L':
		ILLwrite_lp_state_append (line, " <= ");
		ILLwrite_lp_state_append_number (line, lp->rhs[i]);
		break;
	case 'E':
		ILLwrite_lp_state_append (line, " = ");
		ILLwrite_lp_state_append_number (line, lp->rhs[i]);
		break;
	case 'R':
		ILL_FAILtrue (!lp->rangeval, "RANGE constraints without values\n");
		EGlpNumInitVar (ntmp);
		ILLwrite_lp_state_append (line, " >= ");
		ILLwrite_lp_state_append_number (line, lp->rhs[i]);

		ILLwrite_lp_state_append (line, " \t\\ RANGE (");
		ILLwrite_lp_state_append_number (line, lp->rhs[i]);
		ILLwrite_lp_state_append (line, ", ");
		EGlpNumCopySum (ntmp, lp->rhs[i], lp->rangeval[i]);
		ILLwrite_lp_state_append_number (line, ntmp);
		ILLwrite_lp_state_append (line, ")");
		ILLprint_report (lp, "%s\n", line->buf);

		write_the_expr (lp, line, NULL, lprows, i,
										colnames, colInRow, colCoef, lp->nstruct);
		ILLwrite_lp_state_append (line, " <= ");
		ILLwrite_lp_state_append_number (line, ntmp);
		EGlpNumClearVar (ntmp);
		break;
	default:
		ILL_FAILtrue (1, "Unknown row sense\n");
	}

	ILLprint_report (lp, "%s\n", line->buf);
CLEANUP:
	EG_RETURN (rval);
}

static int write_bounds (
	ILLlpdata * lp,
	char **colnames)
{
	int ri, i, rval = 0;
	int prtLower, prtUpper;
	ILLwrite_lp_state l, *line = &l;

	ILL_FAILtrue (lp->lower == NULL || lp->upper == NULL,
								"Should not call write_bounds when lower or upper are NULL");
	ri = ILLraw_first_nondefault_bound (lp);
	if (ri != lp->nstruct)
	{
		ILLprint_report (lp, "Bounds\n");
		ILLwrite_lp_state_init (line, " ");
		ILLwrite_lp_state_save_start (line);

		for (; ri < lp->nstruct; ri++)
		{
			ILLwrite_lp_state_start (line);
			i = lp->structmap[ri];
			if (EGlpNumIsEqqual (lp->lower[i], lp->upper[i]))
			{
				ILLwrite_lp_state_append (line, " ");
				ILLwrite_lp_state_append (line, colnames[ri]);
				ILLwrite_lp_state_append (line, " = ");
				ILLwrite_lp_state_append_number (line, lp->upper[i]);
				ILLprint_report (lp, "%s\n", line->buf);
				continue;
			}
			if ((EGlpNumIsEqqual (lp->lower[i], ILL_MINDOUBLE)) &&
					(EGlpNumIsEqqual (lp->upper[i], ILL_MAXDOUBLE)))
			{
				ILLwrite_lp_state_append (line, colnames[ri]);
				ILLwrite_lp_state_append (line, " free");
				ILLprint_report (lp, "%s\n", line->buf);
				continue;
			}
			prtLower = !ILLraw_default_lower (lp, i);
			prtUpper = !ILLraw_default_upper (lp, i, ri);
			if (prtLower || prtUpper)
			{
				if (prtLower)
				{
					ILLwrite_lp_state_append_number (line, lp->lower[i]);
					ILLwrite_lp_state_append (line, " <= ");
				}
				if (prtLower || prtUpper)
				{
					ILLwrite_lp_state_append (line, colnames[ri]);
				}
				if (prtUpper)
				{
					ILLwrite_lp_state_append (line, " <= ");
					ILLwrite_lp_state_append_number (line, lp->upper[i]);
				}
				ILLprint_report (lp, "%s\n", line->buf);
			}
		}
	}
CLEANUP:
	EG_RETURN (rval);
}

static void write_intvars (
	ILLlpdata * lp,
	char **colnames)
{
	ILLwrite_lp_state ln, *line = &ln;
	int var, j;

	ILLprint_report (lp, "Integer\n");
	ILLwrite_lp_state_init (line, " ");
	ILLwrite_lp_state_save_start (line);

	for (j = 0, var = 0; j < lp->nstruct; j++)
	{
		if (lp->intmarker[j])
		{
			if (var > 0)
			{
				ILLwrite_lp_state_append (line, " ");
			}
			ILLwrite_lp_state_append (line, colnames[j]);
			var++;
			if (line->total >= LINE_LEN)
			{
				ILLprint_report (lp, "%s\n", line->buf);
				ILLwrite_lp_state_init (line, " ");
				var = 0;
			}
		}
	}
	if (var > 0)
	{
		ILLprint_report (lp, "%s\n", line->buf);
	}
}

/* ------------------------------------------------------------ */
/* fix up names that are numbers, i.e. prefix with x X x_ or X_ */

/* 
 * redefine names that start with [0-9]; i.e. give a prefix of 
 *         "x" | "X" | "x_" | "X_"  
 *         or if all these are already taken prefix with "X"<number>
 *         make sure names contain solely the characters: 
 *            [a-zA-Z0-9] and ! " # $ % & ( ) / , . ; ? @ _ ` ' { } | ~
 *         rename names with 'bad' chars to <x X x_ X_><number>
 */
static int fix_names (
	qserror_collector * collector,
	char **names,
	int nnames,
	const char *extra,
	int pref,
	char ***newnames)
{
	ILLsymboltab symt, *symtab = NULL;
	int rval = 0, i, j, n, ind, hit;
	char **n_names = NULL;
	const char *old_name;
	char buf[ILL_namebufsize];
	char p1[2], p2[3];

	p1[0] = pref;
	p1[1] = '\0';
	p2[0] = pref;
	p2[1] = '_';
	p2[2] = '\0';

	ILL_SAFE_MALLOC (n_names, nnames + 1, char *);

	for (i = 0; i < nnames; i++)
	{
		n_names[i] = (char *) NULL;
	}

	for (i = 0; i <= nnames; i++)
	{
		if (i == nnames)
		{
			if (extra == NULL)
				break;
			old_name = extra;
		}
		else
			old_name = names[i];

		n = strlen (old_name);
		strcpy (buf, old_name);
		if (!ILLis_lp_name_char (buf[0], 1))
		{
			sprintf (buf, "%d", i);
		}
		else
		{
			for (j = 1; j < n; j++)
			{
				if (!ILLis_lp_name_char (buf[j], j))
				{
					sprintf (buf, "%d", i);
					break;
				}
			}
		}

		if (!ILLis_lp_name_char (buf[0], 0))
		{
			if (symtab == NULL)
			{
				symtab = &symt;
				ILLsymboltab_init (symtab);
				ILLsymboltab_create (symtab, nnames + 1);
				for (j = 0; j < nnames; j++)
				{
					ILLsymboltab_register (symtab, names[j], -1, &ind, &hit);
					ILL_FAILfalse (ind == j, "ind == j");
				}
				if (extra != NULL)
					ILLsymboltab_register (symtab, extra, -1, &ind, &hit);
			}
			rval = ILLsymboltab_uname (symtab, buf, p1, p2);
			CHECKRVALG (rval, CLEANUP);
			rval = ILLsymboltab_rename (symtab, i, buf);
			CHECKRVALG (rval, CLEANUP);

			ILL_UTIL_STR (n_names[i], buf);
			ILLdata_warn (collector,
										"\"%s\" is not a valid name in LP format; %s\"%s\".",
										old_name, "renaiming to ", buf);
		}
		else
		{
			ILL_UTIL_STR (n_names[i], old_name);
		}
	}

CLEANUP:
	if (symtab != NULL)
	{
		ILLsymboltab_free (symtab);
	}
	*newnames = n_names;
	EG_RETURN (rval);
}

/* 
 * end ILLlpdata_lpwrite 
 * ---------------------------------------------------------------------- */

int ILLread_lp (
	qsline_reader * file,
	const char *fname,
	rawlpdata * lp)
{
/* file format: 
 *     optional problem name (this is in addition to the bix format) 
 *     min or max 
 *     objective fct 
 *     constraints (mandatory !) 
 *     optional bounds 
 *     optional integer (also new) 
 *
 * as opposed to the official bix-format: 
 *     no blanks in variable names 
 *     there may be several bound defs on the same line; 
 *     bound definitions may cross line boundaries
 *     constraints that have no name get generated names:
 *             if constraint i has not name we take the first name from the 
 *             following list that is not in use:
 *                  c<i>,  C<i>, or c<i>_0, c<i>_1, c<i>_2, ...
 */
	int rval = 0;
	ILLread_lp_state lpstate, *state = &lpstate;

	const char *bnds[3], *integer[3], *end[2];

	bnds[0] = "BOUNDS";
	bnds[1] = "BOUND";
	bnds[2] = NULL;
	integer[0] = "INTEGER";
	integer[1] = "INT";
	integer[2] = NULL;
	end[0] = "END";
	end[1] = NULL;

	rval = ILLread_lp_state_init (state, file, fname, 0);
	CHECKRVALG (rval, CLEANUP);

	ILLinit_rawlpdata (lp, file->error_collector);
	rval = ILLsymboltab_create (&lp->rowtab, 100) ||
		ILLsymboltab_create (&lp->coltab, 100);
	CHECKRVALG (rval, CLEANUP);

	if (ILLread_lp_state_next_field (state))
	{
		rval = ILLlp_error (state, "Empty file.\n");
	}
	if (rval == 0)
		rval = read_problem_name (state, lp);
	if (rval == 0)
		rval = read_minmax (state, lp);
	if (rval == 0)
		rval = read_objective (state, lp);
	if (rval == 0)
		rval = read_constraints (state, lp, 1);
	if ((rval == 0) && (lp->ncols == 0 || lp->nrows == 0))
	{
		rval = ILLlp_error (state,
												"Problem must contain at least one %s.\n",
												"non empty constraint");
	}
	CHECKRVALG (rval, CLEANUP);

	if (ILLread_lp_state_keyword (state, bnds) == 0)
	{
		rval = read_bounds (state, lp);
	}
	CHECKRVALG (rval, CLEANUP);

	if (ILLread_lp_state_keyword (state, integer) == 0)
	{
		rval = read_integer (state, lp);
	}
	CHECKRVALG (rval, CLEANUP);

	rval = ILLread_lp_state_keyword (state, end);
	if (rval != 0)
	{
		if (state->eof)
		{
			rval = ILLlp_error (state, "Missing \"End\" at end of file.\n");
		}
		else
		{
			rval = ILLlp_error (state, "\"%s\" unknown keyword\n", state->field);
		}
	}
	if (rval == 0)
	{
		rval = ILLraw_fill_in_rownames (lp) || ILLraw_fill_in_bounds (lp);
	}

CLEANUP:
	EGlpNumClearVar (lpstate.bound_val);
	EG_RETURN (rval);
}

static int read_problem_name (
	ILLread_lp_state * state,
	rawlpdata * lp)
{
	int rval = 0;

	if (!state->fieldOnFirstCol)
	{
		rval = ILLlp_error (state,
												"Keyword \"%s\" not at beginning of line.\n",
												state->field);
	}
	if (!ILLutil_strcasecmp (state->field, "PROBLEM") ||
			!ILLutil_strcasecmp (state->field, "PROB"))
	{
		if (ILLread_lp_state_next_field (state) != 0)
		{
			rval = ILLlp_error (state, "No Problem name field.\n");
		}
		else
		{
			ILL_IFFREE (lp->name, char);

			ILL_UTIL_STR (lp->name, state->field);
			ILL_IFTRACE ("ProblemName: %s\n", state->field);
			(void) ILLread_lp_state_next_field (state);
		}
	}
CLEANUP:
	EG_RETURN (rval);
}

static int read_minmax (
	ILLread_lp_state * state,
	rawlpdata * lp)
{
	int rval = 0;

	if (!state->fieldOnFirstCol)
	{
		rval = ILLlp_error (state,
												"Keyword \"%s\" not at beginning of line.\n",
												state->field);
	}
	if (!ILLutil_strcasecmp (state->field, "MAX") ||
			!ILLutil_strcasecmp (state->field, "MAXIMUM") ||
			!ILLutil_strcasecmp (state->field, "MAXIMIZE"))
	{
		lp->objsense = ILL_MAX;
	}
	else
	{
		if (!ILLutil_strcasecmp (state->field, "MIN") ||
				!ILLutil_strcasecmp (state->field, "MINIMUM") ||
				!ILLutil_strcasecmp (state->field, "MINIMIZE"))
		{
			lp->objsense = ILL_MIN;
		}
		else
		{
			ILLread_lp_state_prev_field (state);
			rval = ILLlp_error (state, "Expecting \"%s\" or \"%s\" keyword.\n",
													"Minimize", "Maximize");
		}
	}
	EG_RETURN (rval);
}

int ILLread_constraint_expr (
	ILLread_lp_state * state,
	rawlpdata * lp,
	int rowind,
	int allowNew)
{
	int rval = 0;
	char firstTerm, haveCoef;
	const char *name;
	EGlpNum_t sign, coef;
	EGlpNum_t ntmp;

	EGlpNumInitVar (ntmp);
	EGlpNumInitVar (sign);
	EGlpNumInitVar (coef);

	firstTerm = 1;
	while (1)
	{
		if (ILLread_lp_state_sign (state, &sign) != 0)
		{
			if (!firstTerm)
			{
				break;									/* we've ssen at least one term, 
																 * this is the constraint's end */
			}
		}
		haveCoef = ILLread_lp_state_possible_coef (state, &coef, oneLpNum);
		if (ILLread_lp_state_next_var (state) == 0)
		{
			EGlpNumCopy (ntmp, coef);
			EGlpNumMultTo (ntmp, sign);
			rval = add_var (lp, state, ntmp, rowind, allowNew);
			CHECKRVALG (rval, CLEANUP);
		}
		else
		{
			if (haveCoef == 0)
			{
				return ILLlp_error (state, "Coefficient without variable.\n");
			}
			else
			{
				break;
			}
		}
		firstTerm = 0;
	}
CLEANUP:
	if ((rval == 0) && firstTerm)
	{
		name = ILLraw_rowname (lp, rowind);
		if (name != NULL)
		{
			ILLlp_warn (state,
									"No terms in constraint expression for \"%s\".\n", name);
		}
		else
		{
			ILLlp_warn (state, "No terms in constraint expression.\n");
		}
	}
	EGlpNumClearVar (ntmp);
	EGlpNumClearVar (sign);
	EGlpNumClearVar (coef);
	EG_RETURN (rval);
}

static int read_objective (
	ILLread_lp_state * state,
	rawlpdata * lp)
{
	int rval = 0;
	char objname[ILL_namebufsize];
	char *name;

	ILL_FAILfalse (lp->nrows == 0, "objective should be first row");
	ILLread_lp_state_skip_blanks (state, 1);
	if (ILLread_lp_state_has_colon (state))
	{
		if (ILLread_lp_state_next_var (state) != 0)
		{
			rval = ILLlp_error (state, "Bad objective function name.\n");
		}
		name = state->field;
		if (rval == 0)
		{
			if (ILLread_lp_state_colon (state) != 0)
			{
				rval = ILLlp_error (state, "':' must follow constraint row name.\n");
			}
		}
	}
	else
	{
		name = NULL;
	}

	if (rval == 0)
	{
		ILL_FAILfalse (lp->rowtab.tablesize == 0,
									 "objective row is first in symbol tab");
		if (name == NULL)
		{
			strcpy (objname, "obj");
			ILLlp_warn (state, "Empty obj name; using \"%s\".\n", objname);
		}
		else
		{
			strcpy (objname, name);
		}
		rval = ILLraw_add_row (lp, objname, 'N', zeroLpNum);
		lp->objindex = lp->nrows - 1;
		CHECKRVALG (rval, CLEANUP);
		rval = ILLread_constraint_expr (state, lp, lp->objindex, 1);
	}
CLEANUP:
	EG_RETURN (rval);
}

int ILLread_constraint_name (
	ILLread_lp_state * state,
	char **rowname)
{
	int rval = 0;

	*rowname = NULL;

	/* if there is a ':' on the line: look for constraint row name */
	if (ILLread_lp_state_has_colon (state))
	{
		if (ILLread_lp_state_next_var (state) != 0)
		{
			rval = ILLlp_error (state, "Bad constraint row name.\n");
		}
		else
		{
			*rowname = state->field;
			if (ILLread_lp_state_colon (state) != 0)
			{
				rval = ILLlp_error (state, "':' must follow constraint row name.\n");
			}
		}
	}
	return rval;
}

int ILLread_one_constraint (
	ILLread_lp_state * state,
	const char *rowname,
	rawlpdata * lp,
	int allowNewCols)
{
	int rval = 0;
	int rowind;
	char sense;
	EGlpNum_t d;

	EGlpNumInitVar (d);

	if ((rowname != NULL) &&
			(ILLsymboltab_lookup (&lp->rowtab, rowname, &rowind) == 0))
	{
		rval = ILLlp_error (state, "Repeated row name \"%s\".\n", rowname);
		CHECKRVALG (rval, CLEANUP);
	}
	rowind = lp->nrows;
	rval = rval || ILLraw_add_row (lp, rowname, 'N', zeroLpNum);

	rval = rval || ILLread_constraint_expr (state, lp, rowind, allowNewCols);
	rval = rval || ILLread_lp_state_sense (state);
	sense = state->sense_val;
	if (rval == 0)
	{
		rval = ILLread_lp_state_value (state, &d);
		if (rval)
		{
			(void) ILLlp_error (state, "No right hand side value in constraint.\n");
		}
	}
	if (rval == 0)
	{
		lp->rowsense[rowind] = sense;
		EGlpNumCopy (lp->rhs[rowind], d);
		ILL_IFTRACE ("SENSE \"%s\": %c %f\n",
								 ILLraw_rowname (lp, rowind), sense, EGlpNumToLf (d));
	}
CLEANUP:
	EGlpNumClearVar (d);
	EG_RETURN (rval);
}

static int read_constraints (
	ILLread_lp_state * state,
	rawlpdata * lp,
	int allowNewCols)
{
	int rval = 0;
	char *rowname = NULL;

	if (ILLcheck_subject_to (state) != 0)
	{
		return ILLlp_error (state, "Constraint section expected.\n");
	}
	while (rval == 0)
	{
		rval = ILLread_constraint_name (state, &rowname);
		if (rval == 0)
		{
			rval = ILLread_one_constraint (state, rowname, lp, allowNewCols);
		}
		if (rval == 0)
		{
			if (ILLread_lp_state_next_constraint (state) != 0)
			{
				break;
			}
		}
	}
	ILLread_lp_state_next_field (state);
	EG_RETURN (rval);
}

/* 
 * return -2 iff next is not a variable and not a keyword
 * return -1 iff next is a keyword and !mustHave 
 * return 1  iff unknown column name or mustHave and keyword 
 * return 0  for success
 */
static int read_colname (
	ILLread_lp_state * state,
	ILLsymboltab * coltab,
	int mustHave)
{
	int rval = 0;
	int colind = ILL_SYM_NOINDEX;

	rval = ILLread_lp_state_next_var (state);
	if (mustHave && (rval != 0))
	{
		return ILLlp_error (state, "Expecting a column name.\n");
	}
	if (rval != 0)
	{
		return (rval == -1) ? rval : -2;
	}
	if (ILLsymboltab_lookup (coltab, state->field, &colind))
	{
		ILLread_lp_state_prev_field (state);
		return ILLlp_error (state, "\"%s\" is not a column name.\n", state->field);
	}
	state->column_index = colind;
	return 0;
}

static int read_integer (
	ILLread_lp_state * state,
	rawlpdata * lp)
{
	int rval = 0;
	ILLsymboltab *coltab = &lp->coltab;

	ILL_FAILfalse (lp->intmarker, "Programming error");

	while ((rval = read_colname (state, coltab, 0)) == 0)
	{
		ILL_FAILtrue (state->column_index == ILL_SYM_NOINDEX, "Programming error");
		lp->intmarker[state->column_index] = 1;
	}
CLEANUP:
	if (rval == -1)
	{															/* last try for a colname gave us a keyword */
		rval = 0;
	}
	else
	{
		rval = ILLlp_error (state, "Expecting a column name.");
	}
	ILLread_lp_state_next_field (state);
	EG_RETURN (rval);
}

static int read_bounds (
	ILLread_lp_state * state,
	rawlpdata * lp)
{
	int rval = 0;
	int colind, haveBound;
	char sense;
	const char *msg;
	ILLsymboltab *coltab;

	ILLraw_init_bounds (lp);
	coltab = &lp->coltab;

	while (1)
	{
		colind = -1;
		haveBound = 0;
		if (ILLread_lp_state_possible_bound_value (state))
		{
			/* this must be for a lower bound */
			ILLtest_lp_state_bound_sense (state);
			if (state->sense_val != 'L')
			{
				rval = ILLlp_error (state, "Expecting \"<=\".\n");
				break;
			}
			rval = read_colname (state, coltab, 1);
			if (rval != 0)
			{
				break;
			}
			colind = state->column_index;
			/* add lower bound value */
			msg = ILLraw_set_lowerBound (lp, colind, state->bound_val);
			ILLlp_warn (state, msg);
			haveBound = 1;
		}
		if (colind == -1)
		{
			rval = read_colname (state, coltab, 0);
			colind = state->column_index;
			if (rval != 0)
			{
				if (rval == -1)
				{
					rval = 0;							/* found a keyword and that's OK */
				}
				else if (rval == -2)
				{
					rval = ILLlp_error (state, "Expecting a column name.\n");
				}
				break;
			}
		}
		ILL_FAILtrue (colind == -1, "must have a valid colname");
		ILLtest_lp_state_bound_sense (state);
		if (state->sense_val != ' ')
		{
			sense = state->sense_val;
			if ((sense != 'L') && (sense != 'E'))
			{
				rval = ILLlp_error (state, "Expecting \"<=\" or \"=\".\n");
				break;
			}
			if (ILLread_lp_state_possible_bound_value (state))
			{
				if (sense == 'E')
				{
					msg = ILLraw_set_fixedBound (lp, colind, state->bound_val);
				}
				else
				{
					msg = ILLraw_set_upperBound (lp, colind, state->bound_val);
				}
				ILLlp_warn (state, msg);
				haveBound = 1;
			}
			else
			{
				rval = ILLlp_error (state, "Expecting bound value.\n");
				break;
			}
		}
		else
		{
			if (ILLtest_lp_state_next_is (state, "FREE"))
			{
				msg = ILLraw_set_unbound (lp, colind);
				ILLlp_warn (state, msg);
				haveBound = 1;
			}
			else
			{
				if (!haveBound)
				{
					rval = ILLlp_error (state, "Not a bound expression.\n");
					break;
				}
			}
		}
		ILL_IFTRACE ("BOUNDS: %f <= %s <= %f\n", EGlpNumToLf (lp->lower[colind]),
								 ILLraw_colname (lp, colind), EGlpNumToLf (lp->upper[colind]));
	}
	ILLread_lp_state_next_field (state);
CLEANUP:
	EG_RETURN (rval);
}

static int add_var (
	rawlpdata * lp,
	ILLread_lp_state * state,
	EGlpNum_t coef,
	int row,
	int allowNew)
{
	char *var = state->field;
	int rval = 0;
	int colind;

	if (ILLsymboltab_lookup (&lp->coltab, var, &colind) != 0)
	{
		if (!allowNew)
		{
			rval = ILLlp_error (state, "Unknown col name \"%s\".\n", var);
		}
		CHECKRVALG (rval, CLEANUP);
		rval = ILLraw_add_col (lp, var, 0 /* not an integer var */ );
		colind = lp->ncols - 1;
		CHECKRVALG (rval, CLEANUP);
	}
	ILL_IFTRACE ("add_var: \"%s\" coef=%f row=%s\n",
							 var, EGlpNumToLf (coef), ILLraw_rowname (lp, row));
	rval = ILLraw_add_col_coef (lp, colind, row, coef);
CLEANUP:
	EG_RETURN (rval);
}
