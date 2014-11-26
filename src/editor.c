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

/* RCS_INFO = "$RCSfile: editor.c,v $ $Revision: 1.2 $ $Date: 2003/11/05 16:49:52 $"; */

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <stdlib.h>
#include <string.h>

#include "qs_config.h"

#include "eg_lpnum.h"
#include "eg_io.h"

#include "qsopt.h"
#include "lpdata.h"
#include "qstruct.h"
#include "qsopt.h"
#include "editor.h"
#include "readline.h"
#include "rawlp.h"
#include "stddefs.h"						/* for MAX */
#include "read_lp.h"
#include "lp.h"
#include "lib.h"
#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

static int TRACE = 0;

#define ILL_BREAK_BODY_IF(rval) if (rval != 0) goto CLEANUP
#define ILL_BREAK_BODY goto CLEANUP

static int transpose ( rawlpdata * lp);
static int pull_info_from_p ( QSdata * p, rawlpdata * lp);
static void add_row ( QSdata * p, rawlpdata * lp, ILLread_lp_state * state);

/* static int new_row(QSdata *p, rawlpdata *lp, ILLread_lp_state *state); */
static void del_row ( QSdata * p, rawlpdata * lp, ILLread_lp_state * state);

static void add_col ( QSdata * p, rawlpdata * lp, ILLread_lp_state * state);
static void del_col ( QSdata * p, rawlpdata * lp, ILLread_lp_state * state);

#define NONE -1
#define QS_EXIT 0
#define ROW 1
#define COL 2
#define PLP 3
#define PRTX 4
#define SOLVE 5
#define PMPS 6
#define HELP 7
#define DEL 8
#define NEW 9
#define ADD 10
#define PRIMAL 11
#define DUAL 12
#define NCOMMAND 13
static const char *commands[NCOMMAND + 1];
static char hasSubCmd[NCOMMAND + 1];

void ILLeditor_init (
	void)
{
	commands[QS_EXIT] = "QS_EXIT";
	commands[ROW] = "ROW";
	commands[COL] = "COL";
	commands[PLP] = "LP";
	commands[PMPS] = "MPS";
	commands[SOLVE] = "SOLVE";
	commands[PRTX] = "PRT";
	commands[HELP] = "HELP";
	commands[ADD] = "ADD";
	commands[DEL] = "DEL";
	commands[NEW] = "NEW";
	commands[PRIMAL] = "PRIMAL";
	commands[DUAL] = "DUAL";
	commands[NCOMMAND] = NULL;

	hasSubCmd[QS_EXIT] = 0;
	hasSubCmd[ROW] = 1;
	hasSubCmd[COL] = 1;
	hasSubCmd[PLP] = 0;
	hasSubCmd[PMPS] = 0;
	hasSubCmd[SOLVE] = 1;
	hasSubCmd[PRTX] = 0;
	hasSubCmd[HELP] = 0;
	hasSubCmd[ADD] = 1;
	hasSubCmd[DEL] = 1;
	hasSubCmd[NEW] = 1;
	hasSubCmd[PRIMAL] = 1;
	hasSubCmd[DUAL] = 1;
	hasSubCmd[NCOMMAND] = 0;
}

static void ILLeditor_help_cmd (
	int cmd,
	int subcmd);

static void ILLeditor_help (
	void)
{
	ILLeditor_help_cmd (ROW, ADD);
	/* ILLeditor_help_cmd(ROW, NEW);  */
	ILLeditor_help_cmd (ROW, DEL);
	ILLeditor_help_cmd (COL, ADD);
	ILLeditor_help_cmd (COL, DEL);
	ILLeditor_help_cmd (SOLVE, NONE);
	ILLeditor_help_cmd (PRTX, NONE);
	ILLeditor_help_cmd (PLP, NONE);
	ILLeditor_help_cmd (PMPS, NONE);
	ILLeditor_help_cmd (QS_EXIT, NONE);
	ILLeditor_help_cmd (HELP, NONE);
}

static void ILLeditor_help_cmd (
	int cmd,
	int subcmd)
{
	if (cmd == ROW && subcmd == ADD)
		fprintf (stdout, "%s ADD:\t%s.\n",
						 commands[ROW], "add a row; enter in LP format");
	if (cmd == COL && subcmd == ADD)
		fprintf (stdout, "%s ADD:\t%s.\n",
						 commands[COL], "add a col; enter in LP format");
	/* if (cmd == ROW && subcmd == NEW) 
	 * fprintf(stdout, "%s NEW:\t%s.\n", 
	 * commands[ROW], "new row; enter rowname: sense rhs");  
	 */
	if (cmd == ROW && subcmd == DEL)
		fprintf (stdout, "%s DEL:\t%s.\n",
						 commands[ROW], "delete a row; give rowname");
	if (cmd == COL && subcmd == DEL)
		fprintf (stdout, "%s DEL:\t%s.\n",
						 commands[COL], "delete a col; give colname");
	if (cmd == SOLVE)
		fprintf (stdout, "%s:\t%s.\n", commands[SOLVE], "solve problem");
	if (cmd == PRTX)
		fprintf (stdout, "%s:\t%s.\n",
						 commands[PRTX], "print variable values for optimal solution");
	if (cmd == PLP)
		fprintf (stdout, "%s [file]:\t%s.\n",
						 commands[PLP], "print problem in LP format to file or stdout");
	if (cmd == PMPS)
		fprintf (stdout, "%s [file]:\t%s.\n",
						 commands[PMPS], "print problem in MPS format to file or stdout");
	if (cmd == QS_EXIT)
		fprintf (stdout, "%s:\t%s.\n", commands[QS_EXIT], "QS_EXIT");
	if (cmd == HELP)
		fprintf (stdout, "%s:\t%s.\n", commands[HELP], "print this help");
}

static void getCmd (
	ILLread_lp_state * state,
	int *cmd,
	int *subcmd)
{
	const char *cmd_str, *subcmd_str;
	int tmp;

	*cmd = ILLutil_index (commands, state->field);
	*subcmd = -1;
	if (hasSubCmd[*cmd] && (ILLread_lp_state_next_field_on_line (state) == 0))
	{
		*subcmd = ILLutil_index (commands, state->field);
		if ((*subcmd == ROW) || (*subcmd == COL) || (*subcmd == SOLVE))
		{
			ILL_SWAP (*subcmd, *cmd, tmp);
		}
	}
	cmd_str = (*cmd >= 0) ? commands[*cmd] : "???";
	subcmd_str = (*subcmd >= 0) ? commands[*subcmd] : "???";
	ILL_IFTRACE ("cmd = %s, subcmd = %s\n", cmd_str, subcmd_str);
}

void ILLeditor (
	QSdata * p)
{
	rawlpdata raw, *lp = &raw;
	int cmd, subcmd, tval, rval = 0;
	ILLread_lp_state lpstate, *state = &lpstate;
	qsline_reader *reader;

	ILL_IFTRACE ("ILLeditor\n");

	reader = ILLline_reader_new ((qsread_line_fct) fgets, stdin);
	rval = ILLread_lp_state_init (state, reader, "STDIN", 1);
	rval = rval || pull_info_from_p (p, lp);
	ILL_BREAK_BODY_IF (rval);

	while (ILLread_lp_state_next_field (state) == 0)
	{
		getCmd (state, &cmd, &subcmd);
		switch (cmd)
		{
		case QS_EXIT:
			ILL_BREAK_BODY;

		case ROW:
			{
				switch (subcmd)
				{
				case ADD:
					add_row (p, lp, state);
					break;
					/* case NEW: rval = new_row(p, lp, state); break; */
				case DEL:
					del_row (p, lp, state);
					break;
				default:
					ILLeditor_help ();
					break;
				}
				break;
			}
		case COL:
			{
				switch (subcmd)
				{
				case ADD:
					add_col (p, lp, state);
					break;
				case DEL:
					del_col (p, lp, state);
					break;
				default:
					ILLeditor_help ();
					break;
				}
				break;
			}

		case SOLVE:
			{
				if (subcmd == PRIMAL)
				{
					(void) ILLeditor_solve (p, PRIMAL_SIMPLEX);
				}
				else if (subcmd == DUAL)
				{
					(void) ILLeditor_solve (p, DUAL_SIMPLEX);
				}
				else
				{
					ILLeditor_help ();
				}
				break;
			}

		case PRTX:
			{
				EGioFile_t*lout = EGioOpenFILE(stdout);
				if ((rval = ILLlib_print_x (lout, p->lp, 0, 0, 1)))
				{
					fprintf (stdout, "The problem may not be feasible.\n");
				}
				EGioClose(lout);
				break;
			}

		case PLP:
		case PMPS:
			{
				if (ILLread_lp_state_next_field_on_line (state) == 0)
				{
					if (cmd == PMPS)
					{
						tval = QSwrite_prob (p, state->field, "MPS");
					}
					else
					{
						tval = QSwrite_prob (p, state->field, "LP");
					}
					if (tval)
					{
						fprintf (stdout, "Could not write problem to \"%s\".\n",
										 state->field);
					}
					else
					{
						fprintf (stdout, "Saved to \"%s\".\n", state->field);
					}
				}
				else
				{
					if (cmd == PMPS)
					{
						(void) QSwrite_prob_file (p, stdout, "MPS");
					}
					else
					{
						(void) QSwrite_prob_file (p, stdout, "LP");
					}
				}
				break;
			}

		case NONE:
			fprintf (stdout, "Unknown command: %s\n", state->field);
		default:
			ILLeditor_help ();
			break;
		}
		fflush (stdout);
		ILLread_lp_state_next_line (state);
	}
CLEANUP:
	ILLline_reader_free (reader);
	ILLfree_rawlpdata (lp);
}

int ILLeditor_solve (
	QSdata * p,
	int salgo)
{
	int rval = 0;
	int status = 0;
	EGlpNum_t val;

	EGlpNumInitVar (val);

	if (salgo == PRIMAL_SIMPLEX)
	{
		rval = QSopt_primal (p, &status);
	}
	else
	{
		rval = QSopt_dual (p, &status);
	}
	ILL_BREAK_BODY_IF (rval);
	rval = QSget_objval (p, &val);
	if (p->simplex_display)
		if (rval == 0)
		{
			fprintf (stdout, "LP Value: %.6f, status %d\n", EGlpNumToLf (val),
							 status);
			fflush (stdout);
		}
CLEANUP:
	EGlpNumClearVar (val);
	ILL_RESULT (rval, "ILLeditor_solve");
}


static int pull_info_from_p (
	QSdata * p,
	rawlpdata * lp)
{
	int i, rval = 0;
	ILLlpdata *qslp = p->lp->O;
	int nrows, ncols;

	ILLinit_rawlpdata (lp, NULL);
	rval = ILLsymboltab_create (&lp->rowtab, 100) ||
		ILLsymboltab_create (&lp->coltab, 100);
	ILL_BREAK_BODY_IF (rval);

	nrows = qslp->nrows;
	ncols = qslp->nstruct;
	/* add rows to lp */
	ILLraw_add_row (lp, qslp->objname, 'N', zeroLpNum);
	for (i = 0; i < nrows; i++)
	{
		ILL_FAILfalse (qslp->rownames[i] != NULL, "should have no NULL names");
		ILLraw_add_row (lp, qslp->rownames[i], qslp->sense[i], qslp->rhs[i]);
	}

	/* add cols to coltab and lp */
	for (i = 0; i < ncols; i++)
	{
		ILL_FAILfalse (qslp->colnames[i] != NULL, "should have no NULL names");
		ILLraw_add_col (lp, qslp->colnames[i],
										(qslp->intmarker) ? qslp->intmarker[i] : 0);
	}
CLEANUP:
	ILL_RETURN (rval, "pull_info_from_p");
}

static int transpose (
	rawlpdata * lp)
{
	int rval = 0;
	int tmp;
	ILLsymboltab tmptab;

	tmp = QSMAX (lp->nrows, lp->ncols);
	if (tmp >= lp->sensesize)
	{
		lp->sensesize *= 1.3;
		lp->sensesize += 1000;
		if (lp->sensesize < tmp + 1)
			lp->sensesize = tmp + 1;
		lp->rowsense = EGrealloc (lp->rowsense, sizeof (char) * lp->sensesize);
		//rval = ILLutil_reallocrus_scale ((void **) &lp->rowsense,
		//                                 &lp->sensesize, tmp + 1,
		//                                 1.3, sizeof (char));
		//ILL_CLEANUP_IF (rval);
	}
	if (tmp >= lp->rhssize)
	{
		lp->rhssize *= 1.3;
		lp->rhssize += 1000;
		if (lp->rhssize < tmp + 1)
			lp->rhssize = tmp + 1;
		EGlpNumReallocArray (&(lp->rhs), lp->rhssize);
		//lp->rhs = EGrealloc(lp->rhs, sizeof(double)*lp->rhssize);
		//rval = ILLutil_reallocrus_scale ((void **) &lp->rhs,
		//                                 &lp->sensesize, tmp + 1,
		//                                 1.3, sizeof (double));
		//ILL_CLEANUP_IF (rval);
	}
	ILL_SWAP (lp->nrows, lp->ncols, tmp);
	ILL_SWAP (lp->rowtab, lp->coltab, tmptab);
	ILL_RETURN (rval, "transpose");
}

static char *get_row_col_name (
	QSdata * p,
	rawlpdata * lp,
	ILLread_lp_state * state,
	int doRow)
{
	int rval = 0;
	int ind;
	char *rname, *thename = NULL;
	char buf[ILL_namebufsize];
	ILLsymboltab *tab = (doRow) ? &lp->rowtab : &lp->coltab;
	int id = (doRow) ? lp->nrows : lp->ncols;

	id--;													/* in rawlpdata obj counts as a row */

	rval = ILLread_constraint_name (state, &rname);
	ILL_BREAK_BODY_IF (rval);

	if (rname == NULL)
	{
		ILLlib_findName (p->qslp, doRow /* forRow */ , rname, id, buf);
		ILL_UTIL_STR (thename, buf);
	}
	else
	{
		ILL_UTIL_STR (thename, rname);
	}
	if (ILLsymboltab_lookup (tab, thename, &ind) == 0)
	{
		rval = ILLlp_error (state, "\"%s\" already exists.", thename);
	}
CLEANUP:
	if (rval != 0)
	{
		ILL_IFFREE (thename, char);
	}
	return thename;
}

static int fill_matrix (
	rawlpdata * lp,
	ILLread_lp_state * state,
	ILLmatrix * m,
	EGlpNum_t * obj,
	int n)
{
	int i, cnt, rval = 0;
	colptr *cp;
	EGlpNum_t val;
	int newCol = (obj != NULL);

	EGlpNumInitVar (val);

	/* rely on fact that objective has rowindex 0 */

	m->matrows = lp->nrows;
	m->matcols = 1;
	m->matval = EGlpNumAllocArray (lp->ncols);
	ILL_SAFE_MALLOC (m->matind, lp->ncols, int);
	ILL_SAFE_MALLOC (m->matbeg, 1, int);
	ILL_SAFE_MALLOC (m->matcnt, 1, int);

	m->matsize = lp->ncols;
	m->matbeg[0] = 0;
	m->matcnt[0] = 0;
	for (i = 0; i < lp->ncols; i++)
	{
		cnt = 0;
		EGlpNumZero (val);
		for (cp = lp->cols[i]; cp != NULL; cp = cp->next)
		{
			ILL_FAILfalse (cp->this_val == n, "n should be the only row around");
			if (EGlpNumIsNeqqZero (cp->coef))
			{
				EGlpNumAddTo (val, cp->coef);
				cnt++;
			}
		}
		if (cnt > 1)
		{
			ILLlp_warn (state, "Multiple coefficients for \"%s\".",
									ILLraw_colname (lp, i));
		}
		if (EGlpNumIsNeqqZero (val))
		{
			if ((i - newCol) >= 0)
			{
				EGlpNumCopy (m->matval[m->matcnt[0]], val);
				m->matind[m->matcnt[0]] = i - newCol;
				m->matcnt[0]++;
			}
			else
			{
				EGlpNumCopy (obj[0], val);
			}
		}
	}
	if (m->matcnt[0] == 0)
	{
		rval = ILLlp_error (state, "There are no non zero coefficients.");
	}
CLEANUP:
	EGlpNumClearVar (val);
	ILL_RESULT (rval, "fill_matrix");
}

static void add_row (
	QSdata * p,
	rawlpdata * lp,
	ILLread_lp_state * state)
{
	int rval = 0;
	int n;
	char *name;
	ILLmatrix m;
	char sense[1];

	ILLmatrix_init (&m);
	n = lp->nrows;
	name = get_row_col_name (p, lp, state, 1 /*doRow */ );

	if (name == NULL)
	{
		rval = 1;
	}
	else
	{
		rval = ILLread_one_constraint (state, name, lp, 0);

		/* adds row name to lp->rowtab the checks constraint expression  */
		if (rval != 0)
		{
			/* failed because of error in expression => 
			 * must remove name from symbol table */
			fprintf (stdout, "Incorrect expression.\n");
		}
		else
		{
			ILL_FAILfalse (lp->nrows == (n + 1), "Should have one row");
			ILL_IFTRACE ("ADDING row %s.\n", name);

			sense[0] = lp->rowsense[n];

			rval = fill_matrix (lp, state, &m, NULL, n);
			ILL_BREAK_BODY_IF (rval);

			QSadd_rows (p, 1, m.matcnt, m.matbeg, m.matind, m.matval,
									&(lp->rhs[n]), sense, (const char **) &name);
		}
	}
CLEANUP:
	ILLmatrix_free (&m);
	if (name != NULL)
	{
		if (rval != 0)
			ILLsymboltab_delete (&lp->rowtab, name);
		ILL_IFFREE (name, char);
	}
	if (rval != 0)
	{
		lp->nrows = n;
	}
	ILLraw_clear_matrix (lp);
}

static void add_col (
	QSdata * p,
	rawlpdata * lp,
	ILLread_lp_state * state)
{
	int rval = 0;
	int n;
	char *name[1];
	int transposed = 1;
	ILLmatrix matrix, *m = &matrix;
	EGlpNum_t obj[1], lower[1], upper[2];

	EGlpNumInitVar (*obj);
	EGlpNumInitVar (*lower);
	EGlpNumInitVar (upper[0]);
	EGlpNumInitVar (upper[1]);

	n = lp->ncols;
	ILLmatrix_init (m);
	name[0] = get_row_col_name (p, lp, state, 0 /*doRow */ );
	rval = (name[0] == NULL);
	ILL_BREAK_BODY_IF (rval);

	transposed = !transpose (lp);
	rval = ILLread_one_constraint (state, name[0], lp, 0);

	/* adds row name to lp->rowtab the checks constraint expression  */
	if (rval != 0)
	{
		/* failed because of error in expression => 
		 * must remove name from symbol table */
		fprintf (stdout, "Incorrect expression.\n");
	}
	else
	{
		ILL_FAILfalse (lp->nrows == (n + 1), "Should have one row");

		rval = fill_matrix (lp, state, m, obj, n);
		ILL_BREAK_BODY_IF (rval);

		fprintf (stdout, "lower ");
		rval = ILLread_lp_state_next_line (state) ||
			ILLread_lp_state_value (state, &(lower[0]));
		ILL_BREAK_BODY_IF (rval);

		fprintf (stdout, "upper ");
		rval = ILLread_lp_state_next_line (state) ||
			ILLread_lp_state_value (state, &(upper[0]));
		ILL_BREAK_BODY_IF (rval);

		ILL_IFTRACE ("ADDING col %s.\n", name[0]);

		QSadd_cols (p, 1, m->matcnt, m->matbeg, m->matind, m->matval,
								obj, lower, upper, (const char **) name);

	}
CLEANUP:
	ILLmatrix_free (m);
	if (name[0] != NULL)
	{
		if (rval != 0)
			ILLsymboltab_delete (&lp->rowtab, name[0]);
		ILL_IFFREE (name[0], char);
	}
	if (rval != 0)
	{
		lp->nrows = n;
	}
	ILLraw_clear_matrix (lp);
	if (transposed)
		transpose (lp);
	ILL_IFFREE (name[0], char);

	EGlpNumClearVar (*obj);
	EGlpNumClearVar (*lower);
	EGlpNumClearVar (upper[0]);
	EGlpNumClearVar (upper[1]);
}

#if 0
#ifndef JAVA_PORT
static void new_row (
	QSdata * p,
	rawlpdata * lp,
	ILLread_lp_state * state)
{
	int rval = 0;
	char *rowname = NULL, *rname = NULL;
	char sense;
	double d;
	int ind, hit;

	rval = ILLread_constraint_name (state, &rname);
	if (rname == NULL)
	{
		rval = 1;
		ILLeditor_help_cmd (ROW, NEW);
	}
	ILL_BREAK_BODY_IF (rval);

	ILLsymboltab_lookup (&lp->rowtab, rname, &ind);
	if (ind != ILL_SYM_NOINDEX)
	{
		rval = ILLlp_error (state, "\"%s\" is already defined.\n", rname);
		ILL_BREAK_BODY_IF (rval);
	}
	ILL_UTIL_STR (rowname, rname);

	rval = ILLread_lp_state_sense (state);
	sense = state->sense_val;
	ILL_BREAK_BODY_IF (rval);

	rval = ILLread_lp_state_value (state, &d);
	ILL_BREAK_BODY_IF (rval);

	rval = QSnew_row (p, d, sense, rowname);
	if (rval != 0)
	{
		fprintf (stderr, "could not add row\n");
	}
	else
	{
		ILLsymboltab_register (&lp->rowtab, rname, &ind, &hit);
	}
CLEANUP:
	ILL_IFFREE (rowname, char);
}
#endif
#endif

static int del_row_or_col (
	QSdata * p,
	rawlpdata * lp,
	ILLread_lp_state * state,
	int isRow)
{
	int i[1], rval = 0;
	char **names = (isRow) ? p->qslp->rownames : p->qslp->colnames;
	int nnames = (isRow) ? p->qslp->nrows : p->qslp->nstruct;
	ILLsymboltab *tab = (isRow) ? &lp->rowtab : &lp->coltab;

	rval = ILLread_lp_state_next_field_on_line (state);
	ILL_BREAK_BODY_IF (rval);

	i[0] = ILLutil_array_index (names, nnames, state->field);
	if (i[0] >= 0)
	{
		rval = (isRow) ? QSdelete_rows (p, 1, i) : QSdelete_cols (p, 1, i);
		if (rval == 0)
		{
			ILLsymboltab_delete (tab, state->field);
		}
	}
	else
	{
		rval = ILLlp_error (state, "\"%s\" is not defined.\n", state->field);
	}

CLEANUP:
	ILL_RESULT (rval, "del_row_or_col");
}

static void del_row (
	QSdata * p,
	rawlpdata * lp,
	ILLread_lp_state * state)
{
	int rval = del_row_or_col (p, lp, state, 1);

	if (rval == 0)
	{
		lp->nrows--;
	}
}

static void del_col (
	QSdata * p,
	rawlpdata * lp,
	ILLread_lp_state * state)
{
	int rval = del_row_or_col (p, lp, state, 0);

	if (rval == 0)
	{
		lp->ncols--;
	}
}
