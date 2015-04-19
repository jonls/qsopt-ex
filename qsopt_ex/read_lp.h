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

/*  RCS_INFO = "$RCSfile: read_lp_state.h,v $ $Revision: 1.2 $ $Date: 2003/11/05 16:57:39 $"; */
#ifndef EGLPNUM_TYPENAME_READ_LP_STATE_H
#define EGLPNUM_TYPENAME_READ_LP_STATE_H

/****************************************************************************/
/*                                                                          */
/*               Routines to support Reading LP Files                       */
/*                                                                          */
/****************************************************************************/

/* 
 * -) anything after '\' is comment 
 * -) variables consist of a-z A-Z 0-9!"#$%(),;.?@_`'{}|~ 
 *    don't start with a digit or '.'
 */

#include "eg_lpnum.h"
#include "symtab.h"

#include "readline_EGLPNUM_TYPENAME.h"


typedef struct EGLPNUM_TYPENAME_ILLread_lp_state
{
	EGLPNUM_TYPENAME_qsline_reader *file;
	const char *file_name;
	char *p;
	EGLPNUM_TYPE bound_val;
	int interactive;
	int line_num;
	int column_index;
	char realline[ILL_namebufsize];
	char line[ILL_namebufsize];
	char field[ILL_namebufsize + 1];
	char fieldOnFirstCol;
	char eof;
	char sense_val;
}
EGLPNUM_TYPENAME_ILLread_lp_state;

extern int EGLPNUM_TYPENAME_ILLread_lp_state_init (
	EGLPNUM_TYPENAME_ILLread_lp_state * state,
	EGLPNUM_TYPENAME_qsline_reader * file,
	const char *fname,
	int interactve);
extern int EGLPNUM_TYPENAME_ILLread_lp_state_next_line (
	EGLPNUM_TYPENAME_ILLread_lp_state * state);
extern int EGLPNUM_TYPENAME_ILLread_lp_state_next_var (
	EGLPNUM_TYPENAME_ILLread_lp_state * state);
extern int EGLPNUM_TYPENAME_ILLread_lp_state_keyword (
	EGLPNUM_TYPENAME_ILLread_lp_state * state,
	const char **kwd);
extern int EGLPNUM_TYPENAME_ILLread_lp_state_bad_keyword (
	EGLPNUM_TYPENAME_ILLread_lp_state * state);
extern int EGLPNUM_TYPENAME_ILLtest_lp_state_keyword (
	EGLPNUM_TYPENAME_ILLread_lp_state * state,
	const char *kwd[]);
extern int EGLPNUM_TYPENAME_ILLread_lp_state_next_field (
	EGLPNUM_TYPENAME_ILLread_lp_state * state);
extern int EGLPNUM_TYPENAME_ILLread_lp_state_next_field_on_line (
	EGLPNUM_TYPENAME_ILLread_lp_state * state);
extern void EGLPNUM_TYPENAME_ILLread_lp_state_prev_field (
	EGLPNUM_TYPENAME_ILLread_lp_state * state);
extern int EGLPNUM_TYPENAME_ILLread_lp_state_sign (
	EGLPNUM_TYPENAME_ILLread_lp_state * state,
	EGLPNUM_TYPE * sign);
extern int EGLPNUM_TYPENAME_ILLread_lp_state_possible_coef (
	EGLPNUM_TYPENAME_ILLread_lp_state * state,
	EGLPNUM_TYPE * coef,
	const EGLPNUM_TYPE defValue);

																				/* returns 1 iff found a number 
																				 * otherwise 0 */
extern int EGLPNUM_TYPENAME_ILLread_lp_state_possible_bound_value (
	EGLPNUM_TYPENAME_ILLread_lp_state * state);

																							 /* returns 1 iff found a number 
																							  * otherwise 0 */
extern int EGLPNUM_TYPENAME_ILLread_lp_state_colon (
	EGLPNUM_TYPENAME_ILLread_lp_state * state);
extern int EGLPNUM_TYPENAME_ILLread_lp_state_has_colon (
	EGLPNUM_TYPENAME_ILLread_lp_state * state);
extern int EGLPNUM_TYPENAME_ILLread_lp_statxe_has_colon (
	EGLPNUM_TYPENAME_ILLread_lp_state * state);
extern int EGLPNUM_TYPENAME_ILLread_lp_state_next_constraint (
	EGLPNUM_TYPENAME_ILLread_lp_state * state);
extern int EGLPNUM_TYPENAME_ILLread_lp_state_sense (
	EGLPNUM_TYPENAME_ILLread_lp_state * state);
extern int EGLPNUM_TYPENAME_ILLtest_lp_state_sense (
	EGLPNUM_TYPENAME_ILLread_lp_state * state,
	int all);
extern void EGLPNUM_TYPENAME_ILLtest_lp_state_bound_sense (
	EGLPNUM_TYPENAME_ILLread_lp_state * state);
extern int EGLPNUM_TYPENAME_ILLread_lp_state_value (
	EGLPNUM_TYPENAME_ILLread_lp_state * state,
	EGLPNUM_TYPE * d);
extern int EGLPNUM_TYPENAME_ILLtest_lp_state_next_is (
	EGLPNUM_TYPENAME_ILLread_lp_state * state,
	const char *str);
extern int EGLPNUM_TYPENAME_ILLread_lp_state_skip_blanks (
	EGLPNUM_TYPENAME_ILLread_lp_state * state,
	int wrapLines);

extern int EGLPNUM_TYPENAME_ILLcheck_subject_to (
	EGLPNUM_TYPENAME_ILLread_lp_state * state);

/*---------------------------------------------------------------------------*/
/* errors and warnings 
 */
extern int EGLPNUM_TYPENAME_ILLlp_error (
	EGLPNUM_TYPENAME_ILLread_lp_state * state,
	const char *format,
	...);
extern void EGLPNUM_TYPENAME_ILLlp_warn (
	EGLPNUM_TYPENAME_ILLread_lp_state * state,
	const char *format,
	...);

/*---------------------------------------------------------------------------*/
/* shared with read_mps_state.c 
 */
extern int EGLPNUM_TYPENAME_ILLget_value (
	char *line,
	EGLPNUM_TYPE * coef);

#endif
