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

/* RCSINFO $Id: lp_EGLPNUM_TYPENAME.h,v 1.2 2003/11/05 16:57:39 meven Exp $ */
#ifndef EGLPNUM_TYPENAME_LP_H
#define EGLPNUM_TYPENAME_LP_H

#include "readline_EGLPNUM_TYPENAME.h"

/****************************************************************************/
/*                                                                          */
/*               Routines to support Reading and Writing LP Files           */
/*                                                                          */
/****************************************************************************/
/* 
 * -) anything after '\' is comment 
 * -) Problem is optional and comes first
 * -) Minimize, Maximize, ... comes after Problem
 * -) variables consist of a-z A-Z 0-9!"#$%(),;.?@_`'{}|~ 
 *    don't start with a digit or '.'
 */

#include "lpdata_EGLPNUM_TYPENAME.h"
#include "rawlp_EGLPNUM_TYPENAME.h"
#include "read_lp_EGLPNUM_TYPENAME.h"
#include "write_lp_EGLPNUM_TYPENAME.h"

extern int EGLPNUM_TYPENAME_ILLread_lp (
	EGLPNUM_TYPENAME_qsline_reader * file,
	const char *fname,
	EGLPNUM_TYPENAME_rawlpdata * lp);
extern int EGLPNUM_TYPENAME_ILLwrite_lp (
	EGLPNUM_TYPENAME_ILLlpdata * l,
	EGLPNUM_TYPENAME_qserror_collector * collector);

			/* write using current lp->reporter */
extern int EGLPNUM_TYPENAME_ILLis_lp_name_char (
	int c,
	int pos);

extern int EGLPNUM_TYPENAME_ILLread_constraint_name (
	EGLPNUM_TYPENAME_ILLread_lp_state * state,
	char **rowname);
extern int EGLPNUM_TYPENAME_ILLread_one_constraint (
	EGLPNUM_TYPENAME_ILLread_lp_state * state,
	const char *rowname,
	EGLPNUM_TYPENAME_rawlpdata * lp,
	int allowNewColsAddRow);
extern int EGLPNUM_TYPENAME_ILLread_constraint_expr (
	EGLPNUM_TYPENAME_ILLread_lp_state * state,
	EGLPNUM_TYPENAME_rawlpdata * lp,
	int rowind,
	int allowNew);

#endif
