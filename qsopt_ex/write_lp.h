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

/*  RCS_INFO = "$RCSfile: wr_lp.h,v $ $Revision: 1.2 $ $Date: 2003/11/05 16:57:39 $"; */
#ifndef EGLPNUM_TYPENAME_WRITE_LP_STATE_H
#define EGLPNUM_TYPENAME_WRITE_LP_STATE_H

#include "eg_lpnum.h"
#include "symtab.h"

/****************************************************************************/
/*                                                                          */
/*               Routines to support writing of LP files                    */
/*                                                                          */
/****************************************************************************/

/* 
 * -) anything after '\' is comment 
 * -) variables consist of a-z A-Z 0-9!"#$%(),;.?@_`'{}|~ 
 *    don't start with a digit or '.'
 */

typedef struct EGLPNUM_TYPENAME_ILLwrite_lp_state
{
	char buf[ILL_namebufsize];
	char *p;
	int startlen;
	int total;
}
EGLPNUM_TYPENAME_ILLwrite_lp_state;

extern void EGLPNUM_TYPENAME_ILLwrite_lp_state_init (
	EGLPNUM_TYPENAME_ILLwrite_lp_state * line,
	const char *str);
extern void EGLPNUM_TYPENAME_ILLwrite_lp_state_append (
	EGLPNUM_TYPENAME_ILLwrite_lp_state * line,
	const char *str);
extern void EGLPNUM_TYPENAME_ILLwrite_lp_state_append_coef (
	EGLPNUM_TYPENAME_ILLwrite_lp_state * line,
	EGLPNUM_TYPE v,
	int cnt);

			/* append number sign ('+', '-') iff cnt > 0 or v < 0.0  
			 * append number iff v != 1.0, v != -1.0  
			 */
extern void EGLPNUM_TYPENAME_ILLwrite_lp_state_append_number (
	EGLPNUM_TYPENAME_ILLwrite_lp_state * line,
	EGLPNUM_TYPE v);
extern void EGLPNUM_TYPENAME_ILLwrite_lp_state_save_start (
	EGLPNUM_TYPENAME_ILLwrite_lp_state * line);
extern void EGLPNUM_TYPENAME_ILLwrite_lp_state_start (
	EGLPNUM_TYPENAME_ILLwrite_lp_state * line);

#endif
