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

/* $RCSfile: rd_mps.h,v $ $Revision: 1.2 $ $Date: 2003/11/05 16:57:39 $ */
#ifndef EGLPNUM_TYPENAME_READ_MPS_STATE_H
#define EGLPNUM_TYPENAME_READ_MPS_STATE_H

#include "symtab.h"
#include "basicdefs.h"

#include "readline_EGLPNUM_TYPENAME.h"


typedef struct EGLPNUM_TYPENAME_ILLread_mps_state_struct
{
	int section[ILL_MPS_N_SECTIONS];
	ILLmps_section active;
	const char *file_name;
	EGLPNUM_TYPENAME_qsline_reader *file;
	unsigned int line_num;
	unsigned int field_num;				/* number of successfully read fields on line */
	int intvar;
	int sosvar;
	char line[ILL_namebufsize];
	char key[ILL_namebufsize];
	char field[ILL_namebufsize];
	char *obj;
	char *p;											/* ptr to next 'unread' character */
}
EGLPNUM_TYPENAME_ILLread_mps_state;

extern int EGLPNUM_TYPENAME_ILLmps_state_init (
	EGLPNUM_TYPENAME_ILLread_mps_state * state,
	EGLPNUM_TYPENAME_qsline_reader * file,
	const char *fname);
extern void EGLPNUM_TYPENAME_ILLmps_state_clear (
	EGLPNUM_TYPENAME_ILLread_mps_state * state);
extern int EGLPNUM_TYPENAME_ILLmps_set_section (
	EGLPNUM_TYPENAME_ILLread_mps_state * state,
	const ILLmps_section sec);

extern int EGLPNUM_TYPENAME_ILLmps_next_line (
	EGLPNUM_TYPENAME_ILLread_mps_state * state);
extern int EGLPNUM_TYPENAME_ILLmps_next_field (
	EGLPNUM_TYPENAME_ILLread_mps_state * state);
extern int EGLPNUM_TYPENAME_ILLmps_next_coef (
	EGLPNUM_TYPENAME_ILLread_mps_state * state,
	EGLPNUM_TYPE * coef);
extern int EGLPNUM_TYPENAME_ILLmps_next_bound (
	EGLPNUM_TYPENAME_ILLread_mps_state * state,
	EGLPNUM_TYPE * coef);
extern void EGLPNUM_TYPENAME_ILLmps_check_end_of_line (
	EGLPNUM_TYPENAME_ILLread_mps_state * state);
extern void EGLPNUM_TYPENAME_ILLmps_set_end_of_line (
	EGLPNUM_TYPENAME_ILLread_mps_state * state);

extern int EGLPNUM_TYPENAME_ILLmps_int_sos_mode (
	EGLPNUM_TYPENAME_ILLread_mps_state * state);

extern const char *EGLPNUM_TYPENAME_ILLmps_possibly_blank_name (
	const char *field,
	EGLPNUM_TYPENAME_ILLread_mps_state * state,
	ILLsymboltab * tab);
extern int EGLPNUM_TYPENAME_ILLmps_empty_key (
	EGLPNUM_TYPENAME_ILLread_mps_state * state);
extern int EGLPNUM_TYPENAME_ILLmps_empty_field (
	EGLPNUM_TYPENAME_ILLread_mps_state * state);

extern int EGLPNUM_TYPENAME_ILLmps_error (
	EGLPNUM_TYPENAME_ILLread_mps_state * state,
	const char *format,
	...);
extern void EGLPNUM_TYPENAME_ILLmps_warn (
	EGLPNUM_TYPENAME_ILLread_mps_state * state,
	const char *format,
	...);

#endif
