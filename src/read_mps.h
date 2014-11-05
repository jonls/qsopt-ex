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
#ifndef READ_MPS_STATE_H
#define READ_MPS_STATE_H

#include "iqsutil.h"

#include "mps.h"

typedef struct ILLread_mps_state_struct
{
	int section[ILL_MPS_N_SECTIONS];
	ILLmps_section active;
	const char *file_name;
	qsline_reader *file;
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
ILLread_mps_state;

extern int ILLmps_state_init (
	ILLread_mps_state * state,
	qsline_reader * file,
	const char *fname);
extern void ILLmps_state_clear (
	ILLread_mps_state * state);
extern int ILLmps_set_section (
	ILLread_mps_state * state,
	const ILLmps_section sec);

extern int ILLmps_next_line (
	ILLread_mps_state * state);
extern int ILLmps_next_field (
	ILLread_mps_state * state);
extern int ILLmps_next_coef (
	ILLread_mps_state * state,
	EGlpNum_t * coef);
extern int ILLmps_next_bound (
	ILLread_mps_state * state,
	EGlpNum_t * coef);
extern void ILLmps_check_end_of_line (
	ILLread_mps_state * state);
extern void ILLmps_set_end_of_line (
	ILLread_mps_state * state);

extern int ILLmps_int_sos_mode (
	ILLread_mps_state * state);

extern const char *ILLmps_possibly_blank_name (
	const char *field,
	ILLread_mps_state * state,
	ILLsymboltab * tab);
extern int ILLmps_empty_key (
	ILLread_mps_state * state);
extern int ILLmps_empty_field (
	ILLread_mps_state * state);

extern int ILLmps_error (
	ILLread_mps_state * state,
	const char *format,
	...);
extern void ILLmps_warn (
	ILLread_mps_state * state,
	const char *format,
	...);

#endif
