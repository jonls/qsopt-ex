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

/* RCSINFO $Id: mps_EGLPNUM_TYPENAME.h,v 1.2 2003/11/05 16:57:39 meven Exp $ */
#ifndef EGLPNUM_TYPENAME_MPS_H
#define EGLPNUM_TYPENAME_MPS_H

#include "readline_EGLPNUM_TYPENAME.h"
#include "format_EGLPNUM_TYPENAME.h"

/****************************************************************************/
/*                                                                          */
/*              Routines to support Reading and Writing MPS Files           */
/*                                                                          */
/****************************************************************************/
#include "basicdefs.h"
extern const char *EGLPNUM_TYPENAME_ILLmps_section_name[ILL_MPS_N_SECTIONS + 2];

#include "lpdata_EGLPNUM_TYPENAME.h"
#include "rawlp_EGLPNUM_TYPENAME.h"
#include "read_mps_EGLPNUM_TYPENAME.h"

extern int EGLPNUM_TYPENAME_ILLread_mps (
	EGLPNUM_TYPENAME_qsline_reader * file,
	const char *filename,
	EGLPNUM_TYPENAME_rawlpdata * lp);

extern int EGLPNUM_TYPENAME_ILLwrite_mps (
	EGLPNUM_TYPENAME_ILLlpdata * lp,
	EGLPNUM_TYPENAME_qserror_collector * collector);

				/* use lp->reporter for output */

#endif
