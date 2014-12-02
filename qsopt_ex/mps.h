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

/* RCSINFO $Id: mps.h,v 1.2 2003/11/05 16:57:39 meven Exp $ */
#ifndef MPS_H
#define MPS_H

#include "readline.h"
#include "format.h"

/****************************************************************************/
/*                                                                          */
/*              Routines to support Reading and Writing MPS Files           */
/*                                                                          */
/****************************************************************************/
#include "basicdefs.h"
extern const char *ILLmps_section_name[ILL_MPS_N_SECTIONS + 2];

#include "lpdata.h"
#include "rawlp.h"
#include "read_mps.h"

extern int ILLread_mps (
	qsline_reader * file,
	const char *filename,
	rawlpdata * lp);

extern int ILLwrite_mps (
	ILLlpdata * lp,
	qserror_collector * collector);

				/* use lp->reporter for output */

#endif
