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

/* RCSINFO $Id: format_EGLPNUM_TYPENAME.h,v 1.3 2003/11/05 16:59:48 meven Exp $ */
#ifndef EGLPNUM_TYPENAME_QS_FORMAT_ERROR_H
#define EGLPNUM_TYPENAME_QS_FORMAT_ERROR_H

#include <stdio.h>

#include "qsopt_EGLPNUM_TYPENAME.h"
#include "eg_io.h"

/****************************************************************************/
/*
   The LP/MPS readers, writers, 
       EGLPNUM_TYPENAME_ILLrawlpdata_to_lpdata, and 
   use EGLPNUM_TYPENAME_ILLformat_error to report problems with their input iff
       the line reader used in reading the problem  or 
       the  EGLPNUM_TYPENAME_qserror_collector pointer passed to EGLPNUM_TYPENAME_ILLwrite_lp_file
   is not NULL.

   The QSgui code uses this feature to collect EGLPNUM_TYPENAME_qsformat_error instances 
   which it uses after reading is done to insert error messages into the 
   input window. 
*/
/****************************************************************************/

/* 
for error type USE: 
          QS_DATA_ERROR			
          QS_DATA_WARN			
          QS_MPS_FORMAT_ERROR		
          QS_MPS_FORMAT_WARN		
          QS_LP_FORMAT_ERROR		
          QS_LP_FORMAT_WARN		
          QS_LP_OBJ_WARN			
          QS_GENERIC_ERROR		
*/

typedef struct EGLPNUM_TYPENAME_qsformat_error
{
	char *desc;
	char *theLine;
	struct EGLPNUM_TYPENAME_qsformat_error *next;
	int type;
	int lineNumber;								/* 1 based line counting */
	int at;
}
EGLPNUM_TYPENAME_qsformat_error;

extern int EGLPNUM_TYPENAME_ILLformat_error_create (
	EGLPNUM_TYPENAME_qsformat_error * error,
	int mode,
	const char *desc,
	int lineNum,
	const char *theLine,
	int atPos);
extern void EGLPNUM_TYPENAME_ILLformat_error_delete (
	EGLPNUM_TYPENAME_qsformat_error * error);

extern void EGLPNUM_TYPENAME_ILLformat_error_print (
	EGioFile_t * out,
	EGLPNUM_TYPENAME_qsformat_error * e);



/*****************************************************************************
 * collecting error messages 
 * either with defining own qsad_error_fct and corresponding data structure 
 * or by using predefined EGLPNUM_TYPENAME_ILLadd_error_to_memory fct with EGLPNUM_TYPENAME_qserror_memory
 */

typedef int (
	*EGLPNUM_TYPENAME_qsadd_error_fct) (
	void *dest,
	const EGLPNUM_TYPENAME_qsformat_error * error);

typedef struct EGLPNUM_TYPENAME_qserror_collector
{
	EGLPNUM_TYPENAME_qsadd_error_fct add_error;
	void *dest;
}
EGLPNUM_TYPENAME_qserror_collector;

typedef struct EGLPNUM_TYPENAME_qserror_memory
{
	unsigned int nerror;
	EGLPNUM_TYPENAME_qsformat_error *error_list;
	char has_error[QS_INPUT_NERROR];
	char hasErrorLines;
}
EGLPNUM_TYPENAME_qserror_memory;


extern EGLPNUM_TYPENAME_qserror_collector *EGLPNUM_TYPENAME_ILLerror_collector_new (
	EGLPNUM_TYPENAME_qsadd_error_fct fct,
	void *dest);

EGLPNUM_TYPENAME_qserror_collector *EGLPNUM_TYPENAME_ILLerror_memory_collector_new (
	EGLPNUM_TYPENAME_qserror_memory * dest);

extern void EGLPNUM_TYPENAME_ILLerror_collector_free (
	EGLPNUM_TYPENAME_qserror_collector * c);

#define EGLPNUM_TYPENAME_ILLformat_error(collector, error)  \
	((collector)->add_error((collector)->dest, error))


extern int EGLPNUM_TYPENAME_ILLadd_error_to_memory (
	void *dest,
	const EGLPNUM_TYPENAME_qsformat_error * error);

extern EGLPNUM_TYPENAME_qserror_memory *EGLPNUM_TYPENAME_ILLerror_memory_create (
	int takeErrorLines);
extern void EGLPNUM_TYPENAME_ILLerror_memory_free (
	EGLPNUM_TYPENAME_qserror_memory * mem);

#endif
