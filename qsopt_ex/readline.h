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

/* RCSINFO $Id: rdline.h,v 1.2 2003/11/05 16:57:39 meven Exp $ */
#ifndef EGLPNUM_TYPENAME_LINE_READER_FILE_H
#define EGLPNUM_TYPENAME_LINE_READER_FILE_H

#include "qsopt_EGLPNUM_TYPENAME.h"
#include "format_EGLPNUM_TYPENAME.h"

/* #ifdef _WINDOWS */
typedef char *(
	*EGLPNUM_TYPENAME_qsread_line_fct) (
	char *s,
	int size,
	void *src);

typedef struct EGLPNUM_TYPENAME_qsline_reader
{
	EGLPNUM_TYPENAME_qsread_line_fct read_line_fct;
	void *data_src;
	struct EGLPNUM_TYPENAME_qserror_collector *error_collector;
}
EGLPNUM_TYPENAME_qsline_reader;

EGLPNUM_TYPENAME_qsline_reader *EGLPNUM_TYPENAME_ILLline_reader_new (
	EGLPNUM_TYPENAME_qsread_line_fct fct,
	void *data_src);
void EGLPNUM_TYPENAME_ILLline_reader_free (
	EGLPNUM_TYPENAME_qsline_reader * reader);

#define EGLPNUM_TYPENAME_ILLline_reader_get(s, size, reader)  \
	(reader)->read_line_fct(s, size, (reader)->data_src)
									 /* used by parsers to retrieve next input line */
/* #else  
 *
 * typedef FILE EGLPNUM_TYPENAME_qsline_reader; 
 *
 * #define EGLPNUM_TYPENAME_ILLline_reader_new(fct, data) ((FILE*) (data)) 
 * #define EGLPNUM_TYPENAME_ILLline_reader_free(reader)  
 * #define EGLPNUM_TYPENAME_ILLline_reader_get(s, size, reader) fgets(s,size,reader) 
 * #endif 
 */

#endif
