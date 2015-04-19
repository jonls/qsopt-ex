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

/* RCSINFO $Id: lpdata_EGLPNUM_TYPENAME.h,v 1.4 2003/11/05 17:00:56 meven Exp $ */
#ifndef EGLPNUM_TYPENAME_ILL_LPDATA_H
#define EGLPNUM_TYPENAME_ILL_LPDATA_H

#include "eg_lpnum.h"
#include "reporter.h"
#include "symtab.h"

#include "qstruct_EGLPNUM_TYPENAME.h"
#include "readline_EGLPNUM_TYPENAME.h"
#include "format_EGLPNUM_TYPENAME.h"
#include "dstruct_EGLPNUM_TYPENAME.h"

extern EGLPNUM_TYPE EGLPNUM_TYPENAME_ILL_MAXDOUBLE;	/*  1e150 */
extern EGLPNUM_TYPE EGLPNUM_TYPENAME_ILL_MINDOUBLE;	/* -1e150 */

#define EGLPNUM_TYPENAME_ILL_MAXINT    (2147483647)	/* this is equal to 2^31-1 */
#define EGLPNUM_TYPENAME_ILL_MIN       (1)				/* Must be same as QS_MIN */
#define EGLPNUM_TYPENAME_ILL_MAX       (-1)			/* Must be same as QS_MAX */

/*  Setting Alg in Presolve  */

#define EGLPNUM_TYPENAME_ILL_PRE_SCALE           1
#define EGLPNUM_TYPENAME_ILL_PRE_FIXED           2
#define EGLPNUM_TYPENAME_ILL_PRE_SINGLE_ROW      4
#define EGLPNUM_TYPENAME_ILL_PRE_FORCING         8
#define EGLPNUM_TYPENAME_ILL_PRE_SINGLE_COL     16
#define EGLPNUM_TYPENAME_ILL_PRE_DUPLICATE_ROW  32
#define EGLPNUM_TYPENAME_ILL_PRE_DUPLICATE_COL  64
#define EGLPNUM_TYPENAME_ILL_PRE_EMPTY_COL     128
#define EGLPNUM_TYPENAME_ILL_PRE_ALL (EGLPNUM_TYPENAME_ILL_PRE_SCALE | EGLPNUM_TYPENAME_ILL_PRE_FIXED | EGLPNUM_TYPENAME_ILL_PRE_SINGLE_ROW           \
                    EGLPNUM_TYPENAME_ILL_PRE_FORCING | EGLPNUM_TYPENAME_ILL_PRE_SINGLE_COL | EGLPNUM_TYPENAME_ILL_PRE_DUPLICATE_ROW \
                   EGLPNUM_TYPENAME_ILL_PRE_DUPLICATE_COL | EGLPNUM_TYPENAME_ILL_PRE_EMPTY_COL)
#define EGLPNUM_TYPENAME_ILL_PRE_SIMPLE (EGLPNUM_TYPENAME_ILL_PRE_FIXED | EGLPNUM_TYPENAME_ILL_PRE_EMPTY_COL)

typedef struct EGLPNUM_TYPENAME_ILLlpdata
{																/* Complete LP data filled in by mpsread.  */
	int nrows;
	int ncols;
	int nstruct;									/* Not including logicals.                 */
	int nzcount;
	int rowsize;									/* Length of row arrays.                   */
	int colsize;									/* Length of col arrays.                   */
	int structsize;								/* Length of intmarker, structmap,         */
	/* colnames                                */
	int objsense;
	char *sense;									/* Original sense, not after logicals.     */
	EGLPNUM_TYPE *obj;
	EGLPNUM_TYPE *rhs;
	EGLPNUM_TYPE *rangeval;
	EGLPNUM_TYPE *lower;
	EGLPNUM_TYPE *upper;
	EGLPNUM_TYPENAME_ILLmatrix A;									/* The coef matrix.                        */
	struct EGLPNUM_TYPENAME_ILLlp_rows *rA;				/* Coef matrix in row form.                */

	char **rownames;
	ILLsymboltab rowtab;					/* contains rownames in no particular order */
	char *objname;								/* if colname is not NULL it is entered into 
																 * the rowtab, see reader fcts in lp.c, mps.c*/

	char **colnames;							/* columns of struct variables */
	ILLsymboltab coltab;					/* contains colnames in no particular order */

	char *probname;
	char *intmarker;
	int *structmap;								/* Indices of structural variables         */
	int *rowmap;									/* Indices of logical and range variables  */
	struct EGLPNUM_TYPENAME_ILLlp_basis *basis;
	struct EGLPNUM_TYPENAME_ILLlp_predata *presolve;
	struct EGLPNUM_TYPENAME_ILLlp_sinfo *sinfo;

	 /**************************************************************************/
	/* these fields are currently only set by mps.c reader fcts               */
	 /**************************************************************************/
	EGLPNUM_TYPENAME_ILLmatrix sos;								/* columns are the sets, rows are the  
																 * problem's structural variables
																 * coefficients are the weights */

	char *sos_type;								/* type of each set */
	int *is_sos_mem;							/* for each structural variable contains 
																 *    -1 == not a set member
																 *     i == member of sos set i 
																 *          where 0 <= i < sos.matcols */
	char *refrowname;							/* name of reference row */
	int refind;										/* index of reference row 
																 *     -1 if refrow was a free row 
																 *          and weights are found only in the 
																 *          sos matrix 
																 *     index >=0 if refrow is also a lp-row */

	 /**************************************************************************
    * EGLPNUM_TYPENAME_QSset_reporter initializes reporter 
    **************************************************************************/
	qsstring_reporter reporter;		/* used from within ILL fcts 
																 * to report feedback */
}
EGLPNUM_TYPENAME_ILLlpdata;

typedef struct EGLPNUM_TYPENAME_ILLlp_basis
{
	int nstruct;
	int nrows;
	int rownorms_size;
	int colnorms_size;
	char *cstat;
	char *rstat;
	EGLPNUM_TYPE *rownorms;
	EGLPNUM_TYPE *colnorms;
}
EGLPNUM_TYPENAME_ILLlp_basis;

typedef struct EGLPNUM_TYPENAME_ILLlp_cache
{
	int nstruct;
	int nrows;
	int status;
	EGLPNUM_TYPE val;
	EGLPNUM_TYPE *x;
	EGLPNUM_TYPE *pi;
	EGLPNUM_TYPE *rc;
	EGLPNUM_TYPE *slack;
}
EGLPNUM_TYPENAME_ILLlp_cache;

typedef struct EGLPNUM_TYPENAME_ILLlp_sinfo
{																/* LP info returned by presolve            */
	int ncols;
	int nrows;
	int nzcount;
	int rowsize;
	int colsize;
	int objsense;

	EGLPNUM_TYPE *obj;
	EGLPNUM_TYPE *rhs;
	EGLPNUM_TYPE *lower;
	EGLPNUM_TYPE *upper;

	EGLPNUM_TYPENAME_ILLmatrix A;

	char **colnames;							/* Just for debugging - not updated */
}
EGLPNUM_TYPENAME_ILLlp_sinfo;

typedef struct EGLPNUM_TYPENAME_ILLlp_preline
{
	EGLPNUM_TYPE rhs;
	EGLPNUM_TYPE obj;
	EGLPNUM_TYPE lower;
	EGLPNUM_TYPE upper;
	int count;
	int *ind;
	int row_or_col;								/* 0 is row, 1 is col */
	EGLPNUM_TYPE *val;
}
EGLPNUM_TYPENAME_ILLlp_preline;

typedef struct EGLPNUM_TYPENAME_ILLlp_preop
{
	int ptype;
	int rowindex;
	int colindex;
	EGLPNUM_TYPENAME_ILLlp_preline line;
}
EGLPNUM_TYPENAME_ILLlp_preop;

typedef struct EGLPNUM_TYPENAME_ILLlp_predata
{																/* Data needed in un-presolve.            */
	int opcount;
	int opsize;
	EGLPNUM_TYPENAME_ILLlp_preop *oplist;
	int r_nrows;
	int r_ncols;
	int *colmap;
	int *rowmap;
	EGLPNUM_TYPE *rowscale;
	EGLPNUM_TYPE *colscale;
	EGLPNUM_TYPE *colfixval;
	EGLPNUM_TYPE *rowfixval;
}
EGLPNUM_TYPENAME_ILLlp_predata;

typedef struct EGLPNUM_TYPENAME_ILLlp_rows
{
	int *rowbeg;
	int *rowcnt;
	int *rowind;
	EGLPNUM_TYPE *rowval;
}
EGLPNUM_TYPENAME_ILLlp_rows;


/****************************************************************************/
/*                                                                          */
/*                             lpdata.c                                     */
/*                                                                          */
/****************************************************************************/

struct EGLPNUM_TYPENAME_qsdata *EGLPNUM_TYPENAME_ILLread (
	EGLPNUM_TYPENAME_qsline_reader * file,
	const char *fname,
	int isMps);
void EGLPNUM_TYPENAME_ILLstart (
	void);							/**< initialize EGLPNUM_TYPENAME_ILL_MAXDOUBLE and other 

													 constants, this funtion should be callef AFTER 
													 EGlpNumStart() */
void EGLPNUM_TYPENAME_ILLend (
	void);						/**< free any internal data asociated with variable 

												 precision numbers */
void EGLPNUM_TYPENAME_ILLchange_precision (
	void);												/**< This function re-compute the internal 

																		 variables precision to the (previously 
																		 set) EGLPNUM_PRECISION value (done with 
																		 EGlpNumSetPrecision) */
void EGLPNUM_TYPENAME_ILLlpdata_init (
	EGLPNUM_TYPENAME_ILLlpdata * lp);
void EGLPNUM_TYPENAME_ILLlpdata_free (
	EGLPNUM_TYPENAME_ILLlpdata * lp);
void EGLPNUM_TYPENAME_ILLlp_basis_init (
	EGLPNUM_TYPENAME_ILLlp_basis * B);
void EGLPNUM_TYPENAME_ILLlp_basis_free (
	EGLPNUM_TYPENAME_ILLlp_basis * B);
void EGLPNUM_TYPENAME_ILLlp_cache_init (
	EGLPNUM_TYPENAME_ILLlp_cache * C);
void EGLPNUM_TYPENAME_ILLlp_cache_free (
	EGLPNUM_TYPENAME_ILLlp_cache * C);
int EGLPNUM_TYPENAME_ILLlp_basis_alloc (
	EGLPNUM_TYPENAME_ILLlp_basis * B,
	int ncols,
	int nrows);
int EGLPNUM_TYPENAME_ILLlp_cache_alloc (
	EGLPNUM_TYPENAME_ILLlp_cache * C,
	int ncols,
	int nrows);

int EGLPNUM_TYPENAME_ILLlp_rows_init (
	EGLPNUM_TYPENAME_ILLlp_rows * lp_rows,
	EGLPNUM_TYPENAME_ILLlpdata * lp,
	int include_logicals);
void EGLPNUM_TYPENAME_ILLlp_rows_clear (
	EGLPNUM_TYPENAME_ILLlp_rows * lp_rows);
int EGLPNUM_TYPENAME_ILLprint_report (
	EGLPNUM_TYPENAME_ILLlpdata * lp,
	const char *format,
	...);

							/* print to lp->reporter */

/****************************************************************************/
/*                                                                          */
/*                             presolve.c                                   */
/*                                                                          */
/****************************************************************************/

void EGLPNUM_TYPENAME_ILLlp_sinfo_init (
	EGLPNUM_TYPENAME_ILLlp_sinfo * sinfo),
  EGLPNUM_TYPENAME_ILLlp_sinfo_free (
	EGLPNUM_TYPENAME_ILLlp_sinfo * sinfo),
  EGLPNUM_TYPENAME_ILLlp_predata_init (
	EGLPNUM_TYPENAME_ILLlp_predata * pre),
  EGLPNUM_TYPENAME_ILLlp_predata_free (
	EGLPNUM_TYPENAME_ILLlp_predata * pre);

int EGLPNUM_TYPENAME_ILLlp_add_logicals (
	EGLPNUM_TYPENAME_ILLlpdata * lp),
  EGLPNUM_TYPENAME_ILLlp_scale (
	EGLPNUM_TYPENAME_ILLlpdata * lp),
  EGLPNUM_TYPENAME_ILLlp_presolve (
	EGLPNUM_TYPENAME_ILLlpdata * lp,
	int pre_types);

/* ========================================================================= */
/* if non-zero, then internal data has been initialized, and there is some
 * memory allocated, if zero, no internal memory has been allocated
 * (or it has been freed) */
extern int EGLPNUM_TYPENAME___QSEX_SETUP;

#endif /* __ILL_LPDATA_H */
