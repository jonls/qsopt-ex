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

/* RCSINFO $Id: lpdata.h,v 1.4 2003/11/05 17:00:56 meven Exp $ */
#ifndef ILL_LPDATA_H
#define ILL_LPDATA_H

#include "qstruct.h"
#include "iqsutil.h"
#include "readline.h"
#include "reporter.h"
#include "format.h"
#include "dstruct.h"

extern EGlpNum_t ILL_MAXDOUBLE;	/*  1e150 */
extern EGlpNum_t ILL_MINDOUBLE;	/* -1e150 */

#define ILL_MAXINT    (2147483647)	/* this is equal to 2^31-1 */
#define ILL_MIN       (1)				/* Must be same as QS_MIN */
#define ILL_MAX       (-1)			/* Must be same as QS_MAX */

/*  Setting Alg in Presolve  */

#define ILL_PRE_SCALE           1
#define ILL_PRE_FIXED           2
#define ILL_PRE_SINGLE_ROW      4
#define ILL_PRE_FORCING         8
#define ILL_PRE_SINGLE_COL     16
#define ILL_PRE_DUPLICATE_ROW  32
#define ILL_PRE_DUPLICATE_COL  64
#define ILL_PRE_EMPTY_COL     128
#define ILL_PRE_ALL (ILL_PRE_SCALE | ILL_PRE_FIXED | ILL_PRE_SINGLE_ROW           \
                    ILL_PRE_FORCING | ILL_PRE_SINGLE_COL | ILL_PRE_DUPLICATE_ROW \
                   ILL_PRE_DUPLICATE_COL | ILL_PRE_EMPTY_COL)
#define ILL_PRE_SIMPLE (ILL_PRE_FIXED | ILL_PRE_EMPTY_COL)

typedef struct ILLlpdata
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
	EGlpNum_t *obj;
	EGlpNum_t *rhs;
	EGlpNum_t *rangeval;
	EGlpNum_t *lower;
	EGlpNum_t *upper;
	ILLmatrix A;									/* The coef matrix.                        */
	struct ILLlp_rows *rA;				/* Coef matrix in row form.                */

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
	struct ILLlp_basis *basis;
	struct ILLlp_predata *presolve;
	struct ILLlp_sinfo *sinfo;

	 /**************************************************************************/
	/* these fields are currently only set by mps.c reader fcts               */
	 /**************************************************************************/
	ILLmatrix sos;								/* columns are the sets, rows are the  
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
    * QSset_reporter initializes reporter 
    **************************************************************************/
	qsstring_reporter reporter;		/* used from within ILL fcts 
																 * to report feedback */
}
ILLlpdata;

typedef struct ILLlp_basis
{
	int nstruct;
	int nrows;
	int rownorms_size;
	int colnorms_size;
	char *cstat;
	char *rstat;
	EGlpNum_t *rownorms;
	EGlpNum_t *colnorms;
}
ILLlp_basis;

typedef struct ILLlp_cache
{
	int nstruct;
	int nrows;
	int status;
	EGlpNum_t val;
	EGlpNum_t *x;
	EGlpNum_t *pi;
	EGlpNum_t *rc;
	EGlpNum_t *slack;
}
ILLlp_cache;

typedef struct ILLlp_sinfo
{																/* LP info returned by presolve            */
	int ncols;
	int nrows;
	int nzcount;
	int rowsize;
	int colsize;
	int objsense;

	EGlpNum_t *obj;
	EGlpNum_t *rhs;
	EGlpNum_t *lower;
	EGlpNum_t *upper;

	ILLmatrix A;

	char **colnames;							/* Just for debugging - not updated */
}
ILLlp_sinfo;

typedef struct ILLlp_preline
{
	EGlpNum_t rhs;
	EGlpNum_t obj;
	EGlpNum_t lower;
	EGlpNum_t upper;
	int count;
	int *ind;
	int row_or_col;								/* 0 is row, 1 is col */
	EGlpNum_t *val;
}
ILLlp_preline;

typedef struct ILLlp_preop
{
	int ptype;
	int rowindex;
	int colindex;
	ILLlp_preline line;
}
ILLlp_preop;

typedef struct ILLlp_predata
{																/* Data needed in un-presolve.            */
	int opcount;
	int opsize;
	ILLlp_preop *oplist;
	int r_nrows;
	int r_ncols;
	int *colmap;
	int *rowmap;
	EGlpNum_t *rowscale;
	EGlpNum_t *colscale;
	EGlpNum_t *colfixval;
	EGlpNum_t *rowfixval;
}
ILLlp_predata;

typedef struct ILLlp_rows
{
	int *rowbeg;
	int *rowcnt;
	int *rowind;
	EGlpNum_t *rowval;
}
ILLlp_rows;


/****************************************************************************/
/*                                                                          */
/*                             lpdata.c                                     */
/*                                                                          */
/****************************************************************************/

struct qsdata *ILLread (
	qsline_reader * file,
	const char *fname,
	int isMps);
void ILLstart (
	void);							/**< initialize ILL_MAXDOUBLE and other 

													 constants, this funtion should be callef AFTER 
													 EGlpNumStart() */
void ILLend (
	void);						/**< free any internal data asociated with variable 

												 precision numbers */
void ILLchange_precision (
	void);												/**< This function re-compute the internal 

																		 variables precision to the (previously 
																		 set) EGLPNUM_PRECISION value (done with 
																		 EGlpNumSetPrecision) */
void ILLlpdata_init (
	ILLlpdata * lp);
void ILLlpdata_free (
	ILLlpdata * lp);
void ILLlp_basis_init (
	ILLlp_basis * B);
void ILLlp_basis_free (
	ILLlp_basis * B);
void ILLlp_cache_init (
	ILLlp_cache * C);
void ILLlp_cache_free (
	ILLlp_cache * C);
int ILLlp_basis_alloc (
	ILLlp_basis * B,
	int ncols,
	int nrows);
int ILLlp_cache_alloc (
	ILLlp_cache * C,
	int ncols,
	int nrows);

int ILLlp_rows_init (
	ILLlp_rows * lp_rows,
	ILLlpdata * lp,
	int include_logicals);
void ILLlp_rows_clear (
	ILLlp_rows * lp_rows);
int ILLprint_report (
	ILLlpdata * lp,
	const char *format,
	...);

							/* print to lp->reporter */

/****************************************************************************/
/*                                                                          */
/*                             presolve.c                                   */
/*                                                                          */
/****************************************************************************/

void ILLlp_sinfo_init (
	ILLlp_sinfo * sinfo),
  ILLlp_sinfo_free (
	ILLlp_sinfo * sinfo),
  ILLlp_predata_init (
	ILLlp_predata * pre),
  ILLlp_predata_free (
	ILLlp_predata * pre);

int ILLlp_add_logicals (
	ILLlpdata * lp),
  ILLlp_scale (
	ILLlpdata * lp),
  ILLlp_presolve (
	ILLlpdata * lp,
	int pre_types);

/* ========================================================================= */
/* if non-zero, then internal data has been initialized, and there is some
 * memory allocated, if zero, no internal memory has been allocated
 * (or it has been freed) */
extern int __QSEX_SETUP;

#endif /* __ILL_LPDATA_H */
