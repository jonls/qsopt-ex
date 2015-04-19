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

/*  RCS_INFO = "$RCSfile: rawlp_EGLPNUM_TYPENAME.h,v $ $Revision: 1.3 $ $Date: 2003/11/05 16:57:39 $"; */
#ifndef EGLPNUM_TYPENAME___ILL_RAWLP_H_
#define EGLPNUM_TYPENAME___ILL_RAWLP_H_

/****************************************************************************/
/* DataStructure and Routines                                               */
/*          to deal with raw lp information as read from mps or lp files    */
/*          support scanning of input                                       */
/*          error reporting                                                 */
/****************************************************************************/

#include "allocrus.h"
#include "eg_lpnum.h"
#include "symtab.h"
#include "trace.h"

#include "lpdata_EGLPNUM_TYPENAME.h"
#include "format_EGLPNUM_TYPENAME.h"
#include "lpdefs_EGLPNUM_TYPENAME.h"

#define EGLPNUM_TYPENAME_ILL_ISBLANK(p) \
             (((*(p))==' '||(*(p))=='\t'||(*(p))=='\r'||(*(p))=='\f') ? 1 : 0)

/* 
 * we rely on ILLsymboltab property:
 *   the ith name added can be retrieved by ILLsymboltab_get(table, i) 
 *   as long as we never delete names from the symbol table 
 */
typedef struct EGLPNUM_TYPENAME_rawlpdata
{
	char *name;

	char *rhsname;
	char *rangesname;
	char *boundsname;

	int objsense;									/* maximize or minimize */
	int objindex;									/* index of objective row */

	int nrows;										/* number of rows in problem */
	ILLsymboltab rowtab;					/* ILLsymboltab_get(rowtab, i) name of ith row */

	int sensesize;								/* size of rowsense */
	char *rowsense;								/* rowsense[i] snese of row[i] */

	char *rhsind;									/* rhsind[i] == 1 we saw an rhs for row[i] */
	/* size is nrows */
	int rhssize;									/* size of rhs array */
	EGLPNUM_TYPE *rhs;								/* rhs values for rows; size is nrows */
	char *rangesind;							/* ranges[i] == 1 we saw a range def for row[i] */
	struct EGLPNUM_TYPENAME_colptr *ranges;				/* list of range values */

	int ncols;										/* number of cols in problem */
	ILLsymboltab coltab;					/* ILLsymboltab_get(coltab, i) name of ith col */
	int colsize;									/* size of cols array */
	struct EGLPNUM_TYPENAME_colptr **cols;

	char *lbind;									/* lbind[i] == 1  we saw a lower bound for col[i] */
	char *ubind;									/* ubind[i] == 1  we saw a upper bound for col[i] */
	EGLPNUM_TYPE *lower;							/* lower[i] = lower bound for col[i] */
	EGLPNUM_TYPE *upper;							/* upper[i] = upper bound for col[i] */

	int intsize;									/* size of intmarker array */
	char *intmarker;							/* intmarker[i] == 1  col[i] is an int var */

	/* sos information is tranfered into EGLPNUM_TYPENAME_ILLmatrix lpdata->sos */
	char *refrow;									/* name of reference row */
	int refrowind;								/* index of refrow or -1  */

	int is_sos_size;							/* size of is_sos_member array */
	int *is_sos_member;						/* for each col contains either               
																 *     -1 == no sos memeber 
																 *     i  == member of set #i */

	int nsos_member;							/* total number of sos set members */
	int sos_weight_size;					/* size of sos_weight array */
	EGLPNUM_TYPE *sos_weight;				/* sos set elem i has weight of sos_weight[i] 
																 * value comes from refrow coeficients */
	int sos_col_size;							/* size of sos_col array */
	int *sos_col;									/* sos elem i is column sos_col[i] */

	int nsos;											/* number of sos sets */
	int sos_setsize;							/* size of sosset array */
	struct EGLPNUM_TYPENAME_sosptr *sos_set;				/* type, size, first element of sos sets 
																 * first is index into sos_weight and sos_col 
																 * arrays */
	EGLPNUM_TYPENAME_qserror_collector *error_collector;
	ILLptrworld ptrworld;
}
EGLPNUM_TYPENAME_rawlpdata;

typedef struct EGLPNUM_TYPENAME_colptr
{
	EGLPNUM_TYPE coef;
	struct EGLPNUM_TYPENAME_colptr *next;
	int this_val;											/* row index */
}
EGLPNUM_TYPENAME_colptr;
extern EGLPNUM_TYPENAME_colptr *EGLPNUM_TYPENAME_ILLcolptralloc (
	ILLptrworld * p);

typedef struct EGLPNUM_TYPENAME_sosptr
{
	int nelem;										/* number of set elements */
	int first;										/* index of first set element in sosmemeber */
	char type;										/* set type */
}
EGLPNUM_TYPENAME_sosptr;
extern const int EGLPNUM_TYPENAME_ILL_SOS_TYPE1;
extern const int EGLPNUM_TYPENAME_ILL_SOS_TYPE2;

extern void EGLPNUM_TYPENAME_ILLinit_rawlpdata (
	EGLPNUM_TYPENAME_rawlpdata * lp,
	EGLPNUM_TYPENAME_qserror_collector * collector);
extern void EGLPNUM_TYPENAME_ILLfree_rawlpdata (
	EGLPNUM_TYPENAME_rawlpdata * lp);
extern void EGLPNUM_TYPENAME_ILLraw_clear_matrix (
	EGLPNUM_TYPENAME_rawlpdata * lp);

extern const char *EGLPNUM_TYPENAME_ILLraw_rowname (
	EGLPNUM_TYPENAME_rawlpdata * lp,
	int i);
extern const char *EGLPNUM_TYPENAME_ILLraw_colname (
	EGLPNUM_TYPENAME_rawlpdata * lp,
	int i);

extern int EGLPNUM_TYPENAME_ILLraw_add_col (
	EGLPNUM_TYPENAME_rawlpdata * lp,
	const char *name,
	int intmarker);
extern int EGLPNUM_TYPENAME_ILLraw_add_row (
	EGLPNUM_TYPENAME_rawlpdata * lp,
	const char *name,
	int sense,
	const EGLPNUM_TYPE rhs);

extern int EGLPNUM_TYPENAME_ILLraw_add_col_coef (
	EGLPNUM_TYPENAME_rawlpdata * lp,
	int colind,
	int rowind,
	EGLPNUM_TYPE coef);

extern int EGLPNUM_TYPENAME_ILLraw_init_ranges (
	EGLPNUM_TYPENAME_rawlpdata * lp);
extern int EGLPNUM_TYPENAME_ILLraw_init_rhs (
	EGLPNUM_TYPENAME_rawlpdata * lp);

extern int EGLPNUM_TYPENAME_ILLraw_add_ranges_coef (
	EGLPNUM_TYPENAME_rawlpdata * lp,
	int rowind,
	EGLPNUM_TYPE coef);


extern int EGLPNUM_TYPENAME_ILLraw_add_sos (
	EGLPNUM_TYPENAME_rawlpdata * lp,
	int sos_type);

																								/* add empty set with type */
extern int EGLPNUM_TYPENAME_ILLraw_add_sos_member (
	EGLPNUM_TYPENAME_rawlpdata * lp,
	int colind);

																								/* add col to last set */
extern int EGLPNUM_TYPENAME_ILLraw_is_mem_other_sos (
	EGLPNUM_TYPENAME_rawlpdata * lp,
	int colind);

extern int EGLPNUM_TYPENAME_ILLraw_set_rhs_name (
	EGLPNUM_TYPENAME_rawlpdata * lp,
	const char *name,
	int *skip);
extern int EGLPNUM_TYPENAME_ILLraw_set_bounds_name (
	EGLPNUM_TYPENAME_rawlpdata * lp,
	const char *name,
	int *skip);
extern int EGLPNUM_TYPENAME_ILLraw_set_ranges_name (
	EGLPNUM_TYPENAME_rawlpdata * lp,
	const char *name,
	int *skip);
extern void EGLPNUM_TYPENAME_ILLprint_rawlpdata (
	EGLPNUM_TYPENAME_rawlpdata * lp);

extern char *EGLPNUM_TYPENAME_ILLraw_unique_name (
	ILLsymboltab * tab,
	char *prefix,
	int i);
extern int EGLPNUM_TYPENAME_ILLraw_fill_in_rownames (
	EGLPNUM_TYPENAME_rawlpdata * lp);

extern int EGLPNUM_TYPENAME_ILLraw_init_bounds (
	EGLPNUM_TYPENAME_rawlpdata * lp);

extern const char *EGLPNUM_TYPENAME_ILLraw_set_lowerBound (
	EGLPNUM_TYPENAME_rawlpdata * lp,
	int i,
	EGLPNUM_TYPE bnd);
extern const char *EGLPNUM_TYPENAME_ILLraw_set_upperBound (
	EGLPNUM_TYPENAME_rawlpdata * lp,
	int i,
	EGLPNUM_TYPE bnd);
extern const char *EGLPNUM_TYPENAME_ILLraw_set_fixedBound (
	EGLPNUM_TYPENAME_rawlpdata * lp,
	int i,
	EGLPNUM_TYPE bnd);
extern const char *EGLPNUM_TYPENAME_ILLraw_set_binaryBound (
	EGLPNUM_TYPENAME_rawlpdata * lp,
	int i);
extern const char *EGLPNUM_TYPENAME_ILLraw_set_unbound (
	EGLPNUM_TYPENAME_rawlpdata * lp,
	int colind);
extern int EGLPNUM_TYPENAME_ILLraw_fill_in_bounds (
	EGLPNUM_TYPENAME_rawlpdata * lp);

extern int EGLPNUM_TYPENAME_ILLraw_first_nondefault_bound (
	EGLPNUM_TYPENAME_ILLlpdata * lp);
extern int EGLPNUM_TYPENAME_ILLraw_default_lower (
	EGLPNUM_TYPENAME_ILLlpdata * lp,
	int i);
extern int EGLPNUM_TYPENAME_ILLraw_default_upper (
	EGLPNUM_TYPENAME_ILLlpdata * lp,
	int i,
	int ri);

extern int EGLPNUM_TYPENAME_ILLrawlpdata_to_lpdata (
	EGLPNUM_TYPENAME_rawlpdata * raw,
	EGLPNUM_TYPENAME_ILLlpdata * lp);

extern int EGLPNUM_TYPENAME_ILLdata_error (
	EGLPNUM_TYPENAME_qserror_collector * collector,
	const char *format,
	...);
extern void EGLPNUM_TYPENAME_ILLdata_warn (
	EGLPNUM_TYPENAME_qserror_collector * collector,
	const char *format,
	...);

#endif
