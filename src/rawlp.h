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

/*  RCS_INFO = "$RCSfile: rawlp.h,v $ $Revision: 1.3 $ $Date: 2003/11/05 16:57:39 $"; */
#ifndef __ILL_RAWLP_H_
#define __ILL_RAWLP_H_

/****************************************************************************/
/* DataStructure and Routines                                               */
/*          to deal with raw lp information as read from mps or lp files    */
/*          support scanning of input                                       */
/*          error reporting                                                 */
/****************************************************************************/

#include "trace.h"
#include "lpdata.h"
#include "iqsutil.h"
#include "format.h"
#include "lpdefs.h"

#define ILL_ISBLANK(p) \
             (((*(p))==' '||(*(p))=='\t'||(*(p))=='\r'||(*(p))=='\f') ? 1 : 0)

/* 
 * we rely on ILLsymboltab property:
 *   the ith name added can be retrieved by ILLsymboltab_get(table, i) 
 *   as long as we never delete names from the symbol table 
 */
typedef struct rawlpdata
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
	EGlpNum_t *rhs;								/* rhs values for rows; size is nrows */
	char *rangesind;							/* ranges[i] == 1 we saw a range def for row[i] */
	struct colptr *ranges;				/* list of range values */

	int ncols;										/* number of cols in problem */
	ILLsymboltab coltab;					/* ILLsymboltab_get(coltab, i) name of ith col */
	int colsize;									/* size of cols array */
	struct colptr **cols;

	char *lbind;									/* lbind[i] == 1  we saw a lower bound for col[i] */
	char *ubind;									/* ubind[i] == 1  we saw a upper bound for col[i] */
	EGlpNum_t *lower;							/* lower[i] = lower bound for col[i] */
	EGlpNum_t *upper;							/* upper[i] = upper bound for col[i] */

	int intsize;									/* size of intmarker array */
	char *intmarker;							/* intmarker[i] == 1  col[i] is an int var */

	/* sos information is tranfered into ILLmatrix lpdata->sos */
	char *refrow;									/* name of reference row */
	int refrowind;								/* index of refrow or -1  */

	int is_sos_size;							/* size of is_sos_member array */
	int *is_sos_member;						/* for each col contains either               
																 *     -1 == no sos memeber 
																 *     i  == member of set #i */

	int nsos_member;							/* total number of sos set members */
	int sos_weight_size;					/* size of sos_weight array */
	EGlpNum_t *sos_weight;				/* sos set elem i has weight of sos_weight[i] 
																 * value comes from refrow coeficients */
	int sos_col_size;							/* size of sos_col array */
	int *sos_col;									/* sos elem i is column sos_col[i] */

	int nsos;											/* number of sos sets */
	int sos_setsize;							/* size of sosset array */
	struct sosptr *sos_set;				/* type, size, first element of sos sets 
																 * first is index into sos_weight and sos_col 
																 * arrays */
	qserror_collector *error_collector;
	ILLptrworld ptrworld;
}
rawlpdata;

typedef struct colptr
{
	EGlpNum_t coef;
	struct colptr *next;
	int this_val;											/* row index */
}
colptr;
extern colptr *ILLcolptralloc (
	ILLptrworld * p);

typedef struct sosptr
{
	int nelem;										/* number of set elements */
	int first;										/* index of first set element in sosmemeber */
	char type;										/* set type */
}
sosptr;
extern const int ILL_SOS_TYPE1;
extern const int ILL_SOS_TYPE2;

extern void ILLinit_rawlpdata (
	rawlpdata * lp,
	qserror_collector * collector);
extern void ILLfree_rawlpdata (
	rawlpdata * lp);
extern void ILLraw_clear_matrix (
	rawlpdata * lp);

extern const char *ILLraw_rowname (
	rawlpdata * lp,
	int i);
extern const char *ILLraw_colname (
	rawlpdata * lp,
	int i);

extern int ILLraw_add_col (
	rawlpdata * lp,
	const char *name,
	int intmarker);
extern int ILLraw_add_row (
	rawlpdata * lp,
	const char *name,
	int sense,
	const EGlpNum_t rhs);

extern int ILLraw_add_col_coef (
	rawlpdata * lp,
	int colind,
	int rowind,
	EGlpNum_t coef);

extern int ILLraw_init_ranges (
	rawlpdata * lp);
extern int ILLraw_init_rhs (
	rawlpdata * lp);

extern int ILLraw_add_ranges_coef (
	rawlpdata * lp,
	int rowind,
	EGlpNum_t coef);


extern int ILLraw_add_sos (
	rawlpdata * lp,
	int sos_type);

																								/* add empty set with type */
extern int ILLraw_add_sos_member (
	rawlpdata * lp,
	int colind);

																								/* add col to last set */
extern int ILLraw_is_mem_other_sos (
	rawlpdata * lp,
	int colind);

extern int ILLraw_set_rhs_name (
	rawlpdata * lp,
	const char *name,
	int *skip);
extern int ILLraw_set_bounds_name (
	rawlpdata * lp,
	const char *name,
	int *skip);
extern int ILLraw_set_ranges_name (
	rawlpdata * lp,
	const char *name,
	int *skip);
extern void ILLprint_rawlpdata (
	rawlpdata * lp);

extern char *ILLraw_unique_name (
	ILLsymboltab * tab,
	char *prefix,
	int i);
extern int ILLraw_fill_in_rownames (
	rawlpdata * lp);

extern int ILLraw_init_bounds (
	rawlpdata * lp);

extern const char *ILLraw_set_lowerBound (
	rawlpdata * lp,
	int i,
	EGlpNum_t bnd);
extern const char *ILLraw_set_upperBound (
	rawlpdata * lp,
	int i,
	EGlpNum_t bnd);
extern const char *ILLraw_set_fixedBound (
	rawlpdata * lp,
	int i,
	EGlpNum_t bnd);
extern const char *ILLraw_set_binaryBound (
	rawlpdata * lp,
	int i);
extern const char *ILLraw_set_unbound (
	rawlpdata * lp,
	int colind);
extern int ILLraw_fill_in_bounds (
	rawlpdata * lp);

extern int ILLraw_first_nondefault_bound (
	ILLlpdata * lp);
extern int ILLraw_default_lower (
	ILLlpdata * lp,
	int i);
extern int ILLraw_default_upper (
	ILLlpdata * lp,
	int i,
	int ri);

extern int ILLrawlpdata_to_lpdata (
	rawlpdata * raw,
	ILLlpdata * lp);

extern int ILLdata_error (
	qserror_collector * collector,
	const char *format,
	...);
extern void ILLdata_warn (
	qserror_collector * collector,
	const char *format,
	...);

#endif
