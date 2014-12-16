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

/*  $RCSfile: qsopt_EGLPNUM_TYPENAME.h,v $ $Revision: 1.3 $ $Date: 2003/11/05 16:57:39 $" */
#ifndef EGLPNUM_TYPENAME___QS_QSOPT_H
#define EGLPNUM_TYPENAME___QS_QSOPT_H

#include <stdlib.h>
#include <stdio.h>

#include <gmp.h>

#ifdef WIN32

#ifdef QSLIB_EXPORTS
#define EGLPNUM_TYPENAME_QSLIB_INTERFACE __declspec(dllexport)
#else
#define EGLPNUM_TYPENAME_QSLIB_INTERFACE __declspec(dllimport)
#endif

#else
#define EGLPNUM_TYPENAME_QSLIB_INTERFACE extern
#endif

#ifdef WIN32
typedef struct EGLPNUM_TYPENAME_QSLIB_INTERFACE EGLPNUM_TYPENAME_qsdata *EGLPNUM_TYPENAME_QSprob;
typedef struct EGLPNUM_TYPENAME_QSLIB_INTERFACE qsbasis *EGLPNUM_TYPENAME_QSbas;
#else
typedef struct EGLPNUM_TYPENAME_qsdata *EGLPNUM_TYPENAME_QSprob;
typedef struct qsbasis *EGLPNUM_TYPENAME_QSbas;
#endif

/****************************************************************************/
/*                                                                          */
/*                 PARAMETERS TO SPECIFY OBJECTIVE SENSE                    */
/*                                                                          */
/****************************************************************************/
#include "basicdefs.h"
/*
#define QS_LP_PRIMAL_FEASIBLE   11
#define QS_LP_PRIMAL_INFEASIBLE 12
#define QS_LP_PRIMAL_UNBOUNDED  13
#define QS_LP_DUAL_FEASIBLE     14
#define QS_LP_DUAL_INFEASIBLE   15
#define QS_LP_DUAL_UNBOUNDED    16
*/

/****************************************************************************/
/*                                                                          */
/*                      QSopt Library Functions                             */
/*                                                                          */
/****************************************************************************/
#ifdef  __cplusplus
extern "C"
{
#endif

#ifdef WIN32
/* 
 *  in WINDOWS we make 
 *     EGLPNUM_TYPENAME_solver_main/EGLPNUM_TYPENAME_reader_main part of DLL
 */
EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_solver_main ( int argc, char **argv);
EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_reader_main ( int argc, char **argv);
#endif

EGLPNUM_TYPENAME_QSLIB_INTERFACE void EGLPNUM_TYPENAME_QSfree ( void *ptr),
		EGLPNUM_TYPENAME_QSfree_prob ( EGLPNUM_TYPENAME_QSprob p),
		EGLPNUM_TYPENAME_QSfree_basis ( EGLPNUM_TYPENAME_QSbas B),
	  EGLPNUM_TYPENAME_QSset_precision ( const unsigned prec),/**< set the precision for floating 
																								 point numbers to the given 
																								 number of bits */
		EGLPNUM_TYPENAME_QSstart ( void),/**< whe we use non native numbers, we need to make 
												 some initializations before operating with the
												 library */
	  EGLPNUM_TYPENAME_QSend ( void);	/**< just to free any internal static data needed by
												 the variable precision numbers */

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSopt_primal ( EGLPNUM_TYPENAME_QSprob p, int *status),
		EGLPNUM_TYPENAME_QSopt_dual ( EGLPNUM_TYPENAME_QSprob p, int *status),
		EGLPNUM_TYPENAME_QSopt_pivotin_col ( EGLPNUM_TYPENAME_QSprob p, int ccnt, int *clist),
		EGLPNUM_TYPENAME_QSopt_pivotin_row ( EGLPNUM_TYPENAME_QSprob p, int rcnt, int *rlist),
		EGLPNUM_TYPENAME_QSopt_strongbranch ( EGLPNUM_TYPENAME_QSprob p, int ncand, int *candidatelist,
			EGLPNUM_TYPE * xlist, EGLPNUM_TYPE * down_vals, EGLPNUM_TYPE * up_vals,
			int iterations, EGLPNUM_TYPE objbound),
		EGLPNUM_TYPENAME_QSchange_objsense ( EGLPNUM_TYPENAME_QSprob p, int newsense),
		EGLPNUM_TYPENAME_QSget_objsense ( EGLPNUM_TYPENAME_QSprob p, int *newsense),
		EGLPNUM_TYPENAME_QSnew_col ( EGLPNUM_TYPENAME_QSprob p,const EGLPNUM_TYPE obj,const EGLPNUM_TYPE lower,const EGLPNUM_TYPE upper,
			const char *name),
		EGLPNUM_TYPENAME_QSadd_cols ( EGLPNUM_TYPENAME_QSprob p, int num, int *cmatcnt, int *cmatbeg, int *cmatind,
			EGLPNUM_TYPE * cmatval, EGLPNUM_TYPE * obj, EGLPNUM_TYPE * lower,
			EGLPNUM_TYPE * upper, const char **names),
		EGLPNUM_TYPENAME_QSadd_col ( EGLPNUM_TYPENAME_QSprob p, int cnt, int *cmatind, EGLPNUM_TYPE * cmatval,
			EGLPNUM_TYPE obj, EGLPNUM_TYPE lower, EGLPNUM_TYPE upper, const char *name),
		EGLPNUM_TYPENAME_QSnew_row ( EGLPNUM_TYPENAME_QSprob p,const EGLPNUM_TYPE rhs, int sense, const char *name),
		EGLPNUM_TYPENAME_QSadd_ranged_rows ( EGLPNUM_TYPENAME_QSprob p, int num, int *rmatcnt, int *rmatbeg, 
			int *rmatind,const EGLPNUM_TYPE * rmatval,const EGLPNUM_TYPE * rhs, char *sense,
			const EGLPNUM_TYPE* range, const char **names),
		EGLPNUM_TYPENAME_QSadd_ranged_row ( EGLPNUM_TYPENAME_QSprob p, int cnt, int *rmatind,const EGLPNUM_TYPE * rmatval,
			const EGLPNUM_TYPE * rhs, int sense,const EGLPNUM_TYPE * range, const char *name),
		EGLPNUM_TYPENAME_QSadd_rows ( EGLPNUM_TYPENAME_QSprob p, int num, int *rmatcnt, int *rmatbeg, int *rmatind,
			const EGLPNUM_TYPE * rmatval,const EGLPNUM_TYPE * rhs, char *sense, const char **names),
		EGLPNUM_TYPENAME_QSadd_row ( EGLPNUM_TYPENAME_QSprob p, int cnt, int *rmatind,const EGLPNUM_TYPE * rmatval,
			const EGLPNUM_TYPE * rhs, int sense, const char *name),
		EGLPNUM_TYPENAME_QSdelete_rows ( EGLPNUM_TYPENAME_QSprob p, int num, int *dellist),
		EGLPNUM_TYPENAME_QSdelete_row ( EGLPNUM_TYPENAME_QSprob p, int rowindex),
		EGLPNUM_TYPENAME_QSdelete_setrows ( EGLPNUM_TYPENAME_QSprob p, int *flags),
		EGLPNUM_TYPENAME_QSdelete_named_row ( EGLPNUM_TYPENAME_QSprob p, const char *rowname),
		EGLPNUM_TYPENAME_QSdelete_named_rows_list ( EGLPNUM_TYPENAME_QSprob p, int num, const char **rownames),
		EGLPNUM_TYPENAME_QSdelete_cols ( EGLPNUM_TYPENAME_QSprob p, int num, int *dellist),
		EGLPNUM_TYPENAME_QSdelete_col ( EGLPNUM_TYPENAME_QSprob p, int colindex),
		EGLPNUM_TYPENAME_QSdelete_setcols ( EGLPNUM_TYPENAME_QSprob p, int *flags),
		EGLPNUM_TYPENAME_QSdelete_named_column ( EGLPNUM_TYPENAME_QSprob p, const char *colname),
		EGLPNUM_TYPENAME_QSdelete_named_columns_list ( EGLPNUM_TYPENAME_QSprob p, int num, const char **colnames),
		EGLPNUM_TYPENAME_QSchange_senses ( EGLPNUM_TYPENAME_QSprob p, int num, int *rowlist, char *sense),
		EGLPNUM_TYPENAME_QSchange_sense ( EGLPNUM_TYPENAME_QSprob p, int rowindex, int sense),
		EGLPNUM_TYPENAME_QSchange_coef ( EGLPNUM_TYPENAME_QSprob p, int rowindex, int colindex, EGLPNUM_TYPE coef),
		EGLPNUM_TYPENAME_QSchange_objcoef ( EGLPNUM_TYPENAME_QSprob p, int indx, EGLPNUM_TYPE coef),
		EGLPNUM_TYPENAME_QSchange_rhscoef ( EGLPNUM_TYPENAME_QSprob p, int indx, EGLPNUM_TYPE coef),
		EGLPNUM_TYPENAME_QSchange_range(EGLPNUM_TYPENAME_QSprob p, int rowindex, EGLPNUM_TYPE range),
		EGLPNUM_TYPENAME_QSchange_bounds ( EGLPNUM_TYPENAME_QSprob p, int num, int *collist, char *lu, 
			const EGLPNUM_TYPE * bounds),
		EGLPNUM_TYPENAME_QSchange_bound ( EGLPNUM_TYPENAME_QSprob p, int indx, int lu,const EGLPNUM_TYPE bound),
		EGLPNUM_TYPENAME_QSload_basis ( EGLPNUM_TYPENAME_QSprob p, EGLPNUM_TYPENAME_QSbas B),
		EGLPNUM_TYPENAME_QSread_and_load_basis ( EGLPNUM_TYPENAME_QSprob p, const char *filename),
		EGLPNUM_TYPENAME_QSload_basis_array ( EGLPNUM_TYPENAME_QSprob p, char *cstat, char *rstat),
		EGLPNUM_TYPENAME_QSload_basis_and_row_norms_array ( EGLPNUM_TYPENAME_QSprob p, char *cstat, char *rstat, 
			EGLPNUM_TYPE * rownorms),
		EGLPNUM_TYPENAME_QSget_basis_array ( EGLPNUM_TYPENAME_QSprob p, char *cstat, char *rstat),
		EGLPNUM_TYPENAME_QSget_basis_and_row_norms_array ( EGLPNUM_TYPENAME_QSprob p, char *cstat, char *rstat,
			EGLPNUM_TYPE * rownorms),
		EGLPNUM_TYPENAME_QSget_binv_row ( EGLPNUM_TYPENAME_QSprob p, int indx, EGLPNUM_TYPE * binvrow),
		EGLPNUM_TYPENAME_QSget_tableau_row ( EGLPNUM_TYPENAME_QSprob p, int indx, EGLPNUM_TYPE * tableaurow),
		EGLPNUM_TYPENAME_QSget_basis_order ( EGLPNUM_TYPENAME_QSprob p, int *basorder), 
    EGLPNUM_TYPENAME_QSget_coef (EGLPNUM_TYPENAME_QSprob p, int rowindex, int colindex, EGLPNUM_TYPE*coef),
		EGLPNUM_TYPENAME_QSget_status ( EGLPNUM_TYPENAME_QSprob p, int *status),
		EGLPNUM_TYPENAME_QSget_solution ( EGLPNUM_TYPENAME_QSprob p, EGLPNUM_TYPE * value, EGLPNUM_TYPE * x,
			EGLPNUM_TYPE * pi, EGLPNUM_TYPE * slack, EGLPNUM_TYPE * rc),
		EGLPNUM_TYPENAME_QSget_objval ( EGLPNUM_TYPENAME_QSprob p, EGLPNUM_TYPE * value),
		EGLPNUM_TYPENAME_QSget_pi_array ( EGLPNUM_TYPENAME_QSprob p, EGLPNUM_TYPE * pi),
		EGLPNUM_TYPENAME_QSget_rc_array ( EGLPNUM_TYPENAME_QSprob p, EGLPNUM_TYPE * rc),
		EGLPNUM_TYPENAME_QSget_x_array ( EGLPNUM_TYPENAME_QSprob p, EGLPNUM_TYPE * x),
		EGLPNUM_TYPENAME_QSget_slack_array ( EGLPNUM_TYPENAME_QSprob p, EGLPNUM_TYPE * slack),
		EGLPNUM_TYPENAME_QSget_infeas_array ( EGLPNUM_TYPENAME_QSprob p, EGLPNUM_TYPE * pi),
		EGLPNUM_TYPENAME_QSget_colcount ( EGLPNUM_TYPENAME_QSprob p),
		EGLPNUM_TYPENAME_QSget_rowcount ( EGLPNUM_TYPENAME_QSprob p),
		EGLPNUM_TYPENAME_QSget_nzcount ( EGLPNUM_TYPENAME_QSprob p),
		EGLPNUM_TYPENAME_QSget_obj_list(EGLPNUM_TYPENAME_QSprob p, int num, int*collist, EGLPNUM_TYPE*obj),
		EGLPNUM_TYPENAME_QSget_obj ( EGLPNUM_TYPENAME_QSprob p, EGLPNUM_TYPE * obj),
		EGLPNUM_TYPENAME_QSget_rhs ( EGLPNUM_TYPENAME_QSprob p, EGLPNUM_TYPE * rhs),
		EGLPNUM_TYPENAME_QSget_ranged_rows_list ( EGLPNUM_TYPENAME_QSprob p, int num, int *rowlist, int **rowcnt,
			int **rowbeg, int **rowind, EGLPNUM_TYPE ** rowval, EGLPNUM_TYPE ** rhs,
			char **sense, EGLPNUM_TYPE **range, char ***names),
		EGLPNUM_TYPENAME_QSget_ranged_rows ( EGLPNUM_TYPENAME_QSprob p, int **rowcnt, int **rowbeg, int **rowind,
			EGLPNUM_TYPE ** rowval, EGLPNUM_TYPE ** rhs, char **sense, 
			EGLPNUM_TYPE ** range, char ***names),
		EGLPNUM_TYPENAME_QSget_senses ( EGLPNUM_TYPENAME_QSprob p, char*senses),
		EGLPNUM_TYPENAME_QSget_rows_list ( EGLPNUM_TYPENAME_QSprob p, int num, int *rowlist, int **rowcnt,
			int **rowbeg, int **rowind, EGLPNUM_TYPE ** rowval, EGLPNUM_TYPE ** rhs,
			char **sense, char ***names),
		EGLPNUM_TYPENAME_QSget_rows ( EGLPNUM_TYPENAME_QSprob p, int **rowcnt, int **rowbeg, int **rowind,
			EGLPNUM_TYPE ** rowval, EGLPNUM_TYPE ** rhs, char **sense, char ***names),
	  EGLPNUM_TYPENAME_QSget_columns_list ( EGLPNUM_TYPENAME_QSprob p, int num, int *collist, int **colcnt,
			int **colbeg, int **colind, EGLPNUM_TYPE ** colval, EGLPNUM_TYPE ** obj,
			EGLPNUM_TYPE ** lower, EGLPNUM_TYPE ** upper, char ***names),
		EGLPNUM_TYPENAME_QSget_columns ( EGLPNUM_TYPENAME_QSprob p, int **colcnt, int **colbeg, int **colind,
			EGLPNUM_TYPE ** colval, EGLPNUM_TYPE ** obj, EGLPNUM_TYPE ** lower,
			EGLPNUM_TYPE ** upper, char ***names),
		EGLPNUM_TYPENAME_QSget_rownames ( EGLPNUM_TYPENAME_QSprob p, char **rownames),
		EGLPNUM_TYPENAME_QSget_colnames ( EGLPNUM_TYPENAME_QSprob p, char **colnames),
		EGLPNUM_TYPENAME_QSget_bound ( EGLPNUM_TYPENAME_QSprob p, int colindex, int lu, EGLPNUM_TYPE * bound),
		EGLPNUM_TYPENAME_QSget_bounds ( EGLPNUM_TYPENAME_QSprob p, EGLPNUM_TYPE * lower, EGLPNUM_TYPE * upper),
		EGLPNUM_TYPENAME_QSget_bounds_list(EGLPNUM_TYPENAME_QSprob p, int num, int*collist, EGLPNUM_TYPE*lb,
			EGLPNUM_TYPE*ub),
		EGLPNUM_TYPENAME_QSget_intflags ( EGLPNUM_TYPENAME_QSprob p, int *intflags),
		EGLPNUM_TYPENAME_QSget_intcount ( EGLPNUM_TYPENAME_QSprob p, int *count),
		EGLPNUM_TYPENAME_QSget_column_index ( EGLPNUM_TYPENAME_QSprob p, const char *name, int *colindex),
		EGLPNUM_TYPENAME_QSget_row_index ( EGLPNUM_TYPENAME_QSprob p, const char *name, int *rowindex),
		EGLPNUM_TYPENAME_QSget_named_x ( EGLPNUM_TYPENAME_QSprob p, const char *colname, EGLPNUM_TYPE * val),
		EGLPNUM_TYPENAME_QSget_named_rc ( EGLPNUM_TYPENAME_QSprob p, const char *colname, EGLPNUM_TYPE * val),
		EGLPNUM_TYPENAME_QSget_named_pi ( EGLPNUM_TYPENAME_QSprob p, const char *rowname, EGLPNUM_TYPE * val),
		EGLPNUM_TYPENAME_QSget_named_slack ( EGLPNUM_TYPENAME_QSprob p, const char *rowname, EGLPNUM_TYPE * val),
		EGLPNUM_TYPENAME_QScompute_row_norms ( EGLPNUM_TYPENAME_QSprob p),
		EGLPNUM_TYPENAME_QSwrite_prob ( EGLPNUM_TYPENAME_QSprob p, const char *filename, const char *filetype),
		EGLPNUM_TYPENAME_QSwrite_prob_file ( EGLPNUM_TYPENAME_QSprob p, FILE * file, const char *filetype),
		EGLPNUM_TYPENAME_QSwrite_basis ( EGLPNUM_TYPENAME_QSprob p, EGLPNUM_TYPENAME_QSbas B, const char *filename),
		EGLPNUM_TYPENAME_QStest_row_norms ( EGLPNUM_TYPENAME_QSprob p),
		EGLPNUM_TYPENAME_QSget_itcnt(EGLPNUM_TYPENAME_QSprob p, int *pI_iter, int *pII_iter, int *dI_iter,
			int *dII_iter, int *tot_iter),
		EGLPNUM_TYPENAME_QSset_param ( EGLPNUM_TYPENAME_QSprob p, int whichparam, int newvalue),
		EGLPNUM_TYPENAME_QSset_param_EGlpNum ( EGLPNUM_TYPENAME_QSprob p, int whichparam, EGLPNUM_TYPE newvalue),
		EGLPNUM_TYPENAME_QSget_param ( EGLPNUM_TYPENAME_QSprob p, int whichparam, int *value),
		EGLPNUM_TYPENAME_QSget_param_EGlpNum ( EGLPNUM_TYPENAME_QSprob p, int whichparam, EGLPNUM_TYPE * value); 

EGLPNUM_TYPENAME_QSLIB_INTERFACE char *EGLPNUM_TYPENAME_QSget_probname ( EGLPNUM_TYPENAME_QSprob p);
EGLPNUM_TYPENAME_QSLIB_INTERFACE char *EGLPNUM_TYPENAME_QSget_objname ( EGLPNUM_TYPENAME_QSprob p);
EGLPNUM_TYPENAME_QSLIB_INTERFACE char *EGLPNUM_TYPENAME_QSversion ( void);

EGLPNUM_TYPENAME_QSLIB_INTERFACE EGLPNUM_TYPENAME_QSprob EGLPNUM_TYPENAME_QScreate_prob ( const char *name, int objsense),
		EGLPNUM_TYPENAME_QSread_prob ( const char *filename, const char *filetype),
		EGLPNUM_TYPENAME_QSload_prob ( const char *probname, int ncols, int nrows, int *cmatcnt,
			int *cmatbeg, int *cmatind, EGLPNUM_TYPE * cmatval, int objsense,
			EGLPNUM_TYPE * obj, EGLPNUM_TYPE * rhs, char *sense, EGLPNUM_TYPE * lower,
			EGLPNUM_TYPE * upper, const char **colnames, const char **rownames),
		EGLPNUM_TYPENAME_QScopy_prob ( EGLPNUM_TYPENAME_QSprob p, const char *newname);

EGLPNUM_TYPENAME_QSLIB_INTERFACE EGLPNUM_TYPENAME_QSbas EGLPNUM_TYPENAME_QSget_basis ( EGLPNUM_TYPENAME_QSprob p),
		EGLPNUM_TYPENAME_QSread_basis ( EGLPNUM_TYPENAME_QSprob p, const char *filename);

#ifdef  __cplusplus
}
#endif

/****************************************************************************
 *
 * This is the undocumented part of the QSlib interface 
 *
 ****************************************************************************/
/* 
 * functions to facilitate line by line reading from other sources than 
 * files from within MPS/LP parsers  
 * 
 * functions to facilitate the collection of error information instead of 
 * having the parsers print messages to stderr
 *                              by mps/lp format writers
 * 
 * a problem's reporter is used by the solver code to provide important 
 * feedback/progress information
 */

#ifdef WIN32
typedef struct EGLPNUM_TYPENAME_QSLIB_INTERFACE EGLPNUM_TYPENAME_qsline_reader *EGLPNUM_TYPENAME_QSline_reader;
typedef struct EGLPNUM_TYPENAME_QSLIB_INTERFACE EGLPNUM_TYPENAME_qsformat_error *EGLPNUM_TYPENAME_QSformat_error;
typedef struct EGLPNUM_TYPENAME_QSLIB_INTERFACE EGLPNUM_TYPENAME_qserror_collector *EGLPNUM_TYPENAME_QSerror_collector;
typedef struct EGLPNUM_TYPENAME_QSLIB_INTERFACE EGLPNUM_TYPENAME_qserror_memory *EGLPNUM_TYPENAME_QSerror_memory;
#else
typedef struct EGLPNUM_TYPENAME_qsline_reader *EGLPNUM_TYPENAME_QSline_reader;
typedef struct EGLPNUM_TYPENAME_qsformat_error *EGLPNUM_TYPENAME_QSformat_error;
typedef struct EGLPNUM_TYPENAME_qserror_collector *EGLPNUM_TYPENAME_QSerror_collector;
typedef struct EGLPNUM_TYPENAME_qserror_memory *EGLPNUM_TYPENAME_QSerror_memory;
#endif

#ifdef  __cplusplus
extern "C"
{
#endif

EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSgrab_cache (EGLPNUM_TYPENAME_QSprob p, int status);

	EGLPNUM_TYPENAME_QSLIB_INTERFACE const char *EGLPNUM_TYPENAME_QSformat_error_type_string (
	int tp);

	EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSerror_get_type (
	EGLPNUM_TYPENAME_QSformat_error error);
	EGLPNUM_TYPENAME_QSLIB_INTERFACE const char *EGLPNUM_TYPENAME_QSerror_get_desc (
	EGLPNUM_TYPENAME_QSformat_error error);
	EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSerror_get_line_number (
	EGLPNUM_TYPENAME_QSformat_error error);
	EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSerror_get_pos (
	EGLPNUM_TYPENAME_QSformat_error error);
	EGLPNUM_TYPENAME_QSLIB_INTERFACE const char *EGLPNUM_TYPENAME_QSerror_get_line (
	EGLPNUM_TYPENAME_QSformat_error error);
	EGLPNUM_TYPENAME_QSLIB_INTERFACE void EGLPNUM_TYPENAME_QSerror_print (
	FILE * f,
	EGLPNUM_TYPENAME_QSformat_error error);

	EGLPNUM_TYPENAME_QSLIB_INTERFACE EGLPNUM_TYPENAME_QSerror_collector EGLPNUM_TYPENAME_QSerror_collector_new (
	void *fct,
	void *dest);
	EGLPNUM_TYPENAME_QSLIB_INTERFACE EGLPNUM_TYPENAME_QSerror_collector EGLPNUM_TYPENAME_QSerror_memory_collector_new (
	EGLPNUM_TYPENAME_QSerror_memory mem);
	EGLPNUM_TYPENAME_QSLIB_INTERFACE void EGLPNUM_TYPENAME_QSerror_collector_free (
	EGLPNUM_TYPENAME_QSerror_collector c);

/****************************************************************************
 * line reader 
 */
	EGLPNUM_TYPENAME_QSLIB_INTERFACE EGLPNUM_TYPENAME_QSline_reader EGLPNUM_TYPENAME_QSline_reader_new (
	void *fct,
	void *data_src);
	/* reader->read_line_fct defaults to fgets */

	EGLPNUM_TYPENAME_QSLIB_INTERFACE void EGLPNUM_TYPENAME_QSline_reader_free (
	EGLPNUM_TYPENAME_QSline_reader reader);

	EGLPNUM_TYPENAME_QSLIB_INTERFACE void EGLPNUM_TYPENAME_QSline_reader_set_error_collector (
	EGLPNUM_TYPENAME_QSline_reader reader,
	EGLPNUM_TYPENAME_QSerror_collector collector);

	EGLPNUM_TYPENAME_QSLIB_INTERFACE char *EGLPNUM_TYPENAME_QSline_reader_get (
	EGLPNUM_TYPENAME_QSline_reader reader,
	char *s,
	int size);

	EGLPNUM_TYPENAME_QSLIB_INTERFACE EGLPNUM_TYPENAME_QSprob EGLPNUM_TYPENAME_QSget_prob (
	EGLPNUM_TYPENAME_QSline_reader reader,
	const char *probname,
	const char *filetype);
	/* the MPS and LP parsers uses the fct from reader 
	 * to get to next input line */


/****************************************************************************
 * error memory 
 */
	EGLPNUM_TYPENAME_QSLIB_INTERFACE EGLPNUM_TYPENAME_QSerror_memory EGLPNUM_TYPENAME_QSerror_memory_create (
	int takeErrorLines);
	EGLPNUM_TYPENAME_QSLIB_INTERFACE void EGLPNUM_TYPENAME_QSerror_memory_free (
	EGLPNUM_TYPENAME_QSerror_memory mem);

	EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSerror_memory_get_nof (
	EGLPNUM_TYPENAME_QSerror_memory mem,
	int error_type);
	EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSerror_memory_get_nerrors (
	EGLPNUM_TYPENAME_QSerror_memory mem);

	EGLPNUM_TYPENAME_QSLIB_INTERFACE EGLPNUM_TYPENAME_QSformat_error EGLPNUM_TYPENAME_QSerror_memory_get_last_error (
	EGLPNUM_TYPENAME_QSerror_memory mem);
	EGLPNUM_TYPENAME_QSLIB_INTERFACE EGLPNUM_TYPENAME_QSformat_error EGLPNUM_TYPENAME_QSerror_memory_get_prev_error (
	EGLPNUM_TYPENAME_QSformat_error e);

/**************************************************************************** 
 * reporter for solver feedback 
 */
	EGLPNUM_TYPENAME_QSLIB_INTERFACE void EGLPNUM_TYPENAME_QSset_reporter (
	EGLPNUM_TYPENAME_QSprob prob,
	int iterskip,
	void *fct,
	void *dest);

	EGLPNUM_TYPENAME_QSLIB_INTERFACE int EGLPNUM_TYPENAME_QSreport_prob (
	EGLPNUM_TYPENAME_QSprob p,
	const char *filetype,
	EGLPNUM_TYPENAME_QSerror_collector c);

#ifdef  __cplusplus
}
#endif
#endif													/* EGLPNUM_TYPENAME___QS_QSOPT_H */
