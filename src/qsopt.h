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

/*  $RCSfile: qsopt.h,v $ $Revision: 1.3 $ $Date: 2003/11/05 16:57:39 $" */
#ifndef __QS_QSOPT_H
#define __QS_QSOPT_H

#include <stdio.h>
#include "qs_config.h"

#ifdef WIN32

#ifdef QSLIB_EXPORTS
#define QSLIB_INTERFACE __declspec(dllexport)
#else
#define QSLIB_INTERFACE __declspec(dllimport)
#endif

#else
#define QSLIB_INTERFACE extern
#endif

#ifdef WIN32
typedef struct QSLIB_INTERFACE qsdata *QSprob;
typedef struct QSLIB_INTERFACE qsbasis *QSbas;
#else
typedef struct qsdata *QSprob;
typedef struct qsbasis *QSbas;
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
 *     solver_main/reader_main part of DLL
 */
QSLIB_INTERFACE int solver_main ( int argc, char **argv);
QSLIB_INTERFACE int reader_main ( int argc, char **argv);
#endif

QSLIB_INTERFACE void QSfree ( void *ptr),
		QSfree_prob ( QSprob p),
		QSfree_basis ( QSbas B),
	  QSset_precision ( const unsigned prec),/**< set the precision for floating 
																								 point numbers to the given 
																								 number of bits */
		QSstart ( void),/**< whe we use non native numbers, we need to make 
												 some initializations before operating with the
												 library */
	  QSend ( void);	/**< just to free any internal static data needed by
												 the variable precision numbers */

QSLIB_INTERFACE int QSopt_primal ( QSprob p, int *status),
		QSopt_dual ( QSprob p, int *status),
		QSopt_pivotin_col ( QSprob p, int ccnt, int *clist),
		QSopt_pivotin_row ( QSprob p, int rcnt, int *rlist),
		QSopt_strongbranch ( QSprob p, int ncand, int *candidatelist,
			EGlpNum_t * xlist, EGlpNum_t * down_vals, EGlpNum_t * up_vals,
			int iterations, EGlpNum_t objbound),
		QSchange_objsense ( QSprob p, int newsense),
		QSget_objsense ( QSprob p, int *newsense),
		QSnew_col ( QSprob p,const EGlpNum_t obj,const EGlpNum_t lower,const EGlpNum_t upper,
			const char *name),
		QSadd_cols ( QSprob p, int num, int *cmatcnt, int *cmatbeg, int *cmatind,
			EGlpNum_t * cmatval, EGlpNum_t * obj, EGlpNum_t * lower,
			EGlpNum_t * upper, const char **names),
		QSadd_col ( QSprob p, int cnt, int *cmatind, EGlpNum_t * cmatval,
			EGlpNum_t obj, EGlpNum_t lower, EGlpNum_t upper, const char *name),
		QSnew_row ( QSprob p,const EGlpNum_t rhs, int sense, const char *name),
		QSadd_ranged_rows ( QSprob p, int num, int *rmatcnt, int *rmatbeg, 
			int *rmatind,const EGlpNum_t * rmatval,const EGlpNum_t * rhs, char *sense,
			const EGlpNum_t* range, const char **names),
		QSadd_ranged_row ( QSprob p, int cnt, int *rmatind,const EGlpNum_t * rmatval,
			const EGlpNum_t * rhs, int sense,const EGlpNum_t * range, const char *name),
		QSadd_rows ( QSprob p, int num, int *rmatcnt, int *rmatbeg, int *rmatind,
			const EGlpNum_t * rmatval,const EGlpNum_t * rhs, char *sense, const char **names),
		QSadd_row ( QSprob p, int cnt, int *rmatind,const EGlpNum_t * rmatval,
			const EGlpNum_t * rhs, int sense, const char *name),
		QSdelete_rows ( QSprob p, int num, int *dellist),
		QSdelete_row ( QSprob p, int rowindex),
		QSdelete_setrows ( QSprob p, int *flags),
		QSdelete_named_row ( QSprob p, const char *rowname),
		QSdelete_named_rows_list ( QSprob p, int num, const char **rownames),
		QSdelete_cols ( QSprob p, int num, int *dellist),
		QSdelete_col ( QSprob p, int colindex),
		QSdelete_setcols ( QSprob p, int *flags),
		QSdelete_named_column ( QSprob p, const char *colname),
		QSdelete_named_columns_list ( QSprob p, int num, const char **colnames),
		QSchange_senses ( QSprob p, int num, int *rowlist, char *sense),
		QSchange_sense ( QSprob p, int rowindex, int sense),
		QSchange_coef ( QSprob p, int rowindex, int colindex, EGlpNum_t coef),
		QSchange_objcoef ( QSprob p, int indx, EGlpNum_t coef),
		QSchange_rhscoef ( QSprob p, int indx, EGlpNum_t coef),
		QSchange_range(QSprob p, int rowindex, EGlpNum_t range),
		QSchange_bounds ( QSprob p, int num, int *collist, char *lu, 
			const EGlpNum_t * bounds),
		QSchange_bound ( QSprob p, int indx, int lu,const EGlpNum_t bound),
		QSload_basis ( QSprob p, QSbas B),
		QSread_and_load_basis ( QSprob p, const char *filename),
		QSload_basis_array ( QSprob p, char *cstat, char *rstat),
		QSload_basis_and_row_norms_array ( QSprob p, char *cstat, char *rstat, 
			EGlpNum_t * rownorms),
		QSget_basis_array ( QSprob p, char *cstat, char *rstat),
		QSget_basis_and_row_norms_array ( QSprob p, char *cstat, char *rstat,
			EGlpNum_t * rownorms),
		QSget_binv_row ( QSprob p, int indx, EGlpNum_t * binvrow),
		QSget_tableau_row ( QSprob p, int indx, EGlpNum_t * tableaurow),
		QSget_basis_order ( QSprob p, int *basorder), 
    QSget_coef (QSprob p, int rowindex, int colindex, EGlpNum_t*coef),
		QSget_status ( QSprob p, int *status),
		QSget_solution ( QSprob p, EGlpNum_t * value, EGlpNum_t * x,
			EGlpNum_t * pi, EGlpNum_t * slack, EGlpNum_t * rc),
		QSget_objval ( QSprob p, EGlpNum_t * value),
		QSget_pi_array ( QSprob p, EGlpNum_t * pi),
		QSget_rc_array ( QSprob p, EGlpNum_t * rc),
		QSget_x_array ( QSprob p, EGlpNum_t * x),
		QSget_slack_array ( QSprob p, EGlpNum_t * slack),
		QSget_infeas_array ( QSprob p, EGlpNum_t * pi),
		QSget_colcount ( QSprob p),
		QSget_rowcount ( QSprob p),
		QSget_nzcount ( QSprob p),
		QSget_obj_list(QSprob p, int num, int*collist, EGlpNum_t*obj),
		QSget_obj ( QSprob p, EGlpNum_t * obj),
		QSget_rhs ( QSprob p, EGlpNum_t * rhs),
		QSget_ranged_rows_list ( QSprob p, int num, int *rowlist, int **rowcnt,
			int **rowbeg, int **rowind, EGlpNum_t ** rowval, EGlpNum_t ** rhs,
			char **sense, EGlpNum_t **range, char ***names),
		QSget_ranged_rows ( QSprob p, int **rowcnt, int **rowbeg, int **rowind,
			EGlpNum_t ** rowval, EGlpNum_t ** rhs, char **sense, 
			EGlpNum_t ** range, char ***names),
		QSget_senses ( QSprob p, char*senses),
		QSget_rows_list ( QSprob p, int num, int *rowlist, int **rowcnt,
			int **rowbeg, int **rowind, EGlpNum_t ** rowval, EGlpNum_t ** rhs,
			char **sense, char ***names),
		QSget_rows ( QSprob p, int **rowcnt, int **rowbeg, int **rowind,
			EGlpNum_t ** rowval, EGlpNum_t ** rhs, char **sense, char ***names),
	  QSget_columns_list ( QSprob p, int num, int *collist, int **colcnt,
			int **colbeg, int **colind, EGlpNum_t ** colval, EGlpNum_t ** obj,
			EGlpNum_t ** lower, EGlpNum_t ** upper, char ***names),
		QSget_columns ( QSprob p, int **colcnt, int **colbeg, int **colind,
			EGlpNum_t ** colval, EGlpNum_t ** obj, EGlpNum_t ** lower,
			EGlpNum_t ** upper, char ***names),
		QSget_rownames ( QSprob p, char **rownames),
		QSget_colnames ( QSprob p, char **colnames),
		QSget_bound ( QSprob p, int colindex, int lu, EGlpNum_t * bound),
		QSget_bounds ( QSprob p, EGlpNum_t * lower, EGlpNum_t * upper),
		QSget_bounds_list(QSprob p, int num, int*collist, EGlpNum_t*lb,
			EGlpNum_t*ub),
		QSget_intflags ( QSprob p, int *intflags),
		QSget_intcount ( QSprob p, int *count),
		QSget_column_index ( QSprob p, const char *name, int *colindex),
		QSget_row_index ( QSprob p, const char *name, int *rowindex),
		QSget_named_x ( QSprob p, const char *colname, EGlpNum_t * val),
		QSget_named_rc ( QSprob p, const char *colname, EGlpNum_t * val),
		QSget_named_pi ( QSprob p, const char *rowname, EGlpNum_t * val),
		QSget_named_slack ( QSprob p, const char *rowname, EGlpNum_t * val),
		QScompute_row_norms ( QSprob p),
		QSwrite_prob ( QSprob p, const char *filename, const char *filetype),
		QSwrite_prob_file ( QSprob p, FILE * file, const char *filetype),
		QSwrite_basis ( QSprob p, QSbas B, const char *filename),
		QStest_row_norms ( QSprob p),
		QSget_itcnt(QSprob p, int *pI_iter, int *pII_iter, int *dI_iter,
			int *dII_iter, int *tot_iter),
		QSset_param ( QSprob p, int whichparam, int newvalue),
		QSset_param_EGlpNum ( QSprob p, int whichparam, EGlpNum_t newvalue),
		QSget_param ( QSprob p, int whichparam, int *value),
		QSget_param_EGlpNum ( QSprob p, int whichparam, EGlpNum_t * value); 

QSLIB_INTERFACE char *QSget_probname ( QSprob p);
QSLIB_INTERFACE char *QSget_objname ( QSprob p);
QSLIB_INTERFACE char *QSversion ( void);

QSLIB_INTERFACE QSprob QScreate_prob ( const char *name, int objsense),
		QSread_prob ( const char *filename, const char *filetype),
		QSload_prob ( const char *probname, int ncols, int nrows, int *cmatcnt,
			int *cmatbeg, int *cmatind, EGlpNum_t * cmatval, int objsense,
			EGlpNum_t * obj, EGlpNum_t * rhs, char *sense, EGlpNum_t * lower,
			EGlpNum_t * upper, const char **colnames, const char **rownames),
		QScopy_prob ( QSprob p, const char *newname);

QSLIB_INTERFACE QSbas QSget_basis ( QSprob p),
		QSread_basis ( QSprob p, const char *filename);

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
typedef struct QSLIB_INTERFACE qsline_reader *QSline_reader;
typedef struct QSLIB_INTERFACE qsformat_error *QSformat_error;
typedef struct QSLIB_INTERFACE qserror_collector *QSerror_collector;
typedef struct QSLIB_INTERFACE qserror_memory *QSerror_memory;
#else
typedef struct qsline_reader *QSline_reader;
typedef struct qsformat_error *QSformat_error;
typedef struct qserror_collector *QSerror_collector;
typedef struct qserror_memory *QSerror_memory;
#endif

#ifdef  __cplusplus
extern "C"
{
#endif
	QSLIB_INTERFACE const char *QSformat_error_type_string (
	int tp);

	QSLIB_INTERFACE int QSerror_get_type (
	QSformat_error error);
	QSLIB_INTERFACE const char *QSerror_get_desc (
	QSformat_error error);
	QSLIB_INTERFACE int QSerror_get_line_number (
	QSformat_error error);
	QSLIB_INTERFACE int QSerror_get_pos (
	QSformat_error error);
	QSLIB_INTERFACE const char *QSerror_get_line (
	QSformat_error error);
	QSLIB_INTERFACE void QSerror_print (
	FILE * f,
	QSformat_error error);

	QSLIB_INTERFACE QSerror_collector QSerror_collector_new (
	void *fct,
	void *dest);
	QSLIB_INTERFACE QSerror_collector QSerror_memory_collector_new (
	QSerror_memory mem);
	QSLIB_INTERFACE void QSerror_collector_free (
	QSerror_collector c);

/****************************************************************************
 * line reader 
 */
	QSLIB_INTERFACE QSline_reader QSline_reader_new (
	void *fct,
	void *data_src);
	/* reader->read_line_fct defaults to fgets */

	QSLIB_INTERFACE void QSline_reader_free (
	QSline_reader reader);

	QSLIB_INTERFACE void QSline_reader_set_error_collector (
	QSline_reader reader,
	QSerror_collector collector);

	QSLIB_INTERFACE char *QSline_reader_get (
	QSline_reader reader,
	char *s,
	int size);

	QSLIB_INTERFACE QSprob QSget_prob (
	QSline_reader reader,
	const char *probname,
	const char *filetype);
	/* the MPS and LP parsers uses the fct from reader 
	 * to get to next input line */


/****************************************************************************
 * error memory 
 */
	QSLIB_INTERFACE QSerror_memory QSerror_memory_create (
	int takeErrorLines);
	QSLIB_INTERFACE void QSerror_memory_free (
	QSerror_memory mem);

	QSLIB_INTERFACE int QSerror_memory_get_nof (
	QSerror_memory mem,
	int error_type);
	QSLIB_INTERFACE int QSerror_memory_get_nerrors (
	QSerror_memory mem);

	QSLIB_INTERFACE QSformat_error QSerror_memory_get_last_error (
	QSerror_memory mem);
	QSLIB_INTERFACE QSformat_error QSerror_memory_get_prev_error (
	QSformat_error e);

/**************************************************************************** 
 * reporter for solver feedback 
 */
	QSLIB_INTERFACE void QSset_reporter (
	QSprob prob,
	int iterskip,
	void *fct,
	void *dest);

	QSLIB_INTERFACE int QSreport_prob (
	QSprob p,
	const char *filetype,
	QSerror_collector c);

#ifdef  __cplusplus
}
#endif
#endif													/* __QS_QSOPT_H */
