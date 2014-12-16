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

/* RCSINFO $Id: lib_EGLPNUM_TYPENAME.h,v 1.4 2003/11/05 17:00:26 meven Exp $ */
#ifndef EGLPNUM_TYPENAME_ILL_LIB_H
#define EGLPNUM_TYPENAME_ILL_LIB_H

#include "lpdefs_EGLPNUM_TYPENAME.h"
#include "lpdata_EGLPNUM_TYPENAME.h"
#include "price_EGLPNUM_TYPENAME.h"
#include "basicdefs.h"

/****************************************************************************/
/*                                                                          */
/*                   Return Status for EGLPNUM_TYPENAME_ILLlib_optimize                      */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/*                               lib.c                                      */
/*                                                                          */
/****************************************************************************/

struct itcnt_t;

int EGLPNUM_TYPENAME_ILLlib_optimize ( EGLPNUM_TYPENAME_lpinfo * lp, EGLPNUM_TYPENAME_ILLlp_basis * B, EGLPNUM_TYPENAME_price_info * pinf, int algo,
			int *status, int simplex_display, struct itcnt_t*itcnt),
		EGLPNUM_TYPENAME_ILLlib_cache_solution ( EGLPNUM_TYPENAME_lpinfo * lp, EGLPNUM_TYPENAME_ILLlp_cache * C), 
		EGLPNUM_TYPENAME_ILLlib_solution ( EGLPNUM_TYPENAME_lpinfo * lp, EGLPNUM_TYPENAME_ILLlp_cache * C, EGLPNUM_TYPE * val, 
			EGLPNUM_TYPE * x, EGLPNUM_TYPE * pi, EGLPNUM_TYPE * slack, EGLPNUM_TYPE * rc),
		EGLPNUM_TYPENAME_ILLlib_get_x ( EGLPNUM_TYPENAME_lpinfo * lp, EGLPNUM_TYPENAME_ILLlp_cache * C, EGLPNUM_TYPE * x),
		EGLPNUM_TYPENAME_ILLlib_get_slack ( EGLPNUM_TYPENAME_lpinfo * lp, EGLPNUM_TYPENAME_ILLlp_cache * C, EGLPNUM_TYPE * slack),
		EGLPNUM_TYPENAME_ILLlib_objval ( EGLPNUM_TYPENAME_lpinfo * lp, EGLPNUM_TYPENAME_ILLlp_cache * C, EGLPNUM_TYPE * val),
		EGLPNUM_TYPENAME_ILLlib_tableau ( EGLPNUM_TYPENAME_lpinfo * lp, int row, EGLPNUM_TYPE * binv, EGLPNUM_TYPE * tabrow),
		EGLPNUM_TYPENAME_ILLlib_basis_order ( EGLPNUM_TYPENAME_lpinfo * lp, int *header),
		EGLPNUM_TYPENAME_ILLlib_newrow ( EGLPNUM_TYPENAME_lpinfo * lp, EGLPNUM_TYPENAME_ILLlp_basis * B,const EGLPNUM_TYPE rhs, int sense,
			const EGLPNUM_TYPE range, const char *name),
		EGLPNUM_TYPENAME_ILLlib_newrows ( EGLPNUM_TYPENAME_lpinfo * lp, EGLPNUM_TYPENAME_ILLlp_basis * B, int num,const EGLPNUM_TYPE * rhs,
			char *sense, const EGLPNUM_TYPE * range, const char **names),
		EGLPNUM_TYPENAME_ILLlib_addrow ( EGLPNUM_TYPENAME_lpinfo * lp, EGLPNUM_TYPENAME_ILLlp_basis * B, int cnt, int *ind,
			const EGLPNUM_TYPE * val, const EGLPNUM_TYPE rhs, int sense,const EGLPNUM_TYPE range,
			const char *rowname),
		EGLPNUM_TYPENAME_ILLlib_addrows ( EGLPNUM_TYPENAME_lpinfo * lp, EGLPNUM_TYPENAME_ILLlp_basis * B, int num, int *rmatcnt, 
			int *rmatbeg, int *rmatind,const EGLPNUM_TYPE * rmatval,const EGLPNUM_TYPE * rhs, 
			char *sense, const EGLPNUM_TYPE * range, const char **names, int *nofactor),
		EGLPNUM_TYPENAME_ILLlib_delrows ( EGLPNUM_TYPENAME_lpinfo * lp, EGLPNUM_TYPENAME_ILLlp_basis * B, EGLPNUM_TYPENAME_ILLlp_cache * C, int num, 
			int *dellist, int *basis_ok, int *cache_ok), 
		EGLPNUM_TYPENAME_ILLlib_newcol ( EGLPNUM_TYPENAME_lpinfo * lp, EGLPNUM_TYPENAME_ILLlp_basis * B,const EGLPNUM_TYPE obj, 
			const EGLPNUM_TYPE lower,const EGLPNUM_TYPE upper, const char *name, int factorok), 
		EGLPNUM_TYPENAME_ILLlib_newcols ( EGLPNUM_TYPENAME_lpinfo * lp, EGLPNUM_TYPENAME_ILLlp_basis * B, int num, EGLPNUM_TYPE * obj,
			EGLPNUM_TYPE * lower, EGLPNUM_TYPE * upper, const char **names, int factorok),
		EGLPNUM_TYPENAME_ILLlib_addcol ( EGLPNUM_TYPENAME_lpinfo * lp, EGLPNUM_TYPENAME_ILLlp_basis * B, int cnt, int *ind,
			EGLPNUM_TYPE * val,const EGLPNUM_TYPE obj,const EGLPNUM_TYPE lower,const EGLPNUM_TYPE upper, 
			const char *name, int factorok),
		EGLPNUM_TYPENAME_ILLlib_addcols ( EGLPNUM_TYPENAME_lpinfo * lp, EGLPNUM_TYPENAME_ILLlp_basis * B, int num, int *cmatcnt,
			int *cmatbeg, int *cmatind, EGLPNUM_TYPE * cmatval, EGLPNUM_TYPE * obj,
			EGLPNUM_TYPE * lower, EGLPNUM_TYPE * upper, const char **names, int factorok),
		EGLPNUM_TYPENAME_ILLlib_delcols ( EGLPNUM_TYPENAME_lpinfo * lp, EGLPNUM_TYPENAME_ILLlp_basis * B, int num, int *dellist,
			int *basis_ok),
		EGLPNUM_TYPENAME_ILLlib_chgcoef ( EGLPNUM_TYPENAME_lpinfo * lp, int rowindex, int colindex, EGLPNUM_TYPE coef), 
		EGLPNUM_TYPENAME_ILLlib_getcoef (EGLPNUM_TYPENAME_lpinfo *lp, int rowindex, int colindex, EGLPNUM_TYPE* coef),
		EGLPNUM_TYPENAME_ILLlib_chgrange (EGLPNUM_TYPENAME_lpinfo *lp, int indx, EGLPNUM_TYPE coef),
		EGLPNUM_TYPENAME_ILLlib_chgsense ( EGLPNUM_TYPENAME_lpinfo * lp, int num, int *rowlist, char *sense),
    EGLPNUM_TYPENAME_ILLlib_getsenses (EGLPNUM_TYPENAME_lpinfo *lp, char *senses),
		EGLPNUM_TYPENAME_ILLlib_getrows ( EGLPNUM_TYPENAME_lpinfo * lp, int num, int *rowlist, int **rowcnt,
			int **rowbeg, int **rowind, EGLPNUM_TYPE ** rowval, EGLPNUM_TYPE ** rhs,
			char **sense, EGLPNUM_TYPE ** range, char ***names),
		EGLPNUM_TYPENAME_ILLlib_getcols ( EGLPNUM_TYPENAME_lpinfo * lp, int num, int *collist, int **colcnt,
			int **colbeg, int **colind, EGLPNUM_TYPE ** colval, EGLPNUM_TYPE ** obj,
			EGLPNUM_TYPE ** lower, EGLPNUM_TYPE ** upper, char ***names),
		EGLPNUM_TYPENAME_ILLlib_getobj ( EGLPNUM_TYPENAME_lpinfo * lp, EGLPNUM_TYPE * obj),
		EGLPNUM_TYPENAME_ILLlib_getobj_list (EGLPNUM_TYPENAME_lpinfo *lp, int num, int* collist, EGLPNUM_TYPE* obj),
		EGLPNUM_TYPENAME_ILLlib_chgobj ( EGLPNUM_TYPENAME_lpinfo * lp, int indx, EGLPNUM_TYPE coef),
		EGLPNUM_TYPENAME_ILLlib_getrhs ( EGLPNUM_TYPENAME_lpinfo * lp, EGLPNUM_TYPE * rhs),
		EGLPNUM_TYPENAME_ILLlib_chgrhs ( EGLPNUM_TYPENAME_lpinfo * lp, int indx, EGLPNUM_TYPE coef),
		EGLPNUM_TYPENAME_ILLlib_getintflags ( EGLPNUM_TYPENAME_lpinfo * lp, int *intflags),
		EGLPNUM_TYPENAME_ILLlib_rownames ( EGLPNUM_TYPENAME_lpinfo * lp, char **rownames),
		EGLPNUM_TYPENAME_ILLlib_colnames ( EGLPNUM_TYPENAME_lpinfo * lp, char **colnames),
		EGLPNUM_TYPENAME_ILLlib_colindex ( EGLPNUM_TYPENAME_lpinfo * lp, const char *name, int *colindex),
		EGLPNUM_TYPENAME_ILLlib_rowindex ( EGLPNUM_TYPENAME_lpinfo * lp, const char *name, int *rowindex),
		EGLPNUM_TYPENAME_ILLlib_chgbnd ( EGLPNUM_TYPENAME_lpinfo * lp, int indx, int lu,const EGLPNUM_TYPE bnd),
		EGLPNUM_TYPENAME_ILLlib_chgbnds ( EGLPNUM_TYPENAME_lpinfo * lp, int cnt, int *indx, char *lu, const EGLPNUM_TYPE * bnd),
		EGLPNUM_TYPENAME_ILLlib_getbnd ( EGLPNUM_TYPENAME_lpinfo * lp, int indx, int lu, EGLPNUM_TYPE * bnd),
		EGLPNUM_TYPENAME_ILLlib_getbnds ( EGLPNUM_TYPENAME_lpinfo * lp, EGLPNUM_TYPE * lower, EGLPNUM_TYPE * upper),
		EGLPNUM_TYPENAME_ILLlib_getbnds_list ( EGLPNUM_TYPENAME_lpinfo *lp, int num, int*collist, EGLPNUM_TYPE *lower,
			EGLPNUM_TYPE *upper),
		EGLPNUM_TYPENAME_ILLlib_strongbranch ( EGLPNUM_TYPENAME_lpinfo * lp, EGLPNUM_TYPENAME_price_info * pinf, int *candidatelist,
			int ncand, EGLPNUM_TYPE * xlist, EGLPNUM_TYPE * downpen, EGLPNUM_TYPE * uppen,
			int iterations, EGLPNUM_TYPE objbound, struct itcnt_t*itcnt),
		EGLPNUM_TYPENAME_ILLlib_getbasis ( EGLPNUM_TYPENAME_lpinfo * lp, char *cstat, char *rstat),
		EGLPNUM_TYPENAME_ILLlib_loadbasis ( EGLPNUM_TYPENAME_ILLlp_basis * B, int nstruct, int nrows, char *cstat,
			char *rstat),
		EGLPNUM_TYPENAME_ILLlib_readbasis ( EGLPNUM_TYPENAME_lpinfo * lp, EGLPNUM_TYPENAME_ILLlp_basis * B, const char *fname),
		EGLPNUM_TYPENAME_ILLlib_writebasis ( EGLPNUM_TYPENAME_lpinfo * lp, EGLPNUM_TYPENAME_ILLlp_basis * B, const char *fname),
		EGLPNUM_TYPENAME_ILLlib_getrownorms ( EGLPNUM_TYPENAME_lpinfo * lp, EGLPNUM_TYPENAME_price_info * pinf, EGLPNUM_TYPE * rownorms),
		EGLPNUM_TYPENAME_ILLlib_loadrownorms ( EGLPNUM_TYPENAME_lpinfo * lp, EGLPNUM_TYPENAME_price_info * pinf, EGLPNUM_TYPE * rownorms),
		EGLPNUM_TYPENAME_ILLlib_recompute_rownorms ( EGLPNUM_TYPENAME_lpinfo * lp, EGLPNUM_TYPENAME_price_info * pinf),
		EGLPNUM_TYPENAME_ILLlib_iter ( EGLPNUM_TYPENAME_lpinfo * lp),
		EGLPNUM_TYPENAME_ILLlib_print_x ( EGioFile_t * fd, EGLPNUM_TYPENAME_lpinfo * lp, EGLPNUM_TYPENAME_ILLlp_cache * C, EGLPNUM_TYPE * x,
			int nonZerosOnly),
		EGLPNUM_TYPENAME_ILLwrite_lp_file ( EGLPNUM_TYPENAME_ILLlpdata * lp, EGioFile_t * eout, EGLPNUM_TYPENAME_qserror_collector * c);


extern int EGLPNUM_TYPENAME_ILLlib_findName (
	EGLPNUM_TYPENAME_ILLlpdata * qslp,
	int forRow,
	const char *name,
	int id,
	char buf[ILL_namebufsize]);

/****************************************************************************/
/*                                                                          */
/*                           presolve.c                                     */
/*                                                                          */
/****************************************************************************/

int EGLPNUM_TYPENAME_ILLpresolve_add_logicals (
	EGLPNUM_TYPENAME_ILLlpdata * lp);


/****************************************************************************/
/*                                                                          */
/*                            binary.c                                      */
/*                                                                          */
/****************************************************************************/

int EGLPNUM_TYPENAME_ILLmip_binary_dfs (
	EGLPNUM_TYPENAME_lpinfo * lp);

#endif /* EGLPNUM_TYPENAME_ILL_LIB_H */
