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

/* RCSINFO $Id: lib.h,v 1.4 2003/11/05 17:00:26 meven Exp $ */
#ifndef ILL_LIB_H
#define ILL_LIB_H

#include "lpdefs.h"
#include "lpdata.h"
#include "price.h"
#include "basicdefs.h"

/****************************************************************************/
/*                                                                          */
/*                   Return Status for ILLlib_optimize                      */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/*                               lib.c                                      */
/*                                                                          */
/****************************************************************************/

struct itcnt_t;

int ILLlib_optimize ( lpinfo * lp, ILLlp_basis * B, price_info * pinf, int algo,
			int *status, int simplex_display, struct itcnt_t*itcnt),
		ILLlib_cache_solution ( lpinfo * lp, ILLlp_cache * C), 
		ILLlib_solution ( lpinfo * lp, ILLlp_cache * C, EGlpNum_t * val, 
			EGlpNum_t * x, EGlpNum_t * pi, EGlpNum_t * slack, EGlpNum_t * rc),
		ILLlib_get_x ( lpinfo * lp, ILLlp_cache * C, EGlpNum_t * x),
		ILLlib_get_slack ( lpinfo * lp, ILLlp_cache * C, EGlpNum_t * slack),
		ILLlib_objval ( lpinfo * lp, ILLlp_cache * C, EGlpNum_t * val),
		ILLlib_tableau ( lpinfo * lp, int row, EGlpNum_t * binv, EGlpNum_t * tabrow),
		ILLlib_basis_order ( lpinfo * lp, int *header),
		ILLlib_newrow ( lpinfo * lp, ILLlp_basis * B,const EGlpNum_t rhs, int sense,
			const EGlpNum_t range, const char *name),
		ILLlib_newrows ( lpinfo * lp, ILLlp_basis * B, int num,const EGlpNum_t * rhs,
			char *sense, const EGlpNum_t * range, const char **names),
		ILLlib_addrow ( lpinfo * lp, ILLlp_basis * B, int cnt, int *ind,
			const EGlpNum_t * val, const EGlpNum_t rhs, int sense,const EGlpNum_t range,
			const char *rowname),
		ILLlib_addrows ( lpinfo * lp, ILLlp_basis * B, int num, int *rmatcnt, 
			int *rmatbeg, int *rmatind,const EGlpNum_t * rmatval,const EGlpNum_t * rhs, 
			char *sense, const EGlpNum_t * range, const char **names, int *nofactor),
		ILLlib_delrows ( lpinfo * lp, ILLlp_basis * B, ILLlp_cache * C, int num, 
			int *dellist, int *basis_ok, int *cache_ok), 
		ILLlib_newcol ( lpinfo * lp, ILLlp_basis * B,const EGlpNum_t obj, 
			const EGlpNum_t lower,const EGlpNum_t upper, const char *name, int factorok), 
		ILLlib_newcols ( lpinfo * lp, ILLlp_basis * B, int num, EGlpNum_t * obj,
			EGlpNum_t * lower, EGlpNum_t * upper, const char **names, int factorok),
		ILLlib_addcol ( lpinfo * lp, ILLlp_basis * B, int cnt, int *ind,
			EGlpNum_t * val,const EGlpNum_t obj,const EGlpNum_t lower,const EGlpNum_t upper, 
			const char *name, int factorok),
		ILLlib_addcols ( lpinfo * lp, ILLlp_basis * B, int num, int *cmatcnt,
			int *cmatbeg, int *cmatind, EGlpNum_t * cmatval, EGlpNum_t * obj,
			EGlpNum_t * lower, EGlpNum_t * upper, const char **names, int factorok),
		ILLlib_delcols ( lpinfo * lp, ILLlp_basis * B, int num, int *dellist,
			int *basis_ok),
		ILLlib_chgcoef ( lpinfo * lp, int rowindex, int colindex, EGlpNum_t coef), 
		ILLlib_getcoef (lpinfo *lp, int rowindex, int colindex, EGlpNum_t* coef),
		ILLlib_chgrange (lpinfo *lp, int indx, EGlpNum_t coef),
		ILLlib_chgsense ( lpinfo * lp, int num, int *rowlist, char *sense),
    ILLlib_getsenses (lpinfo *lp, char *senses),
		ILLlib_getrows ( lpinfo * lp, int num, int *rowlist, int **rowcnt,
			int **rowbeg, int **rowind, EGlpNum_t ** rowval, EGlpNum_t ** rhs,
			char **sense, EGlpNum_t ** range, char ***names),
		ILLlib_getcols ( lpinfo * lp, int num, int *collist, int **colcnt,
			int **colbeg, int **colind, EGlpNum_t ** colval, EGlpNum_t ** obj,
			EGlpNum_t ** lower, EGlpNum_t ** upper, char ***names),
		ILLlib_getobj ( lpinfo * lp, EGlpNum_t * obj),
		ILLlib_getobj_list (lpinfo *lp, int num, int* collist, EGlpNum_t* obj),
		ILLlib_chgobj ( lpinfo * lp, int indx, EGlpNum_t coef),
		ILLlib_getrhs ( lpinfo * lp, EGlpNum_t * rhs),
		ILLlib_chgrhs ( lpinfo * lp, int indx, EGlpNum_t coef),
		ILLlib_getintflags ( lpinfo * lp, int *intflags),
		ILLlib_rownames ( lpinfo * lp, char **rownames),
		ILLlib_colnames ( lpinfo * lp, char **colnames),
		ILLlib_colindex ( lpinfo * lp, const char *name, int *colindex),
		ILLlib_rowindex ( lpinfo * lp, const char *name, int *rowindex),
		ILLlib_chgbnd ( lpinfo * lp, int indx, int lu,const EGlpNum_t bnd),
		ILLlib_chgbnds ( lpinfo * lp, int cnt, int *indx, char *lu, const EGlpNum_t * bnd),
		ILLlib_getbnd ( lpinfo * lp, int indx, int lu, EGlpNum_t * bnd),
		ILLlib_getbnds ( lpinfo * lp, EGlpNum_t * lower, EGlpNum_t * upper),
		ILLlib_getbnds_list ( lpinfo *lp, int num, int*collist, EGlpNum_t *lower,
			EGlpNum_t *upper),
		ILLlib_strongbranch ( lpinfo * lp, price_info * pinf, int *candidatelist,
			int ncand, EGlpNum_t * xlist, EGlpNum_t * downpen, EGlpNum_t * uppen,
			int iterations, EGlpNum_t objbound, struct itcnt_t*itcnt),
		ILLlib_getbasis ( lpinfo * lp, char *cstat, char *rstat),
		ILLlib_loadbasis ( ILLlp_basis * B, int nstruct, int nrows, char *cstat,
			char *rstat),
		ILLlib_readbasis ( lpinfo * lp, ILLlp_basis * B, const char *fname),
		ILLlib_writebasis ( lpinfo * lp, ILLlp_basis * B, const char *fname),
		ILLlib_getrownorms ( lpinfo * lp, price_info * pinf, EGlpNum_t * rownorms),
		ILLlib_loadrownorms ( lpinfo * lp, price_info * pinf, EGlpNum_t * rownorms),
		ILLlib_recompute_rownorms ( lpinfo * lp, price_info * pinf),
		ILLlib_iter ( lpinfo * lp),
		ILLlib_print_x ( EGioFile_t * fd, lpinfo * lp, ILLlp_cache * C, EGlpNum_t * x,
			int nonZerosOnly),
		ILLwrite_lp_file ( ILLlpdata * lp, EGioFile_t * eout, qserror_collector * c);


extern int ILLlib_findName (
	ILLlpdata * qslp,
	int forRow,
	const char *name,
	int id,
	char buf[ILL_namebufsize]);

/****************************************************************************/
/*                                                                          */
/*                           presolve.c                                     */
/*                                                                          */
/****************************************************************************/

int ILLpresolve_add_logicals (
	ILLlpdata * lp);


/****************************************************************************/
/*                                                                          */
/*                            binary.c                                      */
/*                                                                          */
/****************************************************************************/

int ILLmip_binary_dfs (
	lpinfo * lp);

#endif /* ILL_LIB_H */
