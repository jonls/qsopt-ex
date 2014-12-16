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

/* RCSINFO $Id: lpdefs_EGLPNUM_TYPENAME.h,v 1.3 2003/11/05 16:57:39 meven Exp $ */
#ifndef EGLPNUM_TYPENAME___QS_LPDEFS_H
#define EGLPNUM_TYPENAME___QS_LPDEFS_H

#include "urandom.h"

#include "qsopt_EGLPNUM_TYPENAME.h"
#include "lpdata_EGLPNUM_TYPENAME.h"
#include "factor_EGLPNUM_TYPENAME.h"

/* infinity and negative infinity */
#define EGLPNUM_TYPENAME_INFTY  EGLPNUM_TYPENAME_ILL_MAXDOUBLE
#define EGLPNUM_TYPENAME_NINFTY EGLPNUM_TYPENAME_ILL_MINDOUBLE

#include "basicdefs.h"
/* tolerances, these are initialized in EGLPNUM_TYPENAME_ILLstart, file lpdata.c */
/* these three constants are defined in lpdata.c */
extern EGLPNUM_TYPE EGLPNUM_TYPENAME_PARAM_IBASIS_RPIVOT;	/*       0.98 */
extern EGLPNUM_TYPE EGLPNUM_TYPENAME_PARAM_IBASIS_RTRIANG;/*       0.01 */
extern EGLPNUM_TYPE EGLPNUM_TYPENAME_PARAM_MIN_DNORM;			/*      1e-24 */
extern EGLPNUM_TYPE EGLPNUM_TYPENAME_PFEAS_TOLER;					/*       1e-6 */
extern EGLPNUM_TYPE EGLPNUM_TYPENAME_BD_TOLER;						/*       1e-7 */
extern EGLPNUM_TYPE EGLPNUM_TYPENAME_DFEAS_TOLER;					/*       1e-6 */
extern EGLPNUM_TYPE EGLPNUM_TYPENAME_PIVOT_TOLER;					/*      1e-10 */
extern EGLPNUM_TYPE EGLPNUM_TYPENAME_SZERO_TOLER;					/*      1e-15 */
extern EGLPNUM_TYPE EGLPNUM_TYPENAME_PIVZ_TOLER;					/*      1e-12 */
extern EGLPNUM_TYPE EGLPNUM_TYPENAME_OBJBND_TOLER;				/*       1e-2 */
extern EGLPNUM_TYPE EGLPNUM_TYPENAME_DBNDPIV_TOLER;				/*       1e-3 */
extern EGLPNUM_TYPE EGLPNUM_TYPENAME_DBNDPIV_RATIO;				/*       1e-2 */
extern EGLPNUM_TYPE EGLPNUM_TYPENAME_ALTPIV_TOLER;				/*       1e-8 */
//extern EGLPNUM_TYPE DJZERO_TOLER;/*             1e-8 */
extern EGLPNUM_TYPE EGLPNUM_TYPENAME_PROGRESS_ZERO;				/*       1e-7 */
extern EGLPNUM_TYPE EGLPNUM_TYPENAME_PROGRESS_THRESH;			/*       1e-5 */
extern EGLPNUM_TYPE EGLPNUM_TYPENAME_CB_EPS;							/*      0.001 */
extern EGLPNUM_TYPE EGLPNUM_TYPENAME_CB_INF_RATIO;				/*       10.0 */
extern EGLPNUM_TYPE EGLPNUM_TYPENAME_CB_PRI_RLIMIT;				/*       0.25 */

/* structure for statistics */
typedef struct
{
	int ynz_cnt;									/* nz in entering columns */
	int num_y;
	EGLPNUM_TYPE y_ravg;							/* weighted avg. of current & prior y */
	int znz_cnt;									/* nz in ith row of B^{-1}, ie z_i */
	int num_z;
	EGLPNUM_TYPE z_ravg;							/* weighted avg. of current & prior z */
	int zanz_cnt;									/* nz in z^TA */
	int num_za;
	EGLPNUM_TYPE za_ravg;						/* weighted avg. of current & prior za */
	int pnorm_cnt;								/* nz in columns for primal norms */
	int dnorm_cnt;								/* nz in rows for dual norms */
	int pinz_cnt;									/* nz in phase II pi (solve) */
	int num_pi;										/* # of pi solves */
	int pi1nz_cnt;								/* nz in phase I pi (solve) */
	int num_pi1;									/* # of phase I pi solves */
	int upnz_cnt;									/* nz in ftran update vector */
	int num_up;										/* # of ftran_updates */
	int pupv_cnt;									/* nz in primal steep updates */
	int dupv_cnt;									/* nz in dual steep updates */

	int start_slacks;							/* # slacks in beginning */
	int final_slacks;							/* # slacks in the end */
	int start_art;								/* # arts in beginning */
	int final_art;								/* # arts in the end */

	int pI_iter;									/* primal phase I iterations */
	int pII_iter;
	int dI_iter;									/* dual phase I iterations */
	int dII_iter;
	int tot_iter;

	int pivpI[10];								/* sizes of pivots */
	int pivpII[10];
	int pivdI[10];
	int pivdII[10];
}
EGLPNUM_TYPENAME_count_struct;

/* structure for tolerances */
typedef struct
{
	EGLPNUM_TYPE pfeas_tol;
	EGLPNUM_TYPE dfeas_tol;
	EGLPNUM_TYPE pivot_tol;
	EGLPNUM_TYPE szero_tol;
	EGLPNUM_TYPE ip_tol;							/* inner primal & dual feas toler */
	EGLPNUM_TYPE id_tol;
}
EGLPNUM_TYPENAME_tol_struct;

/* bound information */
typedef struct EGLPNUM_TYPENAME_bndinfo
{
	EGLPNUM_TYPE pbound;
	EGLPNUM_TYPE cbound;
	int btype;
	int varnum;
	struct EGLPNUM_TYPENAME_bndinfo *next;
}
EGLPNUM_TYPENAME_bndinfo;

/* bound information */
typedef struct EGLPNUM_TYPENAME_coefinfo
{
	EGLPNUM_TYPE pcoef;
	EGLPNUM_TYPE ccoef;
	int varnum;
	struct EGLPNUM_TYPENAME_coefinfo *next;
}
EGLPNUM_TYPENAME_coefinfo;

/* feasibility info */
typedef struct EGLPNUM_TYPENAME_feas_info
{
	int pstatus;
	int dstatus;
	EGLPNUM_TYPE totinfeas;
}
EGLPNUM_TYPENAME_feas_info;

typedef struct EGLPNUM_TYPENAME_lp_status_info
{
	char optimal;
	char primal_feasible;
	char primal_infeasible;
	char primal_unbounded;
	char dual_feasible;
	char dual_infeasible;
	char dual_unbounded;
	char padd;
}
EGLPNUM_TYPENAME_lp_status_info;

typedef struct EGLPNUM_TYPENAME_pI_uinfo
{
	int tctr;
	int i;
	int *perm;
	int *ix;
	int fs;
	EGLPNUM_TYPE piv;
	EGLPNUM_TYPE *t;
	EGLPNUM_TYPE dty;
	EGLPNUM_TYPE c_obj;
	EGLPNUM_TYPE tz;
}
EGLPNUM_TYPENAME_pI_uinfo;

extern void EGLPNUM_TYPENAME_ILLlp_status_info_init (
	EGLPNUM_TYPENAME_lp_status_info * ls);

/* structure for local lp information
 * contains lp obj values - status - dimensions - input data -
 * solution vecs - basis info - update vecs - work vecs - bound changes -
 * tolerances - time info - statistics 
 */
typedef struct EGLPNUM_TYPENAME_lpinfo
{

	EGLPNUM_TYPE objval;							/* obj info */
	EGLPNUM_TYPE pobjval;						/* intermediate status info */
	EGLPNUM_TYPE dobjval;
	EGLPNUM_TYPE pinfeas;
	EGLPNUM_TYPE dinfeas;
	EGLPNUM_TYPE objbound;
	EGLPNUM_TYPENAME_lp_status_info probstat;			/* final status */
	EGLPNUM_TYPENAME_lp_status_info basisstat;			/* final status */
	int nrows;										/* input info follows; given in col format */
	int ncols;
	int *matcnt;
	int *matbeg;
	int *matind;
	EGLPNUM_TYPE *matval;
	int matfree;
	int matsize;
	EGLPNUM_TYPE *bz;
	EGLPNUM_TYPE *lz;
	EGLPNUM_TYPE *uz;
	EGLPNUM_TYPE *cz;
	int localrows;								/* set to 1 if these are created locally */
	int *rowcnt;									/* row info follows, copy of col info */
	int *rowbeg;
	int *rowind;
	EGLPNUM_TYPE *rowval;

	EGLPNUM_TYPE *xbz;								/* output info x, pi, reduced cost */
	EGLPNUM_TYPE *piz;
	EGLPNUM_TYPE *dz;
	EGLPNUM_TYPE *pIxbz;							/* output info (phase I) x, pi, reduced cost */
	EGLPNUM_TYPE *pIpiz;
	EGLPNUM_TYPE *pIdz;

	int final_phase;							/* final phase, inf & unboundedness info */
	int infub_ix;

	int basisid;									/* basis and variable info follows */
	int nnbasic;
	int *baz;
	int *nbaz;
	int *vstat;
	int *vindex;
	int fbasisid;
	EGLPNUM_TYPENAME_factor_work *f;
	int *vtype;										/* internal var info */
	char *vclass;									/* structural or logical */

	EGLPNUM_TYPENAME_svector zz;										/* local EGLPNUM_TYPENAME_ILLfactor_update vectors z, yj, za */
	EGLPNUM_TYPENAME_svector yjz;
	EGLPNUM_TYPENAME_svector zA;
	EGLPNUM_TYPENAME_svector work;									/* local work vector */
	EGLPNUM_TYPENAME_svector srhs;									/* local vectors for lin. eq. solves */
	EGLPNUM_TYPENAME_svector ssoln;
	int *iwork;										/* local work vector */
	EGLPNUM_TYPENAME_pI_uinfo upd;									/* phase I update info */
	int *bfeas;										/* primal and dual infeasibility info */
	int *dfeas;

	EGLPNUM_TYPENAME_tol_struct *tol;							/* tolerances */
	EGLPNUM_TYPENAME_count_struct *cnts;						/* counts */
	int nbchange;									/* # bound shifts */
	int ncchange;									/* # obj coef shifts */
	EGLPNUM_TYPENAME_bndinfo *bchanges;						/* list of bound shifts */
	EGLPNUM_TYPENAME_coefinfo *cchanges;						/* list of coef shifts */
	int pIratio;									/* ratio tests */
	int pIIratio;
	int dIratio;
	int dIIratio;

	int maxiter;
	int iterskip;
	double maxtime;
	double starttime;
	struct EGLPNUM_TYPENAME_ILLlpdata *O;
	ILLrandstate rstate;

}
EGLPNUM_TYPENAME_lpinfo;

/* pricing structures */
typedef struct
{
	int ninit;
	EGLPNUM_TYPE *norms;
	int *refframe;
}
EGLPNUM_TYPENAME_p_devex_info;

typedef struct
{
	EGLPNUM_TYPE *norms;
}
EGLPNUM_TYPENAME_p_steep_info;

typedef struct
{
	int k;
	int cgroup;
	int ngroups;
	int *gstart;
	int *gshift;
	int *gsize;
	int bsize;
	int *bucket;
	int *perm;
	EGLPNUM_TYPE *infeas;
}
EGLPNUM_TYPENAME_mpart_info;

typedef struct
{
	int ninit;
	EGLPNUM_TYPE *norms;
	int *refframe;
}
EGLPNUM_TYPENAME_d_devex_info;

typedef struct
{
	EGLPNUM_TYPE *norms;
}
EGLPNUM_TYPENAME_d_steep_info;

/* pricing information */
typedef struct EGLPNUM_TYPENAME_price_info
{
	int p_strategy;
	int d_strategy;
	int pI_price;
	int pII_price;
	int dI_price;
	int dII_price;
	int cur_price;
	EGLPNUM_TYPE *p_scaleinf;
	EGLPNUM_TYPE *d_scaleinf;
	EGLPNUM_TYPENAME_p_devex_info pdinfo;
	EGLPNUM_TYPENAME_p_steep_info psinfo;
	EGLPNUM_TYPENAME_mpart_info pmpinfo;
	EGLPNUM_TYPENAME_d_devex_info ddinfo;
	EGLPNUM_TYPENAME_d_steep_info dsinfo;
	EGLPNUM_TYPENAME_mpart_info dmpinfo;
	EGLPNUM_TYPENAME_heap h;
	EGLPNUM_TYPE htrigger;
	int hineff;
	char init;
}
EGLPNUM_TYPENAME_price_info;

#endif /* EGLPNUM_TYPENAME___QS_LPDEFS_H */
