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

/* RCSINFO $Id: lpdefs.h,v 1.3 2003/11/05 16:57:39 meven Exp $ */
#ifndef __QS_LPDEFS_H
#define __QS_LPDEFS_H

#include "urandom.h"

#include "qsopt.h"
#include "lpdata.h"
#include "factor.h"

/* infinity and negative infinity */
#define INFTY  ILL_MAXDOUBLE
#define NINFTY ILL_MINDOUBLE

#include "basicdefs.h"
/* tolerances, these are initialized in ILLstart, file lpdata.c */
//#if EGLPNUM_TYPE != DBL_TYPE && EGLPNUM_TYPE != LDBL_TYPE
/* these three constants are defined in lpdata.c */
extern EGlpNum_t PARAM_IBASIS_RPIVOT;	/*       0.98 */
extern EGlpNum_t PARAM_IBASIS_RTRIANG;/*       0.01 */
extern EGlpNum_t PARAM_MIN_DNORM;			/*      1e-24 */
extern EGlpNum_t PFEAS_TOLER;					/*       1e-6 */
extern EGlpNum_t BD_TOLER;						/*       1e-7 */
extern EGlpNum_t DFEAS_TOLER;					/*       1e-6 */
extern EGlpNum_t PIVOT_TOLER;					/*      1e-10 */
extern EGlpNum_t SZERO_TOLER;					/*      1e-15 */
extern EGlpNum_t PIVZ_TOLER;					/*      1e-12 */
extern EGlpNum_t OBJBND_TOLER;				/*       1e-2 */
extern EGlpNum_t DBNDPIV_TOLER;				/*       1e-3 */
extern EGlpNum_t DBNDPIV_RATIO;				/*       1e-2 */
extern EGlpNum_t ALTPIV_TOLER;				/*       1e-8 */
//extern EGlpNum_t DJZERO_TOLER;/*             1e-8 */
extern EGlpNum_t PROGRESS_ZERO;				/*       1e-7 */
extern EGlpNum_t PROGRESS_THRESH;			/*       1e-5 */
extern EGlpNum_t CB_EPS;							/*      0.001 */
extern EGlpNum_t CB_INF_RATIO;				/*       10.0 */
extern EGlpNum_t CB_PRI_RLIMIT;				/*       0.25 */

/* structure for statistics */
typedef struct
{
	int ynz_cnt;									/* nz in entering columns */
	int num_y;
	EGlpNum_t y_ravg;							/* weighted avg. of current & prior y */
	int znz_cnt;									/* nz in ith row of B^{-1}, ie z_i */
	int num_z;
	EGlpNum_t z_ravg;							/* weighted avg. of current & prior z */
	int zanz_cnt;									/* nz in z^TA */
	int num_za;
	EGlpNum_t za_ravg;						/* weighted avg. of current & prior za */
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
count_struct;

/* structure for tolerances */
typedef struct
{
	EGlpNum_t pfeas_tol;
	EGlpNum_t dfeas_tol;
	EGlpNum_t pivot_tol;
	EGlpNum_t szero_tol;
	EGlpNum_t ip_tol;							/* inner primal & dual feas toler */
	EGlpNum_t id_tol;
}
tol_struct;

/* bound information */
typedef struct bndinfo
{
	EGlpNum_t pbound;
	EGlpNum_t cbound;
	int btype;
	int varnum;
	struct bndinfo *next;
}
bndinfo;

/* bound information */
typedef struct coefinfo
{
	EGlpNum_t pcoef;
	EGlpNum_t ccoef;
	int varnum;
	struct coefinfo *next;
}
coefinfo;

/* feasibility info */
typedef struct feas_info
{
	int pstatus;
	int dstatus;
	EGlpNum_t totinfeas;
}
feas_info;

typedef struct lp_status_info
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
lp_status_info;

typedef struct pI_uinfo
{
	int tctr;
	int i;
	int *perm;
	int *ix;
	int fs;
	EGlpNum_t piv;
	EGlpNum_t *t;
	EGlpNum_t dty;
	EGlpNum_t c_obj;
	EGlpNum_t tz;
}
pI_uinfo;

extern void ILLlp_status_info_init (
	lp_status_info * ls);

/* structure for local lp information
 * contains lp obj values - status - dimensions - input data -
 * solution vecs - basis info - update vecs - work vecs - bound changes -
 * tolerances - time info - statistics 
 */
typedef struct lpinfo
{

	EGlpNum_t objval;							/* obj info */
	EGlpNum_t pobjval;						/* intermediate status info */
	EGlpNum_t dobjval;
	EGlpNum_t pinfeas;
	EGlpNum_t dinfeas;
	EGlpNum_t objbound;
	lp_status_info probstat;			/* final status */
	lp_status_info basisstat;			/* final status */
	int nrows;										/* input info follows; given in col format */
	int ncols;
	int *matcnt;
	int *matbeg;
	int *matind;
	EGlpNum_t *matval;
	int matfree;
	int matsize;
	EGlpNum_t *bz;
	EGlpNum_t *lz;
	EGlpNum_t *uz;
	EGlpNum_t *cz;
	int localrows;								/* set to 1 if these are created locally */
	int *rowcnt;									/* row info follows, copy of col info */
	int *rowbeg;
	int *rowind;
	EGlpNum_t *rowval;

	EGlpNum_t *xbz;								/* output info x, pi, reduced cost */
	EGlpNum_t *piz;
	EGlpNum_t *dz;
	EGlpNum_t *pIxbz;							/* output info (phase I) x, pi, reduced cost */
	EGlpNum_t *pIpiz;
	EGlpNum_t *pIdz;

	int final_phase;							/* final phase, inf & unboundedness info */
	int infub_ix;

	int basisid;									/* basis and variable info follows */
	int nnbasic;
	int *baz;
	int *nbaz;
	int *vstat;
	int *vindex;
	int fbasisid;
	factor_work *f;
	int *vtype;										/* internal var info */
	char *vclass;									/* structural or logical */

	svector zz;										/* local ILLfactor_update vectors z, yj, za */
	svector yjz;
	svector zA;
	svector work;									/* local work vector */
	svector srhs;									/* local vectors for lin. eq. solves */
	svector ssoln;
	int *iwork;										/* local work vector */
	pI_uinfo upd;									/* phase I update info */
	int *bfeas;										/* primal and dual infeasibility info */
	int *dfeas;

	tol_struct *tol;							/* tolerances */
	count_struct *cnts;						/* counts */
	int nbchange;									/* # bound shifts */
	int ncchange;									/* # obj coef shifts */
	bndinfo *bchanges;						/* list of bound shifts */
	coefinfo *cchanges;						/* list of coef shifts */
	int pIratio;									/* ratio tests */
	int pIIratio;
	int dIratio;
	int dIIratio;

	int maxiter;
	int iterskip;
	double maxtime;
	double starttime;
	struct ILLlpdata *O;
	ILLrandstate rstate;

}
lpinfo;

/* pricing structures */
typedef struct
{
	int ninit;
	EGlpNum_t *norms;
	int *refframe;
}
p_devex_info;

typedef struct
{
	EGlpNum_t *norms;
}
p_steep_info;

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
	EGlpNum_t *infeas;
}
mpart_info;

typedef struct
{
	int ninit;
	EGlpNum_t *norms;
	int *refframe;
}
d_devex_info;

typedef struct
{
	EGlpNum_t *norms;
}
d_steep_info;

/* pricing information */
typedef struct price_info
{
	int p_strategy;
	int d_strategy;
	int pI_price;
	int pII_price;
	int dI_price;
	int dII_price;
	int cur_price;
	EGlpNum_t *p_scaleinf;
	EGlpNum_t *d_scaleinf;
	p_devex_info pdinfo;
	p_steep_info psinfo;
	mpart_info pmpinfo;
	d_devex_info ddinfo;
	d_steep_info dsinfo;
	mpart_info dmpinfo;
	heap h;
	EGlpNum_t htrigger;
	int hineff;
	char init;
}
price_info;

#endif /* __QS_LPDEFS_H */
