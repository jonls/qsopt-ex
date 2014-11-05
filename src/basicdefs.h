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
#ifndef __BASICDEFS__
#define __BASICDEFS__

/* storage type */
#define DENSE          0
#define SPARSE         1

/* type of vector */
#define ROW_SOLVE      1
#define COLUMN_SOLVE   2

/* direction of change in non-basic var */
#define VINCREASE      1
#define VDECREASE      2

/* status of variables */
#define STAT_BASIC     1
#define STAT_UPPER     2
#define STAT_LOWER     3
#define STAT_ZERO      4

#define BOUND_LOWER    1
#define BOUND_UPPER    2

/* type of variables */
#define VARTIFICIAL    1
#define VFIXED         2
#define VFREE          4
#define VUPPER         8
#define VLOWER         16
#define VBOUNDED       32

/* class of variables */
#define CLASS_STRUCT   0
#define CLASS_LOGICAL  1

/* algo */
#define PRIMAL_SIMPLEX 1
#define DUAL_SIMPLEX   2
#define PRIMAL_OR_DUAL 3

/* phase */
#define PRIMAL_PHASEI  1
#define PRIMAL_PHASEII 2
#define DUAL_PHASEI    3
#define DUAL_PHASEII   4

/* number of phases */
#define PHASEI         1
#define PHASEII        2

/* type of pricing (all vars or some) */
#define COMPLETE_PRICING   1
#define PARTIAL_PRICING    2
#define MULTI_PART_PRICING 3

/* default pricing */

#define QS_DEFAULT_PRICE_PI  QS_PRICE_PSTEEP
#define QS_DEFAULT_PRICE_PII QS_PRICE_PSTEEP
#define QS_DEFAULT_PRICE_DI  QS_PRICE_DSTEEP
#define QS_DEFAULT_PRICE_DII QS_PRICE_DSTEEP

/* lp sol status */
#define ILL_LP_SOLVED            1
#define ILL_LP_UNSOLVED          2
#define ILL_MAX_ITER             3
#define ILL_MAX_TIME             4
#define ILL_BND_REACHED          5
#define ILL_PPHASEI_ERROR        6
#define ILL_PPHASEII_ERROR       7
#define ILL_DPHASEI_ERROR        8
#define ILL_DPHASEII_ERROR       9
#define ILL_LP_ABORTED          10

/* basis status */
#define OPTIMAL                  1
#define NONOPTIMAL               2
#define PRIMAL_FEASIBLE          3
#define PRIMAL_INFEASIBLE        4
#define PRIMAL_UNBOUNDED         5
#define DUAL_FEASIBLE            7
#define DUAL_INFEASIBLE          8
#define DUAL_UNBOUNDED           9

/* type of ratio test */
#define RATIOTEST_NORMAL 1
#define RATIOTEST_HARRIS 2

/* control parameters */
#define PARAM_PRATIOTESTS        10
#define PARAM_DRATIOTESTS        20
#define PARAM_PRIMAL_REFACTORGAP 50
#define PARAM_PRIMAL_RESOLVEGAP  25
#define PARAM_DUAL_REFACTORGAP   100
#define PARAM_DUAL_RESOLVEGAP    25
#define PARAM_MAX_NOSOLVE        500
#define PARAM_MAX_NOPROG         300
#define PARAM_NOPROG_FACTOR      15

/* numerical parameters */
#define PARAM_BSHIFT             10
#define PARAM_CSHIFT             10

/* general constants */
#define PARAM_HEAP_UTRIGGER      10
#define PARAM_HEAP_RATIO         4.0

/* errors */
#define E_GENERAL_ERROR          1
#define E_INV_LINSOLVE_OPTION    2
#define E_NO_MEMORY              3
#define E_INVALID_OPTION         4
#define E_NULL_ARGUMENT          5
#define E_SIMPLEX_ERROR          6
#define E_BASIS_SINGULAR         7

#ifndef __QS_BASIS__
#define __QS_BASIS__
typedef struct qsbasis
{
	int nstruct;
	int nrows;
	char *cstat;
	char *rstat;
}
QSbasis;
#endif

typedef struct itcnt_t
{
	int pI_iter;
	int pII_iter;
	int dI_iter;
	int dII_iter;
	int tot_iter;
} itcnt_t;

#ifndef QS_DEFINITIONS
#define QS_DEFINITIONS
#define QS_MIN       (1)
#define QS_MAX       (-1)

/****************************************************************************/
/*                                                                          */
/*                 PARAMETERS THAT CAN BE SET BY setparam                   */
/*                                                                          */
/****************************************************************************/


#define QS_PARAM_PRIMAL_PRICING    0
#define QS_PARAM_DUAL_PRICING      2
#define QS_PARAM_SIMPLEX_DISPLAY   4
#define QS_PARAM_SIMPLEX_MAX_ITERATIONS 5
#define QS_PARAM_SIMPLEX_MAX_TIME  6
#define QS_PARAM_SIMPLEX_SCALING   7
#define QS_PARAM_OBJULIM           8
#define QS_PARAM_OBJLLIM           9


/****************************************************************************/
/*                                                                          */
/*                     VALUES FOR PRICING PARAMETERS                        */
/*                                                                          */
/****************************************************************************/

#define QS_PRICE_PDANTZIG 1
#define QS_PRICE_PDEVEX 2
#define QS_PRICE_PSTEEP 3
#define QS_PRICE_PMULTPARTIAL 4

#define QS_PRICE_DDANTZIG 6
#define QS_PRICE_DSTEEP 7
#define QS_PRICE_DMULTPARTIAL 8
#define QS_PRICE_DDEVEX 9


/****************************************************************************/
/*                                                                          */
/*                         VALUES FOR BASIS STATUS                          */
/*                                                                          */
/****************************************************************************/


#define QS_COL_BSTAT_LOWER     '0'
#define QS_COL_BSTAT_BASIC     '1'
#define QS_COL_BSTAT_UPPER     '2'
#define QS_COL_BSTAT_FREE      '3'

#define QS_ROW_BSTAT_LOWER     '0'
#define QS_ROW_BSTAT_BASIC     '1'
#define QS_ROW_BSTAT_UPPER     '2'


/****************************************************************************/
/*                                                                          */
/*         Return Status for dbl_QSopt_primal, dbl_QSopt_dual, dbl_QSget_status         */
/*                                                                          */
/****************************************************************************/

#define QS_LP_OPTIMAL           1
#define QS_LP_INFEASIBLE        2
#define QS_LP_UNBOUNDED         3
#define QS_LP_ITER_LIMIT        4
#define QS_LP_TIME_LIMIT        5
#define QS_LP_UNSOLVED          6
#define QS_LP_ABORTED		7
#define QS_LP_NUMERR            8
#define QS_LP_OBJ_LIMIT         9
#define QS_LP_MODIFIED        100
#define QS_LP_CHANGE_PREC			1024
#endif



/** @brief If set to one, them we allow to re-start the simplex algorithm due to
 * numerical issues */
#define DO_NUMER 0
/** @brief If set to one, then we allow to re-start simplex due to singular
 * basis */
#define DO_SINGULAR 0

/** @brief Factor for wich we change tolerances each time we have to resume
 * simplex */
#define SIMPLEX_FACTOR 5U
#define DENSE_PI 0
#define DENSE_PIIPI 0
#define DENSE_NORM 0
#define SIMPLEX_DEBUG 0


/* possible values of nextstep */
#define SIMPLEX_CONTINUE   1
#define SIMPLEX_TERMINATE  2
#define SIMPLEX_RESUME     3

/* reason for resuming simplex */
#define SIMPLEX_RESUME_SING     1
#define SIMPLEX_RESUME_UNSHIFT  2
#define SIMPLEX_RESUME_NUMER    3

/* values for newphase */
#define SIMPLEX_PHASE_RECOMP  1
#define SIMPLEX_PHASE_NEW     2

#define SIMPLEX_PIVOTINROW 1
#define SIMPLEX_PIVOTINCOL 2
#define SIMPLEX_MAX_RESTART 4
#define SIMPLEX_MAX_PIVOT_FAIL 300


#define FALSE 0
#define TRUE  1
#define QS_FACTOR_MAX_K        1
#define QS_FACTOR_P            2
#define QS_FACTOR_ETAMAX       3
#define QS_FACTOR_FZERO_TOL    4
#define QS_FACTOR_SZERO_TOL    5
#define QS_FACTOR_PARTIAL_TOL  6
#define QS_FACTOR_UR_SPACE_MUL 7
#define QS_FACTOR_UC_SPACE_MUL 8
#define QS_FACTOR_LC_SPACE_MUL 9
#define QS_FACTOR_LR_SPACE_MUL 10
#define QS_FACTOR_ER_SPACE_MUL 11
#define QS_FACTOR_GROW_MUL     12
#define QS_FACTOR_MAXMULT      13
#define QS_FACTOR_MINMULT      14
#define QS_FACTOR_UPDMAXMULT   15
#define QS_FACTOR_DENSE_FRACT  16
#define QS_FACTOR_DENSE_MIN    17
#define E_CHECK_FAILED 6
#define E_NO_PIVOT 7
#define E_FACTOR_BLOWUP 8
#define E_UPDATE_NOSPACE 9
#define E_UPDATE_SINGULAR_ROW 10
#define E_UPDATE_SINGULAR_COL 11
#define E_SING_NO_DATA 12
#define E_SINGULAR_INTERNAL 13
#define SPARSE_FACTOR 0.05
#define CNT_YNZ           1			/* nz in entering columns */
#define CNT_ZNZ           2			/* nz in ith row of B^{-1}, ie z_i */
#define CNT_ZANZ          3			/* nz in ith row of B^{-1}, ie z_i */
#define CNT_PINZ          4			/* nz in phase II pi (solve) */
#define CNT_P1PINZ        5			/* nz in phase I pi (solve) */
#define CNT_UPNZ          6			/* nz in ftran_updates */
#define CNT_PPHASE1ITER   7			/* primal phase I iterations */
#define CNT_PPHASE2ITER   8
#define CNT_DPHASE1ITER   9			/* dual phase I iterations */
#define CNT_DPHASE2ITER   10
#define CNT_PIPIV         11
#define CNT_PIIPIV        12
#define CNT_DIPIV         13
#define CNT_DIIPIV        14
#define CNT_YRAVG         15
#define CNT_ZARAVG        16

#define ROW_PIVOT         0
#define COL_PIVOT         1

#define ILL_LP_OPTIMAL           1
#define ILL_LP_NONOPTIMAL        2
#define ILL_LP_PRIMAL_FEASIBLE   3
#define ILL_LP_PRIMAL_INFEASIBLE 4
#define ILL_LP_PRIMAL_UNBOUNDED  5
#define ILL_LP_DUAL_FEASIBLE     6
#define ILL_LP_DUAL_INFEASIBLE   7
#define ILL_LP_DUAL_UNBOUNDED    8

typedef enum
{ ILL_MPS_NAME, ILL_MPS_OBJSENSE, ILL_MPS_OBJNAME,
	ILL_MPS_ROWS, ILL_MPS_COLS, ILL_MPS_RHS, ILL_MPS_RANGES,
	ILL_MPS_BOUNDS, ILL_MPS_REFROW, ILL_MPS_ENDATA,
	ILL_MPS_NONE
}
ILLmps_section;

#define ILL_MPS_N_SECTIONS ILL_MPS_NONE

/*define as > 0 if heap is to be used */
#define USEHEAP 1

/*result of pricing */
#define PRICE_OPTIMAL 1
#define PRICE_NONOPTIMAL  2

/*type of pricing */
#define ROW_PRICING 1
#define COL_PRICING 2


/****************************************************************************
 * error collection
 */
#define QS_DATA_ERROR			0
#define QS_DATA_WARN			1
#define QS_MPS_FORMAT_ERROR		2
#define QS_MPS_FORMAT_WARN		3
#define QS_LP_FORMAT_ERROR		4
#define QS_LP_FORMAT_WARN		5
#define QS_INPUT_NERROR        8

/* defs for phase I ratio test */
#define BBOUND    1
#define BATOLOWER 2
#define BATOUPPER 3
#define BBTOLOWER 4
#define BBTOUPPER 5
#define BSKIP     6

/* result of ratio test */
#define RATIO_UNBOUNDED 1
#define RATIO_NOBCHANGE 2
#define RATIO_BCHANGE   3
#define RATIO_FAILED    4
#define RATIO_NEGATIVE  5

/* warning level */
#define QSE_WLVL 10000






#endif
