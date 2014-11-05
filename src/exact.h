/* ========================================================================= */
/* ESolver "Exact Mixed Integer Linear Solver" provides some basic structures 
 * and algorithms commons in solving MIP's
 *
 * Copyright (C) 2005 Daniel Espinoza.
 * 
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation; either version 2.1 of the License, or (at your
 * option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public 
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA 
 * */
/* ========================================================================= */
#ifndef __EXACT_H__
#define __EXACT_H__
#include "qs_config.h"

#include "eg_fp.h"
#include "eg_lpnum.h"
#include "eg_io.h"
#include "eg_timer.h"

#include "symtab.h"

#include "dbl_basis.h"
#include "dbl_dheaps_i.h"
#include "dbl_dstruct.h"
#include "dbl_factor.h"
#include "dbl_format.h"
#include "dbl_lpdata.h"
#include "dbl_lpdefs.h"
#include "dbl_mps.h"
#include "dbl_price.h"
#include "dbl_priority.h"
#include "dbl_qsopt.h"
#include "dbl_qstruct.h"
#include "dbl_ratio.h"
#include "dbl_rawlp.h"
#include "dbl_readline.h"
#include "dbl_read_lp.h"
#include "dbl_read_mps.h"
#include "dbl_simplex.h"
#include "dbl_write_lp.h"
#include "dbl_lib.h"
#include "dbl_editor.h"

#include "mpq_lpdata.h"
#include "mpq_lpdefs.h"
#include "mpq_basis.h"
#include "mpq_dheaps_i.h"
#include "mpq_dstruct.h"
#include "mpq_factor.h"
#include "mpq_format.h"
#include "mpq_mps.h"
#include "mpq_price.h"
#include "mpq_priority.h"
#include "mpq_qsopt.h"
#include "mpq_qstruct.h"
#include "mpq_ratio.h"
#include "mpq_rawlp.h"
#include "mpq_readline.h"
#include "mpq_read_lp.h"
#include "mpq_read_mps.h"
#include "mpq_simplex.h"
#include "mpq_write_lp.h"
#include "mpq_lib.h"
#include "mpq_editor.h"

#include "fp20_basis.h"
#include "fp20_dheaps_i.h"
#include "fp20_dstruct.h"
#include "fp20_factor.h"
#include "fp20_format.h"
#include "fp20_lpdata.h"
#include "fp20_lpdefs.h"
#include "fp20_mps.h"
#include "fp20_price.h"
#include "fp20_priority.h"
#include "fp20_qsopt.h"
#include "fp20_qstruct.h"
#include "fp20_ratio.h"
#include "fp20_rawlp.h"
#include "fp20_readline.h"
#include "fp20_read_lp.h"
#include "fp20_read_mps.h"
#include "fp20_simplex.h"
#include "fp20_write_lp.h"
#include "fp20_lib.h"
#include "fp20_editor.h"

#include "mpf_basis.h"
#include "mpf_dheaps_i.h"
#include "mpf_dstruct.h"
#include "mpf_factor.h"
#include "mpf_format.h"
#include "mpf_lpdata.h"
#include "mpf_lpdefs.h"
#include "mpf_mps.h"
#include "mpf_price.h"
#include "mpf_priority.h"
#include "mpf_qsopt.h"
#include "mpf_qstruct.h"
#include "mpf_ratio.h"
#include "mpf_rawlp.h"
#include "mpf_readline.h"
#include "mpf_read_lp.h"
#include "mpf_read_mps.h"
#include "mpf_simplex.h"
#include "mpf_write_lp.h"
#include "mpf_lib.h"
#include "mpf_editor.h"

#if ENABLE_LONG_DOUBLE
#include "ldbl_basis.h"
#include "ldbl_dheaps_i.h"
#include "ldbl_dstruct.h"
#include "ldbl_factor.h"
#include "ldbl_format.h"
#include "ldbl_lpdata.h"
#include "ldbl_lpdefs.h"
#include "ldbl_mps.h"
#include "ldbl_price.h"
#include "ldbl_priority.h"
#include "ldbl_qsopt.h"
#include "ldbl_qstruct.h"
#include "ldbl_ratio.h"
#include "ldbl_rawlp.h"
#include "ldbl_readline.h"
#include "ldbl_read_lp.h"
#include "ldbl_read_mps.h"
#include "ldbl_simplex.h"
#include "ldbl_write_lp.h"
#include "ldbl_lib.h"
#include "ldbl_editor.h"
#endif

#ifdef HAVE_SOFTFLOAT
#include "float128_basis.h"
#include "float128_dheaps_i.h"
#include "float128_dstruct.h"
#include "float128_factor.h"
#include "float128_format.h"
#include "float128_lpdata.h"
#include "float128_lpdefs.h"
#include "float128_mps.h"
#include "float128_price.h"
#include "float128_priority.h"
#include "float128_qsopt.h"
#include "float128_qstruct.h"
#include "float128_ratio.h"
#include "float128_rawlp.h"
#include "float128_readline.h"
#include "float128_read_lp.h"
#include "float128_read_mps.h"
#include "float128_simplex.h"
#include "float128_write_lp.h"
#include "float128_lib.h"
#include "float128_editor.h"
#endif

#include "eg_exutil.h"

/* ========================================================================= */
/** @defgroup Esolver Esolver
 * Here we define an interface to solve LP's (#QSexact_solver) and MIP's 
 * exactly. 
 * @par History:
 * Revision 0.1
 * - 2005-11-14
 * 						- Fix handling of infeasibility testing, the problem is that
 * 						sometimes, the QSget_infeas_array only work after calling primal
 * 						simplex, so, if we are doing dual, ans finish with infeasibility
 * 						status, the call will fail, the fix is to call primal simples on
 * 						those cases before calling the infeasibility proof.
 * - 2005-10-05
 * 						- If one of the floating point approximations fail, keep going
 * 						to the next floating point approximation. A floating point
 * 						approximation may fail because the basis is singular within the
 * 						used precision.
 * - 2005-09-29
 * 						- If the plain double approximation return unbounded we re-try
 * 						in extended precision, but when the extended precision LP solver
 * 						return with QS_LP_UNBOUNDED status, we just give-up and return
 * 						QS_LP_UNBOUNDED status.
 * - 2005-08-17
 * 						- Improve reliability of optimality test.
 * - 2005-07-07
 * 						- Load optimal soplution into the cache.
 * 						- Change the behavior of QSopt, when he wants to re-start simplex,
 * 							we instead increase the precision of the numbers, this behavior
 * 							is managed by DO_NUMER and DO_SINGULAR.
 * - 2005-05-31
 * 						- If the status of the ending call is not optimal, we don't load
 * 						the previous basis, this is because in some examples doing so lead
 * 						to bad behavior of the overall code.
 * - 2005-05-11
 * 						- First definition and implementation
 *
 * */
/** @file
 * @ingroup Esolver */
/** @addtogroup Esolver */
/** @{ */
/* ========================================================================= */

/* ========================================================================= */
/** @brief If enabled, save the last problem proved to be optimal, and its
 * solution. */
#define QSEXACT_SAVE_OPTIMAL 0

/* ========================================================================= */
/** @brief If enabled, save the intermediate problems created by the functions
 * #QScopy_prob_mpq_dbl and #QScopy_prob_mpq_mpf */
#define QSEXACT_SAVE_INT 0

/* ========================================================================= */
/** @brief Copy an exact problem (mpq_QSdata) to a regular double version of the
 * problem (dbl_QSdata) */
dbl_QSdata *QScopy_prob_mpq_dbl (mpq_QSdata * p,
																 const char *newname);

/* ========================================================================= */
/** @brief Copy an exact problem (mpq_QSdata) to a regular double version of the
 * problem (dbl_QSdata) */
mpf_QSdata *QScopy_prob_mpq_mpf (mpq_QSdata * p,
																 const char *newname);

/* ========================================================================= */
/** @brief Test if a given primal/dual solution is feasible and has the same
 * objective value.
 * @param p original problem.
 * @param p_sol primal solution candidate.
 * @param d_sol dual solution candidate.
 * @param basis Basis for wich the current primal/dual vector is a solution.
 * @return one if the given primal/dual solution is optimal, zero otherwise. 
 * @par Description:
 * The input problem has the form \f[ \begin{array}{l}\min cx\\
  s.t. \begin{array}{lcl}Ax&=&b\end{array}\\
  l\leq x\leq u\end{array} \f]
 * where some of the bounds can be \f$\infty\f$ or \f$-\infty\f$. Note that from
 * this the dual problem is allways feasible (we treat \f$\infty\f$ as a
 * suitable large number) because it looks like 
 * \f[ \begin{array}{l}\max by + d_uu-d_ll\\
  s.t. \begin{array}{lcl}A^ty-Id_l+Id_u & =& c \end{array}\\
  d_u,d_l\geq0\end{array} \f] thus we just need to check primal 
 * feasibility and complementary slackness (just to be sure we also check that
 * both dual and primal objective values coincide.
 *
 * If the optimality test is true (i.e. the basis and the given solution, wich
 * might have been modified inside the function) then this function store the
 * optimal solution into the cache of the problem.
 *
 * @note We assume that p_sol and d_sol have the right size for the problem. and
 * moreover, we assume that the problem already has the logical variables added
 * in (to transform it into standard form), this allow us to fix somewhat the
 * primal vector solution to try to get an optimality certificate.
 * */
int QSexact_optimal_test (mpq_QSdata * p,
													mpq_t * p_sol,
													mpq_t * d_sol,
													QSbasis * basis);

/* ========================================================================= */
/** @brief Print into a file the optimal solution.
 * @param p original problem.
 * @param out_f file where to write the solution.
 * @return zero on success, non-zero otherwise.
 * */
int QSexact_print_sol (mpq_QSdata * p,
											 EGioFile_t * out_f);

/* ========================================================================= */
/** @brief Check if the given dual vector is a proof of infeasibility for the
 * given exact problem. 
 * @param p pointer to the problem data structure.
 * @param d_sol array of length at least nrows with the suposed proof of
 * infeasibility.
 * @return zero if the given dual vector is indeed a proof of infeasibility for
 * the problem, non zero otherwise.
 * @par Description:
 * Note that for infeasibility, we just need to proof that the problem 
 \f[ \begin{array}{ll} \min & 0\\ s.t. & Ax = b\\ & l\leq x\leq b\\ 
 \end{array} \f]
 * is infeasible, but it's dual is
 \f[ \begin{array}{ll} \max & by - ud_u + ld_l\\ s.t. & A^ty +Id_l - Id_u = 0\\ 
 & d_u,d_l\geq0\\ \end{array} \f]
 * wich is always feasible (provided \f$y\geq0\f$ (set \f$ (y,d_u,d_l)=0\f$), 
 * and thus we just need to check whether the objective value is \f$\neq 0\f$ 
 * and we have a proof of infeasibility for the primal. That's what this 
 * function perform as a test.
 * */
int QSexact_infeasible_test (mpq_QSdata * p,
														 mpq_t * d_sol);

/* ========================================================================= */
/** @brief create a copy of a mpq_t array into a double array.
 * @param array mpq_t array from where we will create the values. */
#define QScopy_array_mpq_dbl(array) ({ \
	mpq_t*__larray = (array);\
	register unsigned __lsz = __EGlpNumArraySize(__larray);\
	double*__lres = dbl_EGlpNumAllocArray(__lsz);\
	while(__lsz--)\
	{\
		if(mpq_equal(__larray[__lsz],mpq_ILL_MAXDOUBLE))\
			__lres[__lsz] = dbl_ILL_MAXDOUBLE;\
		else if(mpq_equal(__larray[__lsz],mpq_ILL_MINDOUBLE))\
			__lres[__lsz] = dbl_ILL_MINDOUBLE;\
		else __lres[__lsz] = mpq_get_d(__larray[__lsz]);\
	}\
	__lres;})

/* ========================================================================= */
/** @brief create a copy of a mpq_t array into a mpf_t array.
 * @param array mpq_t array from where we will create the values. */
#define QScopy_array_mpq_mpf(array) ({ \
	mpq_t*__larray = (array);\
	register unsigned __lsz = __EGlpNumArraySize(__larray);\
	mpf_t*__lres = mpf_EGlpNumAllocArray(__lsz);\
	while(__lsz--)\
	{\
		if(mpq_equal(__larray[__lsz],mpq_ILL_MAXDOUBLE))\
			mpf_set(__lres[__lsz], mpf_ILL_MAXDOUBLE);\
		else if(mpq_equal(__larray[__lsz],mpq_ILL_MINDOUBLE))\
			mpf_set(__lres[__lsz], mpf_ILL_MINDOUBLE);\
		else mpf_set_q(__lres[__lsz],__larray[__lsz]);\
	}\
	__lres;})

/* ========================================================================= */
/** @brief create a copy of a double array into mpq_t array.
 * @param array original array of double values (note that this array must have
 * been allocated with dbl_EGlpNumAllocArray for this function to work). */
#define QScopy_array_dbl_mpq(array) ({ \
	double*__larray = (array);\
	register unsigned __lsz = __EGlpNumArraySize(__larray);\
	mpq_t*__lres = mpq_EGlpNumAllocArray(__lsz);\
	while(__lsz--)\
	{\
		if(__larray[__lsz] == dbl_ILL_MAXDOUBLE)\
			mpq_set(__lres[__lsz],mpq_ILL_MAXDOUBLE);\
		else if(__larray[__lsz] == dbl_ILL_MINDOUBLE)\
			mpq_set(__lres[__lsz],mpq_ILL_MINDOUBLE);\
		else mpq_EGlpNumSet(__lres[__lsz],__larray[__lsz]);\
	}\
	__lres;})

/* ========================================================================= */
/** @brief create a copy of a mpf_t array into mpq_t array.
 * @param array original array of double values (note that this array must have
 * been allocated with __EGlpNumAllocArray for this function to work). */
#define QScopy_array_mpf_mpq(array) ({ \
	mpf_t*__larray = (array);\
	register unsigned __lsz = __EGlpNumArraySize(__larray);\
	mpq_t*__lres = mpq_EGlpNumAllocArray(__lsz);\
	while(__lsz--)\
	{\
		if(mpf_cmp(__larray[__lsz],mpf_ILL_MAXDOUBLE)==0)\
			mpq_set(__lres[__lsz],mpq_ILL_MAXDOUBLE);\
		else if(mpf_cmp(__larray[__lsz],mpf_ILL_MINDOUBLE)==0)\
			mpq_set(__lres[__lsz],mpq_ILL_MINDOUBLE);\
		mpq_set_f(__lres[__lsz],__larray[__lsz]);\
	}\
	__lres;})

/* ========================================================================= */
/** @brief Write a given row from the LP into the given stream, in exact
 * arithmetic */
void QSexact_write_row (EGioFile_t * out_f,
												mpq_ILLlpdata * lp,
												int row);

/* ========================================================================= */
/** @brief Set the number of bits to use with mpf_t type numbers and change all
 * internal constants as needed. */
#define QSexact_set_precision(precision) mpf_QSset_precision(precision)

#ifndef QS_EXACT_MAX_ITER
/* ========================================================================= */
/** @brief This constant define the maximum number of try's for the exact solver
 * with mpf_t numbers while incrementing the precision */
#define QS_EXACT_MAX_ITER 12
#endif

/* ========================================================================= */
/** @brief test whether given basis is primal and dual feasible in rational arithmetic. 
 * @param p_mpq   the problem data.
 * @param basis   basis to be tested.
 * @param result  where to store whether given basis is primal and dual feasible.
 * @param msg_lvl message level.
 */
int QSexact_basis_optimalstatus(
   mpq_QSdata * p_mpq,
   QSbasis* basis,
   char* result,
   const int msg_lvl
   );

/* ========================================================================= */
/** @brief test whether given basis is dual feasible in rational arithmetic. 
 * @param p_mpq   the problem data.
 * @param basis   basis to be tested.
 * @param result  where to store whether given basis is dual feasible.
 * @param dobjval where to store dual solution value in case of dual feasibility (if not NULL).
 * @param msg_lvl message level.
 */
int QSexact_basis_dualstatus(
   mpq_QSdata * p_mpq,
   QSbasis* basis,
   char* result,
   mpq_t* dobjval,
   const int msg_lvl
   );

/* ========================================================================= */
/** @brief test whether given basis is dual feasible in rational arithmetic. 
 * if wanted it will first directly test the corresponding approximate dual and primal solution 
 * (corrected via dual variables for bounds and primal variables for slacks if possible) for optimality
 * before performing the dual feasibility test on the more expensive exact basic solution. 
 * @param p_mpq   the problem data.
 * @param basis   basis to be tested.
 * @param useprestep whether to directly test approximate primal and dual solution first.
 * @param dbl_p_sol  approximate primal solution to use in prestep 
 *                   (NULL in order to compute it by dual simplex in double precision with given starting basis).
 * @param dbl_d_sol  approximate dual solution to use in prestep 
 *                   (NULL in order to compute it by dual simplex in double precision with given starting basis).
 * @param result  where to store whether given basis is dual feasible.
 * @param dobjval where to store dual solution value in case of dual feasibility (if not NULL).
 * @param msg_lvl message level.
 */
int QSexact_verify (
   mpq_QSdata * p_mpq,
   QSbasis* basis,
   int useprestep,
   double* dbl_p_sol,
   double* dbl_d_sol,
   char* result,
   mpq_t* dobjval,
   const int msg_lvl
   );

/* ========================================================================= */
/** @brief Given an mpq_QSdata problem, solve it exactly.
 * @param x if not null, we store here the primal solution to the 
 * problem (if it exist).
 * @param y if not null, we store here the dual solution to the
 * problem, 
 * @param p_mpq problem to solve exactly.
 * @param status pointer to the integer where we will return the status
 * of the problem, either optimal, infeasible, or unbounded (we could also 
 * return time out).
 * @param simplexalgo whether to use primal or dual simplex while solving
 * to optimality the problem.
 * @param basis if not null, use the given basis to start the
 * iteration of simplex, and store here the optimal basis (if found).
 * @return zero on success, non-zero otherwise. */
int QSexact_solver (mpq_QSdata * p_mpq,
										mpq_t * const x,
										mpq_t * const y,
										QSbasis * const basis,
										int simplexalgo,
										int *status);

/* ========================================================================= */
/** @brief Initializator for global data, this is needed mainly for defining
 * constants in extended floating point precision and for rational precision.
 * This call should be done BEFORE any mpq_xxx mpf_xxx QSxx EGxx call */
extern void QSexactStart(void);
/* ========================================================================= */
/** @brief This function must be called at the end of the program to free all
 * internal data used in the QSexact structures, once this function is called
 * any operation on EGxxx mpq_xxx mpf_xx QSxx may fail.
 * */
extern void QSexactClear(void);
/* ========================================================================= */
/** @brief indicate if the global data needed for QSexact has been initialized,
 * if zero, initialization routine should be called. This is provided to allow
 * syncronization between libraries */
extern int __QSexact_setup;
/** @} */
/* ========================================================================= */
/* end of exact.h */
#endif

