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
#ifndef __EG_EX_UTIL__
#define __EG_EX_UTIL__

#include <stdlib.h>

#include "eg_lpnum.h"
#include "eg_io.h"

#include "qstruct_mpq.h"
#include "lpdata_mpq.h"

/** @file 
 * @ingroup Esolver */
/** @addtogroup Esolver */
/** @{ */
/* ========================================================================= */
/** @brief status for variables that are logicals */
#define EX_STATUS_LGC 1U

/* ========================================================================= */
/** @brief status for variables that are structural */
#define EX_STATUS_STR 2U

/* ========================================================================= */
/** @brief status for variables that are integer */
#define EX_STATUS_INT 4U

/* ========================================================================= */
/** @brief status for variables that bounded from below */
#define EX_STATUS_LB 8U

/* ========================================================================= */
/** @brief status for variables that bounded from above */
#define EX_STATUS_UB 16U

/* ========================================================================= */
/** @brief status for integer variable fixed at its upper bound. This status is
 * not defined while calling cut callbacks. */
#define EX_STATUS_FIX_UB 32U

/* ========================================================================= */
/** @brief status for integer variable fixed at its lower bound. This status is
 * not defined while calling cut callbacks. */
#define EX_STATUS_FIX_LB 64U

/* ========================================================================= */
/** @brief status for integer variable among our set of selected integer
 * variables. This status is not defined while calling cut callbacks. */
#define EX_STATUS_BESTFRAC 128U

/* ========================================================================= */
/** @brief given a cut ax <=> b, write it in integer form ,i.e. set all a,b to
 * integer in such a way that a_i and b are all relativelly prime. 
 * @param n size of the a vector.
 * @param a RHS of the inequality
 * @param b LHS of the inequality.
 * @param maxabs return the maximum absolute value among all resulting a_i,b.
 * @return zero on success, non-zero otherwise.
 * */
int EXutilIntegralize (const unsigned n,
											 mpq_t * const  a,
											 mpq_t b,
											 mpq_t maxabs);

/* ========================================================================= */
/** @brief verbosity level */
#define EX_UTIL_VERBOSE 100

/* ========================================================================= */
/** @brief given two vectors (a,b) and (v,w), compute its inner product and
 * store it into rop.
 * @param dim dimmension of the a and v part of the vector.
 * @param a first part of the first vector.
 * @param b second part of the first vector.
 * @param v first part of the second vector.
 * @param w second part of the second vector.
 * @param rop where we store the result.
 * @note we may take a == v and/or b == w.
 * */
#define EXutilInnProd(dim,a,b,v,w,rop) do{\
	register unsigned __EXuti = (dim);\
	mpq_mul(rop,b,w);\
	while(__EXuti--) mpq_EGlpNumAddInnProdTo(rop,(a)[__EXuti],(v)[__EXuti]);\
	} while(0);

/* ========================================================================= */
/** @brief Compute the number of non-zeros in a given vector. 
 * @param dim size of the a vector.
 * @param a vector where we are operating.
 * @param rop where to return the number of non-zeros (it should be an integer
 * variable, not a pointer to such a variable) */
#define EXutilNzSz(dim,a,rop) do{\
	register unsigned __EXuti = (dim);\
	rop = 0;\
	while(__EXuti--) if(mpz_cmp_ui(mpq_numref((a)[__EXuti]),0UL)) rop++;}while(0)

/* ========================================================================= */
/** @brief Compute the L_1 norm of a given vector. 
 * @param dim size of the a vector.
 * @param a vector where we are operating.
 * @param rop where to retirn the L1 norm value */
#define EXutilL1Norm(dim,a,rop) do{\
	mpq_t*const  __EXuta = (a);\
	mpq_t __qtmp__;\
	register unsigned __EXuti = (dim);\
	mpq_init(__qtmp__);\
	mpq_set_ui(rop,0UL,1UL);\
	while(__EXuti--){\
		mpq_abs(__qtmp__,__EXuta[__EXuti]);\
		mpq_add(rop,rop,__qtmp__);}mpq_clear(__qtmp__);}while(0)

/* ========================================================================= */
/** @brief given a cut ax <=> b, write it in normalized form ,i.e. set all a,b 
 * to integer in such a way that a_i and b are all relativelly prime, and
 * divide them all over the
 * maximum such (a_i,b) (so that the infinity norm of (a,b) is one. 
 * @param n size of the a vector.
 * @param a RHS of the inequality
 * @param b LHS of the inequality.
 * @return zero on success, non-zero otherwise.
 * */
int EXutilSimplify (const unsigned n,
										mpq_t * const  a,
										mpq_t b);

/* ========================================================================= */
/** @brief asign to the first number the fractional part of the second, i.e.
 * 	\f$ rop = op1 - \lfloor op1 \rfloor \f$.
 * 	@param rop where to retunr our result.
 * 	@param op1 number for wich we want to compute its fractional part */
#define mpq_FracPart(rop, op1) do{\
	mpz_fdiv_r (mpq_numref (rop), mpq_numref (op1), mpq_denref (op1));\
	mpz_set (mpq_denref (rop), mpq_denref (op1));\
	mpq_canonicalize((rop));\
} while(0)

/* ========================================================================= */
/** @brief test if the given number is integer.
 * @param op number to test.
 * @return one if the nummber is integer, zero- otherwise. */
#define mpq_IsInteger(op) ({\
	(mpz_cmp(mpq_denref(op),mpz_oneLpNum)==0);})

/* ========================================================================= */
/** @brief round to \f$-\infty\f$ the given number to the closest fraction 
 * of the form \f$a/2^{exp}\f$ from bellow.
 * @param op number to round.
 * @param exp exponent to use in the fraction */
#define mpq_FroundExp(op,exp) do{\
			mpz_t __ztmp__;\
			mpz_init(__ztmp__);\
			mpz_mul_2exp(__ztmp__,mpq_numref((op)),(exp));\
			mpz_fdiv_q(mpq_numref((op)),__ztmp__,mpq_denref((op)));\
			mpz_mul_2exp(mpq_denref((op)),mpz_oneLpNum,(exp));\
			mpq_canonicalize((op));mpz_clear(__ztmp__);}while(0)

/* ========================================================================= */
/** @brief round to \f$+\infty\f$ the given number to the closest fraction 
 * of the form \f$a/2^{exp}\f$ from above.
 * @param op number to round.
 * @param exp exponent to use in the fraction */
#define mpq_CroundExp(op,exp) do{\
			mpz_t __ztmp__;\
			mpz_init(__ztmp__);\
			mpz_mul_2exp(__ztmp__,mpq_numref((op)),(exp));\
			mpz_cdiv_q(mpq_numref((op)),__ztmp__,mpq_denref((op)));\
			mpz_mul_2exp(mpq_denref((op)),mpz_oneLpNum,(exp));\
			mpq_canonicalize((op));mpz_clear(__ztmp__);}while(0)

/* ========================================================================= */
/** @brief round to \f$0\f$ the given number to the closest fraction 
 * of the form \f$a/2^{exp}\f$ towards zero.
 * @param op number to round.
 * @param exp exponent to use in the fraction */
#define mpq_TroundExp(op,exp) do{\
			mpz_t __ztmp__;\
			mpz_init(__ztmp__);\
			mpz_mul_2exp(__ztmp__,mpq_numref((op)),(exp));\
			mpz_tdiv_q(mpq_numref((op)),__ztmp__,mpq_denref((op)));\
			mpz_mul_2exp(mpq_denref((op)),mpz_oneLpNum,(exp));\
			mpq_canonicalize((op));mpz_clear(__ztmp__);}while(0)

/* ========================================================================= */
/** @brief compute the gomory coefficient of the variable given the original
 * coefficient, the multiplier, and all relevant information.
 * @param rop Where we return the gomory coefficient. (it should be different
 * location than the original coefficient).
 * @param coef original coefficient in the equality that we are using to
 * derive the gomory cut.
 * @param is_int this is either zero (indicating that the variable asociated
 * with this coefficient is continuous) or one (indicating that the variable
 * asociated with this coeffcient is integer).
 * @param bound indicate if this variable was complemented to its lower bound
 * (then L) or to its upper bound (then U), any other value will generate an
 * error.
 * @param cut_mlt multiplier to use to the coefficient (and thus the effective
 * coefficient would be coef*cut_mlt).
 * @param b_frac fractional part of the RHS in the equation (after
 * complementing variables). */
void mpq_GomoryCoeff (mpq_t rop,
											mpq_t coef,
											unsigned const is_int,
											int const bound,
											unsigned const cut_mlt,
											mpq_t b_frac);

/* ========================================================================= */
/** @brief Given a vector in QSopt external form, and a row description of the
 * related LP, re-write the vector using only real variables, we do that by
 * substracting the equation defining the logical variable multiplied by the
 * coefficient of the logical variable in the vector to the vector.
 * @param act_lp lp where we are working.
 * @param vector vector of length at least nrows + nstruct where we want to
 * replace all logical coefficients.
 * @param lprows row description of the given LP.
 * @param rhs if we look at vector,rhs as an inequality, then we eliminate the
 * slack coefficient form the inequality as a whole.
 * @return zero on success, non-zero otherwise.
 * */
int EXutilExpandLogicals (mpq_QSdata * const act_lp,
													mpq_t *  const vector,
													mpq_t rhs,
													mpq_ILLlp_rows * const lprows);

/* ========================================================================= */
/** @brief Approximate <A HREF=http://mathworld.wolfram.com/ContinuedFraction.html TARGET=_top>using continued fractions method</A> a given rational 
 * \f$\frac{a}{b} \f$ with another rational \f$\frac{a'}{b'}\f$ that satisfy 
 * that \f$ b' < max_den^2 \f$ and also 
 * \f$|\frac{a}{b} - \frac{a'}{b'}|\leq\frac1{max_den^2}\f$.
 * @param ori original coefficient that we want to represent as a/b with b <=
 * max_den^2. 
 * @param dest we return here the resulting number. 
 * @param max_den maximum allowed denominator in the new representation. */
void EXutilApproximate (mpq_t dest,
												mpq_t ori,
												unsigned const max_den);

/* ========================================================================= */
/** @brief Overestimate the given coefficient by another rational that is
 * representble with denominators not bigger than max_den^2.
 * @param ori original coefficient that we want to represent as a/b with b <=
 * max_den^2. 
 * @param dest we return here the resulting number, note that we always
 * insure that the returned value is bigger than the original value. 
 * @param max_den maximum allowed denominator in the new representation. */
void EXutilOverEstimate (mpq_t dest,
												 mpq_t ori,
												 unsigned const max_den);

/* ========================================================================= */
/** @brief Given an inequality, we try to re-write so that no denominator is
 * bigger than the square of the given number, and ensuring validity. for 
 * coefficients that can't be `nacified' we leave them intact. the process 
 * imply adding multiples of the bounds on variables, and at the end, nicify 
 * the rhs of the inequality.
 * @param max_den the square of this value is the maximum denominator allowed.
 * @param a hand side of the inequality.
 * @param sense sense of the inequality, it should be either 'L' or 'G'.
 * @param b right hand side of the inequality.
 * @param act_prob LP from where we draw the bounds on the variables.
 * @param var_stat status (as defined in #EXmipinfo_t::var_stat) for all
 * variables in the LP, in the internal QSopt ordering. 
 * @note The length of the a vector is at least nstruct, and we assume that
 * entry a[k] corresnpond to the coefficient associated with the k-th
 * structural variable inside. */
void EXutilNicefy (mpq_QSdata * const act_prob,
									 const unsigned char *const var_stat,
									 const unsigned max_den,
									 mpq_t *  a,
									 mpq_t b,
									 int const sense);

/* ========================================================================= */
/** @brief given a variable in internal number, return a pointer to its name.
 * @param iid internal ordering number.
 * @param QSlp pointer to the mpq_QSdata structure containing the LP.
 * @param QSinv_map pointer to an array containing the inverse map from internal
 * numbering to external numbering as in #EXmipinfo_t::inv_map.
 * @return pointer to its name.
 * @note If the variable is a slack variable, it return the name of the
 * inequality. */ 
#define EXutilIidToStr(iid,QSlp,QSinv_map) ({\
	mpq_QSdata*const __EXlp = (QSlp);\
	const int __EXeid = (QSinv_map)[(iid)];\
	(__EXeid >= __EXlp->qslp->nstruct) ? __EXlp->qslp->rownames[__EXeid - __EXlp->qslp->nstruct] : __EXlp->qslp->colnames[__EXeid];})

/* ========================================================================= */
/** @brief Initialize the static variables at start-up */
extern void EXutilDoInit (void);
/* ========================================================================= */
/** @brief Clear all memory related to the static variables */
extern void EXutilDoClear (void);
/* ========================================================================= */
/** @} */
/* end eg_exutil.h */
#endif
