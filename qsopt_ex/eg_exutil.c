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
/** @file 
 * @ingroup Esolver */
/** @addtogroup Esolver */
/** @{ */
#include "eg_exutil.h"

#include "logging-private.h"

/* ========================================================================= */
/** @name EXutilStatics
 * Variables asociated with the #EXutilApproximate function, we use them
 * as static to save some time in intialization */
/*@{*/
/** @brief Array of integers used in the continued fraction method */
static mpz_t Z[7];
/**@brief rational remainder used in the continued fraction method */
static mpq_t cvl;
/* ========================================================================= */
/** @brief Initialize the static variables at start-up */
void EXutilDoInit (void)
{
	unsigned __EXui;
	EGlpNumStart();
	mpq_init (cvl);
	for (__EXui = 7; __EXui--;)
		mpz_init (Z[__EXui]);
}

/* ========================================================================= */
/** @brief Clear all memory related to the static variables */
void EXutilDoClear (void)
{
	unsigned __EXui;
	mpq_clear (cvl);
	for (__EXui = 7; __EXui--;)
		mpz_clear (Z[__EXui]);
	EGlpNumClear();

}

/*@}*/

/* ========================================================================= */
int EXutilIntegralize (const unsigned n,
											 mpq_t * const  a,
											 mpq_t b,
											 mpq_t maxabs)
{
	mpz_t lcm,
	  gcd;
	register unsigned int i;
	mpz_init (lcm);
	mpz_init (gcd);
	mpz_set (lcm, mpq_denref (b));
	mpz_set (gcd, mpq_numref (b));
	if (mpz_cmp_ui (gcd, 0UL) == 0)
		mpz_set_ui (gcd, 1UL);
	/* compute the greatest common divisor ammong the numerator of a_i and b,
	 * and the least common multiple ammong the denominators of a_i and b */
	for (i = n; i--;)
	{
		mpz_lcm (lcm, lcm, mpq_denref (a[i]));
		mpz_gcd (gcd, gcd, mpq_numref (a[i]));
	}
	/* divide everything by lcm/gcd */
	mpz_mul (mpq_numref (b), mpq_numref (b), lcm);
	mpz_mul (mpq_denref (b), mpq_denref (b), gcd);
	mpq_canonicalize (b);
	mpz_abs (mpq_numref (maxabs), mpq_numref (b));
	for (i = n; i--;)
	{
		mpz_mul (mpq_denref (a[i]), mpq_denref (a[i]), gcd);
		mpz_mul (mpq_numref (a[i]), mpq_numref (a[i]), lcm);
		mpq_canonicalize (a[i]);
		if (mpz_cmpabs (mpq_numref (maxabs), mpq_numref (a[i])) < 0)
			mpz_abs (mpq_numref (maxabs), mpq_numref (a[i]));
	}
	/* ending */
	mpz_set_ui (mpq_denref (maxabs), 1UL);
	mpz_clear (gcd);
	mpz_clear (lcm);
	return 0;
}

/* ========================================================================= */
int EXutilSimplify (const unsigned n,
										mpq_t * const  a,
										mpq_t b)
{
	register unsigned i;
	mpq_t maxabs;
	mpq_init (maxabs);
	EXutilIntegralize (n, a, b, maxabs);
	/* normalize the cut so that |(a,b)|_inf == 1 */
	if (mpz_cmp_ui (mpq_numref (maxabs), 0UL))
	{
		mpq_div (b, b, maxabs);
		for (i = n; i--;)
			mpq_div (a[i], a[i], maxabs);
	}
	/* ending */
	mpq_clear (maxabs);
	return 0;
}

/* ========================================================================= */
void mpq_GomoryCoeff (mpq_t rop,
											mpq_t coef,
											unsigned const is_int,
											int const bound,
											unsigned const cut_mlt,
											mpq_t b_frac)
{
	mpq_t fj;
	mpq_init (fj);
	mpq_set_ui (rop, 0UL,1UL);
	/* if the variable is integer */
	if (is_int)
	{
		if (mpq_IsInteger (coef))
		{
			mpq_set (rop, coef);
			mpq_EGlpNumMultUiTo (rop, cut_mlt);
			mpq_EGlpNumMultTo (rop, b_frac);
		}
		/* if the variable is integer, but the coefficient id fractional */
		else
		{
			mpq_set (fj, coef);
			mpq_EGlpNumMultUiTo (fj, cut_mlt);
			mpq_FracPart (fj, fj);
			/* if the variable is complemented to its lower bound */
			if (bound == 'L')
			{
				if (mpq_cmp (fj, b_frac) <= 0)
				{
					mpq_set (rop, coef);
					mpq_EGlpNumMultUiTo (rop, cut_mlt);
					mpq_EGlpNumFloor (rop, rop);
					mpq_EGlpNumMultTo (rop, b_frac);
					mpq_EGlpNumAddTo (rop, fj);
				}
				else
				{
					mpq_set (rop, coef);
					mpq_EGlpNumMultUiTo (rop, cut_mlt);
					mpq_EGlpNumCeil (rop, rop);
					mpq_EGlpNumMultTo (rop, b_frac);
				}
			}
			/* if the variable is complemented to its upper bound */
			else
			{
				mpq_EGlpNumSubTo (fj, mpq_oneLpNum);
				mpq_neg (fj, fj);
				if (mpq_cmp (fj, b_frac) <= 0)
				{
					mpq_set (rop, coef);
					mpq_EGlpNumMultUiTo (rop, cut_mlt);
					mpq_EGlpNumCeil (rop, rop);
					mpq_EGlpNumMultTo (rop, b_frac);
					mpq_EGlpNumSubTo (rop, fj);
				}
				else
				{
					mpq_set (rop, coef);
					mpq_EGlpNumMultUiTo (rop, cut_mlt);
					mpq_EGlpNumFloor (rop, rop);
					mpq_EGlpNumMultTo (rop, b_frac);
				}
			}
		}
	}
	/* if the variable is continuous */
	else
	{
		/* if the variable is complemented to its lower bound */
		if (bound == 'L')
		{
			if (mpq_cmp_ui (coef, 0UL,1UL) > 0)
			{
				mpq_set (rop, coef);
				mpq_EGlpNumMultUiTo (rop, cut_mlt);
			}
		}
		/* if the variable is complemented to its upper bound */
		else
		{
			if (mpq_cmp_ui (coef, 0UL,1UL) < 0)
			{
				mpq_set (rop, coef);
				mpq_EGlpNumMultUiTo (rop, cut_mlt);
			}
		}
	}
	/* done */
	mpq_clear (fj);
}

/* ========================================================================= */
int EXutilExpandLogicals (mpq_QSdata * const act_lp,
													mpq_t *  const vector,
													mpq_t b,
													mpq_ILLlp_rows * const lprows)
{
	const int n_rows = act_lp->qslp->nrows;
	const int n_struct = act_lp->qslp->nstruct;
	int const *const rowmap = act_lp->qslp->rowmap;
	mpq_t * const rhs = act_lp->qslp->rhs;
	mpq_ILLmatrix *const A = &(act_lp->qslp->A);
	int rowbeg;
	int rowcnt;
	int *rowind;
	mpq_t * rowval;
	register int i,
	  k;
	for (i = n_rows; i--;)
	{
		/* convert the vector */
		if (mpz_cmp_ui (mpq_numref (vector[i + n_struct]), 0UL))
		{
			rowbeg = lprows->rowbeg[i];
			rowcnt = lprows->rowcnt[i];
			rowind = lprows->rowind + rowbeg;
			rowval = lprows->rowval + rowbeg;
			/* we use slack again as the multiplier that we need to add the
			 * row to the current vector so as to make dissapear the slack */
			mpq_neg (cvl, vector[i + n_struct]);
			mpq_div (cvl, cvl, A->matval[A->matbeg[rowmap[i]]]);
			MESSAGE (EX_UTIL_VERBOSE + 100, "Replacing constraint %s with multiple"
							 " %lf from integer part", act_lp->qslp->rownames[i],
							 mpq_get_d (cvl));
			mpq_EGlpNumAddInnProdTo (b, rhs[i], cvl);
			for (k = rowcnt; k--;)
				mpq_EGlpNumAddInnProdTo (vector[rowind[k]], cvl, rowval[k]);
			mpq_set_ui (vector[i + n_struct], 0UL, 1UL);
		}
	}
	return 0;
}

/* ========================================================================= */
void EXutilApproximate (mpq_t var,
												mpq_t ori,
												unsigned const max_den)
{
	/* local variables */
	unsigned lsng = mpz_cmp_ui (mpq_numref (ori), 0UL) < 0 ? 1U : 0U;
	int i;
	mpq_t __lpnum__;
	mpq_init(__lpnum__);
	/* check if the given number is zero, if so, set to zero var and return */
	if (mpz_cmp_ui (mpq_numref (ori), 0UL) == 0)
	{
		return;
	}
	/* if not, then we have some work to do */
	/* now we initialize the internal numbers */
	mpq_abs (cvl, ori);
	for (i = 7; i--;)
		mpz_set_ui (Z[i], 0UL);
	mpz_set_ui (Z[0], 1UL);
	mpz_set_ui (Z[4], 1UL);
	mpz_fdiv_q (Z[1], mpq_numref (cvl), mpq_denref (cvl));
	mpq_set_z (__lpnum__, Z[1]);
	mpq_sub (cvl, cvl, __lpnum__);
	/* now we loop until the next t's is more than mpf_eps */
	/* the formula is 
	 * p_i = t_i*p_{i-1} + p_{i-2}, and 
	 * q_i = t_i*q_{i-1} + q_{i-2} 
	 * note that |x-p_i/q_i|<1/q_i^2
	 * for us t_i = Z[6], and the current number is either [0,1,2] in the Z
	 * array, we use those popsitions ciclicly, and use the four position as a
	 * temporary number, Z+4 is used to store q's, at the beginning i = 1. */
	while (1)
	{
		if (mpq_cmp_ui (cvl, 1UL, (unsigned long)max_den) < 0 || (mpz_cmp_ui (Z[4], (unsigned long)max_den) > 0))
		{
			mpz_set (mpq_denref (var), Z[4]);
			mpz_set (mpq_numref (var), Z[1]);
			break;
		}
		/* first run */
		mpq_inv (cvl, cvl);
		mpz_fdiv_q (Z[6], mpq_numref (cvl), mpq_denref (cvl));
		mpq_set_z (__lpnum__, Z[6]);
		mpq_sub (cvl, cvl, __lpnum__);
		mpz_set (Z[2], Z[0]);
		mpz_addmul (Z[2], Z[1], Z[6]);
		mpz_set (Z[5], Z[3]);
		mpz_addmul (Z[5], Z[4], Z[6]);
		if (mpq_cmp_ui (cvl, 1UL, (unsigned long)max_den) < 0 || (mpz_cmp_ui (Z[5], (unsigned long)max_den) > 0))
		{
			mpz_set (mpq_denref (var), Z[5]);
			mpz_set (mpq_numref (var), Z[2]);
			break;
		}
		/* second run */
		mpq_inv (cvl, cvl);
		mpz_fdiv_q (Z[6], mpq_numref (cvl), mpq_denref (cvl));
		mpq_set_z (__lpnum__, Z[6]);
		mpq_sub (cvl, cvl, __lpnum__);
		mpz_set (Z[0], Z[1]);
		mpz_addmul (Z[0], Z[2], Z[6]);
		mpz_set (Z[3], Z[4]);
		mpz_addmul (Z[3], Z[5], Z[6]);
		if (mpq_cmp_ui (cvl, 1UL, (unsigned long)max_den) < 0 || (mpz_cmp_ui (Z[3], (unsigned long)max_den) > 0))
		{
			mpz_set (mpq_denref (var), Z[3]);
			mpz_set (mpq_numref (var), Z[0]);
			break;
		}
		/* third run */
		mpq_inv (cvl, cvl);
		mpz_fdiv_q (Z[6], mpq_numref (cvl), mpq_denref (cvl));
		mpq_set_z (__lpnum__, Z[6]);
		mpq_sub (cvl, cvl, __lpnum__);
		mpz_set (Z[1], Z[2]);
		mpz_addmul (Z[1], Z[0], Z[6]);
		mpz_set (Z[4], Z[5]);
		mpz_addmul (Z[4], Z[3], Z[6]);
	}
	/* ending */
	mpq_canonicalize (var);
	if (lsng)
		mpq_neg (var, var);
	/* clean-up */
	mpq_clear(__lpnum__);
	return;
}

/* ========================================================================= */
void EXutilOverEstimate (mpq_t var,
												 mpq_t ori,
												 unsigned const max_den)
{
	EXutilApproximate (var, ori, max_den);
	/* check if var is < ori, if so, we must add one to the numerator */
	if (mpq_cmp (ori, var) > 0)
	{
		mpq_set_ui (cvl, 1UL, (unsigned long)(max_den * max_den));
		mpq_add (var, var, cvl);
		EXIT (mpq_cmp (ori, var) > 0, "Imposible!");
	}
	return;
}

/* ========================================================================= */
void EXutilNicefy (mpq_QSdata * const act_prob,
									 const unsigned char *const var_stat,
									 const unsigned max_den,
									 mpq_t *  a,
									 mpq_t b,
									 int const sense)
{
	const unsigned square = max_den * max_den;
	const int nstruct = act_prob->qslp->nstruct;
	const int *const  structmap = act_prob->qslp->structmap;
	mpq_t *const  lower = act_prob->qslp->lower;
	mpq_t *const  upper = act_prob->qslp->upper;
	mpq_t num1,
	  num2;
	register int i;
	int colid = 0;
	int sign = 0;
	unsigned cur_stat = 0;
	mpq_init (num1);
	mpq_init (num2);
	if (sense != 'L' && sense != 'G')
		return;
	/* we internally assume that the inequality is of the form ax >= b */
	if (sense == 'L')
	{
		mpq_neg (b, b);
		for (i = nstruct; i--;)
			mpq_neg (a[i], a[i]);
	}
	/* now we first approximate each coefficient */
	for (i = nstruct; i--;)
	{
		colid = structmap[i];
		cur_stat = var_stat[colid];
		/* if the variables is not bounded, we can't nicefy the coefficient */
		if ((cur_stat & (EX_STATUS_UB | EX_STATUS_LB)) == 0)
			continue;
		/* if the numerator is less than square, there is nothing to do */
		if (mpz_cmp_ui (mpq_denref (a[i]), (unsigned long)square) <= 0)
			continue;
		EXutilApproximate (num1, a[i], max_den);
		mpq_sub (num2, num1, a[i]);
		sign = mpz_cmp_ui (mpq_numref (num2), 0UL);
		/* depending on the side that the approximation land we see what we have to
		 * do, first case, num1 > a[i] */
	REDO:
		if (sign > 0)
		{
			/* if the variable is bounded bellow, we just update the RHS and set the
			 * coefficient */
			if (cur_stat & EX_STATUS_LB)
			{
				EXIT ((mpq_EGlpNumIsEqqual (lower[colid], mpq_ILL_MINDOUBLE)),
							"Imposible");
				mpq_EGlpNumAddInnProdTo (b, num2, lower[colid]);
				mpq_set (a[i], num1);
			}
			/* otherwise, we can't approximate by above, but we have to approximate
			 * by bellow */
			else
			{
				mpq_set_ui (num2, 1UL, (unsigned long)square);
				mpq_sub (num1, num1, num2);
				mpq_sub (num2, num1, a[i]);
				sign = mpz_cmp_ui (mpq_numref (num2), 0UL);
				goto REDO;
			}
		}
		/* otherwise, we have that num1 < a[i] */
		else if (sign < 0)
		{
			/* if the variable is bounded by above, we just update the RHS and set
			 * the coefficient */
			if (cur_stat & EX_STATUS_UB)
			{
				EXIT ((mpq_EGlpNumIsEqqual (upper[colid], mpq_ILL_MAXDOUBLE)),
							"Imposible");
				mpq_EGlpNumAddInnProdTo (b, num2, upper[colid]);
				mpq_set (a[i], num1);
			}
			/* otherwise, we can't approximate by bellow, but we have to approximate
			 * by above. */
			else
			{
				mpq_set_ui (num2, 1UL, (unsigned long)square);
				mpq_add (num1, num1, num2);
				mpq_sub (num2, num1, a[i]);
				sign = mpz_cmp_ui (mpq_numref (num2), 0UL);
				goto REDO;
			}
		}
	}
	/* now we round the RHS, we do this by adding the constraint 0 >= -1 */
	EXutilApproximate (num1, b, max_den);
	if (mpq_cmp (num1, b) > 0)
	{
		mpq_set_ui (num2, 1UL, (unsigned long)square);
		mpq_sub (num1, num1, num2);
	}
	mpq_set (b, num1);
	/* before ending, we return the constraint to its normal form */
	if (sense == 'L')
	{
		mpq_neg (b, b);
		for (i = nstruct; i--;)
			mpq_neg (a[i], a[i]);
	}
	/* ending */
	mpq_clear (num2);
	mpq_clear (num1);
	return;
}

/* ========================================================================= */
/** @} */
/* end eg_exutil.c */
