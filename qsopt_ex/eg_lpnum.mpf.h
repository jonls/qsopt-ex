/* EGlib "Efficient General Library" provides some basic structures and
 * algorithms commons in many optimization algorithms.
 *
 * Copyright (C) 2005 Daniel Espinoza and Marcos Goycoolea.
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
#ifndef __EG_LPNUM_MPF__
#define __EG_LPNUM_MPF__

#include <stdlib.h>
#include <string.h>

#include <gmp.h>

#include "eg_lpnum.h"

/** @file
 * @ingroup EGlpNum */
/** @addtogroup EGlpNum */
/** @{ */
/* ========================================================================= */
/** @brief This is the smallest difference (with the current precision) that can
 * be distinguished betwen 1.0 and it's clossest representable number, in some
 * sense it is the absolute minimum epsilon for comparisons */
extern mpf_t mpf_eps;

/* ========================================================================= */
/** extern definitions of constaants for different set-ups */
extern const mpf_t __zeroLpNum_mpf__;
extern const mpf_t __oneLpNum_mpf__;
extern const mpf_t __MaxLpNum_mpf__;
extern const mpf_t __MinLpNum_mpf__;
#define mpf_zeroLpNum __zeroLpNum_mpf__
#define mpf_oneLpNum  __oneLpNum_mpf__
#define mpf_epsLpNum  mpf_eps
#define mpf_MaxLpNum  __MaxLpNum_mpf__
#define mpf_MinLpNum  __MinLpNum_mpf__

/* ========================================================================= */
/** @brief Given a double exp, compute \f$ e^{exp} \f$ and store it in the given
 * mpf_t number.
 * @param exp double exponent to be used.
 * @param num mpf_t number where to store the result.
 * */
#define mpf_EGlpNumEpow(num,exp) ({\
	unsigned int __i = 0;\
	int __lsgn = (exp)<0 ? 1:0;\
	mpf_t __ntmp,__res,__lexp,__err;\
	mpf_init_set_d(__lexp,exp);\
	if(__lsgn) mpf_neg(__lexp,__lexp);\
	mpf_init_set_ui(__ntmp,(unsigned long int)1);\
	mpf_init_set_ui(__res,(unsigned long int)1);\
	mpf_init(__err);\
	mpf_div(__err,__ntmp,__res);\
	while(mpf_cmp(__err,mpf_eps)>0)\
	{\
		mpf_mul(__ntmp,__ntmp,__lexp);\
		mpf_div_ui(__ntmp,__ntmp,(unsigned long int)(++__i));\
		mpf_add(__res,__res,__ntmp);\
		mpf_div(__err,__ntmp,__res);\
	}\
	if(__lsgn) mpf_ui_div(num,(unsigned long int)1,__res);\
	else mpf_set(num,__res);\
	mpf_clear(__ntmp);\
	mpf_clear(__res);\
	mpf_clear(__err);\
	mpf_clear(__lexp);})

/* ========================================================================= */
/** @brief Read from a string a number and store it in the given mpf_t, 
 * @return the number of chars readed from the input string */
#define mpf_EGlpNumReadStr(a,str) ({\
	int __i =0;\
	char __lpstr__[4096];\
	mp_exp_t __lexp;\
	mpf_set_str(a,str,10);\
	mpf_get_str(__lpstr__,&__lexp,10,(size_t)0,a);\
	__i = strlen(__lpstr__);\
	__i;})

/* ========================================================================= */
/** @brief given a mpf_t, write it to a string (to be allocated internally),
 * and return it. */
#define mpf_EGlpNumGetStr(a) ({\
	char *__out= 0;\
	mp_exp_t __lexp = 0;\
	size_t __pos = 0,__lo = (mpf_cmp_ui(a,(unsigned long)0) < 0) ? 1:0;\
	char __lpstr__[4096];\
	mpf_get_str(__lpstr__,&__lexp,10,(size_t)25,a);\
	__pos = strlen(__lpstr__) + ((size_t)(__lo + 2));\
	__out = EGsMalloc(char,__pos+((size_t)15));\
	if(__lo) sprintf(__out,"-0.");\
	else sprintf(__out,"0.");\
	sprintf(__out+__lo+2,"%s",__lpstr__+__lo);\
	if(__pos == 2) __out[1] = '\0';\
	else if(__lexp != 0)\
	{\
		__out[__pos-__lo] = 'e';\
		snprintf(__out+__pos+1-__lo,(size_t)(12-__lo),"%d",(int)__lexp);\
	}\
	__out;})

/* ========================================================================= */
/** @brief given an array of type mpf_t, free it, if the pointer is NULL
 * nothing happen. */
#define mpf_EGlpNumFreeArray(ea) ({\
	size_t __sz = __EGlpNumArraySize(ea);\
	mpf_t* __ptr__ = (ea);\
	while(__sz--) mpf_clear(__ptr__[__sz]);\
	__EGlpNumFreeArray(ea);})


/* ========================================================================= */
/** @brief Reallocate and initialize (if needed) 'size' elements of type 
 * mpf_t and return it, if no more memory, exit(1) */
#define mpf_EGlpNumReallocArray(lptr, lsize) ({ \
	mpf_t** __ptr__ = (lptr); \
	size_t __sz__ = (lsize); \
	size_t *__ntmp__ = (size_t *) *__ptr__; \
	size_t __psz__; \
	/* if no memory allocated before we just call the regular allocator */ \
	if (!*__ptr__) *__ptr__ = mpf_EGlpNumAllocArray (__sz__); \
	else \
	{ \
		/* first check that the previous size is not larger than the current */ \
		__ntmp__--; \
		__psz__ = __ntmp__[0]; \
		if (__psz__ < __sz__) \
		{ \
			/* now we have to do the reallocation */ \
			*__ptr__ = (mpf_t *) __ntmp__; \
			*__ptr__ = EGrealloc(*__ptr__, sizeof(mpf_t) * __sz__ +sizeof(size_t));\
			__ntmp__ = (size_t *) *__ptr__; \
			__ntmp__[0] = __sz__; \
			__ntmp__++; \
			*__ptr__ = (mpf_t *) __ntmp__; \
			for (; __psz__ < __sz__; __psz__++) mpf_init ((*__ptr__)[__psz__]); \
		} \
	} \
})

/* ========================================================================= */
/** @brief Allocate and initialize (if needed) 'size' elements of type mpf_t
 * and return it, if no more memory, exit(1) */
#define mpf_EGlpNumAllocArray(size) ({\
	size_t __i__ = (size);\
	mpf_t *__res = __EGlpNumAllocArray(mpf_t,__i__);\
	while(__i__--) mpf_init(__res[__i__]);\
	__res;})

/* ========================================================================= */
/** @brief set the given number pointer, set its value to the given double.
 * @param var mpf_t where we will store the double value.
 * @param dbl double value to be stored in 'var'.
 * @par Description:
 * This function is intended to set initial values to variables; note that the
 * double is a number and not a pointer to that value, be carefull with this
 * detail. Also, due to implementation details this function can't deal with
 * numbers above 1e158 or smaller than 1e-158. Note also that if the number is
 * writen in the form \f$x=\bar{x}\cdot 2^e\f$ with \f$0.5<|\bar{x}|<1\f$, 
 * then \f$\left|x-\frac{p}{q}\right|<2^{e-64}\f$.
 * */
#define mpf_EGlpNumSet(var, dbl) mpf_set_d(var,(double)(dbl))

/* ========================================================================= */
/** @brief Stores in the first number the ceil value of the second number, i.e.
 * EGlpNumCeil(a,b) <==> a= ceil(b) */
#define mpf_EGlpNumCeil(a, b) mpf_ceil(a,b)

/* ========================================================================= */
/** @brief Stores in the first number the floor value of the second number, i.e.
 * EGlpNumFloor(a,b) <==> a= floor(b) */
#define mpf_EGlpNumFloor(a, b) mpf_floor(a,b)

/* ========================================================================= */
/** @brief store the (multiplicative) inverse of a number to itself, i.e.
 * implement a = 1/a.
 * @param a the number to be inverted. */
#define mpf_EGlpNumInv(a) mpf_ui_div(a,(unsigned long int)1,a)

/* ========================================================================= */
/** @brief Compare if two numbers are equal within a maximum __error.
 * @param a mpf_t first number to compare.
 * @param b mpf_t second number to compare.
 * @return int one in success, zero oterwise.
 * @par Description:
 * Given two numbers 'a','b' return 1 if a == b, otherwise it return 0
 * */
#define mpf_EGlpNumIsEqqual(a,b) (mpf_cmp(a,b) == 0)

/* ========================================================================= */
/** @brief Compare if two numbers are equal within a maximum __error.
 * @param a mpf_t first number to compare.
 * @param b mpf_t second number to compare.
 * @param __error mpf_t maximum difference allowed between both
 * numbers.
 * @return int one in success, zero oterwise.
 * @par Description:
 * Given two numbers 'a','b' and a tolerance '__error',
 * return 1 if |a-b|<= __error, otherwise it return 0.
 * */
#define mpf_EGlpNumIsEqual(a,b,__error) ({\
	mpf_t __lpnum__;int __res__=0;mpf_init(__lpnum__);\
	mpf_sub (__lpnum__, a, b);\
	mpf_abs (__lpnum__, __lpnum__);\
	__res__=(mpf_cmp (__lpnum__, __error) <= 0);\
	mpf_clear(__lpnum__);\
	__res__;\
})

#define mpf_EGlpNumIsNeq(a,b,__error) ({\
	mpf_t __lpnum__;int __res__=0;mpf_init(__lpnum__);\
	mpf_sub (__lpnum__, a, b);\
	mpf_abs (__lpnum__, __lpnum__);\
	__res__=(mpf_cmp (__lpnum__, __error) > 0);\
	mpf_clear(__lpnum__);\
	__res__;\
})

#define mpf_EGlpNumIsNeqZero(a,__error) ({\
	mpf_t __lpnum__;int __res__=0;mpf_init(__lpnum__);\
	mpf_abs (__lpnum__, a);\
	__res__=(mpf_cmp (__lpnum__, __error) > 0);\
	mpf_clear(__lpnum__);\
	__res__;\
})

#define mpf_EGlpNumIsNeqqZero(a)     	(mpf_sgn(a))
#define mpf_EGlpNumIsNeqq(a,b)        (mpf_cmp(a,b)!=0)

/* ========================================================================= */
/** @brief test if the first number is bigger to the second number
 * @param a mpf_t the first number.
 * @param b mpf_t the second number
 * @return int one if success, zero otherwise.
 * @par Description:
 * Given two numbers 'a' and 'b', return one if a < b, zero
 * otherwise.
 * */
#define mpf_EGlpNumIsLess(a,b) (mpf_cmp(a,b) < 0)

/* ========================================================================= */
/** @brief test if a given number is greater than zero
 * @param a number to test
 * @return int one if success, zero otherwise.
 * */
#define mpf_EGlpNumIsGreatZero(a) (mpf_sgn(a) > 0)

/* ========================================================================= */
/** @brief test if a given number is less than zero
 * @param a number to test
 * @return int one if success, zero otherwise.
 * */
#define mpf_EGlpNumIsLessZero(a) (mpf_sgn(a) < 0)

/* ========================================================================= */
/** @brief test if the sum of the first two numbers is less thatn the third
 * number.
 * @param a mpf_t the first number.
 * @param b mpf_t the second number
 * @param c mpf_t the third number
 * @return int one if success, zero otherwise.
 * @par Description:
 * Given a,b, and c, return nonzero if (a + b < c), zero toherwise.
 * */
#define mpf_EGlpNumIsSumLess(a, b, c) ({\
	mpf_t __lpnum__;int __res__=0;mpf_init(__lpnum__);\
	mpf_add (__lpnum__, a, b);\
	__res__=(mpf_cmp (__lpnum__, c) < 0);\
	mpf_clear(__lpnum__);\
	__res__;\
})

/* ========================================================================= */
/** @brief test if the diference of the first two numbers is less thatn the 
 * third number.
 * @param a mpf_t the first number.
 * @param b mpf_t the second number
 * @param c mpf_t the third number
 * @return int one if success, zero otherwise.
 * @par Description:
 * Given a,b, and c, return nonzero if (a - b < c), zero toherwise.
 * */
#define mpf_EGlpNumIsDiffLess(a, b, c) ({\
	mpf_t __lpnum__;int __res__=0;mpf_init(__lpnum__);\
	mpf_sub (__lpnum__, a, b);\
	__res__=(mpf_cmp (__lpnum__, c) < 0);\
	mpf_clear(__lpnum__);\
	__res__;\
})

/* ========================================================================= */
/** @brief test if the first number is bigger to the second number
 * @param a mpf_t the first number.
 * @param b double the second number
 * @return int one if success, zero otherwise.
 * @par Description:
 * Given two numbers 'a' and 'b', return one if a < b, zero
 * otherwise.
 * */
#define mpf_EGlpNumIsLessDbl(a,b) (mpf_cmp_d(a,((double)(b))) < 0)

/* ========================================================================= */
/** @brief test if the first number is bigger to the second number
 * @param a mpf_t the first number.
 * @param b double the second number
 * @return int one if success, zero otherwise.
 * @par Description:
 * Given two numbers 'a' and 'b', return one if a > b, zero
 * otherwise.
 * */
#define mpf_EGlpNumIsGreaDbl(a,b) (mpf_cmp_d(a,((double)(b))) > 0)

/* ========================================================================= */
/** @brief test if the first number is bigger to the second number
 * @param a mpf_t the first number.
 * @param b mpf_t the second number
 * @return int one if success, zero otherwise.
 * @par Description:
 * Given two numbers 'a' and 'b', return one if a <= b, zero
 * otherwise.
 * */
#define mpf_EGlpNumIsLeq(a,b) (mpf_cmp(a,b) <= 0)

/* ========================================================================= */
/** @brief copy the value of the second number to the first.
 * @param a mpf_t source number (it won't change value).
 * @param b mpf_t source number (it won't change value).
 * @param c mpf_t denominator of the difference (it won't change value).
 * @param d mpf_t where to store the value .
 * @par Description:
 * Set @f$a = \frac{b - c}{d} @f$ */
#define mpf_EGlpNumCopyDiffRatio(a, b, c, d) ({\
	mpf_sub (a, b, c);\
	mpf_div (a, a, d);\
})

/* ========================================================================= */
/** @brief copy the value of the second number to the first.
 * @param a mpf_t source number (it won't change value).
 * @param b mpf_t source number (it won't change value).
 * @param dest mpf_t where to store the value stored in 'orig'.
 * @par Description:
 * Set dest = a - b */
#define mpf_EGlpNumCopyDiff(dest,a,b) mpf_sub(dest,a,b)

/* ========================================================================= */
/** @brief copy the value of the sum of the second and third parameter
 * @param a mpf_t source number (it won't change value).
 * @param b mpf_t source number (it won't change value).
 * @param dest mpf_t where to store the sum.
 * @par Description:
 * Set dest = a + b */
#define mpf_EGlpNumCopySum(dest,a,b) mpf_add(dest,a,b)

/* ========================================================================= */
/** @brief copy the value of the second number to the first.
 * @param orig mpf_t source number (it won't change value).
 * @param dest mpf_t where to store the value stored in 'orig'.
 * @par Description:
 * Given two numbers copy the values in 'orig', into 'dest'.
 * */
#define mpf_EGlpNumCopy(dest,orig) mpf_set(dest,orig)

/* ========================================================================= */
/** @brief change the fist number to the maximum between itself and the 
 * absolute value of the second.
 * @param orig mpf_t source number (it won't change value).
 * @param dest mpf_t where to store the value stored in 'orig'.
 * @par Description:
 * implement dest = max(dest,abs(orig))
 * */
#define mpf_EGlpNumSetToMaxAbs(dest, orig) ({\
	mpf_t __lpnum__;mpf_init(__lpnum__);\
	mpf_abs (__lpnum__, orig);\
	if (mpf_cmp (dest, __lpnum__) < 0) mpf_set (dest, __lpnum__);\
	mpf_clear(__lpnum__);})

#define mpf_EGlpNumSetToMinAbs(dest, orig) ({\
	mpf_t __lpnum__;mpf_init(__lpnum__);\
	mpf_abs (__lpnum__, orig);\
	if (mpf_cmp (dest, __lpnum__) > 0) mpf_set (dest, __lpnum__);\
	mpf_clear(__lpnum__);})

/* ========================================================================= */
/** @brief copy the square of the second argument, divided by the third 
 * argument into the first argument.
 * @param dest mpf_t where to store the result
 * @param orig mpf_t second parameter
 * @param den mpf_t third parameter
 * @par Description:
 * compute dest = (orig*orig)/den
 * */
#define mpf_EGlpNumCopySqrOver(dest, orig, den) ({\
	mpf_mul (dest, orig, orig);\
	mpf_div (dest, dest, den);\
})

/* ========================================================================= */
/** @brief copy the value of the absolute value of the second parameter to the 
 * first parameter.
 * @param orig mpf_t source number (it won't change value).
 * @param dest mpf_t where to store the absolute value stored
 * in 'orig'.
 * @par Description:
 * Given a number 'orig', copy its absolute value to 'dest'. i.e.
 * dest = |orig|
 * */
#define mpf_EGlpNumCopyAbs(dest,orig) mpf_abs(dest,orig)

/* ========================================================================= */
/** @brief copy minus the value of the second parameter to the 
 * first parameter.
 * @param orig mpf_t the source number (it won't change value).
 * @param dest mpf_t where to store minus the value stored
 * in 'orig'.
 * @par Description:
 * Given a number 'orig', copy minus the value to 'dest'. i.e.
 * dest = -orig
 * */
#define mpf_EGlpNumCopyNeg(dest,orig) mpf_neg(dest,orig)

/* ========================================================================= */
/** @brief Set des = op1/op2.
 * @param dest mpf_t where we will store the result.
 * @param op1 mpf_t numerator of the fraction (possibly non an integer)
 * @param op2 mpf_t denominator of the fraction (possibly non an integer)
 * @par Description:
 *  Set des = op1/op2
 * */
#define mpf_EGlpNumCopyFrac(dest,op1,op2) mpf_div(dest,op1,op2)

/* ========================================================================= */
/** @brief copy the first 'size' values in the second array to the first array.
 * @param orig mpf_t* pointer to the array from where we will copy the
 * values (it won't change value).
 * @param dest mpf_t* pointer to where to store the first 'size' values 
 * stored in 'orig'.
 * @param size unsigned int specifying how many values of 'orig' will be copied
 * onto 'dest'
 * @par Description:
 * This function is provided to (possible) make fast copies of arrays of
 * numbers, the arrays should be of length at least 'size', and the resulting
 * copy is absolutely independent froom the original, any change in one vale of
 * one array won't change values on the other array.
 * */
#define mpf_EGlpNumCopyArray(dest,orig,size) {\
	register unsigned int __i__ = (size);\
	for(;__i__--;)\
	{\
		mpf_set(dest[__i__],orig[__i__]);\
	}\
}

/* ========================================================================= */
/** @brief Sub to a given number the product of two numbers.
 * @param a mpf_t the number that we are going to Sub to.
 * @param b mpf_t value to be multiplyed.
 * @param c mpf_t value to be multiplyed.
 * @par Description:
 * This function implements a = a - b*c, and clearly don't change the value
 * stored in 'b' nor in 'c'.
 * */
#define mpf_EGlpNumSubInnProdTo(a, b, c) ({\
	mpf_t __lpnum__;mpf_init(__lpnum__);\
	mpf_mul (__lpnum__, b, c);\
	mpf_sub (a, a, __lpnum__);\
	mpf_clear(__lpnum__);\
})

/* ========================================================================= */
/** @brief Add to a given number the product of two numbers.
 * @param a mpf_t the number that we are going to add to.
 * @param b mpf_t value to be multiplyed.
 * @param c mpf_t value to be multiplyed.
 * @par Description:
 * This function implements a = a + b*c, and clearly don't change the value
 * stored in 'b' nor in 'c'.
 * */
#define mpf_EGlpNumAddInnProdTo(a, b, c) ({\
	mpf_t __lpnum__;mpf_init(__lpnum__);\
	mpf_mul (__lpnum__, b, c);\
	mpf_add (a, a, __lpnum__);\
	mpf_clear(__lpnum__);\
})

/* ========================================================================= */
/** @brief Substract to a given number the value of the second number.
 * @param a mpf_t the number that we are going to substract to.
 * @param b unsigned int value to be substracted to 'a'.
 * @par Description:
 * This function implements a = a - b, and clearly don't change the value
 * stored in 'b'.
 * */
#define mpf_EGlpNumSubUiTo(a,b) mpf_sub_ui(a,a,((unsigned long)(b)))

/* ========================================================================= */
/** @brief Add to a given number the value of the second number.
 * @param a mpf_t the number that we are going to add to.
 * @param b unsigned int value to be added to 'a'.
 * @par Description:
 * This function implements a = a + b, and clearly don't change the value
 * stored in 'b'.
 * */
#define mpf_EGlpNumAddUiTo(a,b) mpf_add_ui(a,a,((unsigned long)(b)))

/* ========================================================================= */
/** @brief Add to a given number the value of the second number.
 * @param a mpf_t the number that we are going to add to.
 * @param b mpf_t value to be added to 'a'.
 * @par Description:
 * This function implements a = a + b, and clearly don't change the value
 * stored in 'b'.
 * */
#define mpf_EGlpNumAddTo(a,b) mpf_add(a,a,b)

/* ========================================================================= */
/** @brief Substract to a given number the value of the second number.
 * @param a mpf_t the number that we are going to substract
 * from.
 * @param b mpf_t value to be substracted to 'a'.
 * @par Description:
 * This function implements a = a - b, and clearly don't change the value
 * stored in 'b'.
 * */
#define mpf_EGlpNumSubTo(a,b) mpf_sub(a,a,b)

/* ========================================================================= */
/** @brief Multiply a given number by the value of the second number.
 * @param a mpf_t the number that we are going to multiply by
 * the second number and store the result.
 * @param b mpf_t value to be multyply to 'a'.
 * @par Description:
 * This function implements a = a * b, and clearly don't change the value
 * stored in 'b'.
 * */
#define mpf_EGlpNumMultTo(a,b) mpf_mul(a,a,b)

/* ========================================================================= */
/** @brief Divide a given number by the value of the second number.
 * @param a mpf_t the number that we are going to divide by
 * the second number and store the result.
 * @param b mpf_t value to be divide to 'a'.
 * @par Description:
 * This function implements a = a / b, and clearly don't change the value
 * stored in 'b'.
 * */
#define mpf_EGlpNumDivTo(a,b) mpf_div(a,a,b)

/* ========================================================================= */
/** @brief Divide a given number by the value of the second number.
 * @param a mpf_t the number that we are going to divide by
 * the second number and store the result.
 * @param b unsigned int value to be divided to 'a'.
 * @par Description:
 * This function implements a = a / b, and don't change the value
 * stored in 'b'.
 * */
#define mpf_EGlpNumDivUiTo(a,b) mpf_div_ui(a,a,((unsigned long)(b)))

/* ========================================================================= */
/** @brief Multiply a given number by the value of the second number.
 * @param a mpf_t the number that we are going to multiply by
 * the second number and store the result.
 * @param b unsigned int value to be multyply to 'a'.
 * @par Description:
 * This function implements a = a * b, and clearly don't change the value
 * stored in 'b'.
 * */
#define mpf_EGlpNumMultUiTo(a,b) mpf_mul_ui(a,a,((unsigned long)(b)))

/* ========================================================================= */
/** @brief Reset the value of the pointed number to zero.
 * @param a mpf_t the value to be set to zero.
 * @par Descrpition:
 * Reset a to zero, i.e. implements a = 0;
 * */
#define mpf_EGlpNumZero(a) mpf_set_ui(a,(unsigned long)0)

/* ========================================================================= */
/** @brief Reset the value of the pointed number to one.
 * @param a mpf_t value to be set to one.
 * @par Descrpition:
 * Reset a to zero, i.e. implements a = 1;
 * */
#define mpf_EGlpNumOne(a) mpf_set_ui(a,(unsigned long)1)

/* ========================================================================= */
/** @brief Change the sign of the number.
 * @param a mpf_t number we will change sign.
 * @par Descrpition:
 * Change the sign of the given number, i.e. implements a = -a
 * */
#define mpf_EGlpNumSign(a) mpf_neg(a,a)

/* ========================================================================= */
/** @brief return the closest double value of the given pointer number.
 * @param a mpf_t number that we will be transformed to double.
 * @return double the closest double representation of the given number.
 * par Description:
 * return the double number closest in value to the value stored in a.
 * */
#define mpf_EGlpNumToLf(a) mpf_get_d(a)

/* ========================================================================= */
/** @brief initialize the internal memory of a given variable */
#define mpf_EGlpNumInitVar(a) mpf_init(a)

/* ========================================================================= */
/** @brief free the internal memory of a given variable */
#define mpf_EGlpNumClearVar(a) mpf_clear(a)

/* ========================================================================= */
/** @} */
#endif
