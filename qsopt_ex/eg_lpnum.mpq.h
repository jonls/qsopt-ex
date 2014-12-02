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
#ifndef __EG_LPNUM_MPQ__
#define __EG_LPNUM_MPQ__

#include <stdlib.h>
#include <string.h>

#include <gmp.h>

#include "eg_lpnum.h"

/** @file
 * @brief Interface for rational implementation of the EGlpNum_t type.
 * @par History
 * - 2006-02-01
 * 					- Add verbosity flag for continued fraction conversions.
 * @ingroup EGlpNum */
/** @addtogroup EGlpNum */
/** @{ */
/* ========================================================================= */
/** extern definitions of constaants for different set-ups */
extern const mpq_t __zeroLpNum_mpq__;
extern const mpq_t __oneLpNum_mpq__;
extern const mpq_t __MaxLpNum_mpq__;
extern const mpq_t __MinLpNum_mpq__;
#define mpq_zeroLpNum __zeroLpNum_mpq__
#define mpq_oneLpNum  __oneLpNum_mpq__
#define mpq_epsLpNum  __zeroLpNum_mpq__
#define mpq_MaxLpNum  __MaxLpNum_mpq__
#define mpq_MinLpNum  __MinLpNum_mpq__

/* ========================================================================= */
/** @brief This function read a number in float form and store it in an mpq_t
 * variable, returning how many chars read to create the number, the twist is
 * that it does an 'exact' transformation, in the sense that 0.33333333 will be
 * stored as 33333333/100000000.
 * @param str input string.
 * @param var variable where we will store the number as rational.
 * @return number of reade chars.
 * @par Descriptiom:
 * If the input string doesn't contain a number, 'var' will be set to zero, and
 * the number of readed chars will be zero, we assume that there are no leading
 * empty spaces, (nor tabs), no eschape characters, and trailing zeros are
 * considered as readed, but no the following spaces. 
 * @note
 * This function will only read number in decimal base, in a future release we
 * may include a more general reader. */
int mpq_EGlpNumReadStrXc (mpq_t var,
													char const *str);

/* ========================================================================= */
/** @brief Read from a string a number and store it in the given mpq_t, return
 * the number of chars readed from the input string */
#define mpq_EGlpNumReadStr(a,str) mpq_EGlpNumReadStrXc(a,str)

/* ========================================================================= */
/** @brief given a mpq_t, write it to a string (to be allocated internally), 
 * and return it. */
#define mpq_EGlpNumGetStr(a) ({\
	const size_t __sz = mpz_sizeinbase(mpq_numref(a),10) + mpz_sizeinbase(mpq_denref(a),10) + 3;\
	char *__str=EGsMalloc(char,__sz);\
	mpq_get_str(__str,10,a);})

/* ========================================================================= */
/** @brief given an array of type mpq_t, free it, if the pointer is NULL
 * nothing happen. */
#define mpq_EGlpNumFreeArray(ea) ({\
	size_t __sz = __EGlpNumArraySize(ea);\
	mpq_t* __ptr__ = (ea);\
	while(__sz--) mpq_clear(__ptr__[__sz]);\
	__EGlpNumFreeArray(ea);})


/* ========================================================================= */
/** @brief Reallocate and initialize (if needed) 'size' elements of type 
 * mpq_t and return it, if no more memory, exit(1) */
#define mpq_EGlpNumReallocArray(lptr, lsize) ({\
	mpq_t **__ptr__ = (lptr);\
	size_t *__ntmp__ = (size_t *) *__ptr__, __sz__ = (lsize);\
	size_t __psz__;\
	/* if no memory allocated before we just call the regular allocator */\
	if (!*__ptr__)\
		*__ptr__ = mpq_EGlpNumAllocArray (__sz__);\
	else\
	{\
		/* first check that the previous size is not larger than the current */\
		__ntmp__--;\
		__psz__ = __ntmp__[0];\
		if (__psz__ < __sz__)\
		{\
			/* now we have to do the reallocation */\
			*__ptr__ = (mpq_t *) __ntmp__;\
			*__ptr__ = EGrealloc(*__ptr__,sizeof(mpq_t) * __sz__ + sizeof(size_t));\
			__ntmp__ = (size_t *) *__ptr__;\
			__ntmp__[0] = __sz__;\
			__ntmp__++;\
			*__ptr__ = (mpq_t *) __ntmp__;\
			for (; __psz__ < __sz__; __psz__++) mpq_init ((*__ptr__)[__psz__]);\
		}\
	}\
})

/* ========================================================================= */
/** @brief Allocate and initialize (if needed) 'size' elements of type mpq_t
 * and return it, if no more memory, exit(1) */
#define mpq_EGlpNumAllocArray(size) ({\
	size_t __i__ = (size);\
	mpq_t *__res = __EGlpNumAllocArray(mpq_t,__i__);\
	while(__i__--) mpq_init(__res[__i__]);\
	__res;})

/* ========================================================================= */
/** @brief set the given rational number , to the value to the value of the 
 * given mpf_t, this conversion is done using the continuous fraction method.
 * @param var mpq_t where we will store the value.
 * @param flt mpf_t value to be stored in 'var'.
 * @par Description:
 * This function is intended to set initial values to variables. Note also 
 * that if the number is
 * writen in the form \f$x=\bar{x}\cdot 2^e\f$ with \f$0.5<|\bar{x}|<1\f$, 
 * then \f$\left|x-\frac{p}{q}\right|<2^{e-EGLPNUM_PRECISION}\f$.
 * */
void mpq_EGlpNumSet_mpf (mpq_t var,
												 mpf_t flt);

/* ========================================================================= */
/** @brief set the given number pointer, set its value to the given double.
 * @param var mpq_t where we will store the double value.
 * @param dbl double value to be stored in 'var'.
 * @par Description:
 * This function is intended to set initial values to variables; note that the
 * double is a number and not a pointer to that value, be carefull with this
 * detail. Also, due to implementation details this function can't deal with
 * numbers above 1e158 or smaller than 1e-158. Note also that if the number is
 * writen in the form \f$x=\bar{x}\cdot 2^e\f$ with \f$0.5<|\bar{x}|<1\f$, 
 * then \f$\left|x-\frac{p}{q}\right|<2^{e-64}\f$.
 * */
void mpq_EGlpNumSet (mpq_t var,
										 const double dbl);

/* ========================================================================= */
/** @brief Stores in the first number the ceil value of the second number, i.e.
 * EGlpNumCeil(a,b) <==> a= ceil(b) */
#define mpq_EGlpNumCeil(a, b) ({\
	mpz_cdiv_q (mpq_numref (a), mpq_numref (b), mpq_denref (b));\
	mpz_set_ui (mpq_denref (a), (unsigned long int)1);\
})

/* ========================================================================= */
/** @brief Stores in the first number the floor value of the second number, i.e.
 * EGlpNumFloor(a,b) <==> a= floor(b) */
#define mpq_EGlpNumFloor(a, b) ({\
	mpz_fdiv_q (mpq_numref (a), mpq_numref (b), mpq_denref (b));\
	mpz_set_ui (mpq_denref (a), (unsigned long int)1);\
})

/* ========================================================================= */
/** @brief store the (multiplicative) inverse of a number to itself, i.e.
 * implement a = 1/a.
 * @param a the number to be inverted. */
#define mpq_EGlpNumInv(a) mpq_inv(a,a)

/* ========================================================================= */
/** @brief Compare if two numbers are equal within a maximum error.
 * @param a mpq_t first number to compare.
 * @param b mpq_t second number to compare.
 * @return int one in success, zero oterwise.
 * @par Description:
 * Given two numbers 'a','b' return 1 if a == b, otherwise it return 0
 * */
#define mpq_EGlpNumIsEqqual(a,b) (mpq_equal(a,b))

/* ========================================================================= */
/** @brief Compare if two numbers are equal within a maximum error.
 * @param a mpq_t first number to compare.
 * @param b mpq_t second number to compare.
 * @param error mpq_t maximum difference allowed between both
 * numbers.
 * @return int one in success, zero oterwise.
 * @par Description:
 * Given two numbers 'a','b' and a tolerance 'error',
 * return 1 if |a-b|<= error, otherwise it return 0.
 * */
#define mpq_EGlpNumIsEqual(a,b,error) (mpq_equal(a,b))
#define mpq_EGlpNumIsNeq(a,b,error)   (!(mpq_equal(a,b)))
#define mpq_EGlpNumIsNeqq(a,b)        (!(mpq_equal(a,b)))
#define mpq_EGlpNumIsNeqZero(a,error) (mpq_sgn(a))
#define mpq_EGlpNumIsNeqqZero(a)      (mpq_sgn(a))

/* ========================================================================= */
/** @brief test if the first number is bigger to the second number
 * @param a mpq_t the first number.
 * @param b mpq_t the second number
 * @return int one if success, zero otherwise.
 * @par Description:
 * Given two numbers 'a' and 'b', return one if a < b, zero
 * otherwise.
 * */
#define mpq_EGlpNumIsLess(a,b) (mpq_cmp(a,b) < 0)

/* ========================================================================= */
/** @brief test if a number is greater than zero
 * @param a number to test
 * @return int one if success, zero otherwise.
 * */
#define mpq_EGlpNumIsGreatZero(a) (mpq_sgn(a) > 0)

/* ========================================================================= */
/** @brief test if a number is less than zero
 * @param a number to test
 * @return int one if success, zero otherwise.
 * */
#define mpq_EGlpNumIsLessZero(a) (mpq_sgn(a) < 0)

/* ========================================================================= */
/** @brief test if the sum of the first two numbers is less thatn the third
 * number.
 * @param a mpq_t the first number.
 * @param b mpq_t the second number
 * @param c mpq_t the third number
 * @return int one if success, zero otherwise.
 * @par Description:
 * Given a,b, and c, return nonzero if (a + b < c), zero toherwise.
 * */
#define mpq_EGlpNumIsSumLess(a, b, c) ({\
	mpq_t __lpnum__;int __res__=0;mpq_init(__lpnum__);\
	mpq_add(__lpnum__,a,b);\
	__res__=(mpq_cmp (__lpnum__, c) < 0);\
	mpq_clear(__lpnum__);\
	__res__;\
})

/* ========================================================================= */
/** @brief test if the diference of the first two numbers is less thatn the 
 * third number.
 * @param a mpq_t the first number.
 * @param b mpq_t the second number
 * @param c mpq_t the third number
 * @return int one if success, zero otherwise.
 * @par Description:
 * Given a,b, and c, return nonzero if (a - b < c), zero toherwise.
 * */
#define mpq_EGlpNumIsDiffLess(a, b, c) ({\
	mpq_t __lpnum__;int __res__=0;mpq_init(__lpnum__);\
	mpq_sub (__lpnum__, a, b);\
	__res__=(mpq_cmp (__lpnum__, c) < 0);\
	mpq_clear(__lpnum__);\
	__res__;\
})

/* ========================================================================= */
/** @brief test if the first number is bigger to the second number
 * @param a mpq_t the first number.
 * @param b double the second number
 * @return int one if success, zero otherwise.
 * @par Description:
 * Given two numbers 'a' and 'b', return one if a < b, zero
 * otherwise.
 * */
#define mpq_EGlpNumIsLessDbl(a,b) (mpq_get_d(a) < b)

/* ========================================================================= */
/** @brief test if the first number is bigger to the second number
 * @param a mpq_t the first number.
 * @param b double the second number
 * @return int one if success, zero otherwise.
 * @par Description:
 * Given two numbers 'a' and 'b', return one if a > b, zero
 * otherwise.
 * */
#define mpq_EGlpNumIsGreaDbl(a,b) (mpq_get_d(a) > b)

/* ========================================================================= */
/** @brief test if the first number is bigger to the second number
 * @param a mpq_t the first number.
 * @param b mpq_t the second number
 * @return int one if success, zero otherwise.
 * @par Description:
 * Given two numbers 'a' and 'b', return one if a <= b, zero
 * otherwise.
 * */
#define mpq_EGlpNumIsLeq(a,b) (mpq_cmp(a,b) <= 0)

/* ========================================================================= */
/** @brief copy the value of the second number to the first.
 * @param a mpq_t source number (it won't change value).
 * @param b mpq_t source number (it won't change value).
 * @param c mpq_t denominator of the difference (it won't change value).
 * @param d mpq_t where to store the value .
 * @par Description:
 * Set @f$a = \frac{b - c}{d} @f$ */
#define mpq_EGlpNumCopyDiffRatio(a, b, c, d) ({\
	mpq_sub (a, b, c);\
	mpq_div (a, a, d);\
})

/* ========================================================================= */
/** @brief copy the value of the second number to the first.
 * @param a mpq_t source number (it won't change value).
 * @param b mpq_t source number (it won't change value).
 * @param dest mpq_t where to store the value stored in 'orig'.
 * @par Description:
 * Set dest = a - b */
#define mpq_EGlpNumCopyDiff(dest,a,b) mpq_sub(dest,a,b)

/* ========================================================================= */
/** @brief copy the value of the sum of the second and third parameter
 * @param a mpq_t source number (it won't change value).
 * @param b mpq_t source number (it won't change value).
 * @param dest mpq_t where to store the sum.
 * @par Description:
 * Set dest = a + b */
#define mpq_EGlpNumCopySum(dest,a,b) mpq_add(dest,a,b)

/* ========================================================================= */
/** @brief copy the value of the second number to the first.
 * @param orig mpq_t source number (it won't change value).
 * @param dest mpq_t where to store the value stored in 'orig'.
 * @par Description:
 * Given two numbers copy the values in 'orig', into 'dest'.
 * */
#define mpq_EGlpNumCopy(dest,orig) mpq_set(dest,orig)

/* ========================================================================= */
/** @brief change the fist number to the maximum between itself and the 
 * absolute value of the second.
 * @param orig mpq_t source number (it won't change value).
 * @param dest mpq_t where to store the value stored in 'orig'.
 * @par Description:
 * implement dest = max(dest,abs(orig))
 * */
#define mpq_EGlpNumSetToMaxAbs(dest, orig) ({\
	mpq_t __lpnum__;mpq_init(__lpnum__);\
	mpq_abs (__lpnum__, orig);\
	if (mpq_cmp (dest, __lpnum__) < 0)\
		mpq_set (dest, __lpnum__);\
	mpq_clear(__lpnum__);\
})

#define mpq_EGlpNumSetToMinAbs(dest, orig) ({\
	mpq_t __lpnum__;mpq_init(__lpnum__);\
	mpq_abs (__lpnum__, orig);\
	if (mpq_cmp (dest, __lpnum__) > 0)\
		mpq_set (dest, __lpnum__);\
	mpq_clear(__lpnum__);\
})

/* ========================================================================= */
/** @brief copy the square of the second argument, divided by the third 
 * argument into the first argument.
 * @param dest mpq_t where to store the result
 * @param orig mpq_t second parameter
 * @param den mpq_t third parameter
 * @par Description:
 * compute dest = (orig*orig)/den
 * */
#define mpq_EGlpNumCopySqrOver(dest, orig, den) ({\
	mpq_mul (dest, orig, orig);\
	mpq_div (dest, dest, den);\
})

/* ========================================================================= */
/** @brief copy the value of the absolute value of the second parameter to the 
 * first parameter.
 * @param orig mpq_t source number (it won't change value).
 * @param dest mpq_t where to store the absolute value stored
 * in 'orig'.
 * @par Description:
 * Given a number 'orig', copy its absolute value to 'dest'. i.e.
 * dest = |orig|
 * */
#define mpq_EGlpNumCopyAbs(dest,orig) mpq_abs(dest,orig)

/* ========================================================================= */
/** @brief copy minus the value of the second parameter to the 
 * first parameter.
 * @param orig mpq_t the source number (it won't change value).
 * @param dest mpq_t where to store minus the value stored
 * in 'orig'.
 * @par Description:
 * Given a number 'orig', copy minus the value to 'dest'. i.e.
 * dest = -orig
 * */
#define mpq_EGlpNumCopyNeg(dest,orig) mpq_neg(dest,orig)

/* ========================================================================= */
/** @brief Set des = op1/op2.
 * @param dest mpq_t where we will store the result.
 * @param op1 mpq_t numerator of the fraction (possibly non an integer)
 * @param op2 mpq_t denominator of the fraction (possibly non an integer)
 * @par Description:
 *  Set des = op1/op2
 * */
#define mpq_EGlpNumCopyFrac(dest,op1,op2) mpq_div(dest,op1,op2)

/* ========================================================================= */
/** @brief copy the first 'size' values in the second array to the first array.
 * @param orig mpq_t* pointer to the array from where we will copy the
 * values (it won't change value).
 * @param dest mpq_t* pointer to where to store the first 'size' values 
 * stored in 'orig'.
 * @param size unsigned int specifying how many values of 'orig' will be copied
 * onto 'dest'
 * @par Description:
 * This function is provided to (possible) make fast copies of arrays of
 * numbers, the arrays should be of length at least 'size', and the resulting
 * copy is absolutely independent froom the original, any change in one vale of
 * one array won't change values on the other array.
 * */
#define mpq_EGlpNumCopyArray(dest,orig,size) {\
	register unsigned int __i__ = size;\
	for(;__i__--;)\
	{\
		mpq_set(dest[__i__],orig[__i__]);\
	}\
}

/* ========================================================================= */
/** @brief Sub to a given number the product of two numbers.
 * @param a mpq_t the number that we are going to Sub to.
 * @param b mpq_t value to be multiplyed.
 * @param c mpq_t value to be multiplyed.
 * @par Description:
 * This function implements a = a - b*c, and clearly don't change the value
 * stored in 'b' nor in 'c'.
 * */
#define mpq_EGlpNumSubInnProdTo(a, b, c) ({\
	mpq_t __lpnum__;mpq_init(__lpnum__);\
	mpq_mul (__lpnum__, b, c);\
	mpq_sub (a, a, __lpnum__);\
	mpq_clear(__lpnum__);\
})

/* ========================================================================= */
/** @brief Add to a given number the product of two numbers.
 * @param a mpq_t the number that we are going to add to.
 * @param b mpq_t value to be multiplyed.
 * @param c mpq_t value to be multiplyed.
 * @par Description:
 * This function implements a = a + b*c, and clearly don't change the value
 * stored in 'b' nor in 'c'.
 * */
#define mpq_EGlpNumAddInnProdTo(a, b, c) ({\
	mpq_t __lpnum__;mpq_init(__lpnum__);\
	mpq_mul (__lpnum__, b, c);\
	mpq_add (a, a, __lpnum__);\
	mpq_clear(__lpnum__);\
})

/* ========================================================================= */
/** @brief Substract to a given number the value of the second number.
 * @param a mpq_t the number that we are going to substract to.
 * @param b unsigned int value to be substracted to 'a'.
 * @par Description:
 * This function implements a = a - b, and clearly don't change the value
 * stored in 'b'.
 * */
#define mpq_EGlpNumSubUiTo(a,b) mpz_submul_ui(mpq_numref(a),mpq_denref(a),(unsigned long int)b)

/* ========================================================================= */
/** @brief Add to a given number the value of the second number.
 * @param a mpq_t the number that we are going to add to.
 * @param b unsigned int value to be added to 'a'.
 * @par Description:
 * This function implements a = a + b, and clearly don't change the value
 * stored in 'b'.
 * */
#define mpq_EGlpNumAddUiTo(a,b) mpz_addmul_ui(mpq_numref(a),mpq_denref(a),(unsigned long int)b)

/* ========================================================================= */
/** @brief Add to a given number the value of the second number.
 * @param a mpq_t the number that we are going to add to.
 * @param b mpq_t value to be added to 'a'.
 * @par Description:
 * This function implements a = a + b, and clearly don't change the value
 * stored in 'b'.
 * */
#define mpq_EGlpNumAddTo(a,b) mpq_add(a,a,b)

/* ========================================================================= */
/** @brief Substract to a given number the value of the second number.
 * @param a mpq_t the number that we are going to substract
 * from.
 * @param b mpq_t value to be substracted to 'a'.
 * @par Description:
 * This function implements a = a - b, and clearly don't change the value
 * stored in 'b'.
 * */
#define mpq_EGlpNumSubTo(a,b) mpq_sub(a,a,b)

/* ========================================================================= */
/** @brief Multiply a given number by the value of the second number.
 * @param a mpq_t the number that we are going to multiply by
 * the second number and store the result.
 * @param b mpq_t value to be multyply to 'a'.
 * @par Description:
 * This function implements a = a * b, and clearly don't change the value
 * stored in 'b'.
 * */
#define mpq_EGlpNumMultTo(a,b) mpq_mul(a,a,b)

/* ========================================================================= */
/** @brief Divide a given number by the value of the second number.
 * @param a mpq_t the number that we are going to divide by
 * the second number and store the result.
 * @param b mpq_t value to be divide to 'a'.
 * @par Description:
 * This function implements a = a / b, and clearly don't change the value
 * stored in 'b'.
 * */
#define mpq_EGlpNumDivTo(a,b) mpq_div(a,a,b)

/* ========================================================================= */
/** @brief Divide a given number by the value of the second number.
 * @param a mpq_t the number that we are going to divide by
 * the second number and store the result.
 * @param b unsigned int value to be divided to 'a'.
 * @par Description:
 * This function implements a = a / b, and don't change the value
 * stored in 'b'.
 * */
#define mpq_EGlpNumDivUiTo(a,b) do{mpz_mul_ui(mpq_denref(a),mpq_denref(a),(unsigned long)b);mpq_canonicalize(a);}while(0)

/* ========================================================================= */
/** @brief Multiply a given number by the value of the second number.
 * @param a mpq_t the number that we are going to multiply by
 * the second number and store the result.
 * @param b unsigned int value to be multyply to 'a'.
 * @par Description:
 * This function implements a = a * b, and clearly don't change the value
 * stored in 'b'.
 * */
#define mpq_EGlpNumMultUiTo(a,b) do{mpz_mul_ui(mpq_numref(a),mpq_numref(a),(unsigned long int)b);mpq_canonicalize(a);}while(0)

/* ========================================================================= */
/** @brief Reset the value of the pointed number to zero.
 * @param a mpq_t the value to be set to zero.
 * @par Descrpition:
 * Reset a to zero, i.e. implements a = 0;
 * */
#define mpq_EGlpNumZero(a) mpq_set_ui(a,(unsigned long)0,(unsigned long)1)

/* ========================================================================= */
/** @brief Reset the value of the pointed number to one.
 * @param a mpq_t value to be set to one.
 * @par Descrpition:
 * Reset a to zero, i.e. implements a = 1;
 * */
#define mpq_EGlpNumOne(a) mpq_set_ui(a,(unsigned long)1,(unsigned long)1)

/* ========================================================================= */
/** @brief Change the sign of the number.
 * @param a mpq_t number we will change sign.
 * @par Descrpition:
 * Change the sign of the given number, i.e. implements a = -a
 * */
#define mpq_EGlpNumSign(a) mpq_neg(a,a)

/* ========================================================================= */
/** @brief return the closest double value of the given pointer number.
 * @param a mpq_t number that we will be transformed to double.
 * @return double the closest double representation of the given number.
 * par Description:
 * return the double number closest in value to the value stored in a.
 * */
#define mpq_EGlpNumToLf(a) mpq_get_d(a)

/* ========================================================================= */
/** @brief initialize the internal memory of a given variable */
#define mpq_EGlpNumInitVar(a) mpq_init(a)

/* ========================================================================= */
/** @brief free the internal memory of a given variable */
#define mpq_EGlpNumClearVar(a) mpq_clear(a)

/* ========================================================================= */
/** @} */
#endif
