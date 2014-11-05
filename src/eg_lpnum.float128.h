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
#ifndef __EG_LPNUM_FLOAT128__
#define __EG_LPNUM_FLOAT128__
#ifdef HAVE_SOFTFLOAT
#if HAVE_SOFTFLOAT
#include "softfloat.h"
#include "eg_lpnum.h"
/** @file
 * @ingroup EGlpNum */
/** @addtogroup EGlpNum */
/** @{ */
/* ========================================================================= */
/** @brief This is the smallest difference (with the current precision) that can
 * be distinguished betwen 1.0 and it's clossest representable number, in some
 * sense it is the absolute minimum epsilon for comparisons */
extern float128 float128_eps;

/* ========================================================================= */
/** @brief union to conver C doubles to float64 types as provided by SoftFloat.
 * This is needed because the default conversion does not apply. */
typedef union
{
	float64 fl;
	double dbl;
}
f64dbl_t;

/* ========================================================================= */
/** extern definitions of constants for different set-ups */
extern const float128 __zeroLpNum_float128__;
extern const float128 __oneLpNum_float128__;
extern const float128 __MaxLpNum_float128__;
extern const float128 __MinLpNum_float128__;
#define float128_zeroLpNum __zeroLpNum_float128__
#define float128_oneLpNum  __oneLpNum_float128__
#define float128_epsLpNum  float128_eps
#define float128_MaxLpNum	__MaxLpNum_float128__
#define float128_MinLpNum	__MinLpNum_float128__

/* ========================================================================= */
/** @brief Given a double exp, compute \f$ e^{exp} \f$ and store it in the given
 * float128 number.
 * @param exp double exponent to be used.
 * @param num float128 number where to store the result.
 * */
#define float128_EGlpNumEpow(num,exp) ({\
	unsigned int __i = 0;\
	const f64dbl_t __exp = {.dbl = (exp)};\
	int __lsgn = __exp.dbl <0 ? 1:0;\
	float128 __ntmp,__res,__lexp,__err,__fint;\
	__lext = float64_to_float128(__exp.fl);\
	__lexp.high &= 0x7fffffffffffffffLL;\
	__ntmp = int32_to_float128(1);\
	__res = int32_to_float128(1);\
	__err = float128_div(__ntmp,__res);\
	while(float128_lt(float128_eps,__err))\
	{\
		__ntmp = float128_mul(__ntmp,__lexp);\
		__fint = int32_to_float128(++__i);\
		__ntmp = float128_div(__ntmp,__fint);\
		__res = float128_add(__res,__ntmp);\
		__err = float128_div(__ntmp,__res);\
	}\
	if(__lsgn) num = float128_div(float128_oneLpNum,__res);\
	else num = __res;})

/* ========================================================================= */
/** @brief Read from a string a number and store it in the given float128, 
 * @return the number of chars readed from the input string */
#define float128_EGlpNumReadStr(a,str) ({\
	int __i =0;\
	f64dbl_t __tmp;\
	sscanf(str,"%lf%n",&(__tmp.dbl),&__i);\
	(a) = float64_to_float128(__tmp.fl);\
	__i;})

/* ========================================================================= */
/** @brief given a float128, write it to a string (to be allocated internally), 
 * and return it. */
#define float128_EGlpNumGetStr(a) ({\
	const f64dbl_t __tmp = {.fl = float128_to_float64(a)};\
	char *__str=0;\
	size_t __i = snprintf(__str,(size_t)0,"%.7lg",__tmp.dbl);\
	__str = EGsMalloc(char,__i+1);\
	snprintf(__str,__i+1,"%.7lg",__tmp.dbl);\
	__str;})

/* ========================================================================= */
/** @brief given an array of type float128, free it, if the pointer is NULL
 * nothing happen. */
#define float128_EGlpNumFreeArray(ea) __EGlpNumFreeArray(ea)

/* ========================================================================= */
/** @brief Reallocate and initialize (if needed) 'size' elements of type 
 * float128 and return it, if no more memory, exit(1) */
#define float128_EGlpNumReallocArray(lptr, lsize) ({ \
	float128** __ptr__ = (lptr); \
	size_t __sz__ = (lsize); \
	size_t *__ntmp__ = (size_t *) *__ptr__; \
	size_t __psz__; \
	/* if no memory allocated before we just call the regular allocator */ \
	if (!*__ptr__) *__ptr__ = float128_EGlpNumAllocArray (__sz__); \
	else \
	{ \
		/* first check that the previous size is not larger than the current */ \
		__ntmp__--; \
		__psz__ = __ntmp__[0]; \
		if (__psz__ < __sz__) \
		{ \
			/* now we have to do the reallocation */ \
			*__ptr__ = (float128 *) __ntmp__; \
			*__ptr__ = EGrealloc(*__ptr__, sizeof(float128) * __sz__ +sizeof(size_t));\
			__ntmp__ = (size_t *) *__ptr__; \
			__ntmp__[0] = __sz__; \
			__ntmp__++; \
			*__ptr__ = (float128 *) __ntmp__; \
			for (; __psz__ < __sz__; __psz__++) (*__ptr__)[__psz__] = (float128){0,0}; \
		} \
	} \
})

/* ========================================================================= */
/** @brief Allocate and initialize (if needed) 'size' elements of type float128
 * and return it, if no more memory, exit(1) */
#define float128_EGlpNumAllocArray(size) __EGlpNumAllocArray(float128,(size))

/* ========================================================================= */
/** @brief set the given number pointer, set its value to the given double.
 * @param var float128 where we will store the double value.
 * @param edbl double value to be stored in 'var'.
 * @par Description:
 * This function is intended to set initial values to variables; note that the
 * double is a number and not a pointer to that value, be carefull with this
 * detail. Also, due to implementation details this function can't deal with
 * numbers above 1e158 or smaller than 1e-158. Note also that if the number is
 * writen in the form \f$x=\bar{x}\cdot 2^e\f$ with \f$0.5<|\bar{x}|<1\f$, 
 * then \f$\left|x-\frac{p}{q}\right|<2^{e-64}\f$.
 * */
#define float128_EGlpNumSet(var, edbl) ({\
	f64dbl_t __tmp = {.dbl = (edbl)};\
	var = float64_to_float128(__tmp.fl);})

/* ========================================================================= */
/** @brief Stores in the first number the ceil value of the second number, i.e.
 * EGlpNumCeil(a,b) <==> a= ceil(b) */
#define float128_EGlpNumCeil(a, b) do{\
	(a) = float128_round_to_int(b);\
	if(float128_lt(a,b)) (a) = float128_add((a),float128_oneLpNum);} while(0)

/* ========================================================================= */
/** @brief Stores in the first number the floor value of the second number, i.e.
 * EGlpNumFloor(a,b) <==> a= floor(b) */
#define float128_EGlpNumFloor(a, b) do{\
	(a) = float128_round_to_int(b);\
	if(float128_lt(b,a)) (a) = float128_sub((a),float128_oneLpNum);} while(0)

/* ========================================================================= */
/** @brief store the (multiplicative) inverse of a number to itself, i.e.
 * implement a = 1/a.
 * @param a the number to be inverted. */
#define float128_EGlpNumInv(a) ((a) = float128_div(float128_oneLpNum,a))

/* ========================================================================= */
/** @brief Compare if two numbers are equal within a maximum error.
 * @param a float128 first number to compare.
 * @param b float128 second number to compare.
 * @return int one in success, zero oterwise.
 * @par Description:
 * Given two numbers 'a','b' return 1 if a == b, otherwise it return 0
 * */
#define float128_EGlpNumIsEqqual(a,b) (float128_eq(a,b))

/* ========================================================================= */
/** @brief Compare if two numbers are equal within a maximum error.
 * @param a float128 first number to compare.
 * @param b float128 second number to compare.
 * @param error float128 maximum difference allowed between both
 * numbers.
 * @return int one in success, zero oterwise.
 * @par Description:
 * Given two numbers 'a','b' and a tolerance 'error',
 * return 1 if |a-b|<= error, otherwise it return 0.
 * */
#define float128_EGlpNumIsEqual(a,b,error) ({\
	float128 __lpnum__= float128_sub (a, b);\
	__lpnum__.high &= 0x7fffffffffffffffLL;\
	float128_le(__lpnum__, error);\
})

#define float128_EGlpNumIsNeq(a,b,error) ({\
	float128 __lpnum__= float128_sub(a, b);\
	__lpnum__.high &= 0x7fffffffffffffffLL;\
	float128_lt(error,__lpnum__);})

#define float128_EGlpNumIsNeqZero(a,error) ({\
	float128 __lpnum__= (a);\
	__lpnum__.high &= 0x7fffffffffffffffLL;\
	float128_lt(error,__lpnum__);\
})

#define float128_EGlpNumIsNeqqZero(a)     	(!float128_eq(a,float128_zeroLpNum))
#define float128_EGlpNumIsNeqq(a,b)        (!float128_eq(a,b))

/* ========================================================================= */
/** @brief test if the first number is bigger to the second number
 * @param a float128 the first number.
 * @param b float128 the second number
 * @return int one if success, zero otherwise.
 * @par Description:
 * Given two numbers 'a' and 'b', return one if a < b, zero
 * otherwise.
 * */
#define float128_EGlpNumIsLess(a,b) (float128_lt(a,b))

/* ========================================================================= */
/** @brief test if the first number is greater than zero
 * @param a number to compare.
 * @return int one if success, zero otherwise.
 * */
#define float128_EGlpNumIsGreatZero(a) (!float128_le(a,float128_zeroLpNum))

/* ========================================================================= */
/** @brief test if the first number is less than zero
 * @param a number to compare.
 * @return int one if success, zero otherwise.
 * */
#define float128_EGlpNumIsLessZero(a) (float128_lt(a,float128_zeroLpNum))

/* ========================================================================= */
/** @brief test if the sum of the first two numbers is less thatn the third
 * number.
 * @param a float128 the first number.
 * @param b float128 the second number
 * @param c float128 the third number
 * @return int one if success, zero otherwise.
 * @par Description:
 * Given a,b, and c, return nonzero if (a + b < c), zero toherwise.
 * */
#define float128_EGlpNumIsSumLess(a, b, c) ({\
	float128 __lpnum__= float128_add( a, b);\
	float128_lt(__lpnum__, c);\
})

/* ========================================================================= */
/** @brief test if the diference of the first two numbers is less thatn the 
 * third number.
 * @param a float128 the first number.
 * @param b float128 the second number
 * @param c float128 the third number
 * @return int one if success, zero otherwise.
 * @par Description:
 * Given a,b, and c, return nonzero if (a - b < c), zero toherwise.
 * */
#define float128_EGlpNumIsDiffLess(a, b, c) ({\
	float128 __lpnum__= float128_sub (a, b);\
	float128_lt(__lpnum__, c);\
})

/* ========================================================================= */
/** @brief test if the first number is bigger to the second number
 * @param a float128 the first number.
 * @param b double the second number
 * @return int one if success, zero otherwise.
 * @par Description:
 * Given two numbers 'a' and 'b', return one if a < b, zero
 * otherwise.
 * */
#define float128_EGlpNumIsLessDbl(a,b) ({\
	f64dbl_t __tmp = {.dbl = (b)};\
	float128 __lpnum__= float64_to_float128(__tmp.fl);\
	float128_lt(a,__lpnum__);})

/* ========================================================================= */
/** @brief test if the first number is bigger to the second number
 * @param a float128 the first number.
 * @param b double the second number
 * @return int one if success, zero otherwise.
 * @par Description:
 * Given two numbers 'a' and 'b', return one if a > b, zero
 * otherwise.
 * */
#define float128_EGlpNumIsGreaDbl(a,b) ({\
	f64dbl_t __tmp = {.dbl = (b)};\
	float128 __lpnum__= float64_to_float128(__tmp.fl);\
	float128_lt(__lpnum__,a);})

/* ========================================================================= */
/** @brief test if the first number is bigger to the second number
 * @param a float128 the first number.
 * @param b float128 the second number
 * @return int one if success, zero otherwise.
 * @par Description:
 * Given two numbers 'a' and 'b', return one if a <= b, zero
 * otherwise.
 * */
#define float128_EGlpNumIsLeq(a,b) float128_le(a,b)

/* ========================================================================= */
/** @brief copy the value of the second number to the first.
 * @param a float128 source number (it won't change value).
 * @param b float128 source number (it won't change value).
 * @param c float128 denominator of the difference (it won't change value).
 * @param d float128 where to store the value .
 * @par Description:
 * Set @f$a = \frac{b - c}{d} @f$ */
#define float128_EGlpNumCopyDiffRatio(a, b, c, d) ({\
	float128 __lpnum__= float128_sub (b, c);\
	a = float128_div (__lpnum__, d);\
})

/* ========================================================================= */
/** @brief copy the value of the second number to the first.
 * @param a float128 source number (it won't change value).
 * @param b float128 source number (it won't change value).
 * @param dest float128 where to store the value stored in 'orig'.
 * @par Description:
 * Set dest = a - b */
#define float128_EGlpNumCopyDiff(dest,a,b) (dest = float128_sub(a,b))

/* ========================================================================= */
/** @brief copy the value of the sum of the second and third parameter
 * @param a float128 source number (it won't change value).
 * @param b float128 source number (it won't change value).
 * @param dest float128 where to store the sum.
 * @par Description:
 * Set dest = a + b */
#define float128_EGlpNumCopySum(dest,a,b) (dest = float128_add(a,b))

/* ========================================================================= */
/** @brief copy the value of the second number to the first.
 * @param orig float128 source number (it won't change value).
 * @param dest float128 where to store the value stored in 'orig'.
 * @par Description:
 * Given two numbers copy the values in 'orig', into 'dest'.
 * */
#define float128_EGlpNumCopy(dest,orig) (dest=orig)

/* ========================================================================= */
/** @brief change the fist number to the maximum between itself and the 
 * absolute value of the second.
 * @param orig float128 source number (it won't change value).
 * @param dest float128 where to store the value stored in 'orig'.
 * @par Description:
 * implement dest = max(dest,abs(orig))
 * */
#define float128_EGlpNumSetToMaxAbs(dest, orig) ({\
	float128 __lpnum__= (orig);\
	__lpnum__.high &= 0x7fffffffffffffffLL;\
	if (float128_lt(dest, __lpnum__)) (dest) = __lpnum__;})

#define float128_EGlpNumSetToMinAbs(dest, orig) ({\
	float128 __lpnum__=(orig);\
	__lpnum__.high &= 0x7fffffffffffffffLL;\
	if (float128_lt( __lpnum__, dest)) (dest) = __lpnum__;})

/* ========================================================================= */
/** @brief copy the square of the second argument, divided by the third 
 * argument into the first argument.
 * @param dest float128 where to store the result
 * @param orig float128 second parameter
 * @param den float128 third parameter
 * @par Description:
 * compute dest = (orig*orig)/den
 * */
#define float128_EGlpNumCopySqrOver(dest, orig, den) ({\
	float128 __lpnum__= float128_mul (orig, orig);\
	(dest) = float128_div ( __lpnum__, den);\
})

/* ========================================================================= */
/** @brief copy the value of the absolute value of the second parameter to the 
 * first parameter.
 * @param orig float128 source number (it won't change value).
 * @param dest float128 where to store the absolute value stored
 * in 'orig'.
 * @par Description:
 * Given a number 'orig', copy its absolute value to 'dest'. i.e.
 * dest = |orig|
 * */
#define float128_EGlpNumCopyAbs(dest,orig) ({\
	(dest) = (orig);\
	(dest).high &= 0x7fffffffffffffffLL;})

/* ========================================================================= */
/** @brief copy minus the value of the second parameter to the 
 * first parameter.
 * @param orig float128 the source number (it won't change value).
 * @param dest float128 where to store minus the value stored
 * in 'orig'.
 * @par Description:
 * Given a number 'orig', copy minus the value to 'dest'. i.e.
 * dest = -orig
 * */
#define float128_EGlpNumCopyNeg(dest,orig) ({\
	(dest) = (orig);\
	(dest).high ^= 0x8000000000000000LL;})

/* ========================================================================= */
/** @brief Set des = op1/op2.
 * @param dest float128 where we will store the result.
 * @param op1 float128 numerator of the fraction (possibly non an integer)
 * @param op2 float128 denominator of the fraction (possibly non an integer)
 * @par Description:
 *  Set des = op1/op2
 * */
#define float128_EGlpNumCopyFrac(dest,op1,op2) ((dest) = float128_div(op1,op2))

/* ========================================================================= */
/** @brief copy the first 'size' values in the second array to the first array.
 * @param orig float128* pointer to the array from where we will copy the
 * values (it won't change value).
 * @param dest float128* pointer to where to store the first 'size' values 
 * stored in 'orig'.
 * @param size unsigned int specifying how many values of 'orig' will be copied
 * onto 'dest'
 * @par Description:
 * This function is provided to (possible) make fast copies of arrays of
 * numbers, the arrays should be of length at least 'size', and the resulting
 * copy is absolutely independent froom the original, any change in one vale of
 * one array won't change values on the other array.
 * */
#define float128_EGlpNumCopyArray(dest,orig,size) memcpy(dest,orig,sizeof(float128)*(size))

/* ========================================================================= */
/** @brief Sub to a given number the product of two numbers.
 * @param a float128 the number that we are going to Sub to.
 * @param b float128 value to be multiplyed.
 * @param c float128 value to be multiplyed.
 * @par Description:
 * This function implements a = a - b*c, and clearly don't change the value
 * stored in 'b' nor in 'c'.
 * */
#define float128_EGlpNumSubInnProdTo(a, b, c) ({\
	float128 __lpnum__= float128_mul(b, c);\
	(a) = float128_sub(a, __lpnum__);\
})

/* ========================================================================= */
/** @brief Add to a given number the product of two numbers.
 * @param a float128 the number that we are going to add to.
 * @param b float128 value to be multiplyed.
 * @param c float128 value to be multiplyed.
 * @par Description:
 * This function implements a = a + b*c, and clearly don't change the value
 * stored in 'b' nor in 'c'.
 * */
#define float128_EGlpNumAddInnProdTo(a, b, c) ({\
	float128 __lpnum__= float128_mul ( b, c);\
	(a) = float128_add (a, __lpnum__);\
})

/* ========================================================================= */
/** @brief Substract to a given number the value of the second number.
 * @param a float128 the number that we are going to substract to.
 * @param b unsigned int value to be substracted to 'a'.
 * @par Description:
 * This function implements a = a - b, and clearly don't change the value
 * stored in 'b'.
 * */
#define float128_EGlpNumSubUiTo(a,b) ({\
	float128 __lpnum__= int64_to_float128((long long)(b));\
	(a) = float128_sub(a,__lpnum__);})

/* ========================================================================= */
/** @brief Add to a given number the value of the second number.
 * @param a float128 the number that we are going to add to.
 * @param b unsigned int value to be added to 'a'.
 * @par Description:
 * This function implements a = a + b, and clearly don't change the value
 * stored in 'b'.
 * */
#define float128_EGlpNumAddUiTo(a,b) ({\
	float128 __lpnum__= int64_to_float128((long long)(b));\
	(a) = float128_add(a,__lpnum__);})

/* ========================================================================= */
/** @brief Add to a given number the value of the second number.
 * @param a float128 the number that we are going to add to.
 * @param b float128 value to be added to 'a'.
 * @par Description:
 * This function implements a = a + b, and clearly don't change the value
 * stored in 'b'.
 * */
#define float128_EGlpNumAddTo(a,b) ((a) = float128_add(a,b))

/* ========================================================================= */
/** @brief Substract to a given number the value of the second number.
 * @param a float128 the number that we are going to substract
 * from.
 * @param b float128 value to be substracted to 'a'.
 * @par Description:
 * This function implements a = a - b, and clearly don't change the value
 * stored in 'b'.
 * */
#define float128_EGlpNumSubTo(a,b) ((a) = float128_sub(a,b))

/* ========================================================================= */
/** @brief Multiply a given number by the value of the second number.
 * @param a float128 the number that we are going to multiply by
 * the second number and store the result.
 * @param b float128 value to be multyply to 'a'.
 * @par Description:
 * This function implements a = a * b, and clearly don't change the value
 * stored in 'b'.
 * */
#define float128_EGlpNumMultTo(a,b) ((a) = float128_mul(a,b))

/* ========================================================================= */
/** @brief Divide a given number by the value of the second number.
 * @param a float128 the number that we are going to divide by
 * the second number and store the result.
 * @param b float128 value to be divide to 'a'.
 * @par Description:
 * This function implements a = a / b, and clearly don't change the value
 * stored in 'b'.
 * */
#define float128_EGlpNumDivTo(a,b) ((a) = float128_div(a,b))

/* ========================================================================= */
/** @brief Divide a given number by the value of the second number.
 * @param a float128 the number that we are going to divide by
 * the second number and store the result.
 * @param b unsigned int value to be divided to 'a'.
 * @par Description:
 * This function implements a = a / b, and don't change the value
 * stored in 'b'.
 * */
#define float128_EGlpNumDivUiTo(a,b) ({\
	float128 __lpnum__= int64_to_float128((long long)(b));\
	(a) = float128_div(a,__lpnum__);})

/* ========================================================================= */
/** @brief Multiply a given number by the value of the second number.
 * @param a float128 the number that we are going to multiply by
 * the second number and store the result.
 * @param b unsigned int value to be multyply to 'a'.
 * @par Description:
 * This function implements a = a * b, and clearly don't change the value
 * stored in 'b'.
 * */
#define float128_EGlpNumMultUiTo(a,b) ({\
	float128 __lpnum__= int64_to_float128((long long)(b));\
	(a) = float128_mul(a,__lpnum__);})

/* ========================================================================= */
/** @brief Reset the value of the pointed number to zero.
 * @param a float128 the value to be set to zero.
 * @par Descrpition:
 * Reset a to zero, i.e. implements a = 0;
 * */
#define float128_EGlpNumZero(a) ((a) = (float128){0,0})

/* ========================================================================= */
/** @brief Reset the value of the pointed number to one.
 * @param a float128 value to be set to one.
 * @par Descrpition:
 * Reset a to zero, i.e. implements a = 1;
 * */
#define float128_EGlpNumOne(a) ((a) = (float128){.high = 0x3fff000000000000LL, .low = 0})

/* ========================================================================= */
/** @brief Change the sign of the number.
 * @param a float128 number we will change sign.
 * @par Descrpition:
 * Change the sign of the given number, i.e. implements a = -a
 * */
#define float128_EGlpNumSign(a) ((a).high ^= 0x8000000000000000LL)

/* ========================================================================= */
/** @brief return the closest double value of the given pointer number.
 * @param a float128 number that we will be transformed to double.
 * @return double the closest double representation of the given number.
 * par Description:
 * return the double number closest in value to the value stored in a.
 * */
#define float128_EGlpNumToLf(a) ({\
	const f64dbl_t __tmp = {.fl = float128_to_float64(a)};\
	__tmp.dbl;})

/* ========================================================================= */
/** @brief initialize the internal memory of a given variable */
#define float128_EGlpNumInitVar(a) (a = (float128){0,0})

/* ========================================================================= */
/** @brief free the internal memory of a given variable */
#define float128_EGlpNumClearVar(a)

/* ========================================================================= */
/** @} */
#endif
#endif
#endif
