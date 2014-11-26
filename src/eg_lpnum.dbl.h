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
#ifndef __EG_LPNUM_DBL__
#define __EG_LPNUM_DBL__

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "eg_lpnum.h"

/** @file
 * @ingroup EGlpNum */
/** @addtogroup EGlpNum */
/** @{ */
/* ========================================================================= */
/** extern definitions of constaants for different set-ups */
#define dbl_zeroLpNum 0.0
#define dbl_oneLpNum  1.0
#define dbl_epsLpNum  DBL_EPSILON
#define dbl_MaxLpNum	DBL_MAX
#define dbl_MinLpNum	DBL_MIN

/* ========================================================================= */
/** @brief Read from a string a number and store it in the given double, 
 * @return the number of chars readed from the input string */
#define dbl_EGlpNumReadStr(a,str) ({\
	int __i =0;\
	sscanf(str,"%lf%n",&(a),&__i);\
	__i;})

/* ========================================================================= */
/** @brief given a double, write it to a string (to be allocated internally), 
 * and return it. */
#define dbl_EGlpNumGetStr(a) ({\
	char *__str=0;\
	int __i=snprintf(__str,(size_t)0,"%.7lg",a);\
	__str = EGsMalloc(char,__i+1);\
	snprintf(__str,(size_t)(__i+1),"%.7lg",a);\
	__str;})

/* ========================================================================= */
/** @brief given an array of type double, free it, if the pointer is NULL
 * nothing happen. */
#define dbl_EGlpNumFreeArray(ea) __EGlpNumFreeArray(ea)

/* ========================================================================= */
/** @brief Reallocate and initialize (if needed) 'size' elements of type 
 * EGlpNum_t and return it, if no more memory, exit(1) */
#define dbl_EGlpNumReallocArray(lptr, lsize) ({\
	size_t __csz = (lsize), *__usp = 0;\
	size_t __psz = __EGlpNumArraySize(*lptr);\
	double** __ptr__ = (lptr);\
	if (!__psz) *__ptr__ = dbl_EGlpNumAllocArray (__csz); \
	else if (__psz < __csz) {\
		__usp = (size_t*)(*__ptr__);\
		__usp--;\
		__usp = EGrealloc(__usp, sizeof(double)*__csz+sizeof(size_t));\
		__usp[0] = __csz;\
		*__ptr__ = (double*)(__usp+1);\
		memset((*__ptr__)+__psz,0,sizeof(double)*(__csz-__psz));\
	}\
	*__ptr__;})

/* ========================================================================= */
/** @brief Allocate and initialize (if needed) 'size' elements of type double
 * and return it, if no more memory, exit(1) */
#define dbl_EGlpNumAllocArray(size) __EGlpNumAllocArray(double,size)

/* ========================================================================= */
/** @brief set the given number pointer, set its value to the given double.
 * @param var double where we will store the double value.
 * @param dbl double value to be stored in 'var'.
 * @par Description:
 * This function is intended to set initial values to variables; note that the
 * double is a number and not a pointer to that value, be carefull with this
 * detail. Also, due to implementation details this function can't deal with
 * numbers above 1e158 or smaller than 1e-158. Note also that if the number is
 * writen in the form \f$x=\bar{x}\cdot 2^e\f$ with \f$0.5<|\bar{x}|<1\f$, 
 * then \f$\left|x-\frac{p}{q}\right|<2^{e-64}\f$.
 * */
#define dbl_EGlpNumSet(var, dbl) ((var) = (dbl))

/* ========================================================================= */
/** @brief Stores in the first number the ceil value of the second number, i.e.
 * EGlpNumCeil(a,b) <==> a= ceil(b) */
#define dbl_EGlpNumCeil(a, b) ((a) = ceil(b))

/* ========================================================================= */
/** @brief Stores in the first number the floor value of the second number, i.e.
 * EGlpNumFloor(a,b) <==> a= floor(b) */
#define dbl_EGlpNumFloor(a, b) ((a) = floor(b))

/* ========================================================================= */
/** @brief store the (multiplicative) inverse of a number to itself, i.e.
 * implement a = 1/a.
 * @param a the number to be inverted. */
#define dbl_EGlpNumInv(a) ((a) = 1.0/(a))

/* ========================================================================= */
/** @brief Compare if two numbers are equal within a maximum error.
 * @param a EGlpNum_t first number to compare.
 * @param b EGlpNum_t second number to compare.
 * @return int one in success, zero oterwise.
 * @par Description:
 * Given two numbers 'a','b' return 1 if a == b, otherwise it return 0
 * */
#define dbl_EGlpNumIsEqqual(a,b) ((a) == (b))

/* ========================================================================= */
/** @brief Compare if two numbers are equal within a maximum error.
 * @param a EGlpNum_t first number to compare.
 * @param b EGlpNum_t second number to compare.
 * @param error EGlpNum_t maximum difference allowed between both
 * numbers.
 * @return int one in success, zero oterwise.
 * @par Description:
 * Given two numbers 'a','b' and a tolerance 'error',
 * return 1 if |a-b|<= error, otherwise it return 0.
 * */
#define dbl_EGlpNumIsEqual(a,b,error) (fabs((a)-(b)) <= (error))
#define dbl_EGlpNumIsNeq(a,b,error) (((a)-(b) > (error)) || ((b)-(a) > (error)))
#define dbl_EGlpNumIsNeqq(a,b)  ((a) != (b))
#define dbl_EGlpNumIsNeqZero(a,error) (((a) > (error)) || (-(a) > (error)))
#define dbl_EGlpNumIsNeqqZero(a)     	((a) != 0.0)

/* ========================================================================= */
/** @brief test if the first number is greater than zero
 * @param a EGlpNum_t number to test.
 * @return int one if success, zero otherwise.
 * */
#define dbl_EGlpNumIsGreatZero(a) (a > 0.0)

/* ========================================================================= */
/** @brief test if the first number is less than zero
 * @param a EGlpNum_t number to test.
 * @return int one if success, zero otherwise.
 * */
#define dbl_EGlpNumIsLessZero(a) (a < 0.0)

/* ========================================================================= */
/** @brief test if the first number is bigger to the second number
 * @param a EGlpNum_t the first number.
 * @param b EGlpNum_t the second number
 * @return int one if success, zero otherwise.
 * @par Description:
 * Given two numbers 'a' and 'b', return one if a < b, zero
 * otherwise.
 * */
#define dbl_EGlpNumIsLess(a,b) (a < b)

/* ========================================================================= */
/** @brief test if the sum of the first two numbers is less thatn the third
 * number.
 * @param a EGlpNum_t the first number.
 * @param b EGlpNum_t the second number
 * @param c EGlpNum_t the third number
 * @return int one if success, zero otherwise.
 * @par Description:
 * Given a,b, and c, return nonzero if (a + b < c), zero toherwise.
 * */
#define dbl_EGlpNumIsSumLess(a, b, c) ((a) + (b) < (c))

/* ========================================================================= */
/** @brief test if the diference of the first two numbers is less thatn the 
 * third number.
 * @param a EGlpNum_t the first number.
 * @param b EGlpNum_t the second number
 * @param c EGlpNum_t the third number
 * @return int one if success, zero otherwise.
 * @par Description:
 * Given a,b, and c, return nonzero if (a - b < c), zero toherwise.
 * */
#define dbl_EGlpNumIsDiffLess(a, b, c) ((a) - (b) < (c))

/* ========================================================================= */
/** @brief test if the first number is bigger to the second number
 * @param a EGlpNum_t the first number.
 * @param b double the second number
 * @return int one if success, zero otherwise.
 * @par Description:
 * Given two numbers 'a' and 'b', return one if a < b, zero
 * otherwise.
 * */
#define dbl_EGlpNumIsLessDbl(a,b) ((a) < (b))

/* ========================================================================= */
/** @brief test if the first number is bigger to the second number
 * @param a EGlpNum_t the first number.
 * @param b double the second number
 * @return int one if success, zero otherwise.
 * @par Description:
 * Given two numbers 'a' and 'b', return one if a > b, zero
 * otherwise.
 * */
#define dbl_EGlpNumIsGreaDbl(a,b) ((a) > (b))

/* ========================================================================= */
/** @brief test if the first number is bigger to the second number
 * @param a EGlpNum_t the first number.
 * @param b EGlpNum_t the second number
 * @return int one if success, zero otherwise.
 * @par Description:
 * Given two numbers 'a' and 'b', return one if a <= b, zero
 * otherwise.
 * */
#define dbl_EGlpNumIsLeq(a,b) ((a) <= (b))

/* ========================================================================= */
/** @brief copy the value of the second number to the first.
 * @param a EGlpNum_t source number (it won't change value).
 * @param b EGlpNum_t source number (it won't change value).
 * @param den EGlpNum_t denominator of the difference (it won't change value).
 * @param dest EGlpNum_t where to store the value .
 * @par Description:
 * Set dest = (a - b) / den */
#define dbl_EGlpNumCopyDiffRatio(dest,a, b, den) ((dest) = ((a) - (b)) / (den))

/* ========================================================================= */
/** @brief copy the value of the second number to the first.
 * @param a EGlpNum_t source number (it won't change value).
 * @param b EGlpNum_t source number (it won't change value).
 * @param dest EGlpNum_t where to store the value stored in 'orig'.
 * @par Description:
 * Set dest = a - b */
#define dbl_EGlpNumCopyDiff(dest,a,b) ((dest) = (a) - (b))

/* ========================================================================= */
/** @brief copy the value of the sum of the second and third parameter
 * @param a EGlpNum_t source number (it won't change value).
 * @param b EGlpNum_t source number (it won't change value).
 * @param dest EGlpNum_t where to store the sum.
 * @par Description:
 * Set dest = a + b */
#define dbl_EGlpNumCopySum(dest,a,b) ((dest) = (a) + (b))

/* ========================================================================= */
/** @brief copy the value of the second number to the first.
 * @param orig EGlpNum_t source number (it won't change value).
 * @param dest EGlpNum_t where to store the value stored in 'orig'.
 * @par Description:
 * Given two numbers copy the values in 'orig', into 'dest'.
 * */
#define dbl_EGlpNumCopy(dest,orig) ((dest) = (orig))

/* ========================================================================= */
/** @brief change the fist number to the maximum between itself and the 
 * absolute value of the second.
 * @param orig EGlpNum_t source number (it won't change value).
 * @param dest EGlpNum_t where to store the value stored in 'orig'.
 * @par Description:
 * implement dest = max(dest,abs(orig))
 * */
#define dbl_EGlpNumSetToMaxAbs(dest, orig) if((dest) < fabs(orig)) \
																							(dest) = fabs(orig)
#define dbl_EGlpNumSetToMinAbs(dest, orig) if((dest) > fabs(orig)) \
																							(dest) = fabs(orig)

/* ========================================================================= */
/** @brief copy the square of the second argument, divided by the third 
 * argument into the first argument.
 * @param dest EGlpNum_t where to store the result
 * @param orig EGlpNum_t second parameter
 * @param den EGlpNum_t third parameter
 * @par Description:
 * compute dest = (orig*orig)/den
 * */
#define dbl_EGlpNumCopySqrOver(dest, orig, den) ((dest) = (orig)*(orig)/(den))

/* ========================================================================= */
/** @brief copy the value of the absolute value of the second parameter to the 
 * first parameter.
 * @param orig EGlpNum_t source number (it won't change value).
 * @param dest EGlpNum_t where to store the absolute value stored
 * in 'orig'.
 * @par Description:
 * Given a number 'orig', copy its absolute value to 'dest'. i.e.
 * dest = |orig|
 * */
#define dbl_EGlpNumCopyAbs(dest,orig) ((dest) = fabs(orig))

/* ========================================================================= */
/** @brief copy minus the value of the second parameter to the 
 * first parameter.
 * @param orig EGlpNum_t the source number (it won't change value).
 * @param dest EGlpNum_t where to store minus the value stored
 * in 'orig'.
 * @par Description:
 * Given a number 'orig', copy minus the value to 'dest'. i.e.
 * dest = -orig
 * */
#define dbl_EGlpNumCopyNeg(dest,orig) ((dest) = -(orig))

/* ========================================================================= */
/** @brief Set des = op1/op2.
 * @param dest EGlpNum_t where we will store the result.
 * @param op1 EGlpNum_t numerator of the fraction (possibly non an integer)
 * @param op2 EGlpNum_t denominator of the fraction (possibly non an integer)
 * @par Description:
 *  Set des = op1/op2
 * */
#define dbl_EGlpNumCopyFrac(dest,op1,op2) ((dest) = (op1)/(op2))

/* ========================================================================= */
/** @brief copy the first 'size' values in the second array to the first array.
 * @param orig EGlpNum_t* pointer to the array from where we will copy the
 * values (it won't change value).
 * @param dest EGlpNum_t* pointer to where to store the first 'size' values 
 * stored in 'orig'.
 * @param size unsigned int specifying how many values of 'orig' will be copied
 * onto 'dest'
 * @par Description:
 * This function is provided to (possible) make fast copies of arrays of
 * numbers, the arrays should be of length at least 'size', and the resulting
 * copy is absolutely independent froom the original, any change in one vale of
 * one array won't change values on the other array.
 * */
#define dbl_EGlpNumCopyArray(dest,orig,size) memcpy(dest,orig,sizeof(double)*(size))

/* ========================================================================= */
/** @brief Sub to a given number the product of two numbers.
 * @param a EGlpNum_t the number that we are going to Sub to.
 * @param b EGlpNum_t value to be multiplyed.
 * @param c EGlpNum_t value to be multiplyed.
 * @par Description:
 * This function implements a = a - b*c, and clearly don't change the value
 * stored in 'b' nor in 'c'.
 * */
#define dbl_EGlpNumSubInnProdTo(a, b, c) ((a) -= (b)*(c))

/* ========================================================================= */
/** @brief Add to a given number the product of two numbers.
 * @param a EGlpNum_t the number that we are going to add to.
 * @param b EGlpNum_t value to be multiplyed.
 * @param c EGlpNum_t value to be multiplyed.
 * @par Description:
 * This function implements a = a + b*c, and clearly don't change the value
 * stored in 'b' nor in 'c'.
 * */
#define dbl_EGlpNumAddInnProdTo(a, b, c) ((a) += (b)*(c))

/* ========================================================================= */
/** @brief Substract to a given number the value of the second number.
 * @param a EGlpNum_t the number that we are going to substract to.
 * @param b unsigned int value to be substracted to 'a'.
 * @par Description:
 * This function implements a = a - b, and clearly don't change the value
 * stored in 'b'.
 * */
#define dbl_EGlpNumSubUiTo(a,b) ((a) -= (b))

/* ========================================================================= */
/** @brief Add to a given number the value of the second number.
 * @param a EGlpNum_t the number that we are going to add to.
 * @param b unsigned int value to be added to 'a'.
 * @par Description:
 * This function implements a = a + b, and clearly don't change the value
 * stored in 'b'.
 * */
#define dbl_EGlpNumAddUiTo(a,b) ((a) += (b))

/* ========================================================================= */
/** @brief Add to a given number the value of the second number.
 * @param a EGlpNum_t the number that we are going to add to.
 * @param b EGlpNum_t value to be added to 'a'.
 * @par Description:
 * This function implements a = a + b, and clearly don't change the value
 * stored in 'b'.
 * */
#define dbl_EGlpNumAddTo(a,b) ((a) += (b))

/* ========================================================================= */
/** @brief Substract to a given number the value of the second number.
 * @param a EGlpNum_t the number that we are going to substract
 * from.
 * @param b EGlpNum_t value to be substracted to 'a'.
 * @par Description:
 * This function implements a = a - b, and clearly don't change the value
 * stored in 'b'.
 * */
#define dbl_EGlpNumSubTo(a,b) ((a) -= (b))

/* ========================================================================= */
/** @brief Multiply a given number by the value of the second number.
 * @param a EGlpNum_t the number that we are going to multiply by
 * the second number and store the result.
 * @param b EGlpNum_t value to be multyply to 'a'.
 * @par Description:
 * This function implements a = a * b, and clearly don't change the value
 * stored in 'b'.
 * */
#define dbl_EGlpNumMultTo(a,b) ((a) *= (b))

/* ========================================================================= */
/** @brief Divide a given number by the value of the second number.
 * @param a EGlpNum_t the number that we are going to divide by
 * the second number and store the result.
 * @param b EGlpNum_t value to be divide to 'a'.
 * @par Description:
 * This function implements a = a / b, and clearly don't change the value
 * stored in 'b'.
 * */
#define dbl_EGlpNumDivTo(a,b) ((a) /= (b))

/* ========================================================================= */
/** @brief Divide a given number by the value of the second number.
 * @param a EGlpNum_t the number that we are going to divide by
 * the second number and store the result.
 * @param b unsigned int value to be divided to 'a'.
 * @par Description:
 * This function implements a = a / b, and don't change the value
 * stored in 'b'.
 * */
#define dbl_EGlpNumDivUiTo(a,b) ((a) /= (b))

/* ========================================================================= */
/** @brief Multiply a given number by the value of the second number.
 * @param a EGlpNum_t the number that we are going to multiply by
 * the second number and store the result.
 * @param b unsigned int value to be multyply to 'a'.
 * @par Description:
 * This function implements a = a * b, and clearly don't change the value
 * stored in 'b'.
 * */
#define dbl_EGlpNumMultUiTo(a,b) ((a) *= (b))

/* ========================================================================= */
/** @brief Reset the value of the pointed number to zero.
 * @param a EGlpNum_t the value to be set to zero.
 * @par Descrpition:
 * Reset a to zero, i.e. implements a = 0;
 * */
#define dbl_EGlpNumZero(a) ((a) = 0.0)

/* ========================================================================= */
/** @brief Reset the value of the pointed number to one.
 * @param a EGlpNum_t value to be set to one.
 * @par Descrpition:
 * Reset a to zero, i.e. implements a = 1;
 * */
#define dbl_EGlpNumOne(a) ((a) = 1.0)

/* ========================================================================= */
/** @brief Change the sign of the number.
 * @param a EGlpNum_t number we will change sign.
 * @par Descrpition:
 * Change the sign of the given number, i.e. implements a = -a
 * */
#define dbl_EGlpNumSign(a) ((a) = -(a))

/* ========================================================================= */
/** @brief return the closest double value of the given pointer number.
 * @param a EGlpNum_t number that we will be transformed to double.
 * @return double the closest double representation of the given number.
 * par Description:
 * return the double number closest in value to the value stored in a.
 * */
#define dbl_EGlpNumToLf(a) ((double)a)

/* ========================================================================= */
/** @brief initialize the internal memory of a given variable */
#define dbl_EGlpNumInitVar(a) dbl_EGlpNumZero(a)

/* ========================================================================= */
/** @brief free the internal memory of a given variable */
#define dbl_EGlpNumClearVar(a)

/* ========================================================================= */
/** @} */
#endif
