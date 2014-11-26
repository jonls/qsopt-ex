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
#ifndef __EG_LPNUM_H__
#define __EG_LPNUM_H__
/* ========================================================================= */
/** @defgroup EGlpNum EGlpNum
 *
 * Here we define a common interface to handle numbers in general, the idea is
 * to be able to work with infinite precicion numbers, plain doubles, floats,
 * integers, or fixed point numbers, without actually making different codes
 * for each of those types, rather we preffer to fix that at compyle time.
 *
 * @par History:
 * Revision 1.1.0
 * 	- 2013-04-18
 * 						- Add uint32_EGutilPermSort and uint32_EGutilPermSort2
 * 	- 2010-09-02
 * 						- Add support for int32 templates, streamline the template
 * 						implementation too
 * 						- Add MaxLpNum and MinLpNum keywords
 *  We start doing the migration to gmp and 'natural' types, this means that we
 *  drop support for EGrat, this allow us to drop the requirement to use
 *  pointers, and instead we can just call the functions with the original
 *  parameters, still we have to be VERY carefull regarding changing
 *  local/external copies.
 *  - 2007-10-08
 *  					- Move EGswap, EGabs, EGmin and Egmax to eg_numutil.h
 *  - 2005-10-31
 *  					- Add EGswap to swap elements of any predefined type.
 *  - 2005-08-31
 *  					- Add EGmin and EGmax for built in types (i.e. for types where
 *  					the < comparison works as we want).
 *  - 2005-08-16
 *  					- Streamline mpq_EGlpNumGetStr
 *  					- Minor Fixes for zeroLpNum
 *  - 2005-07-29
 *  					- Add EGabs definition.
 *  - 2005-07-24
 *  					- Split eg_lpnum.h into different headers for each type of
 *  						suported numbers.
 *  					- Deprecate EGlpNumCOmpUFrac
 *  - 2005-05-26
 *  					- Add epsLpNum
 *  - 2005-05-17
 *  					- Add mpq_EGlpNumReadStrXc(mpq_t,__EGPstr__)
 *  - 2005-05-16
 *  					- Add mpq_EGlpNumSet_mpf(mpq,mpf)
 *  - 2005-05-12
 *  					- Add mpf_EGlpNumEpow(num,power)
 *  					- Add function to change precision of the numbers on the fly.
 *  					- Add EGlpNumReadStr to set a number from a given input string.
 *  					- Add EGlpNumGetStr to get (hopefully) the exact representation 
 *  						of the given input as string.
 *  - 2005-05-03
 *  					- Change the structure of the header so it provides an interface
 *  					in which a program can use all types of numbers at the same time,
 *  					this implies that we must define a start-up and clean-up function
 *  					that would initialize all constants for all numbers, and change
 *  					the naming scheme accordingly to support this.
 *  					- Deprecate EGlpNumInitArray
 *  					- Deprecate EGlpNumFreeIntArray
 *  					- Change all static-inline definitions to Statement Exprs style.
 *  - 2005-04-25
 *  					- Add EGlpNumIsSumLess(__EGPa__,__EGPb__,__EGPc__)
 *  					- Add EGlpNumIsDiffLess(__EGPa__,__EGPb__,__EGPc__)
 *  - 2005-04-13
 *  					- Add EGlpNumCopyDiffRatio(__EGPa__,__EGPb__,__EGPc__,d)
 *  					- Add EGlpNumIsEqqual(__EGPa__,__EGPb__)
 *  					- Add EGlpNumIsEqual(__EGPa__,__EGPb__,__EGPerr__)
 *  					- Add EGlpNumCopy(__EGPa__,__EGPb__)
 *  					- Add EGlpNumIsLess(__EGPa__,__EGPb__)
 *  					- Add EGlpNumToLf(__EGPa__)
 *  					- Add EGlpNumCopyDiff(__EGPa__,__EGPb__,__EGPc__)
 *  					- Add EGlpNumCopyAbs(__EGPa__,__EGPb__)
 *  					- Add EGlpNumSubTo(__EGPa__,__EGPb__)
 *  					- Add EGlpNumAddTo(__EGPa__,__EGPb__)
 *  					- Add EGlpNumDivTo(__EGPa__,__EGPb__)
 *  					- Add EGlpNumMultTo(__EGPa__,__EGPb__)
 *  					- Add EGlpNumZero(__EGPa__)
 *  					- Add EGlpNumOne(__EGPa__)
 *  					- Add EGlpNumAddInnProdTo(__EGPa__,__EGPb__,__EGPc__)
 *  					- Add EGlpNumSubInnProdTo(__EGPa__,__EGPb__,__EGPc__)
 *  					- Add EGlpNumSign(__EGPa__)
 *  					- Add EGlpNumCopyNeg(__EGPa__,__EGPb__)
 *  					- Add EGlpNumDivUiTo(__EGPa__,__EGPb__)
 *  					- Add EGlpNumMultUiTo(__EGPa__,__EGPb__)
 *  					- Add EGlpNumIsLeq(__EGPa__,__EGPb__)
 *  					- Add EGlpNumCopySqrOver(__EGPa__,__EGPb__,__EGPc__)
 *  					- Add EGlpNumSet(__EGPa__,__EGPb__)
 *  					- Add EGlpNumCopyFrac(__EGPa__,__EGPb__,__EGPc__)
 *  					- Add EGlpNumAddUiTo(__EGPa__,__EGPb__)
 *  					- Add EGlpNumSubUiTo(__EGPa__,__EGPb__)
 *  					- Add EGlpNumCopySum(__EGPa__,__EGPb__,__EGPc__)
 *  					- Add EGlpNumInv(__EGPa__)
 *  					- Add EGlpNumFloor(__EGPa__)
 *  					- Add EGlpNumCeil(__EGPa__)
 *  					- Add EGlpNumIsLessDbl(__EGPa__,__EGPb__)
 *  					- Add EGlpNumIsGreaDbl(__EGPa__,__EGPb__)
 *  					- Add EGlpNumSetToMaxAbs(__EGPa__,__EGPb__)
 *  					- Add EGlpNumAllocArray(__EGPszb__)
 *  					- Add EGlpNumFreeArray(__EGPa__)
 *  					- Add EGlpNumReallocArray(__EGPa__,__EGPszb__)
 *  					- Add EGlpNumInitVar(__EGPa__)
 *  					- Add EGlpNumClearVar(__EGPa__)
 *  					- Add EGlpNumIsNeqq(__EGPa__,__EGPb__)
 *  					- Add EGlpNumIsNeq(__EGPa__,__EGPb__,__EGPc__)
 *  					- Add EGlpNumIsNeqqZero(__EGPa__)
 *  					- Add EGlpNumIsNeqZero(__EGPa__,__EGPb__)
 *  					- Add EGlpNumSetToMinAbs(__EGPa__,__EGPb__)
 * Revision 0.0.1
 * - 2004-07-15
 * 						- Add support for GNU_MP_F types
 * - 2004-07-14
 * 						- Add support for GNU_MP_Q types
 * - 2004-07-12
 * 						- Add support for EG-rationals
 * - 2004-06-21
 * 						- First Implementation/Definition
 * */

/** @file
 * @ingroup EGlpNum */
/** @addtogroup EGlpNum */
/** @{ */
/* ========================================================================= */

#include <stdlib.h>
#include <stdint.h>

#include <gmp.h>

/* ========================================================================= */
/** @name Number Types Definitions:
 * Define (as its name suggest) an internal identifier for the given 
 * type. this definitions are provided to select different types of data at 
 * compile  time, thus allowing us to provide limited template support. */
/* @{ */
/** C double type. */
#define DBL_TYPE 0
/** C float type. */
#define FLT_TYPE 1
/** C int type. */
#define INT_TYPE 2
/** EGlib #EGfp10_t type, this is an implementation of fixed precision
 * arithmetic with 10 bits for fractional representation. */
#define FP10_TYPE 3
/** EGlib #EGfp20_t type, this is an implementation of fixed precision
 * arithmetic with 20 bits for fractional representation. */
#define FP20_TYPE 4
/** EGlib #EGfp28_t type, this is an implementation of fixed precision
 * arithmetic with 28 bits for fractional representation. */
#define FP28_TYPE 5
/** EGlib #EGfp25_t type, this is an implementation of fixed precision
 * arithmetic with 25 bits for fractional representation. */
#define FP25_TYPE 6
/** GNU_MP library mpz_t type */
#define GNU_MP_Z 8
/** GNU_MP library mpq_t type */
#define GNU_MP_Q 9
/** GNU_MP library mpf_t type */
#define GNU_MP_F 10
/** C long double type */
#define LDBL_TYPE 11
/** C long long int type */
#define LLINT_TYPE 12
/** SoftFloat 128-bit floating point numbner */
#define FLOAT128_TYPE 13
/** int32_t integers */
#define INT32_TYPE 14
/* @} */
/* ========================================================================= */

/* We have no mpz specific header so we include these here */
extern const mpz_t __zeroLpNum_mpz__;
extern const mpz_t __oneLpNum_mpz__;
extern const mpz_t __MaxLpNum_mpz__;
extern const mpz_t __MinLpNum_mpz__;
#define mpz_zeroLpNum __zeroLpNum_mpz__
#define mpz_oneLpNum  __oneLpNum_mpz__
#define mpz_epsLpNum  __zeroLpNum_mpz__
#define mpz_MaxLpNum  __MaxLpNum_mpz__
#define mpz_MinLpNum  __MinLpNum_mpz__

#include "eg_lpnum.dbl.h"
#include "eg_lpnum.mpq.h"
#include "eg_lpnum.mpf.h"
#include "eg_macros.h"
#include "eg_mem.h"
#include "eg_nummacros.h"
/* ========================================================================= */
/** @brief Debugging verbosity messages deped on the value of DEBUG (defined in
 * eg_configure.h) and on the value of EGLPNUM_DEBUGL macro defined here.
* */
#define EGLPNUM_DEBUGL 100

/* ========================================================================= */
/** @brief Set the default number of __BITS__ used in the precision of the
 * float point numbers (mpf_t), a normal double use up to 56-64 bits., the 
 * default precision is set to 128 */
extern unsigned long int EGLPNUM_PRECISION;

/* ========================================================================= */
/** @brief Change the default precision for mpf_t numbers. */
void EGlpNumSetPrecision (const unsigned prec);

/* ========================================================================= */
/** @brief Allocate an array of a given type and store (sizeof(size_t) bytes 
 * before the actual array) the size of the allocated array. 
 * @param __type the type of the array to be returned.
 * @param __size the length of the array to be returned, note that it can be
 * zero, in wich case no memory allocation is made and NULL is returned. */
#define __EGlpNumAllocArray(__type,__size) ({\
	size_t __sz = (__size);\
	size_t *__utmp = __sz ? (size_t*) EGmalloc (sizeof(__type) * __sz + sizeof(size_t)) : 0;\
	if(__sz) __utmp[0] = __sz;\
	(__type*)(__sz ? (__utmp+1):0);})

/* ========================================================================= */
/** @brief Given an array allocated with __EGlpNumAllocArray, return the size of
 * the given array, if the array is null, return zero. 
 * @param __array the array from where we need the size. */
/* ========================================================================= */
#define __EGlpNumArraySize(__array) ({\
	size_t *__utmp = (size_t*)(__array);\
	if(__utmp) __utmp--;\
	__utmp ? __utmp[0]:0;})

/* ========================================================================= */
/** @brief, given an array allocated by __EGlpNumAllocArray, free the allocated
 * memory.
 * @param __array the array to be freed, it can be null. The given array will
 * always pooint to NULL when this function is done.
 * */
/* ========================================================================= */
#define __EGlpNumFreeArray(__array) ({\
	size_t *__utmp = (size_t*)(__array);\
	if(__utmp) free (__utmp-1);\
	(__array) = 0;})

/* ========================================================================= */
/** @brief indicate if the global data needed for EGlpNum has been initialized,
 * if zero, initialization routine should be called. This is provided to allow
 * syncronization between libraries */
extern int __EGlpNum_setup;
/* ========================================================================= */
/** @brief initialization routine for global data. This function is called as a
 * constructor, but calling it twice won't cause any problems, it is provided
 * to ensure that all EGlpnum globals are initialized at the beggining 
 * and in case they where not (__EGlpNMum_setup==0), then call the initializator */
extern void EGlpNumStart(void);
/* ========================================================================= */
/** @brief This function must be called at the end of the program to free all
 * internal data used in the EGlpNum_t structures, once this function is called
 * any operation on EGlpNum_t types may fail.
 * */
extern void EGlpNumClear(void);
/* ========================================================================= */
/** @brief provided for backwards compatibility */
#define EGlpNumExit EGlpNumClear
/** @}*/
/* ========================================================================= */
/** @name Common EGlpNum Interface functions:
 * Here we define the basic functions needed to declare when implementing a
 * number template. */
/** @{ */
#ifdef UNUSED_INTERFACE
/* ========================================================================= */
/** extern definitions of constants for different set-ups */
#define zeroLpNum 
#define oneLpNum 
#define epsLpNum
#define MaxLpNum 
#define MinLpNum
/* ========================================================================= */
/** @brief Read from a string a number and store it in the given number, 
 * @return the number of chars readed from the input string */
#define EGlpNumReadStr(__EGPa__,__EGPstr__) 
/* ========================================================================= */
/** @brief given a int, write it to a string (to be allocated internally), 
 * and return it. */
#define EGlpNumGetStr(__EGPa__) 
/* ========================================================================= */
/** @brief given an array of type int, free it, if the pointer is NULL
 * nothing happen. */
#define EGlpNumFreeArray(__EGParr__)
/* ========================================================================= */
/** @brief Reallocate and initialize (if needed) 'size' elements of type 
 * EGlpNum_t and return it, if no more memory, exit(1) */
#define EGlpNumReallocArray(__EGPptr__,__EGPsza__) 
/* ========================================================================= */
/** @brief Allocate and initialize (if needed) 'size' elements of type int
 * and return it, if no more memory, exit(1) */
#define EGlpNumAllocArray(__EGPszb__)
/* ========================================================================= */
/** @brief set the given number, set its value to the given double.
 * @param var where we will store the int value.
 * @param dbl_var value to be stored in 'var'.
 * @par Description:
 * This function is intended to set initial values to variables; note that the
 * int is a number and not a pointer to that value, be carefull with this
 * detail. */
#define EGlpNumSet(__EGPnum__,__EGPdnum__)
/* ========================================================================= */
/** @brief Stores in the first number the ceil value of the second number, i.e.
 * EGlpNumCeil(__EGPa__,__EGPb__) <==> a= ceil(__EGPb__) */
#define EGlpNumCeil(__EGPa__,__EGPb__)
/* ========================================================================= */
/** @brief Stores in the first number the floor value of the second number, i.e.
 * EGlpNumFloor(__EGPa__,__EGPb__) <==> a= floor(__EGPb__) */
#define EGlpNumFloor(__EGPa__,__EGPb__)
/* ========================================================================= */
/** @brief store the (multiplicative) inverse of a number to itself, i.e.
 * implement a = 1/a.
 * @param a the number to be inverted. */
#define EGlpNumInv(__EGPa__)
/* ========================================================================= */
/** @brief Compare if two numbers are equal within a maximum error.
 * @param a EGlpNum_t first number to compare.
 * @param b EGlpNum_t second number to compare.
 * @return int one in success, zero oterwise.
 * @par Description:
 * Given two numbers 'a','b' return 1 if a == b, otherwise it return 0
 * */
#define EGlpNumIsEqqual(__EGPa__,__EGPb__)
/* ========================================================================= */
/** @brief Compare if two numbers are equal within a maximum error.
 * @param a EGlpNum_t first number to compare.
 * @param b EGlpNum_t second number to compare.
 * @param error EGlpNum_t maximum difference allowed between both
 * numbers.
 * @return int one in success, zero oterwise.
 * @par Description:
 * Given two numbers 'a','b' and a tolerance 'error',
 * return 1 if |a-b| <= error, otherwise it return 0.
 * */
#define EGlpNumIsEqual(__EGPa__,__EGPb__,__EGPerr__)
/* return 1 if |a-b| > error */
#define EGlpNumIsNeq(__EGPa__,__EGPb__,__EGPerr__)
/* return 1 if a == b */
#define EGlpNumIsNeqq(__EGPa__,__EGPb__)
/* return 1 if |a| <= error */
#define EGlpNumIsNeqZero(__EGPa__,__EGPerr__)
/* return 1 if a!= 0 */
#define EGlpNumIsNeqqZero(__EGPa__)
/* ========================================================================= */
/** @brief test if the first number is bigger to the second number
 * @param a EGlpNum_t the first number.
 * @param b EGlpNum_t the second number
 * @return int one if success, zero otherwise.
 * @par Description:
 * Given two numbers 'a' and 'b', return one if a < b, zero
 * otherwise.
 * */
#define EGlpNumIsLess(__EGPa__,__EGPb__)
/* ========================================================================= */
/** @brief test if a numer is greater than zero
 * @param a number to test
 * @return int one if success, zero otherwise.
 * */
#define EGlpNumIsGreatZero(__EGPa__)

/* ========================================================================= */
/** @brief test if a numer is less than zero
 * @param a number to test
 * @return int one if success, zero otherwise.
 * */
#define EGlpNumIsLessZero(__EGPa__)
/* ========================================================================= */
/** @brief test if the sum of the first two numbers is less thatn the third
 * number.
 * @param a EGlpNum_t the first number.
 * @param b EGlpNum_t the second number
 * @param c EGlpNum_t the third number
 * @return int one if success, zero otherwise.
 * @par Description:
 * Given a,__EGPb__, and c, return nonzero if (a + b < c), zero toherwise.
 * */
#define EGlpNumIsSumLess(__EGPa__,__EGPb__,__EGPc__)
/* ========================================================================= */
/** @brief test if the diference of the first two numbers is less thatn the 
 * third number.
 * @param a EGlpNum_t the first number.
 * @param b EGlpNum_t the second number
 * @param c EGlpNum_t the third number
 * @return int one if success, zero otherwise.
 * @par Description:
 * Given a,__EGPb__, and c, return nonzero if (a - b < c), zero toherwise.
 * */
#define EGlpNumIsDiffLess(__EGPa__,__EGPb__,__EGPc__)
/* ========================================================================= */
/** @brief test if the first number is bigger to the second number
 * @param a EGlpNum_t the first number.
 * @param b int the second number
 * @return int one if success, zero otherwise.
 * @par Description:
 * Given two numbers 'a' and 'b', return one if a < b, zero
 * otherwise.
 * */
#define EGlpNumIsLessDbl(__EGPa__,__EGPb__)
/* ========================================================================= */
/** @brief test if the first number is bigger to the second number
 * @param a EGlpNum_t the first number.
 * @param b int the second number
 * @return int one if success, zero otherwise.
 * @par Description:
 * Given two numbers 'a' and 'b', return one if a > b, zero
 * otherwise.
 * */
#define EGlpNumIsGreaDbl(__EGPa__,__EGPb__)
/* ========================================================================= */
/** @brief test if the first number is bigger to the second number
 * @param a EGlpNum_t the first number.
 * @param b EGlpNum_t the second number
 * @return int one if success, zero otherwise.
 * @par Description:
 * Given two numbers 'a' and 'b', return one if a <= b, zero
 * otherwise.
 * */
#define EGlpNumIsLeq(__EGPa__,__EGPb__)
/* ========================================================================= */
/** @brief copy the value of the second number to the first.
 * @param a EGlpNum_t source number (it won't change value).
 * @param b EGlpNum_t source number (it won't change value).
 * @param __EGden__ EGlpNum_t denominator of the difference (it won't change value).
 * @param __EGdest__ EGlpNum_t where to store the value .
 * @par Description:
 * Set __EGdest__ = (a - b) / __EGden__ */
#define EGlpNumCopyDiffRatio(__EGdest__,__EGPa__,__EGPb__,__EGden__)
/* ========================================================================= */
/** @brief copy the value of the second number to the first.
 * @param a EGlpNum_t source number (it won't change value).
 * @param b EGlpNum_t source number (it won't change value).
 * @param __EGdest__ EGlpNum_t where to store the value stored in '__EGorig__'.
 * @par Description:
 * Set __EGdest__ = a - b */
#define EGlpNumCopyDiff(__EGdest__,__EGPa__,__EGPb__)
/* ========================================================================= */
/** @brief copy the value of the sum of the second and third parameter
 * @param a EGlpNum_t source number (it won't change value).
 * @param b EGlpNum_t source number (it won't change value).
 * @param __EGdest__ EGlpNum_t where to store the sum.
 * @par Description:
 * Set __EGdest__ = a + b */
#define EGlpNumCopySum(__EGdest__,__EGPa__,__EGPb__)
/* ========================================================================= */
/** @brief copy the value of the second number to the first.
 * @param __EGorig__ EGlpNum_t source number (it won't change value).
 * @param __EGdest__ EGlpNum_t where to store the value stored in '__EGorig__'.
 * @par Description:
 * Given two numbers copy the values in '__EGorig__', into '__EGdest__'.
 * */
#define EGlpNumCopy(__EGdest__,__EGorig__)
/* ========================================================================= */
/** @brief change the fist number to the maximum between itself and the 
 * absolute value of the second.
 * @param __EGorig__ EGlpNum_t source number (it won't change value).
 * @param __EGdest__ EGlpNum_t where to store the value stored in '__EGorig__'.
 * @par Description:
 * implement __EGdest__ = max(__EGdest__,abs(__EGorig__))
 * */
#define EGlpNumSetToMaxAbs(__EGdest__,__EGorig__) 
#define EGlpNumSetToMinAbs(__EGdest__,__EGorig__)
/* ========================================================================= */
/** @brief copy the square of the second argument, divided by the third 
 * argument into the first argument.
 * @param __EGdest__ EGlpNum_t where to store the result
 * @param __EGorig__ EGlpNum_t second parameter
 * @param __EGden__ EGlpNum_t third parameter
 * @par Description:
 * compute __EGdest__ = (__EGorig__*__EGorig__)/__EGden__
 * */
#define EGlpNumCopySqrOver(__EGdest__,__EGorig__,__EGden__)
/* ========================================================================= */
/** @brief copy the value of the absolute value of the second parameter to the 
 * first parameter.
 * @param __EGorig__ EGlpNum_t source number (it won't change value).
 * @param __EGdest__ EGlpNum_t where to store the absolute value stored
 * in '__EGorig__'.
 * @par Description:
 * Given a number '__EGorig__', copy its absolute value to '__EGdest__'. i.e.
 * __EGdest__ = |__EGorig__|
 * */
#define EGlpNumCopyAbs(__EGdest__,__EGorig__)
/* ========================================================================= */
/** @brief copy minus the value of the second parameter to the 
 * first parameter.
 * @param __EGorig__ EGlpNum_t the source number (it won't change value).
 * @param __EGdest__ EGlpNum_t where to store minus the value stored
 * in '__EGorig__'.
 * @par Description:
 * Given a number '__EGorig__', copy minus the value to '__EGdest__'. i.e.
 * __EGdest__ = -__EGorig__
 * */
#define EGlpNumCopyNeg(__EGdest__,__EGorig__)
/* ========================================================================= */
/** @brief Set des = __EGPa__/__EGPb__.
 * @param __EGdest__ EGlpNum_t where we will store the result.
 * @param __EGPa__ EGlpNum_t numerator of the fraction (possibly non an integer)
 * @param __EGPb__ EGlpNum_t denominator of the fraction (possibly non an integer)
 * @par Description:
 *  Set __EGdest__ = __EGPa__/__EGPb__
 * */
#define EGlpNumCopyFrac(__EGdest__,__EGPa__,__EGPb__)
/* ========================================================================= */
/** @brief copy the first 'size' values in the second array to the first array.
 * @param __EGorig__ EGlpNum_t* pointer to the array from where we will copy the
 * values (it won't change value).
 * @param __EGdest__ EGlpNum_t* pointer to where to store the first 'size' values 
 * stored in '__EGorig__'.
 * @param size unsigned int specifying how many values of '__EGorig__' will be copied
 * onto '__EGdest__'
 * @par Description:
 * This function is provided to (possible) make fast copies of arrays of
 * numbers, the arrays should be of length at least 'size', and the resulting
 * copy is absolutely independent froom the original, any change in one vale of
 * one array won't change values on the other array.
 * */
#define EGlpNumCopyArray(__EGdest__,__EGorig__,__EGPszb__)
/* ========================================================================= */
/** @brief Sub to a given number the product of two numbers.
 * @param a EGlpNum_t the number that we are going to Sub to.
 * @param b EGlpNum_t value to be multiplyed.
 * @param c EGlpNum_t value to be multiplyed.
 * @par Description:
 * This function implements a = a - b*c, and clearly don't change the value
 * stored in 'b' nor in 'c'.
 * */
#define EGlpNumSubInnProdTo(__EGPa__,__EGPb__,__EGPc__) 
/* ========================================================================= */
/** @brief Add to a given number the product of two numbers.
 * @param a EGlpNum_t the number that we are going to add to.
 * @param b EGlpNum_t value to be multiplyed.
 * @param c EGlpNum_t value to be multiplyed.
 * @par Description:
 * This function implements a = a + b*c, and clearly don't change the value
 * stored in 'b' nor in 'c'.
 * */
#define EGlpNumAddInnProdTo(__EGPa__,__EGPb__,__EGPc__)
/* ========================================================================= */
/** @brief Substract to a given number the value of the second number.
 * @param a EGlpNum_t the number that we are going to substract to.
 * @param b unsigned int value to be substracted to 'a'.
 * @par Description:
 * This function implements a = a - b, and clearly don't change the value
 * stored in 'b'.
 * */
#define EGlpNumSubUiTo(__EGPa__,__EGPb__)
/* ========================================================================= */
/** @brief Add to a given number the value of the second number.
 * @param a EGlpNum_t the number that we are going to add to.
 * @param b unsigned int value to be added to 'a'.
 * @par Description:
 * This function implements a = a + b, and clearly don't change the value
 * stored in 'b'.
 * */
#define EGlpNumAddUiTo(__EGPa__,__EGPb__)
/* ========================================================================= */
/** @brief Add to a given number the value of the second number.
 * @param a EGlpNum_t the number that we are going to add to.
 * @param b EGlpNum_t value to be added to 'a'.
 * @par Description:
 * This function implements a = a + b, and clearly don't change the value
 * stored in 'b'.
 * */
#define EGlpNumAddTo(__EGPa__,__EGPb__)
/* ========================================================================= */
/** @brief Substract to a given number the value of the second number.
 * @param a EGlpNum_t the number that we are going to substract
 * from.
 * @param b EGlpNum_t value to be substracted to 'a'.
 * @par Description:
 * This function implements a = a - b, and clearly don't change the value
 * stored in 'b'.
 * */
#define EGlpNumSubTo(__EGPa__,__EGPb__)
/* ========================================================================= */
/** @brief Multiply a given number by the value of the second number.
 * @param a EGlpNum_t the number that we are going to multiply by
 * the second number and store the result.
 * @param b EGlpNum_t value to be multyply to 'a'.
 * @par Description:
 * This function implements a = a * b, and clearly don't change the value
 * stored in 'b'.
 * */
#define EGlpNumMultTo(__EGPa__,__EGPb__)
/* ========================================================================= */
/** @brief Divide a given number by the value of the second number.
 * @param a EGlpNum_t the number that we are going to divide by
 * the second number and store the result.
 * @param b EGlpNum_t value to be divide to 'a'.
 * @par Description:
 * This function implements a = a / b, and clearly don't change the value
 * stored in 'b'.
 * */
#define EGlpNumDivTo(__EGPa__,__EGPb__)
/* ========================================================================= */
/** @brief Divide a given number by the value of the second number.
 * @param a EGlpNum_t the number that we are going to divide by
 * the second number and store the result.
 * @param b unsigned int value to be divided to 'a'.
 * @par Description:
 * This function implements a = a / b, and don't change the value
 * stored in 'b'.
 * */
#define EGlpNumDivUiTo(__EGPa__,__EGPb__)
/* ========================================================================= */
/** @brief Multiply a given number by the value of the second number.
 * @param a EGlpNum_t the number that we are going to multiply by
 * the second number and store the result.
 * @param b unsigned int value to be multyply to 'a'.
 * @par Description:
 * This function implements a = a * b, and clearly don't change the value
 * stored in 'b'.
 * */
#define EGlpNumMultUiTo(__EGPa__,__EGPb__)
/* ========================================================================= */
/** @brief Reset the value of the pointed number to zero.
 * @param a EGlpNum_t the value to be set to zero.
 * @par Descrpition:
 * Reset a to zero, i.e. implements a = 0;
 * */
#define EGlpNumZero(__EGPa__)
/* ========================================================================= */
/** @brief Reset the value of the pointed number to one.
 * @param a EGlpNum_t value to be set to one.
 * @par Descrpition:
 * Reset a to zero, i.e. implements a = 1;
 * */
#define EGlpNumOne(__EGPa__)
/* ========================================================================= */
/** @brief Change the sign of the number.
 * @param a EGlpNum_t number we will change sign.
 * @par Descrpition:
 * Change the sign of the given number, i.e. implements a = -a
 * */
#define EGlpNumSign(__EGPa__)
/* ========================================================================= */
/** @brief return the closest double value of the given pointer number.
 * @param a EGlpNum_t number that we will be transformed to int.
 * @return int the closest int representation of the given number.
 * par Description:
 * return the int number closest in value to the value stored in a.
 * */
#define EGlpNumToLf(__EGPa__)
/* ========================================================================= */
/** @brief initialize the internal memory of a given variable and set its
 * initial value to zero */
#define EGlpNumInitVar(__EGPa__)

/* ========================================================================= */
/** @brief free the internal memory of a given variable */
#define EGlpNumClearVar(__EGPa__)
#endif
 /** @} */
/* ========================================================================= */
/** @brief Sort (in increasing order) a sub-set of entries in an array using 
 * quicksort, by permutating the order of the elements in the subset rather 
 * than in the whole original array.
 * @param sz length of the permutation array.
 * @param perm array of indices of elements that we want to sort.
 * @param elem array (of length at least max(perm[k]:k=0,...,sz-1)) containing
 * the elements to be sorted.
 * @note The array of elements is not changed by this function.
 * @note This code is based in concorde's implementation of
 * permutation-quick-sort.
 * */
void uint32_EGutilPermSort (const size_t sz,
										 int *const perm,
										 const uint32_t * const elem);

/* ========================================================================= */
/** @brief Sort (in decreasing order) a sub-set of entries in an array using 
 * quicksort, by permutating the order of the elements in the subset rather 
 * than in the whole original array.
 * @param sz length of the permutation array.
 * @param perm array of indices of elements that we want to sort.
 * @param elem array (of length at least max(perm[k]:k=0,...,sz-1)) containing
 * the elements to be sorted.
 * @note The array of elements is not changed by this function.
 * @note This code is based in concorde's implementation of
 * permutation-quick-sort.
 * */
void uint32_EGutilPermSort2 (const size_t sz,
										 int*const perm,
										 const uint32_t*const elem);


/* ========================================================================= */
#endif
