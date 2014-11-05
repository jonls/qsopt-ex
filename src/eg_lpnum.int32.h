/* EGlib "Efficient General Library" provides some basic structures and
 * algorithms commons in many optimization algorithms.
 *
 * Copyright (C) 2005 Daniel Espinoza and Marcos Goycoolea.
 * 
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation; either version 2.1 of the License,or (at your
 * option) any later version.
 *
 * This library is distributed in the hope that it will be useful,but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public 
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not,write to the Free Software Foundation,
 * Inc.,51 Franklin St,Fifth Floor,Boston,MA  02110-1301  USA 
 * */
#ifndef __EG_LPNUM_INT32__
#define __EG_LPNUM_INT32__
#include "qs_config.h"
#include "eg_lpnum.h"
/** @file
 * @ingroup EGlpNum */
/** @addtogroup EGlpNum */
/** @{ */
/*=========================================================================*/
/** extern definitions of constaants for different set-ups */
#define int32_zeroLpNum 0
#define int32_oneLpNum  1
#define int32_epsLpNum  0
#define int32_MaxLpNum INT32_MAX
#define int32_MinLpNum INT32_MIN

/*=========================================================================*/
/** @brief Read from a string a number and store it in the given number,
 * @return the number of chars readed from the input string */
#define int32_EGlpNumReadStr(__EGPa__,__EGstr__) ({\
	int __i=0;\
	sscanf(__EGstr__,"%d%n",&(__EGPa__),&__i);\
	__i;})

/*=========================================================================*/
/** @brief given a number,write it to a string (to be allocated internally),
 * and return it. */
#define int32_EGlpNumGetStr(__EGPa__) ({\
	char *__str=0;\
	size_t __i=snprintf(__str,(size_t)0,"%d",__EGPa__);\
	__str=EGsMalloc(char,__i+1);\
	snprintf(__str,__i+1,"%d",__EGPa__);\
	__str;})

/*=========================================================================*/
/** @brief given an array of type number,free it,if the pointer is NULL
 * nothing happen. */
#define int32_EGlpNumFreeArray(__EGea__) __EGlpNumFreeArray(__EGea__)

/*=========================================================================*/
/** @brief Reallocate and initialize (if needed) '__EGlsz__' elements of type 
 * int32_t and return it,if no more memory,exit(1) */
#define int32_EGlpNumReallocArray(__EGpone__,__EGlsz__) ({\
	size_t __csz=(__EGlsz__),*__usp=0;\
	size_t __psz=__EGlpNumArraySize(*(__EGpone__));\
	int32_t** __ptr__=(__EGpone__);\
	if (!__psz) *__ptr__=int32_EGlpNumAllocArray (__csz); \
	else if (__psz < __csz) {\
		__usp=(size_t*)(*__ptr__);\
		__usp--;\
		__usp=EGrealloc(__usp,sizeof(int32_t)*__csz+sizeof(size_t));\
		__usp[0]=__csz;\
		*__ptr__=(int32_t*)(__usp+1);\
		memset((*__ptr__)+__psz,0,sizeof(int32_t)*(__csz-__psz));\
	}\
	*__ptr__;})

/*=========================================================================*/
/** @brief Allocate and initialize (if needed) '__EGsz__' elements of type int32_t
 * and return it,if no more memory,exit(1) */
#define int32_EGlpNumAllocArray(__EGsz__) __EGlpNumAllocArray(int32_t,__EGsz__)

/*=========================================================================*/
/** @brief set the given number pointer,set its value to the given double.
 * @param (__EGvone__) double where we will store the double value.
 * @param __EGdbl__ double value to be stored in '(__EGvone__)'.
 * @par Description:
 * This function is intended to set initial values to variables; note that the
 * double is (__EGPa__) number and not (__EGPa__) pointer to that value,be carefull with this
 * detail. Also,due to implementation details this function can't deal with
 * numbers above 1e158 or smaller than 1e-158. Note also that if the number is
 * writen in the form \f$x=\bar{x}\cdot 2^e\f$ with \f$0.5<|\bar{x}|<1\f$,
 * then \f$\left|x-\frac{p}{q}\right|<2^{e-64}\f$.
 * */
#define int32_EGlpNumSet(__EGvone__,__EGdbl__) ((__EGvone__)=(int32_t)(__EGdbl__))

/*=========================================================================*/
/** @brief Stores in the first number the ceil value of the second number,i.e.
 * EGlpNumCeil(__EGPa__,__EGPb__) <==> (__EGPa__)=ceil(__EGPb__) */
#define int32_EGlpNumCeil(__EGPa__,__EGPb__) (__EGPa__=__EGPb__)

/*=========================================================================*/
/** @brief Stores in the first number the floor value of the second number,i.e.
 * EGlpNumFloor(__EGPa__,__EGPb__) <==> (__EGPa__)=floor(__EGPb__) */
#define int32_EGlpNumFloor(__EGPa__,__EGPb__) (__EGPa__=__EGPb__)

/*=========================================================================*/
/** @brief store the (multiplicative) inverse of (__EGPa__) number to itself,i.e.
 * implement (__EGPa__)=1/(__EGPa__).
 * @param (__EGPa__) the number to be inverted. */
#define int32_EGlpNumInv(__EGPa__) ((__EGPa__)=0)

/*=========================================================================*/
/** @brief Compare if two numbers are equal within (__EGPa__) maximum (__EGPerr__).
 * @param (__EGPa__) EGlpNum_t first number to compare.
 * @param (__EGPb__) EGlpNum_t second number to compare.
 * @return int one in success,zero oterwise.
 * @par Description:
 * Given two numbers '(__EGPa__)','(__EGPb__)' return 1 if (__EGPa__)==(__EGPb__),otherwise it return 0
 * */
#define int32_EGlpNumIsEqqual(__EGPa__,__EGPb__) ((__EGPa__)==(__EGPb__))

/*=========================================================================*/
/** @brief Compare if two numbers are equal within (__EGPa__) maximum (__EGPerr__).
 * @param (__EGPa__) EGlpNum_t first number to compare.
 * @param (__EGPb__) EGlpNum_t second number to compare.
 * @param (__EGPerr__) EGlpNum_t maximum difference allowed between both
 * numbers.
 * @return int one in success,zero oterwise.
 * @par Description:
 * Given two numbers '(__EGPa__)','(__EGPb__)' and (__EGPa__) tolerance '(__EGPerr__)',
 * return 1 if |(__EGPa__)-(__EGPb__)|<=(__EGPerr__),otherwise it return 0.
 * */
#define int32_EGlpNumIsEqual(__EGPa__,__EGPb__,__EGPerr__) (EGabs((__EGPa__)-(__EGPb__))<=(__EGPerr__))
#define int32_EGlpNumIsNeq(__EGPa__,__EGPb__,__EGPerr__) (EGabs((__EGPa__)-(__EGPb__))>(__EGPerr__))
#define int32_EGlpNumIsNeqq(__EGPa__,__EGPb__)  ((__EGPa__)!=(__EGPb__))
#define int32_EGlpNumIsNeqZero(__EGPa__,__EGPerr__) (((__EGPa__)>(__EGPerr__))||(-(__EGPa__) > (__EGPerr__)))
#define int32_EGlpNumIsNeqqZero(__EGPa__)     	(__EGPa__)

/*=========================================================================*/
/** @brief test if the first number is greater than zero
 * @param (__EGPa__) EGlpNum_t number to test.
 * @return int one if success,zero otherwise.
 * */
#define int32_EGlpNumIsGreatZero(__EGPa__) ((__EGPa__)>0)

/*=========================================================================*/
/** @brief test if the first number is less than zero
 * @param (__EGPa__) EGlpNum_t number to test.
 * @return int one if success,zero otherwise.
 * */
#define int32_EGlpNumIsLessZero(__EGPa__) ((__EGPa__)<0)

/*=========================================================================*/
/** @brief test if the first number is bigger to the second number
 * @param (__EGPa__) EGlpNum_t the first number.
 * @param (__EGPb__) EGlpNum_t the second number
 * @return int one if success,zero otherwise.
 * @par Description:
 * Given two numbers '(__EGPa__)' and '(__EGPb__)',return one if (__EGPa__) < (__EGPb__),zero
 * otherwise.
 * */
#define int32_EGlpNumIsLess(__EGPa__,__EGPb__) ((__EGPa__)<(__EGPb__))

/*=========================================================================*/
/** @brief test if the sum of the first two numbers is less thatn the third
 * number.
 * @param (__EGPa__) EGlpNum_t the first number.
 * @param (__EGPb__) EGlpNum_t the second number
 * @param (__EGPc__) EGlpNum_t the third number
 * @return int one if success,zero otherwise.
 * @par Description:
 * Given (__EGPa__),__EGPb__,and (__EGPc__),return nonzero if ((__EGPa__) + (__EGPb__) < (__EGPc__)),zero toherwise.
 * */
#define int32_EGlpNumIsSumLess(__EGPa__,__EGPb__,__EGPc__) ((__EGPa__)+(__EGPb__)<(__EGPc__))

/*=========================================================================*/
/** @brief test if the diference of the first two numbers is less thatn the 
 * third number.
 * @param (__EGPa__) EGlpNum_t the first number.
 * @param (__EGPb__) EGlpNum_t the second number
 * @param (__EGPc__) EGlpNum_t the third number
 * @return int one if success,zero otherwise.
 * @par Description:
 * Given (__EGPa__),__EGPb__,and (__EGPc__),return nonzero if ((__EGPa__) - (__EGPb__) < (__EGPc__)),zero toherwise.
 * */
#define int32_EGlpNumIsDiffLess(__EGPa__,__EGPb__,__EGPc__) ((__EGPa__)-(__EGPb__)<(__EGPc__))

/*=========================================================================*/
/** @brief test if the first number is bigger to the second number
 * @param (__EGPa__) EGlpNum_t the first number.
 * @param (__EGPb__) double the second number
 * @return int one if success,zero otherwise.
 * @par Description:
 * Given two numbers '(__EGPa__)' and '(__EGPb__)',return one if (__EGPa__) < (__EGPb__),zero
 * otherwise.
 * */
#define int32_EGlpNumIsLessDbl(__EGPa__,__EGPb__) (((double)(__EGPa__)) < (__EGPb__))

/*=========================================================================*/
/** @brief test if the first number is bigger to the second number
 * @param (__EGPa__) EGlpNum_t the first number.
 * @param (__EGPb__) double the second number
 * @return int one if success,zero otherwise.
 * @par Description:
 * Given two numbers '(__EGPa__)' and '(__EGPb__)',return one if (__EGPa__) > (__EGPb__),zero
 * otherwise.
 * */
#define int32_EGlpNumIsGreaDbl(__EGPa__,__EGPb__) (((double)(__EGPa__)) > (__EGPb__))

/*=========================================================================*/
/** @brief test if the first number is bigger to the second number
 * @param (__EGPa__) EGlpNum_t the first number.
 * @param (__EGPb__) EGlpNum_t the second number
 * @return int one if success,zero otherwise.
 * @par Description:
 * Given two numbers '(__EGPa__)' and '(__EGPb__)',return one if (__EGPa__) <=(__EGPb__),zero
 * otherwise.
 * */
#define int32_EGlpNumIsLeq(__EGPa__,__EGPb__) ((__EGPa__)<=(__EGPb__))

/*=========================================================================*/
/** @brief copy the value of the second number to the first.
 * @param (__EGPa__) EGlpNum_t source number (it won't change value).
 * @param (__EGPb__) EGlpNum_t source number (it won't change value).
 * @param (__EGPdn__) EGlpNum_t denominator of the difference (it won't change value).
 * @param (__EGPd__) EGlpNum_t where to store the value .
 * @par Description:
 * Set (__EGPd__)=((__EGPa__) - (__EGPb__)) / (__EGPdn__) */
#define int32_EGlpNumCopyDiffRatio(__EGPd__,__EGPa__,__EGPb__,__EGPdn__) ((__EGPd__)=((__EGPa__)-(__EGPb__))/(__EGPdn__))

/*=========================================================================*/
/** @brief copy the value of the second number to the first.
 * @param (__EGPa__) EGlpNum_t source number (it won't change value).
 * @param (__EGPb__) EGlpNum_t source number (it won't change value).
 * @param (__EGPd__) EGlpNum_t where to store the value stored in '(__EGPo__)'.
 * @par Description:
 * Set (__EGPd__)=(__EGPa__) - (__EGPb__) */
#define int32_EGlpNumCopyDiff(__EGPd__,__EGPa__,__EGPb__) ((__EGPd__)=((__EGPa__)-(__EGPb__)))

/*=========================================================================*/
/** @brief copy the value of the sum of the second and third parameter
 * @param (__EGPa__) EGlpNum_t source number (it won't change value).
 * @param (__EGPb__) EGlpNum_t source number (it won't change value).
 * @param (__EGPd__) EGlpNum_t where to store the sum.
 * @par Description:
 * Set (__EGPd__)=(__EGPa__) + (__EGPb__) */
#define int32_EGlpNumCopySum(__EGPd__,__EGPa__,__EGPb__) ((__EGPd__)=((__EGPa__)+(__EGPb__)))

/*=========================================================================*/
/** @brief copy the value of the second number to the first.
 * @param (__EGPo__) EGlpNum_t source number (it won't change value).
 * @param (__EGPd__) EGlpNum_t where to store the value stored in '(__EGPo__)'.
 * @par Description:
 * Given two numbers copy the values in '(__EGPo__)',into '(__EGPd__)'.
 * */
#define int32_EGlpNumCopy(__EGPd__,__EGPo__) ((__EGPd__)=(__EGPo__))

/*=========================================================================*/
/** @brief change the fist number to the maximum between itself and the 
 * absolute value of the second.
 * @param (__EGPo__) EGlpNum_t source number (it won't change value).
 * @param (__EGPd__) EGlpNum_t where to store the value stored in '(__EGPo__)'.
 * @par Description:
 * implement (__EGPd__)=max(__EGPd__,abs(__EGPo__))
 * */
#define int32_EGlpNumSetToMaxAbs(__EGPd__,__EGPo__) ({\
	const int32_t __EGcv__=(__EGPo__)<0?(__EGPo__):-(__EGPo__);\
	__EGPd__=(__EGcv__)>(__EGPd__)?(__EGcv__):(__EGPd__);})
#define int32_EGlpNumSetToMinAbs(__EGPd__,__EGPo__) ({\
	const int32_t __EGcv__=(__EGPo__)<0?(__EGPo__):-(__EGPo__);\
	__EGPd__=(__EGcv__)<(__EGPd__)?(__EGcv__):(__EGPd__);})

/*=========================================================================*/
/** @brief copy the square of the second argument,divided by the third 
 * argument into the first argument.
 * @param (__EGPd__) EGlpNum_t where to store the result
 * @param (__EGPo__) EGlpNum_t second parameter
 * @param (__EGPdn__) EGlpNum_t third parameter
 * @par Description:
 * compute (__EGPd__)=((__EGPo__)*(__EGPo__))/(__EGPdn__)
 * */
#define int32_EGlpNumCopySqrOver(__EGPd__,__EGPo__,__EGPdn__) ((__EGPd__)=(__EGPo__)*(__EGPo__)/(__EGPdn__))

/*=========================================================================*/
/** @brief copy the value of the absolute value of the second parameter to the 
 * first parameter.
 * @param (__EGPo__) EGlpNum_t source number (it won't change value).
 * @param (__EGPd__) EGlpNum_t where to store the absolute value stored
 * in '(__EGPo__)'.
 * @par Description:
 * Given (__EGPa__) number '(__EGPo__)',copy its absolute value to '(__EGPd__)'. i.e.
 * (__EGPd__)=|(__EGPo__)|
 * */
#define int32_EGlpNumCopyAbs(__EGPd__,__EGPo__) ((__EGPd__)=((__EGPo__)<0?-(__EGPo__):(__EGPo__)))

/*=========================================================================*/
/** @brief copy minus the value of the second parameter to the 
 * first parameter.
 * @param (__EGPo__) EGlpNum_t the source number (it won't change value).
 * @param (__EGPd__) EGlpNum_t where to store minus the value stored
 * in '(__EGPo__)'.
 * @par Description:
 * Given (__EGPa__) number '(__EGPo__)',copy minus the value to '(__EGPd__)'. i.e.
 * (__EGPd__)=-(__EGPo__)
 * */
#define int32_EGlpNumCopyNeg(__EGPd__,__EGPo__) ((__EGPd__)=-(__EGPo__))

/*=========================================================================*/
/** @brief Set des=(__EGPone__)/(__EGPtwo__).
 * @param (__EGPd__) EGlpNum_t where we will store the result.
 * @param (__EGPone__) EGlpNum_t numerator of the fraction (possibly non an integer)
 * @param (__EGPtwo__) EGlpNum_t denominator of the fraction (possibly non an integer)
 * @par Description:
 *  Set des=(__EGPone__)/(__EGPtwo__)
 * */
#define int32_EGlpNumCopyFrac(__EGPd__,__EGPone__,__EGPtwo__) ((__EGPd__)=(__EGPone__)/(__EGPtwo__))

/*=========================================================================*/
/** @brief copy the first '(__EGsz__)' values in the second array to the first array.
 * @param (__EGPo__) EGlpNum_t* pointer to the array from where we will copy the
 * values (it won't change value).
 * @param (__EGPd__) EGlpNum_t* pointer to where to store the first '(__EGsz__)' values 
 * stored in '(__EGPo__)'.
 * @param (__EGsz__) unsigned int specifying how many values of '(__EGPo__)' will be copied
 * onto '(__EGPd__)'
 * @par Description:
 * This function is provided to (possible) make fast copies of arrays of
 * numbers,the arrays should be of length at least '(__EGsz__)',and the resulting
 * copy is absolutely independent froom the original,any change in one vale of
 * one array won't change values on the other array.
 * */
#define int32_EGlpNumCopyArray(__EGPd__,__EGPo__,__EGsz__) memcpy(__EGPd__,__EGPo__,sizeof(int32_t)*(__EGsz__))

/*=========================================================================*/
/** @brief Sub to (__EGPa__) given number the product of two numbers.
 * @param (__EGPa__) EGlpNum_t the number that we are going to Sub to.
 * @param (__EGPb__) EGlpNum_t value to be multiplyed.
 * @param (__EGPc__) EGlpNum_t value to be multiplyed.
 * @par Description:
 * This function implements (__EGPa__)=(__EGPa__) - (__EGPb__)*(__EGPc__),and clearly don't change the value
 * stored in '(__EGPb__)' nor in '(__EGPc__)'.
 * */
#define int32_EGlpNumSubInnProdTo(__EGPa__,__EGPb__,__EGPc__) ((__EGPa__)-=(__EGPb__)*(__EGPc__))

/*=========================================================================*/
/** @brief Add to (__EGPa__) given number the product of two numbers.
 * @param (__EGPa__) EGlpNum_t the number that we are going to add to.
 * @param (__EGPb__) EGlpNum_t value to be multiplyed.
 * @param (__EGPc__) EGlpNum_t value to be multiplyed.
 * @par Description:
 * This function implements (__EGPa__)=(__EGPa__) + (__EGPb__)*(__EGPc__),and clearly don't change the value
 * stored in '(__EGPb__)' nor in '(__EGPc__)'.
 * */
#define int32_EGlpNumAddInnProdTo(__EGPa__,__EGPb__,__EGPc__) ((__EGPa__)+=(__EGPb__)*(__EGPc__))

/*=========================================================================*/
/** @brief Substract to (__EGPa__) given number the value of the second number.
 * @param (__EGPa__) EGlpNum_t the number that we are going to substract to.
 * @param (__EGPb__) unsigned int value to be substracted to '(__EGPa__)'.
 * @par Description:
 * This function implements (__EGPa__)=(__EGPa__) - (__EGPb__),and clearly don't change the value
 * stored in '(__EGPb__)'.
 * */
#define int32_EGlpNumSubUiTo(__EGPa__,__EGPb__) ((__EGPa__)-=((int32_t)(__EGPb__)))

/*=========================================================================*/
/** @brief Add to (__EGPa__) given number the value of the second number.
 * @param (__EGPa__) EGlpNum_t the number that we are going to add to.
 * @param (__EGPb__) unsigned int value to be added to '(__EGPa__)'.
 * @par Description:
 * This function implements (__EGPa__)=(__EGPa__) + (__EGPb__),and clearly don't change the value
 * stored in '(__EGPb__)'.
 * */
#define int32_EGlpNumAddUiTo(__EGPa__,__EGPb__) ((__EGPa__)+=((int32_t)(__EGPb__)))

/*=========================================================================*/
/** @brief Add to (__EGPa__) given number the value of the second number.
 * @param (__EGPa__) EGlpNum_t the number that we are going to add to.
 * @param (__EGPb__) EGlpNum_t value to be added to '(__EGPa__)'.
 * @par Description:
 * This function implements (__EGPa__)=(__EGPa__) + (__EGPb__),and clearly don't change the value
 * stored in '(__EGPb__)'.
 * */
#define int32_EGlpNumAddTo(__EGPa__,__EGPb__) ((__EGPa__)+=(__EGPb__))

/*=========================================================================*/
/** @brief Substract to (__EGPa__) given number the value of the second number.
 * @param (__EGPa__) EGlpNum_t the number that we are going to substract
 * from.
 * @param (__EGPb__) EGlpNum_t value to be substracted to '(__EGPa__)'.
 * @par Description:
 * This function implements (__EGPa__)=(__EGPa__) - (__EGPb__),and clearly don't change the value
 * stored in '(__EGPb__)'.
 * */
#define int32_EGlpNumSubTo(__EGPa__,__EGPb__) ((__EGPa__)-=(__EGPb__))

/*=========================================================================*/
/** @brief Multiply (__EGPa__) given number by the value of the second number.
 * @param (__EGPa__) EGlpNum_t the number that we are going to multiply by
 * the second number and store the result.
 * @param (__EGPb__) EGlpNum_t value to be multyply to '(__EGPa__)'.
 * @par Description:
 * This function implements (__EGPa__)=(__EGPa__) * (__EGPb__),and clearly don't change the value
 * stored in '(__EGPb__)'.
 * */
#define int32_EGlpNumMultTo(__EGPa__,__EGPb__) ((__EGPa__)*=(__EGPb__))

/*=========================================================================*/
/** @brief Divide (__EGPa__) given number by the value of the second number.
 * @param (__EGPa__) EGlpNum_t the number that we are going to divide by
 * the second number and store the result.
 * @param (__EGPb__) EGlpNum_t value to be divide to '(__EGPa__)'.
 * @par Description:
 * This function implements (__EGPa__)=(__EGPa__) / (__EGPb__),and clearly don't change the value
 * stored in '(__EGPb__)'.
 * */
#define int32_EGlpNumDivTo(__EGPa__,__EGPb__) ((__EGPa__)/=(__EGPb__))

/*=========================================================================*/
/** @brief Divide (__EGPa__) given number by the value of the second number.
 * @param (__EGPa__) EGlpNum_t the number that we are going to divide by
 * the second number and store the result.
 * @param (__EGPb__) unsigned int value to be divided to '(__EGPa__)'.
 * @par Description:
 * This function implements (__EGPa__)=(__EGPa__) / (__EGPb__),and don't change the value
 * stored in '(__EGPb__)'.
 * */
#define int32_EGlpNumDivUiTo(__EGPa__,__EGPb__) ((__EGPa__)/=((int32_t)(__EGPb__)))

/*=========================================================================*/
/** @brief Multiply (__EGPa__) given number by the value of the second number.
 * @param (__EGPa__) EGlpNum_t the number that we are going to multiply by
 * the second number and store the result.
 * @param (__EGPb__) unsigned int value to be multyply to '(__EGPa__)'.
 * @par Description:
 * This function implements (__EGPa__)=(__EGPa__) * (__EGPb__),and clearly don't change the value
 * stored in '(__EGPb__)'.
 * */
#define int32_EGlpNumMultUiTo(__EGPa__,__EGPb__) ((__EGPa__)*=((int32_t)(__EGPb__)))

/*=========================================================================*/
/** @brief Reset the value of the pointed number to zero.
 * @param (__EGPa__) EGlpNum_t the value to be set to zero.
 * @par Descrpition:
 * Reset (__EGPa__) to zero,i.e. implements (__EGPa__)=0;
 * */
#define int32_EGlpNumZero(__EGPa__) ((__EGPa__)=0)

/*=========================================================================*/
/** @brief Reset the value of the pointed number to one.
 * @param (__EGPa__) EGlpNum_t value to be set to one.
 * @par Descrpition:
 * Reset (__EGPa__) to zero,i.e. implements (__EGPa__)=1;
 * */
#define int32_EGlpNumOne(__EGPa__) ((__EGPa__)=1)

/*=========================================================================*/
/** @brief Change the sign of the number.
 * @param (__EGPa__) EGlpNum_t number we will change sign.
 * @par Descrpition:
 * Change the sign of the given number,i.e. implements (__EGPa__)=-(__EGPa__)
 * */
#define int32_EGlpNumSign(__EGPa__) ((__EGPa__)=-(__EGPa__))

/*=========================================================================*/
/** @brief return the closest double value of the given pointer number.
 * @param (__EGPa__) EGlpNum_t number that we will be transformed to double.
 * @return double the closest double representation of the given number.
 * par Description:
 * return the double number closest in value to the value stored in (__EGPa__).
 * */
#define int32_EGlpNumToLf(__EGPa__) ((double)(__EGPa__))

/*=========================================================================*/
/** @brief initialize the internal memory of (__EGPa__) given variable */
#define int32_EGlpNumInitVar(__EGPa__) int32_EGlpNumZero(__EGPa__)

/*=========================================================================*/
/** @brief free the internal memory of (__EGPa__) given variable */
#define int32_EGlpNumClearVar(__EGPa__)

/* ========================================================================= */
/** @} */
#endif
