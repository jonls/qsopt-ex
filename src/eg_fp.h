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
 * You should have received __A copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA 
 * */
/* ========================================================================= */
/* This header include definitions of macros to work with fixed point
 * representation of decimal numbers.
 *
 * - 2003-06-26
 * 					- First Implementation
 * */
/* ========================================================================= */
#ifndef __EG_FP_H__
#define __EG_FP_H__
#include "qs_config.h"
#include "eg_macros.h"
/* ========================================================================= */
/* EG_FP options */
/* ========================================================================= */
/* this enable check of overflows in the EGfp_t operations, this means that if
 * we perform an overflow operation, or work with overflow numbers, the resoult
 * will have the overflow bit set. */
#ifndef __EGFP_CHECK_OVERFLOW__
#define __EGFP_CHECK_OVERFLOW__ 0
#endif

/* this enable (on top of overflow check) varbose messages when overflow is
 * detected */
#ifndef __EGFP_CHECK_VERBOSE__
#define __EGFP_CHECK_VERBOSE__ 0
#endif


/* ========================================================================= */
/* local definitions, this implies acuracy */
/* ========================================================================= */

/* we choose the size_t integers because this is the usual natural length for
 * the machine, thus should be fastest choice. */
typedef long EGfp10_t;
typedef long EGfp20_t;
typedef long EGfp25_t;
typedef long EGfp28_t;

/* define how many bits there are in our structure */
#define EGFP_BIT CHAR_BIT*sizeof(long)

/* we reserve the second to last bit as an overflow mark, we are assuming that
 * the sign bit is the last one... this may be ot true in some machines, but
 * should be according to the STDC99 */
#define EGFP_OFBIT (((long)1)<<(EGFP_BIT-2))

/* this is the sign bit of the representation */
#define EGFP_SGBIT (1<<(EGFP_BIT-1))

/* this is how many bits we use for fractional representation */
#define EGFP_FRBIT10 10
#define EGFP_FRBIT20 20
#define EGFP_FRBIT25 25
#define EGFP_FRBIT28 28

/* this is how many bits we use for the integer parrt, we compute it as the
 * resoult of substract the other's fields to the total space */
#define EGFP_INBIT10 (EGFP_BIT-2-EGFP_FRBIT10)
#define EGFP_INBIT20 (EGFP_BIT-2-EGFP_FRBIT20)
#define EGFP_INBIT25 (EGFP_BIT-2-EGFP_FRBIT25)
#define EGFP_INBIT28 (EGFP_BIT-2-EGFP_FRBIT28)

/* now we define the maximum representable number */
#define EGFP_MAX10 ((double)((((long)1)<<30)-1))/((long)1<<EGFP_FRBIT10)
#define EGFP_MAX20 ((double)((((long)1)<<30)-1))/((long)1<<EGFP_FRBIT20)
#define EGFP_MAX25 ((double)((((long)1)<<30)-1))/((long)1<<EGFP_FRBIT25)
#define EGFP_MAX28 ((double)((((long)1)<<30)-1))/((long)1<<EGFP_FRBIT28)

/* now the minimum representable value */
#define EGFP_MIN10 (-EGFP_MAX10)
#define EGFP_MIN20 (-EGFP_MAX20)
#define EGFP_MIN25 (-EGFP_MAX25)
#define EGFP_MIN28 (-EGFP_MAX28)

/* this is the conversion factor used to transform __A real type to our fixed
 * point representation. */
#define EGFP_FACTOR10 (1LL<<EGFP_FRBIT10)
#define EGFP_FACTOR20 (1LL<<EGFP_FRBIT20)
#define EGFP_FACTOR25 (1LL<<EGFP_FRBIT25)
#define EGFP_FACTOR28 (1LL<<EGFP_FRBIT28)

/* the difference between one and the nearest representable number. Note that in
 * our case this number may be seen as the discretization for the numbers that
 * can be represented */
#define EGFP10_EPSILON (((double)1)/(EGFP_FACTOR10))
#define EGFP20_EPSILON (((double)1)/(EGFP_FACTOR20))
#define EGFP25_EPSILON (((double)1)/(EGFP_FACTOR25))
#define EGFP28_EPSILON (((double)1)/(EGFP_FACTOR28))

/* ========================================================================= */
/* now we define some function to manipulate this numbers */
/* ========================================================================= */

/* this check and maintain the overflow bit if the checking was enable, see the
 * eg_config.h to see the default value */
#if __EGFP_CHECK_OVERFLOW__
/* chek if __A FP_T number is overflow */
#define __EGFP_ABS__(__X) (__X>0?__X:-__X)
#define __EGFP_SGN__(__X) (__X<0?-1:1)
#if __EGFP_CHECK_VERBOSE__
#define EGFP_OCHK(__X) (((__EGFP_ABS__(__X))&EGFP_OFBIT)?fprintf(stderr,"Overflow in %s:%d\n",__FILE__,__LINE__),EGFP_OFBIT&(__EGFP_ABS__(__X)):0)
#else
#define EGFP_OCHK(__X) ((__EGFP_ABS__(__X))&EGFP_OFBIT)
#endif
/* define range check */
#define EGFP_RCHK10(__X) (((EGfp10_t)(((__X)>EGFP_MAX10)||((__X)<EGFP_MIN10)))<<(EGFP_BIT-2))
#define EGFP_RCHK20(__X) (((EGfp20_t)(((__X)>EGFP_MAX20)||((__X)<EGFP_MIN20)))<<(EGFP_BIT-2))
#define EGFP_RCHK25(__X) (((EGfp25_t)(((__X)>EGFP_MAX25)||((__X)<EGFP_MIN25)))<<(EGFP_BIT-2))
#define EGFP_RCHK28(__X) (((EGfp28_t)(((__X)>EGFP_MAX28)||((__X)<EGFP_MIN28)))<<(EGFP_BIT-2))
#else
#define EGFP_RCHK10(__X) 0
#define EGFP_RCHK20(__X) 0
#define EGFP_RCHK25(__X) 0
#define EGFP_RCHK28(__X) 0
#define EGFP_OCHK(__X) 0
#define __EGFP_SGN__(__X) 1
#endif

/* this function check if __A EG_fp_t is overflown */
#define EGfpCheckOverflow10(__X) EGFP_OCHK(__X)
#define EGfpCheckOverflow20(__X) EGFP_OCHK(__X)
#define EGfpCheckOverflow25(__X) EGFP_OCHK(__X)
#define EGfpCheckOverflow28(__X) EGFP_OCHK(__X)

/* convert __A float to our fixed point representation */
#define ftofp10(__F) (__EGFP_SGN__(__F)*(((EGfp10_t)((__F*__EGFP_SGN__(__F))*EGFP_FACTOR10))|EGFP_RCHK10(__F)))
#define ftofp20(__F) (__EGFP_SGN__(__F)*(((EGfp20_t)((__F*__EGFP_SGN__(__F))*EGFP_FACTOR20))|EGFP_RCHK20(__F)))
#define ftofp25(__F) (__EGFP_SGN__(__F)*(((EGfp25_t)((__F*__EGFP_SGN__(__F))*EGFP_FACTOR25))|EGFP_RCHK25(__F)))
#define ftofp28(__F) (__EGFP_SGN__(__F)*(((EGfp28_t)((__F*__EGFP_SGN__(__F))*EGFP_FACTOR28))|EGFP_RCHK28(__F)))

/* convert __A double to our fixed point representation */
#define lftofp10(__F) (__EGFP_SGN__(__F)*(((EGfp10_t)((__F*__EGFP_SGN__(__F))*EGFP_FACTOR10))|EGFP_RCHK10(__F)))
#define lftofp20(__F) (__EGFP_SGN__(__F)*(((EGfp20_t)((__F*__EGFP_SGN__(__F))*EGFP_FACTOR20))|EGFP_RCHK20(__F)))
#define lftofp25(__F) (__EGFP_SGN__(__F)*(((EGfp25_t)((__F*__EGFP_SGN__(__F))*EGFP_FACTOR25))|EGFP_RCHK25(__F)))
#define lftofp28(__F) (__EGFP_SGN__(__F)*(((EGfp28_t)((__F*__EGFP_SGN__(__F))*EGFP_FACTOR28))|EGFP_RCHK28(__F)))

#ifdef HAVE_LONG_DOUBLE
/* convert __A long double to our fixed point representation */
#define llftofp10(__F) (__EGFP_SGN__(__F)*(((EGfp10_t)((__F*__EGFP_SGN__(__F))*EGFP_FACTOR10))|EGFP_RCHK10(__F)))
#define llftofp20(__F) (__EGFP_SGN__(__F)*(((EGfp20_t)((__F*__EGFP_SGN__(__F))*EGFP_FACTOR20))|EGFP_RCHK20(__F)))
#define llftofp25(__F) (__EGFP_SGN__(__F)*(((EGfp25_t)((__F*__EGFP_SGN__(__F))*EGFP_FACTOR25))|EGFP_RCHK25(__F)))
#define llftofp28(__F) (__EGFP_SGN__(__F)*(((EGfp28_t)((__F*__EGFP_SGN__(__F))*EGFP_FACTOR28))|EGFP_RCHK28(__F)))
#endif

/* convert an integer to our fixed point representation */
#define itofp10(__F) (__EGFP_SGN__(__F)*(((EGfp10_t)((__F*__EGFP_SGN__(__F))*EGFP_FACTOR10))|EGFP_RCHK10(__F)))
#define itofp20(__F) (__EGFP_SGN__(__F)*(((EGfp20_t)((__F*__EGFP_SGN__(__F))*EGFP_FACTOR20))|EGFP_RCHK20(__F)))
#define itofp25(__F) (__EGFP_SGN__(__F)*(((EGfp25_t)((__F*__EGFP_SGN__(__F))*EGFP_FACTOR25))|EGFP_RCHK25(__F)))
#define itofp28(__F) (__EGFP_SGN__(__F)*(((EGfp28_t)((__F*__EGFP_SGN__(__F))*EGFP_FACTOR28))|EGFP_RCHK28(__F)))

/* convert an integer to our fixed point representation */
#define ltofp10(__F) (__EGFP_SGN__(__F)*(((EGfp10_t)((__F*__EGFP_SGN__(__F))*EGFP_FACTOR10))|EGFP_RCHK10(__F)))
#define ltofp20(__F) (__EGFP_SGN__(__F)*(((EGfp20_t)((__F*__EGFP_SGN__(__F))*EGFP_FACTOR20))|EGFP_RCHK20(__F)))
#define ltofp25(__F) (__EGFP_SGN__(__F)*(((EGfp25_t)((__F*__EGFP_SGN__(__F))*EGFP_FACTOR25))|EGFP_RCHK25(__F)))
#define ltofp28(__F) (__EGFP_SGN__(__F)*(((EGfp28_t)((__F*__EGFP_SGN__(__F))*EGFP_FACTOR28))|EGFP_RCHK28(__F)))

/* convert an integer to our fixed point representation */
#define lltofp10(__F) (__EGFP_SGN__(__F)*(((EGfp10_t)((__F*__EGFP_SGN__(__F))*EGFP_FACTOR10))|EGFP_RCHK10(__F)))
#define lltofp20(__F) (__EGFP_SGN__(__F)*(((EGfp20_t)((__F*__EGFP_SGN__(__F))*EGFP_FACTOR20))|EGFP_RCHK20(__F)))
#define lltofp25(__F) (__EGFP_SGN__(__F)*(((EGfp25_t)((__F*__EGFP_SGN__(__F))*EGFP_FACTOR25))|EGFP_RCHK25(__F)))
#define lltofp28(__F) (__EGFP_SGN__(__F)*(((EGfp28_t)((__F*__EGFP_SGN__(__F))*EGFP_FACTOR28))|EGFP_RCHK28(__F)))

/* convert anEGfp_t to double */
#define fptolf10(__F) (((double)__F)/EGFP_FACTOR10)
#define fptolf20(__F) (((double)__F)/EGFP_FACTOR20)
#define fptolf25(__F) (((double)__F)/EGFP_FACTOR25)
#define fptolf28(__F) (((double)__F)/EGFP_FACTOR28)

/* convert anEGfp_t to float */
#define fptof10(__F) (((float)__F)/EGFP_FACTOR10)
#define fptof20(__F) (((float)__F)/EGFP_FACTOR20)
#define fptof25(__F) (((float)__F)/EGFP_FACTOR25)
#define fptof28(__F) (((float)__F)/EGFP_FACTOR28)

/* convert anEGfp_t to int */
#define fptoi10(__F) (((int)__F)/EGFP_FACTOR10)
#define fptoi20(__F) (((int)__F)/EGFP_FACTOR20)
#define fptoi25(__F) (((int)__F)/EGFP_FACTOR25)
#define fptoi28(__F) (((int)__F)/EGFP_FACTOR28)

/* convert anEGfp_t to long */
#define fptol10(__F) (((long)__F)/EGFP_FACTOR10)
#define fptol20(__F) (((long)__F)/EGFP_FACTOR20)
#define fptol25(__F) (((long)__F)/EGFP_FACTOR25)
#define fptol28(__F) (((long)__F)/EGFP_FACTOR28)

/* convert anEGfp_t to long long */
#define fptoll10(__F) (((long long)__F)/EGFP_FACTOR10)
#define fptoll20(__F) (((long long)__F)/EGFP_FACTOR20)
#define fptoll25(__F) (((long long)__F)/EGFP_FACTOR25)
#define fptoll28(__F) (((long long)__F)/EGFP_FACTOR28)

/* this function add to FP numbers */
#if __EGFP_CHECK_OVERFLOW__
#define EGfpAdd(__A,__B) ((EGFP_OCHK(__A)|EGFP_OCHK(__B)|EGFP_OCHK((__A+__B))|((__A+__B)*__EGFP_SGN__(__A+__B)))*__EGFP_SGN__(__A+__B))
#else
#define EGfpAdd(__A,__B) (__A+__B)
#endif
#define EGfpAdd10(__A,__B) EGfpAdd(__A,__B)
#define EGfpAdd20(__A,__B) EGfpAdd(__A,__B)
#define EGfpAdd25(__A,__B) EGfpAdd(__A,__B)
#define EGfpAdd28(__A,__B) EGfpAdd(__A,__B)

/* this function substract to FP numbers */
#if __EGFP_CHECK_OVERFLOW__
#define EGfpSub(__A,__B) ((EGFP_OCHK(__A)|EGFP_OCHK(__B)|EGFP_OCHK((__A-__B))|((__A-__B)*__EGFP_SGN__(__A-__B)))*__EGFP_SGN__(__A-__B))
#else
#define EGfpSub(__A,__B) (__A-__B)
#endif
#define EGfpSub10(__A,__B) EGfpSub(__A,__B)
#define EGfpSub20(__A,__B) EGfpSub(__A,__B)
#define EGfpSub25(__A,__B) EGfpSub(__A,__B)
#define EGfpSub28(__A,__B) EGfpSub(__A,__B)

/* this define multiplication of FP numbers */
#define __EGfpIntMul10__(__A,__B) ((long)(((long long)(__A))*((long long)(__B))>>EGFP_FRBIT10))
#define __EGfpIntMul20__(__A,__B) ((long)(((long long)(__A))*((long long)(__B))>>EGFP_FRBIT20))
#define __EGfpIntMul25__(__A,__B) ((long)(((long long)(__A))*((long long)(__B))>>EGFP_FRBIT25))
#define __EGfpIntMul28__(__A,__B) ((long)(((long long)(__A))*((long long)(__B))>>EGFP_FRBIT28))
#if __EGFP_CHECK_OVERFLOW__
#define EGfpMul10(__A,__B) ((\
		(EGFP_OCHK(__A))|\
		(EGFP_OCHK(__B))|\
		(EGFP_RCHK10((((double)(__A))*((double)(__B)))/(EGFP_FACTOR10*EGFP_FACTOR10)))|\
		(__EGfpIntMul10__(__A,__B)*__EGFP_SGN__(__EGfpIntMul10__(__A,__B))))*__EGFP_SGN__(__EGfpIntMul10__(__A,__B)))
#define EGfpMul20(__A,__B) ((\
		(EGFP_OCHK(__A))|\
		(EGFP_OCHK(__B))|\
		(EGFP_RCHK20((((double)(__A))*((double)(__B)))/(EGFP_FACTOR20*EGFP_FACTOR20)))|\
		(__EGfpIntMul20__(__A,__B)*__EGFP_SGN__(__EGfpIntMul20__(__A,__B))))*__EGFP_SGN__(__EGfpIntMul20__(__A,__B)))
#define EGfpMul25(__A,__B) ((\
		(EGFP_OCHK(__A))|\
		(EGFP_OCHK(__B))|\
		(EGFP_RCHK25((((double)(__A))*((double)(__B)))/(EGFP_FACTOR25*EGFP_FACTOR25)))|\
		(__EGfpIntMul25__(__A,__B)*__EGFP_SGN__(__EGfpIntMul25__(__A,__B))))*__EGFP_SGN__(__EGfpIntMul25__(__A,__B)))
#define EGfpMul28(__A,__B) ((\
		(EGFP_OCHK(__A))|\
		(EGFP_OCHK(__B))|\
		(EGFP_RCHK28((((double)(__A))*((double)(__B)))/(EGFP_FACTOR28*EGFP_FACTOR28)))|\
		(__EGfpIntMul28__(__A,__B)*__EGFP_SGN__(__EGfpIntMul28__(__A,__B))))*__EGFP_SGN__(__EGfpIntMul28__(__A,__B)))
#else
#define EGfpMul10(__A,__B) __EGfpIntMul10__(__A,__B)
#define EGfpMul20(__A,__B) __EGfpIntMul20__(__A,__B)
#define EGfpMul25(__A,__B) __EGfpIntMul25__(__A,__B)
#define EGfpMul28(__A,__B) __EGfpIntMul28__(__A,__B)
#endif

/* this define divition of FP numbers */
#define __EGfpIntDiv10__(__A,__B) ((EGfp10_t)((((long long)(__A))<<EGFP_FRBIT10)/(__B)))
#define __EGfpIntDiv20__(__A,__B) ((EGfp20_t)((((long long)(__A))<<EGFP_FRBIT20)/(__B)))
#define __EGfpIntDiv25__(__A,__B) ((EGfp25_t)((((long long)(__A))<<EGFP_FRBIT25)/(__B)))
#define __EGfpIntDiv28__(__A,__B) ((EGfp28_t)((((long long)(__A))<<EGFP_FRBIT28)/(__B)))
#if __EGFP_CHECK_OVERFLOW__
#define EGfpDiv10(__A,__B) ((\
		(EGFP_OCHK(__A))|\
		(EGFP_OCHK(__B))|\
		(EGFP_RCHK10(((((((double)(__A))*EGFP_FACTOR10)/(__B)))/EGFP_FACTOR10)))|\
		(__EGfpIntDiv10__(__A,__B)*__EGFP_SGN__(__EGfpIntDiv10__(__A,__B))))*__EGFP_SGN__(__EGfpIntDiv10__(__A,__B)))
#define EGfpDiv20(__A,__B) ((\
		(EGFP_OCHK(__A))|\
		(EGFP_OCHK(__B))|\
		(EGFP_RCHK20(((((((double)(__A))*EGFP_FACTOR20)/(__B)))/EGFP_FACTOR20)))|\
		(__EGfpIntDiv20__(__A,__B)*__EGFP_SGN__(__EGfpIntDiv20__(__A,__B))))*__EGFP_SGN__(__EGfpIntDiv20__(__A,__B)))
#define EGfpDiv25(__A,__B) ((\
		(EGFP_OCHK(__A))|\
		(EGFP_OCHK(__B))|\
		(EGFP_RCHK25(((((((double)(__A))*EGFP_FACTOR25)/(__B)))/EGFP_FACTOR25)))|\
		(__EGfpIntDiv25__(__A,__B)*__EGFP_SGN__(__EGfpIntDiv25__(__A,__B))))*__EGFP_SGN__(__EGfpIntDiv25__(__A,__B)))
#define EGfpDiv28(__A,__B) ((\
		(EGFP_OCHK(__A))|\
		(EGFP_OCHK(__B))|\
		(EGFP_RCHK28(((((((double)(__A))*EGFP_FACTOR28)/(__B)))/EGFP_FACTOR28)))|\
		(__EGfpIntDiv28__(__A,__B)*__EGFP_SGN__(__EGfpIntDiv28__(__A,__B))))*__EGFP_SGN__(__EGfpIntDiv28__(__A,__B)))
#else
#define EGfpDiv10(__A,__B) __EGfpIntDiv10__(__A,__B)
#define EGfpDiv20(__A,__B) __EGfpIntDiv20__(__A,__B)
#define EGfpDiv25(__A,__B) __EGfpIntDiv25__(__A,__B)
#define EGfpDiv28(__A,__B) __EGfpIntDiv28__(__A,__B)
#endif

/* this defines the change of sign of __A FP number */
#define EGfpMinus(__A) (__A*-1)
#define EGfpMinus10(__A) EGfpMinus(__A)
#define EGfpMinus20(__A) EGfpMinus(__A)
#define EGfpMinus25(__A) EGfpMinus(__A)
#define EGfpMinus28(__A) EGfpMinus(__A)

/* end of eg_fp.h */
#endif
