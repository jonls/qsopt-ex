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

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#include <gmp.h>

#include "qs_config.h"
#include "logging-private.h"

#include "eg_memslab.h"
#include "eg_nummacros.h"

/** @file
 * @ingroup EGlpNum */
/** @addtogroup EGlpNum */
/** @{ */
/* ========================================================================= */
/** @brief This constant define the acuracy required while
 * converting doubles to rationals, a good number is 1e-5. More exactly, we 
 * stop the continued fraction method whenever the next e_i-[e_i] computed is
 * less than EGLPNUM_MINEPS. Note that this value can't be smaller than
 * 1/ULONG_MAX, otherwise we will have problems in the confertion step. */
#ifndef EGLPNUM_MINEPS
#define EGLPNUM_MINEPS 0x1ep-20
#else
#if EGLPNUM_MINEPS < 3e-10
#undef EGLPNUM_MINEPS
#define EGLPNUM_MINEPS 3e-10
#endif
#endif


/* ========================================================================= */
/** @brief if non-zero, use slab-pool allocator for GMP, otherwise, use malloc/
 * realloc / free, */
#ifndef EG_LPNUM_MEMSLAB
#define EG_LPNUM_MEMSLAB 1
#endif
/* ========================================================================= */
/** @brief type-dependant constants and helper numbers @{ */
mpz_t __zeroLpNum_mpz__;
mpz_t __oneLpNum_mpz__;
mpz_t __MaxLpNum_mpz__;
mpz_t __MinLpNum_mpz__;
mpq_t __zeroLpNum_mpq__;
mpq_t __oneLpNum_mpq__;
mpq_t __MaxLpNum_mpq__;
mpq_t __MinLpNum_mpq__;
mpf_t __zeroLpNum_mpf__;
mpf_t __MaxLpNum_mpf__;
mpf_t __MinLpNum_mpf__;
mpf_t __oneLpNum_mpf__;
mpf_t mpf_eps;
unsigned long int EGLPNUM_PRECISION = 128;

/** @} */

/* ========================================================================= */
static int __EGlpNum_setup=0;
/* ========================================================================= */
/** @name data to handle memory allocations within gmp */
/** @{*/
#define __GMP_MEM_VERBOSE 100
#define __GMP_MEM_STATS__ 0
#define __GMP_MEM_MAX__ 256
static const uint8_t _EGgmpPlTable[257] = \
	{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,	/* pool for 0-16 bytes */
		1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,	/* pool for 17-32 bytes */
		2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,	/* pool for 33-64 bytes */
		2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
		3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,	/* pool for 65-128 bytes */
		3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
		3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
		3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
		4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,	/* pool for 129-255 bytes */
		4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
		4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
		4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
		4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
		4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
		4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
		4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4};
static const size_t _EGgmpPlSz[5] = {16,32,64,128,256};
#define __GMP_MEM_NPOOL__ 5
EGmemSlabPool_t EGgmpPl[__GMP_MEM_NPOOL__];
#if __GMP_MEM_STATS__
static size_t __alloc_sz[__GMP_MEM_NPOOL__] = {0,0,0,0,0};
static size_t __rlloc_sz[__GMP_MEM_NPOOL__] = {0,0,0,0,0};
static size_t __totmem = 0;
static size_t __maxmem = 0;
static size_t __nallocs = 0;
static size_t __nrllocs = 0;
static size_t __nalarge = 0;
static size_t __nrlarge = 0;
static size_t __laaverage=0;
static size_t __lraverage=0;
#endif
/** @}*/
/* ========================================================================= */
#if __GMP_MEM_STATS__
#define _TACC(__m) do{\
	__totmem += __m;\
	if(__totmem > __maxmem) __maxmem = __totmem;}while(0)
#else
#define _TACC(__m)
#endif
/* ========================================================================= */
#if __GMP_MEM_STATS__
#define _AACC(__a) do{\
	__nallocs++;\
	if((__a) <= __GMP_MEM_MAX__)\
		__alloc_sz[_EGgmpPlTable[__a]]++;\
	else{\
		__nalarge++;\
		__laaverage += (__a);}}while(0)
#else
#define _AACC(__a)
#endif
/* ========================================================================= */
/** @brief dummy malloc function for gmp, at this stage is used for 
 * creating an account of the allocation */
static void* __EGgmp_malloc(size_t sz)
{
	void*ptr = 0;
	_TACC(sz);
	_AACC(sz);
	if(sz <= __GMP_MEM_MAX__)
	{
		ptr = EGmemSlabPoolAlloc(EGgmpPl+_EGgmpPlTable[sz]);
		MESSAGE(__GMP_MEM_VERBOSE,"alloc %p [%zd]",ptr, sz);
		return ptr;
		/*memset(ptr,0,sz);*/
	}
	else
	{
		ptr = malloc(sz);
		if(!ptr) EXIT(1,"No more memory");
		MESSAGE(__GMP_MEM_VERBOSE,"alloc %p [%zd]",ptr, sz);
		return ptr;
	}
	/*QSlog("allocating %p [%zd]",ptr,sz);*/
	/*return ptr;*/
}
/* ========================================================================= */
#if __GMP_MEM_STATS__
#define _RACC(__a) do{\
	__nrllocs++;\
	if(__a<__GMP_MEM_MAX__)\
		__rlloc_sz[_EGgmpPlTable[__a]]++;\
	else{\
		__nrlarge++;\
		__lraverage += (__a);}}while(0)
#else
#define _RACC(__a)
#endif
/* ========================================================================= */
/** @brief dummy realloc for gmp */
static void* __EGgmp_realloc(void*ptr,size_t osz,size_t nsz)
{
	const size_t a1 = nsz > __GMP_MEM_MAX__ ? __GMP_MEM_NPOOL__ : _EGgmpPlTable[nsz];
	const size_t a2 = osz > __GMP_MEM_MAX__ ? __GMP_MEM_NPOOL__ : _EGgmpPlTable[osz];
	const size_t msz = nsz > osz ? osz : nsz;
	void*rptr=0;
	#if __GMP_MEM_STATS__
	_RACC(nsz);
	__totmem -= osz;
	_TACC(nsz);
	#endif
	if(a1 < __GMP_MEM_NPOOL__)
	{
		if(a2 == a1)
		{
			rptr = ptr;
			MESSAGE(__GMP_MEM_VERBOSE,"realloc %p [%zd] to %p [%zd]",ptr, osz, rptr, nsz);
			return rptr;
		}
		else
		{
			rptr = EGmemSlabPoolAlloc(EGgmpPl+a1);
			memcpy(rptr,ptr,msz);
			if(a2 < __GMP_MEM_NPOOL__)
			{
				MESSAGE(__GMP_MEM_VERBOSE,"realloc %p [%zd] to %p [%zd]",ptr, osz, rptr, nsz);
				EGmemSlabPoolFree(ptr);
				return rptr;
			}
			else
			{
				MESSAGE(__GMP_MEM_VERBOSE,"realloc %p [%zd] to %p [%zd]",ptr, osz, rptr, nsz);
				EGfree(ptr);
				return rptr;
			}
		}
	}
	else if(a2 < __GMP_MEM_NPOOL__)
	{
		rptr = EGmalloc(nsz);
		memcpy(rptr,ptr,msz);
		MESSAGE(__GMP_MEM_VERBOSE,"realloc %p [%zd] to %p [%zd]",ptr, osz, rptr, nsz);
		EGmemSlabPoolFree(ptr);
		return rptr;
	}
	else
	{
		rptr = EGrealloc(ptr,nsz);
		MESSAGE(__GMP_MEM_VERBOSE,"realloc %p [%zd] to %p [%zd]",ptr, osz, rptr, nsz);
		return rptr;
	}
	/*QSlog("Re-allocating %p [%zd] to %p [%zd]",ptr,osz,rptr,nsz);*/
	/*return rptr;*/
}
/* ========================================================================= */
/** @brief dummy free for gmp */
static void __EGgmp_free(void*ptr,size_t sz)
{
	const size_t a = sz > __GMP_MEM_MAX__ ? __GMP_MEM_NPOOL__ : _EGgmpPlTable[sz];
	#if __GMP_MEM_STATS__
	__totmem -= sz;
	#endif
	MESSAGE(__GMP_MEM_VERBOSE,"freeing %p",ptr);
	if(a<__GMP_MEM_NPOOL__)
	{
		EGmemSlabPoolFree(ptr);
		return;
	}
	else
	{
		EGfree(ptr);
		return;
	}
}

/* ========================================================================= */
/*void EGlpNumStart(void) __attribute__ ((constructor));*/
void EGlpNumStart(void)
{
	int rval=0;
	register int i;
	if(__EGlpNum_setup) return;
	if(EG_LPNUM_MEMSLAB)
	{
		for( i = __GMP_MEM_NPOOL__ ; i-- ; )
		{
			EGmemSlabPoolInit(EGgmpPl+i,_EGgmpPlSz[i],0,0);
			rval = EGmemSlabPoolSetParam(EGgmpPl+i,EG_MSLBP_FREEFREE,0);
			EXIT(rval,"Unknown error");
		}
		mp_set_memory_functions(__EGgmp_malloc, __EGgmp_realloc, __EGgmp_free);
	}

	mpf_set_default_prec (EGLPNUM_PRECISION);
	mpz_init (__zeroLpNum_mpz__);
	mpz_init (__oneLpNum_mpz__);
	mpz_init (__MaxLpNum_mpz__);
	mpz_init (__MinLpNum_mpz__);
	mpz_set_ui (__zeroLpNum_mpz__, (unsigned long int)0);
	mpz_set_ui (__oneLpNum_mpz__, (unsigned long int)1);
	mpq_init (__MaxLpNum_mpq__);
	mpq_init (__MinLpNum_mpq__);
	mpf_init (__MaxLpNum_mpf__);
	mpf_init (__MinLpNum_mpf__);
	mpf_init (__zeroLpNum_mpf__);
	mpf_init (__oneLpNum_mpf__);
	mpf_set_ui(__MaxLpNum_mpf__,1UL);
	mpf_set_si(__MinLpNum_mpf__,-1L);
	mpf_mul_2exp(__MaxLpNum_mpf__,__MaxLpNum_mpf__,4096);
	/*mpf_mul_2exp(__MaxLpNum_mpf__,__MaxLpNum_mpf__,ULONG_MAX);*/
	mpf_mul_2exp(__MinLpNum_mpf__,__MinLpNum_mpf__,4096);
	/*mpf_mul_2exp(__MinLpNum_mpf__,__MinLpNum_mpf__,ULONG_MAX);*/
	mpq_set_f(__MaxLpNum_mpq__,__MaxLpNum_mpf__);
	mpq_set_f(__MinLpNum_mpq__,__MinLpNum_mpf__);
	mpz_set_f(__MaxLpNum_mpz__,__MaxLpNum_mpf__);
	mpz_set_f(__MinLpNum_mpz__,__MinLpNum_mpf__);
	mpf_set_ui (__oneLpNum_mpf__, (unsigned long int)1);
	mpf_set_ui (__zeroLpNum_mpf__, (unsigned long int)0);
	mpf_init_set_ui (mpf_eps, (unsigned long int)1);
	mpf_div_2exp (mpf_eps, mpf_eps,  (unsigned long int)(EGLPNUM_PRECISION - 1));
	mpq_init (__zeroLpNum_mpq__);
	mpq_init (__oneLpNum_mpq__);
	mpq_set_ui (__oneLpNum_mpq__, (unsigned long int)1, (unsigned long int)1);
	mpq_set_ui (__zeroLpNum_mpq__, (unsigned long int)0, (unsigned long int)1);
	__EGlpNum_setup=1;
}

/* ========================================================================= */
void EGlpNumSetPrecision (const unsigned prec)
{
	EGLPNUM_PRECISION = prec;
	mpf_set_default_prec (EGLPNUM_PRECISION);
	mpf_clear (mpf_eps);
	mpf_init_set_ui (mpf_eps, (unsigned long int)1);
	mpf_div_2exp (mpf_eps, mpf_eps,  (unsigned long int)(EGLPNUM_PRECISION - 1));
}

/* ========================================================================= */
/*void EGlpNumExit(void) __attribute__ ((destructor));*/
void EGlpNumClear(void)
{
	#if __GMP_MEM_STATS__
	const char mc[5][3] = {"b ","Kb","Mb","Gb","Tb"};
	#endif
	int i;
	if(!__EGlpNum_setup) return;
	mpf_clear (__zeroLpNum_mpf__);
	mpf_clear (__oneLpNum_mpf__);
	mpf_clear (__MaxLpNum_mpf__);
	mpf_clear (__MinLpNum_mpf__);
	mpf_clear (mpf_eps);
	mpq_clear (__zeroLpNum_mpq__);
	mpq_clear (__oneLpNum_mpq__);
	mpq_clear (__MinLpNum_mpq__);
	mpq_clear (__MaxLpNum_mpq__);
	mpz_clear (__zeroLpNum_mpz__);
	mpz_clear (__oneLpNum_mpz__);
	mpz_clear (__MaxLpNum_mpz__);
	mpz_clear (__MinLpNum_mpz__);
	if(EG_LPNUM_MEMSLAB)
	{
		mp_set_memory_functions(0, 0, 0);
		for(i = __GMP_MEM_NPOOL__ ; i-- ; )
		{
			EGmemSlabPoolClear(EGgmpPl+i);
		}
		#if __GMP_MEM_STATS__
		QSlog("GMP alloc statistics:");
		for( i = 0 ; i < __GMP_MEM_NPOOL__ ; i++)
		{
			if(__maxmem > 1024*1024) __maxmem= (__maxmem+1023)/1024;
			else break;
		}
		QSlog("\tmaximum memory allocated      : %8.3lf %s",
								((double)__maxmem)/1024, mc[i+1]);
		QSlog("\tmalloc calls                  : %11zd",__nallocs);
		QSlog("\trealloc calls                 : %11zd",__nrllocs);
		QSlog("\tsmall size allocs-reallocs:");
		for( i = 0 ; i < __GMP_MEM_NPOOL__ ; i++)
		{
			if(__alloc_sz[i] || __rlloc_sz[i])
				QSlog("\t%4d %11zd (%5.2lf%%) %11zd (%5.2lf%%)", 
										_EGgmpPlSz[i], __alloc_sz[i], 
										100.0*((double)__alloc_sz[i])/(__nallocs+__nrllocs),
										__rlloc_sz[i], 
										100.0*((double)__rlloc_sz[i])/(__nrllocs+__nallocs));
		}
		if(__nalarge)
			QSlog("\tlarge size allocs (cals/avg)  : %11zd %11zd",
									__nalarge,__laaverage/__nalarge);
		if(__nrlarge)
			QSlog("\tlarge size reallocs (cals/avg): %11zd %11zd",
									__nrlarge,__lraverage/__nalarge);
		#endif
		QSlog("Disabling EG-GMP mempool");
	}
	__EGlpNum_setup=0;
}

/* ========================================================================= */
void mpq_EGlpNumSet_mpf (mpq_t var,
												 mpf_t flt)
{
	/* local variables */
	unsigned long int __lsgn = mpf_cmp_ui (flt, (unsigned long int)0) < 0 ? (unsigned long int)1 : (unsigned long int)0;
	mpz_t __utmp,
	  __z[7],
	  max_den;
	long int __lexp = 0;
	int i;
	unsigned long int uexp;
	mpf_t __cvl,__lpnum__;
	mpf_init(__lpnum__);
	/* check if the given number is zero, if so, set to zero var and return */
	if (mpf_cmp_ui (flt, (unsigned long int)0) == 0)
	{
		mpq_set_ui (var, (unsigned long int)0, (unsigned long int)1);
		return;
	}
	/* if not, then we have some work to do */
	/* now we initialize the internal numbers */
	mpf_init (__cvl);
	mpf_abs (__cvl, flt);
	mpz_init_set_ui (__utmp, (unsigned long int)0);
	for (i = 7; i--;)
		mpz_init_set_ui (__z[i], (unsigned long int)0);
	mpz_set_ui (__z[0], (unsigned long int)1);
	mpz_set_ui (__z[4], (unsigned long int)1);
	/* max_den is the maximum denominator that we want to see, this number should
	 * be sligtly larger than the square root of 2^EGLPNUM_PRECISION */
	mpz_init_set_ui (max_den, (unsigned long int)1);
	mpz_mul_2exp (max_den, max_den, EGLPNUM_PRECISION >> 1);
	/* first we compute the exponent stored in the limbs */
	/* now we compute the 2^n part needed to set this number between 0.5 and 1 */
	mpf_get_d_2exp(&__lexp,__cvl);
	if (__lexp < 0)
	{
		uexp = (unsigned long int)(-__lexp);
		mpf_mul_2exp (__cvl, __cvl, (unsigned long int) uexp);
	}
	else
	{
		uexp = (unsigned long int)(__lexp);
		mpf_div_2exp (__cvl, __cvl, (unsigned long int) uexp);
	}
	/* now we loop until the next t's is more than mpf_eps */
	/* the formula is 
	 * p_i = t_i*p_{i-1} + p_{i-2}, and 
	 * q_i = t_i*q_{i-1} + q_{i-2} 
	 * note that |x-p_i/q_i|<1/q_i^2
	 * for us t_i = __utmp, and the current number is either [0,1,2] in the __z
	 * array, we use those popsitions ciclicly, and use the four position as a
	 * temporary number, __z+4 is used to store q's, at the beginning i = 1. */
	while (1)
	{
		if (mpf_cmp (__cvl, mpf_eps) < 0 || (mpz_cmp (__z[4], max_den) > 0))
		{
			mpz_set (mpq_denref (var), __z[4]);
			mpz_set (mpq_numref (var), __z[1]);
			break;
		}
		/* first run */
		mpf_ui_div (__cvl, (unsigned long int)1, __cvl);
		mpz_set_f (__utmp, __cvl);
		mpf_set_z (__lpnum__, __utmp);
		mpf_sub (__cvl, __cvl, __lpnum__);
		mpz_set (__z[6], __utmp);
		mpz_set (__z[2], __z[0]);
		mpz_addmul (__z[2], __z[1], __z[6]);
		mpz_set (__z[5], __z[3]);
		mpz_addmul (__z[5], __z[4], __z[6]);
		if (mpf_cmp (__cvl, mpf_eps) < 0 || (mpz_cmp (__z[5], max_den) > 0))
		{
			mpz_set (mpq_denref (var), __z[5]);
			mpz_set (mpq_numref (var), __z[2]);
			break;
		}
		/* second run */
		mpf_ui_div (__cvl, (unsigned long int)1, __cvl);
		mpz_set_f (__utmp, __cvl);
		mpf_set_z (__lpnum__, __utmp);
		mpf_sub (__cvl, __cvl, __lpnum__);
		mpz_set (__z[6], __utmp);
		mpz_set (__z[0], __z[1]);
		mpz_addmul (__z[0], __z[2], __z[6]);
		mpz_set (__z[3], __z[4]);
		mpz_addmul (__z[3], __z[5], __z[6]);
		if (mpf_cmp (__cvl, mpf_eps) < 0 || (mpz_cmp (__z[3], max_den) > 0))
		{
			mpz_set (mpq_denref (var), __z[3]);
			mpz_set (mpq_numref (var), __z[0]);
			break;
		}
		/* third run */
		mpf_ui_div (__cvl, (unsigned long int)1, __cvl);
		mpz_set_f (__utmp, __cvl);
		mpf_set_z (__lpnum__, __utmp);
		mpf_sub (__cvl, __cvl, __lpnum__);
		mpz_set (__z[6], __utmp);
		mpz_set (__z[1], __z[2]);
		mpz_addmul (__z[1], __z[0], __z[6]);
		mpz_set (__z[4], __z[5]);
		mpz_addmul (__z[4], __z[3], __z[6]);
	}

	/* ending */
	mpq_canonicalize (var);
	if (__lsgn)
		mpq_neg (var, var);
	if (__lexp > 0)
		mpq_mul_2exp (var, var, (unsigned long int) __lexp);
	if (__lexp < 0)
		mpq_div_2exp (var, var, (unsigned long int) (-__lexp));
	for (i = 7; i--;)
		mpz_clear (__z[i]);
	mpf_clear (__cvl);
	mpz_clear (max_den);
	mpz_clear (__utmp);
	mpf_clear(__lpnum__);
	return;
}
/* ========================================================================= */
/** @brief verbosity of continued fraction conversion method */
#ifndef MPQ_VERBOSE_CNT_FRAC
#define MPQ_VERBOSE_CNT_FRAC 1000
#endif
/* ========================================================================= */
void mpq_EGlpNumSet (mpq_t var,
										 double const dbl)
{
	/* local variables */
	double __dbl = dbl;
	unsigned __lsgn = __dbl > 0.0 ? 0U: 1U;
	unsigned long __utmp = 0;
	int __lexp = 0;
	double __cvl = __dbl = fabs (__dbl);
	/* we use the first three numbers for p, and the last three numbers for q */
	/* first check that the dbl is not zero */
	if (__dbl < 1e-151)
	{
		mpq_set_ui (var, (unsigned long int)0, (unsigned long int)1);
		__lsgn = 0;
	}
	else if (__dbl > 1e151)
		mpq_set_d (var, __dbl);
	else
	{
		/* now we initialize the integer numbers */
		mpz_t __z[7];
		for (__utmp = 7; __utmp--;)
			mpz_init (__z[__utmp]);
		mpz_set_ui (__z[0], (unsigned long int)1);
		mpz_set_ui (__z[4], (unsigned long int)1);
		/* now we compute the 2^n part needed to set this number between 0 and 1 */
#define __HI_EXP(x,e,v,lv) {if( x >=v ){ e = e + lv; x /= v;}}
		if (__cvl > 1)
		{
			__HI_EXP (__cvl, __lexp,
								115792089237316195423570985008687907853269984665640564039457584007913129639936.0,
								256);
			__HI_EXP (__cvl, __lexp, 340282366920938463463374607431768211456.0, 128);
			__HI_EXP (__cvl, __lexp, 18446744073709551616.0, 64);
			__HI_EXP (__cvl, __lexp, 4294967296.0, 32);
			__HI_EXP (__cvl, __lexp, 65536.0, 16);
			__HI_EXP (__cvl, __lexp, 256.0, 8);
			__HI_EXP (__cvl, __lexp, 16.0, 4);
			__HI_EXP (__cvl, __lexp, 4.0, 2);
			__HI_EXP (__cvl, __lexp, 2.0, 1);
#undef __HI_EXP
		}
		else if (__cvl < 0.5)
		{
#define __LO_EXP(x,e,v,lv) {if( x < 1/v ) { e = e - lv; x *= v;}}
			__LO_EXP (__cvl, __lexp,
								115792089237316195423570985008687907853269984665640564039457584007913129639936.0,
								256);
			__LO_EXP (__cvl, __lexp, 340282366920938463463374607431768211456.0, 128);
			__LO_EXP (__cvl, __lexp, 18446744073709551616.0, 64);
			__LO_EXP (__cvl, __lexp, 4294967296.0, 32);
			__LO_EXP (__cvl, __lexp, 65536.0, 16);
			__LO_EXP (__cvl, __lexp, 256.0, 8);
			__LO_EXP (__cvl, __lexp, 16.0, 4);
			__LO_EXP (__cvl, __lexp, 4.0, 2);
			__LO_EXP (__cvl, __lexp, 2.0, 1);
#undef __LO_EXP
		}
		/* now we loop until the next t's is more than EGLPNUM_MINEPS */
		/* the formula is 
		 * p_i = t_i*p_{i-1} + p_{i-2}, and 
		 * q_i = t_i*q_{i-1} + q_{i-2} 
		 * note that |x-p_i/q_i|<1/q_i^2
		 * for us t_i = __utmp, and the current number is either [0,1,2] in the __z
		 * array, we use those popsitions ciclicly, and use the four position as a
		 * temporary number, __z+4 is used to store q's, at the beginning i = 1. */
		while (1)
		{
			if (__cvl < EGLPNUM_MINEPS || (mpz_cmp_ui (__z[4], (unsigned long int)0xfffffff) > 0))
			{
				mpz_set (mpq_denref (var), __z[4]);
				mpz_set (mpq_numref (var), __z[1]);
				break;
			}
			MESSAGE(MPQ_VERBOSE_CNT_FRAC,"cur approximate %ld/%ld, error %10.7lg", 
							mpz_get_si(__z[1]), mpz_get_si(__z[4]), __cvl);
			/* first run */
			__cvl = 1 / __cvl;
			__dbl = floor(__cvl);
			__utmp = (unsigned long)__dbl;
			__cvl -= __dbl;
			mpz_set_ui (__z[6], __utmp);
			mpz_set (__z[2], __z[0]);
			mpz_addmul (__z[2], __z[1], __z[6]);
			mpz_set (__z[5], __z[3]);
			mpz_addmul (__z[5], __z[4], __z[6]);
			if (__cvl < EGLPNUM_MINEPS || (mpz_cmp_ui (__z[5], (unsigned long int)0xfffffff) > 0))
			{
				mpz_set (mpq_denref (var), __z[5]);
				mpz_set (mpq_numref (var), __z[2]);
				break;
			}
			MESSAGE(MPQ_VERBOSE_CNT_FRAC,"cur approximate %ld/%ld, error %10.7lg", 
							mpz_get_si(__z[2]), mpz_get_si(__z[5]), __cvl);
			/* second run */
			__cvl = 1 / __cvl;
			__dbl = floor(__cvl);
			__utmp = (unsigned long)(__dbl);
			__cvl -= __dbl;
			mpz_set_ui (__z[6], __utmp);
			mpz_set (__z[0], __z[1]);
			mpz_addmul (__z[0], __z[2], __z[6]);
			mpz_set (__z[3], __z[4]);
			mpz_addmul (__z[3], __z[5], __z[6]);
			if (__cvl < EGLPNUM_MINEPS || (mpz_cmp_ui (__z[3], (unsigned long int)0xfffffff) > 0))
			{
				mpz_set (mpq_denref (var), __z[3]);
				mpz_set (mpq_numref (var), __z[0]);
				break;
			}
			MESSAGE(MPQ_VERBOSE_CNT_FRAC,"cur approximate %ld/%ld, error %10.7lg", 
							mpz_get_si(__z[0]), mpz_get_si(__z[3]), __cvl);
			/* third run */
			__cvl = 1 / __cvl;
			__dbl = floor(__cvl);
			__utmp = (unsigned long)(__dbl);
			__cvl -= __dbl;
			mpz_set_ui (__z[6], __utmp);
			mpz_set (__z[1], __z[2]);
			mpz_addmul (__z[1], __z[0], __z[6]);
			mpz_set (__z[4], __z[5]);
			mpz_addmul (__z[4], __z[3], __z[6]);
		}
		for (__utmp = 7; __utmp--;)
			mpz_clear (__z[__utmp]);
	}
	/* ending */
	mpq_canonicalize (var);
	if (__lsgn)
		mpq_neg (var, var);
	if (__lexp > 0)
		mpq_mul_2exp (var, var, (unsigned long int) __lexp);
	if (__lexp < 0)
		mpq_div_2exp (var, var, (unsigned long int) (-__lexp));
	return;
}

/* ========================================================================= */
int mpz_EGlpNumReadStr (mpz_t var,
												char const *str)
{
	/* local variables */
	char unsigned a_sgn = 1;
	char unsigned sgn = 0;
	char c = 0;
	int n_char = 0;
	/* now we read the string */
	c = str[n_char];
	mpz_set_ui (var, (unsigned long int)0);
	while ((('0' <= c) && (c <= '9')) ||	/* allow to read digits */
				 (a_sgn && (c == '+' || c == '-')) /* allow sign for exponent */ )
	{
		switch (c)
		{
		case '0':
		case '1':
		case '2':
		case '3':
		case '4':
		case '5':
		case '6':
		case '7':
		case '8':
		case '9':
			mpz_mul_ui (var, var, (unsigned long int)10);
			mpz_add_ui (var, var,  (unsigned long int)(c - '0'));
			a_sgn = 0;
			break;
		case '-':
			sgn = 1;
		case '+':
			a_sgn = 0;
			break;
		}
		/* advance the reading character */
		c = str[++n_char];
	}
	if (sgn)
		mpz_neg (var, var);
	return n_char;
}

/* ========================================================================= */
int mpq_EGlpNumReadStrXc (mpq_t var,
													char const *str)
{
	/* local variables */
	char unsigned a_dot = 1,
	  a_exp = 0,
	  a_exp_sgn = 0,
	  a_sgn = 1,
	  a_div = 1;
	char c = 0;
	int l_exp = 0,
	  sgn = 0,
	  exp_sgn = 0;
	int n_char = 0,
	  n_dig = 0,
	  cn = 0;
	mpq_t den[2];
	mpq_init (den[0]);
	mpq_init (den[1]);
	mpq_set_ui (den[1], (unsigned long int)1, (unsigned long int)1);
	mpq_set_ui (den[0], (unsigned long int)0, (unsigned long int)1);

	/* now we read the string */
	c = str[n_char];
	while ((('0' <= c) && (c <= '9')) ||	/* allow to read digits */
				 (a_dot && (c == '.')) ||	/* allow to read a dot point */
				 (a_exp && (c == 'e' || c == 'E')) ||	/* allow an exponent marker */
				 (a_sgn && (c == '+' || c == '-')) ||	/* allow a number sign */
				 (a_div && c == '/') ||	/* allow the division sign */
				 (a_exp_sgn && (c == '+' || c == '-')) /* allow sign for exponent */ )
	{
		switch (c)
		{
		case '0':
		case '1':
		case '2':
		case '3':
		case '4':
		case '5':
		case '6':
		case '7':
		case '8':
		case '9':
			/* if we haven't read the exponent then the digits bellongs to the mantisa
			 * */
			if (a_exp || n_dig == 0)
			{
				if (!a_dot)
					mpz_mul_ui (mpq_denref (den[cn]), mpq_denref (den[cn]), (unsigned long int)10);
				mpz_mul_ui (mpq_numref (den[cn]), mpq_numref (den[cn]), (unsigned long int)10);
				mpz_add_ui (mpq_numref (den[cn]), mpq_numref (den[cn]),
										 (unsigned long int)(c - '0'));
				n_dig++;
				a_exp = 1;
			}
			/* otherwise, if we have read the exponent, the digits should go to the
			 * exponent */
			else
			{
				l_exp = 10 * l_exp + c - '0';
				a_exp_sgn = 0;
			}
			a_sgn = 0;
			break;
		case '.':
			a_sgn = 0;
			a_dot = 0;
			a_sgn = 0;
			break;
		case '-':
			if (a_sgn)
				sgn = 1;
			else
				exp_sgn = 1;
		case '+':
			if (a_sgn)
				a_sgn = 0;
			if (a_exp_sgn)
				a_exp_sgn = 0;
			break;
		case 'e':
		case 'E':
			a_sgn = 0;
			a_exp = 0;
			a_exp_sgn = 1;
			break;
		case '/':
			if (exp_sgn)
				l_exp = -l_exp;
			if (l_exp > 0)
				while (l_exp--)
					mpz_mul_ui (mpq_numref (den[0]), mpq_numref (den[0]), (unsigned long int)10);
			else if (l_exp < 0)
			{
				l_exp = -l_exp;
				while (l_exp--)
					mpz_mul_ui (mpq_denref (den[0]), mpq_denref (den[0]), (unsigned long int)10);
			}
			if (sgn)
				mpz_neg (mpq_numref (den[0]), mpq_numref (den[0]));
			mpq_canonicalize (den[0]);
			mpq_set_ui (den[1], (unsigned long int)0, (unsigned long int)1);
			sgn = 0;
			exp_sgn = 0;
			l_exp = 0;
			a_div = 0;
			n_dig = 0;
			a_dot = 1;
			a_exp = 0;
			a_exp_sgn = 0;
			a_sgn = 1;
			cn = 1;
			break;
		}
		/* advance the reading character */
		c = str[++n_char];
	}
	if (n_char)
	{
		/* now expand the exponent of the denominator */
		if (exp_sgn)
			l_exp = -l_exp;
		if (l_exp > 0)
			while (l_exp--)
				mpz_mul_ui (mpq_numref (den[cn]), mpq_numref (den[cn]), (unsigned long int)10);
		else if (l_exp < 0)
		{
			l_exp = -l_exp;
			while (l_exp--)
				mpz_mul_ui (mpq_denref (den[cn]), mpq_denref (den[cn]), (unsigned long int)10);
		}
		/* check the sign of the whole number */
		if (sgn)
			mpz_neg (mpq_numref (den[cn]), mpq_numref (den[cn]));
		/* ending */
		mpq_canonicalize (den[0]);
		mpq_canonicalize (den[1]);
		mpq_div (var, den[0], den[1]);
	}
	mpq_clear (den[0]);
	mpq_clear (den[1]);
	return n_char;
}

/* ========================================================================= */
void uint32_EGutilPermSort (const size_t sz,
										 int *const perm,
										 const uint32_t * const elem)
{
	size_t i,
	  j;
	int temp;
	uint32_t t;
	if (sz <= 1)
		return;

	EGswap (perm[0], perm[(sz - 1) / 2], temp);
	i = 0;
	j = sz;
	(t= elem[perm[0]]);
	for (;;)
	{
		do
			i++;
		while (i < sz && (elem[perm[i]]< t));
		do
			j--;
		while (j && (t< elem[perm[j]]));
		if (j < i)
			break;
		EGswap (perm[i], perm[j], temp);
	}
	EGswap (perm[0], perm[j], temp);
	uint32_EGutilPermSort (j, perm, elem);
	uint32_EGutilPermSort (sz - i, perm + i, elem);
}

/* ========================================================================= */
void uint32_EGutilPermSort2 (const size_t sz,
										 int *const perm,
										 const uint32_t * const elem)
{
	size_t i,
	  j;
	int temp;
	uint32_t t;
	if (sz <= 1)
		return;

	EGswap (perm[0], perm[(sz - 1) / 2], temp);
	i = 0;
	j = sz;
	(t= elem[perm[0]]);
	for (;;)
	{
		do
			i++;
		while (i < sz && (t< elem[perm[i]]));
		do
			j--;
		while ((elem[perm[j]]< t));
		if (j < i)
			break;
		EGswap (perm[i], perm[j], temp);
	}
	EGswap (perm[0], perm[j], temp);
	uint32_EGutilPermSort2 (j, perm, elem);
	uint32_EGutilPermSort2 (sz - i, perm + i, elem);
}

/* ========================================================================= */
/** @} */
