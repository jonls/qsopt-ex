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
#ifndef __EG_NUMMACROS_H__
#define __EG_NUMMACROS_H__

#include "eg_macros.h"

/* ========================================================================= */
/** @defgroup EGlpNumMacros General Number Utilities
 * Here we put some utilities common for number.
 * 
 * @par History:
 * Revision 0.0.2
 *  - 2007-10-08
 *  					- Move EGabs, EGswap, EGmin and EGmax to this file
 * */
/** @{*/
/** @file
 * @brief This file provide the user interface and function definitions 
 * for general number utilities.
 * */

/* ========================================================================= */
/** @brief Given tree numbers N1, N2 and Ntmp, swap values of N1 and N2 using
 * Ntmp as a temporal number. The variables should be of some primitive type of
 * C for this macro to work.
 * @param N1 first number.
 * @param N2 second number.
 * @param Ntmp temporal variable.
 * */
#define EGswap(N1,N2,Ntmp) do{\
	Ntmp = N1;\
	N1 = N2;\
	N2 = Ntmp;} while(0)

/* ========================================================================= */
/** @brief given two variables (of the same type, and of some predefined type)
 * return the maximum value among the two of them. */
#define EGmax(a,b) ({\
	const typeof(a) __EGma = (a);\
	const typeof(b) __EGmb = (b);\
	(__EGma > __EGmb ? __EGma : __EGmb);})

/* ========================================================================= */
/** @brief given two variables (of the same type, and of some predefined type)
 * return the minimum value among the two of them. */
#define EGmin(a,b) ({\
	const typeof(a) __EGma = (a);\
	const typeof(b) __EGmb = (b);\
	(__EGma < __EGmb ? __EGma : __EGmb);})

/* ========================================================================= */
/** @brief a general macro to return the absolute value of the given variable
 * @param var variable whose absolute value we want to compute.
 * @return value of the absolute value of the given variable, note that this
 * macro will only work in built-in types, and will use the default comparison
 * for those internal types. */
#define EGabs(var) ({\
		const typeof(var) __EGav = (var);\
		(__EGav < 0) ? -__EGav : __EGav;})

/* ========================================================================= */
/** @}*/
#endif
