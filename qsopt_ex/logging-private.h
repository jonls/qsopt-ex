/* logging - Private functions for logging/error reporting
 *
 * Copyright (C) 2015  Jon Lund Steffensen <jonlst@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * */

#ifndef QS_LOGGING_PRIVATE_H__
#define QS_LOGGING_PRIVATE_H__

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>


void QSlogv(const char *format, va_list args)
	__attribute__ ((format(__printf__, 1, 0)));
void QSlog(const char *format, ...)
	__attribute__ ((format(__printf__, 1, 0)));


/** @name Code Location Utility:
 * this are utility macros to print information about where we are.
 * @note For some reason __func__ don't work correctly for sun's cc in inline
 * functions, so we don't use in SUN architecture */
/*  @{  */
#define __EG_PRINTLOCF__  QSlog(", in %s (%s:%d)",__func__,__FILE__,__LINE__)
#define __EG_PRINTLOC__      __EG_PRINTLOCF__
#define __EG_PRINTLOC2F__ QSlog("in %s (%s:%d)",__func__,__FILE__,__LINE__)
#define __EG_PRINTLOC2__     __EG_PRINTLOC2F__
/*  @}  */

/** @brief if the given display level is over the given treshold, display the
 * given message and location information */
#define IFMESSAGE(__display,...) do{\
	if((__display)>0){\
		QSlog(__VA_ARGS__);										\
		__EG_PRINTLOC__;}}while(0)

/** @brief if return value is non-zero warn on screen */
#define WARNIF(__L__) do{\
	const int __Wval__= (__L__);\
	if(__Wval__){QSlog("WARNING: In %s (%s:%d) "#__L__" = %d",__func__,__FILE__,__LINE__,__Wval__);}}while(0)
#if DEBUG>=1

/** @brief This macro is to print error messages and to return with value one
 * from the current function, it also print the file and line where this
 * happend, but the condition is looked only if the debug level is at least __L
 * */
#define EXITL(__L,__A,...) ({\
	if(__L<=DEBUG){\
	if(__A){\
		QSlog(__VA_ARGS__);										\
		__EG_PRINTLOC__;\
		_exit(1);}}})

/** @brief This macro is to print error messages and to return with value one
 * from the current function, it also print the file and line where this
 * happend, but the condition is looked only if the debug level is at least __L
 * */
#define TESTL(__L,__A,...) ({\
	if(__L<=DEBUG){\
	if(__A){\
		QSlog(__VA_ARGS__);										\
		__EG_PRINTLOC__;\
		return 1;}}})

/** @brief this macro check if the value of a pointer is not bellow the first
 * 64Kb, if so it return the given value */
#define PTRTEST(__PTR,__RVAL) {\
	if(__PTR) ADVTESTL(0,((size_t)(__PTR)) < (1U<<16),__RVAL, \
								 "%s=%p is not a valid pointer",\
									#__PTR, (void*)(__PTR));}

/** @brief This macro is to print error messages and jump to the given point
 * in the code, it also print the file and line where this happend */
#define TESTGL(__L,__A,__B,...) ({\
	if((__L)<=DEBUG && (__A)){\
		QSlog("ERROR: " __VA_ARGS__);					\
		__EG_PRINTLOC__;\
		goto __B;}})

/** @brief return value macro, if the value is non-zero, write to standard
 * error the returning code and where this happen */
#define EG_RETURN(__A) {\
	const int __RVAL__ = (__A);\
	if(__RVAL__){\
		QSlog("rval %d",__RVAL__);						\
		__EG_PRINTLOC__;}\
	return __RVAL__;}

/** @brief This macro is to print error messages and jump to the given point
 * in the code, it also print the file and line where this happend */
#define TESTG(__A,__B,...) ({\
	if(__A){\
		QSlog("ERROR: " __VA_ARGS__);					\
		__EG_PRINTLOC__;\
		goto __B;}})

/** @brief This macro is to print error messages and jump to the CLEANUP label
 * in the current function, it also print the file and line where this happend */
#define TESTGD(__A,...) ({\
	if(__A){\
		QSlog("ERROR: " __VA_ARGS__);					\
		__EG_PRINTLOC__;\
		goto CLEANUP;}})

/** @brief This macro is to print error messages and to return with value one
 * from the current function, it also print the file and line where this
 * happend */
#define TEST(__A,...) ({\
	if(__A){\
		QSlog(__VA_ARGS__);										\
		__EG_PRINTLOC__;\
		return 1;}})

/** @brief This macro print messages to the screen when the debug level is as
 * big as the first parameter, if the debug level is zero we eliminate the
 * code and reduce it to the empty instruction. */
#define MESSAGEF(__A,...) ({\
	if(__A <= DEBUG ){\
		QSlog(__VA_ARGS__);										\
		__EG_PRINTLOCF__(__F);}})

/** @brief This macro print messages to the screen when the debug level is as
 * big as the first parameter, if the debug level is zero we eliminate the
 * code and reduce it to the empty instruction. */
#define MESSAGE(__A,...) ({\
	if(__A <= DEBUG ){\
		QSlog(__VA_ARGS__);										\
		__EG_PRINTLOC__;}})

/** @brief This macro print messages to the screen when the verbose level is as
 * big as the first parameter, if the verbose level is zero we eliminate the
 * code and reduce it to the empty instruction. */
#if VERBOSE_LEVEL >= 1
#define OUTPUT(__A,...) ({\
	if(__A <= VERBOSE_LEVEL ){\
		QSlog(__VA_ARGS__);}})
#else
#define OUTPUT(__A,...) ;
#endif

/** @brief This macro print messages to the screen when the condition __A is
 * true .if the debug level is one we don't print any warning message. if
 * the debug level is zero we eliminate the code and reduce it to the empty
 * instruction. */
#define WARNINGL(__L,__A,...) ({\
	if((__A)&&(DEBUG>=__L)){\
		QSlog("WARNING: " __VA_ARGS__);\
		__EG_PRINTLOC__;}})

#else
#define TESTL(__L,__A,...) ;
#define EXITL(__L,__A,...) ;
#define TEST(__A,...) ;
#define TESTG(__A,__B,...) ;
#define TESTGL(__A,__B,...) ;
#define MESSAGE(__A,...) ;
#define MESSAGEF(__A,__F,...) ;
#define WARNINGL(__L,__A,...) ;
#define EG_RETURN(__A) do{return (__A);}while(0)
#define PTRTEST(__PTR,__RVAL) ;
#endif

/** @brief This macro is to print error messages and to return with value one
 * from the current function, it also print the file and line where this
 * happend */
#define FTEST(__A,...) ({\
	if(__A){\
		QSlog(__VA_ARGS__);										\
		__EG_PRINTLOC__;\
		return 1;}})

/** @brief This macro is to print error messages and to return with value one
 * from the current function, it also print the file and line where this
 * happend */
#define FTESTG(__A,__B,...) ({\
	if(__A){\
		QSlog("ERROR: " __VA_ARGS__);					\
		__EG_PRINTLOC__;\
		goto __B;}})

/** @brief This macro print messages to the screen when the condition __A is
 * true. */
#define WARNING(__A,...) ({									\
	if(__A){																	\
		QSlog("WARNING: " __VA_ARGS__);		\
		__EG_PRINTLOC__;}})

/** @brief this macro test if a value is non zero, if it is it print where is
 * it and exit EXIT_FAILURE. The idea is to use it to check return values of functions,
 * and the calling function can't return a status, and then we are forced to
 * exit. */
#define EXITRVAL(__A) ({\
	if(__A){\
		__EG_PRINTLOC2__;\
		exit(EXIT_FAILURE);}})

/** @brief This macro is to print error messages and exit the program with
 * code one from the current function, it also print the file and line where
 * this happend */
#define EXIT(__A,...) ({										\
	if(__A){																	\
		QSlog("EXIT: " __VA_ARGS__);			\
		__EG_PRINTLOC__;\
		exit(EXIT_FAILURE);}})

/** @brief this macro test if a value is non zero, if it is it print where is
 * it and return __B. The idea is to use it to check return values of functions
 * */
#define ADVCHECKRVAL(__A,__B) ({\
	if(__A){\
		__EG_PRINTLOC2__;\
		return __B;}})

/** @brief This macro test a condition '__A' when the debug level used at
 * compile time is at least '__L'. If the condition is true, it print the
 * message and return the '__RVAL' value. */
#define ADVTESTL(__L,__A,__RVAL,...) ({\
	if((DEBUG>=__L)&&(__A)){\
		QSlog(__VA_ARGS__);										\
		__EG_PRINTLOC__;\
		return __RVAL;}})

/** @brief this macro test if a value is non zero, if it is it print where is
 * it and return 1. The idea is to use it to check return values of functions
 * */
#define CHECKRVAL(__A) ({\
	if(__A){\
		__EG_PRINTLOC2__;\
		return __A;}})

/** @brief, if a non-zero value is given as an argument, check the errno stored
 * in the system, print the related message, and return the non-zero given
 * parameter, otherwise, do nothing.
 * @param __value if non-zero check systems errors, and return this value
 * */
#define TESTERRNOIF(__value) do{\
	if(__value){\
		const int __EGserrno = errno;\
		QSlog("failed with errno %d, %s",__EGserrno, strerror(__EGserrno)); \
		__EG_PRINTLOC2__;\
		return __value;}}while(0)

/** @brief this function, if the input is non zero, print a message of
 * function, file and line and then goto the second parameter */
#define CHECKRVALG(__A,__B) do{if(__A){__EG_PRINTLOC2__;goto __B;}}while(0)

#endif /* ! QS_LOGGING_PRIVATE_H__ */
