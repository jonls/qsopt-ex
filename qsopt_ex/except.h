/****************************************************************************/
/*                                                                          */
/*  This file is part of QSopt_ex.                                          */
/*                                                                          */
/*  (c) Copyright 2006 by David Applegate, William Cook, Sanjeeb Dash,      */
/*  and Daniel Espinoza                                                     */
/*                                                                          */
/*  Sanjeeb Dash ownership of copyright in QSopt_ex is derived from his     */
/*  copyright in QSopt.                                                     */
/*                                                                          */
/*  This code may be used under the terms of the GNU General Public License */
/*  (Version 2.1 or later) as published by the Free Software Foundation.    */
/*                                                                          */
/*  Alternatively, use is granted for research purposes only.               */
/*                                                                          */
/*  It is your choice of which of these two licenses you are operating      */
/*  under.                                                                  */
/*                                                                          */
/*  We make no guarantees about the correctness or usefulness of this code. */
/*                                                                          */
/****************************************************************************/

/* RCSINFO $Id: except.h,v 1.3 2003/11/05 17:02:10 meven Exp $ */
#ifndef ILL_except
#define ILL_except

#include "allocrus.h"
#include "trace.h"

/* Version 2.4 and later of GCC define a magical variable `__PRETTY_FUNCTION__'
   which contains the name of the function currently being defined.
#  define __DEV_FUNCTION__     __PRETTY_FUNCTION__
   This is broken in G++ before version 2.6.
   C9x has a similar variable called __func__, but prefer the GCC one since
   it demangles C++ function names.  */
# ifdef __GNUC__
#  if __GNUC__ > 2 || (__GNUC__ == 2 \
                       && __GNUC_MINOR__ >= (defined __cplusplus ? 6 : 4))
#   define __DEV_FUNCTION__    __PRETTY_FUNCTION__
#  else
#   define __DEV_FUNCTION__    ((__const char *) 0)
#  endif
# else
#  if defined __STDC_VERSION__ && __STDC_VERSION__ >= 199901L
#   define __DEV_FUNCTION__    __func__
#  else
#   define __DEV_FUNCTION__    ((const char *) 0)
#  endif
# endif


/* put debugger breakpoint here */
extern void ILL_report (
	const char *msg,
	const char *fct,
	const char *file,
	unsigned int line,
	int with_source_info);

/* printed message looks as follows
 *
 * with_source_info == 0:   <msg> + "\n"
 *
 * with_source_info == 1:   if (fct != NULL) 
 *                              <msg> + " in function <fct>\n";
 *                          else 
 *                              <msg> + " in file <file> line <line>\n"; 
 */

#define ILL_GENERAL_ERROR -1
#define ILL_NO_MEMORY 2
#define ILL_NULL_PTR  3

#define ILL_REPORT(msg,with) \
    ILL_report(msg, __DEV_FUNCTION__, __FILE__, __LINE__, with)
#ifdef NDEBUG
#define ILL_REPRT(msg) \
    ILL_report(msg, __DEV_FUNCTION__, __FILE__, __LINE__, 0)
#else
#define ILL_REPRT(msg) \
    ILL_report(msg, __DEV_FUNCTION__, __FILE__, __LINE__, 1)
#endif

#define ILL_RESULT(expr, msg) \
{								\
	if (TRACE > 0) { ILL_RETURN(expr, msg); } \
    return expr; \
}

#define ILL_RETURN_PTR(ptr, msg) \
    { void *ILL_RETURN_p = ptr; \
      if (ILL_RETURN_p == NULL) {  \
        if (TRACE > 0) ILL_REPRT(msg); \
      } \
      return ILL_RETURN_p;  \
    }

#ifdef NDEBUG
#define ILL_RETURN(expr, msg)           \
{										\
	if (expr != 0) {					\
			if (TRACE > 0) ILL_REPRT(msg); \
    }									\
    return expr;						\
}

#else
#define ILL_RETURN(expr, msg)           \
	{									\
		if (expr != 0) {                \
			ILL_REPRT(msg);             \
		}                               \
		ILL_IFTRACE("%s: returning %d\n", __DEV_FUNCTION__, expr); \
		return expr;					\
	}
#endif

#define ILL_CHECKnull(expr, msg) \
    { if ((expr) == NULL)  {  \
         ILL_REPRT(msg); \
         rval = ILL_NULL_PTR; \
         goto CLEANUP;  \
      } }

#define ILL_FAILtrue(expr, msg) \
    { if (expr)  {  \
         ILL_REPRT(msg); \
         rval = ILL_GENERAL_ERROR; \
         goto CLEANUP;  \
      } }

#define ILL_FAILtrue_no_rval(expr, msg) \
    { if (expr)  {  \
         ILL_REPRT(msg); \
         goto CLEANUP;  \
      } }


#define ILL_FAILfalse(expr, msg)  ILL_FAILtrue(!(expr), msg)
#define ILL_FAILfalse_no_rval(expr, msg)  ILL_FAILtrue_no_rval(!(expr), msg)

#define ILL_ERROR(rval, msg)      {									\
									fprintf(stderr, "%s\n", msg);	\
									rval = 1; goto CLEANUP;			\
								  }
#define ILL_CLEANUP_IF(rval)      { if ((rval) != 0) { goto CLEANUP; } }
#define ILL_CLEANUP               goto CLEANUP

#define ILL_SAFE_MALLOC(lhs, n, type) \
    { lhs = ILL_UTIL_SAFE_MALLOC(n, type, lhs); \
      if (lhs == NULL)  {  \
         ILL_REPRT("Out of memory"); \
         rval = ILL_NO_MEMORY; \
         goto CLEANUP;  \
      }}

#define ILL_SAFE_MALLOC_no_rval(lhs, n, type) \
    { lhs = ILL_UTIL_SAFE_MALLOC(n, type, lhs); \
      if (lhs == NULL)  {  \
         ILL_REPRT("Out of memory"); \
         goto CLEANUP;  \
      }}


#define ILL_NEW(ptr, type) ILL_SAFE_MALLOC(ptr, 1, type)
#define ILL_NEW_no_rval(ptr, type) ILL_SAFE_MALLOC_no_rval(ptr, 1, type)

/* we define debugging verbosity for singular basis */
extern int __QS_SB_VERB;
#endif
