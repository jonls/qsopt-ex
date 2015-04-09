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

/* RCSINFO $Id: trace.h,v 1.2 2003/11/05 16:47:22 meven Exp $ */
#ifndef ILL_trace_h
#define ILL_trace_h

#include "logging-private.h"

/* users of these macros must declare a static int TRACE variable */
#ifndef NDEBUG
#define ILL_IFTRACE        if (TRACE) QSlog
#define ILL_IFTRACE2       if (TRACE > 1) QSlog
#define ILL_IFDOTRACE      if (TRACE)
#else
/* the optimizer will take care of this */
#define ILL_IFTRACE        if (0) QSlog
#define ILL_IFTRACE        if (0) QSlog
#define ILL_IFDOTRACE      if (0)
#endif

#endif
