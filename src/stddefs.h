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

/* $RCSfile: stddefs.h,v $ $Revision: 1.3 $ $Date: 2003/11/05 16:57:39 $ */
#ifndef ILL_STDDEFS_H
#define ILL_STDDEFS_H

#define SWAP(x,y,temp) {temp = x; x = y; y = temp;}

#define QSMIN(x,y) ((x) < (y) ? (x) : (y))
#define QSMAX(x,y) ((x) > (y) ? (x) : (y))

#define ABS(x) ((x) >= 0 ? (x) : -(x))

#endif /* ILL_STDDEFS_H */
