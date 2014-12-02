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

#ifndef __ZEIT_H__
#define __ZEIT_H__
/****************************************************************************/
/*                                                                          */
/*                             zeit.c                                       */
/*                                                                          */
/****************************************************************************/

typedef struct ILLutil_timer
{
	double szeit;
	double cum_zeit;
	char name[40];
	int count;
}
ILLutil_timer;


double ILLutil_zeit (
	void),
  ILLutil_real_zeit (
	void),
  ILLutil_stop_timer (
	ILLutil_timer * t,
	int printit),
  ILLutil_total_timer (
	ILLutil_timer * t,
	int printit);


void ILLutil_init_timer (
	ILLutil_timer * t,
	const char *name),
  ILLutil_start_timer (
	ILLutil_timer * t),
  ILLutil_suspend_timer (
	ILLutil_timer * t),
  ILLutil_resume_timer (
	ILLutil_timer * t);



#endif
