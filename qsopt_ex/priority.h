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

#ifndef __PRIORITY_H__
#define __PRIORITY_H__
#include "dheaps_i.h"
/****************************************************************************/
/*                                                                          */
/*                             priority.c                                   */
/*                                                                          */
/****************************************************************************/

typedef struct ILLpriority
{
	ILLdheap heap;
	union ILLpri_data
	{
		void *data;
		int next;
	}
	 *pri_info;
	int space;
	int freelist;
}
ILLpriority;

void ILLutil_priority_free (
	ILLpriority * pri),
  ILLutil_priority_delete (
	ILLpriority * pri,
	int handle),
  ILLutil_priority_changekey (
	ILLpriority * pri,
	int handle,
	EGlpNum_t * newkey),
  ILLutil_priority_findmin (
	ILLpriority * pri,
	EGlpNum_t * keyval,
	void **en),
  ILLutil_priority_deletemin (
	ILLpriority * pri,
	EGlpNum_t * keyval,
	void **en);

int ILLutil_priority_init (
	ILLpriority * pri,
	int k),
  ILLutil_priority_insert (
	ILLpriority * pri,
	void *data,
	EGlpNum_t * keyval,
	int *handle);



#endif
