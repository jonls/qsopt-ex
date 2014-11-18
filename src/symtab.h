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

/* RCSINFO $Id: symtab.h,v 1.2 2003/11/05 16:47:22 meven Exp $ */
#ifndef ILL_SYMTAB_H
#define ILL_SYMTAB_H

#include <stdlib.h>
#include <stdio.h>

/* we allow 256KB of buffer space i.e. 2^17*/
#define ILL_namebufsize 0x20000U

typedef struct ILLsymbolent
{
	int symbol;
	int index;
	int next;
}
ILLsymbolent;

typedef struct ILLsymboltab
{
	int *hashtable;
	ILLsymbolent *nametable;
	char *namelist;
	int tablesize;
	int strsize;
	int hashspace;
	int name_space;
	int strspace;
	int freedchars;
	int the_hash;
	int the_index;
	int the_prev_index;
	int index_ok;
}
ILLsymboltab;

/* 
 * hashtable[stringhash(entry) % hashspace] either NO_INDEX or some hash number
 * nametable[hash number] = { next:    another index for nametable 
 *                            symbol:  index into namelist where string resides
 *                          }
 * tablesize:  number of entries   (tablesize <= name_space)
 * name_space:  length of nametable and indexlist
 * hashspace:  length of hashtable   nextprime(name_space)
 * strsize:    number of chars used in namelist 
 * strspace:   length of namelist
 * indexlist:  LP col/row indices for the table entries
 * indexlist_ok:  1 if column indices in indexlist are up-to-date, 0 otherwise
 *
 * Deletion of entries affects their ordering in symboltab->nametable. 
 * Strings may move around within symboltab->namelist. 
 */


#define ILL_SYM_NOINDEX (-1)
extern void ILLsymboltab_init (
	ILLsymboltab * h),
  ILLsymboltab_free (
	ILLsymboltab * h),
  ILLsymboltab_size (
	const ILLsymboltab * h,
	int *p_size),
  ILLsymboltab_prt (
	FILE * fd,
	ILLsymboltab * h);

extern int ILLsymboltab_create (
	ILLsymboltab * h,
	int init_size),
  ILLsymboltab_copy (
	ILLsymboltab * src,
	ILLsymboltab * dst),
  ILLsymboltab_register (
	ILLsymboltab * h,
	const char *s,
	int itemindex,
	int *p_index,
	int *p_existed),
  ILLsymboltab_lookup (
	ILLsymboltab * h,
	const char *s,
	int *p_index),
  ILLsymboltab_index_ok (
	ILLsymboltab * h),
  ILLsymboltab_index_reset (
	ILLsymboltab * h,
	int icount,
	char **names),
  ILLsymboltab_getindex (
	ILLsymboltab * h,
	const char *name,
	int *hindex),
  ILLsymboltab_contains (
	ILLsymboltab * h,
	const char *s),
  ILLsymboltab_delete (
	ILLsymboltab * h,
	const char *s),
  ILLsymboltab_uname (
	ILLsymboltab * h,
	char name[ILL_namebufsize],
	const char *try_prefix1,
	const char *try_prefix2);

extern void ILLsymboltab_unique_name (
	ILLsymboltab * tab,
	int i,
	const char *pref,
	char uname[ILL_namebufsize]);

extern const char *ILLsymboltab_get (
	const ILLsymboltab * tab,
	int i);
extern int ILLsymboltab_rename (
	ILLsymboltab * h,
	int i,
	const char *new_name);

#endif /* __SYMTAB_H */
