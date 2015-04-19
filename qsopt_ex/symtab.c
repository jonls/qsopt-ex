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

/* RCSINFO $Id: symboltab.c,v 1.2 2003/11/05 16:47:22 meven Exp $ */

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "eg_mem.h"

#include "logging-private.h"

#include "util.h"
#include "trace.h"
#include "except.h"
#include "symtab.h"


static int TRACE = 0;

static unsigned int stringhash (
	const char *key,
	int tsize);
static int look_it_up (
	ILLsymboltab * h,
	const char *s);
static int grow_symboltab (
	ILLsymboltab * h);
static int grow_namelist (
	ILLsymboltab * h);
static int add_string (
	ILLsymboltab * h,
	const char *s,
	int *symbol);
static void delete_from_list (
	ILLsymboltab * h,
	int del_ind,
	int prev_ind,
	int x);

#ifdef TRY_CODE
/* debug support */
static const char *get_str (
	const ILLsymboltab * h,
	int indx);
static void prt_xchain (
	const ILLsymboltab * h,
	int x);
static void prt_chain (
	const ILLsymboltab * h,
	char *s);
#endif


void ILLsymboltab_init (
	ILLsymboltab * h)
{
	h->tablesize = 0;
	h->strsize = 0;
	h->freedchars = 0;
	h->hashspace = 0;
	h->name_space = 0;
	h->strspace = 0;
	h->hashtable = (int *) NULL;
	h->nametable = (ILLsymbolent *) NULL;
	h->namelist = (char *) NULL;
}

void ILLsymboltab_free (
	ILLsymboltab * h)
{
	ILL_IFFREE (h->hashtable, int);

	ILL_IFFREE (h->nametable, ILLsymbolent);
	ILL_IFFREE (h->namelist, char);

	ILLsymboltab_init (h);
}

int ILLsymboltab_create (
	ILLsymboltab * h,
	int init_size)
{
	int rval = 0;
	int i;

	if (init_size <= 0)
		init_size = 1000;

	ILLsymboltab_free (h);

	h->tablesize = 0;
	h->strsize = 0;
	h->freedchars = 0;
	h->name_space = init_size;
	h->hashspace = ILLutil_nextprime (((unsigned) h->name_space));
#ifdef TRY_CODE
	h->strspace = init_size;
#else
	h->strspace = init_size * 5;
#endif
	h->index_ok = 0;

	ILL_SAFE_MALLOC (h->hashtable, h->hashspace, int);

	ILL_SAFE_MALLOC (h->nametable, h->name_space, ILLsymbolent);
	ILL_SAFE_MALLOC (h->namelist, h->strspace, char);

	for (i = 0; i < h->hashspace; i++)
	{
		h->hashtable[i] = ILL_SYM_NOINDEX;
	}

CLEANUP:
	if (rval)
	{
		ILLsymboltab_free (h);
	}
	ILL_RETURN (rval, "ILLsymboltab_create");
}

int ILLsymboltab_copy (
	ILLsymboltab * src,
	ILLsymboltab * dst)
{
	int rval = 0;
	int i;

	ILLsymboltab_free (dst);

	*dst = *src;
	ILL_SAFE_MALLOC (dst->hashtable, dst->hashspace, int);

	ILL_SAFE_MALLOC (dst->nametable, dst->name_space, ILLsymbolent);
	ILL_SAFE_MALLOC (dst->namelist, dst->strspace, char);

	for (i = 0; i < src->hashspace; i++)
	{
		dst->hashtable[i] = src->hashtable[i];
	}
	for (i = 0; i < src->tablesize; i++)
	{
		dst->nametable[i] = src->nametable[i];
	}
	for (i = 0; i < src->strsize; i++)
	{
		dst->namelist[i] = src->namelist[i];
	}

CLEANUP:
	if (rval)
	{
		ILLsymboltab_free (dst);
	}
	ILL_RETURN (rval, "ILLsymboltab_copy");
}


const char *ILLsymboltab_get (
	const ILLsymboltab * h,
	int i)
{
	char *name = NULL;

	ILL_FAILfalse_no_rval ((i >= 0) && (i <= h->tablesize), "Index out of range");
	if (h->nametable[i].symbol != -1)
		name = h->namelist + h->nametable[i].symbol;
CLEANUP:
	return name;
}

int ILLsymboltab_index_ok (
	ILLsymboltab * h)
{
	return (h && (h->index_ok == 1));
}

int ILLsymboltab_index_reset (
	ILLsymboltab * h,
	int icount,
	char **names)
{
	int rval = 0;
	int i, k;

	/* Note may have tablesize 1 larger than icount, due to objname */

	if ((h->tablesize != icount) && (h->tablesize != icount + 1))
	{
		QSlog("symbol table (%d) does not match reset list (%d)",
								h->tablesize, icount);
		rval = 1;
		ILL_CLEANUP;
	}

	for (i = 0; i < icount; i++)
	{
		k = look_it_up (h, (const char *) names[i]);
		if (k)
		{
			QSlog("Symbol %s is not in table", names[i]);
			rval = 1;
			ILL_CLEANUP;
		}
		k = h->the_index;
		h->nametable[k].index = i;
	}

	h->index_ok = 1;

CLEANUP:
	return rval;
}

int ILLsymboltab_getindex (
	ILLsymboltab * h,
	const char *name,
	int *hindex)
{
	int rval = 0;
	int k;

	*hindex = -1;

	if (!h || !h->index_ok)
	{
		QSlog("symbol table index out of date");
		rval = 1;
		ILL_CLEANUP;
	}

	k = look_it_up (h, name);
	if (k)
	{
		QSlog("Symbol %s is not in table", name);
		ILL_CLEANUP;
	}
	k = h->the_index;

	*hindex = h->nametable[k].index;

CLEANUP:

	ILL_RETURN (rval, "ILLsymboltab_getindex");
}

int ILLsymboltab_rename (
	ILLsymboltab * h,
	int i,
	const char *new_name)
{
	int rval = 0;
	int symbol = 0;

	ILL_FAILfalse ((i >= 0) && (i <= h->tablesize), "Index out of range");
	if ((new_name != NULL) && (look_it_up (h, new_name) == 0))
	{
		rval = (i != h->the_index);
		ILL_RETURN (rval, "ILLsymboltab_rename");
	}
	if (h->nametable[i].symbol != -1)
	{
		(void) look_it_up (h, h->namelist + h->nametable[i].symbol);
		ILL_FAILfalse (i == h->the_index, "must find it at i");
		delete_from_list (h, i, h->the_prev_index, h->the_hash);
	}
	if (new_name != NULL)
	{
		rval = add_string (h, new_name, &symbol);
		ILL_CLEANUP_IF (rval);
		h->the_hash = stringhash (new_name, h->hashspace);
		h->nametable[i].symbol = symbol;
		h->nametable[i].next = h->hashtable[h->the_hash];
		h->hashtable[h->the_hash] = i;
	}
	else
	{
		h->nametable[i].symbol = -1;
		h->nametable[i].next = ILL_SYM_NOINDEX;
	}
CLEANUP:
	ILL_RETURN (rval, "ILLsymboltab_rename");
}

int ILLsymboltab_lookup (
	ILLsymboltab * h,
	const char *s,
	int *ind)
{
	int rval = look_it_up (h, s);

	*ind = h->the_index;
	return rval;
}

int ILLsymboltab_contains (
	ILLsymboltab * tab,
	const char *name)
{
	return !look_it_up (tab, name);
}

void ILLsymboltab_size (
	const ILLsymboltab * h,
	int *p_size)
{
	*p_size = h->tablesize;
}

int ILLsymboltab_register (
	ILLsymboltab * h,
	const char *s,
	int itemindex,
	int *the_prev_index,
	int *existed)
{
	ILLsymbolent *nametable = h->nametable;
	int e, symbol = 0;
	int rval = 0;

	if (itemindex < 0)
	{
		h->index_ok = 0;
	}

	h->the_prev_index = ILL_SYM_NOINDEX;
	h->the_index = ILL_SYM_NOINDEX;
	if (s == NULL)
	{
		e = h->tablesize;
		h->the_index = e;
		*existed = 0;
		while (h->tablesize >= h->name_space)
		{
			rval = grow_symboltab (h);
			ILL_CLEANUP_IF (rval);
		}
		h->tablesize++;
		h->nametable[e].symbol = -1;
		h->nametable[e].index = itemindex;
		h->nametable[e].next = ILL_SYM_NOINDEX;
		ILL_IFTRACE ("register: %s NULL entry#=%d\n", (*existed) ? "OLD" : "NEW",
								 e);
	}
	else
	{
		*existed = !look_it_up (h, s);
		if (*existed)
		{
			ILL_IFTRACE ("register: OLD %s entry#=%d hash=%d\n",
									 s, h->the_index, h->the_hash);
			return 0;
		}

		rval = add_string (h, s, &symbol);
		ILL_CLEANUP_IF (rval);

		while (h->tablesize >= h->name_space)
		{
			rval = grow_symboltab (h);
			ILL_CLEANUP_IF (rval);
			h->the_hash = stringhash (s, h->hashspace);	/*hash changes with bigger table */
		}

		nametable = h->nametable;

		e = h->tablesize;
		h->tablesize++;
		h->the_prev_index = e;

		nametable[e].symbol = symbol;
		nametable[e].index = itemindex;
		nametable[e].next = h->hashtable[h->the_hash];
		h->hashtable[h->the_hash] = e;

		ILL_IFTRACE ("register: %s NULL entry#=%d\n", (*existed) ? "OLD" : "NEW",
								 e);
	}

CLEANUP:
	*the_prev_index = h->the_prev_index;
	ILL_RETURN (rval, "ILLsymboltab_register");
}

int ILLsymboltab_delete (
	ILLsymboltab * h,
	const char *s)
{
	int del_ind, rval = 0;
	char *last;

	ILL_FAILtrue (s == NULL, "must give non NULL str");

	rval = look_it_up (h, s);
	del_ind = h->the_index;
	ILL_CLEANUP_IF (rval);				/* was not in table */
	ILL_FAILfalse ((del_ind != ILL_SYM_NOINDEX) &&
								 (h->nametable[del_ind].symbol != -1),
								 "we should have found this non NULL str");
	h->index_ok = 0;
	delete_from_list (h, del_ind, h->the_prev_index, h->the_hash);

	h->tablesize--;
	if (del_ind != h->tablesize)
	{
		if (h->nametable[h->tablesize].symbol != -1)
		{
			last = h->namelist + h->nametable[h->tablesize].symbol;
			rval = look_it_up (h, last);
			ILL_FAILfalse ((rval == 0) && (h->the_index == h->tablesize),
										 "Should find last entry");
			if (h->the_prev_index != ILL_SYM_NOINDEX)
			{
				h->nametable[h->the_prev_index].next = del_ind;
			}
			else
			{
				h->hashtable[h->the_hash] = del_ind;
			}
		}
		h->nametable[del_ind] = h->nametable[h->tablesize];
	}
CLEANUP:
	ILL_RETURN (rval, "ILLsymboltab_delete");
}

void ILLsymboltab_prt (
	FILE * fd,
	ILLsymboltab * h)
{
	char *str;
	int i;

	for (i = 0; i < h->tablesize; i++)
	{
		if (h->nametable[i].symbol == -1)
		{
			fprintf (fd, "%d: NULL nohash\n", i);
		}
		else
		{
			str = h->namelist + h->nametable[i].symbol;
			fprintf (fd, "%d: %s hash=%d\n", i, str, stringhash (str, h->hashspace));
		}
	}
}

static int look_it_up (
	ILLsymboltab * h,
	const char *s)
{
	ILLsymbolent *nametable = h->nametable;
	char *namelist = h->namelist;
	int e;

	if(!h->hashspace) goto CLEANUP;
	ILL_FAILfalse_no_rval (s, "Should never call with NULL string");
	h->the_prev_index = ILL_SYM_NOINDEX;
	h->the_hash = stringhash (s, h->hashspace);
	for (e = h->hashtable[h->the_hash]; e != ILL_SYM_NOINDEX;
			 e = nametable[e].next)
	{
		if (strcmp (namelist + nametable[e].symbol, s) == 0)
		{
			h->the_index = e;
			ILL_IFTRACE ("look_it_up: OLD %s entry#=%d hash=%d\n", s, e, h->the_hash);
			return 0;
		}
		h->the_prev_index = e;
	}
CLEANUP:
	h->the_index = ILL_SYM_NOINDEX;
	ILL_IFTRACE ("look_it_up: NEW %s \n", s);
	return 1;
}


static void delete_from_list (
	ILLsymboltab * h,
	int del_ind,
	int prev_ind,
	int x)
{
	if (prev_ind != ILL_SYM_NOINDEX)
	{
		ILL_FAILtrue_no_rval (h->nametable[prev_ind].symbol == -1,
													"A NULL str with same hash ?");
		h->nametable[prev_ind].next = h->nametable[del_ind].next;
	}
	else
	{
		h->hashtable[x] = h->nametable[del_ind].next;
	}
	h->freedchars += strlen (h->namelist + h->nametable[del_ind].symbol) + 1;
CLEANUP:
	;
}


static int grow_symboltab (
	ILLsymboltab * h)
{
	int newnamespace, newhashspace, *newhash, tablesize;
	ILLsymbolent *newname;
	char *namelist = h->namelist;
	int i;
	unsigned int x;
	int rval = 0;

	newnamespace = h->name_space * 2;
	newhashspace = ILLutil_nextprime (((unsigned) newnamespace));
	// rval = ILLutil_reallocrus_count ((void **) (&h->nametable), newnamespace,
	//                                 sizeof (ILLsymbolent));
	//ILL_CLEANUP_IF (rval);
	h->nametable = EGrealloc (h->nametable, sizeof (ILLsymbolent) * newnamespace);
	newname = h->nametable;

	ILL_SAFE_MALLOC (newhash, newhashspace, int);
	ILL_IFFREE (h->hashtable, int);

	h->hashtable = newhash;

	h->name_space = newnamespace;
	h->hashspace = newhashspace;

	newhash = h->hashtable;
	tablesize = h->tablesize;

	for (i = 0; i < newhashspace; i++)
	{
		newhash[i] = ILL_SYM_NOINDEX;
	}
	for (i = 0; i < tablesize; i++)
	{
		if (newname[i].symbol != -1)
		{
			x = stringhash (namelist + newname[i].symbol, newhashspace);
			newname[i].next = newhash[x];
			newhash[x] = i;
		}
	}
CLEANUP:
	ILL_RETURN (rval, "grow_symboltab");
}

static int grow_namelist (
	ILLsymboltab * h)
{
	int newstrspace, i, j, newsymbol, rval = 0;
	char *newnamelist, *newc;

	if (2 * h->freedchars >= h->strspace)
	{
		/* compact string array */
		ILL_SAFE_MALLOC (newnamelist, h->strspace, char);

		newc = newnamelist;
		for (i = 0; i < h->tablesize; i++)
		{
			if (h->nametable[i].symbol != -1)
			{
				newsymbol = newc - newnamelist;
				for (j = h->nametable[i].symbol; h->namelist[j] != '\0'; j++)
				{
					*newc = h->namelist[j];
					newc++;
				}
				*newc = '\0';
				newc++;
				h->nametable[i].symbol = newsymbol;
			}
		}
		ILL_IFFREE (h->namelist, char);

		h->namelist = newnamelist;
		h->strsize = newc - newnamelist;
		h->freedchars = 0;
	}
	else
	{
		newstrspace = h->strspace * 2;
		h->namelist = EGrealloc (h->namelist, sizeof (char) * newstrspace);
		//rval = ILLutil_reallocrus_count ((void **) &h->namelist, newstrspace,
		//                                 sizeof (char));
		//ILL_CLEANUP_IF (rval);
		h->strspace = newstrspace;
	}
CLEANUP:
	ILL_RETURN (rval, "grow_namelist");
}

static int add_string (
	ILLsymboltab * h,
	const char *s,
	int *symbol)
{
	int l, rval = 0;

	l = strlen (s) + 1;
	while (h->strsize + l > h->strspace)
	{
		rval = grow_namelist (h);
		ILL_CLEANUP_IF (rval);
	}
	strcpy (h->namelist + h->strsize, s);
	*symbol = h->strsize;
	h->strsize += l;
CLEANUP:
	ILL_RETURN (rval, "add_string");
}

static unsigned int stringhash (
	const char *key,
	int tsize)
{
	unsigned int x = 0;

#ifdef TRY_CODE
	while (*key)
	{
		x += *key;
		key++;
	}
#else
	while (*key)
	{
		x = 37 * x + *key;
		key++;
	}
#endif
	return x % tsize;
}

/**************************************************************************/
/* ILLsymboltab_unique_name and its support                                */
/**************************************************************************/
static void make_var (
	char *new_var,
	const char *prefix,
	char *name)
{
	size_t plen = strlen (prefix);
	size_t nlen = strlen (name);
	char *p;

	if (nlen + plen >= ILL_namebufsize)
	{
		nlen = ILL_namebufsize - plen - 1;
	}
	strcpy (new_var, prefix);
	p = new_var + plen;
	strncpy (p, name, nlen + 1);
}

int ILLsymboltab_uname (
	ILLsymboltab * symtab,
	char *name,
	const char *try_prefix1,
	const char *try_prefix2)
{
	int nvars = symtab->tablesize;
	int rval = 0;
	int i, found, numlen;
	const char *try_prefix[3];
	char prefix[ILL_namebufsize];
	char new_pre[ILL_namebufsize];
	char new[ILL_namebufsize];

	ILL_FAILtrue (try_prefix1 == NULL, "try_prefix must not be NULL");
	try_prefix[0] = try_prefix1;
	try_prefix[1] = try_prefix2;
	try_prefix[2] = NULL;
	new[0] = '\0';
	found = 0;
	for (i = 0; (!found) && try_prefix[i]; i++)
	{
		make_var (new, try_prefix[i], name);
		found = !ILLsymboltab_contains (symtab, new);
	}
	if (!found)
	{
		i = 0;
		sprintf (prefix, "%s", try_prefix[0]);
		numlen = (log10 ((double) (symtab->tablesize - 1) * 10)) + 1;
		while (!found)
		{
			ILL_FAILfalse (i <= nvars, "something wrong in find_unique_name");
			make_var (new_pre, prefix, name);
			new_pre[ILL_namebufsize - numlen - 1] = '\0';
			sprintf (new, "%s_%d", new_pre, i);
			found = !ILLsymboltab_contains (symtab, new);
			i++;
		}
	}

CLEANUP:
	strcpy (name, new);
	return rval;
}

void ILLsymboltab_unique_name (
	ILLsymboltab * tab,
	int i,
	const char *pref,
	char uname2[ILL_namebufsize])
{
	int notUnique;

	sprintf (uname2, "%d", i);
	notUnique = ILLsymboltab_uname (tab, uname2, pref, NULL);
	ILL_FAILtrue_no_rval (notUnique, "Programming error");
CLEANUP:
	return;
}

#ifdef TRY_CODE
/* debugging */
static const char *get_str (
	const ILLsymboltab * h,
	int indx)
{
	if (indx < 0 || indx >= h->tablesize)
	{
		return "INDEX OUT OF RANGE";
	}
	if (h->nametable[indx].symbol == -1)
	{
		return "NULL";
	}
	else
	{
		return h->namelist + h->nametable[indx].symbol;
	}
}

static void prt_xchain (
	const ILLsymboltab * h,
	int x)
{
	int e;
	char *str;

	x = x % h->hashspace;
	QSlog("chain hash %d:", x);
	for (e = h->hashtable[x]; e != ILL_SYM_NOINDEX; e = h->nametable[e].next)
	{
		if (h->nametable[e].symbol >= 0)
		{
			str = h->namelist + h->nametable[e].symbol;
			x = stringhash (str, h->hashspace);
			QSlog(" %s(h=%d, e=%d)", str, x, e);
		}
		else
		{
			QSlog(" NULL");
		}
	}
}

static void prt_chain (
	const ILLsymboltab * h,
	char *s)
{
	int e;

	if (s != NULL)
	{
		prt_xchain (h, stringhash (s, h->hashspace));
	}
	else
	{
		for (e = 0; e < h->hashspace; e++)
		{
			if (h->hashtable[e] != ILL_SYM_NOINDEX)
				prt_xchain (h, e);
		}
	}
}
#endif

#ifdef TRY_CODE
int main (
	int ac,
	char **av)
{

	int i, rval, index, pre_exist, nwords;
	const char *prefix[3];
	ILLsymboltab t, *tab = &t;
	char cmd[100], symbol[100], line[256], str[100];
	const char *s;
	int ok;

	TRACE = 1;
	prefix[0] = "C";
	prefix[1] = "c";
	prefix[2] = NULL;
	ILLsymboltab_init (tab);
	ILLsymboltab_create (tab, 1);

	fprintf (stdout, "> ");
	fflush (stdout);
	while (fgets (line, 100, stdin))
	{
		ok = 0;
		symbol[0] = '\0';
		nwords = sscanf (line, "%s%s%s", cmd, symbol, str);
		if (nwords >= 1)
		{
			fprintf (stdout, ":: %s", line);
			if ((nwords >= 2) && strcmp (cmd, "REG") == 0)
			{
				ok = 1;
				if (strcmp (symbol, "NULL") == 0)
				{
					rval = ILLsymboltab_register (tab, NULL, &index, &pre_exist);
				}
				else
				{
					rval = ILLsymboltab_register (tab, symbol, &index, &pre_exist);
				}
			}
			if ((nwords >= 2) && strcmp (cmd, "LOOK") == 0)
			{
				ok = 1;
				rval = ILLsymboltab_register (tab, symbol, &index, &pre_exist);
			}
			if ((nwords >= 2) && strcmp (cmd, "DEL") == 0)
			{
				ok = 1;
				rval = ILLsymboltab_delete (tab, symbol);
			}
			if ((nwords >= 2) && strcmp (cmd, "UNIQUE") == 0)
			{
				ok = 1;
				rval = ILLsymboltab_uname (tab, symbol, "c", "C");
			}
			if ((nwords >= 1) && strcmp (cmd, "PRT") == 0)
			{
				ok = 1;
				ILLsymboltab_prt (stdout, tab);
			}
			if ((nwords >= 2) && strcmp (cmd, "GET") == 0)
			{
				ok = 1;
				i = atoi (symbol);
				s = ILLsymboltab_get (tab, i);
				fprintf (stdout, "%d: %s\n", i, (s != NULL) ? s : "NULL");
			}
			if ((nwords >= 3) && strcmp (cmd, "RENAME") == 0)
			{
				ok = 1;
				i = atoi (symbol);
				if (strcmp (str, "NULL") == 0)
				{
					ILLsymboltab_rename (tab, i, NULL);
				}
				else
				{
					ILLsymboltab_rename (tab, i, str);
				}
			}
			if (strcmp (cmd, "CHAIN") == 0)
			{
				ok = 1;
				if ((nwords == 1) || strcmp (symbol, "NULL") == 0)
				{
					prt_chain (tab, NULL);
					fprintf (stdout,
									 "last %s(%d) strsize/space: %d/%d freedchars: %d\n",
									 get_str (tab, tab->tablesize - 1), tab->tablesize - 1,
									 tab->strsize, tab->strspace, tab->freedchars);
				}
				else
				{
					prt_chain (tab, symbol);
				}
			}
			if ((nwords >= 1) && strcmp (cmd, "INFO") == 0)
			{
				ok = 1;
				fprintf (stdout,
								 "last %s(%d) strsize/space: %d/%d freedchars: %d\n",
								 get_str (tab, tab->tablesize - 1), tab->tablesize - 1,
								 tab->strsize, tab->strspace, tab->freedchars);
			}
			if ((nwords >= 1) && strcmp (cmd, "STR") == 0)
			{
				ok = 1;
				for (i = 0; i < tab->strsize; i++)
				{
					if (tab->namelist[i] == '\0')
						fprintf (stdout, "\\0");
					else
						fprintf (stdout, "%c", tab->namelist[i]);
				}
				fprintf (stdout, "\n");
			}
		}
		if (ok == 0)
		{
			fprintf (stdout, "commands: REG, LOOK, DEL, RENAME, GET, %s\n",
							 "RENAME, UNIQUE, CHAIN, STR, INFO, PRT");
		}
		fprintf (stdout, "> ");
		fflush (stdout);
	}

	return 0;
}

/* sample input for testing symboltab */
#ifdef NEVER
REG aabb REG aabb								// its already there 
  REG abab REG abab							// its already there 
  REG baba REG baba							// its already there 
  REG bbaa REG bbaa							// its already there 
  CHAIN aabb STR								// all with same hash 
  DEL abab											// remove from middle
  CHAIN baba DEL aabb						// remove last
  CHAIN baba DEL bbaa						// remove first
  CHAIN baba DEL baba						// remove the only element 
  CHAIN PRT STR REG 123456789012	// compact name list no entries in table
  STR REG longnametoo INFO REG a REG b DEL 123456789012 DEL longnametoo INFO STR REG more	// compact name list with entries in table
  STR INFO PRT REG NULL					// add NULL str 
  PRT DEL b											// remove with NULL being last entry
  PRT CHAIN DEL a								// remove with NULL somewhere in table
  PRT DEL more CHAIN						// nothing left in table
  REG NULL											// add NULL str 
  REG NULL											// add NULL str 
  REG NULL											// add NULL str 
  REG more PRT CHAIN						// only chain is for more
  RENAME 0 more									// should fail
  RENAME 0 another							// own hash 
  RENAME 0 another							// must fail 
  CHAIN RENAME 1 orem						// existing hash 
  CHAIN RENAME 2 makestrarraygrowevenlongerwiththisid CHAIN PRT INFO STR DEL makestrarraygrowevenlongerwiththisid STR RENAME 3 somemore	// should trigger grownamelist
 
	PRT
	REG omemore
	STR
	CHAIN
	INFO
	UNIQUE a
	UNIQUE a
	REG Ca
	UNIQUE a REG ca UNIQUE a REG Ca_0 REG Ca_1 REG Ca_2 REG Ca_3 REG Ca_4 UNIQUE a
#endif
#endif
