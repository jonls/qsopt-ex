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

/* RCS_INFO = "$RCSfile: dstruct.c,v $ $Revision: 1.2 $ $Date: 2003/11/05 16:49:52 $"; */
static int TRACE = 0;

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "logging-private.h"

#include "allocrus.h"
#include "eg_lpnum.h"
#include "eg_io.h"
#include "except.h"

#include "dstruct_EGLPNUM_TYPENAME.h"
#include "qsopt_EGLPNUM_TYPENAME.h"
#include "lpdefs_EGLPNUM_TYPENAME.h"

/****************************************************************************/
/*                                                                          */
/*                            EGLPNUM_TYPENAME_svector                                       */
/*                                                                          */
/*  Written by:  Applegate, Cook, Dash                                      */
/*  Date:                                                                   */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/****************************************************************************/

void EGLPNUM_TYPENAME_ILLsvector_init (
	EGLPNUM_TYPENAME_svector * s)
{
	s->nzcnt = 0;
	s->indx = 0;
	s->coef = 0;
}

void EGLPNUM_TYPENAME_ILLsvector_free (
	EGLPNUM_TYPENAME_svector * s)
{
	ILL_IFFREE (s->indx, int);

	EGLPNUM_TYPENAME_EGlpNumFreeArray (s->coef);
	s->nzcnt = 0;
}

int EGLPNUM_TYPENAME_ILLsvector_alloc (
	EGLPNUM_TYPENAME_svector * s,
	int nzcnt)
{
	int rval = 0;

	s->nzcnt = nzcnt;
	if (nzcnt == 0)
	{
		s->indx = 0;
		s->coef = 0;
	}
	else
	{
		ILL_SAFE_MALLOC (s->indx, nzcnt, int);

		s->coef = EGLPNUM_TYPENAME_EGlpNumAllocArray (nzcnt);
	}
	return 0;
CLEANUP:
	ILL_IFFREE (s->indx, int);
	EGLPNUM_TYPENAME_EGlpNumFreeArray (s->coef);
	ILL_RETURN (rval, "EGLPNUM_TYPENAME_ILLsvector_alloc");
}

int EGLPNUM_TYPENAME_ILLsvector_copy (
	const EGLPNUM_TYPENAME_svector * s_in,
	EGLPNUM_TYPENAME_svector * s_out)
{
	int i;
	int nzcnt = s_in->nzcnt;
	int rval = 0;

	rval = EGLPNUM_TYPENAME_ILLsvector_alloc (s_out, nzcnt);
	ILL_CLEANUP_IF (rval);
	for (i = 0; i < nzcnt; i++)
	{
		s_out->indx[i] = s_in->indx[i];
		EGLPNUM_TYPENAME_EGlpNumCopy (s_out->coef[i], s_in->coef[i]);
	}

CLEANUP:
	ILL_RETURN (rval, "EGLPNUM_TYPENAME_ILLsvector_copy");
}

/****************************************************************************/
/*                                                                          */
/*                            EGLPNUM_TYPENAME_heap                                          */
/*                                                                          */
/*  Written by:  Applegate, Cook, Dash                                      */
/*  Date:                                                                   */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/****************************************************************************/

#define DEBUG_HEAP 0

#define HEAP_D 3
#define HEAP_UP(x) (((x)-1)/HEAP_D)
#define HEAP_DOWN(x) (((x)*HEAP_D)+1)

static int siftup (
	EGLPNUM_TYPENAME_heap * h,
	int hloc,
	int ix),
  siftdown (
	EGLPNUM_TYPENAME_heap * h,
	int hloc,
	int ix),
  maxchild (
	EGLPNUM_TYPENAME_heap * h,
	int hloc);

static int siftup (
	EGLPNUM_TYPENAME_heap * h,
	int hloc,
	int ix)
{
	int i = hloc;
	int p = HEAP_UP (i);
	EGLPNUM_TYPE val;

	EGLPNUM_TYPENAME_EGlpNumInitVar (val);
	EGLPNUM_TYPENAME_EGlpNumCopy (val, h->key[ix]);

	while (i > 0 && EGLPNUM_TYPENAME_EGlpNumIsLess (h->key[h->entry[p]], val))
	{
		h->entry[i] = h->entry[p];
		h->loc[h->entry[i]] = i;
		i = p;
		p = HEAP_UP (p);
	}
	h->entry[i] = ix;
	h->loc[ix] = i;
	ILL_IFTRACE2 ("%s:%la:%d:%d:%d\n", __func__, EGLPNUM_TYPENAME_EGlpNumToLf (val), hloc, ix, i);
	EGLPNUM_TYPENAME_EGlpNumClearVar (val);
	return i;
}

static int siftdown (
	EGLPNUM_TYPENAME_heap * h,
	int hloc,
	int ix)
{
	int i = hloc;
	int c = maxchild (h, i);
	EGLPNUM_TYPE val;

	EGLPNUM_TYPENAME_EGlpNumInitVar (val);
	EGLPNUM_TYPENAME_EGlpNumCopy (val, h->key[ix]);
	ILL_IFTRACE2 ("%s:%d:%d:%d:%la", __func__, hloc, ix, c, EGLPNUM_TYPENAME_EGlpNumToLf (val));

	while (c != -1 && EGLPNUM_TYPENAME_EGlpNumIsLess (val, h->key[h->entry[c]]))
	{
		h->entry[i] = h->entry[c];
		h->loc[h->entry[i]] = i;
		i = c;
		c = maxchild (h, c);
	}
	h->entry[i] = ix;
	h->loc[ix] = i;
	EGLPNUM_TYPENAME_EGlpNumClearVar (val);
	ILL_IFTRACE2 ("%s:%d:%d\n", __func__, ix, i);
	return i;
}

//extern EGLPNUM_TYPE EGLPNUM_TYPENAME_ILL_MINDOUBLE;
static int maxchild (
	EGLPNUM_TYPENAME_heap * h,
	int hloc)
{
	int i;
	int mc = -1;
	int hmin = HEAP_D * hloc + 1;
	int hmax = HEAP_D * hloc + HEAP_D;
	EGLPNUM_TYPE val;

	EGLPNUM_TYPENAME_EGlpNumInitVar (val);
	EGLPNUM_TYPENAME_EGlpNumCopy (val, EGLPNUM_TYPENAME_ILL_MINDOUBLE);
	ILL_IFTRACE2 (" %s:%d", __func__, hloc);

	for (i = hmin; i <= hmax && i < h->size; i++)
		if (EGLPNUM_TYPENAME_EGlpNumIsLess (val, h->key[h->entry[i]]))
		{
			EGLPNUM_TYPENAME_EGlpNumCopy (val, h->key[h->entry[i]]);
			mc = i;
			ILL_IFTRACE2 (":%d:%la", mc, EGLPNUM_TYPENAME_EGlpNumToLf (val));
		}
	EGLPNUM_TYPENAME_EGlpNumClearVar (val);
	ILL_IFTRACE2 ("\n");
	return mc;
}

#if DEBUG_HEAP > 0

static void printheap (
	EGLPNUM_TYPENAME_heap * h)
{
	int i;

	QSlog("entry (%d): ", h->size);
	for (i = 0; i < h->size; i++)
		QSlog("%d ", h->entry[i]);
	QSlog(" loc: ");
	for (i = 0; i < h->maxsize; i++)
		QSlog("%d ", h->loc[i]);
	QSlog("\n key: ");
	for (i = 0; i < h->maxsize; i++)
		QSlog("%la ", EGLPNUM_TYPENAME_EGlpNumToLf (h->key[i]));
	QSlog("\n key(sorted): ");
	for (i = 0; i < h->size; i++)
		QSlog("%la ", EGLPNUM_TYPENAME_EGlpNumToLf (h->key[h->entry[i]]));
}

static void heapcheck (
	EGLPNUM_TYPENAME_heap * h)
{
	int i, tcnt = 0;

	for (i = 0; i < h->maxsize; i++)
	{
		if (h->loc[i] < -1)
			QSlog("error in EGLPNUM_TYPENAME_heap\n");
		else if (h->loc[i] > -1)
			tcnt++;
	}
	if (tcnt != h->size)
		QSlog("error 3 in EGLPNUM_TYPENAME_heap\n");

	for (i = 0; i < h->size; i++)
	{
		if (h->loc[h->entry[i]] != i)
			QSlog("error 1 in EGLPNUM_TYPENAME_heap\n");
		if (!EGLPNUM_TYPENAME_EGlpNumIsNeqqZero (h->key[h->entry[i]]))
			QSlog("error 2 in EGLPNUM_TYPENAME_heap\n");
		if (EGLPNUM_TYPENAME_EGlpNumIsLess (h->key[h->entry[HEAP_UP (i)]], h->key[h->entry[i]]))
			QSlog("error 4 in EGLPNUM_TYPENAME_heap\n");
	}
}

#endif

void EGLPNUM_TYPENAME_ILLheap_insert (
	EGLPNUM_TYPENAME_heap * const h,
	int const ix)
{
	int i = h->size;

	ILL_IFTRACE ("%s:%d:%la\n", __func__, ix, EGLPNUM_TYPENAME_EGlpNumToLf (h->key[ix]));

	i = siftup (h, i, ix);
	h->size++;

#if DEBUG_HEAP > 0
	heapcheck (h);
#endif
#if DEBUG_HEAP > 1
	printheap (h);
#endif
}

void EGLPNUM_TYPENAME_ILLheap_modify (
	EGLPNUM_TYPENAME_heap * const h,
	int const ix)
{
	int i = h->loc[ix];
	int pi = i;

	ILL_IFTRACE ("%s:%d\n", __func__, ix);

	if (h->loc[ix] == -1)
		return;
	i = siftup (h, i, ix);
	if (pi == i)
		i = siftdown (h, i, ix);

#if DEBUG_HEAP > 0
	heapcheck (h);
#endif
#if DEBUG_HEAP > 1
	printheap (h);
#endif
}

void EGLPNUM_TYPENAME_ILLheap_delete (
	EGLPNUM_TYPENAME_heap * const h,
	int const ix)
{
	int i = h->loc[ix];
	int pi = i;
	int nix = h->entry[h->size - 1];

	ILL_IFTRACE ("%s:%d:%d:%d\n", __func__, ix, nix, pi);

	h->loc[ix] = -1;
	h->size--;
	if (nix == ix)
	{
#if DEBUG_HEAP > 0
		heapcheck (h);
#endif
#if DEBUG_HEAP > 1
		printheap (h);
#endif
		return;
	}

	h->entry[i] = nix;
	h->loc[nix] = i;

	i = siftup (h, i, nix);
	ILL_IFTRACE ("%s:%d:%d:%d:%d\n", __func__, ix, nix, pi, i);
	if (pi == i)
		siftdown (h, i, nix);

#if DEBUG_HEAP > 0
	heapcheck (h);
#endif
#if DEBUG_HEAP > 1
	printheap (h);
#endif
}

int EGLPNUM_TYPENAME_ILLheap_findmin (
	EGLPNUM_TYPENAME_heap * const h)
{
	if (h->hexist == 0 || h->size <= 0)
		return -1;
	return h->entry[0];
}

void EGLPNUM_TYPENAME_ILLheap_init (
	EGLPNUM_TYPENAME_heap * const h)
{
	h->entry = NULL;
	h->loc = NULL;
	h->key = NULL;
	h->hexist = 0;
}

int EGLPNUM_TYPENAME_ILLheap_build (
	EGLPNUM_TYPENAME_heap * const h,
	int const nelems,
	EGLPNUM_TYPE * key)
{
	int rval = 0;
	int i, n = 0;

	ILL_IFTRACE ("%s:%d\n", __func__, nelems);

	h->hexist = 1;
	h->size = 0;
	h->maxsize = nelems;
	h->key = key;
	ILL_SAFE_MALLOC (h->entry, nelems, int);
	ILL_SAFE_MALLOC (h->loc, nelems, int);

	for (i = 0; i < nelems; i++)
	{
		if (EGLPNUM_TYPENAME_EGlpNumIsGreatZero (key[i]))
		{
			h->entry[n] = i;
			h->loc[i] = n;
			n++;
		}
		else
			h->loc[i] = -1;
	}
	h->size = n;
	for (i = n - 1; i >= 0; i--)
	{
		ILL_IFTRACE2 ("insert %la\n", EGLPNUM_TYPENAME_EGlpNumToLf (h->key[h->entry[i]]));
		siftdown (h, i, h->entry[i]);
	}

#if DEBUG_HEAP > 0
	heapcheck (h);
#endif
#if DEBUG_HEAP > 1
	printheap (h);
#endif

CLEANUP:
	if (rval)
		EGLPNUM_TYPENAME_ILLheap_free (h);
	ILL_RETURN (rval, "EGLPNUM_TYPENAME_ILLheap_init");
}

void EGLPNUM_TYPENAME_ILLheap_free (
	EGLPNUM_TYPENAME_heap * const h)
{
	if (h->hexist)
	{
		ILL_IFFREE (h->entry, int);
		ILL_IFFREE (h->loc, int);

		h->hexist = 0;
		h->maxsize = 0;
		h->size = 0;
	}
}


/****************************************************************************/
/*                                                                          */
/*                          matrix                                          */
/*                                                                          */
/*  Written by:  Applegate, Cook, Dash                                      */
/*  Date:                                                                   */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/****************************************************************************/

void EGLPNUM_TYPENAME_ILLmatrix_init (
	EGLPNUM_TYPENAME_ILLmatrix * A)
{
	if (A)
	{
		A->matval = 0;
		A->matcnt = 0;
		A->matbeg = 0;
		A->matind = 0;
		A->matcols = 0;
		A->matcolsize = 0;
		A->matrows = 0;
		A->matsize = 0;
		A->matfree = 0;
	}
}

void EGLPNUM_TYPENAME_ILLmatrix_free (
	EGLPNUM_TYPENAME_ILLmatrix * A)
{
	if (A)
	{
		EGLPNUM_TYPENAME_EGlpNumFreeArray (A->matval);
		ILL_IFFREE (A->matcnt, int);
		ILL_IFFREE (A->matbeg, int);
		ILL_IFFREE (A->matind, int);

		EGLPNUM_TYPENAME_ILLmatrix_init (A);
	}
}

void EGLPNUM_TYPENAME_ILLmatrix_prt (
	EGioFile_t * fd,
	EGLPNUM_TYPENAME_ILLmatrix * A)
{
	int j, k;

	if (A == NULL)
	{
		EGioPrintf (fd, "Matrix %p: empty\n", (void *) A);
	}
	else
	{
		EGioPrintf (fd, "Matrix %p: nrows = %d ncols = %d\n",
						 (void *) A, A->matrows, A->matcols);
		for (j = 0; j < A->matcols; j++)
		{
			EGioPrintf (fd, "col %d: ", j);
			for (k = A->matbeg[j]; k < A->matbeg[j] + A->matcnt[j]; k++)
			{
				EGioPrintf (fd, "row %d=%.3f ", A->matind[k], EGLPNUM_TYPENAME_EGlpNumToLf (A->matval[k]));
			}
			EGioPrintf (fd, "\n");
		}
	}
}
