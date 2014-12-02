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

/* RCS_INFO = "$RCSfile: ftest.c,v $ $Revision: 1.2 $ $Date: 2003/11/05 16:49:52 $"; */

#include "qs_config.h"
#include "config.h"
#include "iqsutil.h"
#include "lpdefs.h"
#include "factor.h"

#define NREP 2
#define MAXITER 10000
#undef  DENSE_READ
//static int TRACE = 0;

#if 0
static int handle_singularity (
	void *sdata,
	int c,
	int r)
{
	fprintf (stderr, "singular basis, replace column %d by row %d logical\n",
					 r, c);
	return 0;
}
#endif

int main (
	int ac,
	char **av)
{
	factor_work f;
	int dim = 0;
	int ncol = 0;
	int nzcnt = 0;
	int *basis = (int *) NULL;
	int *cbeg = (int *) NULL;
	int *clen = (int *) NULL;
	int *cind = (int *) NULL;
	EGlpNum_t *coef = (EGlpNum_t *) NULL;
	int niter = 0;
	char str[4096];
	char*l_argv[1024];
	int l_argc,l_par;
	char cmd[MAXITER];
	int col[MAXITER];
	svector a[MAXITER];
	svector x;
	svector upd;
	int rval;
	int i;
	int j;
	int rep;
	int nz;
	EGioFile_t *fin = 0;
	char c;
	int iter;
	int nsing, *psrow, *pscol;
	double szeit;
	double factor_szeit;
	double ftran_szeit;
	double btran_szeit;
	double update_szeit;
	double factor_zeit = 0.0;
	double ftran_zeit = 0.0;
	double btran_zeit = 0.0;
	double update_zeit = 0.0;
	double dtmp;
	QSexactStart();
	ILLfactor_init_factor_work (&f);
	for (i = 0; i < MAXITER; i++)
	{
		ILLsvector_init (&a[i]);
	}
	ILLsvector_init (&x);
	ILLsvector_init (&upd);

	szeit = ILLutil_zeit ();

	if (ac > 1)
	{
		fin = EGioOpen (av[1], "r");
		if (!fin)
		{
			perror (av[1]);
			fprintf (stderr, "Unable to open %s for input\n", av[1]);
			return 1;
		}
	}
	else
	{
		fin = EGioOpenFILE(stdin);
	}
	do{
		EGioGets(str,4096,fin);
		EGioNParse(str,1024," ","%#",&i,l_argv);
	}while(i==0 && !EGioEof(fin));
	WARNING (i != 3, "Couldn't read dimension, columns and nonzeros, only read"
					 " %d elements", i);
	dim = atoi(l_argv[0]);
	ncol = atoi(l_argv[1]);
	nzcnt = atoi(l_argv[2]);
	fprintf (stderr, "Using dimension %d columns %d and nonzero %d\n", dim,
					 ncol, nzcnt);

	ILL_SAFE_MALLOC (basis, dim, int);
	ILL_SAFE_MALLOC (cbeg, ncol, int);
	ILL_SAFE_MALLOC (clen, ncol, int);
	ILL_SAFE_MALLOC (cind, nzcnt, int);

	coef = EGlpNumAllocArray (nzcnt);

	do{
		EGioGets(str,4096,fin);
		EGioNParse(str,1024," ","%#",&i,l_argv);
	}while(i==0 && !EGioEof(fin));
	TESTG((rval=(i!=dim)), CLEANUP, "Wrong format!");
	for (i = 0; i < dim; i++)
	{
		basis[i] = atoi(l_argv[i]);
	}

	nz = 0;
	for (i = 0; i < ncol; i++)
	{
		cbeg[i] = nz;
		do
		{
			EGioGets(str,4096,fin);
			EGioNParse(str,1024," ","%#",&l_argc,l_argv);
		}while(l_argc==0 && !EGioEof(fin));
		TESTG((rval=(l_argc<1)), CLEANUP, "Wrong format!");
		clen[i] = l_argv[0];
		TESTG((rval=(l_argc!= 2*clen[i]+1)), CLEANUP, "Wrong format or too long lines!");
		for (j = 0; j < clen[i]; j++)
		{
			cind[nz] = atoi(l_argv[2*i+1]);
			dtmp = strtod(l_argv[2*i+2],0);
			EGlpNumSet (coef[nz], dtmp);
			nz++;
		}
	}

	rval = ILLsvector_alloc (&x, dim);
	ILL_CLEANUP_IF (rval);
	rval = ILLsvector_alloc (&upd, dim);
	ILL_CLEANUP_IF (rval);

	while (niter < MAXITER)
	{
		do
		{
			EGioGets(str,4096,fin);
			EGioNParse(str,1024," ","%#",&l_argc,l_argv);
		}while(l_argc==0 && !EGioEof(fin));
		l_par = 0;
		if(l_argc==0) break;
		TESTG((rval=(l_argc<2)), CLEANUP,"Wrong Format");
		c = l_argv[l_par++][0];
		cmd[niter] = c;

		switch (c)
		{
		case 'u':
			col[niter] = atoi(l_argv[l_par++]);
		case 'F':
		case 'f':
		case 'b':
#ifdef DENSE_READ
			nzcnt = 0;
			TESTG((rval=(l_argc!=dim+l_par)), CLEANUP, "Wrng format or too long lines!");
			for (i = 0; i < dim; i++)
			{
				dtmp = strtod(l_argv[l_par++],0);
				if (dtmp != 0.0)
				{
					x.indx[nzcnt] = i;
					EGlpNumSet (x.coef[nzcnt], dtmp);
					nzcnt++;
				}
			}
#else
			nzcnt = atoi(l_argv[l_par++]);
			TESTG((rval=(l_argc!=2*nzcnt+l_par)), CLEANUP, "Wrong format or too log lines!");
			for (i = 0; i < nzcnt; i++)
			{
				x.indx[i] = atoi(l_argv[l_par++]);
				dtmp = strtod(l_argv[l_par++],0);
				EGlpNumSet (x.coef[i], dtmp);

			}
#endif
			x.nzcnt = nzcnt;
			rval = ILLsvector_copy (&x, &a[niter]);
			ILL_CLEANUP_IF (rval);
			break;
		}
		niter++;
	}

	printf ("Matrix and iterations read in %.2f seconds\n",
					ILLutil_zeit () - szeit);
	fflush (stdout);

	szeit = ILLutil_zeit ();

	for (rep = 0; rep < NREP; rep++)
	{

		factor_szeit = ILLutil_zeit ();

#if 0
		ILLfactor_init_factor_work (&f);
#endif

		rval = ILLfactor_create_factor_work (&f, dim);
		ILL_CLEANUP_IF (rval);

		rval =
			ILLfactor (&f, basis, cbeg, clen, cind, coef, &nsing, &psrow, &pscol);
		ILL_CLEANUP_IF (rval);

		factor_zeit += ILLutil_zeit () - factor_szeit;

		for (iter = 0; iter < niter; iter++)
		{
			switch (cmd[iter])
			{
			case 'f':
				ftran_szeit = ILLutil_zeit ();
				ILLfactor_ftran (&f, &a[iter], &x);
				ftran_zeit += ILLutil_zeit () - ftran_szeit;
				break;
			case 'F':
				ftran_szeit = ILLutil_zeit ();
				ILLfactor_ftran_update (&f, &a[iter], &upd, &x);
				ftran_zeit += ILLutil_zeit () - ftran_szeit;
				break;
			case 'b':
				btran_szeit = ILLutil_zeit ();
				ILLfactor_btran (&f, &a[iter], &x);
				btran_zeit += ILLutil_zeit () - btran_szeit;
				break;
			case 'u':
				update_szeit = ILLutil_zeit ();
				ILLfactor_update (&f, &a[iter], col[iter], &nz);
				update_zeit += ILLutil_zeit () - update_szeit;
				break;
			}
		}
		ILLfactor_free_factor_work (&f);
	}
	printf ("%d reps of %d steps finished in %.2f seconds\n", rep, niter,
					ILLutil_zeit () - szeit);
	printf ("factor: %.2f\n", factor_zeit);
	printf ("ftran:  %.2f\n", ftran_zeit);
	printf ("btran:  %.2f\n", btran_zeit);
	printf ("update: %.2f\n", update_zeit);
	fflush (stdout);

	rval = 0;

CLEANUP:
	if (fin)
	{
		EGioClose (fin);
	}
	ILLfactor_free_factor_work (&f);
	ILL_IFFREE (basis, int);
	ILL_IFFREE (cbeg, int);
	ILL_IFFREE (clen, int);
	ILL_IFFREE (cind, int);

	EGlpNumFreeArray (coef);
	for (i = 0; i < niter; i++)
	{
		ILLsvector_free (&a[i]);
	}
	ILLsvector_free (&upd);
	ILLsvector_free (&x);
	QSexactClear();
	return rval;
}
