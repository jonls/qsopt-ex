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

/* RCS_INFO = "$RCSfile: reader.c,v $ $Revision: 1.2 $ $Date: 2003/11/05 16:49:52 $"; */

#include "reader.h"
#include "qstruct.h"
#include "qsopt.h"

static int TRACE = 0;
static int TEST_ERROR_COLLECTOR = 0;
static int TEST_ERROR_MEMORY = 0;

static char *fname = (char *) NULL;
static char *out_lp = (char *) NULL;
static char *out_mps = (char *) NULL;
static int lpfile = 0;
static int stats = 0;

static void usage (
	char *s);

static int parseargs (
	int ac,
	char **av);

static int add_error (
	void *dest,
	QSformat_error error)
{
	const char *type = "Error";
	const char *line;
	int tp, i, at;
	FILE *out = (FILE *) dest;

	at = QSerror_get_pos (error);
	tp = QSerror_get_type (error);
	type = QSformat_error_type_string (tp);

	fprintf (out, "ADD: ");
	fprintf (out, "%s  line %d pos %d\n",
					 type, QSerror_get_line_number (error), at);
	line = QSerror_get_line (error);
	if (line != NULL)
	{
		fprintf (out, "LINE %s", line);
		if (at >= 0)
		{
			fprintf (out, ".....");
			for (i = 0; i <= (at - 1); i++)
			{
				if (line[i] == '\t')
				{
					fputc ('\t', out);
				}
				else
				{
					fputc ('.', out);
				}
			}
			fprintf (out, "^\n");
		}
	}
	else
	{
		fprintf (out, "NO LINE\n");
	}

	fprintf (out, "MSG: %s\n", QSerror_get_desc (error));
	return 0;
}

QSLIB_INTERFACE int reader_main (
	int ac,
	char **av)
{
	int rval = 0;
	int rvalmps = 0;
	int rvallp = 0;
	QSdata *p = NULL;
	ILLutil_timer timer_read;
	ILLutil_timer timer_write;
	FILE *fin = NULL;
	QSline_reader reader = NULL;
	QSerror_collector collector = NULL;
	QSerror_memory error_mem = NULL;
	QSformat_error e = NULL;

	if (parseargs (ac, av))
		goto CLEANUP;

	ILLutil_init_timer (&timer_read, "READER_READ");
	ILLutil_start_timer (&timer_read);
	if (TEST_ERROR_COLLECTOR || TEST_ERROR_MEMORY)
	{
		fin = fopen (fname, "r");
		reader = QSline_reader_new ((void *) fgets, fin);
		if (TEST_ERROR_COLLECTOR)
		{
			collector = QSerror_collector_new ((void *) add_error, stderr);
		}
		if (TEST_ERROR_MEMORY)
		{
			error_mem = QSerror_memory_create (0);
			ILL_CHECKnull (error_mem, "Could not make error memory");
			collector = QSerror_memory_collector_new (error_mem);
		}
		if (fin == NULL)
		{
			fprintf (stderr, "Can't open \"%s\" for reading.\n", fname);
		}
		rval = (fin == NULL) || (reader == NULL) || (collector == NULL);
		ILL_CLEANUP_IF (rval);

		QSline_reader_set_error_collector (reader, collector);
		p = QSget_prob (reader, fname, (lpfile == 0) ? "MPS" : "LP");

		if (TEST_ERROR_MEMORY)
		{
			int n = QSerror_memory_get_nerrors (error_mem);

			fprintf (stderr, "#error %d\n", n);
			for (e = QSerror_memory_get_last_error (error_mem);
					 e != NULL; e = QSerror_memory_get_prev_error (e))
			{
				QSerror_print (stderr, e);
			}
		}
	}
	else
	{
		if (lpfile == 0)
		{
			p = QSread_prob (fname, "MPS");
		}
		else
		{
			p = QSread_prob (fname, "LP");
		}
	}
	ILL_IFTRACE ("QSread_prob %s\n", fname);
	ILLutil_stop_timer (&timer_read, 1);
	rval = (p == NULL);
	fprintf (stdout,
					 "read \"%s\" %s.\n", fname, (rval == 0) ? "succeeded" : "failed");
	ILL_CLEANUP_IF (rval);

	if (stats)
	{
		fprintf (stdout, "\n");
		fprintf (stdout, "The problem \"%s\" has %d rows and %d cols.\n",
						 QSget_probname (p), QSget_rowcount (p), QSget_colcount (p));
	}
	ILLutil_init_timer (&timer_write, "READER_WRITE");
	ILLutil_start_timer (&timer_write);
	if (out_mps)
	{
		fprintf (stdout, "\n");
		rvalmps = QSwrite_prob (p, out_mps, "MPS");
		fprintf (stdout, "write \"%s\" %s.\n", out_mps,
						 (rvalmps == 0) ? "succeeded" : "failed");
	}
	if (out_lp)
	{
		fprintf (stdout, "\n");
		rvallp = QSwrite_prob (p, out_lp, "LP");
		fprintf (stdout, "write \"%s\" %s.\n", out_lp,
						 (rvallp == 0) ? "succeeded" : "failed");
	}
	ILLutil_stop_timer (&timer_write, 1);

	rval = rvalmps || rvallp;
	ILL_CLEANUP_IF (rval);

CLEANUP:

	if (p != NULL)
		QSfree_prob (p);
	if (fin != NULL)
		fclose (fin);
	if (reader != NULL)
		QSline_reader_free (reader);
	if (collector != NULL)
		QSerror_collector_free (collector);
	if (error_mem != NULL)
		QSerror_memory_free (error_mem);

	return rval;									/* main return */
}

static void usage (
	char *s)
{
	fprintf (stderr, "Usage: %s [- below -] prob_file\n", s);
	fprintf (stderr, "   -L    input file is in lp format (default: mps)\n");
	fprintf (stderr, "   -s    print information about problem to stdout\n");
	fprintf (stderr, "   -l f  write lp in LP format to file f\n");
	fprintf (stderr, "   -m f  write lp in MPS format to file f\n");
#if 0
	fprintf (stderr, "   -E    test error_memory\n");
	fprintf (stderr, "   -e    test error_collector\n");
#endif
}

static int parseargs (
	int ac,
	char **av)
{
	int c;
	int boptind = 1;
	char *boptarg = (char *) NULL;

	while ((c =
					ILLutil_bix_getopt (ac, av, "Ll:m:tseE", &boptind, &boptarg)) != EOF)
		switch (c)
		{
		case 'L':
			lpfile = 1;
			break;
		case 'l':
			out_lp = boptarg;
			break;
		case 'm':
			out_mps = boptarg;
			break;
		case 's':
			stats = 1;
			break;
		case 't':
			TRACE++;
			break;
#if 0
		case 'E':
			TEST_ERROR_MEMORY++;
			break;
		case 'e':
			TEST_ERROR_COLLECTOR++;
			break;
#endif
		case '?':
		default:
			usage (av[0]);
			return 1;
		}

	if (boptind != (ac - 1))
	{
		usage (av[0]);
		return 1;
	}

	fname = av[boptind++];
	return 0;
}


#ifndef WIN32
int main (
	int ac,
	char **av)
{
	return reader_main (ac, av);
}
#endif
