/* EGlib "Efficient General Library" provides some basic structures and
 * algorithms commons in many optimization algorithms.
 *
 * Copyright (C) 2005 Daniel Espinoza and Marcos Goycoolea.
 * 
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation; either version 2.1 of the License, or (at your
 * option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public 
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA 
 * */
/** @file 
 * @ingroup EGeHeap */
/** @addtogroup EGeHeap */
/** @{ */
/* ========================================================================= */
/** @brief This program reads a list of double values from a text file and 
 * uses a heap to sort them. It also allows the user to change the value of 
 * one of the read numbers after they have been placed in the heap. The 
 * purpose of this program is to illustrate the use of the @ref EGeHeap 
 * structure and its associated functions. 
 * 
 * The input file format is: 
 * n
 * Value_0
 * Value_1
 * Value_2
 * Value_3
 *   ...  
 * Value_n 
 * */
/* ========================================================================= */
#include "qs_config.h"
#include "QSopt_ex.h"

/* ========================================================================= */
/** @brief display the usage message for this program */
void eheap_usage (char *program)
{
	fprintf (stdout, "Usage: %s [options]\n", program);
	fprintf (stdout, "Options:\n");
	fprintf (stdout, "     -d n   'd' value.\n");
	fprintf (stdout, "     -f n   file name.\n");
	fprintf (stdout, "     -c n   item whose value to change.\n");
	fprintf (stdout, "     -v n   new item value (0 by default).\n");

}

/* ========================================================================= */
/** @brief parse the input argumenbts for the program */
int eheap_parseargs (int argc,
										 char **argv,
										 unsigned int *d,
										 unsigned int *ch,
										 EGlpNum_t * v,
										 char **file_name)
{
	int c;
	while ((c = getopt (argc, argv, "f:d:c:v:")) != EOF)
	{
		switch (c)
		{
		case 'f':
			*file_name = optarg;
			break;
		case 'd':
			*d = atoi (optarg);
			break;
		case 'c':
			*ch = atoi (optarg);
			break;
		case 'v':
			EGlpNumReadStr (*v, optarg);
			break;
		default:
			eheap_usage (argv[0]);
			return 1;
		}
	}
	/* reporting the options */
	if (!*file_name)
	{
		eheap_usage (argv[0]);
		return 1;
	}
	fprintf (stdout, "Parsed Options:\n");
	fprintf (stdout, "input         : %s\n", *file_name);
	fprintf (stdout, "d             : %u\n", *d);
	if (*ch != UINT_MAX)
	{
		fprintf (stdout, "c             : %u\n", *ch);
		fprintf (stdout, "v             : %lf\n", EGlpNumToLf (*v));
	}
	return 0;
}

/* ========================================================================= */
/** @brief main function */
int main (int argc,
					char **argv)
{

	int rval = 0;
	unsigned int i,
	  c = UINT_MAX,
	  d = 2;
	char *file_name = 0,
	  str1[1024];
	FILE *file;
	unsigned int nval;
	EGlpNum_t v;
	EGeHeap_t my_heap;
	EGeHeapCn_t *all_cn = 0,
	 *ccn;
	EGlpNumStart();
	/* set signal and limits */
	EGsigSet(rval,CLEANUP);
	EGsetLimits(3600.0,4294967295UL);
	EGlpNumInitVar (v);
	EGlpNumZero (v);
	EGeHeapInit (&my_heap);
	rval = EGeHeapCheck (&my_heap);
	CHECKRVALG (rval, CLEANUP);

	/* Parse command line input to get 'file name' and 'd'. */
	rval = eheap_parseargs (argc, argv, &d, &c, &v, &file_name);
	CHECKRVALG (rval, CLEANUP);
	EGeHeapChangeD (&my_heap, d);
	rval = EGeHeapCheck (&my_heap);
	CHECKRVALG (rval, CLEANUP);

	/* Read the objects to sort from the file */
	file = fopen (file_name, "r");
	TEST (!file, "unable to open file %s", file_name);
	fscanf (file, "%u", &nval);
	all_cn = EGsMalloc (EGeHeapCn_t, nval);
	EGeHeapResize (&my_heap, nval);
	rval = EGeHeapCheck (&my_heap);
	CHECKRVALG (rval, CLEANUP);
	IFMESSAGE (1, "Inserting %u elements into the heap", nval);
	for (i = 0; i < nval; i++)
	{
		EGeHeapCnInit (all_cn + i);
		fscanf (file, "%s", str1);
		EGlpNumReadStr (all_cn[i].val, str1);
		IFMESSAGE (1, "Adding value (%s,%lg) to the heap", str1,
						 EGlpNumToLf (all_cn[i].val));
		rval = EGeHeapAdd (&my_heap, all_cn + i);
		CHECKRVALG (rval, CLEANUP);
		rval = EGeHeapCheck (&my_heap);
		CHECKRVALG (rval, CLEANUP);
	}
	fclose (file);

	/* Check if change value is in range */
	TESTG (c != UINT_MAX && c >= nval, CLEANUP,
				 "Change item (%u) is out of range (only %u objects)", c, nval);

	/* Popping the values from the heap */
	fprintf (stderr, "\nRemoving:\n\n");
	for (i = 0; i < nval; i++)
	{
		ccn = EGeHeapGetMin (&my_heap);
		IFMESSAGE (1, "%u: item %zd : %lg", i, ccn - all_cn,
						 EGlpNumToLf (ccn->val));
		EGeHeapDel (&my_heap, ccn);
		rval = EGeHeapCheck (&my_heap);
		CHECKRVALG (rval, CLEANUP);
	}

	if (c == UINT_MAX)
		goto CLEANUP;

	/* Re-inserting the values into the heap */
	fprintf (stderr, "\nRe-Inserting.\n\n");
	for (i = 0; i < nval; i++)
	{
		rval = EGeHeapAdd (&my_heap, all_cn + i);
		CHECKRVALG (rval, CLEANUP);
		rval = EGeHeapCheck (&my_heap);
		CHECKRVALG (rval, CLEANUP);
	}

	/* Changing value of an item */
	fprintf (stderr, "Changing value of item %u from %lf to %lf.\n", c,
					 EGlpNumToLf (all_cn[c].val), EGlpNumToLf (v));
	rval = EGeHeapChangeVal (&my_heap, all_cn + c, v);
	CHECKRVALG (rval, CLEANUP);
	rval = EGeHeapCheck (&my_heap);
	CHECKRVALG (rval, CLEANUP);

	/* Popping the values from the heap */
	fprintf (stderr, "\nRemoving:\n\n");
	for (i = 0; i < nval; i++)
	{
		//ccn = EGeHeapGetMin (&my_heap);
		ccn = my_heap.cn[my_heap.sz-1];
		IFMESSAGE (1, "%u: item %zd : %lf", i, ccn - all_cn,
						 EGlpNumToLf (ccn->val));
		EGeHeapDel (&my_heap, ccn);
		rval = EGeHeapCheck (&my_heap);
		CHECKRVALG (rval, CLEANUP);
	}

	/* Liberating allocated memory */
CLEANUP:
	if (all_cn) EGfree (all_cn);
	EGeHeapClear (&my_heap);
	EGlpNumClearVar (v);
	EGlpNumClear();
	return rval;
}
/* ========================================================================= */
/** @} */
