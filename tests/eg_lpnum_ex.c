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
 * @ingroup EGlpNum
 * */
/** @addtogroup EGlpNum */
/** @{ */
#include "eg_lpnum.h"
/* ========================================================================= */
/** @brief Tester program for EGlpNum_t structure and functions
 * @return zero on success, non-zero otherwise 
 * @par Description:
 * Perform various tests on EGlpNum_t and their functions 
 * */
int main (int argc,
					char **argv)
{
	/* local variables */
	EGlpNum_t ntmp[4];
	double dtmp[5];
	char *strnum1,
	 *strnum2,
	 *strnum3;
	int rval=0;
	#ifdef HAVE_LIBGMP
	int n_char = 0;
	mpq_t qnum;
	#endif
	EGlpNumStart();
	/* set signal and limits */
	EGsigSet(rval,CLEANUP);
	EGsetLimits(3600.0,4294967295UL);
	#ifdef HAVE_LIBGMP
	EGlpNumSetPrecision (128);
	mpq_init (qnum);
	#endif
	EGlpNumInitVar (ntmp[0]);
	EGlpNumInitVar (ntmp[1]);
	EGlpNumInitVar (ntmp[2]);
	EGlpNumInitVar (ntmp[3]);
	/* the input should have at least two parameters, namelly the number where we
	 * will work on */
	if (argc < 3)
	{
		fprintf (stderr,
						 "usage: %s num1 num2\n\tWhere num1 and num2 are numbers in"
						 " the number format (either a/b or regular doubles)\n", argv[0]);
		exit (1);
	}

	/* we ask for two numbers and perform some basic operations and compare
	 * aganist double arithmetic */
	#ifdef HAVE_LIBGMP
	n_char = mpq_EGlpNumReadStrXc (qnum, argv[1]);
	#endif
	EGlpNumReadStr (ntmp[0], argv[1]);
	EGlpNumReadStr (ntmp[1], argv[2]);
	dtmp[0] = EGlpNumToLf (ntmp[0]);
	dtmp[1] = EGlpNumToLf (ntmp[1]);

	/* convert numbers */
	strnum1 = EGlpNumGetStr (ntmp[0]);
	strnum2 = EGlpNumGetStr (ntmp[1]);
	fprintf (stderr, "You Input %s (%lg) and %s (%lg)\n", strnum1, dtmp[0],
					 strnum2, dtmp[1]);
	#ifdef HAVE_LIBGMP
	strnum3 = mpq_EGlpNumGetStr (qnum);
	fprintf (stderr, "Your first number represented as exact rational is %s, "
					 "readed with %d chars\n", strnum3, n_char);
	free (strnum3);
	mpq_EGlpNumSet (qnum, dtmp[0]);
	strnum3 = mpq_EGlpNumGetStr (qnum);
	fprintf (stderr, "Your first number represented as continuous fraction "
					 "is %s, readed with %d chars\n", strnum3, n_char);
	free (strnum3);
	#endif

	/* internal constants */
	strnum3 = EGlpNumGetStr (oneLpNum);
	fprintf (stderr, "1.0 = %s\n", strnum3);
	EGfree (strnum3);
	strnum3 = EGlpNumGetStr (zeroLpNum);
	fprintf (stderr, "0.0 = %s\n", strnum3);
	EGfree (strnum3);
	strnum3 = EGlpNumGetStr (epsLpNum);
	fprintf (stderr, "eps = %s\n", strnum3);
	EGfree (strnum3);
	strnum3 = EGlpNumGetStr (MaxLpNum);
	fprintf (stderr, "Max = %s\n", strnum3);
	EGfree (strnum3);
	strnum3 = EGlpNumGetStr (MinLpNum);
	fprintf (stderr, "Min = %s\n", strnum3);
	EGfree (strnum3);

	/* copying functions */
	EGlpNumCopy (ntmp[2], ntmp[0]);
	strnum3 = EGlpNumGetStr (ntmp[2]);
	fprintf (stderr, "%s = %s (%lg)\n", strnum3, strnum1, dtmp[0]);
	EGfree (strnum3);
	EGlpNumCopyDiff (ntmp[2], ntmp[0], ntmp[1]);
	strnum3 = EGlpNumGetStr (ntmp[2]);
	fprintf (stderr, "%s = %s - %s (%lg)\n", strnum3, strnum1, strnum2,
					 dtmp[0] - dtmp[1]);
	EGfree (strnum3);
	EGlpNumCopyDiffRatio (ntmp[2], ntmp[0], ntmp[1], ntmp[0]);
	strnum3 = EGlpNumGetStr (ntmp[2]);
	fprintf (stderr, "%s = (%s  - %s)/%s (%lg)\n", strnum3, strnum1, strnum2,
					 strnum1, (dtmp[0] - dtmp[1]) / dtmp[0]);
	EGfree (strnum3);
	EGlpNumCopySum (ntmp[2], ntmp[0], ntmp[1]);
	strnum3 = EGlpNumGetStr (ntmp[2]);
	fprintf (stderr, "%s = %s + %s (%lg)\n", strnum3, strnum1, strnum2,
					 dtmp[0] + dtmp[1]);
	EGfree (strnum3);
	EGlpNumCopySqrOver (ntmp[2], ntmp[1], ntmp[0]);
	strnum3 = EGlpNumGetStr (ntmp[2]);
	fprintf (stderr, "%s = %s^2/%s (%lg)\n", strnum3, strnum2, strnum1,
					 dtmp[1] * dtmp[1] / dtmp[0]);
	EGfree (strnum3);
	EGlpNumCopyAbs (ntmp[2], ntmp[0]);
	strnum3 = EGlpNumGetStr (ntmp[2]);
	fprintf (stderr, "%s = |%s| (%lg)\n", strnum3, strnum1, fabs (dtmp[0]));
	EGfree (strnum3);
	EGlpNumCopyNeg (ntmp[2], ntmp[0]);
	strnum3 = EGlpNumGetStr (ntmp[2]);
	fprintf (stderr, "%s = -1*%s (%lg)\n", strnum3, strnum1, -dtmp[0]);
	EGfree (strnum3);
	EGlpNumCopyFrac (ntmp[2], ntmp[1], ntmp[0]);
	strnum3 = EGlpNumGetStr (ntmp[2]);
	fprintf (stderr, "%s = %s/%s (%lg)\n", strnum3, strnum2, strnum1,
					 dtmp[1] / dtmp[0]);
	EGfree (strnum3);

	/* add */
	EGlpNumCopy (ntmp[2], ntmp[0]);
	EGlpNumAddTo (ntmp[2], ntmp[1]);
	strnum3 = EGlpNumGetStr (ntmp[2]);
	fprintf (stderr, "%s + %s = %s (%lg)\n", strnum1, strnum2, strnum3,
					 dtmp[0] + dtmp[1]);
	EGfree (strnum3);
	EGlpNumCopy (ntmp[2], ntmp[0]);
	EGlpNumAddUiTo (ntmp[2], 0xffU);
	strnum3 = EGlpNumGetStr (ntmp[2]);
	fprintf (stderr, "%s + %u = %s (%lg)\n", strnum1, 0xffU, strnum3,
					 dtmp[0] + 0xffU);
	EGfree (strnum3);

	/* substract */
	EGlpNumCopy (ntmp[2], ntmp[0]);
	EGlpNumSubTo (ntmp[2], ntmp[1]);
	strnum3 = EGlpNumGetStr (ntmp[2]);
	fprintf (stderr, "%s - %s = %s (%lg)\n", strnum1, strnum2, strnum3,
					 dtmp[0] - dtmp[1]);
	EGfree (strnum3);
	EGlpNumCopy (ntmp[2], ntmp[0]);
	EGlpNumSubUiTo (ntmp[2], 0xffU);
	strnum3 = EGlpNumGetStr (ntmp[2]);
	fprintf (stderr, "%s - %u = %s (%lg)\n", strnum1, 0xffU, strnum3,
					 dtmp[0] - 0xffU);
	EGfree (strnum3);

	/* multiply */
	EGlpNumCopy (ntmp[2], ntmp[0]);
	EGlpNumMultTo (ntmp[2], ntmp[1]);
	strnum3 = EGlpNumGetStr (ntmp[2]);
	fprintf (stderr, "%s * %s = %s (%lg)\n", strnum1, strnum2, strnum3,
					 dtmp[0] * dtmp[1]);
	EGfree (strnum3);

	/* multiply unsigned */
	EGlpNumCopy (ntmp[2], ntmp[0]);
	EGlpNumMultUiTo (ntmp[2], 13);
	strnum3 = EGlpNumGetStr (ntmp[2]);
	fprintf (stderr, "%s * 13 = %s (%lg)\n", strnum1, strnum3, dtmp[0] * 13);
	EGfree (strnum3);

	/* inverse */
	EGlpNumCopy (ntmp[2], ntmp[0]);
	EGlpNumInv (ntmp[2]);
	strnum3 = EGlpNumGetStr (ntmp[2]);
	fprintf (stderr, "1/%s = %s (%lg)\n", strnum1, strnum3, 1.0 / dtmp[0]);
	EGfree (strnum3);

	/* floor and ceil */
	EGlpNumFloor (ntmp[2], ntmp[0]);
	strnum3 = EGlpNumGetStr (ntmp[2]);
	fprintf (stderr, "floor(%s) = %s (%lg)\n", strnum1, strnum3, floor (dtmp[0]));
	EGfree (strnum3);
	EGlpNumCeil (ntmp[2], ntmp[0]);
	strnum3 = EGlpNumGetStr (ntmp[2]);
	fprintf (stderr, "ceil(%s) = %s (%lg)\n", strnum1, strnum3, ceil (dtmp[0]));
	EGfree (strnum3);

	/* negative and inner products */
	EGlpNumCopy (ntmp[2], ntmp[0]);
	EGlpNumSign (ntmp[2]);
	strnum3 = EGlpNumGetStr (ntmp[2]);
	fprintf (stderr, "-(%s) = %s (%lg)\n", strnum1, strnum3, -dtmp[0]);
	EGfree (strnum3);
	EGlpNumOne (ntmp[2]);
	EGlpNumAddInnProdTo (ntmp[2], ntmp[0], ntmp[1]);
	strnum3 = EGlpNumGetStr (ntmp[2]);
	fprintf (stderr, "1.0 + %s*%s = %s (%lg)\n", strnum1, strnum2, strnum3,
					 1.0 + dtmp[1] * dtmp[0]);
	EGfree (strnum3);
	EGlpNumOne (ntmp[2]);
	EGlpNumSubInnProdTo (ntmp[2], ntmp[0], ntmp[1]);
	strnum3 = EGlpNumGetStr (ntmp[2]);
	fprintf (stderr, "1.0 - %s*%s = %s (%lg)\n", strnum1, strnum2, strnum3,
					 1.0 - dtmp[1] * dtmp[0]);
	EGfree (strnum3);

	/* divide */
	EGlpNumCopy (ntmp[2], ntmp[0]);
	EGlpNumDivTo (ntmp[2], ntmp[1]);
	strnum3 = EGlpNumGetStr (ntmp[2]);
	fprintf (stderr, "%s / %s = %s (%lg)\n", strnum1, strnum2, strnum3,
					 dtmp[0] / dtmp[1]);
	EGfree (strnum3);
	EGlpNumCopy (ntmp[2], ntmp[0]);
	EGlpNumDivUiTo (ntmp[2], 0xffU);
	strnum3 = EGlpNumGetStr (ntmp[2]);
	fprintf (stderr, "%s / %u = %s (%lg)\n", strnum1, 0xffU, strnum3,
					 dtmp[0] / 0xffU);
	EGfree (strnum3);
	EGfree (strnum1);
	EGfree (strnum2);

#ifdef HAVE_LIBGMP
	/* test transformation to rationals */
	{
		mpq_t n1, n2;
		mpf_t f1, f2;
		mpq_init (n1);
		mpq_init (n2);
		mpf_init (f1);
		mpf_init (f2);
		mpf_EGlpNumSet(f1,EGlpNumToLf(ntmp[0]));
		mpf_EGlpNumSet(f2,EGlpNumToLf(ntmp[1]));
		mpq_EGlpNumSet_mpf (n1, f1);
		mpq_EGlpNumSet_mpf (n2, f2);
		strnum1 = mpq_EGlpNumGetStr (n1);
		strnum2 = mpq_EGlpNumGetStr (n2);
		fprintf (stderr, "Your input in rational was:\n\t(%10.7lf) %s\n\t(%10.7lf)"
				" %s\n", EGlpNumToLf (ntmp[0]), strnum1, EGlpNumToLf (ntmp[1]),
				strnum2);
		EGfree (strnum1);
		EGfree (strnum2);

		/* test natural exponentiation */
		mpf_EGlpNumEpow (f1, -38.81624211135693732736499880);
		mpf_set_ui (f2, (unsigned long)1);
		mpf_div_2exp (f2, f2, (unsigned long)56);
		strnum1 = mpf_EGlpNumGetStr (f1);
		strnum2 = mpf_EGlpNumGetStr (f2);
		fprintf (stderr, "2^-56 = %s ~ %s\n", strnum1, strnum2);
		EGfree (strnum1);
		EGfree (strnum2);
		mpq_clear (n1);
		mpq_clear (n2);
		mpf_clear (f1);
		mpf_clear (f2);
	}
#endif

	/* ending */
	CLEANUP:
	EGlpNumClearVar (ntmp[0]);
	EGlpNumClearVar (ntmp[1]);
	EGlpNumClearVar (ntmp[2]);
	EGlpNumClearVar (ntmp[3]);
	#ifdef HAVE_LIBGMP
	mpq_clear (qnum);
	#endif
	EGlpNumClear();
	return 0;
}

/* ========================================================================= */
/** @} */
