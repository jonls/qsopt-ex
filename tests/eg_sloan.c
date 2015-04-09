
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include "qs_config.h"
#include "logging-private.h"

#include "QSopt_ex.h"
/* ========================================================================= */
/** @name Static Variables
 * Set of static variables for this program. */
/*@{*/
/** @brief range of values for first OA level */
static unsigned s1 = 3;
/** @brief range of values for second OA level */
static unsigned s2 = 3;
/** @brief number of columns for first OA level */
static unsigned k1 = 3;
/** @brief number of columns for second OA level */
static unsigned k2 = 3;
/** @brief use (or not) scaling of the resulting LP */
static int use_scaling = 1;
/** @brief use dlouble floating point output */
static int use_double = 0;
/** @brief file name output. */
const char *out_file = "output.lp";
/** @brief Strength level required */
static unsigned t = 2;
/** @brief temporal string space */
char strtmp[1024];
/*@}*/

/* ========================================================================= */
/** @brief Display options to the screen */
static inline void sl_usage (char const *const s)
{
	fprintf (stderr,
					 "This programs compute bounds for mixed-level orthogonal arrays (OA) of different strengths for two levels\n");
	fprintf (stderr, "Usage: %s [options]\n", s);
	fprintf (stderr, "    -a n   range of values for the first OA level (>=2)\n");
	fprintf (stderr,
					 "    -b n   range of values for the second OA level (>=2)\n");
	fprintf (stderr,
					 "    -c n   number of columns for the first OA level (>=2)\n");
	fprintf (stderr,
					 "    -e n   number of columns for the second OA level (>=2)\n");
	fprintf (stderr,
					 "    -d n   if n > 0, write the LP in double precition arithmetic, otherwise, write it in rational form\n");
	fprintf (stderr, "    -o f   filename where to write the output LP\n");
	fprintf (stderr,
					 "    -s n   if n > 0 scale the resulting LP, otherwise present the unscaled LP\n");
	fprintf (stderr, "    -t n   required strength of the OA (>=1)\n");
}

/* ========================================================================= */
/** @brief Display options to the screen */
static inline int sl_parseargs (int argc,
																char **argv)
{
	int c;
	while ((c = getopt (argc, argv, "a:b:c:e:s:d:o:t:")) != EOF)
	{
		switch (c)
		{
		case 'a':
			s1 = atoi (optarg);
			break;
		case 'b':
			s2 = atoi (optarg);
			break;
		case 'c':
			k1 = atoi (optarg);
			break;
		case 'e':
			k2 = atoi (optarg);
			break;
		case 's':
			use_scaling = atoi (optarg);
			break;
		case 'd':
			use_double = atoi (optarg);
			break;
		case 'o':
			out_file = optarg;
			break;
		case 't':
			t = atoi (optarg);
			break;
		default:
			sl_usage (argv[0]);
			return 1;
		}
	}
	if (s1 < 2 || s2 < 2 || k1 < 1 || k2 < 1 || t < 1)
	{
		sl_usage (argv[0]);
		return 1;
	}
	fprintf (stderr, "Running %s\nOptions:\n", argv[0]);
	fprintf (stderr, "\tfirst level columns = %d\n", k1);
	fprintf (stderr, "\tsecond level columns = %d\n", k2);
	fprintf (stderr, "\tfirst level range = %d\n", s1);
	fprintf (stderr, "\tsecond level range = %d\n", s2);
	fprintf (stderr, "\tt = %d\n", t);
	fprintf (stderr, "\t%s scaling\n", use_scaling ? "using" : "not using");
	fprintf (stderr, "\t%s output\n", use_double ? "double" : "rational");
	fprintf (stderr, "\toutput file : %s\n", out_file);
	return 0;
}

/* ========================================================================= */
/** @brief compute the krawtchouk ppolynomial, defined as
 * \f[P^s_j(x,m)=\sum\limits_{i=0}^j(-1)^i(s-1)^{j-i}\binom{x}{i}\binom{m-x}{j-i}\f]
 * @param x the \f$x\f$ parameter of the function.
 * @param s the \f$s\f$ parameter of the function.
 * @param j the \f$j\f$ parameter of the function.
 * @param m the \f$m\f$ parameter of the function.
 * @param rop where to store the result */
void Kpoly (unsigned x,
						unsigned s,
						unsigned j,
						unsigned m,
						mpz_t rop)
{
	register unsigned i = j + 1;
	mpz_t z1,
	  z2;
	mpz_init (z1);
	mpz_init (z2);
	mpz_set_ui (rop, (unsigned long)0);
	EXIT (m < x, "m < x! impossible!");
	EXIT (s < 1, "s < 1! impossible!");
	while (i--)
	{
		mpz_bin_uiui (z1, (unsigned long)x, (unsigned long)i);
		mpz_set (z2, z1);
		mpz_bin_uiui (z1, (unsigned long)(m - x), (unsigned long)(j - i));
		mpz_mul (z2, z2, z1);
		mpz_ui_pow_ui (z1, (unsigned long)(s - 1), (unsigned long)(j - i));
		mpz_mul (z2, z2, z1);
		if (i & 1U)
			mpz_sub (rop, rop, z2);
		else
			mpz_add (rop, rop, z2);
	}
	mpz_clear (z1);
	mpz_clear (z2);
}

/* ========================================================================= */
/** @brief main function, here we build the LP, execute the options and exit */
int main (int argc,
					char **argv)
{
	int rval = 0;
	mpq_t v1,
	  v2;
	mpq_QSdata *p_mpq = 0;
	dbl_QSdata *p_dbl = 0;
	register unsigned i,
	  j,
	  k,
	  l;
	/* parse input */
	QSopt_ex_version();
	QSexactStart();
	rval = sl_parseargs (argc, argv);
	if (rval)
		return rval;
	mpq_init (v1);
	mpq_init (v2);
	/* create the problem with the appropriate number of variables and
	 * constraints */
	snprintf (strtmp, (size_t)1023, "OA_SL_%d-%d_%d-%d_t-%d", s1, k1, s2, k2, t);
	p_mpq = mpq_QScreate_prob (strtmp, QS_MIN);
	for (i = 0; i <= k1; i++)
		for (j = 0; j <= k2; j++)
		{
			snprintf (strtmp, (size_t)1023, "V_%d_%d", i, j);
			rval =
				mpq_QSnew_col (p_mpq, mpq_oneLpNum,
											 (i + j == 0) ? mpq_oneLpNum : mpq_zeroLpNum,
											 mpq_ILL_MAXDOUBLE, strtmp);
			CHECKRVALG (rval, CLEANUP);
			strtmp[0] = 'C';
			rval = mpq_QSnew_row (p_mpq, mpq_zeroLpNum,
											 ((i + j >= 1) && (i + j <= t)) ? 'E' : 'G', strtmp);
			CHECKRVALG (rval, CLEANUP);
		}
	/* now set the coefficients */
	for (i = 0; i <= k1; i++)
		for (j = 0; j <= k2; j++)
			for (k = 0; k <= k1; k++)
				for (l = 0; l <= k2; l++)
				{
					Kpoly (k, s1, i, k1, mpq_numref (v2));
					Kpoly (l, s2, j, k2, mpq_numref (v1));
					mpq_mul (v1, v1, v2);
					if (mpz_cmp_ui (mpq_numref (v1), (unsigned long)0))
					{
						rval =
							mpq_QSchange_coef (p_mpq, ((int)(j + i * (k2 + 1))), (int)(l + k * (k2 + 1)), v1);
						CHECKRVALG (rval, CLEANUP);
					}
				}
	/* now, if we are using scaling, we scale the non-zeros */
	if (use_scaling)
	{
		mpq_set_ui (v1, (unsigned long)1, (unsigned long)1);
		EXutilSimplify ((unsigned)(p_mpq->qslp->A.matsize), p_mpq->qslp->A.matval, v1);
		mpq_div (v2, mpq_oneLpNum, v1);
		fprintf (stderr, "Scale factor %lf\n", mpq_get_d (v2));
	}
	/* now we save the LP */
	if (use_double)
	{
		snprintf (strtmp, (size_t)1023, "OA_SL_%d-%d_%d-%d_t-%d", s1, k1, s2, k2, t);
		p_dbl = QScopy_prob_mpq_dbl (p_mpq, strtmp);
		rval = dbl_QSwrite_prob (p_dbl, out_file, "LP");
		CHECKRVALG (rval, CLEANUP);
	}
	else
	{
		rval = mpq_QSwrite_prob (p_mpq, out_file, "LP");
		CHECKRVALG (rval, CLEANUP);
	}

	/* ending */
CLEANUP:
	if (p_mpq)
		mpq_QSfree_prob (p_mpq);
	if (p_dbl)
		dbl_QSfree_prob (p_dbl);
	mpq_clear (v1);
	mpq_clear (v2);
	QSexactClear();
	return rval;
}
