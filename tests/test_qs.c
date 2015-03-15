/*
 * This file is part of QSopt_ex.
 *
 * (c) Copyright 2006 by David Applegate, William Cook, Sanjeeb Dash,
 * and Daniel Espinoza
 * (c) Copyright 2015  Jon Lund Steffensen <jonlst@gmail.com>
 *
 * Sanjeeb Dash ownership of copyright in QSopt_ex is derived from his
 * copyright in QSopt.
 *
 * This code may be used under the terms of the GNU General Public License
 * (Version 2.1 or later) as published by the Free Software Foundation.
 *
 * Alternatively, use is granted for research purposes only.
 *
 * It is your choice of which of these two licenses you are operating
 * under.
 *
 * We make no guarantees about the correctness or usefulness of this code.
 */

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#include "QSopt_ex.h"

/* Number of tests. Must be kept up-to-date manually */
#define TEST_COUNT  8

/*                                                      */
/*  Using QSload_prob to load the following LP problem  */
/*       Maximize  3.0x + 2.0y + 4.0z                   */
/*       Subject to                                     */
/*                 3.0x + 2.0y + 1.0z <= 12.0           */
/*                 5.0x + 1.0y         = 10.0           */
/*                 x >= 2.0                             */
/*                 y free                               */
/*                 1.0 <= z <= 10.0                     */
/*                                                      */

static int load_test_problem (mpq_QSprob *p)
{
    int i, rval = 0;
    int cmatcnt[3] = { 2, 2, 1 };
    int cmatbeg[3] = { 0, 2, 4 };
    int cmatind[5] = { 0, 1, 0, 1, 0 };
    char sense[2] = { 'L', 'E' };
    const char *colnames[3] = { "x", "y", "z" };
    const char *rownames[2] = { "c1", "c2"};
    mpq_t cmatval[5];
    mpq_t obj[3];
    mpq_t rhs[2];
    mpq_t lower[3];
    mpq_t upper[3];

    for (i = 0; i < 5; i++) mpq_init (cmatval[i]);
    mpq_set_d (cmatval[0], 3.0);
    mpq_set_d (cmatval[1], 5.0);
    mpq_set_d (cmatval[2], 2.0);
    mpq_set_d (cmatval[3], 1.0);
    mpq_set_d (cmatval[4], 1.0);

    for (i = 0; i < 3; i++) mpq_init (obj[i]);
    mpq_set_d (obj[0], 3.0);
    mpq_set_d (obj[1], 2.0);
    mpq_set_d (obj[2], 4.0);

    for (i = 0; i < 2; i++) mpq_init (rhs[i]);
    mpq_set_d (rhs[0], 12.0);
    mpq_set_d (rhs[1], 10.0);

    for (i = 0; i < 3; i++) mpq_init (lower[i]);
    mpq_set_d (lower[0], 2.0);
    mpq_set (lower[1], mpq_ILL_MINDOUBLE);
    mpq_set_d (lower[2], 1.0);

    for (i = 0; i < 3; i++) mpq_init (upper[i]);
    mpq_set (upper[0], mpq_ILL_MAXDOUBLE);
    mpq_set (upper[1], mpq_ILL_MAXDOUBLE);
    mpq_set_d (upper[2], 10.0);

    /*  CPXcopylpwnames  */

    *p = mpq_QSload_prob ("small", 3, 2, cmatcnt, cmatbeg, cmatind, cmatval,
                      QS_MAX, obj, rhs, sense, lower, upper, colnames,
                      rownames);

    if (*p == (mpq_QSprob) NULL) {
        fprintf (stderr, "Unable to load the LP problem\n");
        rval = 1;  goto CLEANUP;
    }

CLEANUP:

    for (i = 0; i < 5; i++) mpq_clear (cmatval[i]);
    for (i = 0; i < 3; i++) mpq_clear (obj[i]);
    for (i = 0; i < 2; i++) mpq_clear (rhs[i]);
    for (i = 0; i < 3; i++) mpq_clear (lower[i]);
    for (i = 0; i < 3; i++) mpq_clear (upper[i]);

    return rval;
}

int main (int ac, char **av)
{
    int rval = 0, status = 0;
    int i, nrows = 0, ncols = 0;
    mpq_t value;
    mpq_t *x = (mpq_t *) NULL;
    mpq_t *y = (mpq_t *) NULL;
    mpq_QSprob p = (mpq_QSprob) NULL;
    char **colnames = (char **) NULL;
    char **rownames = (char **) NULL;

    QSexactStart();

    mpq_init (value);

    if (ac > 1) {
        printf ("Usage: %s\n", *av);
        rval = EXIT_FAILURE; goto CLEANUP;
    }

    /*  Must call QSexact_set_precision before any other QS function       */
    /*  Cannot use mpq_ILL_MAXDOUBLE or mpq_ILL_MINDOUBLE before this call */

    QSexact_set_precision (128);   

    /* Print number of tests */
    printf("1..%u\n", TEST_COUNT);

    /* Test 1 - load test problem */
    rval = load_test_problem (&p);
    if (rval) {
        printf("not ok 1 - Unable to load the LP\n");
        goto CLEANUP;
    }

    printf("ok 1 - Loaded test LP problem\n");

    /* Test 2 - get row count */
    nrows = mpq_QSget_rowcount (p);
    if (nrows != 2) {
        printf("not ok 2 - Invalid number of constraints: %d\n", nrows);
    } else {
        printf("ok 2 - Got number of constraints\n");
    }

    /* Test 3 - get col count */
    ncols = mpq_QSget_colcount (p);
    if (ncols != 3) {
        printf("not ok 3 - Invalid number of variables: %d\n", ncols);
    } else {
        printf("ok 3 - Got number of variables\n");
    }

    /*  Test 4 - Solve problem  */
    rval = QSexact_solver (p, NULL, NULL, NULL, DUAL_SIMPLEX, &status);
    if (rval) {
        printf("not ok 4 - Solver failed\n");
        goto CLEANUP;
    }
    if (status != QS_LP_OPTIMAL) {
        printf("not ok 4 - Did not find an optimal solution.\n");
        goto CLEANUP;
    }

    printf("ok 4 - Solution was found\n");

    /* Test 5 - Get objective value  */
    rval = mpq_QSget_objval (p, &value);
    if (rval) {
        printf("not ok 5 - Could not get obj value, error code %d\n", rval);
    } else {
        if (mpq_cmp_ui(value, 42, 1) != 0) {
            printf("not ok 5 - Unexpected obj value: %.6f\n", mpq_get_d(value));
        } else {
            printf("ok 5 - The correct objective value was obtained\n");
        }
    }

    /* Test 6 - Get variable values */
    x = (mpq_t *) malloc (ncols * sizeof (mpq_t));
    for (i = 0; i < ncols; i++) mpq_init (x[i]);
    rval = mpq_QSget_x_array (p, x);
    if (rval) {
        printf("not ok 6 - Could not get x-vector, error code %d\n", rval);
    } else {
        colnames = (char **) malloc (ncols * sizeof (char *));
        mpq_QSget_colnames (p, colnames);
        printf("ok 6 - Obtained variable values and names\n");
    }

    /* Test 7 - Get dual values */
    y = (mpq_t *) malloc (nrows * sizeof (mpq_t));
    for (i = 0; i < nrows; i++) mpq_init (y[i]);
    rval = mpq_QSget_pi_array (p, y);
    if (rval) {
        printf("not ok 7 - Could not get dual values, error code %d\n", rval);
    } else {
        rownames = (char **) malloc (nrows * sizeof (char *));
        mpq_QSget_rownames (p, rownames);
        printf("ok 7 - Obtained dual values and contraint names\n");
    }

    /* Test 8 - Write problem out to /dev/null */
    rval = mpq_QSwrite_prob (p, "/dev/null", "LP");
    if (rval) {
        printf("not ok 8 - Could not write the LP, error code %d\n", rval);
    } else {
        printf("ok 8 - LP written to /dev/null\n");
    }

CLEANUP:

    if (p) mpq_QSfree_prob (p);  /*  CPXfreeprob  */

    mpq_clear (value);
    if (x) {
        for (i = 0; i < ncols; i++) mpq_clear (x[i]);
        free (x);
    }
    if (y) {
        for (i = 0; i < nrows; i++) mpq_clear (y[i]);
        free (y);
    }
    if (rownames) {
        for (i = 0; i < nrows; i++) if (rownames[i]) free (rownames[i]);
        free (rownames);
    }
    if (colnames) {
        for (i = 0; i < ncols; i++) if (colnames[i]) free (colnames[i]);
        free (colnames);
    }

    QSexactClear();
    return rval;
}
