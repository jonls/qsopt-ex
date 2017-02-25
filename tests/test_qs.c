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

typedef void test_func(int test_id);


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
static int load_test_problem(mpq_QSprob *p)
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

    if (*p == NULL) {
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

/* Load and solve test problem. */
static int solve_test_problem(mpq_QSprob *p, int *status)
{
    int rval = 0;

    rval = load_test_problem(p);
    if (rval) return rval;

    rval = QSexact_solver(*p, NULL, NULL, NULL, DUAL_SIMPLEX, status);
    if (rval) return rval;

    return rval;
}

static void test_new_column(int test_id)
{
    mpq_t objective, lower, upper;
    mpq_init(objective);
    mpq_init(lower);
    mpq_init(upper);

    mpq_set_d(objective, 4.0);
    mpq_set_d(lower, -5.0);
    mpq_set_d(upper, 75.5);

    const char *name = "test_var";

    mpq_QSprob p = mpq_QScreate_prob("test", QS_MAX);
    if (p == NULL) {
        printf("not ok %i - Unable to create LP problem\n", test_id);
        goto CLEANUP;
    }

    /* Test create new column */
    int rval = mpq_QSnew_col(p, objective, lower, upper, name);
    if (rval) {
        printf("not ok %i - Failed to create variable\n", test_id);
    } else {
        printf("ok %i - New variable created\n", test_id);
    }

CLEANUP:
    mpq_clear(objective);
    mpq_clear(lower);
    mpq_clear(upper);
    if (p) mpq_QSfree_prob(p);
}

static void test_new_column_without_name(int test_id)
{
    mpq_t objective, lower, upper;
    mpq_init(objective);
    mpq_init(lower);
    mpq_init(upper);

    mpq_set_d(objective, 4.0);
    mpq_set_d(lower, -5.0);
    mpq_set_d(upper, 75.5);

    mpq_QSprob p = mpq_QScreate_prob("test", QS_MAX);
    if (p == NULL) {
        printf("not ok %i - Unable to create LP problem\n", test_id);
        goto CLEANUP;
    }

    /* Test create new column */
    int rval = mpq_QSnew_col(p, objective, lower, upper, NULL);
    if (rval) {
        printf("not ok %i - Failed to create variable\n", test_id);
    } else {
        printf("ok %i - New variable created without name\n", test_id);
    }

CLEANUP:
    mpq_clear(objective);
    mpq_clear(lower);
    mpq_clear(upper);
    if (p) mpq_QSfree_prob(p);
}

static void test_delete_column(int test_id)
{
    mpq_QSprob p = NULL;

    int rval = load_test_problem(&p);
    if (rval) {
        printf("not ok %i - Unable to load the LP\n", test_id);
        goto CLEANUP;
    }

    /* Test delete column */
    rval = mpq_QSdelete_col(p, 0);
    if (rval) {
        printf("not ok %i - Failed to delete variable\n", test_id);
    } else {
        printf("ok %i - Deleted varable\n", test_id);
    }

CLEANUP:
    if (p) mpq_QSfree_prob(p);
}

static void test_delete_invalid_column(int test_id)
{
    mpq_QSprob p = NULL;

    int rval = load_test_problem(&p);
    if (rval) {
        printf("not ok %i - Unable to load the LP\n", test_id);
        goto CLEANUP;
    }

    /* Test delete column */
    rval = mpq_QSdelete_col(p, 10);
    if (rval) {
        printf("ok %i - Error reported when deleting invalid variable\n",
               test_id);
    } else {
        printf("not ok %i - No error reported\n", test_id);
    }

CLEANUP:
    if (p) mpq_QSfree_prob(p);
}

static void test_delete_column_in_empty_problem(int test_id)
{
    mpq_QSprob p = mpq_QScreate_prob("test", QS_MAX);
    if (p == NULL) {
        printf("not ok %i - Unable to create LP problem\n", test_id);
        goto CLEANUP;
    }

    /* Test delete column */
    int rval = mpq_QSdelete_col(p, 0);
    if (rval) {
        printf("ok %i - Error reported when deleting invalid variable\n",
               test_id);
    } else {
        printf("not ok %i - No error reported\n", test_id);
    }

CLEANUP:
    if (p) mpq_QSfree_prob(p);
}

static void test_delete_row(int test_id)
{
    mpq_QSprob p = NULL;

    int rval = load_test_problem(&p);
    if (rval) {
        printf("not ok %i - Unable to load the LP\n", test_id);
        goto CLEANUP;
    }

    /* Test delete row */
    rval = mpq_QSdelete_row(p, 1);
    if (rval) {
        printf("not ok %i - Failed to delete constraint\n", test_id);
    } else {
        printf("ok %i - Deleted constraint\n", test_id);
    }

CLEANUP:
    if (p) mpq_QSfree_prob(p);
}

static void test_delete_invalid_row(int test_id)
{
    mpq_QSprob p = NULL;

    int rval = load_test_problem(&p);
    if (rval) {
        printf("not ok %i - Unable to load the LP\n", test_id);
        goto CLEANUP;
    }

    /* Test delete row */
    rval = mpq_QSdelete_row(p, 5);
    if (rval) {
        printf("ok %i - Error reported when deleting invalid constraint\n",
               test_id);
    } else {
        printf("not ok %i - No error reported\n", test_id);
    }

CLEANUP:
    if (p) mpq_QSfree_prob(p);
}

static void test_delete_row_in_empty_problem(int test_id)
{
    mpq_QSprob p = mpq_QScreate_prob("test", QS_MAX);
    if (p == NULL) {
        printf("not ok %i - Unable to create LP problem\n", test_id);
        goto CLEANUP;
    }

    /* Test delete row */
    int rval = mpq_QSdelete_row(p, 0);
    if (rval) {
        printf("ok %i - Error reported when deleting invalid constraint\n",
               test_id);
    } else {
        printf("not ok %i - No error reported\n", test_id);
    }

CLEANUP:
    if (p) mpq_QSfree_prob(p);
}

static void test_row_count(int test_id)
{
    mpq_QSprob p = NULL;

    int rval = load_test_problem(&p);
    if (rval) {
        printf("not ok %i - Unable to load the LP\n", test_id);
        goto CLEANUP;
    }

    /* Test get row count */
    int nrows = mpq_QSget_rowcount (p);
    if (nrows != 2) {
        printf("not ok %i - Invalid number of constraints: %d\n",
               test_id, nrows);
    } else {
        printf("ok %i - Got number of constraints\n", test_id);
    }

CLEANUP:
    if (p) mpq_QSfree_prob(p);
}

static void test_col_count(int test_id)
{
    mpq_QSprob p = NULL;

    int rval = load_test_problem(&p);
    if (rval) {
        printf("not ok %i - Unable to load the LP\n", test_id);
        goto CLEANUP;
    }

    /* Test get col count */
    int ncols = mpq_QSget_colcount(p);
    if (ncols != 3) {
        printf("not ok %i - Invalid number of variables: %d\n",
               test_id, ncols);
    } else {
        printf("ok %i - Got number of variables\n", test_id);
    }

CLEANUP:
    if (p) mpq_QSfree_prob(p);
}

static void test_solution_is_optimal(int test_id)
{
    mpq_QSprob p = NULL;
    int status = 0;

    int rval = solve_test_problem(&p, &status);
    if (rval) {
        printf("not ok %i - Unable to solve the LP\n", test_id);
        goto CLEANUP;
    }

    /* Test status */
    if (status != QS_LP_OPTIMAL) {
        printf("not ok %i - Did not find an optimal solution.\n", test_id);
        goto CLEANUP;
    }

    printf("ok %i - Solution was found\n", test_id);

CLEANUP:
    if (p) mpq_QSfree_prob(p);
}

static void test_solution_objective(int test_id)
{
    mpq_QSprob p = NULL;
    int status = 0;
    int rval = 0;

    mpq_t value;
    mpq_init(value);

    rval = solve_test_problem(&p, &status);
    if (rval) {
        printf("not ok %i - Unable to solve the LP\n", test_id);
        goto CLEANUP;
    }

    if (status != QS_LP_OPTIMAL) {
        printf("not ok %i - Did not find an optimal solution.\n", test_id);
        goto CLEANUP;
    }

    /* Test objective */
    rval = mpq_QSget_objval(p, &value);
    if (rval) {
        printf("not ok %i - Could not get obj value, error code %d\n",
               test_id, rval);
    } else {
        if (mpq_cmp_ui(value, 42, 1) != 0) {
            printf("not ok %i - Unexpected obj value: %.6f\n", test_id,
                   mpq_get_d(value));
        } else {
            printf("ok %i - The correct objective value was obtained\n",
                   test_id);
        }
    }

CLEANUP:
    if (p) mpq_QSfree_prob(p);
    mpq_clear(value);
}

static void test_solution_get_variables(int test_id)
{
    mpq_QSprob p = NULL;
    mpq_t *x = NULL;
    int status = 0;
    int rval = 0;
    int i;
    int ncols = 3;

    rval = solve_test_problem(&p, &status);
    if (rval) {
        printf("not ok %i - Unable to solve the LP\n", test_id);
        goto CLEANUP;
    }

    if (status != QS_LP_OPTIMAL) {
        printf("not ok %i - Did not find an optimal solution.\n", test_id);
        goto CLEANUP;
    }

    /* Test get variables */
    x = malloc(ncols * sizeof(mpq_t));
    for (i = 0; i < ncols; i++) mpq_init(x[i]);
    rval = mpq_QSget_x_array(p, x);
    if (rval) {
        printf("not ok %i - Could not get x-vector, error code %d\n",
               test_id, rval);
    } else {
        char **colnames = malloc(ncols * sizeof(char *));
        mpq_QSget_colnames(p, colnames);
        printf("ok %i - Obtained variable values and names\n", test_id);
        free(colnames);
    }

CLEANUP:
    if (x) {
        for (i = 0; i < ncols; i++) mpq_clear(x[i]);
        free(x);
    }
    if (p) mpq_QSfree_prob(p);
}

static void test_solution_get_dual_values(int test_id)
{
    mpq_QSprob p = NULL;
    mpq_t *y = NULL;
    int status = 0;
    int rval = 0;
    int i;
    int nrows = 2;

    rval = solve_test_problem(&p, &status);
    if (rval) {
        printf("not ok %i - Unable to solve the LP\n", test_id);
        goto CLEANUP;
    }

    if (status != QS_LP_OPTIMAL) {
        printf("not ok %i - Did not find an optimal solution.\n", test_id);
        goto CLEANUP;
    }

    /* Test get dual values */
    y = malloc(nrows * sizeof (mpq_t));
    for (i = 0; i < nrows; i++) mpq_init(y[i]);
    rval = mpq_QSget_pi_array(p, y);
    if (rval) {
        printf("not ok %i - Could not get dual values, error code %d\n",
               test_id, rval);
    } else {
        char **rownames = malloc(nrows * sizeof(char *));
        mpq_QSget_rownames(p, rownames);
        printf("ok %i - Obtained dual values and contraint names\n", test_id);
        free(rownames);
    }

CLEANUP:
    if (y) {
        for (i = 0; i < nrows; i++) mpq_clear(y[i]);
        free(y);
    }
    if (p) mpq_QSfree_prob(p);
}

static void test_write_problem_to_file(int test_id)
{
    mpq_QSprob p = NULL;

    int rval = load_test_problem(&p);
    if (rval) {
        printf("not ok %i - Unable to load the LP\n", test_id);
        goto CLEANUP;
    }

    /* Test write to file */
    rval = mpq_QSwrite_prob (p, "/dev/null", "LP");
    if (rval) {
        printf("not ok %i - Could not write the LP, error code %d\n",
               test_id, rval);
    } else {
        printf("ok %i - LP written to /dev/null\n", test_id);
    }

CLEANUP:
    if (p) mpq_QSfree_prob(p);
}

int main(int ac, char **av)
{
    int rval = 0;

    QSexactStart();

    if (ac > 1) {
        printf ("Usage: %s\n", *av);
        rval = EXIT_FAILURE; goto CLEANUP;
    }

    /*  Must call QSexact_set_precision before any other QS function       */
    /*  Cannot use mpq_ILL_MAXDOUBLE or mpq_ILL_MINDOUBLE before this call */

    QSexact_set_precision(128);

    static test_func* test_functions[] = {
        test_new_column,
        test_new_column_without_name,
        test_delete_row,
        test_delete_invalid_row,
        test_delete_row_in_empty_problem,
        test_delete_column,
        test_row_count,
        test_col_count,
        test_solution_is_optimal,
        test_solution_objective,
        test_solution_get_variables,
        test_solution_get_dual_values,
        test_write_problem_to_file
    };

    /* Print number of tests */
    int count = sizeof(test_functions) / sizeof(test_functions[0]);
    printf("1..%i\n", count);

    /* Run tests */
    int i;
    for (i = 0; i < count; i++) {
        int test_id = i + 1;
        test_functions[i](test_id);
    }

CLEANUP:

    QSexactClear();
    return rval;
}
