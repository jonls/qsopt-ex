
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <libgen.h>

#include <gmp.h>

#include "qs_config.h"
#include "QSopt_ex.h"

#define max(x,y) ((x) > (y) ? (x) : (y))

mpq_QSdata *try_read(const char *base, const char *name)
{
    mpq_QSdata *result = NULL;
    size_t size = 0;
    size_t new_size;
    char *path = NULL;

    /* Try to open file using the given name as MPS */
    new_size = snprintf(NULL, 0, "%s/%s", base, name);
    if (new_size > size) {
        size = new_size;
        path = realloc(path, size + 1);
        if (path == NULL) abort();
    }
    sprintf(path, "%s/%s", base, name);

    result = mpq_QSread_prob(path, "MPS");
    if (result != NULL) {
        free(path);
        return result;
    }

    /* Try to open file with ".mps" extension */
    new_size = snprintf(NULL, 0, "%s/%s.mps", base, name);
    if (new_size > size) {
        size = new_size;
        path = realloc(path, size + 1);
        if (path == NULL) abort();
    }
    sprintf(path, "%s/%s.mps", base, name);

    result = mpq_QSread_prob(path, "MPS");
    if (result != NULL) {
        free(path);
        return result;
    }

    /* Try to open file with ".lp" extension */
    new_size = snprintf(NULL, 0, "%s/%s.lp", base, name);
    if (new_size > size) {
        size = new_size;
        path = realloc(path, size + 1);
        if (path == NULL) abort();
    }
    sprintf(path, "%s/%s.lp", base, name);

    result = mpq_QSread_prob(path, "LP");
    if (result != NULL) {
        free(path);
        return result;
    }

    free(path);
    return NULL;
}

int main(int argc, char *argv[])
{
    /* Initialize */
    QSopt_ex_version();
    QSexactStart();
    QSexact_set_precision(128);

    if (argc < 2) {
        fprintf(stderr, "Usage: %s FILE\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    /* Open file containing problem descriptions */
    char buffer[512];
    FILE *f = fopen(argv[1], "r");
    if (f == NULL) {
        perror("open");
        exit(EXIT_FAILURE);
    }

    /* Problem files will be loaded relative to this path */
    char *base_path = dirname(argv[1]);
    if (base_path == NULL) {
        perror("dirname");
        exit(EXIT_FAILURE);
    }

    base_path = strdup(base_path);
    if (base_path == NULL) {
        perror("strdup");
        exit(EXIT_FAILURE);
    }

    int problem_count = 0;
    int problem_array_size = 0;
    char **problem_names = NULL;
    mpq_t *problem_ref = NULL;

    int line = 0;
    while (fgets(buffer, 512, f) != NULL) {
        line += 1;

        if (buffer[0] == '\0' || buffer[0] == '#') continue;

        char *value_str = strchr(buffer, '\t');
        if (value_str == NULL) {
            fprintf(stderr, "Error parsing line %d: No value found!\n", line);
            exit(EXIT_FAILURE);
        }

        /* Separate name and value */
        *value_str = '\0';
        value_str += 1;

        if (problem_count >= problem_array_size) {
            problem_array_size = 2 * max(problem_array_size, 1);
            problem_names = realloc(problem_names, problem_array_size*sizeof(char*));
            mpq_EGlpNumReallocArray(&problem_ref, problem_array_size);
        }

        problem_names[problem_count] = strdup(buffer);
        if (problem_names[problem_count] == NULL) {
            perror("strdup");
            exit(EXIT_FAILURE);
        }
        mpq_EGlpNumReadStr(problem_ref[problem_count], value_str);
        problem_count += 1;
    }

    if (ferror(f)) {
        fprintf(stderr, "Error reading file!\n");
        exit(EXIT_FAILURE);
    }

    fclose(f);

    if (problem_count < 1) {
        fprintf(stderr, "No problems found!\n");
    }

    printf("1..%d\n", problem_count);

    for (int i = 0; i < problem_count; i++) {
        int r;
        int status;
        int prob_id = i+1;

        mpq_QSdata *prob = try_read(base_path, problem_names[i]);
        if (prob == NULL) {
            printf("not ok %d - Failed to read problem %s\n", prob_id, problem_names[i]);
            continue;
        }

        r = mpq_QSset_param(prob, QS_PARAM_SIMPLEX_DISPLAY, 1);
        if (r != 0) {
            printf("not ok %d - Failed to set parameter: SIMPLEX_DISPLAY\n", prob_id);
            mpq_QSfree_prob(prob);
            continue;
        }

        r = mpq_QSset_param(prob, QS_PARAM_PRIMAL_PRICING, QS_PRICE_PSTEEP);
        if (r != 0) {
            printf("not ok %d - Failed to set parameter: PRIMAL_PRICING\n", prob_id);
            mpq_QSfree_prob(prob);
            continue;
        }

        r = mpq_QSset_param(prob, QS_PARAM_DUAL_PRICING, QS_PRICE_DSTEEP);
        if (r != 0) {
            printf("not ok %d - Failed to set parameter: DUAL_PRICING\n", prob_id);
            mpq_QSfree_prob(prob);
            continue;
        }

        r = mpq_QSset_param(prob, QS_PARAM_SIMPLEX_SCALING, 1);
        if (r != 0) {
            printf("not ok %d - Failed to set parameter: SIMPLEX_SCALING\n", prob_id);
            mpq_QSfree_prob(prob);
            continue;
        }

        r = QSexact_solver (prob, NULL, NULL, NULL, PRIMAL_SIMPLEX, &status);
        if (r != 0) {
            printf("not ok %d - Failed to solve problem\n", prob_id);
            mpq_QSfree_prob(prob);
            continue;
        }

        if (status != QS_LP_OPTIMAL) {
            printf("not ok %d - Unable to solve problem\n", prob_id);
            mpq_QSfree_prob(prob);
            continue;
        }

        mpq_t obj;
        mpq_init(obj);
        mpq_QSget_objval(prob, &obj);

        fprintf(stderr, "Expected: %g, got: %g\n",
                mpq_EGlpNumToLf(problem_ref[i]), mpq_EGlpNumToLf(obj));
        mpq_clear(obj);

        printf("ok %d\n", prob_id);

        mpq_QSfree_prob(prob);
    }

    /* Clean up */
    QSexactClear();

    return 0;
}
