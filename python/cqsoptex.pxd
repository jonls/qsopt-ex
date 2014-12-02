
from cgmp cimport mpq_t

cdef extern from 'qsopt_ex/basicdefs.h' nogil:
    # Simplex algorithm
    enum: PRIMAL_SIMPLEX
    enum: DUAL_SIMPLEX
    enum: PRIMAL_OR_DUAL

    # Objective sense
    enum: QS_MIN
    enum: QS_MAX

    # Solution status
    enum: QS_LP_OPTIMAL
    enum: QS_LP_INFEASIBLE
    enum: QS_LP_UNBOUNDED
    enum: QS_LP_ITER_LIMIT
    enum: QS_LP_TIME_LIMIT
    enum: QS_LP_UNSOLVED
    enum: QS_LP_ABORTED
    enum: QS_LP_NUMERR
    enum: QS_LP_OBJ_LIMIT
    enum: QS_LP_MODIFIED
    enum: QS_LP_CHANGE_PREC

    # Parameters that can be set by QSset_param
    enum: QS_PARAM_PRIMAL_PRICING
    enum: QS_PARAM_DUAL_PRICING
    enum: QS_PARAM_SIMPLEX_DISPLAY
    enum: QS_PARAM_SIMPLEX_MAX_ITERATIONS
    enum: QS_PARAM_SIMPLEX_MAX_TIME
    enum: QS_PARAM_SIMPLEX_SCALING
    enum: QS_PARAM_OBJULIM
    enum: QS_PARAM_OBJLLIM

    # Values for primal pricing parameter
    enum: QS_PRICE_PDANTZIG
    enum: QS_PRICE_PDEVEX
    enum: QS_PRICE_PSTEEP
    enum: QS_PRICE_PMULTPARTIAL

    # Values for dual pricing parameter
    enum: QS_PRICE_DDANTZIG
    enum: QS_PRICE_DSTEEP
    enum: QS_PRICE_DMULTPARTIAL
    enum: QS_PRICE_DDEVEX

    ctypedef struct QSbasis:
        pass

cdef extern from 'qsopt_ex/lpdefs_mpq.h' nogil:
    const mpq_t mpq_INFTY
    const mpq_t mpq_NINFTY

cdef extern from 'qsopt_ex/qstruct_mpq.h' nogil:
    ctypedef struct mpq_QSdata:
        pass

cdef extern from 'qsopt_ex/QSopt_ex_version.h' nogil:
    void QSopt_ex_version()

cdef extern from 'qsopt_ex/qsopt_mpq.h' nogil:
    mpq_QSdata* mpq_QScreate_prob(const char* name, int objsense)
    void mpq_QSfree_prob(mpq_QSdata* problem)

    int mpq_QSget_param(mpq_QSdata* problem, int param, int* value)
    int mpq_QSset_param(mpq_QSdata* problem, int param, int value)

    int mpq_QSget_objval(mpq_QSdata* problem, mpq_t* value)
    int mpq_QSget_x_array(mpq_QSdata* problem, mpq_t* value)
    int mpq_QSget_solution(mpq_QSdata* problem, mpq_t* value, mpq_t* x, mpq_t* pi, mpq_t* slack, mpq_t* rc)
    int mpq_QSget_status(mpq_QSdata* problem, int* status)

    int mpq_QSget_colcount(mpq_QSdata* problem)
    int mpq_QSget_rowcount(mpq_QSdata* problem)

    int mpq_QSget_column_index(mpq_QSdata* problem, const char* name, int* index)

    int mpq_QSchange_objsense(mpq_QSdata* problem, int sense)
    int mpq_QSchange_objcoef(mpq_QSdata* problem, int index, mpq_t value)

    int mpq_QSnew_col(mpq_QSdata* problem, const mpq_t objective, const mpq_t lower, const mpq_t upper, const char* name)
    int mpq_QSadd_row(mpq_QSdata* problem, int count, int* indices, const mpq_t* values, const mpq_t* rhs, char sense, const char* name)

    mpq_QSdata* mpq_QSread_prob(const char* filepath, const char* filetype)
    int mpq_QSwrite_prob(mpq_QSdata* problem, const char* filepath, const char* filetype)

cdef extern from 'qsopt_ex/exact.h' nogil:
    void QSexactStart()
    void QSexactClear()

    int QSexact_solver(mpq_QSdata* problem, mpq_t* const x, mpq_t* const y, QSbasis* const basis, int simplexalgo, int* status)
    void QSexact_set_precision(unsigned int precision)
