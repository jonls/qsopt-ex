
from libc cimport stdlib, limits
cimport cqsoptex
cimport cgmp

import numbers, fractions


# Initialize library
# Python modules never unload so we don't have any place to
# run cqsoptex.QSexactClear()
cqsoptex.QSexactStart()
cqsoptex.QSexact_set_precision(128)


### Constants
# Simplex algorithm
class SimplexAlgorithm(object):
    PRIMAL = cqsoptex.PRIMAL_SIMPLEX
    DUAL = cqsoptex.DUAL_SIMPLEX
    PRIMAL_OR_DUAL = cqsoptex.PRIMAL_OR_DUAL

# Objective sense
class ObjectiveSense(object):
    MINIMIZE = cqsoptex.QS_MIN
    MAXIMIZE = cqsoptex.QS_MAX

# Solution status
class SolutionStatus(object):
    OPTIMAL = cqsoptex.QS_LP_OPTIMAL
    INFEASIBLE = cqsoptex.QS_LP_INFEASIBLE
    UNBOUNDED = cqsoptex.QS_LP_UNBOUNDED
    ITER_LIMIT = cqsoptex.QS_LP_ITER_LIMIT
    TIME_LIMIT = cqsoptex.QS_LP_TIME_LIMIT
    UNSOLVED = cqsoptex.QS_LP_UNSOLVED
    ABORTED = cqsoptex.QS_LP_ABORTED
    NUMERR = cqsoptex.QS_LP_NUMERR
    OBJ_LIMIT = cqsoptex.QS_LP_OBJ_LIMIT
    MODIFIED = cqsoptex.QS_LP_MODIFIED
    CHANGE_PREC = cqsoptex.QS_LP_CHANGE_PREC

# Parameters that can be set by QSset_param
class Parameter(object):
    PRIMAL_PRICING = cqsoptex.QS_PARAM_PRIMAL_PRICING
    DUAL_PRICING = cqsoptex.QS_PARAM_DUAL_PRICING
    SIMPLEX_DISPLAY = cqsoptex.QS_PARAM_SIMPLEX_DISPLAY
    SIMPLEX_MAX_ITERATIONS = cqsoptex.QS_PARAM_SIMPLEX_MAX_ITERATIONS
    SIMPLEX_MAX_TIME = cqsoptex.QS_PARAM_SIMPLEX_MAX_TIME
    SIMPLEX_SCALING = cqsoptex.QS_PARAM_SIMPLEX_SCALING
    OBJULIM = cqsoptex.QS_PARAM_OBJULIM
    OBJLLIM = cqsoptex.QS_PARAM_OBJLLIM

# Values for pricing parameter
class Pricing(object):
    # Primal pricing
    PDANTZIG = cqsoptex.QS_PRICE_PDANTZIG
    PDEVEX = cqsoptex.QS_PRICE_PDEVEX
    PSTEEP = cqsoptex.QS_PRICE_PSTEEP
    PMULTPARTIAL = cqsoptex.QS_PRICE_PMULTPARTIAL

    # Dual pricing
    DDANTZIG = cqsoptex.QS_PRICE_DDANTZIG
    DDEVEX = cqsoptex.QS_PRICE_DDEVEX
    DSTEEP = cqsoptex.QS_PRICE_DSTEEP
    DMULTPARTIAL = cqsoptex.QS_PRICE_DMULTPARTIAL

# Constraint sense
class ConstraintSense(object):
    EQUAL = 'E'
    GREATER = 'G'
    LESS = 'L'


# This function is not yet in the cpython.long pxd file
cdef extern from 'Python.h':
    long PyLong_AsLongAndOverflow(object value, int* overflow)


# Conversion from numeric Python types to GMP types
cdef int mpz_set_pylong(cgmp.mpz_t rop, object value) except -1:
    '''Set mpz_t value from Python int/long object'''

    cdef int overflow
    cdef long long_value
    cdef cgmp.mpz_t z_temp

    if not isinstance(value, (int, long)):
        raise ValueError('Value must be a Python int or long')

    long_value = PyLong_AsLongAndOverflow(value, &overflow)
    if not overflow:
        cgmp.mpz_set_si(rop, long_value)
    else:
        cgmp.mpz_init(z_temp)
        cgmp.mpz_set_ui(rop, 0)
        try:
            sign = -1 if value < 0 else 1
            value = abs(value)
            offset = 0
            while value > 0:
                cgmp.mpz_set_ui(z_temp, value & 0xffff)
                cgmp.mpz_mul_2exp(z_temp, z_temp, offset)
                cgmp.mpz_add(rop, rop, z_temp)
                value >>= 16
                offset += 16

            cgmp.mpz_mul_si(rop, rop, sign)
        finally:
            cgmp.mpz_clear(z_temp)

    return 0

cdef int mpq_set_pyrational(cgmp.mpq_t rop, object value) except -1:
    '''Set mpq_t value from Python Rational object'''

    if not isinstance(value, numbers.Rational):
        raise ValueError('Value must be a Python numbers.Rational')

    mpz_set_pylong(cgmp.mpq_numref(rop), value.numerator)
    mpz_set_pylong(cgmp.mpq_denref(rop), value.denominator)
    cgmp.mpq_canonicalize(rop)

    return 0

cdef int mpq_set_pynumeric(cgmp.mpq_t rop, object value) except -1:
    '''Set mpq_t value from Python numeric object'''
    if not isinstance(value, numbers.Rational):
        # Try to convert to fractional value
        value = fractions.Fraction(value)
    mpq_set_pyrational(rop, value)

    return 0


# Conversion from GMP types to Python numeric types
cdef object pylong_from_mpz(const cgmp.mpz_t value):
    '''Return Python long representing the value'''
    cdef size_t i

    result = 0
    for i in range(cgmp.mpz_size(value), 0, -1):
        result <<= limits.CHAR_BIT*sizeof(cgmp.mp_limb_t)
        result += cgmp.mpz_getlimbn(value, i-1)

    return cgmp.mpz_sgn(value) * result

cdef object pyrational_from_mpq(const cgmp.mpq_t value):
    '''Return Python Rational (long or Fraction) representing the value'''

    num = pylong_from_mpz(cgmp.mpq_numref(value))
    denom = pylong_from_mpz(cgmp.mpq_denref(value))
    if denom == 1:
        return num
    return fractions.Fraction(num, denom)


def print_version():
    cqsoptex.QSopt_ex_version()


# Wrapper around exact LP problem
class ExactProblemError(Exception):
    pass

cdef class ExactProblem:
    '''Representation of exact LP problem'''

    # C API Problem
    cdef cqsoptex.mpq_QSdata* _c_qsdata

    # Solution values as mpq_t
    cdef int _c_sol_nvars
    cdef cqsoptex.mpq_t* _c_sol_x

    def __cinit__(self):
        qsdata = cqsoptex.mpq_QScreate_prob('problem', ObjectiveSense.MAXIMIZE)
        if qsdata is NULL:
            raise MemoryError()
        self._c_qsdata = qsdata

        self._c_sol_nvars = 0
        self._c_sol_x = NULL

    def __dealloc__(self):
        cqsoptex.mpq_QSfree_prob(self._c_qsdata)
        stdlib.free(self._c_sol_x)

    def _invalidate_solution(self):
        '''Free the cached solution values'''
        if self._c_sol_x is not NULL:
            for i in range(self._c_sol_nvars):
                cgmp.mpq_clear(self._c_sol_x[i])
            stdlib.free(self._c_sol_x)
        self._c_sol_nvars = 0
        self._c_sol_x = NULL

    def _index_maybe_string(self, variable):
        '''Get the column index from a variable that is either name (string) or index'''

        cdef int r, colindex

        if isinstance(variable, basestring):
            r = cqsoptex.mpq_QSget_column_index(self._c_qsdata, variable, &colindex)
            if r != 0:
                raise ExactProblemError('An error occured in QSget_column_index()')
            return colindex
        else:
            return int(variable)

    def get_variable_count(self):
        '''Get number of variables defined in the problem'''
        return cqsoptex.mpq_QSget_colcount(self._c_qsdata)

    def get_constraint_count(self):
        '''Get number of constraints defined in the problem'''
        return cqsoptex.mpq_QSget_rowcount(self._c_qsdata)

    def read(self, path, filetype='MPS'):
        '''Read LP problem from the given file'''
        self._invalidate_solution()
        qsdata = cqsoptex.mpq_QSread_prob(path, filetype)
        if qsdata is NULL:
            raise ExactProblemError('An error occured in QSread_prob()')
        cqsoptex.mpq_QSfree_prob(self._c_qsdata)
        self._c_qsdata = qsdata

    def write(self, path, filetype='MPS'):
        '''Write LP problem to the given file'''
        cdef int r = cqsoptex.mpq_QSwrite_prob(self._c_qsdata, path, filetype)
        if r != 0:
            raise ExactProblemError('An error occured in QSwrite_prob()')

    def add_variable(self, objective=None, lower=None, upper=None, name=None):
        '''Add variable to problem'''
        cdef cgmp.mpq_t objective_q, lower_q, upper_q
        cdef const char* n
        cdef int r

        # Remove cache of solution
        self._invalidate_solution()

        cgmp.mpq_init(objective_q)
        cgmp.mpq_init(lower_q)
        cgmp.mpq_init(upper_q)

        try:
            # Set objective
            if objective is not None:
                mpq_set_pynumeric(objective_q, objective)

            # Set lower bound
            if lower is not None:
                mpq_set_pynumeric(lower_q, lower)
            else:
                cgmp.mpq_set(lower_q, cqsoptex.mpq_NINFTY)

            # Set upper bound
            if upper is not None:
                mpq_set_pynumeric(upper_q, upper)
            else:
                cgmp.mpq_set(upper_q, cqsoptex.mpq_INFTY)

            # Set name
            if name is not None:
                n = name
            else:
                n = NULL

            # Create column
            r = cqsoptex.mpq_QSnew_col(self._c_qsdata, objective_q, lower_q, upper_q, n)
            if r != 0:
                raise ExactProblemError('An error occured in QSnew_col()')
        finally:
            cgmp.mpq_clear(objective_q)
            cgmp.mpq_clear(lower_q)
            cgmp.mpq_clear(upper_q)

    def add_linear_constraint(self, sense, values=None, rhs=0, name=None):
        '''Add linear constraint to problem'''
        cdef cgmp.mpq_t rhs_q
        cdef const char* n
        cdef char sense_c
        cdef int r, i

        cdef int count = 0
        cdef int* indices = NULL
        cdef int count_q = 0
        cdef cgmp.mpq_t* values_q = NULL

        # Remove cache of solution
        self._invalidate_solution()

        # There seems to be no way of obtaining the index of the row that
        # was just added. Since we do not know the name either (it could be None
        # or sometimes libqsopt_ex will even change it), we cannot use
        # QSget_row_index() either.
        # Essentially, we cannot refer to the column the was just added so
        # all values must be added in the same call. (Otherwise we could use
        # QSnew_row() and followed by a number of calls to QSchange_coef()).
        # This makes the following code quite complex as we have to build
        # the array of indices and values before the call to QSadd_row() and
        # we have to properly clean up the allocated memory afterwards.

        # Count the number of values in constraint
        if values is None:
            values = ()
        if isinstance(values, dict):
            values = values.iteritems()
        values = list(values)
        count = len(values)

        cgmp.mpq_init(rhs_q)

        try:
            # Allocate memory for values
            indices = <int*>stdlib.malloc(count * sizeof(int))
            if indices is NULL:
                raise MemoryError()

            values_q = <cgmp.mpq_t*>stdlib.malloc(count * sizeof(cgmp.mpq_t))
            if values_q is NULL:
                raise MemoryError()
            count_q = count

            # Initialize mpq_t values
            for i in range(count_q):
                cgmp.mpq_init(values_q[i])

            # Copy variable values from Python numbers
            for i, pair in enumerate(values):
                variable, value = pair

                # Get variable index and value
                indices[i] = self._index_maybe_string(variable)
                mpq_set_pynumeric(values_q[i], value)

            # Set right-hand side
            mpq_set_pynumeric(rhs_q, rhs)

            # First character in Python string is sense symbol
            sense_c = ord(sense[0])

            # Set name
            if name is not None:
                n = name
            else:
                n = NULL

            r = cqsoptex.mpq_QSadd_row(self._c_qsdata, count, indices, values_q, &rhs_q, sense_c, n)
            if r != 0:
                raise ExactProblemError('An error occured in QSnew_row()')
        finally:
            cgmp.mpq_clear(rhs_q)
            stdlib.free(indices)

            for i in range(count_q):
                cgmp.mpq_clear(values_q[i])
            stdlib.free(values_q)

    def set_objective_sense(self, sense):
        '''Set objective sense (i.e. minimize or maximize)'''

        cdef int r = cqsoptex.mpq_QSchange_objsense(self._c_qsdata, sense)
        if r != 0:
            raise ExactProblemError('An error occured in QSchange_objsense()')

    def set_linear_objective(self, values):
        '''Set linear objective values from dict or iterator of pairs'''

        cdef cgmp.mpq_t value_q
        cdef int r

        if isinstance(values, dict):
            values = values.iteritems()

        cgmp.mpq_init(value_q)
        try:
            for variable, value in values:
                index = self._index_maybe_string(variable)
                mpq_set_pynumeric(value_q, value)

                r = cqsoptex.mpq_QSchange_objcoef(self._c_qsdata, index, value_q)
                if r != 0:
                    raise ExactProblemError('An error occured in QSchange_objcoef()')
        finally:
            cgmp.mpq_clear(value_q)

    def solve(self):
        '''Solve problem and return status'''
        cdef int status
        cdef int r, i

        # Free any previous solution values
        self._invalidate_solution()

        r = cqsoptex.QSexact_solver(self._c_qsdata, NULL, NULL, NULL,
                                    cqsoptex.PRIMAL_SIMPLEX, &status)
        if r != 0:
            raise ExactProblemError('An error occured in QSexact_solver()')

        if status == SolutionStatus.OPTIMAL:
            # Allocate space for solution values
            self._c_sol_nvars = cqsoptex.mpq_QSget_colcount(self._c_qsdata)
            self._c_sol_x = <cgmp.mpq_t*>stdlib.malloc(self._c_sol_nvars * sizeof(cgmp.mpq_t))
            if self._c_sol_x is NULL:
                raise MemoryError()

            for i in range(self._c_sol_nvars):
                cgmp.mpq_init(self._c_sol_x[i])

            # Cache solution values
            r = cqsoptex.mpq_QSget_x_array(self._c_qsdata, self._c_sol_x)
            if r != 0:
                raise ExactProblemError('An error occured in QSget_x_array()')

        return status

    def get_status(self):
        '''Get status of problem'''
        cdef int status
        cdef int r = cqsoptex.mpq_QSget_status(self._c_qsdata, &status)
        if r != 0:
            raise ExactProblemError('An error occured in QSget_status()')
        return status

    def get_param(self, param):
        '''Get value of parameter'''
        cdef int value
        cdef int r = cqsoptex.mpq_QSget_param(self._c_qsdata, param, &value)
        if r != 0:
            raise ExactProblemError('An error occured in QSget_param()')
        return value

    def set_param(self, param, value):
        '''Set parameter to value'''
        cdef int r = cqsoptex.mpq_QSset_param(self._c_qsdata, param, value)
        if r != 0:
            raise ExactProblemError('An error occured in QSset_param()')

    def get_value(self, variable):
        '''Get value of variable'''

        cdef int index = self._index_maybe_string(variable)
        if index < 0 or index >= self._c_sol_nvars:
            raise IndexError('Invalid variable index')
        return pyrational_from_mpq(self._c_sol_x[index])

    def get_objective_value(self):
        '''Get value of objective'''
        cdef int r
        cdef cgmp.mpq_t value
        cgmp.mpq_init(value)

        result = None
        try:
            r = cqsoptex.mpq_QSget_objval(self._c_qsdata, &value)
            if r != 0:
                raise ExactProblemError('An error occured in QSget_objval()')
            result = pyrational_from_mpq(value)
        finally:
            cgmp.mpq_clear(value)

        return result
