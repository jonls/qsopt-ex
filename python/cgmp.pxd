
from libc.stdio cimport FILE

cdef extern from 'gmp.h' nogil:
    # Special types
    ctypedef int mp_size_t
    ctypedef int mp_limb_t
    ctypedef int mp_bitcnt_t

    # Integer type
    ctypedef struct mpz_t:
        pass

    void mpz_init(mpz_t x)
    void mpz_clear(mpz_t x)

    void mpz_set(mpz_t rop, const mpz_t op)
    void mpz_set_ui(mpz_t rop, unsigned long int op)
    void mpz_set_si(mpz_t rop, signed long int op)

    void mpz_add(mpz_t rop, const mpz_t op1, const mpz_t op2)
    void mpz_add_ui(mpz_t rop, const mpz_t op1, unsigned long int op2)
    void mpz_mul_si(mpz_t rop, const mpz_t op1, signed long int op2)
    void mpz_mul_2exp(mpz_t rop, const mpz_t op1, mp_bitcnt_t op2)
    void mpz_neg(mpz_t rop, const mpz_t op)
    void mpz_abs(mpz_t rop, const mpz_t op)

    int mpz_sgn(const mpz_t op)

    size_t mpz_size(const mpz_t op)
    mp_limb_t mpz_getlimbn(const mpz_t x, mp_size_t n)

    size_t mpz_out_str(FILE* stream, int base, const mpz_t op)
    size_t mpz_inp_str(mpz_t rop, FILE* stream, int base)

    # Rational type
    ctypedef struct mpq_t:
        pass

    void mpq_canonicalize(mpq_t op)

    void mpq_init(mpq_t x)
    void mpq_clear(mpq_t x)

    void mpq_set(mpq_t rop, const mpq_t op)
    void mpq_set_ui(mpq_t rop, unsigned long int op1, unsigned long int op2)
    void mpq_set_si(mpq_t rop, signed long int op1, unsigned long int op2)

    void mpq_neg(mpq_t negated_operand, const mpq_t operand)

    mpz_t mpq_numref(const mpq_t op)
    mpz_t mpq_denref(const mpq_t op)

    size_t mpq_out_str(FILE* stream, int base, const mpq_t op)
    size_t mpq_inp_str(mpq_t rop, FILE* stream, int base)

    # Floating-point type
    ctypedef struct mpf_t:
        pass
