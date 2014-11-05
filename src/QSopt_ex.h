/* QSopt_ex.h -- QSopt_ex main header file
 * This file is part of QSopt_ex.
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
 *
 * Copyright (c) 2014  Jon Lund Steffensen <jonlst@gmail.com> 
 */

#ifndef QSOPT_EX_H
#define QSOPT_EX_H

#include "eg_elist.h"
#include "eg_fp.h"
#include "eg_io.h"
#include "eg_lpnum.h"
#include "eg_macros.h"
#include "eg_mem.h"
#include "eg_memslab.h"

#include "eg_lpnum.dbl.h"
#include "eg_lpnum.float128.h"
#include "eg_lpnum.fp20.h"
#include "eg_lpnum.int32.h"
#include "eg_lpnum.int.h"
#include "eg_lpnum.ldbl.h"
#include "eg_lpnum.llint.h"
#include "eg_lpnum.mpf.h"
#include "eg_lpnum.mpq.h"
#include "eg_lpnum.mpz.h"

#include "basicdefs.h"
#include "urandom.h"
#include "symtab.h"
#include "reporter.h"
#include "allocrus.h"
#include "bgetopt.h"
#include "zeit.h"
#include "except.h"
#include "qs_config.h"

/* Double template headers */
#include "qstruct_dbl.h"
#include "editor_dbl.h"
#include "dstruct_dbl.h"
#include "factor_dbl.h"
#include "lpdefs_dbl.h"
#include "readline_dbl.h"
#include "lpdata_dbl.h"
#include "basis_dbl.h"
#include "dheaps_i_dbl.h"
#include "qsopt_dbl.h"
#include "format_dbl.h"
#include "rawlp_dbl.h"
#include "mps_dbl.h"
#include "price_dbl.h"
#include "priority_dbl.h"
#include "ratio_dbl.h"
#include "read_lp_dbl.h"
#include "read_mps_dbl.h"
#include "simplex_dbl.h"
#include "write_lp_dbl.h"
#include "lib_dbl.h"
#include "eg_numutil_dbl.h"
#include "eg_eheap_dbl.h"

/* FP20 template headers */
#include "qstruct_fp20.h"
#include "editor_fp20.h"
#include "dstruct_fp20.h"
#include "factor_fp20.h"
#include "lpdefs_fp20.h"
#include "readline_fp20.h"
#include "lpdata_fp20.h"
#include "basis_fp20.h"
#include "dheaps_i_fp20.h"
#include "qsopt_fp20.h"
#include "format_fp20.h"
#include "rawlp_fp20.h"
#include "mps_fp20.h"
#include "price_fp20.h"
#include "priority_fp20.h"
#include "ratio_fp20.h"
#include "read_lp_fp20.h"
#include "read_mps_fp20.h"
#include "simplex_fp20.h"
#include "write_lp_fp20.h"
#include "lib_fp20.h"
#include "eg_numutil_fp20.h"
#include "eg_eheap_fp20.h"

/* MPQ template headers */
#include "qstruct_mpq.h"
#include "editor_mpq.h"
#include "dstruct_mpq.h"
#include "factor_mpq.h"
#include "lpdefs_mpq.h"
#include "readline_mpq.h"
#include "lpdata_mpq.h"
#include "basis_mpq.h"
#include "dheaps_i_mpq.h"
#include "qsopt_mpq.h"
#include "format_mpq.h"
#include "rawlp_mpq.h"
#include "mps_mpq.h"
#include "price_mpq.h"
#include "priority_mpq.h"
#include "ratio_mpq.h"
#include "read_lp_mpq.h"
#include "read_mps_mpq.h"
#include "simplex_mpq.h"
#include "write_lp_mpq.h"
#include "lib_mpq.h"
#include "eg_numutil_mpq.h"
#include "eg_eheap_mpq.h"

/* MPF template headers */
#include "qstruct_mpf.h"
#include "editor_mpf.h"
#include "dstruct_mpf.h"
#include "factor_mpf.h"
#include "lpdefs_mpf.h"
#include "readline_mpf.h"
#include "lpdata_mpf.h"
#include "basis_mpf.h"
#include "dheaps_i_mpf.h"
#include "qsopt_mpf.h"
#include "format_mpf.h"
#include "rawlp_mpf.h"
#include "mps_mpf.h"
#include "price_mpf.h"
#include "priority_mpf.h"
#include "ratio_mpf.h"
#include "read_lp_mpf.h"
#include "read_mps_mpf.h"
#include "simplex_mpf.h"
#include "write_lp_mpf.h"
#include "lib_mpf.h"
#include "eg_numutil_mpf.h"
#include "eg_eheap_mpf.h"

/* Float128 template headers */
#ifdef HAVE_SOFTFLOAT
#include "qstruct_float128.h"
#include "editor_float128.h"
#include "dstruct_float128.h"
#include "factor_float128.h"
#include "lpdefs_float128.h"
#include "readline_float128.h"
#include "lpdata_float128.h"
#include "basis_float128.h"
#include "dheaps_i_float128.h"
#include "qsopt_float128.h"
#include "format_float128.h"
#include "rawlp_float128.h"
#include "mps_float128.h"
#include "price_float128.h"
#include "priority_float128.h"
#include "ratio_float128.h"
#include "read_lp_float128.h"
#include "read_mps_float128.h"
#include "simplex_float128.h"
#include "write_lp_float128.h"
#include "lib_float128.h"
#include "eg_numutil_float128.h"
#include "eg_eheap_float128.h"
#endif

/* Long double template headers */
#if ENABLE_LONG_DOUBLE
#include "qstruct_ldbl.h"
#include "editor_ldbl.h"
#include "dstruct_ldbl.h"
#include "factor_ldbl.h"
#include "lpdefs_ldbl.h"
#include "readline_ldbl.h"
#include "lpdata_ldbl.h"
#include "basis_ldbl.h"
#include "dheaps_i_ldbl.h"
#include "qsopt_ldbl.h"
#include "format_ldbl.h"
#include "rawlp_ldbl.h"
#include "mps_ldbl.h"
#include "price_ldbl.h"
#include "priority_ldbl.h"
#include "ratio_ldbl.h"
#include "read_lp_ldbl.h"
#include "read_mps_ldbl.h"
#include "simplex_ldbl.h"
#include "write_lp_ldbl.h"
#include "lib_ldbl.h"
#include "eg_numutil_ldbl.h"
#include "eg_eheap_ldbl.h"
#endif

#include "exact.h"
#include "eg_exutil.h"

#endif /* !QSOPT_EX_H */
