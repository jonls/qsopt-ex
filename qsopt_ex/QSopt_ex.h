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
 * Copyright (c) 2014-2015  Jon Lund Steffensen <jonlst@gmail.com>
 */

#ifndef QSOPT_EX_H
#define QSOPT_EX_H

#include "eg_elist.h"
#include "eg_io.h"
#include "eg_lpnum.h"
#include "eg_macros.h"
#include "eg_mem.h"
#include "eg_memslab.h"

#include "eg_lpnum.dbl.h"
#include "eg_lpnum.mpf.h"
#include "eg_lpnum.mpq.h"

#include "basicdefs.h"
#include "urandom.h"
#include "symtab.h"
#include "reporter.h"
#include "bgetopt.h"
#include "zeit.h"
#include "logging.h"

#include "QSopt_ex_version.h"

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
#include "price_dbl.h"
#include "priority_dbl.h"
#include "ratio_dbl.h"
#include "read_lp_dbl.h"
#include "read_mps_dbl.h"
#include "simplex_dbl.h"
#include "write_lp_dbl.h"
#include "lib_dbl.h"
#include "eg_numutil_dbl.h"

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
#include "price_mpq.h"
#include "priority_mpq.h"
#include "ratio_mpq.h"
#include "read_lp_mpq.h"
#include "read_mps_mpq.h"
#include "simplex_mpq.h"
#include "write_lp_mpq.h"
#include "lib_mpq.h"
#include "eg_numutil_mpq.h"

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
#include "price_mpf.h"
#include "priority_mpf.h"
#include "ratio_mpf.h"
#include "read_lp_mpf.h"
#include "read_mps_mpf.h"
#include "simplex_mpf.h"
#include "write_lp_mpf.h"
#include "lib_mpf.h"
#include "eg_numutil_mpf.h"

#include "exact.h"
#include "eg_exutil.h"

#endif /* !QSOPT_EX_H */
