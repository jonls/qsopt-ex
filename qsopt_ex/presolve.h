/****************************************************************************/
/*                                                                          */
/*  This file is part of QSopt_ex.                                          */
/*                                                                          */
/*  (c) Copyright 2006 by David Applegate, William Cook, Sanjeeb Dash,      */
/*  and Daniel Espinoza                                                     */
/*                                                                          */
/*  Sanjeeb Dash ownership of copyright in QSopt_ex is derived from his     */
/*  copyright in QSopt.                                                     */
/*                                                                          */
/*  This code may be used under the terms of the GNU General Public License */
/*  (Version 2.1 or later) as published by the Free Software Foundation.    */
/*                                                                          */
/*  Alternatively, use is granted for research purposes only.               */
/*                                                                          */
/*  It is your choice of which of these two licenses you are operating      */
/*  under.                                                                  */
/*                                                                          */
/*  We make no guarantees about the correctness or usefulness of this code. */
/*                                                                          */
/****************************************************************************/

#ifndef EGLPNUM_TYPENAME_PRESOLVE_H
#define EGLPNUM_TYPENAME_PRESOLVE_H

#include "lpdata_EGLPNUM_TYPENAME.h"

void EGLPNUM_TYPENAME_ILLlp_sinfo_init(EGLPNUM_TYPENAME_ILLlp_sinfo * sinfo);
void EGLPNUM_TYPENAME_ILLlp_sinfo_free(EGLPNUM_TYPENAME_ILLlp_sinfo * sinfo);
int EGLPNUM_TYPENAME_ILLlp_sinfo_print(EGLPNUM_TYPENAME_ILLlp_sinfo * s);

void EGLPNUM_TYPENAME_ILLlp_predata_init(EGLPNUM_TYPENAME_ILLlp_predata * pre);
void EGLPNUM_TYPENAME_ILLlp_predata_free(EGLPNUM_TYPENAME_ILLlp_predata * pre);

void EGLPNUM_TYPENAME_ILLlp_preop_init(EGLPNUM_TYPENAME_ILLlp_preop * op);
void EGLPNUM_TYPENAME_ILLlp_preop_free(EGLPNUM_TYPENAME_ILLlp_preop * op);

void EGLPNUM_TYPENAME_ILLlp_preline_init(EGLPNUM_TYPENAME_ILLlp_preline * line);
void EGLPNUM_TYPENAME_ILLlp_preline_free(EGLPNUM_TYPENAME_ILLlp_preline * line);

#endif /* ! EGLPNUM_TYPENAME_PRESOLVE_H */
