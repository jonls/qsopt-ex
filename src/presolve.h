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

#ifndef PRESOLVE_H
#define PRESOLVE_H

#include "lpdata.h"

void ILLlp_sinfo_init(ILLlp_sinfo * sinfo);
void ILLlp_sinfo_free(ILLlp_sinfo * sinfo);
int ILLlp_sinfo_print(ILLlp_sinfo * s);

void ILLlp_predata_init(ILLlp_predata * pre);
void ILLlp_predata_free(ILLlp_predata * pre);

void ILLlp_preop_init(ILLlp_preop * op);
void ILLlp_preop_free(ILLlp_preop * op);

void ILLlp_preline_init(ILLlp_preline * line);
void ILLlp_preline_free(ILLlp_preline * line);

#endif /* ! PRESOLVE_H */
