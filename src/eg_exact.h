/* QSopt-Exact "An exact LP solver"
 *
 * Copyright (C) 2006 Daniel Espinoza.
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
 * */
/* ========================================================================= */
/** @mainpage QSopt-Exact Home Page
 *
 * @section Introduction
 
<P>This is a joint project of <A HREF=http://www.dii.uchile.cl/~daespino TARGET=_top>Daniel Espinoza</A>, <A HREF=http://www.isye.gatech.edu/~wcook TARGET=_top>William Cook</A>, <A HREF=http://www.research.ibm.com/people/s/sanjeebd TARGET=_top>Sanjeeb Dash</A> and <A HREF=http://public.research.att.com/viewPage.cfm?PageID=424 TARGET=_top>David Applegate</A>.
Also, <A HREF=http://www.zib.de/wolter/ TARGET=_top>Kati Wolter</A> has contributed with bug-fixes, bug-reports and special functionality for SCIP exact.</P>

<P>The objective of this software is to provide a solver for Linear Programming (and Integer Programming to a lesser degree)  that  returns true (rational) optimal solutions.</P>

<P>It relies heavilly on the <A HREF=http://www.swox.com/gmp TARGET=_top>GNUMP library</A>, that provides a multiprecision library for both floating point and also rational arithmetic.
Note that if you use a dynamicly linked version of QSopt-Exact, then the GMP library should have been compiled with the option --enable-alloca=malloc-reentrant, this is needed to avoid memory corruption. 
The basis for the LP solver was taken from <A HREF=http://www.isye.gatech.edu/~wcook/qsopt TARGET=_top>QSopt</A>, which is an LP solver based on floating point arithmetic and available for free for research purposes.
A brief description of the implementation and the obtained results can be obtained <A HREF=http://www.dii.uchile.cl/~daespino/files/exact_simplex.pdf >here (pdf)</A>, and a longer description (which is part of my Ph.D. thesis) can be found <A HREF=http://www.dii.uchile.cl/~daespino/files/espinoza_daniel_g_200605_phd.pdf >here</A>
Much of the functionality used in QSopt-Exact comes from <A HREF=http://www.dii.uchile.cl/~daespino/EGlib_doc/main.html TARGET=_self>EGlib</A>. 
</P>

<P>You can see the <A HREF=http://www.dii.uchile.cl/~daespino/QSoptExact_doc/modules.html TARGET=_self>documentation</A> or download the <A HREF=http://www.dii.uchile.cl/~daespino/SOurce/QSoptExact.tar.bz2 TARGET=_top>program source</A>.</P>
<P>We also have made available binaries for <A HREF=http://www.dii.uchile.cl/~daespino/SOurce/QSoptExact-32.tar.bz2 TARGET=_top>linux 32 bit</A> and for <A HREF=http://www.dii.uchile.cl/~daespino/SOurce/QSoptExact-64.tar.bz2 TARGET=_top>linux 64 bit</A>
 
<P>Finally, many thanks to all people that have contributed with bug-reports, comments, and help.</P>

 * */
