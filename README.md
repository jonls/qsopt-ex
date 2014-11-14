
NEEDED SOFTWARE:

- GCC c compiler (gcc), porting to other C compilers have not been attempted 
	(if you succeed with other compilers please let me know).
- gawk (other versions of awk don't work, be advised), it has been tested with 
	GNU awk version 3.1.6
- exuberant ctags (regular ctags don't work), it has been tested with 
	Exuberant Ctags 5.6, Copyright (C) 1996-2004 Darren Hiebert
	Addresses: <dhiebert@users.sourceforge.net>, http://ctags.sourceforge.net
- GNU MP (we have tested the 4.x.x and with 5.0.x version series without problem ), be aware 
	that you should compile and install a version using option 
	--enable-alloca=malloc-reentrant	
	to ensure no memory overwriting problems.
- EGlib you will need version 2.6.20 or later, and is available as a
	subversion repository in 
	https://conexo.dii.uchile.cl/SVN/EGlib/EGlib2/tags/EGlib-2.6.20/
	or you can get the bleading edge version at
	https://conexo.dii.uchile.cl/SVN/EGlib/EGlib2/trunk
	or you can get the full source at
	http://www.dii.uchile.cl/~daespino/
	Under EGlib, note that you should have installed before the previous programs,
	and you should ensure that Makefile.common uses gawk and exuberant-ctags, also 
	you should edit make.conf and enable GMP support and SoftFloat support.
- libz to read/write gz-compresed files
- libbz2 to read/write bz2-compresed files

INSTALLING

	After installing all pre-requisites, ensure that Makefile.common uses 
	gawk and exuberant-ctags, also you should edit make.conf to ensure that 
	proper paths are suplied, an example make.conf.default is provided.

	Set-up path and locations of needed software (see --help for options)

	./configure

	Then just type 

		make
	you will generate several executables. The main solver
	is named esolver or esolver_dyn (depending on compilation options), to force 
	the generation of a static binary  type 
		make -f Makefile.library esolver 

USING IT AS A LIBRARY
	To see an example of how to use this software as a C library, see the file
	src/esolver.c or the file SLoan_LPs/eg_sloan.h, the library is provided in
	shared and static form, they are mamed lib/QSopt_ex.so and lib/QSopt_ex.a 
	and the main include file is called include/QSopt_ex.h
	A simple example showing the basic functions and details of using mpq_t types 
	is shown in src/demo_qs.c

NOTES 
	You could use the release version that is MUCH simpler to compile, se the 
	webpage and https://conexo.dii.uchile.cl/SVN/ESolver/tags/QSopt_ex-2.5.8

COMMENTS

	As usual, no waranties are made about the software, see the LICNECE file.
	For comments or questions send an e-mail to daespino __at__ gmail __dot__ com

ISSUES/KNOWN BUGS

	It seems that reading LPs with maximizing objective functions is broken, check if that is the case for your problem
