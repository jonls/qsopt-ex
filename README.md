
QSopt Exact
===========

[![Build Status](https://travis-ci.org/jonls/redshift.svg?branch=master)](https://travis-ci.org/jonls/redshift)

Exact linear programming solver. This is a fork of QSopt_ex by Daniel
Espinoza et al. version 2.5.10 released under the LGPL 2.1
(http://www.math.uwaterloo.ca/~bico/qsopt/ex/). The authors of
QSopt_ex also granted a free license to use the software for research
purposes but this license does not extend to the changes introduced by
this project.

The goal of this fork is to update the software, and in particular the
build system, to be more friendly. In addition the external
dependencies have been reduced by removing the dependency on EGlib and
GNU awk. The dependencies may be further reduced later.

Dependencies
------------

- C compiler: Tested with GCC and Clang.
- Libtool: To build QSopt_ex as a library.
- Exuberant Ctags: Regular Ctags won't work according to the QSopt_ex
  authors. This program is needed for the custom templating build
  system that was implemented in QSopt_ex. This system may be
  refined in a future update, removing this dependency.
- GNU MP: Tested with 6.0.0. The original authors stated that QSopt_ex
  was tested with the 4.x.x and with 5.0.x version series. The authors
  also noted that GNU MP should be compiled using option
  `--enable-alloca=malloc-reentrant` but this does not seem to be
  required anymore.
- libz: To read/write gz-compresed files.
- libbz2 To read/write bz2-compresed files.

Installing
----------

If you have just cloned the source code with Git, run the `bootstrap`
script to automatically set up the build system.

``` shell
$ ./bootstrap
```

This script calls `autoreconf` and `libtoolize` with the proper
arguments. This will also regenerate the `configure` script. It is
recommended to build out of source directory. This is simply done
by running `configure` from an empty directory.

``` shell
$ mkdir build && cd build
$ ../configure
```

Use `./configure --help` to see available options. Now the test
programs and library can be compiled using `make`. It is possible
to do a parallel build using the `-jX` switch where `X` is the number
of parallel processes.

``` shell
$ make -j4
```

Running the solver
------------------

The exact solver is available though the `esolver` executable. It can be
invoked to solve an LP or MPS format problem.

``` shell
$ ./esolver cycle.mps
```

See `./esolver -h` for more information on command line options.

Using it as a library
---------------------
To see an example of how to use this software as a C library, see the file
`src/esolver.c`.
