
QSopt Exact
===========

[![Build Status](https://travis-ci.org/jonls/qsopt-ex.svg?branch=master)](https://travis-ci.org/jonls/qsopt-ex)

Exact linear programming solver. This is a fork of QSopt_ex, originally
released by Daniel Espinoza _et al._
[version 2.5.10](http://www.math.uwaterloo.ca/~bico/qsopt/ex/) under the
GPL 3 (or later). The authors of QSopt_ex also granted a free license to
use the software for research purposes but this license does not extend
to the changes introduced by this project.

The goal of this fork is to update the software, and in particular the
build system, to be more friendly. In addition the external
dependencies have been reduced by removing the dependency on EGlib,
GNU awk and Exuberant Ctags.

Dependencies
------------

- C compiler: Tested with GCC and Clang.
- Libtool: To build QSopt_ex as a library.
- [GNU MP](https://gmplib.org/): Tested with 6.0.0. The original authors
  stated that QSopt_ex was tested with the 4.x.x and with 5.0.x version
  series. The authors also noted that GNU MP should be compiled using option
  `--enable-alloca=malloc-reentrant` but this does not seem to be required
  anymore.
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

To install the libraries and executables run

``` shell
$ make install
```

This will install into the prefix specified when `configure` was run.

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
To see an example of how to use this software as a C library, see the test
[tests/test_qs.c](tests/test_qs.c) or the program
[esolver/esolver.c](esolver/esolver.c).

Python module
-------------

The Python module has moved to a separate repository at
[jonls/python-qsoptex](https://github.com/jonls/python-qsoptex).
