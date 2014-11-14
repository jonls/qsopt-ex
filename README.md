
QSopt Exact
===========

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

- C compiler: the README in the original QSopt_ex notes that GCC is
  required and that porting to other C compilers had not been
  attempted. This could indicate that QSopt_ex is using GCC specific
  extensions (this has not yet been verified).
- Libtool: To build QSopt_ex as a library.
- Exuberant Ctags: Regular Ctags won't work according to the QSopt_ex
  authors. This program is needed for the custom templating build
  system that was implemented in QSopt_ex. This system may be
  refined in a future update, removing this dependency.
- GNU MP: (original QSopt_ex was tested with the 4.x.x and with 5.0.x
  version series without problem, according to the authors). The
  authors also note that GNU MP should be compiled using option
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
arguments. This will also regenerate the `configure` script.

``` shell
$ ./configure
```

Use `./configure --help` to see available options. Now the test
programs and library can be compiled using

``` shell
$ make
```

Using it as a library
---------------------
To see an example of how to use this software as a C library, see the file
`src/esolver.c`.
