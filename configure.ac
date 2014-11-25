#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.69)
AC_INIT([QSopt_ex], [2.5.10], [jonlst@gmail.com])
AC_CONFIG_SRCDIR([src/esolver.c])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign subdir-objects dist-xz])

m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

AC_USE_SYSTEM_EXTENSIONS

### Versioning
# This is the internal library version information used by libtool.
# IT DOES NOT FOLLOW THE NORMAL RELEASE VERSIONING!
# Instead follow these rules from the libtool documentation:
# 1. Update the version information only immediately before a public release of your software.
# 2. If the library source code has changed at all since the last update, then increment revision.
#    ("c:r:a" becomes "c:r+1:a")
# 3. If any interfaces have been added, removed, or changed since the last update, increment current,
#    and set revision to 0. ("c+1:0:a")
# 4. If any interfaces have been added since the last public release, then increment age.
# 5. If any interfaces have been removed since the last public release, then set age to 0.

m4_define([lt_current],[1])
m4_define([lt_revision],[1])
m4_define([lt_age],[1])
m4_define([lt_version_info],[lt_current:lt_revision:lt_age])
AC_SUBST([LT_VERSION_INFO], [lt_version_info])

### Checks for programs.
AC_PROG_CC
AC_PROG_SED
AC_PROG_LIBTOOL

AM_PATH_PYTHON([2.7], [have_python=yes], [have_python=no])
AC_PATH_PROGS([CYTHON], [cython cython-$PYTHON_VERSION], [:])
AS_IF([test "$CYTHON" != :], [have_cython=yes], [have_cython=no])

### Check for Exuberant ctags
AC_ARG_VAR([XCTAGS], [Exuberant Ctags command])
AS_IF([test -z "$XCTAGS"], [
	AC_PATH_PROGS([XCTAGS], [ctags-exuberant ctags], [:])
	AC_MSG_CHECKING([whether ctags progam is Exuberant Ctags])
	AS_IF([test "$XCTAGS" != :], [
		AS_IF([$XCTAGS --version 2>&1 | grep -q "Exuberant"], [
			have_xctags=yes
			AC_MSG_RESULT([yes])
		], [
			have_xctags=no
			AC_MSG_RESULT([no])
		])
	], [
		have_xctags=no
		AC_MSG_RESULT([no])
	])
], [
  have_xctags=yes
])
AS_IF([test "x$have_xctags" != xyes], [AC_MSG_ERROR([Exuberant Ctags not found, see http://ctags.sourceforge.net/])])

### Check for flags to use for python modules
PYTHON_INCLUDES=""
PYTHON_LDFLAGS=""
AS_IF([test "$PYTHON" != :], [
	PYTHON_CONFIG="$PYTHON-config"
	AS_IF([test -x "$PYTHON_CONFIG"], [
		PYTHON_INCLUDES=`$PYTHON_CONFIG --includes`
		PYTHON_LDFLAGS=`$PYTHON_CONFIG --ldflags`
	], [
		AC_MSG_ERROR([Unable to find python-config])
	])
])
AC_SUBST(PYTHON_INCLUDES)
AC_SUBST(PYTHON_LDFLAGS)

### Check for Python module support
AC_MSG_CHECKING([whether to enable Python module])
AC_ARG_ENABLE([python-module], [AC_HELP_STRING([--enable-python-module],
	[enable Python module (qsoptex)])],
	[enable_python_module=$enableval], [enable_python_module=maybe])
AS_IF([test "x$enable_python_module" != xno], [
	AS_IF([test "x$have_python" = xyes -a "x$have_cython" = xyes], [
		AC_MSG_RESULT([yes])
		enable_python_module=yes
	], [
		AC_MSG_RESULT([missing Python and/or Cython])
		AS_IF([test "x$enable_python_module" = xyes], [
			AC_MSG_ERROR([missing dependencies for Python module support])
		])
		enable_python_module=no
	])
], [
	AC_MSG_RESULT([no])
	enable_python_module=no
])
AM_CONDITIONAL([ENABLE_PYTHON_MODULE], [test "x$enable_python_module" = "xyes"])

### Check for debug mode
AC_ARG_ENABLE([debug],
	[AS_HELP_STRING([--enable-debug], [enable debug mode @<:@default=yes@:>@])],
	[enable_debug=$enable])
AS_IF([test "x$enable_debug" != xno], [
	AC_SUBST([ENABLE_DEBUG],[1])
	AC_DEFINE([DEBUG], [1], [Enable debug mode])
], [
  AC_SUBST([ENABLE_DEBUG],[0])
	AC_DEFINE([DEBUG], [0], [Disable debug mode])
])

### Check for verbosity level
AC_ARG_ENABLE([verbose-level],
	[AS_HELP_STRING([--enable-verbose-level], [fix verbose-level for compile-time routines @<:@default=100@:>@])],
	[],[enable_verbose_level=100])
AC_DEFINE_UNQUOTED([VERBOSE_LEVEL],[$enable_verbose_level],[Compilation-time verbose level])

### Check for libgmp
AC_ARG_VAR([GMP_CFLAGS], [C compiler flags for GMP, overriding automatic detection])
AC_ARG_VAR([GMP_LIBS], [linker flags for GMP, overriding automatic detection])
AS_IF([test -z "$GMP_CFLAGS" -a -z "$GMP_LIBS"], [
	AC_CHECK_LIB([gmp], [__gmpz_init], [
		GMP_CFLAGS=""
		GMP_LIBS="-lgmp"
	], [AC_MSG_ERROR([GNU MP not found, see https://gmplib.org/])])
])

### Check for libz
AC_CHECK_LIB([z], [gzopen], [], AC_MSG_NOTICE([compiling without libz]))

### Check for libbz2
AC_CHECK_LIB([bz2], [BZ2_bzopen], [], AC_MSG_NOTICE([compiling without libbz2]))

### Check for math library
AC_SEARCH_LIBS([pow], [m], [], [
	AC_MSG_ERROR([unable to find library with math functions (pow)])])

### Check for pthread library
AC_CHECK_LIB([pthread], [pthread_mutex_lock])

### Checks for header files.
AC_HEADER_STDC
AC_HEADER_TIME
AC_CHECK_HEADERS([errno.h float.h getopt.h inttypes.h limits.h math.h setjmp.h signal.h stdarg.h stdint.h stdio.h string.h sys/param.h sys/resource.h sys/time.h sys/times.h sys/utsname.h unistd.h])

### Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_TYPEOF
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_INT8_T
AC_TYPE_SIGNAL
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T

### Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_MEMCMP
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_FUNC_VPRINTF
AC_CHECK_FUNCS([floor memset pow sqrt strdup strerror uname posix_memalign sleep getrusage times clock sigaction signal])

### Print output files
AC_CONFIG_FILES([
	Makefile
])
AC_OUTPUT


echo "
 $PACKAGE_NAME $VERSION

    prefix:		${prefix}
    compiler:		${CC}
    cflags:		${CFLAGS}
    ldflags:		${LDFLAGS}

    Python module:	${enable_python_module}
"
