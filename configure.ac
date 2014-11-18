#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.69)
AC_INIT([QSopt_ex], [2.5.10], [daespino@gmail.com])
AC_CONFIG_SRCDIR([src/esolver.c])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign subdir-objects dist-xz])

m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

AC_USE_SYSTEM_EXTENSIONS

### Checks for programs.
AC_PROG_CC
AC_PROG_SED
AC_PROG_LIBTOOL

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
AC_CHECK_LIB([gmp], [__gmpz_init], [],
	[AC_MSG_ERROR([GNU MP not found, see https://gmplib.org/])])

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
"
