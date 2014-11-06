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

# ======================================================================
# Checks for programs.
# ======================================================================
AC_PROG_CC
AC_PROG_SED
AC_PROG_LIBTOOL

# ======================================================================
# Optional packages to use.
# ======================================================================
AC_ARG_WITH(
	[gmp-lib-dir],
	[AS_HELP_STRING([--with-gmp-lib-dir],[specify where to look for -lgmp])],
	[	AC_SUBST([LDFLAGS],["$LDFLAGS -L $with_gmp_lib_dir"]) 
		AC_SUBST([EXTRA_LIBS_DIR],["$EXTRA_LIBS_DIR $with_gmp_lib_dir"])])
AC_ARG_WITH(
	[gmp-include-dir],
	[AS_HELP_STRING([--with-gmp-include-dir],[specify where to look for gmp.h])],
	[	AC_SUBST([CPPFLAGS],["$CPPFLAGS -I $with_gmp_include_dir"])
		AC_SUBST([EXTRA_INCLUDE_DIR],["$EXTRA_INCLUDE_DIR $with_gmp_include_dir"])])
AC_ARG_WITH( [gmp], [AS_HELP_STRING( [--with-gmp], [support gnump templates @<:@default=check@:>@]) ], [], [with_gmp=check])
AS_IF(
	[test "x$with_gmp" != xno],
	[AC_CHECK_LIB(
		[gmp],
		[__gmpz_add],
		[AC_SUBST( [LIBGMP], ["-lgmp"])
		AC_SUBST([ENABLE_LIBGMP],[1])
		AC_DEFINE( [HAVE_LIBGMP], [1], [Define if you have libgmp])
		enable_libgmp=yes
		],[
		AC_SUBST([ENABLE_LIBGMP],[0])
		enable_libgmp=no
		if test "x$with_gmp" != xcheck; then
			AC_MSG_FAILURE([--with-gmp was given, but test for gmp failed]) fi ], ["-lgmp"])])
AM_CONDITIONAL([ENABLE_LIBGMP], [test "x$enable_libgmp" = xyes])

# ======================================================================
# gcc-based debug flags
# ======================================================================
AC_ARG_ENABLE([debug],[AS_HELP_STRING([--enable-debug],[use GCC debug comilation flags and support -ggdb3 @<:@default=yes@:>@])],[enable_debug=$enableval],[enable_debug=yes])
AS_IF(
	[test "x$enable_debug" != xno],
	[AC_SUBST([ENABLE_DEBUG],[1])],
	[AC_SUBST([ENABLE_DEBUG],[0])])

# ======================================================================
# compilation-time debug level
# ======================================================================
AC_ARG_ENABLE([debug-level],[AS_HELP_STRING([--enable-debug-level],[fix debug-level for compile-time routines @<:@default=1@:>@])],[enable_debug_level=$enableval],[enable_debug_level=1])
AC_DEFINE_UNQUOTED([DEBUG],[$enable_debug_level],[Compilation-time debug level])

# ======================================================================
# Check for long double support
# ======================================================================
AC_ARG_ENABLE([long-double],[AS_HELP_STRING([--enable-long-double],[support long double (if present in the system) @<:@default=yes@:>@])],[],[enable_long_double=yes])
AS_IF(
	[test "x$enable_long_double" != xno],
	[AC_DEFINE([ENABLE_LONG_DOUBLE],[1],[Define to one if you want long double support (if present in the system)])
	AC_SUBST([ENABLE_LONG_DOUBLE],[1])],
	[AC_SUBST([ENABLE_LONG_DOUBLE],[0])])

# ======================================================================
# gcc-based optimization flags
# ======================================================================
AC_ARG_ENABLE([optimize],[AS_HELP_STRING([--enable-optimize],[use GCC optimization comilation flags @<:@default=yes@:>@])],[],[enable_optimize=yes])
AS_IF(
	[test "x$enable_optimize" != xno],
	[AC_SUBST([ENABLE_OPTIMIZE],[1])],
	[AC_SUBST([ENABLE_OPTIMIZE],[0])])

# ======================================================================
# check for softfloat support
# ======================================================================
AC_ARG_ENABLE([softfloat],[AS_HELP_STRING([--enable-softfloat],[support software-based 128-bit floating point numbers @<:@default=no@:>@])],[],[enable_softfloat=no])
AS_IF(
	[test "x$enable_softfloat" != xno],
	[AC_DEFINE([HAVE_SOFTFLOAT],[1],[Define to one if you want software-based 128-bit floating point numbers])
	AC_SUBST([ENABLE_SOFTFLOAT],[1])
	enable_softfloat=yes
	],
	[AC_SUBST([ENABLE_SOFTFLOAT],[0])
	enable_softfloat=no])
AM_CONDITIONAL([ENABLE_SOFTFLOAT], [test "x$enable_softfloat" = xyes])

# ======================================================================
# compilation-time verbose level
# ======================================================================
AC_ARG_ENABLE([verbose-level],[AS_HELP_STRING([--enable-verbose-level],[fix verbose-level for compile-time routines @<:@default=100@:>@])],[],[enable_verbose_level=100])
AC_DEFINE_UNQUOTED([VERBOSE_LEVEL],[$enable_verbose_level],[Compilation-time verbose level])

# ======================================================================
# Checks for libraries.
# ======================================================================
AC_CHECK_LIB(z, gzopen, [], AC_MSG_NOTICE([Compiling without libz.]), [])
AC_CHECK_LIB(bz2, BZ2_bzopen, [], AC_MSG_NOTICE([Compiling without libbz2.]), [])

# find math library
AC_CHECK_LIB([m], [pow])

# find nsl library
AC_CHECK_LIB([nsl], [y_bind])

# thread library
AC_CHECK_LIB([pthread], [pthread_mutex_lock])

# resolve library
AC_CHECK_LIB([resolv], [_gethtbyname])

# ======================================================================
# Checking for supported compiler flags (taken from valgrind 3.5.0)
# ======================================================================
# does this compiler support -fforce-addr ?
AC_MSG_CHECKING([if gcc accepts -fforce-addr])
safe_CFLAGS=$CFLAGS
CFLAGS="-fforce-addr"
AC_TRY_COMPILE(,
	[ return 0; ],
	[ FLAG_FFORCE_ADDR="-fforce-addr"; AC_MSG_RESULT([yes]) ],
	[ FLAG_FFORCE_ADDR=""; AC_MSG_RESULT([no]) ])
CFLAGS=$safe_CFLAGS
AC_SUBST(FLAG_FFORCE_ADDR)

# ======================================================================
# does this compiler support -finline-limit=1000000 ?
AC_MSG_CHECKING([if gcc accepts -finline-limit=1000000])
safe_CFLAGS=$CFLAGS
CFLAGS="-finline-limit=1000000"
AC_TRY_COMPILE(,
	[ return 0; ],
	[ FLAG_FINLINE_LIMIT="-finline-limit=1000000"; AC_MSG_RESULT([yes]) ],
	[ FLAG_FINLINE_LIMIT=""; AC_MSG_RESULT([no]) ])
CFLAGS=$safe_CFLAGS
AC_SUBST(FLAG_FINLINE_LIMIT)

# ======================================================================
# does this compiler support -frerun-cse-after-loop ?
AC_MSG_CHECKING([if gcc accepts -frerun-cse-after-loop])
safe_CFLAGS=$CFLAGS
CFLAGS="-frerun-cse-after-loop"
AC_TRY_COMPILE(,
	[ return 0; ],
	[ FLAG_FRERUN_CSE_AFTER_LOOP="-frerun-cse-after-loop"; AC_MSG_RESULT([yes]) ],
	[ FLAG_FRERUN_CSE_AFTER_LOOP=""; AC_MSG_RESULT([no]) ])
CFLAGS=$safe_CFLAGS
AC_SUBST(FLAG_FRERUN_CSE_AFTER_LOOP)

# ======================================================================
# does this compiler support -frerun-loop-opt ?
AC_MSG_CHECKING([if gcc accepts -frerun-loop-opt])
safe_CFLAGS=$CFLAGS
CFLAGS="-frerun-loop-opt"
AC_TRY_COMPILE(,
	[ return 0; ],
	[ FLAG_FRERUN_LOOP_OPT="-frerun-loop-opt"; AC_MSG_RESULT([yes]) ],
	[ FLAG_FRERUN_LOOP_OPT=""; AC_MSG_RESULT([no]) ])
CFLAGS=$safe_CFLAGS
AC_SUBST(FLAG_FRERUN_LOOP_OPT)

# ======================================================================
# does this compiler support -funroll-loops ?
AC_MSG_CHECKING([if gcc accepts -funroll-loops])
safe_CFLAGS=$CFLAGS
CFLAGS="-funroll-loops"
AC_TRY_COMPILE(,
	[ return 0; ],
	[ FLAG_FUNROLL_LOOPS="-funroll-loops"; AC_MSG_RESULT([yes]) ],
	[ FLAG_FUNROLL_LOOPS=""; AC_MSG_RESULT([no]) ])
CFLAGS=$safe_CFLAGS
AC_SUBST(FLAG_FUNROLL_LOOPS)

# ======================================================================
# does this compiler support -ggdb3 ?
AC_MSG_CHECKING([if gcc accepts -ggdb3 -ggdb or -g])
safe_CFLAGS=$CFLAGS
CFLAGS="-ggdb3"
AC_TRY_COMPILE(,
	[ return 0; ],
	[ FLAG_GGDB3="-ggdb3"; AC_MSG_RESULT([-ggdb3]) ],
	[ CFLAGS="-ggdb"
		AC_TRY_COMPILE(,
			[ return 0; ],
			[ FLAG_GGDB3="-ggdb"; AC_MSG_RESULT([-ggdb])],
			[	CFLAGS="-g"
				AC_TRY_COMPILE(,
					[ return 0;],
					[ FLAG_GGDB3="-g"; AC_MSG_RESULT([-g])],
					[ FLAG_GGDB3=""; AC_MSG_RESULT([not supported])])
			])
	])
CFLAGS=$safe_CFLAGS
AC_SUBST(FLAG_GGDB3)

# ======================================================================
# does this compiler support -m32 ?
AC_MSG_CHECKING([if gcc accepts -m32])
safe_CFLAGS=$CFLAGS
CFLAGS="-m32"
AC_TRY_COMPILE(,
	[ return 0; ],
	[ FLAG_M32="-m32"; AC_MSG_RESULT([yes]) ],
	[ FLAG_M32=""; AC_MSG_RESULT([no]) ])
CFLAGS=$safe_CFLAGS
AC_SUBST(FLAG_M32)

# ======================================================================
# does this compiler support -m3dnow ?
AC_MSG_CHECKING([if gcc accepts -m3dnow])
safe_CFLAGS=$CFLAGS
CFLAGS="-m3dnow"
AC_TRY_COMPILE(, 
	[ return 0; ], 
	[ FLAG_M3DNOW="-m3dnow"; AC_MSG_RESULT([yes]) ],
	[ FLAG_M3DNOW=""; AC_MSG_RESULT([no]) ])
CFLAGS=$safe_CFLAGS
AC_SUBST(FLAG_M3DNOW)

# ======================================================================
# does this compiler support -maix32 ?
AC_MSG_CHECKING([if gcc accepts -maix32])
safe_CFLAGS=$CFLAGS
CFLAGS="-maix32"
AC_TRY_COMPILE(, 
	[ return 0; ], 
	[ FLAG_MAIX32="-maix32"; AC_MSG_RESULT([yes]) ],
	[ FLAG_MAIX32=""; AC_MSG_RESULT([no]) ])
CFLAGS=$safe_CFLAGS
AC_SUBST(FLAG_MAIX32)

# ======================================================================
# does this compiler support -m64 ?
AC_MSG_CHECKING([if gcc accepts -m64])
safe_CFLAGS=$CFLAGS
CFLAGS="-m64"
AC_TRY_COMPILE(, 
	[ return 0; ],
	[ FLAG_M64="-m64"; AC_MSG_RESULT([yes]) ],
	[ FLAG_M64=""; AC_MSG_RESULT([no]) ])
CFLAGS=$safe_CFLAGS
AC_SUBST(FLAG_M64)

# ======================================================================
# does this compiler support -maix64 ?
AC_MSG_CHECKING([if gcc accepts -maix64])
safe_CFLAGS=$CFLAGS
CFLAGS="-maix64"
AC_TRY_COMPILE(,
	[ return 0; ],
	[ FLAG_MAIX64="-maix64"; AC_MSG_RESULT([yes]) ],
	[ FLAG_MAIX64=""; AC_MSG_RESULT([no]) ])
CFLAGS=$safe_CFLAGS
AC_SUBST(FLAG_MAIX64)

# ======================================================================
# does this compiler support -march=native ?
AC_MSG_CHECKING([if gcc accepts -march=native])
safe_CFLAGS=$CFLAGS
CFLAGS="-march=native"
AC_TRY_COMPILE(,
	[ return 0; ],
	[ FLAG_MARCH="-march=native"; AC_MSG_RESULT([yes]) ],
	[ FLAG_MARCH=""; AC_MSG_RESULT([no]) ])
CFLAGS=$safe_CFLAGS
AC_SUBST(FLAG_MARCH)

# ======================================================================
# does this compiler support -mfpmath=sse ?
AC_MSG_CHECKING([if gcc accepts -mfpmath=sse])
safe_CFLAGS=$CFLAGS
CFLAGS="-mfpmath=sse"
AC_TRY_COMPILE(,
	[ return 0; ],
	[ FLAG_MFPMATH_SSE="-mfpmath=sse"; AC_MSG_RESULT([yes]) ],
	[ FLAG_MFPMATH_SSE=""; AC_MSG_RESULT([no]) ])
CFLAGS=$safe_CFLAGS
AC_SUBST(FLAG_MFPMATH_SSE)

# ======================================================================
# does this compiler support -mmmx ?
AC_MSG_CHECKING([if gcc accepts -mmmx])
safe_CFLAGS=$CFLAGS
CFLAGS="-mmmx"
AC_TRY_COMPILE(,
	[ return 0; ],
	[ FLAG_MMMX="-mmmx"; AC_MSG_RESULT([yes]) ],
	[ FLAG_MMMX=""; AC_MSG_RESULT([no]) ])
CFLAGS=$safe_CFLAGS
AC_SUBST(FLAG_MMMX)

# ======================================================================
# does this compiler support -mpreferred-stack-boundary=4 ?
AC_MSG_CHECKING([if gcc accepts -mpreferred-stack-boundary])
safe_CFLAGS=$CFLAGS
CFLAGS="-mpreferred-stack-boundary=4"
AC_TRY_COMPILE(, 
	[ return 0; ],
	[ FLAG_MPREFERRED_STACK_BOUNDARY="-mpreferred-stack-boundary=4";
		AC_MSG_RESULT([yes]) ],
	[ FLAG_MPREFERRED_STACK_BOUNDARY=""; AC_MSG_RESULT([no]) ])
CFLAGS=$safe_CFLAGS
AC_SUBST(FLAG_MPREFERRED_STACK_BOUNDARY)

# ======================================================================
# does this compiler support -msse ?
AC_MSG_CHECKING([if gcc accepts -msse])
safe_CFLAGS=$CFLAGS
CFLAGS="-msse"
AC_TRY_COMPILE(, 
	[ return 0; ], 
	[ FLAG_MSSE="-msse"; AC_MSG_RESULT([yes]) ],
	[ FLAG_MSSE=""; AC_MSG_RESULT([no]) ])
CFLAGS=$safe_CFLAGS
AC_SUBST(FLAG_MSSE)

# ======================================================================
# does this compiler support -msse2 ?
AC_MSG_CHECKING([if gcc accepts -msse2])
safe_CFLAGS=$CFLAGS
CFLAGS="-msse2"
AC_TRY_COMPILE(, 
	[ return 0; ], 
	[ FLAG_MSSE2="-msse2"; AC_MSG_RESULT([yes]) ],
	[ FLAG_MSSE2=""; AC_MSG_RESULT([no]) ])
CFLAGS=$safe_CFLAGS
AC_SUBST(FLAG_MSSE2)

# ======================================================================
# does this compiler support -msse3 ?
AC_MSG_CHECKING([if gcc accepts -msse3])
safe_CFLAGS=$CFLAGS
CFLAGS="-msse3"
AC_TRY_COMPILE(, 
	[ return 0; ], 
	[ FLAG_MSSE3="-msse3"; AC_MSG_RESULT([yes]) ],
	[ FLAG_MSSE3=""; AC_MSG_RESULT([no]) ])
CFLAGS=$safe_CFLAGS
AC_SUBST(FLAG_MSSE3)

# ======================================================================
# does this compiler support -msse4 ?
AC_MSG_CHECKING([if gcc accepts -msse4])
safe_CFLAGS=$CFLAGS
CFLAGS="-msse4"
AC_TRY_COMPILE(, 
	[ return 0; ], 
	[ FLAG_MSSE4="-msse4"; AC_MSG_RESULT([yes]) ],
	[ FLAG_MSSE4=""; AC_MSG_RESULT([no]) ])
CFLAGS=$safe_CFLAGS
AC_SUBST(FLAG_MSSE4)

# ======================================================================
# does this compiler support -mtune ?
AC_MSG_CHECKING([if gcc accepts -mtune=native])
safe_CFLAGS=$CFLAGS
CFLAGS="-mtune=native"
AC_TRY_COMPILE(, 
	[ return 0; ], 
	[ FLAG_MTUNE="-mtune=native"; AC_MSG_RESULT([yes]) ],
	[ FLAG_MTUNE=""; AC_MSG_RESULT([no]) ])
CFLAGS=$safe_CFLAGS
AC_SUBST(FLAG_MTUNE)

# ======================================================================
# does this compiler support -O3 ?
AC_MSG_CHECKING([if gcc accepts -O3])
safe_CFLAGS=$CFLAGS
CFLAGS="-O3"
AC_TRY_COMPILE(, 
	[ return 0; ],
	[ FLAG_O3="-O3";
		AC_MSG_RESULT([yes]) ],
	[ FLAG_O3=""; AC_MSG_RESULT([no]) ])
CFLAGS=$safe_CFLAGS
AC_SUBST(FLAG_O3)
# ======================================================================
# does this compiler support -std=gnu89 ?
AC_MSG_CHECKING([if gcc accepts -std=gnu89])
safe_CFLAGS=$CFLAGS
CFLAGS="-std=gnu89"
AC_TRY_COMPILE(, 
	[ return 0; ],
	[ FLAG_STD_GNU89="-std=gnu89";
		AC_MSG_RESULT([yes]) ],
	[ FLAG_STD_GNU89=""; AC_MSG_RESULT([no]) ])
CFLAGS=$safe_CFLAGS
AC_SUBST(FLAG_STD_GNU89)

# ======================================================================
# does this compiler support -Wall ?
AC_MSG_CHECKING([if gcc accepts -Wall])
safe_CFLAGS=$CFLAGS
CFLAGS="-Wall"
AC_TRY_COMPILE(,
	[ return 0; ],
	[	FLAG_WALL="-Wall"; AC_MSG_RESULT([yes]) ],
	[ FLAG_WALL=""; AC_MSG_RESULT([no]) ])
CFLAGS=$safe_CFLAGS
AC_SUBST(FLAG_WALL)

# ======================================================================
# does this compiler support -Wbad-function-cast ?
AC_MSG_CHECKING([if gcc accepts -Wbad-function-cast])
safe_CFLAGS=$CFLAGS
CFLAGS="-Wbad-function-cast"
AC_TRY_COMPILE(,
	[ return 0; ],
	[	FLAG_WBAD_FUNCTION_CAST="-Wbad-function-cast"; AC_MSG_RESULT([yes]) ],
	[ FLAG_WBAD_FUNCTION_CAST=""; AC_MSG_RESULT([no]) ])
CFLAGS=$safe_CFLAGS
AC_SUBST(FLAG_WBAD_FUNCTION_CAST)

# ======================================================================
# does this compiler support -Wcast-equal ?
AC_MSG_CHECKING([if gcc accepts -Wcast-equal])
safe_CFLAGS=$CFLAGS
CFLAGS="-Wcast-equal"
AC_TRY_COMPILE(,
	[ return 0; ],
	[	FLAG_WCAST_EQUAL="-Wundef"; AC_MSG_RESULT([yes]) ],
	[ FLAG_WCAST_EQUAL=""; AC_MSG_RESULT([no]) ])
CFLAGS=$safe_CFLAGS
AC_SUBST(FLAG_WCAST_EQUAL)

# ======================================================================
# does this compiler support -Wdeclaration-after-statement ?
AC_MSG_CHECKING([if gcc accepts -Wdeclaration-after-statement])
safe_CFLAGS=$CFLAGS
CFLAGS="-Wdeclaration-after-statement"
AC_TRY_COMPILE(,
	[ return 0; ],
	[	FLAG_WDECLARATION_AFTER_STATEMENT="-Wdeclaration-after-statement"; AC_MSG_RESULT([yes]) ],
	[ FLAG_WDECLARATION_AFTER_STATEMENT=""; AC_MSG_RESULT([no]) ])
CFLAGS=$safe_CFLAGS
AC_SUBST(FLAG_WDECLARATION_AFTER_STATEMENT)

# ======================================================================
# does this compiler support -Wextra or the older -W ?
AC_MSG_CHECKING([if gcc accepts -Wextra or -W])
safe_CFLAGS=$CFLAGS
CFLAGS="-Wextra"
AC_TRY_COMPILE(,
	[ return 0; ],
	[ AC_SUBST([FLAG_WEXTRA], [-Wextra])
		AC_MSG_RESULT([-Wextra]) ],
	[ CFLAGS="-W"
  	AC_TRY_COMPILE(,
			[ return 0; ],
			[ AC_SUBST([FLAG_WEXTRA], [-W])
				AC_MSG_RESULT([-W]) ],
			[ AC_SUBST([FLAG_WEXTRA], [])
				AC_MSG_RESULT([not supported]) ])
	])
CFLAGS=$safe_CFLAGS

# ======================================================================
# does this compiler support -Wfloat-equal ?
AC_MSG_CHECKING([if gcc accepts -Wfloat-equal])
safe_CFLAGS=$CFLAGS
CFLAGS="-Wfloat-equal"
AC_TRY_COMPILE(,
	[ return 0; ],
	[	FLAG_WFLOAT_EQUAL="-Wfloat-equal"; AC_MSG_RESULT([yes]) ],
	[ FLAG_WFLOAT_EQUAL=""; AC_MSG_RESULT([no]) ])
CFLAGS=$safe_CFLAGS
AC_SUBST(FLAG_WFLOAT_EQUAL)

# ======================================================================
# does this compiler support -Winline ?
AC_MSG_CHECKING([if gcc accepts -Winline])
safe_CFLAGS=$CFLAGS
CFLAGS="-Winline"
AC_TRY_COMPILE(,
	[ return 0; ],
	[	FLAG_WINLINE="-Winline"; AC_MSG_RESULT([yes]) ],
	[ FLAG_WINLINE=""; AC_MSG_RESULT([no]) ])
CFLAGS=$safe_CFLAGS
AC_SUBST(FLAG_WINLINE)

# ======================================================================
# does this compiler support -Wnexted-externs ?
AC_MSG_CHECKING([if gcc accepts -Wnexted-externs])
safe_CFLAGS=$CFLAGS
CFLAGS="-Wnexted-externs"
AC_TRY_COMPILE(,
	[ return 0; ],
	[	FLAG_WNESTED_EXTERNS="-Wnexted-externs"; AC_MSG_RESULT([yes]) ],
	[ FLAG_WNESTED_EXTERNS=""; AC_MSG_RESULT([no]) ])
CFLAGS=$safe_CFLAGS
AC_SUBST(FLAG_WNESTED_EXTERNS)

# ======================================================================
# does this compiler support -Wno-empty-body ?
AC_MSG_CHECKING([if gcc accepts -Wno-empty-body])
safe_CFLAGS=$CFLAGS
CFLAGS="-Wno-empty-body"
AC_TRY_COMPILE(,
	[ return 0; ],
	[	FLAG_WNO_EMPTY_BODY="-Wno-empty-body"; AC_MSG_RESULT([yes]) ],
	[ FLAG_WNO_EMPTY_BODY=""; AC_MSG_RESULT([no]) ])
CFLAGS=$safe_CFLAGS
AC_SUBST(FLAG_WNO_EMPTY_BODY)

# ======================================================================
# does this compiler support -Wno-pointer-sign ?
AC_MSG_CHECKING([if gcc accepts -Wno-pointer-sign])
safe_CFLAGS=$CFLAGS
CFLAGS="-Wno-pointer-sign"
AC_TRY_COMPILE(,
	[ return 0; ],
	[	FLAG_WNO_POINTER_SIGN="-Wno-pointer-sign"; AC_MSG_RESULT([yes]) ],
	[ FLAG_WNO_POINTER_SIGN=""; AC_MSG_RESULT([no]) ])
CFLAGS=$safe_CFLAGS
AC_SUBST(FLAG_WNO_POINTER_SIGN)

# ======================================================================
# does this compiler support -Wno-format-zero-length ?
AC_MSG_CHECKING([if gcc accepts -Wno-format-zero-length])
safe_CFLAGS=$CFLAGS
CFLAGS="-Wno-format-zero-length"
AC_TRY_COMPILE(,
	[ return 0; ],
	[	FLAG_WNO_FORMAT_ZERO_LENGTH="-Wno-format-zero-length"; AC_MSG_RESULT([yes]) ],
	[ FLAG_WNO_FORMAT_ZERO_LENGTH=""; AC_MSG_RESULT([no]) ])
CFLAGS=$safe_CFLAGS
AC_SUBST(FLAG_WNO_FORMAT_ZERO_LENGTH)

# ======================================================================
# does this compiler support -Wno-uninitialized ?
AC_MSG_CHECKING([if gcc accepts -Wno-uninitialized])
safe_CFLAGS=$CFLAGS
CFLAGS="-Wno-uninitialized"
AC_TRY_COMPILE(,
	[ return 0; ],
	[	FLAG_WNO_UNINITIALIZED="-Wno-uninitialized"; AC_MSG_RESULT([yes]) ],
	[ FLAG_WNO_UNINITIALIZED=""; AC_MSG_RESULT([no]) ])
CFLAGS=$safe_CFLAGS
AC_SUBST(FLAG_WNO_UNINITIALIZED)

# ======================================================================
# does this compiler support -Wpointer-arith ?
AC_MSG_CHECKING([if gcc accepts -Wpointer-arith])
safe_CFLAGS=$CFLAGS
CFLAGS="-Wpointer-arith"
AC_TRY_COMPILE(,
	[ return 0; ],
	[	FLAG_WPOINTER_ARITH="-Wpointer-arith"; AC_MSG_RESULT([yes]) ],
	[ FLAG_WPOINTER_ARITH=""; AC_MSG_RESULT([no]) ])
CFLAGS=$safe_CFLAGS
AC_SUBST(FLAG_WPOINTER_ARITH)

# ======================================================================
# does this compiler support -Wshadow ?
AC_MSG_CHECKING([if gcc accepts -Wshadow])
safe_CFLAGS=$CFLAGS
CFLAGS="-Wshadow"
AC_TRY_COMPILE(,
	[ return 0; ],
	[	FLAG_WSHADOW="-Wshadow"; AC_MSG_RESULT([yes]) ],
	[ FLAG_WSHADOW=""; AC_MSG_RESULT([no]) ])
CFLAGS=$safe_CFLAGS
AC_SUBST(FLAG_WSHADOW)

# ======================================================================
# does this compiler support -Wsign-compare ?
AC_MSG_CHECKING([if gcc accepts -Wsign-compare])
safe_CFLAGS=$CFLAGS
CFLAGS="-Wsign-compare"
AC_TRY_COMPILE(,
	[ return 0; ],
	[	FLAG_WSIGN_COMPARE="-Wsign-compare"; AC_MSG_RESULT([yes]) ],
	[ FLAG_WSIGN_COMPARE=""; AC_MSG_RESULT([no]) ])
CFLAGS=$safe_CFLAGS
AC_SUBST(FLAG_WSIGN_COMPARE)

# ======================================================================
# does this compiler support -Wstrict-prototypes ?
AC_MSG_CHECKING([if gcc accepts -Wstrict-prototypes])
safe_CFLAGS=$CFLAGS
CFLAGS="-Wstrict-prototypes"
AC_TRY_COMPILE(,
	[ return 0; ],
	[	FLAG_WSTRICT_PROTOTYPES="-Wstrict-prototypes"; AC_MSG_RESULT([yes]) ],
	[ FLAG_WSTRICT_PROTOTYPES=""; AC_MSG_RESULT([no]) ])
CFLAGS=$safe_CFLAGS
AC_SUBST(FLAG_WSTRICT_PROTOTYPES)

# ======================================================================
# does this compiler support -Wundef ?
AC_MSG_CHECKING([if gcc accepts -Wundef])
safe_CFLAGS=$CFLAGS
CFLAGS="-Wundef"
AC_TRY_COMPILE(,
	[ return 0; ],
	[	FLAG_WUNDEF="-Wundef"; AC_MSG_RESULT([yes]) ],
	[ FLAG_WUNDEF=""; AC_MSG_RESULT([no]) ])
CFLAGS=$safe_CFLAGS
AC_SUBST(FLAG_WUNDEF)

# ======================================================================
# does this compiler support -Wuninitialized ?
AC_MSG_CHECKING([if gcc accepts -Wuninitialized])
safe_CFLAGS=$CFLAGS
CFLAGS="-Wuninitialized"
AC_TRY_COMPILE(,
	[ return 0; ],
	[	FLAG_WUNINITIALIZED="-Wuninitialized"; AC_MSG_RESULT([yes]) ],
	[ FLAG_WUNINITIALIZED=""; AC_MSG_RESULT([no]) ])
CFLAGS=$safe_CFLAGS
AC_SUBST(FLAG_WUNINITIALIZED)

# ======================================================================
# does this compiler support -Wunused-function ?
AC_MSG_CHECKING([if gcc accepts -Wunused-function])
safe_CFLAGS=$CFLAGS
CFLAGS="-Wunused-function"
AC_TRY_COMPILE(,
	[ return 0; ],
	[	FLAG_WUNUSED_FUNCTION="-Wunused-function"; AC_MSG_RESULT([yes]) ],
	[ FLAG_WUNUSED_FUNCTION=""; AC_MSG_RESULT([no]) ])
CFLAGS=$safe_CFLAGS
AC_SUBST(FLAG_WUNUSED_FUNCTION)

# ======================================================================
# does this compiler support -Wunused-label ?
AC_MSG_CHECKING([if gcc accepts -Wunused-label])
safe_CFLAGS=$CFLAGS
CFLAGS="-Wunused-label"
AC_TRY_COMPILE(,
	[ return 0; ],
	[	FLAG_WUNUSED_LABEL="-Wunused-label"; AC_MSG_RESULT([yes]) ],
	[ FLAG_WUNUSED_LABEL=""; AC_MSG_RESULT([no]) ])
CFLAGS=$safe_CFLAGS
AC_SUBST(FLAG_WUNUSED_LABEL)

# ======================================================================
# does this compiler support -Wunused-value ?
AC_MSG_CHECKING([if gcc accepts -Wunused-value])
safe_CFLAGS=$CFLAGS
CFLAGS="-Wunused-value"
AC_TRY_COMPILE(,
	[ return 0; ],
	[	FLAG_WUNUSED_VALUE="-Wunused-value"; AC_MSG_RESULT([yes]) ],
	[ FLAG_WUNUSED_VALUE=""; AC_MSG_RESULT([no]) ])
CFLAGS=$safe_CFLAGS
AC_SUBST(FLAG_WUNUSED_VALUE)

# ======================================================================
# does this compiler support -Wunused-variable ?
AC_MSG_CHECKING([if gcc accepts -Wunused-variable])
safe_CFLAGS=$CFLAGS
CFLAGS="-Wunused-variable"
AC_TRY_COMPILE(,
	[ return 0; ],
	[	FLAG_WUNUSED_VARIABLE="-Wunused-variable"; AC_MSG_RESULT([yes]) ],
	[ FLAG_WUNUSED_VARIABLE=""; AC_MSG_RESULT([no]) ])
CFLAGS=$safe_CFLAGS
AC_SUBST(FLAG_WUNUSED_VARIABLE)

# ======================================================================
# does this compiler support -Wunused-parameter ?
AC_MSG_CHECKING([if gcc accepts -Wunused-parameter])
safe_CFLAGS=$CFLAGS
CFLAGS="-Wunused-parameter"
AC_TRY_COMPILE(,
	[ return 0; ],
	[	FLAG_WUNUSED_PARAMETER="-Wunused-parameter"; AC_MSG_RESULT([yes]) ],
	[ FLAG_WUNUSED_PARAMETER=""; AC_MSG_RESULT([no]) ])
CFLAGS=$safe_CFLAGS
AC_SUBST(FLAG_WUNUSED_PARAMETER)

# ======================================================================
# does this compiler support -Wwrite-strings ?
AC_MSG_CHECKING([if gcc accepts -Wwrite-strings])
safe_CFLAGS=$CFLAGS
CFLAGS="-Wwrite-strings"
AC_TRY_COMPILE(,
	[ return 0; ],
	[	FLAG_WWRITE_STRINGS="-Wwrite-strings"; AC_MSG_RESULT([yes]) ],
	[ FLAG_WWRITE_STRINGS=""; AC_MSG_RESULT([no]) ])
CFLAGS=$safe_CFLAGS
AC_SUBST(FLAG_WWRITE_STRINGS)

# ======================================================================
# Checks for header files.
AC_HEADER_STDC
AC_HEADER_TIME
AC_CHECK_HEADERS(sys/time.h time.h stdio.h sys/types.h sys/stat.h stdlib.h stddef.h string.h strings.h inttypes.h stdint.h unistd.h limits.h errno.h assert.h malloc.h fcntl.h signal.h sys/socket.h netinet/in.h netdb.h sys/resource.h sys/param.h sys/times.h float.h zlib.h bzlib.h stdarg.h gmp.h math.h sys/utsname.h setjmp.h)
AC_DEFUN([EG_INCLUDES],[
		#define _XOPEN_SOURCE 600
		#ifdef TIME_WITH_SYS_TIME
		# include <sys/time.h>
		# include <time.h>
		#else
		# ifdef HAVE_SYS_TIME_H
		# include <sys/time.h>
		# else
		# include <time.h>
		# endif
		#endif
		#ifdef HAVE_STDIO_H
		#include <stdio.h>
		#endif
		#ifdef HAVE_SYS_TYPES_H
		# include <sys/types.h>
		#endif
		#ifdef HAVE_SYS_STAT_H
		# include <sys/stat.h>
		#endif
		#ifdef STDC_HEADERS
		# include <stdlib.h>
		# include <stddef.h>
		#else
		# ifdef HAVE_STDLIB_H
		# include <stdlib.h>
		# endif
		# ifdef HAVE_STDDEF_H
		# include <stddef.h>
		#endif
		#endif
		#ifdef HAVE_STRING_H
		# if !defined STDC_HEADERS && defined HAVE_MEMORY_H
		# include <memory.h>
		# endif
		# include <string.h>
		#endif
		#ifdef HAVE_STRINGS_H
		# include <strings.h>
		#endif
		#ifdef HAVE_INTTYPES_H
		# include <inttypes.h>
		#endif
		#ifdef HAVE_STDINT_H
		# include <stdint.h>
		#endif
		#ifdef HAVE_UNISTD_H
		# include <unistd.h>
		#endif
		#ifdef HAVE_LIMITS_H
		# include <limits.h>
		#endif
		#ifdef HAVE_ERRNO_H
		# include <errno.h>
		#endif
		#ifdef HAVE_ASSERT_H
		# include <assert.h>
		#endif
		#ifdef HAVE_MALLOC_H
		# include <malloc.h>
		#endif
		#ifdef HAVE_FCNTL_H
		# include <fcntl.h>
		#endif
		#ifdef HAVE_SIGNAL_H
		# include <signal.h>
		#endif
		#ifdef HAVE_SYS_SOCKET_H
		# include <sys/socket.h>
		#endif
		#ifdef HAVE_NETINET_IN_H
		# include <netinet/in.h>
		#endif
		#ifdef HAVE_NETDB_H
		# include <netdb.h>
		#endif
		#ifdef HAVE_SYS_RESOURCE_H
		# include <sys/resource.h>
		#endif
		#ifdef HAVE_SYS_PARAM_H
		# include <sys/param.h>
		#endif
		#ifdef HAVE_SYS_TIMES_H
		# include <sys/times.h>
		#endif
		#ifdef HAVE_FLOAT_H
		# include <float.h>
		#endif
		#ifdef HAVE_ZLIB_H
		# include <zlib.h>
		#endif
		#ifdef HAVE_BZLIB_H
		# include <bzlib.h>
		#endif
		#ifdef HAVE_STDARG_H
		# include <stdarg.h>
		#endif
		#ifdef HAVE_GMP_H
		# include <gmp.h>
		#endif
		#ifdef HAVE_MATH_H
		# include <math.h>
		#endif
		#ifdef HAVE_SYS_UTSNAME_H
		# include <sys/utsname.h>
		#endif
		#ifdef HAVE_SETJMP_H
		# include <setjmp.h>
		#endif
		])

# ======================================================================
# Checks for typedefs, structures, and compiler characteristics.
# ======================================================================
AC_C_CONST
AC_C_INLINE
AC_C_VOLATILE
AC_C_RESTRICT
AC_C_TYPEOF
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_INT8_T
AC_TYPE_PID_T
AC_TYPE_SIGNAL
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T
AC_TYPE_LONG_DOUBLE

# Check to enable long double
AS_IF([test "x$ac_cv_type_long_double" != xno], [
	AC_SUBST([HAVE_LONG_DOUBLE],[1])
	enable_ldbl=yes
], [
	AC_SUBST([HAVE_LONG_DOUBLE],[0])
	enable_ldbl=no
])
AM_CONDITIONAL([ENABLE_LONG_DOUBLE], [test "x$enable_ldbl" = xyes])

AC_CHECK_TYPES([gzFile, BZFILE*],[],[],[EG_INCLUDES])

# ======================================================================
# Check declaration of some functions
# ======================================================================
AC_CHECK_DECLS([gzopen, gzeof, gzerror, gzclose, gzgets, gzwrite],[],[],[EG_INCLUDES])
AC_CHECK_DECLS([BZ2_bzopen, BZ2_bzerror, BZ2_bzclose, BZ2_bzread, BZ2_bzwrite],[],[],[EG_INCLUDES])

# ======================================================================
# Checks for library functions.
# ======================================================================
AC_FUNC_MALLOC
AC_FUNC_MEMCMP
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_FUNC_VPRINTF
AC_CHECK_FUNCS([floor gethostbyname memset pow socket sqrt strdup strerror uname posix_memalign sleep getrusage times clock sigaction signal])

# ======================================================================
# Test that GMP is available
# ======================================================================
AS_IF([test \(\("x$with_gmp" = xyes\) -a \("x$ac_cv_have_header_gmp_h" != xyes\)\)],[AC_ERROR([You need gmp.h (usually in gmp-dev package) in order to compile gmp support for EGlib, or specify the correct include path])])

# ======================================================================
# Print output files
# ======================================================================
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
