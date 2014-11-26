
ACLOCAL_AMFLAGS = -I m4

# These source files are used as is for the library
MAIN_SOURCE_FILES = \
	src/allocrus.c \
	src/bgetopt.c \
	src/eg_io.c \
	src/eg_lpnum.c \
	src/except.c \
	src/urandom.c \
	src/zeit.c \
	src/names.c \
	src/symtab.c \
	src/util.c \
	src/exact.c \
	src/reporter.c \
	src/eg_exutil.c \
	src/eg_macros.c \
	src/eg_memslab.c \
	src/sortrus_common.c \
	src/QSopt_ex_version.c

# Header files
MAIN_HEADER_FILES = \
	src/allocrus.h \
	src/bgetopt.h \
	src/eg_io.h \
	src/eg_lpnum.h \
	src/except.h \
	src/urandom.h \
	src/zeit.h \
	src/names.h \
	src/symtab.h \
	src/util.h \
	src/exact.h \
	src/reporter.h \
	src/eg_exutil.h \
	src/eg_macros.h \
	src/eg_memslab.h \
	src/eg_mem.h \
	src/eg_elist.h \
	src/eg_nummacros.h \
	src/eg_lpnum.dbl.h \
	src/eg_lpnum.mpq.h \
	src/eg_lpnum.mpf.h \
	src/eg_timer.h \
	src/basicdefs.h \
	src/stddefs.h \
	src/trace.h \
	src/sortrus_common.h \
	src/QSopt_ex_version.h \
	src/QSopt_ex.h

PRIVATE_HEADER_FILES = \
	src/qs_config.h

# These source files have to be instantiated using the
# type template mechanism
TEMPLATE_SOURCES = \
	src/rawlp.c \
	src/mps.c \
	src/read_mps.c \
	src/lp.c \
	src/write_lp.c \
	src/read_lp.c \
	src/readline.c \
	src/lpdata.c \
	src/presolve.c \
	src/factor.c \
	src/basis.c \
	src/price.c \
	src/dstruct.c \
	src/simplex.c \
	src/fct.c \
	src/ratio.c \
	src/lib.c \
	src/binary.c \
	src/qsopt.c \
	src/sortrus.c \
	src/dheaps_i.c \
	src/priority.c \
	src/editor.c \
	src/format.c \
	src/util.c \
	src/eg_numutil.c

TEMPLATE_HEADER_FILES = \
	src/rawlp.h \
	src/mps.h \
	src/read_mps.h \
	src/lp.h \
	src/write_lp.h \
	src/read_lp.h \
	src/readline.h \
	src/lpdata.h \
	src/presolve.h \
	src/factor.h \
	src/basis.h \
	src/price.h \
	src/dstruct.h \
	src/simplex.h \
	src/fct.h \
	src/ratio.h \
	src/lib.h \
	src/binary.h \
	src/qsopt.h \
	src/sortrus.h \
	src/dheaps_i.h \
	src/priority.h \
	src/editor.h \
	src/format.h \
	src/lpdefs.h \
	src/qstruct.h \
	src/solver.h \
	src/iqsutil.h \
	src/util.h \
	src/eg_numutil.h

# Template names to substitute
TEMPLATE_NAMES = \
	MaxLpNum MinLpNum epsLpNum EGlpNumGetStr EGlpNumReadStr EGlpNumCeil EGlpNumFloor \
        EGlpNumInv zeroLpNum oneLpNum EGlpNumSet EGlpNumIsEqual EGlpNumIsEqqual \
        EGlpNumIsNeq EGlpNumIsNeqq EGlpNumIsNeqZero EGlpNumIsNeqqZero EGlpNumIsLess \
        EGlpNumIsSumLess EGlpNumIsDiffLess EGlpNumIsLessDbl EGlpNnumIsGreaDbl \
        EGlpNumIsLeq EGlpNumCopyDiffRatio EGlpNumCopyDiff EGlpNumCopySum EGlpNumCopy \
        EGlpNumSetToMaxAbs EGlpNumSetToMinAbs EGlpNumCopySqrOver EGlpNumCopyAbs \
        EGlpNumCopyNeg EGlpNumCopyFrac EGlpNumCopyArray EGlpNumSubInnProdTo \
        EGlpNumAddInnProdTo EGlpNumSubUiTo EGlpNumAddUiTo EGlpNumAddTo EGlpNumSubTo \
        EGlpNumMultTo EGlpNumDivTo EGlpNumDivUiTo EGlpNumMultUiTo EGlpNumZero \
        EGlpNumOne EGlpNumSign EGlpNumToLf EGlpNumAllocArray EGlpNumInitVar \
        EGlpNumInitVar EGlpNumClearVar EGlpNumReallocArray EGlpNumFreeArray \
        EGlpNumIsGreaDbl EGlpNumIsLessDbl EGutilPermSort EGlpNumInnProd \
        EGlpNumIsLessZero EGlpNumIsGreatZero

# Template names that cannot be substituted
TEMPLATE_NAMES_GLOBAL = \
	main sscanf sprintf fprintf fprint ungetc perror parseargs \
        memset listen fscanf fflush fclose TRACE connect bind accept \
        time printf getrusage strcasecmp strncasecmp \
        NULL ILL_PTRWORLD_ROUTINES DEBUG ILL_IFTRACE itcnt_t

# Define specific template sources
TEMPLATE_SOURCES_DBL = $(TEMPLATE_SOURCES:.c=_dbl.c)
TEMPLATE_HEADERS_DBL = $(TEMPLATE_HEADER_FILES:.h=_dbl.h)

$(TEMPLATE_SOURCES_DBL): %_dbl.c: %.c template_dbl.sed Makefile
	$(AM_V_GEN)$(MKDIR_P) $(@D) && \
		$(SED) -f template_dbl.sed $< > $@
$(TEMPLATE_HEADERS_DBL): %_dbl.h: %.h template_dbl.sed Makefile
	$(AM_V_GEN)$(MKDIR_P) $(@D) && \
		$(SED) -f template_dbl.sed $< > $@

TEMPLATE_SOURCES_MPQ = $(TEMPLATE_SOURCES:.c=_mpq.c)
TEMPLATE_HEADERS_MPQ = $(TEMPLATE_HEADER_FILES:.h=_mpq.h)

$(TEMPLATE_SOURCES_MPQ): %_mpq.c: %.c template_mpq.sed Makefile
	$(AM_V_GEN)$(MKDIR_P) $(@D) && \
		$(SED) -f template_mpq.sed $< > $@
$(TEMPLATE_HEADERS_MPQ): %_mpq.h: %.h template_mpq.sed Makefile
	$(AM_V_GEN)$(MKDIR_P) $(@D) && \
		$(SED) -f template_mpq.sed $< > $@

TEMPLATE_SOURCES_MPF = $(TEMPLATE_SOURCES:.c=_mpf.c)
TEMPLATE_HEADERS_MPF = $(TEMPLATE_HEADER_FILES:.h=_mpf.h)

$(TEMPLATE_SOURCES_MPF): %_mpf.c: %.c template_mpf.sed Makefile
	$(AM_V_GEN)$(MKDIR_P) $(@D) && \
		$(SED) -f template_mpf.sed $< > $@
$(TEMPLATE_HEADERS_MPF): %_mpf.h: %.h template_mpf.sed Makefile
	$(AM_V_GEN)$(MKDIR_P) $(@D) && \
		$(SED) -f template_mpf.sed $< > $@

# Define complete set of template files
TEMPLATE_FILES = \
	$(TEMPLATE_SOURCES_DBL) $(TEMPLATE_HEADERS_DBL) \
	$(TEMPLATE_SOURCES_MPQ) $(TEMPLATE_HEADERS_MPQ) \
	$(TEMPLATE_SOURCES_MPF) $(TEMPLATE_HEADERS_MPF)

template_%.sed: templatetags Makefile typenames
	$(AM_V_GEN)$(MKDIR_P) $(@D) && \
		type_prefix="$$(echo $(@) | $(SED) -e 's|^template_\(.\+\)\.sed$$|\1|')" && \
		type_name="$$(grep ^$${type_prefix} $(top_srcdir)/typenames | cut -f2)" && \
		echo "# WARNING! This file was autogenerated by Makefile" > $@ && \
		( echo "1i\\" && echo "/* WARNING! This file was autogenerated from template */" && echo ) >> $@ && \
		echo "s|\<EGlpNum_t\>|$${type_name}|g" >> $@ && \
		for name in $$(cat templatetags); do \
			echo "s|\<$${name}\>|$${type_prefix}_&|g" \
		; done >> $@ ; \
		for header in $(TEMPLATE_HEADER_FILES); do \
			header_base="$$(basename $${header})" && \
			echo "s|\<$${header_base}\>|$${header_base%%.h}_$${type_prefix}.h|g" \
		; done >> $@

templatetags: $(TEMPLATE_HEADER_FILES)
	$(AM_V_GEN)$(MKDIR_P) $(@D) && \
		( ( $(XCTAGS) -x --c-kinds=+xp-m --c++-kinds=+xp-m $^ | \
		cut -d ' ' -f 1 ; \
		for name in $(TEMPLATE_NAMES); do echo "$${name}"; done ) | \
		sort | uniq ; \
		for name in $(TEMPLATE_NAMES_GLOBAL); do echo "$${name}"; echo "$${name}"; done ) | \
		sort | uniq -u > $@

# Tell automake that these are generated files
BUILT_SOURCES = $(TEMPLATE_FILES)

# Define tag files
TEMPLATE_SED_FILES = $(wildcard template_*.sed)


AM_CFLAGS = -I$(top_srcdir)/src -I$(top_builddir)/src

# Library files
lib_LTLIBRARIES = libqsopt_ex.la
libqsopt_ex_la_SOURCES = \
	$(MAIN_SOURCE_FILES) $(MAIN_HEADER_FILES) \
	$(PRIVATE_HEADER_FILES)
nodist_libqsopt_ex_la_SOURCES = $(TEMPLATE_FILES)
libqsopt_ex_la_CFLAGS = $(AM_CFLAGS) $(GMP_CFLAGS)
libqsopt_ex_la_LDFLAGS = \
	-no-undefined \
	-version-info $(LT_VERSION_INFO) \
	$(AM_LDFLAGS)
libqsopt_ex_la_LIBADD = $(GMP_LIBS)

# Library header files
headerdir = $(includedir)/qsopt_ex
header_DATA = \
	$(MAIN_HEADER_FILES) $(TEMPLATE_HEADERS_DBL) \
	$(TEMPLATE_HEADERS_MPQ) $(TEMPLATE_HEADERS_MPF)


# Program files
bin_PROGRAMS = esolver
noinst_PROGRAMS = demo_qs eg_sloan

demo_qs_SOURCES = src/demo_qs.c
demo_qs_CFLAGS = $(AM_CFLAGS) $(GMP_CFLAGS)
demo_qs_LDADD = $(GMP_LIBS) libqsopt_ex.la

eg_sloan_SOURCES = src/eg_sloan.c src/eg_sloan.h
eg_sloan_CFLAGS = $(AM_CFLAGS) $(GMP_CFLAGS)
eg_sloan_LDADD = $(GMP_LIBS) libqsopt_ex.la

esolver_SOURCES = src/esolver.c
esolver_CFLAGS = $(AM_CFLAGS) $(GMP_CFLAGS)
esolver_LDADD = $(GMP_LIBS) libqsopt_ex.la

# Python module
MODULE_PYX_FILES = \
	python/qsoptex.pyx

MODULE_PXD_FILES = \
	python/cgmp.pxd \
	python/cqsoptex.pxd

MODULE_C_FILES = $(MODULE_PYX_FILES:.pyx=.c)

$(MODULE_C_FILES): %.c: %.pyx Makefile
	$(AM_V_GEN)$(MKDIR_P) $(@D) && \
		$(CYTHON) -o $@ $<

if ENABLE_PYTHON_MODULE
pyexec_LTLIBRARIES = qsoptex.la
nodist_qsoptex_la_SOURCES = $(MODULE_C_FILES)
qsoptex_la_CFLAGS = $(AM_CFLAGS) $(GMP_CFLAGS) $(PYTHON_INCLUDES)
qsoptex_la_LIBADD = $(GMP_LIBS) libqsopt_ex.la
qsoptex_la_LDFLAGS = -avoid-version -module $(PYTHON_LDFLAGS)
endif


MODULE_TEST_FILES = python/test_qsoptex.py

if ENABLE_PYTHON_MODULE
# Run module unit tests after install
installcheck-local:
	PYTHONPATH="$(pyexecdir):$$PYTHONPATH" \
		$(PYTHON) $(top_srcdir)/python/test_qsoptex.py -v
endif


# Additional files to be distributed. Always run `make distcheck`
# to be sure all necessary files are distributed!
EXTRA_DIST = \
	README.md Doxyfile License.txt \
	$(TEMPLATE_SOURCES) $(TEMPLATE_HEADER_FILES) typenames \
	$(MODULE_PYX_FILES) $(MODULE_PXD_FILES) $(MODULE_TEST_FILES)

# Clean files
CLEANFILES = \
	$(TEMPLATE_FILES) $(TEMPLATE_SED_FILES) templatetags \
	$(MODULE_C_FILES)
