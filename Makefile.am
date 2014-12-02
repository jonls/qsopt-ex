
ACLOCAL_AMFLAGS = -I m4

# These source files are used as is for the library
MAIN_SOURCE_FILES = \
	qsopt_ex/allocrus.c \
	qsopt_ex/bgetopt.c \
	qsopt_ex/eg_io.c \
	qsopt_ex/eg_lpnum.c \
	qsopt_ex/except.c \
	qsopt_ex/urandom.c \
	qsopt_ex/zeit.c \
	qsopt_ex/names.c \
	qsopt_ex/symtab.c \
	qsopt_ex/util.c \
	qsopt_ex/exact.c \
	qsopt_ex/reporter.c \
	qsopt_ex/eg_exutil.c \
	qsopt_ex/eg_macros.c \
	qsopt_ex/eg_memslab.c \
	qsopt_ex/sortrus_common.c \
	qsopt_ex/QSopt_ex_version.c

# Header files
MAIN_HEADER_FILES = \
	qsopt_ex/allocrus.h \
	qsopt_ex/bgetopt.h \
	qsopt_ex/eg_io.h \
	qsopt_ex/eg_lpnum.h \
	qsopt_ex/except.h \
	qsopt_ex/urandom.h \
	qsopt_ex/zeit.h \
	qsopt_ex/names.h \
	qsopt_ex/symtab.h \
	qsopt_ex/util.h \
	qsopt_ex/exact.h \
	qsopt_ex/reporter.h \
	qsopt_ex/eg_exutil.h \
	qsopt_ex/eg_macros.h \
	qsopt_ex/eg_memslab.h \
	qsopt_ex/eg_mem.h \
	qsopt_ex/eg_elist.h \
	qsopt_ex/eg_nummacros.h \
	qsopt_ex/eg_lpnum.dbl.h \
	qsopt_ex/eg_lpnum.mpq.h \
	qsopt_ex/eg_lpnum.mpf.h \
	qsopt_ex/eg_timer.h \
	qsopt_ex/basicdefs.h \
	qsopt_ex/stddefs.h \
	qsopt_ex/trace.h \
	qsopt_ex/sortrus_common.h \
	qsopt_ex/QSopt_ex_version.h \
	qsopt_ex/QSopt_ex.h

PRIVATE_HEADER_FILES = \
	qsopt_ex/qs_config.h

# These source files have to be instantiated using the
# type template mechanism
TEMPLATE_SOURCES = \
	qsopt_ex/rawlp.c \
	qsopt_ex/mps.c \
	qsopt_ex/read_mps.c \
	qsopt_ex/lp.c \
	qsopt_ex/write_lp.c \
	qsopt_ex/read_lp.c \
	qsopt_ex/readline.c \
	qsopt_ex/lpdata.c \
	qsopt_ex/presolve.c \
	qsopt_ex/factor.c \
	qsopt_ex/basis.c \
	qsopt_ex/price.c \
	qsopt_ex/dstruct.c \
	qsopt_ex/simplex.c \
	qsopt_ex/fct.c \
	qsopt_ex/ratio.c \
	qsopt_ex/lib.c \
	qsopt_ex/binary.c \
	qsopt_ex/qsopt.c \
	qsopt_ex/sortrus.c \
	qsopt_ex/dheaps_i.c \
	qsopt_ex/priority.c \
	qsopt_ex/editor.c \
	qsopt_ex/format.c \
	qsopt_ex/util.c \
	qsopt_ex/eg_numutil.c

TEMPLATE_HEADER_FILES = \
	qsopt_ex/rawlp.h \
	qsopt_ex/mps.h \
	qsopt_ex/read_mps.h \
	qsopt_ex/lp.h \
	qsopt_ex/write_lp.h \
	qsopt_ex/read_lp.h \
	qsopt_ex/readline.h \
	qsopt_ex/lpdata.h \
	qsopt_ex/presolve.h \
	qsopt_ex/factor.h \
	qsopt_ex/basis.h \
	qsopt_ex/price.h \
	qsopt_ex/dstruct.h \
	qsopt_ex/simplex.h \
	qsopt_ex/fct.h \
	qsopt_ex/ratio.h \
	qsopt_ex/lib.h \
	qsopt_ex/binary.h \
	qsopt_ex/qsopt.h \
	qsopt_ex/sortrus.h \
	qsopt_ex/dheaps_i.h \
	qsopt_ex/priority.h \
	qsopt_ex/editor.h \
	qsopt_ex/format.h \
	qsopt_ex/lpdefs.h \
	qsopt_ex/qstruct.h \
	qsopt_ex/iqsutil.h \
	qsopt_ex/util.h \
	qsopt_ex/eg_numutil.h

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


AM_CFLAGS = -I$(top_srcdir)/qsopt_ex -I$(top_builddir)/qsopt_ex

# Library files
lib_LTLIBRARIES = libqsopt_ex.la
libqsopt_ex_la_SOURCES = \
	$(MAIN_SOURCE_FILES) $(MAIN_HEADER_FILES) \
	$(PRIVATE_HEADER_FILES)
nodist_libqsopt_ex_la_SOURCES = $(TEMPLATE_FILES)
libqsopt_ex_la_CFLAGS = \
	$(AM_CFLAGS) $(GMP_CFLAGS)
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
bin_PROGRAMS = esolver/esolver
noinst_PROGRAMS = tests/demo_qs tests/eg_sloan

tests_demo_qs_SOURCES = tests/demo_qs.c
tests_demo_qs_CFLAGS = $(AM_CFLAGS) $(GMP_CFLAGS)
tests_demo_qs_LDADD = $(GMP_LIBS) libqsopt_ex.la

tests_eg_sloan_SOURCES = tests/eg_sloan.c tests/eg_sloan.h
tests_eg_sloan_CFLAGS = $(AM_CFLAGS) $(GMP_CFLAGS)
tests_eg_sloan_LDADD = $(GMP_LIBS) libqsopt_ex.la

esolver_esolver_SOURCES = esolver/esolver.c
esolver_esolver_CFLAGS = $(AM_CFLAGS) $(GMP_CFLAGS)
esolver_esolver_LDADD = $(GMP_LIBS) libqsopt_ex.la


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
	README.md NEWS.md Doxyfile License.txt \
	$(TEMPLATE_SOURCES) $(TEMPLATE_HEADER_FILES) typenames \
	$(MODULE_PYX_FILES) $(MODULE_PXD_FILES) $(MODULE_TEST_FILES)

# Clean files
CLEANFILES = \
	$(TEMPLATE_FILES) $(TEMPLATE_SED_FILES) templatetags \
	$(MODULE_C_FILES)
