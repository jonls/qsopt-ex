
ACLOCAL_AMFLAGS = -I m4

# These source files are used as is
MAIN_SOURCES = \
	src/allocrus.c src/allocrus.h \
	src/bgetopt.c src/bgetopt.h \
	src/eg_io.c src/eg_io.h \
	src/eg_lpnum.c src/eg_lpnum.h \
	src/except.c src/except.h \
	src/urandom.c src/urandom.h \
	src/zeit.c src/zeit.h \
	src/names.c src/names.h \
	src/symtab.c src/symtab.h \
	src/util.c src/util.h \
	src/exact.c src/exact.h \
	src/reporter.c src/reporter.h \
	src/eg_exutil.c src/eg_exutil.h \
	src/eg_macros.c src/eg_macros.h \
	src/eg_memslab.c src/eg_memslab.h \
	src/basicdefs.h \
	src/sortrus_common.c src/sortrus_common.h \
	src/QSopt_ex_version.c

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

# Define library sources
LIBRARY_SOURCES = \
	$(MAIN_SOURCES) \
	$(TEMPLATE_SOURCES_DBL) $(TEMPLATE_HEADERS_DBL) \
	$(TEMPLATE_SOURCES_MPQ) $(TEMPLATE_HEADERS_MPQ) \
	$(TEMPLATE_SOURCES_MPF) $(TEMPLATE_HEADERS_MPF)

# Define complete set of template files
TEMPLATE_FILES = \
	$(TEMPLATE_SOURCES_DBL) $(TEMPLATE_HEADERS_DBL) \
	$(TEMPLATE_SOURCES_MPQ) $(TEMPLATE_HEADERS_MPQ) \
	$(TEMPLATE_SOURCES_MPF) $(TEMPLATE_HEADERS_MPF)

template_%.sed: templatetags Makefile typenames
	$(AM_V_GEN)$(MKDIR_P) $(@D) && \
		type_prefix="$$(echo $(@) | $(SED) -e 's|^template_\(.\+\)\.sed$$|\1|')" && \
		type_name="$$(grep ^$${type_prefix} $(top_srcdir)/typenames | cut -f2)" && \
		echo "s|\<EGlpNum_t\>|$${type_name}|g" > $@ && \
		for name in $$(cat templatetags); do \
			echo "s|\<$${name}\>|$${type_prefix}_&|g" \
		; done >> $@ ; \
		for header in $(TEMPLATE_HEADER_FILES); do \
			header_base="$$(basename $${header})" && \
			echo "s|\<$${header_base}\>|$${header_base%%.h}_$${type_prefix}.h|g" \
		; done >> $@

templatetags: $(TEMPLATE_SOURCES) $(TEMPLATE_HEADER_FILES)
	$(AM_V_GEN)$(MKDIR_P) $(@D) && \
		( ( $(CTAGS) -x --c-kinds=+xp-m --c++-kinds=+xp-m $^ | \
		cut -d ' ' -f 1 ; \
		for name in $(TEMPLATE_NAMES); do echo "$${name}"; done ) | \
		sort | uniq ; \
		for name in $(TEMPLATE_NAMES_GLOBAL); do echo "$${name}"; echo "$${name}"; done ) | \
		sort | uniq -u > $@

# Tell automake that these are generated files
BUILT_SOURCES = $(TEMPLATE_FILES)

# Define tag files
TEMPLATE_SED_FILES = $(wildcard template_*.sed)

# Clean files
CLEANFILES = $(TEMPLATE_FILES) $(TEMPLATE_SED_FILES) templatetags


AM_CFLAGS = -I$(top_srcdir)/src -I$(top_builddir)/src

# Library files
lib_LTLIBRARIES = libqsopt_ex.la
libqsopt_ex_la_SOURCES = $(LIBRARY_SOURCES)
libqsopt_ex_la_LIBADD = $(LIBGMP)

# Program files
bin_PROGRAMS = demo_qs eg_sloan esolver

demo_qs_SOURCES = src/demo_qs.c
demo_qs_LDADD = libqsopt_ex.la

eg_sloan_SOURCES = src/eg_sloan.c src/eg_sloan.h
eg_sloan_LDADD = libqsopt_ex.la

esolver_SOURCES = src/esolver.c
esolver_LDADD = libqsopt_ex.la
