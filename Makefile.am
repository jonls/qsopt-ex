
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
	src/eg_eheap.h \
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

$(TEMPLATE_SOURCES_DBL): $(TEMPLATE_SOURCES) Makefile template_dbl.awk
	$(AM_V_GEN)$(MKDIR_P) $(@D) && \
		$(AWK) -f template_dbl.awk $(top_srcdir)/$(@:_dbl.c=.c) > $@
$(TEMPLATE_HEADERS_DBL): $(TEMPLATE_HEADER_FILES) Makefile template_dbl.awk
	$(AM_V_GEN)$(MKDIR_P) $(@D) && \
		$(AWK) -f template_dbl.awk $(top_srcdir)/$(@:_dbl.h=.h) > $@

TEMPLATE_SOURCES_FP20 = $(TEMPLATE_SOURCES:.c=_fp20.c)
TEMPLATE_HEADERS_FP20 = $(TEMPLATE_HEADER_FILES:.h=_fp20.h)

$(TEMPLATE_SOURCES_FP20): $(TEMPLATE_SOURCES) Makefile template_fp20.awk
	$(AM_V_GEN)$(MKDIR_P) $(@D) && \
		$(AWK) -f template_fp20.awk $(top_srcdir)/$(@:_fp20.c=.c) > $@
$(TEMPLATE_HEADERS_FP20): $(TEMPLATE_HEADER_FILES) Makefile template_fp20.awk
	$(AM_V_GEN)$(MKDIR_P) $(@D) && \
		$(AWK) -f template_fp20.awk $(top_srcdir)/$(@:_fp20.h=.h) > $@

TEMPLATE_SOURCES_MPQ = $(TEMPLATE_SOURCES:.c=_mpq.c)
TEMPLATE_HEADERS_MPQ = $(TEMPLATE_HEADER_FILES:.h=_mpq.h)

$(TEMPLATE_SOURCES_MPQ): $(TEMPLATE_SOURCES) Makefile template_mpq.awk
	$(AM_V_GEN)$(MKDIR_P) $(@D) && \
		$(AWK) -f template_mpq.awk $(top_srcdir)/$(@:_mpq.c=.c) > $@
$(TEMPLATE_HEADERS_MPQ): $(TEMPLATE_HEADER_FILES) Makefile template_mpq.awk
	$(AM_V_GEN)$(MKDIR_P) $(@D) && \
		$(AWK) -f template_mpq.awk $(top_srcdir)/$(@:_mpq.h=.h) > $@

TEMPLATE_SOURCES_MPF = $(TEMPLATE_SOURCES:.c=_mpf.c)
TEMPLATE_HEADERS_MPF = $(TEMPLATE_HEADER_FILES:.h=_mpf.h)

$(TEMPLATE_SOURCES_MPF): $(TEMPLATE_SOURCES) Makefile template_mpf.awk
	$(AM_V_GEN)$(MKDIR_P) $(@D) && \
		$(AWK) -f template_mpf.awk $(top_srcdir)/$(@:_mpf.c=.c) > $@
$(TEMPLATE_HEADERS_MPF): $(TEMPLATE_HEADER_FILES) Makefile template_mpf.awk
	$(AM_V_GEN)$(MKDIR_P) $(@D) && \
		$(AWK) -f template_mpf.awk $(top_srcdir)/$(@:_mpf.h=.h) > $@

TEMPLATE_SOURCES_FLOAT128 = $(TEMPLATE_SOURCES:.c=_float128.c)
TEMPLATE_HEADERS_FLOAT128 = $(TEMPLATE_HEADER_FILES:.h=_float128.h)

$(TEMPLATE_SOURCES_FLOAT128): $(TEMPLATE_SOURCES) Makefile template_float128.awk
	$(AM_V_GEN)$(MKDIR_P) $(@D) && \
		$(AWK) -f template_float128.awk $(top_srcdir)/$(@:_float128.c=.c) > $@
$(TEMPLATE_HEADERS_FLOAT128): $(TEMPLATE_HEADER_FILES) Makefile template_float128.awk
	$(AM_V_GEN)$(MKDIR_P) $(@D) && \
		$(AWK) -f template_float128.awk $(top_srcdir)/$(@:_float128.h=.h) > $@

TEMPLATE_SOURCES_LDBL = $(TEMPLATE_SOURCES:.c=_ldbl.c)
TEMPLATE_HEADERS_LDBL = $(TEMPLATE_HEADER_FILES:.h=_ldbl.h)

$(TEMPLATE_SOURCES_LDBL): $(TEMPLATE_SOURCES) Makefile template_ldbl.awk
	$(AM_V_GEN)$(MKDIR_P) $(@D) && \
		$(AWK) -f template_ldbl.awk $(top_srcdir)/$(@:_ldbl.c=.c) > $@
$(TEMPLATE_HEADERS_LDBL): $(TEMPLATE_HEADER_FILES) Makefile template_ldbl.awk
	$(AM_V_GEN)$(MKDIR_P) $(@D) && \
		$(AWK) -f template_ldbl.awk $(top_srcdir)/$(@:_ldbl.h=.h) > $@

# Define complete set of sources
PROGRAM_SOURCES = \
	$(MAIN_SOURCES) \
	$(TEMPLATE_SOURCES_DBL) $(TEMPLATE_HEADERS_DBL) \
	$(TEMPLATE_SOURCES_FP20) $(TEMPLATE_HEADERS_FP20)

if ENABLE_LIBGMP
PROGRAM_SOURCES += \
	$(TEMPLATE_SOURCES_MPQ) $(TEMPLATE_HEADERS_MPQ) \
	$(TEMPLATE_SOURCES_MPF) $(TEMPLATE_HEADERS_MPF)
endif
if ENABLE_SOFTFLOAT
PROGRAM_SOURCES += \
	$(TEMPLATE_SOURCES_FLOAT128) $(TEMPLATE_HEADERS_FLOAT128)
endif
if ENABLE_LONG_DOUBLE
PROGRAM_SOURCES += \
	$(TEMPLATE_SOURCES_LDBL) $(TEMPLATE_HEADERS_LDBL)
endif

# Define complete set of template files
TEMPLATE_FILES = \
	$(TEMPLATE_SOURCES_DBL) $(TEMPLATE_HEADERS_DBL) \
	$(TEMPLATE_SOURCES_FP20) $(TEMPLATE_HEADERS_FP20) \
	$(TEMPLATE_SOURCES_MPQ) $(TEMPLATE_HEADERS_MPQ) \
	$(TEMPLATE_SOURCES_MPF) $(TEMPLATE_HEADERS_MPF) \
	$(TEMPLATE_SOURCES_FLOAT128) $(TEMPLATE_HEADERS_FLOAT128) \
	$(TEMPLATE_SOURCES_LDBL) $(TEMPLATE_HEADERS_LDBL)

template_%.awk: templatetags Makefile typenames
	$(AM_V_GEN)$(MKDIR_P) $(@D) && \
		type_prefix="$$(echo $(@) | $(SED) -e 's|^template_\(.\+\)\.awk$$|\1|')" && \
		type_name="$$(grep ^$${type_prefix} $(top_srcdir)/typenames | cut -f2)" && \
		echo "/\<EGlpNum_t\>/ {gsub(/\<EGlpNum_t\>/,\"$${type_name}\")};" > $@ && \
		for name in $$(cat templatetags); do \
			echo "/\<$${name}\>/ {gsub(/\<$${name}\>/,\"$${type_prefix}_$${name}\")};" \
		; done >> $@ ; \
		for header in $(TEMPLATE_HEADER_FILES); do \
			header_base="$$(basename $${header})" && \
			echo "/\<$${header_base}\>/ {gsub(/\<$${header_base}\>/,\"$${header_base%%.h}_$${type_prefix}.h\")};" \
		; done >> $@ ; \
		echo "{print};" >> $@

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
TEMPLATE_AWK_FILES = $(wildcard template_*.awk)

# Clean files
CLEANFILES = $(TEMPLATE_FILES) $(TEMPLATE_AWK_FILES) templatetags

# Program files
bin_PROGRAMS = demo_qs eg_sloan esolver
AM_CFLAGS = -I$(top_srcdir)/src -I$(top_builddir)/src
LDADD = $(LIBGMP)

demo_qs_SOURCES = $(PROGRAM_SOURCES) src/demo_qs.c
eg_sloan_SOURCES = $(PROGRAM_SOURCES) src/eg_sloan.c src/eg_sloan.h
esolver_SOURCES = $(PROGRAM_SOURCES) src/esolver.c
