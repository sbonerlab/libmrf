EXTRA_DIST = m4
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}
AM_MAKEFLAGS = --no-print-directory

AM_CPPFLAGS = \
	-I${top_srcdir}/src

PC_SED = \
	$(AM_V_GEN)$(MKDIR_P) $(dir $@) && $(SED) \
	-e 's,@VERSION\@,$(VERSION),g' \
	-e 's,@prefix\@,$(prefix),g' \
	-e 's,@exec_prefix\@,$(exec_prefix),g' \
	-e 's,@libdir\@,$(libdir),g' \
	-e 's,@includedir\@,$(includedir),g' \
	< $< > $@ || rm $@

%.pc: %.pc.in Makefile
	$(PC_SED)

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = mrf.pc
EXTRA_DIST += mrf.pc

lib_LTLIBRARIES = libmrf.la
libmrf_la_SOURCES = \
	mrf/mrf.c \
	mrf/mrfUtil.c \
	mrf/sam.c \
	mrf/segmentationUtil.c

libmrf_la_LIBADD = -lbios
nobase_dist_include_HEADERS = \
	mrf/mrf.h \
    mrf/mrfUtil.h \
    mrf/sam.h \
    mrf/segmentationUtil.h
