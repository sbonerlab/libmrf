SUBDIRS = lib
ACLOCAL_AMFLAGS = -I m4
EXTRA_DIST = m4

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = mrf.pc

debug:
	$(MAKE) all "CFLAGS=-g -DDEBUG " $(AM_MAKEFLAGS)