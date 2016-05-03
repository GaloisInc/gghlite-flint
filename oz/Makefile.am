AUTOMAKE_OPTIONS = foreign -Wall

AM_CFLAGS  = ${DEBUG_CFLAGS} -I$(top_srcdir) -fopenmp

lib_LTLIBRARIES=liboz.la

liboz_la_SOURCES = oz.c flint-addons.c util.c sqrt.c invert.c mul.c ntt.c norm.c rem.c
liboz_la_LDFLAGS = -version-info $(OZ_VERSION_INFO) -no-undefined
liboz_la_INCLUDEDIR = $(includedir)/oz
liboz_la_LIBADD = -lgomp

pkgincludesubdir = $(includedir)/oz
pkgincludesub_HEADERS = oz.h flags.h flint-addons.h sqrt.h invert.h mul.h \
	norm.h rem.h ntt.h
noinst_HEADERS = util.h
