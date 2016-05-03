AUTOMAKE_OPTIONS = foreign
AM_CFLAGS = ${DEBUG_CFLAGS} -I$(top_srcdir) -I$(top_srcdir)/dgs -std=c99 -fopenmp

lib_LTLIBRARIES = libdgsl.la

libdgsl_la_SOURCES = dgsl.c gso.c

pkgincludesubdir = $(includedir)/dgsl
pkgincludesub_HEADERS = dgsl.h gso.h

libdgsl_la_LDFLAGS = -version-info $(DGSL_VERSION_INFO) -no-undefined
libdgsl_la_INCLUDEDIR = $(includedir)/dgsl
libdgsl_la_LIBADD = -L$(top_builddir)/dgs -ldgs -L$(top_builddir)/oz -loz
