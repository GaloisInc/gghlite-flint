AUTOMAKE_OPTIONS = foreign
AM_CFLAGS = ${DEBUG_CFLAGS} -I$(top_srcdir) -I$(top_srcdir)/dgs \
            -std=c99 -D_DEFAULT_SOURCE -fopenmp

lib_LTLIBRARIES = libgghlite.la

pkgincludesubdir = $(includedir)/gghlite

include_HEADERS = gghlite.h
pkgincludesub_HEADERS = config.h gghlite-defs.h \
                        gghlite-internals.h \
                        misc.h

libgghlite_la_SOURCES = gghlite.c \
                        gghlite_pk.c \
                        misc.c \
                        lattice_reduction.c \
                        gghlite.h \
                        misc.h \
                        ggh-defs.h \
                        ggh-internals.h \
                        api.c
libgghlite_la_LIBADD = $(top_builddir)/oz/liboz.la \
                       $(top_builddir)/dgs/libdgs.la \
                       $(top_builddir)/dgsl/libdgsl.la

libgghlite_la_LDFLAGS = -version-info $(GGHLITE_VERSION_INFO) -no-undefined
libgghlite_la_INCLUDEDIR = $(includedir)/gghlite
