AUTOMAKE_OPTIONS = foreign
AM_CFLAGS=${DEBUG_CFLAGS} -std=c99

lib_LTLIBRARIES = libdgs.la

libdgs_la_SOURCES = dgs/dgs_bern.c \
                    dgs/dgs_gauss_dp.c \
                    dgs/dgs_gauss_mp.c	

pkgincludesubdir = $(includedir)/dgs
pkgincludesub_HEADERS = dgs/dgs_bern.h \
                        dgs/dgs_gauss.h \
                        dgs/dgs.h \
                        dgs/dgs_misc.h

libdgs_la_LDFLAGS = -version-info $(DGSL_VERSION_INFO) -no-undefined

check_PROGRAMS=test_gauss_z
test_gauss_z_SOURCES=tests/test_gauss_z.c
test_gauss_z_LDFLAGS= -no-install -lm -ldgs
test_gauss_z_CFLAGS=$(AM_CFLAGS)

TESTS = test_gauss_z

clean-local:
	(cd tests; make clean; cd ..)
