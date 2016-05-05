m4_define([version_current], [1])
m4_define([version_revision], [0])
m4_define([version_age], [0])

AC_INIT([gghflint], [version_current.version_revision.version_age])
AC_CONFIG_AUX_DIR([.])

AC_CONFIG_SRCDIR([gghlite])

AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([subdir-objects])

GGHLITE_VERSION_INFO=1:0:0
AC_SUBST(GGHLITE_VERSION_INFO)
OZ_VERSION_INFO=1:0:0
AC_SUBST(OZ_VERSION_INFO)
DGS_VERSION_INFO=1:0:0
AC_SUBST(DGS_VERSION_INFO)
DGSL_VERSION_INFO=1:0:0
AC_SUBST(DGSL_VERSION_INFO)

# support older automake versions
m4_ifdef([AM_SILENT_RULES],
         [AM_SILENT_RULES([yes])],
         [AC_SUBST([AM_DEFAULT_VERBOSITY], [1])])

AC_CONFIG_HEADERS([gghlite/config.h])

m4_pattern_allow([AM_PROG_AR])
AM_PROG_AR()

AC_PROG_LIBTOOL()

AC_PROG_CC_C99()

AS_IF([test "$ac_cv_prog_cc_c99" = "no"], AC_MSG_ERROR([C99 support is required but not found.]))

AC_ARG_ENABLE(debug,        [  --enable-debug          Enable assert() statements for debugging.], [enable_debug=yes])
AC_ARG_ENABLE(heuristics,   [  --enable-heuristics     Enable various heuristic improvements.], [enable_heuristics=yes])

CFLAGS=                         dnl get rid of default -g -O2
COMMON_CFLAGS="-Wall -Wformat -Wformat-security -Wextra -Wunused \
-Wshadow -Wmissing-prototypes -Wfloat-equal -Wpointer-arith -Wcast-align \
-Wstrict-prototypes -Wredundant-decls -Wendif-labels -Wcast-qual \
-std=gnu11 -Wpedantic"

if test "x$enable_debug" = x"yes"; then
  EXTRA_CFLAGS="-O0 -g"
else
  EXTRA_CFLAGS="-O3"
  AC_DEFINE(NDEBUG,1,[Define whether debugging is enabled])
fi
AC_SUBST(COMMON_CFLAGS)
AC_SUBST(EXTRA_CFLAGS)

AX_VALGRIND_CHECK()

# Heuristics
if test "x$enable_heuristics" = x"yes"; then
   AC_DEFINE(GGHLITE_HEURISTICS,1,[Define whether to enable performance improvements])
fi

AC_CHECK_HEADERS([omp.h])

AC_SEARCH_LIBS(mpz_urandomm_aes,aesrand)
if test "x$ac_cv_search_mpz_urandomm_aes" = "xno"; then
  AC_MSG_ERROR([libaesrand not found])
fi
AC_SEARCH_LIBS(mpfr_init2,mpfr)
if test "x$ac_cv_search_mpfr_init2" = "xno"; then
  AC_MSG_ERROR([libmpfr not found])
fi
AC_SEARCH_LIBS(fmpz_init,flint)
if test "x$ac_cv_search_fmpz_init" = "xno"; then
  AC_MSG_ERROR([libflint not found])
fi

AC_CONFIG_FILES([Makefile oz/Makefile dgs/Makefile dgsl/Makefile gghlite/Makefile applications/Makefile tests/Makefile])

AC_OUTPUT
