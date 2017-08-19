#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([vcftools],
        [m4_esyscmd([build-aux/git-version-gen .tarball-version])],
        [https://github.com/vcftools/vcftools/issues])
AC_CONFIG_SRCDIR([src/cpp/vcftools.cpp])
AC_CONFIG_HEADERS([config.h])

# Automake invocation.
AM_INIT_AUTOMAKE([foreign])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_CPP

# Checks for perl.
AC_PATH_PROGS([PERL], [perl] [perl5], [false])
AS_IF([test "x$PERL" = "xfalse"],[
  AC_MSG_ERROR([Perl not found; check your \$PATH.])
])

pmdir_relative_path=`\
$PERL -MConfig \
      -wle '($_ = $Config{installsitelib})
            =~ s!^\Q$Config{siteprefix}/!!; \
            print'`

AC_ARG_WITH(
    [pmdir],
    AS_HELP_STRING(
      [--with-pmdir=DIR],
      [install Perl modules in DIR]),
    [PMDIR=${withval}],
    [PMDIR="$pmdir_relative_path"])

AC_SUBST([PMDIR])

# Checks for libraries.
PKG_CHECK_MODULES([ZLIB], [zlib])

# Checks for header files.
AC_CHECK_HEADERS([arpa/inet.h fcntl.h limits.h netdb.h stdint.h stdlib.h string.h sys/socket.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_INT8_T
AC_TYPE_OFF_T
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT8_T

# Checks for operating system services or capabilities.
AC_SYS_LARGEFILE

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([gethostbyaddr gethostbyname memset pow select socket sqrt strchr strdup strerror strstr strtol])

# Optional features.
AC_ARG_ENABLE([pca],
              AS_HELP_STRING([--enable-pca], [enable PCA feature]),
              [pca=${enableval}],
              [pca=no])

AS_IF([test "x$pca" = "xyes"],[
  AC_CHECK_LIB([lapack], [dgeev_])
])

# Generate output.
AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/cpp/Makefile
                 src/perl/Makefile])
AC_OUTPUT
