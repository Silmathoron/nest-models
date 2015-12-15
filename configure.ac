AC_PREREQ(2.52)
AC_INIT([nngt_module], [0.1], [nest_user@nest-initiative.org])

# Exporting source and build directories requires full path names.
# Thus we have to expand.
# Here, we are in top build dir, since source dir must exist, we can just
# move there and call pwd
if test "x$srcdir" = x ; then
  PKGSRCDIR=`pwd`
else
  PKGSRCDIR=`cd $srcdir && pwd`
fi
PKGBUILDDIR=`pwd`

# If this is not called, install-sh will be put into .. by bootstrap.sh
# moritz, 06-26-06
AC_CONFIG_AUX_DIR(.)

AM_INIT_AUTOMAKE([tar-ustar])

# obtain host system type; HEP 2004-12-20
AC_CANONICAL_HOST

# This is necessary to use AC_CHECK_HEADER
AC_USE_SYSTEM_EXTENSIONS

# ------------------------------------------------------------------------
# Handle options
#
# NOTE: No programs/compilations must be run in this section;
#       otherwise CFLAGS and CXXFLAGS may take on funny default
#       values.
#       HEP 2004-12-20
# ------------------------------------------------------------------------

# nest-config
NEST_CONFIG=`which nest-config`
AC_ARG_WITH(nest,[  --with-nest=script  nest-config script including path],
[
  if test "$withval" != yes; then
    NEST_CONFIG=$withval
  else
    AC_MSG_ERROR([--with-nest-config expects the nest-config script as argument. See README for details.])
  fi
])

# ---------------------------------------------------------------
# Do not build dynamic libraries for user models if on BlueGene
# ---------------------------------------------------------------

configure_bluegene=`$NEST_CONFIG --configure-bluegene`
AM_CONDITIONAL(BUILD_DYNAMIC_USER_MODULES, test x$configure_bluegene = xno)

if test x$configure_bluegene = xyes ; then
  AC_DISABLE_SHARED
  AC_ENABLE_STATIC
else
  AC_DISABLE_STATIC
  AC_ENABLE_SHARED
fi

# -------------------------------------------
# END Handle options
# -------------------------------------------

# does nest-config work
AC_MSG_CHECKING([for nest-config ])
AC_CHECK_FILE($NEST_CONFIG, HAVE_NEST=yes,
              AC_MSG_ERROR([No usable nest-config was found. You may want to use --with-nest-config.]))
AC_MSG_RESULT(found)

# the following will crash if nest-config does not run
# careful, lines below must not break
AC_MSG_CHECKING([for NEST directory information ])
NEST_PREFIX=`$NEST_CONFIG --prefix`
NEST_CPPFLAGS=`$NEST_CONFIG --cflags`
NEST_COMPILER=`$NEST_CONFIG --compiler`
if test $prefix = NONE; then prefix=`$NEST_CONFIG --prefix`; fi
AC_MSG_RESULT($NEST_CPPFLAGS)

# Set the platform-dependent compiler flags based on the canonical
# host string.  These flags are placed in AM_{C,CXX}FLAGS.  If
# {C,CXX}FLAGS are given as environment variables, then they are
# appended to the set of automatically chosen flags.  After
# {C,CXX}FLAGS have been read out, they must be cleared, since
# system-dependent defaults will otherwise be placed into the
# Makefiles.  HEP 2004-12-20.

# Before we can determine the proper compiler flags, we must know
# which compiler we are using.  Since the pertaining AC macros run the
# compiler and set CFLAGS, CXXFLAGS to system-dependent values, we
# need to save command line/enviroment settings of these variables
# first. AC_AIX must run before the compiler is run, so we must run it
# here.
# HEP 2004-12-21

MYMODULE_SAVE_CXXFLAGS=$CXXFLAGS

# Must first check if we are on AIX
AC_AIX

# Check for C++ compiler, looking for the same compiler
# used with NEST
AC_PROG_CXX([ $NEST_COMPILER ])

# the following is makeshift, should have the macro set proper
# MYMODULE_SET_CXXFLAGS
AM_CXXFLAGS=$MYMODULE_SAVE_CXXFLAGS
CXXFLAGS=

## Configure C environment

AC_PROG_LD
AC_PROG_INSTALL

AC_LIBLTDL_CONVENIENCE     ## put libltdl into a convenience library
AC_PROG_LIBTOOL      ## use libtool
AC_CONFIG_SUBDIRS(libltdl) ## also configure subdir containing libltdl


#-- Stuff for GSL -----------------------------------------
# NOTE: MUST come BEFORE switching to AC_LANG_CPLUSPLUS !!!
# HEP 2001-08-07
# NOTE: Adding @GSL_CFLAGS@ to the AM_CPPFLAGS variable
#       in Makefile.am will result in warnings about re-ordering
#       the search path.
# HEP 2002-07-10
# GSL Related HAVE_* variables:
#
#   HAVE_GSL     indicates that some GSL release version is available
#                (i.e., version 1.0 or higher)
#   HAVE_GSL_X_Y indicates that the GSL is version >= X.Y
#
# Currently implemented:
#   HAVE_GSL
#   HAVE_GSL_1_2   # minimum requirement for GSL Random Generators
#   HAVE_GSL_1_11  # fixes bug in ODE solvers, see #313 and #88
#
# Add further HAVE_GSL_MAJOR_MINOR if necessary according to scheme for
# HAVE_GSL_1_2 below.
#
# For each HAVE_GSL_X_Y, define also an automake conditional
# GSL_1_2_AVAILABLE to allow for conditional compilation.  See example
# below, and librandom for an example of how to use it in Makefile.am.
#
# HEP 2002-08-21
AM_PATH_GSL(1.0, ac_have_gsl="yes", ac_have_gsl="no")

if test "$ac_have_gsl" != no ; then
  AC_DEFINE(HAVE_GSL,1, [Is the GNU Science Library available (ver. >= 1.0)?])
  AC_SUBST(HAVE_GSL)
  AC_SUBST(GSL_LIBS)
  AC_SUBST(GSL_CFLAGS)
  AM_PATH_GSL(1.2, ac_have_gsl_1_2="yes", ac_have_gsl_1_2="no")
  if test "x$ac_have_gsl_1_2" = xyes ; then
    AC_DEFINE(HAVE_GSL_1_2, 1, [Is GSL Version >= 1.2?])

    AM_PATH_GSL(1.11, ac_have_gsl_1_11="yes", ac_have_gsl_1_11="no")
    if test "x$ac_have_gsl_1_11" = xyes ; then
      AC_DEFINE(HAVE_GSL_1_11, 1, [Is GSL Version >= 1.11?])
    else
      ac_have_gsl_1_11="no"
    fi
  fi
else
  unset GSL_LIBS
  unset GSL_CFLAGS
  unset GSL_VERSION
  ac_have_gsl_1_2="no"
fi
# AM_CONDITIONAL must not be invoked conditionally
AM_CONDITIONAL(GSL_1_2_AVAILABLE, test x$ac_have_gsl_1_2 = xyes)
AM_CONDITIONAL(GSL_1_11_AVAILABLE, test x$ac_have_gsl_1_11 = xyes)

# C Inline stuff revised, HEP 2003-03-11
SLI_C_INLINE

if test "$ac_cv_c_inline" != no ; then
  AC_DEFINE(HAVE_INLINE,1, [Do we have the inline macro in C?])
  AC_SUBST(HAVE_INLINE)
fi

## --- end GSL STUFF -----------------------------------

#-- Set the language to C++
AC_LANG_CPLUSPLUS

#-- Look for programs needed in the Makefile
AC_PROG_CXXCPP
AM_PROG_LIBTOOL
AC_PATH_PROGS([MAKE],[gmake make],[make])


# ---------------------------------------------------------------
# Configure directories to be built
# ---------------------------------------------------------------

PKGDATADIR=$datadir/$PACKAGE
PKGDOCDIR=$datadir/doc/$PACKAGE

# set up directories from which to build help
# second line replaces space with colon as separator
HELPDIRS="$PKGSRCDIR $PKGSRCDIR/sli"
HELPDIRS=`echo $HELPDIRS | tr " " ":"`

#-- Replace these variables in *.in
AC_SUBST(HAVE_NEST)
AC_SUBST(NEST_CONFIG)
AC_SUBST(NEST_CPPFLAGS)
AC_SUBST(NEST_COMPILER)
AC_SUBST(NEST_PREFIX)
AC_SUBST(HELPDIRS)
AC_SUBST(PKGSRCDIR)
AC_SUBST(PKGBUILDDIR)
AC_SUBST(PKGDATADIR)
AC_SUBST(PKGDOCDIR)
AC_SUBST(KERNEL)
AC_SUBST(HOST)
AC_SUBST(SED)
AC_SUBST(LD)
AC_SUBST(host_os)
AC_SUBST(host_cpu)
AC_SUBST(host_vendor)
AC_SUBST(AS)
AC_SUBST(CXX)
AC_SUBST(AR)
AC_SUBST(ARFLAGS)
AC_SUBST(CXX_AR)
AC_SUBST(AM_CXXFLAGS)
AC_SUBST(AM_CFLAGS)
AC_SUBST(MAKE)
AC_SUBST(MAKE_FLAGS)
AC_SUBST(INCLTDL)
AC_SUBST(LIBLTDL)

AC_CONFIG_HEADER(nngt_module_config.h:nngt_module_config.h.in)
AC_CONFIG_FILES(Makefile)

# -----------------------------------------------
# Create output
# -----------------------------------------------
AC_OUTPUT


# -----------------------------------------------
# Report, after output at end of configure run
# Must come after AC_OUTPUT, so that it is
# displayed after libltdl has been configured
# -----------------------------------------------

echo
echo "-------------------------------------------------------"
echo "NngtModule Configuration Summary"
echo "-------------------------------------------------------"
echo
echo "C++ compiler        : $CXX"
echo "C++ compiler flags  : $AM_CXXFLAGS"
echo "NEST compiler flags : $NEST_CPPFLAGS"

# these variables will still contain '${prefix}'
# we want to have the versions where this is resolved, too:
eval eval eval PKGDOCDIR_AS_CONFIGURED=$PKGDOCDIR
eval eval eval PKGDATADIR_AS_CONFIGURED=$PKGDATADIR
eval eval eval LIBDIR_AS_CONFIGURED=$libdir
echo
echo "-------------------------------------------------------"
echo
echo "You can build and install NngtModule now, using"
echo "  make"
echo "  make install"
echo
echo "NngtModule will be installed to: $LIBDIR_AS_CONFIGURED"
echo
