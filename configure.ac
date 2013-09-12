# Process this file with autoconf to produce a configure script.

AC_INIT([SELEN], [2.9.11], BUG-REPORT-ADDRESS)
AC_CONFIG_MACRO_DIR([m4])

# Checks for programs.
AC_LANG(Fortran)
AC_PROG_FC

# Check for OpenMP
AC_OPENMP

# Check for Fortran specific flags, etc
AC_FC_MODULE_OUTPUT_FLAG

# Compiler must support at least 132 character columns
AC_FC_LINE_LENGTH(132)

# Compiler must support 16-bit or 10-bit quad precision
AX_FORTRAN_QUAD_REAL

# Run script to guess compiler flags
export FC
flags_guess="$SHELL flags.guess"
AC_MSG_NOTICE([running $flags_guess])
flags=`$flags_guess` ||
    AC_MSG_ERROR([$flags_guess failed])
eval $flags
AC_ARG_VAR(FCFLAGS_PLATFORM, [Fortran platform specific compiler flags])
AC_ARG_VAR(LDFLAGS_PLATFORM, [Fortran platform specific linker flags])
AC_ARG_VAR(RUN_COMMANDS, [Fortran runtime commands])

# Check for support programs
AX_PROG_GMT
ACLTX_PROG_PS2PDF([echo "WARNING: GMT scripts will fail"])

AC_CONFIG_FILES([Makefile])

AC_OUTPUT

