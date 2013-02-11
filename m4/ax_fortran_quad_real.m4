# ===========================================================================
#
# SYNOPSIS
#
#   AX_FORTRAN_QUAD_REAL()
#
# DESCRIPTION
#
#   This macro checks if Fortran supports quad precision reals, either at
#   10 bit or 16 bit precision. If it does support quad precision, the macro
#   adds -DSUPPORTS_QUAD_10 or -DSUPPORTS_QUAD_16 to the Fortran flags.
#   If it doesn't support quad precision, a failure occurs.
#
# LICENSE
#
#   Copyright (c) 2013 Eric Heien <emheien@ucdavis.edu>

AC_DEFUN([AX_FORTRAN_QUAD_REAL],[
AC_MSG_CHECKING(whether Fortran supports quad precision reals)
AC_LINK_IFELSE([
				AC_LANG_PROGRAM([],
					[       [integer, parameter :: qp = 16]
       [REAL(qp) :: test = 1.0_qp]
       [test = log10 (test)]])],
				FCFLAGS="$FCFLAGS -DSUPPORTS_QUAD_16"; AC_MSG_RESULT(128 bit REAL(16)),
				AC_LINK_IFELSE([
								AC_LANG_PROGRAM([],
									[       [integer, parameter :: qp = 10]
				       [REAL(qp) :: test = 1.0_qp]
				       [test = log10 (test)]])],
								FCFLAGS="$FCFLAGS -DSUPPORTS_QUAD_10"; AC_MSG_RESULT(80 bit REAL(10)),
								AC_MSG_ERROR(no)))

])

