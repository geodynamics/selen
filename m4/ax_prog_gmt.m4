# ===========================================================================
#
# SYNOPSIS
#
#   AX_PROG_GMT([ACTION-IF-NOT-FOUND])
#
# DESCRIPTION
#
#   This macro searches for the GMT mapping tools software and sets
#   the variable "gmt" to the name of the application or to "no" if
#   not found. If ACTION-IF-NOT-FOUND is not specified, configure will
#   fail when the program is not found.
#
#   Example:
#
#     AX_PROG_GMT()
#     AX_PROG_GMT([gmt_avail="no"])
#
# LICENSE
#
#   Copyright (c) 2013 Eric Heien <emheien@ucdavis.edu>

AC_DEFUN([AX_PROG_GMT],[
AC_CHECK_PROGS(gmt,[gmt GMT],no)
if test "x$gmt" = "xno" ;
then
	ifelse($#,0,[AC_MSG_WARN([GMT not found, all GMT scripts will fail])],
        $1)
fi
])

