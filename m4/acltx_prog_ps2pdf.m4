AC_DEFUN([ACLTX_PROG_PS2PDF],[
AC_ARG_VAR(ps2pdf,[specify default ps2pdf application])
if test "$ac_cv_env_ps2pdf_set" = "set" ; then
    AC_MSG_CHECKING([Checking for ps2pdf])
    ps2pdf="$ac_cv_env_ps2pdf_value";
    AC_MSG_RESULT([$ps2pdf (from parameter)])
else
    AC_CHECK_PROGS(ps2pdf,[ps2pdf14 ps2pdf13 ps2pdf12 ps2pdf pstopdf],no)
fi
if test $ps2pdf = "no" ;
then
        ifelse($#,0,[AC_MSG_ERROR([Unable to find the ps2pdf application])],
        $1)
fi
])

