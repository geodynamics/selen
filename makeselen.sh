# 
# This is file "makeselen.sh" 
#
# Last modified GS 04-07-2008 [INTEL port]
# Revised GS & FC July 2009 - "Varying coastlines"
# Revised GS & FC May 2010 - v 2.9 for g95
# Revised DM Apr 2011 - gfortran/Linux compatibility
# Revised DM July 2011 - config.dat parsing for auto-compiler setup
#
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
# Copyright (C) 2008 Giorgio Spada, Florence Colleoni, and Paolo Stocchi 
#
# This file is part of SELEN. 
#  
# SELEN is free software: you can redistribute it and/or modify it under the 
# terms of the GNU General Public License as published by the Free Software 
# Foundation, either version 3 of the License, or at your option) any later 
# version. 
#
# SELEN is distributed in the hope that it will be useful, but WITHOUT ANY 
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
# details. 
# 
# You should have received a copy of the GNU General Public License along 
# with SELEN.  If not, see <http://www.gnu.org/licenses/>.
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#
#
# ----------------------------------------------------------------
#  This self-explicatory script controls the execution os SELEN- 
#  GS & FC 12-9-2007. Help available from giorgio.spada@gmail.com 
# ----------------------------------------------------------------
#
#
#
# --- Removes the stop file 'stop.dat' from a possible previous crash
if [ -f stop.dat ] 
then 
echo 
echo --- Recovering from previous crash
/bin/rm stop.dat 
fi 
#
echo
date
#
# --- A small welcome message 
echo
echo ========= This is SELEN 2.9.10 ========= 
#
# --- Parse the config.dat file to identify the fortran compiler
SYS=`grep 121 config.dat | awk -F\' '{print $2}'`
case "$SYS" in
   '1' | '3' | '7' )
     FC='gfortran -w'
     ;;
   '2' | '4' | '5' )
     FC='ifort'
     export FORT_FMT_RECL=500
     ;;
   '6' )
     FC='xlf90'
     ;;
   '8' )
     FC='g95'
     ;;     
esac
echo
echo --- Compile command is $FC
#
# --- Compiles the SHTOOLS configuration file 
echo
echo --- Compiling shtools.f90
$FC -c shtools.f90
#  
#
# --- Compiles and executes 'config.f90' according to SELEN settings in 'config.dat'
#
echo --- Compiling config.f90
$FC config.f90 -o config.exe   
#
echo --- Executing config.f90
./config.exe 
#
#
# --- If no error conditions are met, selen.sh is executed.      
if [ ! -f stop.dat ] 
then 
echo --- No errors detected 
echo --- SELEN is launched  
time sh selen.sh
# --- Purging before execution
echo
date
else
echo --- Errors have been detected 
echo --- SELEN is NOT launched 
echo --- For help, please check the log files 
fi 
