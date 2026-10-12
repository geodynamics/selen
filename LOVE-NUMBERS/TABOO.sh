#!/bin/sh
# 
# =================================================
#
# This sequence of commands produces the full set of Love numbers 
# and their combinations that are required for the functionality 
# of SELEN4. 
# 
# GS December 2019 
# DM March 2021
#
# =================================================
#
 CFG_FILE=$1
#
 if [ -z "$CFG_FILE" ]; then
     echo
     echo "ERROR: a configuration file name must be specified."
     echo
     exit 1
 fi
#
 taboo_path=bin
#
# Store the CFG_FILE name 
 echo $CFG_FILE > cfg_name
#
echo 
#---------------------------------
echo " >->-> Executing TABOO..."
#---------------------------------
$taboo_path/taboo.exe  $CFG_FILE
#
# 
echo 
#---------------------------------
echo " >->-> Secular Love numbers ..."
#---------------------------------
$taboo_path/fluid-ln.exe   $CFG_FILE
#
#
echo 
#---------------------------------
echo " >->-> PMTF and RRFs (traditional theory)"
#---------------------------------
$taboo_path/pmtf.exe   $CFG_FILE
#
#
echo 
#---------------------------------
echo " >->-> Gamma functions (traditional theory)"
#---------------------------------
$taboo_path/guv.exe    $CFG_FILE
#
#
#
echo 
#---------------------------------
echo " >->-> PMTF and RRFs (revised theory)"
#---------------------------------
$taboo_path/pmtf-rev.exe   $CFG_FILE
#
#
echo 
#---------------------------------
 echo " >->-> Gamma functions (revised theory)"
#---------------------------------
$taboo_path/guv-rev.exe    $CFG_FILE
#
#
#---------------------------------
 echo ""
 echo " >->-> All done!"
 echo ""
 echo " ** You may want to archive the outputs of this run with the command"
 echo " **"
 echo " **    $ sh pproc.sh <folder_name>"
 echo " **"
 echo " ** and/or to clean up the working directory with the command"
 echo " **" 
 echo " **    $ sh purge.sh"
 echo ""
#---------------------------------
#
