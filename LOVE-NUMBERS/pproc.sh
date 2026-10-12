#!/bin/sh
# 
# =================================================
#
# This sequence of commands constitutes the post-processing phase
# 
# GS December 2019 
# DM March 2021
#
# =================================================
#
#
#
 LABEL=$1
#
 if [ -z "$LABEL" ]; then
     echo
     echo "Usage: pproc.sh <folder_name>"
     echo
     exit 1
 fi

#
echo
echo "POST-PROCESSING..."
#
#LABEL="RUN-NV20-N52-M07" 
#
echo
echo "the label chosen is " $LABEL
echo
if  [ -d "$LABEL" ] 
then
echo WARNING: data exist in $LABEL - they will be overwritten  
echo 
else 
echo Folder $LABEL does not exist - a new one is created 
mkdir $LABEL
fi
# -------------------------------------------------
echo Copying the backup model file   
#
if  [ -f "run-model-parameters.dat" ] 
then 
mv run-model-parameters.dat "./"$LABEL"/run-model-parameters.dat"
fi
# -------------------------------------------------
echo Copying the load Love numbers  
#
if  [ -f "h-load.dat" ] 
then 
mv h-load.dat "./"$LABEL"/h-load.dat"
fi
if  [ -f "k-load.dat" ] 
then 
mv k-load.dat "./"$LABEL"/k-load.dat"
fi
if  [ -f "l-load.dat" ] 
then 
mv l-load.dat "./"$LABEL"/l-load.dat"
fi
# -------------------------------------------------
echo Copying the tidal Love numbers  
#
if  [ -f "h-tide.dat" ] 
then 
mv h-tide.dat "./"$LABEL"/h-tide.dat"
fi
if  [ -f "k-tide.dat" ] 
then 
mv k-tide.dat "./"$LABEL"/k-tide.dat"
fi
if  [ -f "l-tide.dat" ] 
then 
mv l-tide.dat "./"$LABEL"/l-tide.dat"
fi
# -------------------------------------------------
echo Copying the spectrum of relaxation
#
if  [ -f "spectrum.dat" ] 
then 
mv spectrum.dat "./"$LABEL"/spectrum.dat"
fi
# -------------------------------------------------
echo Copying the degree 2 data   
#
if  [ -f "l2-data.dat" ] 
then 
mv l2-data.dat "./"$LABEL"/l2-data.dat"
fi
# -------------------------------------------------
echo Copying the E and Beta functions files  
#
if  [ -f "ebg.dat" ] 
then 
mv ebg.dat "./"$LABEL"/ebg.dat"
fi
if  [ -f "ebu.dat" ] 
then 
mv ebu.dat "./"$LABEL"/ebu.dat"
fi
if  [ -f "ebv.dat" ] 
then 
mv ebv.dat "./"$LABEL"/ebv.dat"
fi
if  [ -f "ebs.dat" ] 
then 
mv ebs.dat "./"$LABEL"/ebs.dat"
fi
# -------------------------------------------------
echo Copying the PMTF files  
#
if  [ -f "PMTF.dat" ] 
then 
mv PMTF.dat "./"$LABEL"/PMTF.dat"
fi
if  [ -f "PMTF-REV.dat" ] 
then 
mv PMTF-REV.dat "./"$LABEL"/PMTF-REV.dat"
fi
# -------------------------------------------------
echo Copying the Gamma functions files  
#
if  [ -f "gamma.dat" ] 
then 
mv gamma.dat "./"$LABEL"/gamma.dat"
fi
if  [ -f "gamma-REV.dat" ] 
then 
mv gamma-REV.dat "./"$LABEL"/gamma-REV.dat"
fi
# -------------------------------------------------
echo Copying miscellanea configuration and log files  
#
echo TABOO configuration file
if  [ -f "cfg_name" ] 
then 
cp `cat cfg_name` "./"$LABEL"/taboo.cfg"
rm cfg_name
fi
echo taboo.log
if  [ -f "taboo.log" ] 
then 
mv taboo.log "./"$LABEL"/taboo.log"
fi
echo the secular fluid k LN file
if  [ -f "./fluid-ln.dat" ] 
then 
mv fluid-ln.dat "./"$LABEL"/fluid-ln.dat"
fi
echo fluid-ln.log
if  [ -f "fluid-ln.log" ] 
then 
mv fluid-ln.log "./"$LABEL"/fluid-ln.log"
fi
echo Good bye and good luck 
echo ""






