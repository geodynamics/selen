 echo
############################################################# 
 date
#
 SLE_SOLVER=./sle.exe
#
 PPR_PROGRAM=./ppr.exe
#
 RUN_NAME=$1
 CONFIG_NAME=$2
#
 DEPOT_NAME=RUN_$RUN_NAME
# 
# 
 echo
 echo ::::::: Executing the script: make_sle.sh 
 echo
 echo +---------------------------------------------+  
 echo +------ Solving the Sea Level Equation -------+
 echo +------------- This is SELEN 4 ---------------+
 echo +--- Copyright G. Spada and D. Melini 2019 ---+
 echo +---------------------------------------------+   
 echo 
#
 if [ -z "$RUN_NAME" ]; then
     echo
     echo "Usage: make_sle.sh <RUN_LABEL> <CONFIG_NAME>"
     echo "       if not explictly given, CONFIG_NAME=config.sle.RUN_LABEL"
     echo
     exit 1
 fi
#
 if [ -z "$CONFIG_NAME" ]; then
     CONFIG_NAME=config.sle.$RUN_NAME
 fi
#
 echo ======= The RUN label is: $RUN_NAME
 echo         Outputs will be stored in: $DEPOT_NAME
 echo
#
 echo ======= The config file is: $CONFIG_NAME
 echo
#
 echo ======= The SLE solver is: $SLE_SOLVER
 echo
 echo ======= The POST-PROCESSING program is: $PPR_PROGRAM
 echo
# 
#
 echo
 echo - - - - Executing program: $SLE_SOLVER  
 $SLE_SOLVER $CONFIG_NAME
 date 
 echo
 echo +----------------------------------+   
 echo +------ Done solving the SLE ------+ 
 echo +----------------------------------+    
 echo 
#############################################################
 echo
 echo +:::::::::::::::::::::::::::::+  
 echo +:::::: Post-processing ::::::+ 
 echo +:::::::::::::::::::::::::::::+   
 echo 
 date 
#
if [ -d $DEPOT_NAME ]
then
echo ; echo ---- WARNING: Data are overwritten in folder $DEPOT_NAME ; echo  
fi 
#
if [ ! -d $DEPOT_NAME ]
then
echo ; echo ---- Folder $DEPOT_NAME is created with subfolders ; echo
fi
mkdir -p $DEPOT_NAME
mkdir -p $DEPOT_NAME/CON 
mkdir -p $DEPOT_NAME/OFU
mkdir -p $DEPOT_NAME/TOP 
mkdir -p $DEPOT_NAME/ICE 
mkdir -p $DEPOT_NAME/FPR
mkdir -p $DEPOT_NAME/RSL
mkdir -p $DEPOT_NAME/TGS
mkdir -p $DEPOT_NAME/GEO
mkdir -p $DEPOT_NAME/STK
mkdir -p $DEPOT_NAME/PMT
mkdir -p $DEPOT_NAME/BIN 
mkdir -p $DEPOT_NAME/LOA
#
#
# //////////////////////////////////////////////////////
#
#
#
 echo - - - - Executing program: $PPR_PROGRAM
 $PPR_PROGRAM $CONFIG_NAME
#
#
# ------------------------------------------------------  CONFIGURATION 
echo
TO=$DEPOT_NAME/CON
echo ---- Copying the CONFIG data into folder: $TO
cp $CONFIG_NAME $TO
#cp make_sle.sh $TO
#cp  $SLE_SOLVER $TO
#cp  $PPR_PROGRAM $TO
#
#
# ------------------------------------------------------  OCEAN FUNCTION
TO=$DEPOT_NAME/OFU
echo ---- Moving the OF data into folder: $TO
if [ -e ice_floating.000.0.dat ]       ; then mv ice_floating*.dat $TO ; fi 
if [ -e ice_grounded_below.000.0.dat ] ; then mv ice_grounded_below*.dat $TO ; fi 
if [ -e ice_grounded_above.000.0.dat ] ; then mv ice_grounded_above*.dat $TO ; fi 
if [ -e continent.000.0.dat ]          ; then mv continent*.dat    $TO ; fi 
if [ -e ocean.000.0.dat ]              ; then mv ocean.*.dat       $TO ; fi 
cp ./GMT_scripts/of-REV4-nn*.gmt $TO 
#
#
# ------------------------------------------------------  TOPOGRAPHY
TO=$DEPOT_NAME/TOP
echo ---- Moving the TOPO data into folder: $TO
if [ -e topo.000.0.dat ]          ; then mv topo*.dat        $TO ; fi 
cp ./GMT_scripts/topo-REV2-nn42.gmt $TO
cp ./GMT_scripts/topo-REV2-nn52.gmt $TO
#
#
#
# ------------------------------------------------------  ICE DATA
TO=$DEPOT_NAME/ICE
echo ---- Moving the ICE data into folder: $TO
if [ -e ice.000.0.dat ]          ; then mv ice.*.dat           $TO ; fi 
cp ./GMT_scripts/ice-REV3-nn42.gmt $TO 
cp ./GMT_scripts/ice-REV3-nn52.gmt $TO 
#
#
#
# ------------------------------------------------------  LOAD FUNCTION
TO=$DEPOT_NAME/LOA
echo ---- Moving the LOAD FUNCTION data into folder: $TO
if [ -e load.000.0.dat ]          ; then mv load.*.dat           $TO ; fi 
#
#
#
# ------------------------------------------------------  FINGERPRINTS
TO=$DEPOT_NAME/FPR
echo ---- Moving the FINGERPRINTS data into folder: $TO
if [ -e sdot.pix ]          ; then mv sdot.pix           $TO ; fi 
if [ -e udot.pix ]          ; then mv udot.pix           $TO ; fi 
if [ -e ndot.pix ]          ; then mv ndot.pix           $TO ; fi 
if [ -e gdot.pix ]          ; then mv gdot.pix           $TO ; fi 
if [ -e ldot.pix ]          ; then mv ldot.pix           $TO ; fi 
if [ -e vdot_east.pix ]     ; then mv vdot_east.pix      $TO ; fi 
if [ -e vdot_sout.pix ]     ; then mv vdot_sout.pix      $TO ; fi 
if [ -e edot.pix ]          ; then mv edot.pix           $TO ; fi 
if [ -e fps-stats.dat ]     ; then mv fps-stats.dat      $TO ; fi 
cp ./GMT_scripts/fps-sle.gmt $TO 
#
#
#
# ------------------------------------------------------  RELATIVE SEA LEVEL
TO=$DEPOT_NAME/RSL
echo ---- Moving the RSL data into folder: $TO
if [ -e rslp-101.dat ]          ; then mv rslp*.dat          $TO ; fi 
if [ -e rsld-101.dat ]          ; then mv rsld*.dat          $TO ; fi 
if [ -e RSL.DAT ]               ; then mv RSL.DAT            $TO ; fi 
cp ./GMT_scripts/rsl-MULTI.gmt $TO
cp ./GMT_scripts/rsl-SINGLE.gmt $TO
#
#
#
# ------------------------------------------------------  TIDE GAUGES
TO=$DEPOT_NAME/TGS
echo ---- Moving the TG data into folder: $TO
if [ -e tg.dat ]            ; then mv tg.dat              $TO ; fi 
#
#
#
# ------------------------------------------------------  GEODETIC VELOCITIES
TO=$DEPOT_NAME/GEO
echo ---- Moving the geodetic velocities into folder: $TO
if [ -e velocities.dat ]    ; then mv velocities.dat      $TO ; fi 
#
#
#
# ------------------------------------------------------  STOKES COEFFICIENTS
TO=$DEPOT_NAME/STK
echo ---- Moving the STOKES COEFF into folder: $TO
if [ -e stokes.dat ]            ; then mv stokes.dat      $TO ; fi 
#
#
# ------------------------------------------------------  POLAR MOTION
TO=$DEPOT_NAME/PMT
echo ---- Moving the POLAR MOTION data into folder: $TO
if [ -e m.dot ]            ; then mv m.dat            $TO ; fi 
if [ -e m.dot ]            ; then mv m.dot            $TO ; fi 
#
 echo
 echo +::::::::::::::::::::::::::::::::::::+   
 echo +:::::: End of post-processing ::::::+ 
 echo +::::::::::::::::::::::::::::::::::::+    
 echo 
 date 

