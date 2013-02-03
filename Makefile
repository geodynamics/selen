#=====================================================================
#
#          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
#          --------------------------------------------------
#
#          Main authors: Dimitri Komatitsch and Jeroen Tromp
#                        Princeton University, USA
#             and University of Pau / CNRS / INRIA, France
# (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
#                            April 2011
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
#=====================================================================

# Makefile.  Generated from Makefile.in by configure.

FC = g95
FCFLAGS = -I${SETUP} -I${B}

FCCOMPILE = ${FC} ${FCFLAGS}

AR = @AR@
ARFLAGS = @ARFLAGS@
RANLIB = @RANLIB@

## compilation directories
# B : build directory
B = .
# E : executables directory
E = ${B}/bin
# O : objects directory
O = ${B}/obj
# S : source file directory
S = ${B}/src
## setup file directory
SETUP = ${B}/tmp
## output file directory
OUTPUT = ${B}/OUTPUT_FILES

#######################################
# solver objects with statically allocated arrays; dependent upon
# values_from_mesher.h

CONFIG_OBJECTS = \
	$O/shtools.o \
	$O/config.o \
	$O/harmonics.o

SLE_OBJECTS = \
	$O/harmonics.o \
	$O/sle.o

# Equivalent sea level
ESL_OBJECTS = \
	$O/harmonics.o \
	$O/esl.o

# Ice sheets contours
MS_OBJECTS = \
	$O/harmonics.o \
	$O/ms.o

# Love numbers tools (TABOO)
TB_OBJECTS = \
	$O/harmonics.o \
	$O/tb.o

# Computation of ice Shape factors and SH dechomposition
SHAPE_FACTORS_OBJECTS = \
	$O/harmonics.o \
	$O/shape_factors.o
SHICE_OBJECTS = \
	$O/harmonics.o \
	$O/shice.o

# Pixelization (i)
PX_OBJECTS = \
	$O/harmonics.o \
	$O/px.o

# Pixelization (ii) 
PXREC_OBJECTS = \
	$O/harmonics.o \
	$O/px_rec.o

# Pixelization partitioning - missing source file px_part.f90
# Copy to local storage - missing source file px_copy.f90
# Parallel wet/dry pixel separation - missing source file px_select.f90

# Spherical harmonics 
SH_OBJECTS = \
	$O/harmonics.o \
	$O/sh.o

# Window function
WNW_OBJECTS = \
	$O/harmonics.o \
	$O/wnw.o


#######################################
####
#### targets
####

# default targets
DEFAULT = \
	config \
	sle \
	px \
	pxrec

default: $(DEFAULT)

all: clean default

backup:
	mkdir -p bak
	cp *f90 *h Makefile bak

bak: backup

req_dirs: bindir objdir tmpdir

bindir:
	mkdir -p bin
	
objdir:
	mkdir -p obj

tmpdir:
	mkdir -p tmp

#######################################
####
#### rules for executables
####

config: req_dirs $(CONFIG_OBJECTS)
	${FCCOMPILE} -o ${E}/config $(CONFIG_OBJECTS)

sle: req_dirs $(SLE_OBJECTS)
	${FCCOMPILE} -o ${E}/sle $(SLE_OBJECTS)

esl: req_dirs $(ESL_OBJECTS)
	${FCCOMPILE} -o ${E}/esl $(ESL_OBJECTS)

ms: req_dirs $(MS_OBJECTS)
	${FCCOMPILE} -o ${E}/ms $(MS_OBJECTS)

px: req_dirs $(PX_OBJECTS)
	${FCCOMPILE} -o ${E}/px $(PX_OBJECTS)

pxrec: req_dirs $(PXREC_OBJECTS)
	${FCCOMPILE} -o ${E}/pxrec $(PXREC_OBJECTS)

tb: req_dirs $(TB_OBJECTS)
	${FCCOMPILE} -o ${E}/tb $(TB_OBJECTS)

clean:
	rm -rf $O/* bin obj tmp

#######################################
####
#### rule for each .o file below
####

$O/config.o: $S/harmonics.f90 $S/config.f90 $S/shtools.f90
	${FCCOMPILE} -c -o $O/config.o $S/config.f90

$O/esl.o: ${SETUP}/data.inc $S/harmonics.f90 $S/esl.f90
	${FCCOMPILE} -c -o $O/esl.o $S/esl.f90

$O/geo.o: ${SETUP}/data.inc $S/harmonics.f90 $S/geo.f90
	${FCCOMPILE} -c -o $O/geo.o $S/geo.f90

$O/gmaps.o: ${SETUP}/data.inc $S/harmonics.f90 $S/gmaps.f90
	${FCCOMPILE} -c -o $O/gmaps.o $S/gmaps.f90

$O/harmonics.o: $S/harmonics.f90
	${FCCOMPILE} -c -o $O/harmonics.o $S/harmonics.f90

$O/ms.o: ${SETUP}/data.inc $S/harmonics.f90 $S/ms.f90
	${FCCOMPILE} -c -o $O/ms.o $S/ms.f90

$O/of_dv.o: ${SETUP}/data.inc $S/harmonics.f90 $S/of_dv.f90
	${FCCOMPILE} -c -o $O/of_dv.o $S/of_dv.f90

$O/px.o: ${SETUP}/data.inc $S/harmonics.f90 $S/px.f90
	${FCCOMPILE} -c -o $O/px.o $S/px.f90

$O/px_rec.o: ${SETUP}/data.inc $S/px_rec.f90
	${FCCOMPILE} -c -o $O/px_rec.o $S/px_rec.f90

$O/rec_ice.o: ${SETUP}/data.inc $S/harmonics.f90 $S/rec_ice.f90
	${FCCOMPILE} -c -o $O/rec_ice.o $S/rec_ice.f90

$O/rec_of.o: ${SETUP}/data.inc $S/harmonics.f90 $S/rec_of.f90
	${FCCOMPILE} -c -o $O/rec_of.o $S/rec_of.f90

$O/rmaps.o: ${SETUP}/data.inc $S/harmonics.f90 $S/rmaps.f90
	${FCCOMPILE} -c -o $O/rmaps.o $S/rmaps.f90

$O/rsl.o: ${SETUP}/data.inc $S/harmonics.f90 $S/rsl.f90
	${FCCOMPILE} -c -o $O/rsl.o $S/rsl.f90

$O/rsl_zones.o: ${SETUP}/data.inc $S/harmonics.f90 $S/rsl_zones.f90
	${FCCOMPILE} -c -o $O/rsl_zones.o $S/rsl_zones.f90

$O/rslc.o: ${SETUP}/data.inc $S/harmonics.f90 $S/rslc.f90
	${FCCOMPILE} -c -o $O/rslc.o $S/rslc.f90

$O/sh.o: ${SETUP}/data.inc $S/harmonics.f90 $S/sh.f90
	${FCCOMPILE} -c -o $O/sh.o $S/sh.f90

$O/sh_of.o: ${SETUP}/data.inc $S/harmonics.f90 $S/sh_of.f90
	${FCCOMPILE} -c -o $O/sh_of.o $S/sh_of.f90

$O/sh_rsl.o: ${SETUP}/data.inc $S/harmonics.f90 $S/sh_rsl.f90
	${FCCOMPILE} -c -o $O/sh_rsl.o $S/sh_rsl.f90

$O/sh_rslc.o: ${SETUP}/data.inc $S/harmonics.f90 $S/sh_rslc.f90
	${FCCOMPILE} -c -o $O/sh_rslc.o $S/sh_rslc.f90

$O/sh_tgauges.o: ${SETUP}/data.inc $S/harmonics.f90 $S/sh_tgauges.f90
	${FCCOMPILE} -c -o $O/sh_tgauges.o $S/sh_tgauges.f90

$O/shape_factors.o: ${SETUP}/data.inc $S/harmonics.f90 $S/shape_factors.f90
	${FCCOMPILE} -c -o $O/shape_factors.o $S/shape_factors.f90

$O/shice.o: ${SETUP}/data.inc $S/harmonics.f90 $S/shice.f90
	${FCCOMPILE} -c -o $O/shice.o $S/shice.f90

$O/shtools.o: $S/harmonics.f90 $S/shtools.f90
	${FCCOMPILE} -c -o $O/shtools.o $S/shtools.f90

$O/sle.o: ${SETUP}/data.inc $S/harmonics.f90 $S/sle.f90
	${FCCOMPILE} -c -o $O/sle.o $S/sle.f90

$O/stokes.o: ${SETUP}/data.inc $S/harmonics.f90 $S/stokes.f90
	${FCCOMPILE} -c -o $O/stokes.o $S/stokes.f90

$O/tb.o: ${SETUP}/data.inc $S/harmonics.f90 $S/tb.F90
	${FCCOMPILE} -c -o $O/tb.o $S/tb.F90

$O/tgauges.o: ${SETUP}/data.inc $S/harmonics.f90 $S/tgauges.f90
	${FCCOMPILE} -c -o $O/tgauges.o $S/tgauges.f90

$O/wnw.o: ${SETUP}/data.inc $S/harmonics.f90 $S/wnw.f90
	${FCCOMPILE} -c -o $O/wnw.o $S/wnw.f90

${SETUP}/data.inc: config config.dat
	(cd tmp ; ../${E}/config config.dat ${SETUP}/data.inc ; cd ..)

