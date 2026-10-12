#
# Build options for gfortran
FC=gfortran
FOPTS=-O
#
# Build options for ifort
#FC=ifort
#FOPTS=-O -no-wrap-margin
#
# Build options for gfortran/debug
#FC=gfortran
#FOPTS=-g -fbacktrace -fcheck=all
#
PMTF_OBJS=utils.o Balanc_Nlib.o Bubble_Sort_QP.o HQR_Nlib.o HQR_Driver.o Swap_QP.o RH.o Read_DEG2.o
#
TABOO_OBJS=config.o taboo.o utils.o
#
FLUID_LN_OBJS=fluid-ln.o utils.o
#
FMLIB_OBJS=fmsave.o fmzm90.o fm.o
#
BIN_DIR=../bin
#
all: ../bin/taboo.exe ../bin/pmtf.exe ../bin/pmtf-rev.exe ../bin/guv.exe ../bin/guv-rev.exe ../bin/fluid-ln.exe
#
../bin/taboo.exe: $(FMLIB_OBJS) $(TABOO_OBJS) | $(BIN_DIR)
	$(FC) $(FOPTS) $^ -o $@
#
../bin/pmtf.exe: config.o pmtf.o $(FMLIB_OBJS) $(PMTF_OBJS) | $(BIN_DIR)
	$(FC) $(FOPTS) $^ -o $@
#
../bin/pmtf-rev.exe: config.o pmtf-rev.o $(FMLIB_OBJS) $(PMTF_OBJS) | $(BIN_DIR)
	$(FC) $(FOPTS) $^ -o $@
#
../bin/guv.exe: $(FMLIB_OBJS) config.o guv.o utils.o RH.o | $(BIN_DIR)
	$(FC) $(FOPTS) $^ -o $@
#
../bin/guv-rev.exe: $(FMLIB_OBJS) config.o guv-rev.o utils.o RH.o | $(BIN_DIR)
	$(FC) $(FOPTS) $^ -o $@
#
../bin/fluid-ln.exe: $(FMLIB_OBJS) config.o $(FLUID_LN_OBJS) | $(BIN_DIR)
	$(FC) $(FOPTS) $^ -o $@
#
$(BIN_DIR):
	mkdir -p $@
#
clean:
	$(RM) ../bin/*.exe *.mod *.o
#
%.o: %.f
	$(FC) $(FOPTS) -c $<
#
%.o: %.f90
	$(FC) $(FOPTS) -c $<
#
%.o: %.f95
	$(FC) $(FOPTS) -c $<
#
#
