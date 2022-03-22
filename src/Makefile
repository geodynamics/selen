#
#
FC=gfortran
FOPTS=-m64 -O
#
CURRENT_REVISION=REV22
#
SLE_SOURCES=sle-$(CURRENT_REVISION).f90 stop_config.f90 mj.f90 lj.f90 dom.f90 j_index.f90 config.f90
#
PPR_SOURCES=pproc-sle-$(CURRENT_REVISION).f90 stop_config.f90 mj.f90 lj.f90 dom.f90 \
            plmbar_mod.f90 plmbar.f90 harmo.f90 j_index.f90 config.f90 \
            y.f90  y_t.f90 y_l.f90 y_tt.f90 y_ll.f90 y_ct.f90 y_cl.f90 y_tl.f90 \
            plmbar_mod_d1.f90 plmbar_d1.f90                                      
#
SH_SOURCES=sha-REV2.f90 sindd.f90 cosdd.f90 plmbar.f90 j_index.f90 plmbar_mod.f90 stop_config.f90 mj.f90
#
all: ../sle.exe ../ppr.exe ../sha.exe
#
clean:
	$(RM) ../sle.exe ../ppr.exe ../sha.exe *.mod *.o
#
../sle.exe: $(SLE_SOURCES) shtools.mod parameters.inc
	$(FC) $(FOPTS) -fopenmp $(SLE_SOURCES) -o ../sle.exe
#
../ppr.exe: $(PPR_SOURCES)
	$(FC) $(FOPTS) $(PPR_SOURCES) -o ../ppr.exe
#
../sha.exe: $(SH_SOURCES)
	$(FC) $(FOPTS) -fopenmp $(SH_SOURCES) -o ../sha.exe
#
shtools.mod: shtools.f90
	$(FC) $(FOPTS) -c shtools.f90
#
