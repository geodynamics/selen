# Create the selen program
bin_PROGRAMS = config selen

# Config requires these source files
config_SOURCES = shtools.f90 config.f90

# Tell how to generate data.inc (requires config and config.dat)
data.inc: config config.dat
	./config

# SELEN requires these source files
selen_SOURCES = data.inc esl.f90 geo.f90 gmaps.f90 harmonics.f90 ms.f90 of_dv.f90 \
		px.f90 px_rec.f90 rec_ice.f90 rec_of.f90 rmaps.f90 rsl.f90 \
		rsl_zones.f90 rslc.f90 sh.f90 sh_of.f90 sh_rsl.f90 sh_rslc.f90 \
		sh_tgauges.f90 shape_factors.f90 shice.f90 shtools.f90 sle.f90 \
		stokes.f90 tb.F90 tgauges.f90 wnw.f90

# Extra files to clean up
clean-local:
	-rm -f data.inc *.grd *.jun *.log *.tmp *.cpt *tmp* junk* *.mod *.old *.bak

