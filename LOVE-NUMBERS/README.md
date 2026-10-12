# SELEN Green's functions calculator

This directory contains utilities based upon the TABOO code (https://github.com/danielemelini/TABOO) which compute SELEN surface and rotational Green's functions for a given user-supplied rheological profile.

To build the code, type the `make` command inside the `src/` directory. 

The SELEN Green's functions calcuator uses the FMLIB multi-precision arithmetic library by D.M. Smith (https://dmsmith.lmu.build/), whose source code is included in the `src/` directory. A copy of the FMLIB license is also included in the `src/FMLIB-LICENSE.txt` file.

# Defining the Earth model

The Earth model is described in a data file. An example model file is listed below:

```
#----------------------------------------------------------------------------------------------------
# DM 26-Nov-2025 / PREM_LT60_NV5_VM5a.dat / From Table 2 of Roy&Peltier 2017
#      r<(m)          r>(m)          rho(kg/m^3)    rig(Pa)       visco(Pa.s)   Thick(km)     Layer
#----------------------------------------------------------------------------------------------------
   6   6.31100E+06    6.37100E+06    3.09544E+03    5.3767E+10    1.0000E+99    6.0000E+01    LT
   5   6.27100E+06    6.31100E+06    3.37472E+03    6.7446E+10    1.0000E+22    4.0000E+01    M1
   4   5.95100E+06    6.27100E+06    3.44828E+03    7.3874E+10    5.0000E+20    3.2000E+02    M2
   3   5.70100E+06    5.95100E+06    3.89413E+03    1.1088E+11    5.0000E+20    2.5000E+02    M3
   2   5.11100E+06    5.70100E+06    4.54824E+03    1.8321E+11    1.5700E+21    5.9000E+02    M4
   1   3.48550E+06    5.11100E+06    5.10252E+03    2.4446E+11    3.2300E+21    1.6255E+03    M5
   0   0.00000E+00    3.48550E+06    1.09626E+04    8.4634E+09    0.0000E+00    3.4855E+03    CO
```

The first four lines are intended as headers describing the contents of the file and will not be used. Each of the remaining lines describes a model layer. Layers are ordered from top to bottom (from the surface to the core). The columns contain the following information:

* Column 1: Layer number (from $N-1$ to $0$, where $N$ is the total number of layers)
* Column 2 & 3: Bottom and top radii of the layer (in m)
* Column 4: Density $\rho$ (in kg/m^3)
* Column 5; Rigidity $\mu$ (i.e., shear modulus, in Pa)
* Column 6 Viscosity $\eta$ (in Pa*s). Its value is ignored for the outermost layer (which is assumed to be elastic, i.e. with infinite viscosity) and for the inner core (which is assumed to be an inviscid fluid, i.e. with zero viscosity)
* Column 7: Thickness of the layer, in km. The value is present for the user's convenience and is not read by the program.
* Column: 8 Layer label. This column is present for user’s convenience and it is not read by the program.

Example model files for rheological profiles VM1, VM2a and VM5a are available in the `MODELS/` folder.

# Configuring the Green's functions calculator

The parameters controlling the generation of the Green's functions are specified in a configuration file. An example configuration file is listed below.

```
!  -------------------------------------------------------
!  This is the configuration file for procedure "TABOO"
!  -------------------------------------------------------
!
! >>>>>>> Maximum harmonic degrees
!
  512
!
! >>>>>>> Multi-precision parameter 
!
  256
!
! >>>>>>> Number of time steps and their length in kyrs
!
  52 
  0.5   
!
! >>>>>>> Number of mantle layers and model file
!
   5 
!
   MODELS/PREM_LT60_NV5_VM5a.dat    
!
! ================ End of configuration file ================
```

Lines starting with an exclamation mark (`!`) are interpreted as comments and are therefore ignored by the program. The parameters in the configuration file are, in order of appearance:
* the maximum harmonic degree (`512` in the example above). 
* the number of significant digits used in the computation (`128` in the example).
* the number of time steps and the time step length in kyrs (`60` and `0.5` in the example above). Note that the time step length shall match the temporal discretization step employed in the ice model and the number of time steps should be at least equal to the number of time steps in the ice model.
* the number $N_V$ of visco-elastic layers in the model (`5` in the example avove).
* the full path to the model file (`MODELS/PREM_LT90_NV5_UM1_TZ1_LM3_ANU.dat`). Note that the number of mantle layers specified in the config file shall be coherent with the total number of layers in the model file ($N=N_V+2$, where $N$ is the total number of layers in the model file and $N_V$ is the number of mantle layers specified in the config file).

Example configuration files for generating the surface and rotational Green's functions for models VM1, VM2a and VM5a are available in the `CONFIGS/` folder.

# Running the Green's functions calculator

To compute the Green's functions for a given configuration, type the following command:

```
sh TABOO.sh config_file
```

where `config_file` is the full path to the configuration file described above. For instance:

```
sh TABOO.sh CONFIGS/taboo_VM5a.cfg
```

will compute GFs for the VM5a viscosity model, suitable for use with the ICE-6G GIA model, according to the parameters specified in the `CONFIGS/taboo_VM5a.cfg` configuration file.

After the program has run, the following files will be made available:

| File| Contents |
| ----| -------- |
| `ebg.dat`  | surface Green's function for geoid height change |
| `ebs.dat`  | surface Green's function for sea level change |
| `ebu.dat`  | surface Green's function for vertical displacement |
| `ebv.dat`  | surface Green's function for horizontal displacement |
| `gamma.dat` | rotational Green's function according to the traditional theory |
| `gamma-REV.dat` | rotational Green's function according to the revised theory |
| `PMTF.dat`   | polar motion transfer function according to the traditional rotational theory |
| `PMTF-REV.dat` | polar motion transfer function according to the revised rotational theory |
| `h-load.dat` | $h_n$ loading Love Numbers in multi-exponential form |
| `l-load.dat` | $l_n$ loading Love Numbers in multi-exponential form |
| `k-load.dat` | $l_n$ loading Love Numbers in multi-exponential form |
| `h-tide.dat` | $h_2$ tidal Love Number in multi-exponential form |
| `l-tide.dat` | $l_2$ tidal Love Number in multi-exponential form |
| `k-tide.dat` | $k_2$ tidal Love Number in multi-exponential form |
| `l2-data.dat`  | degree 2 Love Numbers used to compute the PMTF |
| `fluid-ln.dat` | fluid tidal LNs assuming a zero effective lithosphere thickness |
| `run-model-parameters.dat` | a copy of the Earth model file |
| `spectrum.dat` | relaxation spectrum of the model |
| `taboo.cfg` | a copy of the configuration file |
| `taboo.log` | log file of the TABOO module |

These files can be stored in a folder by typing
```
sh pproc.sh folder_name
```
In that case a folder called `folder_name` will be created and all the output files will be moved inside that folder.

Otherwise, the working directory can be purged with the command
```
sh purge.sh
```
which will delete all the output files listed above.
