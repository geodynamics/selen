# SELEN4
SealEveL EquatioN solver - version 4

## Features

This repository contains SELEN, an open-source Fortran code for the numerical solution of the *Sea Level Equation* (SLE) for a spherical, layered, rotating Earth with viscoelastic rheology.

SELEN implements a *gravitationally* and *topographically* self-consistent SLE and can compute several quantities of interest for the Glacial Isostatic Adjustment (GIA) problem, including Relative Sea-Level (RSL) curves, present-day sealevel rates at tide gauges, surface deformations and perturbations of the gravity field, geodetic *fingerprints* and paleo-topography maps.

SELEN includes portions of the [SHTOOLS library](https://shtools.oca.eu/shtools/) by Mark A. Wieczorek and Matthias Meschede and subroutines from Max Tegmark for the icosahedron-shaped pixelization of the sphere.

SELEN is distributed under the 3-clause BSD license.

## Installation

To install SELEN, just download the repository to a local directory.  Some files in the `DATA/` directory are compressed to reduce disk usage; before running SELEN for the first time they should be expanded by typing `gunzip *.gz` inside the `DATA/` directory.  The SELEN executables can be built by typing `make` in the `src/` directory. For further details, see the SELEN User guide `User_Manual.pdf`.



