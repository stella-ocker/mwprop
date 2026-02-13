# Fortran Version of NE2025

The Fortran source code of NE2025 is provided here for users wishing for speed and/or the highest level of computational fidelity. This directory additionally provides a copy of the calibration data used to fit NE2025, as well as instructions for reading in the data (see cal.NE2025).

-----

## Fortran Installation

The Fortran code is built using gfortran (ensure gcc is up-to-date). These build commands have been tested on Mac OS. Please submit any issues encountered during installation so that we can maintain a record of bug fixes for future users.

    cd src.NE2025
    gfortran -ffixed-line-length-0 -c -o NE2025.o NE2025.f
    gfortran -ffixed-line-length-0 -c -o dmdsm.NE2025.o dmdsm.NE2025.f
    gfortran -ffixed-line-length-0 -c -o density.NE2025.o density.NE2025.f # this may produce several warnings; these are expected and do not affect subsequent steps
    gfortran -ffixed-line-length-0 -c -o neLISM.NE2025.o neLISM.NE2025.f
    gfortran -ffixed-line-length-0 -c -o neclumpN.NE2025.o neclumpN.NE2025.f
    gfortran -ffixed-line-length-0 -c -o nevoidN.NE2025.o nevoidN.NE2025.f
    gfortran -ffixed-line-length-0 -c -o scattering98.o scattering98.f
    gfortran NE2025.o dmdsm.NE2025.o density.NE2025.o neLISM.NE2025.o neclumpN.NE2025.o nevoidN.NE2025.o scattering98.o -Wl,-rpath,/anaconda3/lib -o ../bin.NE2025/NE2025 # modify library path as needed

-----

## Usage

    cd bin.NE2025
    ./NE2025 l   b   DM/D   ndir 
