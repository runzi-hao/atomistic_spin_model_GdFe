# Atomistic Spin Model for FeGd (C++ Implementation)

## Overview
This project aims to reproduce Figure 3 in Radu, I., Vahaplar, K., Stamm, C. et al. 
Nature 472, 205–208 (2011). https://doi.org/10.1038/nature09901
It implements an atomistic spin dynamics model for a binary alloy (Fe-Gd).  
It supports:
- Building FCC lattice neighbors and assigning atom species
- Computing effective magnetic fields (exchange, anisotropy, thermal noise, total)
- Time evolution using LLG integrator (Heun)
- Bulk magnetization reductions
- Input/output via CSV files
The code is modular and designed for research/educational use.

## Project Structure
- lattice.h/.cpp              : Build FCC neighbors, assign species, count FCC sites, invert linear index 
- fields.h/.cpp               : Compute exchange, anisotropy, thermal, total fields
- integrator.h/.cpp           : Time evolution kernel and normalizations
- reductions.h/.cpp           : Compute bulk magnetizations and fields
- io_temperature_csv.h/.cpp   : Read temperature vs time series
- io.h/.cpp                   : Input/output CSV utilities (settings, helper, bulk properties, neighbors, species)
- test.h/.cpp                 : Unit tests (e.g., atom counts)
- init.h                      : Initialize per-site magnetization, wrapper
- io_csv_utils.h              : CSV parsing helpers to trim string, split a line
- math_utils.h                : Vector normalizations, interpolation
- params.h                    : Constants, data types, control/lattice/species parameters, bulk properties
- rng.h                       : Random number generator wrapper
- main.cpp                    : Main driver with time loop
- temperature_series.h        : Currently not used

## Requirements
- C++20 or newer
- CMake 3.16 or newer
- A C++ compiler (tested with GCC)

## Build Instructions (Bash)
# Configure
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
# Build
cmake --build build -j

# Run (example)
./build/atomistic_spin_model_GdFe 

## Input/Output
Input: 
- input.csv with lattice and species parameter specifications
- Temperature: temperature series csv for electron temperature vs time
Output: 
- bulk_values_vs_time.csv with magnetizations and fields vs time
- nearest_neighbors.txt with neighbor sites for each site
- Gd_sites.txt with all Gd sites

## License
This project is licensed under CC BY-NC 4.0 (non-commercial use only). 
See LICENSE file for details.

## Citation
If you use this code for research, please cite appropriately. 
### Plain text
Runzi Hao. 2025. Atomistic Spin Model for FeGd (C++ Implementation).  
GitHub repository: https://github.com/Ryanne87/atomistic_spin_model_GdFe  
Licensed under CC BY-NC 4.0.

