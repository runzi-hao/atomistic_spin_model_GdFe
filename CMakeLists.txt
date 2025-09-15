cmake_minimum_required(VERSION 3.30)
project(atomistic_spin_model_GdFe)

set(CMAKE_CXX_STANDARD 20)

add_executable(atomistic_spin_model_GdFe main.cpp
        integrator.cpp
        integrator.h
        fields.cpp
        fields.h
        params.h
        math_utils.h
        lattice.cpp
        lattice.h
        rng.h
        reductions.cpp
        reductions.h
        io.cpp
        io.h
        temperature_series.h
        io_temperature_csv.cpp
        io_temperature_csv.h
        io_csv_utils.h
        init.h)
