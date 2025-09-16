cmake_minimum_required(VERSION 3.16) # 3.30
project(atomistic_spin_model_GdFe LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Export compile_commands.json (clangd/VS Code)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type specified; default to RelWithDebInfo")
    set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Choose build type")
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
            "Debug" "Release" "RelWithDebInfo" "MinSizeRel")
endif()
message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

# Target
add_executable(atomistic_spin_model_GdFe main.cpp
        integrator.cpp
        integrator.h
        fields.cpp
        fields.h
        lattice.cpp
        lattice.h
        reductions.cpp
        reductions.h
        io.cpp
        io.h
        io_temperature_csv.cpp
        io_temperature_csv.h
        params.h
        math_utils.h
        rng.h
        io_csv_utils.h
        init.h)

target_include_directories(atomistic_spin_model_GdFe PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}")

target_compile_definitions(atomistic_spin_model_GdFe PRIVATE
        $<$<CONFIG:Debug>:DEBUG_BUILD>                    # custom macro
        $<$<CONFIG:RelWithDebInfo>:RELWITHDEBINFO_BUILD>  # custom macro
        $<$<CONFIG:Release>:NDEBUG>
        $<$<CONFIG:MinSizeRel>:NDEBUG MINSIZEREL_BUILD>   # standard, custom
)

# Warnings per compiler
target_compile_options(atomistic_spin_model_GdFe PRIVATE
        $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra -Wpedantic>
        $<$<CXX_COMPILER_ID:MSVC>:/W4 /permissive->)

