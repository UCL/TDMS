# General setup ---------------------------------------------------------------
cmake_minimum_required(VERSION 3.21)

project(tdms LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 14)
include(CTest)

# Allow RPATH on mac
set(CMAKE_MACOSX_RPATH TRUE)

# Append the cmake/ directory to the search path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/targets.cmake)

# and the  include directories for the project
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/matlabio/)


# Command line options for CMake ----------------------------------------------
# Type of derivatives to use
set(DERIVATIVE_TYPE "PS" CACHE STRING "Type of derivatives to compute.
Either Use pseudo spectral (PS) or finite differences (FD).
Options: {FD, PS}")
set_property(CACHE DERIVATIVE_TYPE PROPERTY STRINGS "FD" "PS")
get_property(x CACHE DERIVATIVE_TYPE PROPERTY STRINGS)
if (NOT DERIVATIVE_TYPE IN_LIST x)
    message(FATAL_ERROR "Cannot set DERIVATIVE_TYPE=${DERIVATIVE_TYPE}. Options: ${x}")
endif()

if (${DERIVATIVE_TYPE} STREQUAL FD)
    set(DFLAG "-DFDFLAG")
endif()


# Matlab ----------------------------------------------------------------------
find_package(Matlab REQUIRED MAT_LIBRARY)
include_directories(${Matlab_INCLUDE_DIRS})

# Currently matlab only has x84 libraries on Mac (no arm)
set(CMAKE_OSX_ARCHITECTURES "x86_64")

# Set the RPATH for the executable to find the dynamically linked libraries
get_filename_component(MATLAB_LIB_ROOT "${Matlab_MAT_LIBRARY}" DIRECTORY)
set(CMAKE_INSTALL_RPATH ${MATLAB_LIB_ROOT})


# FFTW3 -----------------------------------------------------------------------
find_package(FFTW REQUIRED)
include_directories(${FFTW_INCLUDE_DIRS})


# OpenMP + libc++ -------------------------------------------------------------
if (CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    find_package(OpenMP_Mac REQUIRED)
    find_package(Lib_CXX_Mac)
else()
    find_package(OpenMP REQUIRED)
endif()


# TDMS target -----------------------------------------------------------------
set(SOURCES
    src/fields/base.cpp
    src/fields/electric.cpp
    src/fields/magnetic.cpp
    src/fields/split.cpp
    src/argument_parser.cpp
    src/fdtd_grid_initialiser.cpp
    src/grid_labels.cpp
    src/interpolate.cpp
    src/iterator.cpp
    src/interpolation_methods.cpp
    src/matrix_collection.cpp
    src/mesh_base.cpp
    src/numeric.cpp
    src/numerical_derivative.cpp
    src/simulation_parameters.cpp
    src/timer.cpp
    src/utils.cpp
    matlabio/matlabio.cpp
    )

if (BUILD_TESTING)
    test_target()
else()
    release_target()
endif()


# Compile options -------------------------------------------------------------
target_compile_options(tdms PUBLIC -DMX_COMPAT_32 -c ${DFLAG} -O3)

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_compile_options(tdms PRIVATE /W4)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(tdms PRIVATE -Wall -Wno-maybe-uninitialized)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    target_compile_options(tdms PRIVATE -Wall)
endif()


# Install ---------------------------------------------------------------------
install(TARGETS tdms)

# Testing ---------------------------------------------------------------------
if (BUILD_TESTING)
    enable_testing()

    add_executable(tdms_tests
            tests/unit/test_openandorder.cpp
            tests/unit/test_fields.cpp
            tests/unit/test_interpolate.cpp)

    target_link_libraries(tdms_tests PRIVATE
            Catch2::Catch2WithMain
            tdms_lib
            )

    add_test(test_all tdms_tests)
endif()
