project(tdms)

# General setup ---------------------------------------------------------------
cmake_minimum_required(VERSION 3.21)

# Allow RPATH on mac
set(CMAKE_MACOSX_RPATH TRUE)

# Append the cmake/ directory to the search path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

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
if (${CMAKE_SYSTEM_NAME} STREQUAL Darwin)
    find_package(OpenMP_Mac REQUIRED)
    find_package(Lib_CXX_Mac)
else()
    find_package(OpenMP REQUIRED)
endif()

# Compile options -------------------------------------------------------------
add_compile_options(-DMX_COMPAT_32 -c ${DFLAG} -O3)

add_executable(tdms)
set_property(TARGET tdms PROPERTY CXX_STANDARD 11)

target_sources(tdms PRIVATE
        matlabio/matlabio.cpp
        src/argument_parser.cpp
        src/fdtd_grid_initialiser.cpp
        src/interpolate.cpp
        src/iterator.cpp
        src/matrix_collection.cpp
        src/mesh_base.cpp
        src/numeric.cpp
        src/numerical_derivative.cpp
        src/openandorder.cpp
        src/utils.cpp
        )

target_link_libraries(tdms
        FFTW::Double
        ${Matlab_MEX_LIBRARY}
        ${Matlab_MX_LIBRARY}
        ${Matlab_MAT_LIBRARY}
        ${LIBCXX_LIBRARY}
        OpenMP::OpenMP_CXX
        )

if(MSVC)
  target_compile_options(tdms PRIVATE /W4)
else()
  target_compile_options(tdms PRIVATE -Wall -Wno-maybe-uninitialized)
endif()

install(TARGETS tdms)
