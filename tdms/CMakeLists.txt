project(tdms)

# General setup ---------------------------------------------------------------
cmake_minimum_required(VERSION 3.22)

# Allow RPATH on mac
set(CMAKE_MACOSX_RPATH 1)

# Append the cmake/ directory to the search path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

# and the  include directories for the project
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/matlabio/)
include_directories(${CMAKE_SOURCE_DIR}/mesh/)


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

# Currently matlab only has x84 libraries on Mac (no arm)
set(CMAKE_OSX_ARCHITECTURES "x86_64")

# Set the RPATH for the executable to find the dynamically linked libraries
get_filename_component(MATLAB_LIB_ROOT "${Matlab_MAT_LIBRARY}" DIRECTORY)
set(CMAKE_INSTALL_RPATH ${MATLAB_LIB_ROOT})
message(STATUS "rpath: ${CMAKE_INSTALL_RPATH}")

# FFTW3 -----------------------------------------------------------------------
find_package(FFTW REQUIRED)
include_directories(${FFTW_INCLUDE_DIRS})
include_directories(${Matlab_INCLUDE_DIRS})


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
target_sources(tdms PRIVATE
        src/openandorder_PSTD.cpp
        src/interpolate.cpp
        src/numeric.cpp
        src/queuelinked.cpp
        src/tdms_iterator.cpp
        src/mesh_base.cpp
        src/numerical_derivative.cpp
        matlabio/matlabio.cpp
        )

target_link_libraries(tdms
        Matlab::mex
        Matlab::mx
        Matlab::mat
        ${FFTW_LIBRARIES}
        ${LIBCXX_LIBRARY}
        OpenMP::OpenMP_CXX
        )

install(TARGETS tdms)
