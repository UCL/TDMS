project(tdms)

# General setup ---------------------------------------------------------------
cmake_minimum_required(VERSION 3.22)

# Append the cmake/ directory to the search path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

# and the  include directories for the project
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/matlabio/)
include_directories(${CMAKE_SOURCE_DIR}/mesh/)


# Command line options for CMake ----------------------------------------------
# Type of I/O to use
set(IO_TYPE "MATLAB" CACHE STRING "Use MATLAB or C++ input/output")
set_property(CACHE IO_TYPE PROPERTY STRINGS "MATLAB" "OWN")
get_property(x CACHE IO_TYPE PROPERTY STRINGS)
if (NOT IO_TYPE IN_LIST x)
    message(FATAL_ERROR "Cannot set IO_TYPE=${IO_TYPE}. Options: ${x}")
endif()

# Type of derivatives to use
set(DERIVATIVE_TYPE "PS" CACHE STRING "Type of derivatives to compute.
Either Use pseudo spectral (PS) or finite differences (FD).
Options: {FD, PS}")
set_property(CACHE DERIVATIVE_TYPE PROPERTY STRINGS "FD" "PS")
get_property(x CACHE DERIVATIVE_TYPE PROPERTY STRINGS)
if (NOT DERIVATIVE_TYPE IN_LIST x)
    message(FATAL_ERROR "Cannot set DERIVATIVE_TYPE=${DERIVATIVE_TYPE}. Options: ${x}")
endif()


# Matlab ----------------------------------------------------------------------
if (${IO_TYPE} STREQUAL MATLAB)
    set(DFLAG "-D MIOFLAG")

    find_package(Matlab REQUIRED MAT_LIBRARY)
    # Currently matlab only has x84 libraries on Mac (no arm)
    set(CMAKE_OSX_ARCHITECTURES "x86_64")
    set(NUMERICAL_DERIVATIVE_SOURCE_FILE "src/numerical_derivative.cpp")
else()
    set(CXX_MATLAB_IO_SOURCE_FILE "matlabio/Mat_io.cpp")
endif()


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
        src/cli.cpp
        src/interpolate.cpp
        src/io.cpp
        src/numeric.cpp
        src/queuelinked.cpp
        src/tdms_iterator.cpp
        matlabio/matlabio.cpp
        mesh/mesh_base.cpp
        ${NUMERICAL_DERIVATIVE_SOURCE_FILE}
        ${CXX_MATLAB_IO_SOURCE_FILE}
        )

target_link_libraries(tdms
        ${Matlab_LIBRARIES}
        ${FFTW_LIBRARIES}
        ${LIBCXX_LIBRARY}
        OpenMP::OpenMP_CXX
        )

install(TARGETS tdms)
