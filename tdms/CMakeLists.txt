project(.)
cmake_minimum_required(VERSION 3.22)

# Append the cmake/ directory to the search path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
        "${CMAKE_SOURCE_DIR}/cmake/")

# and the  include directories for the project
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/matlabio/)
include_directories(${CMAKE_SOURCE_DIR}/mesh/)

# Command line options for CMake ----------------------------------------------
option(MATLAB_IO "Use MATLAB input/output" ON)
option(OWN_IO "Use inbuilt input/output" OFF)
option(USE_PS "Use pseudo spectral" ON)
option(USE_FD "Use finite differences" OFF)

if (${MATLAB_IO} STREQUAL ON)
    set(DFLAG -DMIOFLAG)  # TODO: Something better. Copied from Makefile
endif()

# Matlab ----------------------------------------------------------------------
find_package(Matlab REQUIRED MAT_LIBRARY)
# Currently matlab only has x84 libraries on Mac (no arm) so hard code that arch
set(CMAKE_OSX_ARCHITECTURES "x86_64")


# FFTW3 -----------------------------------------------------------------------
find_package(FFTW REQUIRED)
include_directories(${FFTW_INCLUDE_DIRS})
include_directories(${Matlab_INCLUDE_DIRS})

# OpenMP ----------------------------------------------------------------------
find_package(OpenMP REQUIRED)

add_compile_options(-DMX_COMPAT_32 -c ${DFLAG} -O3)

add_executable(tdms)
target_sources(tdms PRIVATE
        src/cli.cpp
        src/interpolate.cpp
        src/io.cpp
        src/numeric.cpp
        src/numerical_derivative.cpp
        src/queuelinked.cpp
        src/tdms_iterator.cpp
        matlabio/matlabio.cpp
        mesh/mesh.cpp
        mesh/mesh_base.cpp
        )


target_link_libraries(tdms
        ${Matlab_LIBRARIES}
        ${FFTW_LIBRARIES}
        OpenMP::OpenMP_CXX
        )

install(TARGETS tdms)
