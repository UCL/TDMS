# General setup ---------------------------------------------------------------
cmake_minimum_required(VERSION 3.21)

project(tdms LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)

# Allow building with testing, and default to off
option(BUILD_TESTING "" OFF)
include(CTest)

# Allow user to specify git access protocol (default to https)
option(GIT_SSH "" OFF)
if(GIT_SSH)
    set(GITHUB_PREFIX git@github.com:)
else()
    set(GITHUB_PREFIX https://github.com/)
endif()


# Allow RPATH on mac
set(CMAKE_MACOSX_RPATH TRUE)

# Add -fPIC flag
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Append the cmake/ directory to the search path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/targets.cmake)

# and the  include directories for the project
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/matlabio/)


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


# spdlog ----------------------------------------------------------------------
find_package(spdlog NO_CMAKE_PACKAGE_REGISTRY QUIET)
if(NOT spdlog_FOUND)
    include(FetchContent)
    FetchContent_Declare(spdlog
	    GIT_REPOSITORY ${GITHUB_PREFIX}gabime/spdlog
            GIT_TAG v1.10.0
            )
    FetchContent_MakeAvailable(spdlog)
endif()

# Set the active logging level to debug if required
if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug"
        OR "${CMAKE_BUILD_TYPE}" STREQUAL "RelWithDebInfo")
    add_definitions(-DSPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_DEBUG)
elseif("${CMAKE_BUILD_TYPE}" STREQUAL "Test")
    message(FATAL_ERROR "CMAKE_BUILT_TYPE \"Test\" is not recognised.
    Did you mean -DBUILD_TESTING=ON ?")
    add_definitions(-DSPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_INFO)
else()
    add_definitions(-DSPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_INFO)
endif()


# TDMS target -----------------------------------------------------------------
set(SOURCES
    src/fields/base.cpp
    src/fields/electric.cpp
    src/fields/magnetic.cpp
    src/fields/split.cpp
    src/fields/td_field_exporter_2d.cpp
    src/argument_parser.cpp
    src/array_init.cpp
    src/arrays.cpp
    src/dimensions.cpp
    src/fdtd_grid_initialiser.cpp
    src/grid_labels.cpp
    src/interface.cpp
    src/iterator.cpp
    src/interpolation_methods.cpp
    src/matrix_collection.cpp
    src/mesh_base.cpp
    src/numerical_derivative.cpp
    src/shapes.cpp
    src/simulation_parameters.cpp
    src/source.cpp
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
            tests/unit/test_fields.cpp
            tests/unit/test_BLi_vs_cubic_interpolation.cpp
            tests/unit/test_field_interpolation.cpp
            tests/unit/test_interpolation_determination.cpp
            tests/unit/test_interpolation_functions.cpp
            tests/unit/test_numerical_derivative.cpp
            tests/unit/test_openandorder.cpp)

    target_link_libraries(tdms_tests PRIVATE
            Catch2::Catch2WithMain
            tdms_lib
            )

    add_test(test_all tdms_tests)
endif()
