# General setup ---------------------------------------------------------------
cmake_minimum_required(VERSION 3.21)

project(tdms LANGUAGES CXX C)
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
include_directories(${CMAKE_SOURCE_DIR}/include/)
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

# hdf5 ------------------------------------------------------------------------
find_package(HDF5 REQUIRED COMPONENTS CXX)
include_directories(${HDF5_INCLUDE_DIR})

# spdlog ----------------------------------------------------------------------
find_package(spdlog NO_CMAKE_PACKAGE_REGISTRY QUIET)
if(NOT spdlog_FOUND)
    include(FetchContent)
    FetchContent_Declare(spdlog
	    GIT_REPOSITORY ${GITHUB_PREFIX}gabime/spdlog
            GIT_TAG v1.11.0
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

# Compile definintions
if (CMAKE_SOURCE_DIR)
    add_compile_definitions(CMAKE_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
    message(STATUS "Unit tests will use ${CMAKE_SOURCE_DIR} as a reference for finding data files.")
endif()


# TDMS version ----------------------------------------------------------------
# if supplied via the cmake configuration (-DTDMS_VERSION=whatever_I_want) then
# respect that, otherwise construct a version from a git tag or commit hash.
if (NOT TDMS_VERSION)
    include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/version_from_git.cmake)
endif()

if (TDMS_VERSION)
    add_compile_definitions(TDMS_VERSION="${TDMS_VERSION}")
    message(STATUS "Compiling TDMS version: ${TDMS_VERSION}")
endif()


# TDMS target -----------------------------------------------------------------

# wildcard to find all *.cpp files in the src directory except for main.cpp
# (which we don't want for the testing target)
file(GLOB_RECURSE SOURCES "src/*.cpp")
list(FILTER SOURCES EXCLUDE REGEX "src/main.cpp")

# TODO: can delete this line when matlabio is removed.
set(SOURCES ${SOURCES} matlabio/matlabio.cpp)

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


# Code Coverage ---------------------------------------------------------------
add_library(coverage_config INTERFACE)

option(CODE_COVERAGE "Enable coverage reporting" OFF)
if(CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  message("Setting up code coverage.")
  # Add required flags (GCC & LLVM/Clang)
  target_compile_options(coverage_config INTERFACE
    -O0        # no optimization
    -g         # generate debug info
    --coverage # sets all required flags
  )
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.13)
    target_link_options(coverage_config INTERFACE --coverage)
  else()
    target_link_libraries(coverage_config INTERFACE --coverage)
  endif()
endif(CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")


# Testing ---------------------------------------------------------------------
if (BUILD_TESTING)
    enable_testing()

    # find all .cpp files in the unit tests directory (and subdirectories).
    file(GLOB_RECURSE UNIT_TEST_SRC "${CMAKE_SOURCE_DIR}/tests/unit/test_*.cpp")

    add_executable(tdms_tests "${UNIT_TEST_SRC}")

    # add the additional include directories to the tdms_tests target
    target_include_directories(tdms_tests PRIVATE "${CMAKE_SOURCE_DIR}/tests/include")
    # link additional libraries for testing
    target_link_libraries(tdms_tests PRIVATE
            Catch2::Catch2WithMain
            tdms_lib
            coverage_config
            )

    add_test(test_all tdms_tests)
endif()
