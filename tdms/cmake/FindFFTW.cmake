# - Find the FFTW library
#
# Adapted from: https://github.com/egpbos/findFFTW/blob/master/FindFFTW.cmake
#
# Usage:
#   find_package(FFTW [REQUIRED])
#
# It sets the following variables:
#   FFTW_FOUND                  ... true if fftw is found on the system
#   FFTW_LIBRARIES              ... full paths to all found fftw libraries
#   FFTW_[component]_LIB        ... full path to one of the components
#   FFTW_INCLUDE_DIRS           ... fftw include directory paths
#
# The following variables will be checked by the function
#   FFTW_USE_STATIC_LIBS        ... if true, only static libraries are found, otherwise both static and shared.
#   FFTW_ROOT                   ... if set, the libraries are exclusively searched
#                                   under this path


if( NOT FFTW_ROOT AND DEFINED ENV{FFTWDIR} )
    set( FFTW_ROOT $ENV{FFTWDIR} )
endif()

# Check if we can use PkgConfig
find_package(PkgConfig)
if( PKG_CONFIG_FOUND AND NOT FFTW_ROOT )
    pkg_check_modules( PKG_FFTW QUIET "fftw3" )
endif()

#Check whether to search static or dynamic libs
set( CMAKE_FIND_LIBRARY_SUFFIXES_SAV ${CMAKE_FIND_LIBRARY_SUFFIXES} )

if( ${FFTW_USE_STATIC_LIBS} )
    set( CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX} )
else()
    set( CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES_SAV} )
endif()

if()
    find_library(
            FFTW_DOUBLE_LIB
            NAMES "fftw3" libfftw3-3
            PATHS ${FFTW_ROOT}
            PATH_SUFFIXES "lib" "lib64"
            NO_DEFAULT_PATH
    )

    find_path(FFTW_INCLUDE_DIRS
            NAMES "fftw3.h"
            PATHS ${FFTW_ROOT}
            PATH_SUFFIXES "include"
            NO_DEFAULT_PATH
            )

else()

    find_library(
            FFTW_DOUBLE_LIB
            NAMES "fftw3"
            PATHS ${PKG_FFTW_LIBRARY_DIRS} ${LIB_INSTALL_DIR}
    )

    find_path(FFTW_INCLUDE_DIRS
            NAMES "fftw3.h"
            PATHS ${PKG_FFTW_INCLUDE_DIRS} ${INCLUDE_INSTALL_DIR}
            )

endif()

if (NOT FFTW_DOUBLE_LIB)
    message("Failed to find FFTW3. Set FFTW_ROOT with e.g.\n
             cmake .. -DFFTW_ROOT=/path/to/fftw3/\n")
endif()

set(FFTW_LIBRARIES ${FFTW_LIBRARIES} ${FFTW_DOUBLE_LIB})
add_library(FFTW::Double INTERFACE IMPORTED)
set_target_properties(FFTW::Double
        PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${FFTW_INCLUDE_DIRS}"
        INTERFACE_LINK_LIBRARIES "${FFTW_DOUBLE_LIB}"
        )

set( CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES_SAV} )

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(FFTW
        REQUIRED_VARS FFTW_INCLUDE_DIRS
        HANDLE_COMPONENTS
        )

mark_as_advanced(
        FFTW_INCLUDE_DIRS
        FFTW_LIBRARIES
        FFTW_DOUBLE_LIB
)