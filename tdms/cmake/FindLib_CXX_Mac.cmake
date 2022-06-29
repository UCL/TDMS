# Find libc++ on a Mac where a specific directory is searched
# required on an ARM mac targeting x86 and using a non-default compiler

if (${CXX_ROOT})
    find_library(LIBCXX_LIBRARY
            NAMES c++ cxx
            PATHS ${CXX_ROOT}
            PATH_SUFFIXES "lib"
            NO_DEFAULT_PATH
            )

    mark_as_advanced(LIBCXX_LIBRARY)
endif()
