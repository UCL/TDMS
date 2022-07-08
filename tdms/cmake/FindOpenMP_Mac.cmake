# Find openMP specifically on a mac which does not ship with it
# Allows for the setting of OMP_ROOT which should be searched

find_library(OpenMP_CXX_LIBRARY
        NAMES omp
        PATHS ${OMP_ROOT}
        PATH_SUFFIXES "lib"
        NO_DEFAULT_PATH
        )

mark_as_advanced(OpenMP_CXX_LIBRARY)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OpenMP_Mac DEFAULT_MSG
        OpenMP_CXX_LIBRARY)

if (OpenMP_Mac_FOUND)
    set(OpenMP_CXX_LIBRARY ${OpenMP_CXX_LIBRARY})
    set(OpenMP_COMPILE_OPTIONS -Xpreprocessor -fopenmp)

    add_library(OpenMP::OpenMP_CXX SHARED IMPORTED)
    set_target_properties(OpenMP::OpenMP_CXX PROPERTIES
            IMPORTED_LOCATION "${OpenMP_CXX_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${OpenMP_INCLUDE_DIRS}"
            INTERFACE_COMPILE_OPTIONS "${OpenMP_COMPILE_OPTIONS}"
            )
endif()
