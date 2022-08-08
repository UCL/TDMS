function(release_target)
    add_executable(tdms)

    target_sources(tdms PRIVATE
            src/argument_parser.cpp
            src/fdtd_grid_initialiser.cpp
            src/openandorder.cpp
            src/interpolate.cpp
            src/numeric.cpp
            src/iterator.cpp
            src/mesh_base.cpp
            src/numerical_derivative.cpp
            src/utils.cpp
            src/matrix_collection.cpp
            matlabio/matlabio.cpp
            )

    target_link_libraries(tdms
            FFTW::Double
            ${Matlab_MEX_LIBRARY}
            ${Matlab_MX_LIBRARY}
            ${Matlab_MAT_LIBRARY}
            ${LIBCXX_LIBRARY}
            OpenMP::OpenMP_CXX
            )

    # Windows -----------------------------------------------------------------
    if (WIN32)
        # Windows requires the shared dlls to be in the same directory as the
        # target or to be present in the $PATH environment variable, MATLAB
        # requires a whole folder, which is wasteful to copy. Therefore, the
        # MATLABRootDIR/bin/win64 directory must be added to the $PATH
        # environment variable
        add_custom_command(TARGET tdms POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy ${FFTW_ROOT}/bin/fftw3.dll $<TARGET_FILE_DIR:tdms>
                COMMAND_EXPAND_LISTS
                )
    endif(WIN32)

endfunction()

function(test_target)

    if (WIN32)
        message(FATAL_ERROR Cannot build TDMS tests on Windows)
    endif(WIN32)

    # catch2 tests ----------------------------------------------------------------
    Include(FetchContent)

    FetchContent_Declare(
            Catch2
            GIT_REPOSITORY https://github.com/catchorg/Catch2.git
            GIT_TAG        v3.0.1
    )

    FetchContent_MakeAvailable(Catch2)
    list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/contrib)
    include(CTest)
    include(Catch)

    add_library(tdms_lib SHARED)
    target_sources(tdms_lib PUBLIC
            src/argument_parser.cpp
            src/fdtd_grid_initialiser.cpp
            src/interpolate.cpp
            src/iterator.cpp
            src/matrix_collection.cpp
            src/mesh_base.cpp
            src/numeric.cpp
            src/numerical_derivative.cpp
            src/utils.cpp
            matlabio/matlabio.cpp
            )

    add_executable(tdms "src/openandorder.cpp")

    target_link_libraries(tdms_lib LINK_PUBLIC
            FFTW::Double
            ${Matlab_MEX_LIBRARY}
            ${Matlab_MX_LIBRARY}
            ${Matlab_MAT_LIBRARY}
            ${LIBCXX_LIBRARY}
            OpenMP::OpenMP_CXX
            )

    target_link_libraries(tdms tdms_lib)
    target_compile_options(tdms_lib PUBLIC -DMX_COMPAT_32 -c ${DFLAG} -O3)

endfunction()

