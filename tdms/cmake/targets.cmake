function(release_target)
    add_executable(tdms)

    target_sources(tdms PRIVATE ${SOURCES} src/openandorder.cpp)

    target_link_libraries(tdms
            FFTW::Double
            ${Matlab_MEX_LIBRARY}
            ${Matlab_MX_LIBRARY}
            ${Matlab_MAT_LIBRARY}
            ${LIBCXX_LIBRARY}
            OpenMP::OpenMP_CXX
            spdlog::spdlog
            )
endfunction()

function(test_target)

    if (WIN32)
        message(FATAL_ERROR Cannot build TDMS tests on Windows)
    endif(WIN32)

    # catch2 tests ------------------------------------------------------------
    Include(FetchContent)

    FetchContent_Declare(
            Catch2
            GIT_REPOSITORY ${GITHUB_PREFIX}catchorg/Catch2.git
            GIT_TAG        v3.0.1
    )

    FetchContent_MakeAvailable(Catch2)
    list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/contrib)
    include(CTest)
    include(Catch)

    add_library(tdms_lib SHARED)
    target_sources(tdms_lib PUBLIC ${SOURCES})
    add_executable(tdms "src/main.cpp")

    target_link_libraries(tdms_lib LINK_PUBLIC
            FFTW::Double
            ${Matlab_MEX_LIBRARY}
            ${Matlab_MX_LIBRARY}
            ${Matlab_MAT_LIBRARY}
            ${LIBCXX_LIBRARY}
            OpenMP::OpenMP_CXX
            spdlog::spdlog
            )

    target_link_libraries(tdms tdms_lib)
    target_compile_options(tdms_lib PUBLIC
            -c ${DFLAG}
            -O3
            -DSPDLOG_BUILD_SHARED=ON
            -DMX_COMPAT_32)

endfunction()
