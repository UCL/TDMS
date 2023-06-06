#pragma once

#include <complex>
#include <vector>

#include <H5Cpp.h>

#include "globals.h"

/** @brief Functionality that will allow us to read complex data from .mat
 * files, via the HDF5 C++ API. */
namespace hdf5_MATLAB_complex {
/**
 * @brief Structure array matching the compound datatype that MATLAB .mat (v7.3)
 * files are saved as.
 * @details These data-structures are always interleaved float64s (C++ native
 * doubles). The datatype itself has two components, "real" and "imag", which
 * need corresponding members in the compound datatype we define for HDF5 to
 * use.
 */
struct MATLAB_complex {
  double real;//<! Read part
  double imag;//<! Imaginary part

  /**
   * @brief Cast to std::complex<double>
   * @details Utility function, as this allows us to cast to
   * std::complex<double> so that we can actually do mathematical operations and
   * use C++ stdlib functions.
   */
  operator std::complex<double>() const {
    return real + imag * tdms_math_constants::IMAGINARY_UNIT;
  }
};

/**
 * @brief Returns the compound type that HDF5 needs in order to read complex
 * values from .mat files.
 * @details We cannot typedef a compound type, because they must be instantiated
 * and have members assigned. Therefore, this function performs the necessary
 * steps to assemble the compound datatype, and return it so it can be passed to
 * the HDF5 API.
 * @return H5::CompType Compound type representing the MATLAB_complex structure
 */
H5::CompType to_hdf5_CompType();

/** @brief Converts a vector of MATLAB_complex values to a vector of
 * std::complex<doubles>. */
std::vector<std::complex<double>>
to_complex_vector(const std::vector<MATLAB_complex> &matlab_complex);
}// namespace hdf5_MATLAB_complex
