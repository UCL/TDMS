#include <complex>
#include <vector>

#include <H5Cpp.h>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

#include "globals.h"
#include "hdf5_io/hdf5_compound_types.h"
#include "hdf5_io/hdf5_reader.h"
#include "unit_test_utils.h"

using namespace std;
using namespace hdf5_MATLAB_complex;
using Catch::Approx;
using tdms_math_constants::IMAGINARY_UNIT;
using tdms_unit_test_data::struct_testdata;

TEST_CASE("Compound Types: MATLAB_complex") {
  spdlog::info("Compound types: MATLAB_complex");
  HDF5Reader MATFile(struct_testdata);

  SECTION("Check casting to complex<double>") {
    // Cast a standalone MATLAB_complex
    MATLAB_complex c = {3.14156, 2.7182818};
    complex<double> c_via_cast = c;
    REQUIRE(c.real == Approx(c_via_cast.real()));
    REQUIRE(c.imag == Approx(c_via_cast.imag()));

    // Now check the vector conversion function
    const int n_vec_elements = 5;
    vector<MATLAB_complex> c_vec(n_vec_elements);
    vector<complex<double>> c_vec_expected(n_vec_elements);
    // Manually set the expected values
    for (unsigned int i = 0; i < n_vec_elements; i++) {
      c_vec[i] = {i + 0.5, n_vec_elements - 0.5};
      c_vec_expected[i] = {i + 0.5, n_vec_elements - 0.5};
    }
    // Convert via function and check we got the expected values
    vector<complex<double>> c_vec_via_cast = to_complex_vector(c_vec);
    bool all_elements_correct = true;
    for (unsigned int i = 0; i < n_vec_elements; i++) {
      all_elements_correct =
              all_elements_correct &&
              abs(c_vec_via_cast[i] - c_vec_expected[i]) == Approx(0.);
    }
    REQUIRE(all_elements_correct);
  }

  SECTION("Read MATLAB complex array") {
    vector<MATLAB_complex> data_vector;
    MATFile.read("example_struct/complex_22", data_vector, to_hdf5_CompType());

    // Should have read the 2x2 Pauli matrix
    REQUIRE(data_vector.size() == 4);
    // Check that the values are as expected, cast to std::complex types for
    // ease of comparison, previous SECTION checks that this is OK to do.
    vector<complex<double>> read_data = to_complex_vector(data_vector);
    vector<complex<double>> pauli_matrix(4);
    pauli_matrix[0] = 0.;
    pauli_matrix[1] = IMAGINARY_UNIT;
    pauli_matrix[2] = -IMAGINARY_UNIT;
    pauli_matrix[3] = 0.;
    bool all_elements_correct = true;
    for (unsigned int i = 0; i < 4; i++) {
      all_elements_correct = all_elements_correct &&
                             abs(read_data[i] - pauli_matrix[i]) == Approx(0.);
    }
    REQUIRE(all_elements_correct);
  }
}
