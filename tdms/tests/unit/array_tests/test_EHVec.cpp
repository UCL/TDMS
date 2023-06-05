#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

#include "arrays/eh_vec.h"
#include "unit_test_utils.h"

using tdms_tests::is_close;

TEST_CASE("EHVec") {
  const int REAL = 0, IMAG = 1;
  const int n_rows = 4, n_cols = 8;
  EHVec eh;

  // allocate memory
  REQUIRE(!eh.is_allocated());
  eh.allocate(n_rows, n_cols);
  REQUIRE(eh.is_allocated());
  spdlog::info("Completed assignment");

  // check that we an assign fftw_complexes to the elements
  eh[0][0][REAL] = 1.;
  eh[0][0][IMAG] = 0.;
  eh[0][1][REAL] = 0.;
  eh[0][1][IMAG] = 1.;
  fftw_complex fftw_unit{1., 0.};
  fftw_complex fftw_imag_unit{0., 1.};

  bool elements_set_correctly =
          is_close(eh[0][0][REAL], fftw_unit[REAL]) &&
          is_close(eh[0][0][IMAG], fftw_unit[IMAG]) &&
          is_close(eh[0][1][REAL], fftw_imag_unit[REAL]) &&
          is_close(eh[0][1][IMAG], fftw_imag_unit[IMAG]);
  REQUIRE(elements_set_correctly);
}
