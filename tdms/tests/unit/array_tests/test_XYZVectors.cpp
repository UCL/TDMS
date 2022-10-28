/**
 * @file test_XYZVectors.cpp
 * @author William Graham
 * @brief Tests allocation, deallocation, and methods for the XYZVectors class
 */
#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

#include "arrays.h"
#include "globals.h"

const double tol = 1e-8;

TEST_CASE("Test XYZVectors allocation") {
  SPDLOG_INFO("== Testing XYZVectors allocation/deallocation");
  // dimensions for the test-tensor
  const int nx = 4, ny = 8, nz = 16;

  // try creating a XYZVector
  XYZVectors v;
  // check that, although this has been declared, its members still point to nullptrs
  REQUIRE(v.x == nullptr);
  REQUIRE(v.y == nullptr);
  REQUIRE(v.z == nullptr);

  // create some arrays to assign to the members
  double x_vec[nx] = {0.}, y_vec[ny] = {2., 1.}, z_vec[nz] = {4., 3., 2., 1.};
  // and for testing they're pointing to the right values
  double x_tot = 0., y_tot = 0., z_tot = 0.;

  // now let's try assigning
  v.set_ptr(AxialDirection::X, x_vec);
  REQUIRE(v.x != nullptr);
  v.set_ptr(AxialDirection::Y, y_vec);
  REQUIRE(v.y != nullptr);
  v.set_ptr(AxialDirection::Z, z_vec);
  REQUIRE(v.z != nullptr);
  // and manipulating the components
  for (int ix = 0; ix < nx; ix++) { x_tot += v.x[ix]; }
  CHECK(abs(x_tot - 0.) <= tol);
  for (int iy = 0; iy < ny; iy++) { y_tot += v.y[iy]; }
  CHECK(abs(y_tot - 3.) <= tol);
  for (int iz = 0; iz < nz; iz++) { z_tot += v.z[iz]; }
  CHECK(abs(z_tot - 10.) <= tol);

  // in theory, we can also swap the components by reassigning the pointers
  v.set_ptr('x', z_vec);
  x_tot = 0.;
  for (int ix_but_really_iz = 0; ix_but_really_iz < nz; ix_but_really_iz++) {
    x_tot += v.x[ix_but_really_iz];
  }
  CHECK(abs(x_tot - 10.) <= tol);
}
