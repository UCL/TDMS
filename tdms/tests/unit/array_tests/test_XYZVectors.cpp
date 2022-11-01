/**
 * @file test_XYZVectors.cpp
 * @author William Graham
 * @brief Tests allocation, deallocation, and methods for the XYZVectors class
 */
#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>
#include <numeric>

#include "arrays.h"
#include "globals.h"

using std::accumulate;

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
  // now let's try assigning
  v.set_ptr(AxialDirection::X, x_vec);
  REQUIRE(v.x != nullptr);
  v.set_ptr(AxialDirection::Y, y_vec);
  REQUIRE(v.y != nullptr);
  v.set_ptr(AxialDirection::Z, z_vec);
  REQUIRE(v.z != nullptr);
  // and manipulating the components
  double x_tot = accumulate(v.x, v.x + nx, 0.);
  CHECK(abs(x_tot - 0.) <= tol);
  double y_tot = accumulate(v.y, v.y + ny, 0.);
  CHECK(abs(y_tot - 3.) <= tol);
  double z_tot = accumulate(v.z, v.z + nz, 0.);
  CHECK(abs(z_tot - 10.) <= tol);

  // in theory, we can also swap the components by reassigning the pointers
  v.set_ptr('x', z_vec);
  x_tot = accumulate(v.x, v.x + nz, 0.);
  CHECK(abs(x_tot - 10.) <= tol);
}
