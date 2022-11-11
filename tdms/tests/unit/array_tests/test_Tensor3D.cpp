/**
 * @file test_Tensor3D.cpp
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Tests the functionality of the Tensor3D class, which is the building block for several further field classes
 *
 */
#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

#include "arrays.h"

const double tol = 1e-16;

TEST_CASE("Tensor3D") {
    SPDLOG_INFO("== Tensor3D");
    // some mock-up array dimensions
    int n_layers = 4, n_cols = 8, n_rows = 16;

    // default constructor should assign all dimensions to 0, and the tensor itself should be a nullptr
    Tensor3D<double> *default_constructed_tensor = new Tensor3D<double>;
    // should have no elements
    REQUIRE(!(default_constructed_tensor->has_elements()));
    // let us assign some free memory to this tensor via allocate()
    default_constructed_tensor->allocate(n_layers, n_cols, n_rows);
    // we should now "have elements", even though they are unassigned
    REQUIRE(default_constructed_tensor->has_elements());
    // now we can try to zero all the elements in the array, the Frobenius norm of the array should be zero if this was sucessful
    default_constructed_tensor->zero();
    REQUIRE(abs(default_constructed_tensor->frobenius()) < tol);

    // we should also check what happens when we use the overloaded constructor, providing a pre-built pointer to a "Tensor" in memory
    int ***p = (int ***) malloc(n_layers * sizeof(int **));
    for (int k = 0; k < n_layers; k++) {
        p[k] = (int **) malloc(n_cols * sizeof(int *));
        for (int j = 0; j < n_cols; j++) {
            p[k][j] = (int *) malloc(n_rows * sizeof(int));
        }
    }
    // assign some values to this tensor. We'll go with =0 if i+j+k is even, and =1 if odd
    // this gives us 4*8*16/2 = 4^4 = 256 non-zero entries, which are 1, so the analytic norm is 16.
    for (int k=0; k<n_layers; k++) {
        for (int j=0; j<n_cols; j++) {
            for (int i=0; i<n_rows; i++) {
                p[k][j][i] = (i+j+k) % 2;
            }
        }
    }
    // the analytic frobenuis norm
    double target_fro = 16.;
    // now prepare another tensor
    Tensor3D<int> *overloaded_constructed_tensor = new Tensor3D(p, n_layers, n_cols, n_rows);
    // this tensor should be flagged as "having elements", since we provided a pointer in the constructor
    REQUIRE(overloaded_constructed_tensor->has_elements());
    // it should also have the target Frobenius norm, since we preallocated our array via an int***
    REQUIRE(abs(overloaded_constructed_tensor->frobenius()-target_fro)<tol);

    // destruct both Tensors and free memory
    delete default_constructed_tensor;
    delete overloaded_constructed_tensor;
}
