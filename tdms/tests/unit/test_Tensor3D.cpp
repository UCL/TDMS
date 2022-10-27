/**
 * @file test_Tensor3D.cpp
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Tests the functionality of the Tensor3D class, which is the building block for several further field classes
 * 
 */
#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>

#include "field.h"

TEST_CASE("Tensor3D: allocation and deallocation") {
    // default constructor should assign all dimensions to 0, and the tensor itself should be a nullptr
    
}