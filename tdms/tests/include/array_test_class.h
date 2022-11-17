/**
 * @file array_test_class.h
 * @brief Defines the test classes for the objects in arrays.h
 */
#pragma once

/**
 * @brief Abstract container for the tests that classes in arrays.h will have to pass.
 *
 * Private methods are empty, and designed to be overridden in the derived class. However this also means that definitions for these tests can be left out of those subclasses, whilst not changing the run_all_tests() function.
 *
 */
class AbstractArrayTest {
private:
  int empty_dimensions[2] = {0, 1};            //< For initialising empty arrays
  int I_tot = 4, J_tot = 8, K_tot = 4;         //< For mimicing simulation grid dimensions
  int dimensions_2d[2] = {1, 1};               //< For defining 2D MATLAB arrays
  int dimensions_3d[3] = {I_tot, J_tot, K_tot};//< For defining 3D MATLAB arrays
  // we could put some common values, EG empty_dimensions, I_tot, etc here as private values to save redefining everything

  void test_empty_construction(){};//< Tests the behaviour of construction when passed an empty MATLAB array
  void test_correct_construction(){};//< Tests the behaviour of construction methods when passing the expected inputs

public:
  AbstractArrayTest() = default;

  // Runs all the class's private test methods sequentially
  void run_all_class_tests() {
    test_empty_construction();
    test_correct_construction();
  }
};
