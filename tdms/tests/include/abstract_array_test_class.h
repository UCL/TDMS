/**
 * @file array_test_class.h
 * @brief Defines the abstract class for testing objects defined in arrays.h
 */
#pragma once

#include <catch2/catch_test_macros.hpp>
#include <spdlog/spdlog.h>
#include <string>

#include "arrays.h"

/**
 * @brief Abstract container for the tests that classes in arrays.h will have to pass.
 *
 * Private methods are empty, and designed to be overridden in the derived class. However this also means that definitions for these tests can be left out of those subclasses, whilst not changing the run_all_tests() function.
 *
 */
class AbstractArrayTest {
protected:
  int empty_dimensions[2] = {0, 1};   //< For initialising empty arrays
  int I_tot = 4, J_tot = 8, K_tot = 4;//< For mimicing simulation grid dimensions
  int n_numeric_elements =
          8;//< For standardising the number of elements we initialise MATLAB vectors with
  int dimensions_2d[2] = {1, 1};               //< For defining 2D MATLAB arrays
  int dimensions_3d[3] = {I_tot, J_tot, K_tot};//< For defining 3D MATLAB arrays

  // To point to the (top-level) matlab inputs the classes require. Top-level structure ensures mxDestroyArray clears all sub-arrays too
  mxArray *matlab_input = nullptr;
  // Flags if we have created an mxArray manually, so we are forced to destroy to free memory when done
  bool matlab_input_assigned = false;

  /**
   * @brief Destroys the matlab_input array if it is assigned.
   *
   * @returns true If mxDestroyArray was called successfully
   * @returns false matlab_input was not assigned, so no deletion was necessary
   */
  bool destroy_matlab_input() {
    if (matlab_input_assigned) {
      mxDestroyArray(matlab_input);
      matlab_input_assigned = false;
      return true;
    } else {
      return false;
    }
  }
  /**
   * @brief Creates a MATLAB structure array
   *
   * @param n_rows,n_cols Dimensions of the array. This is the number of identical structs in the array, not the size of a particular field
   * @param n_fields Number of fields
   * @param fields Names of the fields
   */
  void create_struct_array(int n_rows, int n_cols, int n_fields, const char *fields[]) {
    matlab_input = mxCreateStructMatrix(n_rows, n_cols, n_fields, fields);
    matlab_input_assigned = true;
  }
  /**
   * @brief Creates a MATLAB structure array
   *
   * @param n_dimensions Number of dimensions of the array.
   * @param dimensions Array or vector, containing the number of elements in each axis of the array
   * @param n_fields Number of fields
   * @param fields Names of the fields
   */
  void create_struct_array(int n_dimensions, int *dimensions, int n_fields, const char *fields[]) {
    matlab_input = mxCreateStructArray(n_dimensions, (const mwSize *) dimensions, n_fields, fields);
    matlab_input_assigned = true;
  }
  /**
   * @brief Creates a 1-by-1 MATLAB structure array
   *
   * @param n_fields Number of fields
   * @param fields Names of the fields
   */
  void create_1by1_struct(int n_fields, const char *fields[]) {
    create_struct_array(1, 1, n_fields, fields);
  }
  /**
   * @brief Creates the MATLAB empty struct (no fields, 0 size)
   *
   */
  void create_empty_struct() {
    create_struct_array(0, 1, 0, {});
  }
  /**
   * @brief Creates a MATLAB numeric array
   *
   * @param n_dims Number of dimensions in the array
   * @param dimensions Array/vector containing the number of elements in each dimension, sequentially
   * @param data_type The MATLAB data type to populate the array with
   * @param complex_flag Whether the data is real or complex
   */
  void create_numeric_array(int n_dims, int *dimensions, mxClassID data_type = mxDOUBLE_CLASS,
                            mxComplexity complex_flag = mxREAL) {
    matlab_input =
            mxCreateNumericArray(n_dims, (const mwSize *) dimensions, data_type, complex_flag);
    matlab_input_assigned = true;
  }

  /* Test functions to be overriden by each class.
  These methods are designed to be overriden by the subclasses specific to each arrays.h class.

  The default behaviour of each test_ method is to set ran_a_test to false - this means that subclasses that skip over certain functionality (EG do not need to test the result of providing an empty MATLAB input) are skipped in when run_all_class_tests() is run, and do not bloat the log.
  If overriden, a test_ method should not alter ran_a_test
  */

  bool ran_a_test = true;
  /**
  * @brief Tests the behaviour of construction when passed an empty MATLAB array
  */
  virtual void test_empty_construction() { ran_a_test = false; }
  /**
   * @brief Tests the behaviour of construction when passed an array of the incorrect dimensions
   */
  virtual void test_wrong_input_dimensions() { ran_a_test = false; }
  /**
   * @brief Tests the behaviour of construction when passed an array of the incorrect type
   */
  virtual void test_wrong_input_type() { ran_a_test = false; }
  /**
   * @brief Tests the behaviour of construction when passed a struct with the wrong number of fields
   */
  virtual void test_incorrect_number_of_fields() { ran_a_test = false; }
  /**
   * @brief Tests the behaviour of construction methods when passed a struct with an incorrect fieldname
   */
  virtual void test_incorrect_fieldname() { ran_a_test = false; }
  /**
   * @brief Tests the behaviour of construction methods when passing the expected inputs
   */
  virtual void test_correct_construction() { ran_a_test = false; }
  /**
   * @brief Tests the behaviour of the initialise method
   */
  virtual void test_initialise_method() { ran_a_test = false; }
  /**
   * @brief Tests any other methods that the class may have
   */
  virtual void test_other_methods() { ran_a_test = false; }

public:
  AbstractArrayTest() = default;

  /**
   * @brief Returns the name of the class, for logging purposes
   */
  virtual std::string get_class_name() { return ""; }

  /**
   * @brief Runs all the unit tests associated to this class
   */
  void run_all_class_tests() {
    // String to print to the log
    std::string logging_string = get_class_name() + ": ";

    // run all tests
    SECTION("Empty construction") {
      logging_string += "Empty construction";
      test_empty_construction();
    }
    SECTION("Incorrect dimension of input") {
      logging_string += "Incorrect dimension of input";
      test_wrong_input_dimensions();
    }
    SECTION("Incorrect type of input") {
      logging_string += "Incorrect type of input";
      test_wrong_input_type();
    }
    SECTION("Incorrect number of fields") {
      logging_string += "Incorrect number of fields";
      test_incorrect_number_of_fields();
    }
    SECTION("Incorrect fieldname") {
      logging_string += "Incorrect fieldname";
      test_incorrect_fieldname();
    }
    SECTION("Correct construction") {
      logging_string += "Correct construction";
      test_correct_construction();
    }
    SECTION("initialise() method") {
      logging_string += "initialise() method";
      test_initialise_method();
    }
    SECTION("Testing class-specific methods") {
      logging_string += "class specific method";
      test_other_methods();
    }

    // suppress output if no test was run
    if (ran_a_test) {
      // tear down is necessary if a test was run
      bool needed_to_destroy = destroy_matlab_input();
      if (needed_to_destroy) {
        logging_string += " (matlab_input destroyed)";
      } else {
        logging_string += " (nothing to tear down)";
      }
      SPDLOG_INFO(logging_string);
    }
  }
};
