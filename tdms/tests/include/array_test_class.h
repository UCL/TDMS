/**
 * @file array_test_class.h
 * @brief Defines the test classes for the objects in arrays.h
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
  int empty_dimensions[2] = {0, 1};            //< For initialising empty arrays
  int I_tot = 4, J_tot = 8, K_tot = 4;         //< For mimicing simulation grid dimensions
  int n_rows = 4, n_cols = 8;//< For mimicing matrix dimensions
  int n_numeric_elements = 8; //< For standardising the number of elements we initialise MATLAB vectors with
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
  void create_empty_struct(int n_fields, const char *fields[]) {
    create_struct_array(0, 1, n_fields, fields);
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

  The default behaviour of each test_ method is to return false - this means that subclasses that skip over certain functionality (EG do not need to test the result of providing an empty MATLAB input) are skipped in when run_all_class_tests() is run, and do not bloat the log.
  If overriden, a test_ method should return true, to indicate to run_all_class_tests() that this test did attempt to run (and might have failed).
  */

  /**
  * @brief Tests the behaviour of construction when passed an empty MATLAB array
  */
  virtual bool test_empty_construction() { return false;}
  /**
   * @brief Tests the behaviour of construction when passed an array of the incorrect dimensions
   */
  virtual bool test_wrong_input_dimensions() { return false;}
  /**
   * @brief Tests the behaviour of construction when passed a struct with the wrong number of fields
   */
  virtual bool test_incorrect_number_of_fields() { return false;}
  /**
   * @brief Tests the behaviour of construction methods when passing the expected inputs
   */
  virtual bool test_correct_construction() { return false;}
  /**
   * @brief Tests the behaviour of the initialise method
   */
  virtual bool test_initialise_method() { return false;}

public:
  AbstractArrayTest() = default;

  /**
   * @brief Returns the name of the class, for logging purposes
   */
  virtual std::string get_class_name() { return "";}

  /**
   * @brief Runs all the unit tests associated to this class
   */
  void run_all_class_tests() {
    // String to print to the log
    std::string logging_string = get_class_name() + ": ";
    // Skips over tests that do not need to be executed
    bool test_executed;

    // run all tests
    SECTION("Empty construction") {
      logging_string += "Empty construction";
      test_executed = test_empty_construction();
    }
    SECTION("Incorrect dimension of input") {
      logging_string += "Incorrect dimension of input";
      test_executed = test_wrong_input_dimensions();
    }
    SECTION("Incorrect number of fields") {
      logging_string += "Incorrect number of fields";
      test_executed = test_incorrect_number_of_fields();
    }
    SECTION("Correct construction") {
      logging_string += "Correct construction";
      test_executed = test_correct_construction();
    }

    // suppress output if no test was run
    if (test_executed) {
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

class CCollectionTest : public AbstractArrayTest {
  private:
    const char *fieldnames[9] = { "Cax", "Cay", "Caz", "Cbx", "Cby", "Cbz", "Ccx", "Ccy", "Ccz" };

    bool test_incorrect_number_of_fields() override;
    bool test_correct_construction() override;

  public:
    std::string get_class_name() override {return "CCollection";}
};

class DCollectionTest : public AbstractArrayTest {
  private:
    const char *fieldnames[6] = {"Dax", "Day", "Daz", "Dbx", "Dby", "Dbz"};

    bool test_incorrect_number_of_fields() override;
    bool test_correct_construction() override;

  public:
    std::string get_class_name() override { return "DCollection"; }
};

class ComplexAmplitudeSampleTest : public AbstractArrayTest {
private:
  const char *fieldnames[2] = {"vertices", "components"};

  bool test_empty_construction() override;
  bool test_correct_construction() override;

public:
  std::string get_class_name() override { return "ComplexAmplitudeSample"; }
};

class DetectorSensitivityArraysTest : public AbstractArrayTest {
private:
  bool test_correct_construction() override;
  bool test_initialise_method() override;

public:
  std::string get_class_name() override { return "DetectorSensitivityArray";}
};
