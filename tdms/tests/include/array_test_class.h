/**
 * @file array_test_class.h
 * @brief Defines the test classes for the objects in arrays.h
 */
#pragma once

#include "abstract_array_test_class.h"

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

class DispersiveMultilayerTest : public AbstractArrayTest {
private:
  const int n_fields = 9;
  const char *fieldnames[9] = {"alpha",   "beta",    "gamma",   "kappa_x", "kappa_y",
                               "kappa_z", "sigma_x", "sigma_y", "sigma_z"};

  bool test_empty_construction() override;
  bool test_wrong_input_type() override;
  bool test_correct_construction() override;

public:
  std::string get_class_name() override { return "DispersiveMultilayerTest"; }
};

// DTilde test methods check the performance of initialise, as this is the de-facto constructor
class DTildeTest : public AbstractArrayTest {
private:
  const int n_fields = 2;
  const char *fieldnames[2] = { "Dx_tilde", "Dy_tilde" };

  bool test_empty_construction() override;
  bool test_wrong_input_type() override;
  bool test_incorrect_number_of_fields() override;
  bool test_correct_construction() override;
  bool test_initialise_method() override;

public:
  std::string get_class_name() override { return ""; }
};

/*
class : public AbstractArrayTest {
private:
  bool test_empty_construction() override;
  bool test_wrong_input_dimensions() override;
  bool test_wrong_input_type() override;
  bool test_incorrect_number_of_fields() override;
  bool test_correct_construction() override;
  bool test_initialise_method() override;

public:
  std::string get_class_name() override { return ""; }
};
*/
