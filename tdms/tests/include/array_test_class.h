/**
 * @file array_test_class.h
 * @brief Defines the test classes for the objects in arrays.h
 */
#pragma once

#include "abstract_array_test_class.h"

class CCollectionTest : public AbstractArrayTest {
  private:
    const char *fieldnames[9] = { "Cax", "Cay", "Caz", "Cbx", "Cby", "Cbz", "Ccx", "Ccy", "Ccz" };

    void test_incorrect_number_of_fields() override;
    void test_correct_construction() override;

  public:
    std::string get_class_name() override {return "CCollection";}
};

class DCollectionTest : public AbstractArrayTest {
  private:
    const char *fieldnames[6] = {"Dax", "Day", "Daz", "Dbx", "Dby", "Dbz"};

    void test_incorrect_number_of_fields() override;
    void test_correct_construction() override;

  public:
    std::string get_class_name() override { return "DCollection"; }
};

class ComplexAmplitudeSampleTest : public AbstractArrayTest {
private:
  const char *fieldnames[2] = {"vertices", "components"};

  void test_empty_construction() override;
  void test_correct_construction() override;

public:
  std::string get_class_name() override { return "ComplexAmplitudeSample"; }
};

class DetectorSensitivityArraysTest : public AbstractArrayTest {
private:
  void test_correct_construction() override;
  void test_initialise_method() override;

public:
  std::string get_class_name() override { return "DetectorSensitivityArray";}
};

class DispersiveMultilayerTest : public AbstractArrayTest {
private:
  const int n_fields = 9;
  const char *fieldnames[9] = {"alpha",   "beta",    "gamma",   "kappa_x", "kappa_y",
                               "kappa_z", "sigma_x", "sigma_y", "sigma_z"};

  void test_empty_construction() override;
  void test_wrong_input_type() override;
  void test_correct_construction() override;

public:
  std::string get_class_name() override { return "DispersiveMultilayer"; }
};

// Test methods check the performance of initialise, as this is the de-facto constructor
class DTildeTest : public AbstractArrayTest {
private:
  const int n_fields = 2;
  const char *fieldnames[2] = { "Dx_tilde", "Dy_tilde" };

  void test_empty_construction() override;
  void test_wrong_input_type() override;
  void test_incorrect_number_of_fields() override;
  void test_correct_construction() override;
  void test_initialise_method() override;

public:
  std::string get_class_name() override { return "DTilde"; }
};

class FieldSampleTest : public AbstractArrayTest {
private:
  const int n_fields = 4;
  const char *fieldnames[4] = {"i", "j", "k", "n"};

  void test_empty_construction() override;
  void test_wrong_input_type() override;
  void test_incorrect_number_of_fields() override;
  void test_correct_construction() override;

public:
  std::string get_class_name() override { return "FieldSample"; }
};

// Test methods check the performance of initialise, as this is the de-facto constructor
class FrequencyVectorsTest : public AbstractArrayTest {
private:
  const int n_fields = 2;
  const char *fieldnames[2] = {"fx_vec", "fy_vec"};

  void test_empty_construction() override;
  void test_wrong_input_type() override;
  void test_incorrect_number_of_fields() override;
  void test_correct_construction() override;
  void test_initialise_method() override;

public:
  std::string get_class_name() override { return ""; }
};
/*
class : public AbstractArrayTest {
private:
  void test_empty_construction() override;
  void test_wrong_input_dimensions() override;
  void test_wrong_input_type() override;
  void test_incorrect_number_of_fields() override;
  void test_correct_construction() override;
  void test_initialise_method() override;

public:
  std::string get_class_name() override { return ""; }
};
*/
