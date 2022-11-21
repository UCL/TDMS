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
  int n_rows = 4, n_cols = 8;

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
  std::string get_class_name() override { return "FrequencyVectors"; }
};

class IncidentFieldTest : public AbstractArrayTest {
private:
  const int n_rows = 6, n_cols = 4, n_layers = 5;
  const int n_fields = 2;
  const char *fieldnames[2] = {"exi", "eyi"};

  void test_empty_construction() override;
  void test_wrong_input_type() override;
  void test_incorrect_number_of_fields() override;
  void test_correct_construction() override;

public:
  std::string get_class_name() override { return "IncidentField"; }
};

class CMaterialTest : public AbstractArrayTest {
private:
  const int n_fields = 9;
  const char *fieldnames[9] = {"Cax", "Cay", "Caz", "Cbx", "Cby", "Cbz", "Ccx", "Ccy", "Ccz"};
  const char *wrong_fieldnames[9] = {"Dax", "Cay", "Daz", "Cbx", "Dby", "Cbz", "Ccx", "Ccy", "Ccz"};

  void test_incorrect_number_of_fields() override;
  void test_incorrect_fieldname() override;
  void test_correct_construction() override;

public:
  std::string get_class_name() override { return "CMaterial"; }
};

class DMaterialTest : public AbstractArrayTest {
private:
  const int n_fields = 6;
  const char *fieldnames[6] = {"Dax", "Day", "Daz", "Dbx", "Dby", "Dbz"};
  const char *wrong_fieldnames[6] = {"Dax", "Cay", "Daz", "Cbx", "Dby", "Cbz"};

  void test_incorrect_number_of_fields() override;
  void test_incorrect_fieldname() override;
  void test_correct_construction() override;

public:
  std::string get_class_name() override { return "DMaterial"; }
};

class MatrixTest : public AbstractArrayTest {
private:
  const int n_rows = 4, n_cols = 8;

  void test_correct_construction() override;
  // allocate() must be tested
  void test_other_methods() override;

public:
  std::string get_class_name() override { return "Matrix"; }
};

class VerticesTest : public AbstractArrayTest {
private:
  const int n_fields = 1;
  const char *fieldnames[1] = {"vertices"};

  void test_correct_construction() override;

public:
  std::string get_class_name() override { return "Vertices"; }
};

class GratingStructureTest : public AbstractArrayTest {
private:
  void test_empty_construction() override;
  void test_wrong_input_dimensions() override;
  void test_correct_construction() override;

public:
  std::string get_class_name() override { return "GratingStructure"; }
};

class PupilTest : public AbstractArrayTest {
private:
  const int n_rows = 4, n_cols = 8;

  void test_empty_construction() override;
  void test_wrong_input_dimensions() override;
  void test_correct_construction() override;

public:
  std::string get_class_name() override { return "Pupil"; }
};

class EHVecTest : public AbstractArrayTest {
private:
  const int n_rows = 4, n_cols = 8;

  // allocate() needs testing
  void test_other_methods() override;

public:
  std::string get_class_name() override { return "EHVec"; }
};

class Tensor3DTest : public AbstractArrayTest {
private:
  const int n_layers = 4, n_cols = 8, n_rows = 16;

  void test_correct_construction() override;
  // allocate() (and) zero(), frobenius()
  void test_other_methods() override;

public:
  std::string get_class_name() override { return "Tensor3D"; }
};

/*
class : public AbstractArrayTest {
private:
  void test_empty_construction() override;
  void test_wrong_input_dimensions() override;
  void test_wrong_input_type() override;
  void test_incorrect_number_of_fields() override;
  void test_incorrect_fieldname() override;
  void test_correct_construction() override;
  void test_initialise_method() override;
  void test_other_methods() override;

public:
  std::string get_class_name() override { return "SETME"; }
};
*/
