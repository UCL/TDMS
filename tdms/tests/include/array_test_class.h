/**
 * @file array_test_class.h
 * @brief Defines the test classes for the objects in arrays.h
 *
 *
 * TEMPLATE FOR MAKING A NEW TEST ARRAY CLASS
 * class : public AbstractArrayTest {
 * private:
 *   void test_empty_construction() override;
 *   void test_wrong_input_dimensions() override;
 *   void test_wrong_input_type() override;
 *   void test_incorrect_number_of_fields() override;
 *   void test_incorrect_fieldname() override;
 *   void test_correct_construction() override;
 *   void test_initialise_method() override;
 *   void test_other_methods() override;
 *
 * public:
 *   std::string get_class_name() override { return "SETME"; }
 * };
 */
#pragma once

#include "abstract_array_test_class.h"

/** @brief Unit tests for CCollection */
class CCollectionTest : public AbstractArrayTest {
private:
  const char *fieldnames[9] = {"Cax", "Cay", "Caz", "Cbx", "Cby",
                               "Cbz", "Ccx", "Ccy", "Ccz"};

  void test_incorrect_number_of_fields() override;
  void test_correct_construction() override;

public:
  std::string get_class_name() override { return "CCollection"; }
};

/** @brief Unit tests for DCollection */
class DCollectionTest : public AbstractArrayTest {
private:
  const char *fieldnames[6] = {"Dax", "Day", "Daz", "Dbx", "Dby", "Dbz"};

  void test_incorrect_number_of_fields() override;
  void test_correct_construction() override;

public:
  std::string get_class_name() override { return "DCollection"; }
};

/** @brief Unit tests for DetectorSensitivityArrays */
class DetectorSensitivityArraysTest : public AbstractArrayTest {
private:
  int n_rows = 4, n_cols = 8;

  void test_correct_construction() override;
  void test_initialise_method() override;

public:
  std::string get_class_name() override { return "DetectorSensitivityArray"; }
};

// Test methods check the performance of initialise, as this is the de-facto
// constructor
class DTildeTest : public AbstractArrayTest {
private:
  const int n_fields = 2;
  const char *fieldnames[2] = {"Dx_tilde", "Dy_tilde"};

  void test_empty_construction() override;
  void test_wrong_input_type() override;
  void test_incorrect_number_of_fields() override;
  void test_correct_construction() override;
  void test_initialise_method() override;

public:
  std::string get_class_name() override { return "DTilde"; }
};

/** @brief Unit tests for FieldSample */
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

/** @brief Unit tests for FullFieldSnapshot */
class FullFieldSnapshotTest : public AbstractArrayTest {
private:
  // multiply_{E,H}_by()
  void test_other_methods() override;

public:
  std::string get_class_name() override { return "FullFieldSnapshot"; }
};

/** @brief Unit tests for IncidentField */
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

/** @brief Unit tests for CMaterial */
class CMaterialTest : public AbstractArrayTest {
private:
  const int n_fields = 9;
  const char *fieldnames[9] = {"Cax", "Cay", "Caz", "Cbx", "Cby",
                               "Cbz", "Ccx", "Ccy", "Ccz"};
  const char *wrong_fieldnames[9] = {"Dax", "Cay", "Daz", "Cbx", "Dby",
                                     "Cbz", "Ccx", "Ccy", "Ccz"};

  void test_incorrect_number_of_fields() override;
  void test_incorrect_fieldname() override;
  void test_correct_construction() override;

public:
  std::string get_class_name() override { return "CMaterial"; }
};

/** @brief Unit tests for DMaterial */
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

/** @brief Unit tests for Matrix */
class MatrixTest : public AbstractArrayTest {
private:
  const int n_rows = 4, n_cols = 8;

  void test_correct_construction() override;
  // allocate() must be tested
  void test_other_methods() override;

public:
  std::string get_class_name() override { return "Matrix"; }
};

/** @brief Unit tests for Vertices */
class VerticesTest : public AbstractArrayTest {
private:
  const int n_fields = 1;
  const char *fieldnames[1] = {"vertices"};

  void test_correct_construction() override;

public:
  std::string get_class_name() override { return "Vertices"; }
};

/** @brief Unit tests for GratingStructure */
class GratingStructureTest : public AbstractArrayTest {
private:
  void test_empty_construction() override;
  void test_wrong_input_dimensions() override;
  void test_correct_construction() override;

public:
  std::string get_class_name() override { return "GratingStructure"; }
};

/** @brief Unit tests for Pupil */
class PupilTest : public AbstractArrayTest {
private:
  const int n_rows = 4, n_cols = 8;

  void test_empty_construction() override;
  void test_wrong_input_dimensions() override;
  void test_correct_construction() override;

public:
  std::string get_class_name() override { return "Pupil"; }
};

/** @brief Unit tests for EHVec */
class EHVecTest : public AbstractArrayTest {
private:
  const int n_rows = 4, n_cols = 8;

  // allocate() needs testing
  void test_other_methods() override;

public:
  std::string get_class_name() override { return "EHVec"; }
};

/** @brief Unit tests for Tensor3D */
class Tensor3DTest : public AbstractArrayTest {
private:
  const int n_layers = 4, n_cols = 8, n_rows = 16;

  void test_correct_construction() override;
  // allocate() (and) zero(), frobenius()
  void test_other_methods() override;

public:
  std::string get_class_name() override { return "Tensor3D"; }
};

/** @brief Unit tests for Vector */
class VectorTest : public AbstractArrayTest {
private:
  void test_correct_construction() override;

public:
  std::string get_class_name() override { return "Vector"; }
};

/** @brief Unit tests for FieldComponentsVector */
class FieldComponentsVectorTest : public AbstractArrayTest {
private:
  const int n_fields = 1;
  const char *fieldnames[1] = {"components"};

  void test_correct_construction() override;
  void test_initialise_method() override;

public:
  std::string get_class_name() override { return "FieldComponentsVector"; }
};

/** @brief Unit tests for FrequencyExtractVector */
class FrequencyExtractVectorTest : public AbstractArrayTest {
private:
  const double omega_an = 1.;

  void test_empty_construction() override;
  void test_wrong_input_dimensions() override;
  void test_correct_construction() override;
  // max() method
  void test_other_methods() override;

public:
  std::string get_class_name() override { return "FrequencyExtractVector"; }
};

/** @brief Unit tests for VertexPhasors */
class VertexPhasorsTest : public AbstractArrayTest {
private:
  const char *fieldnames[2] = {"vertices", "components"};

  void test_empty_construction() override;
  void test_correct_construction() override;

public:
  std::string get_class_name() override { return "VertexPhasors"; }
};

/** @brief Unit tests for XYZTensor3D */
class XYZTensor3DTest : public AbstractArrayTest {
private:
  const int n_layers = 4, n_cols = 8, n_rows = 16;

  void test_correct_construction() override;
  // allocate()
  void test_other_methods() override;

public:
  std::string get_class_name() override { return "XYZTensor3D"; }
};

/** @brief Unit tests for XYZVectors */
class XYZVectorsTest : public AbstractArrayTest {
private:
  const static int n_layers = 4;
  const static int n_cols = 8;
  const static int n_rows = 16;

  void test_correct_construction() override;
  // set_ptr(), all_elements_less_than()
  void test_other_methods() override;

public:
  std::string get_class_name() override { return "XYZVectors"; }
};
