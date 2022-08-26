#pragma once
#include <string>
#include <stdexcept>
#include "matlabio.h"


class XYZVectors {
public:
  double* x = nullptr;
  double* y = nullptr;
  double* z = nullptr;

  /**
   * Default constructor
   */
  XYZVectors() = default;

  /**
   * Set the pointer for one of the vectors in this collection with a name of c
   * @param c Character labeling the vector
   * @param ptr Pointer to assign
   */
  void set_ptr(char c, double* ptr);
};

class MaterialCollection{
protected:
  static void init_xyz_vectors(const mxArray *ptr, XYZVectors &arrays, const std::string &prefix);
};

class CCollectionBase {
public:
  XYZVectors a;
  XYZVectors b;
  XYZVectors c;
};

class CCollection : public CCollectionBase {
private:
  void init_xyz_vectors(const mxArray *ptr, XYZVectors &arrays, const std::string &prefix);

public:
  bool is_multilayer = false;
  bool is_disp_ml = false;

  explicit CCollection(const mxArray *ptr);
};

class CMaterial : public CCollectionBase, MaterialCollection {
public:
  explicit CMaterial(const mxArray *ptr);
};

class DCollectionBase {
public:
  XYZVectors a;
  XYZVectors b;
};

class DCollection: public DCollectionBase{
private:
  static void init_xyz_vectors(const mxArray *ptr, XYZVectors &arrays, const std::string &prefix);

public:
  explicit DCollection(const mxArray *ptr);
};

class DMaterial : public DCollectionBase, MaterialCollection {
public:
  explicit DMaterial(const mxArray *ptr);
};

class DispersiveMultiLayer{
public:
  double* alpha = nullptr;
  double* beta = nullptr;
  double* gamma = nullptr;
  XYZVectors kappa;
  XYZVectors sigma;

public:
  explicit DispersiveMultiLayer(const mxArray *ptr);
};

class GratingStructure{
private:
  int** matrix = nullptr;

public:
  inline int* operator[] (int value) const { return matrix[value]; }

  GratingStructure(const mxArray *ptr, int I_tot);

  bool has_elements(){ return matrix != nullptr; };

  ~GratingStructure();
};

class FrequencyExtractVector{
private:
  double* vector = nullptr; // Internal array
  int n = 0;               // Number of elements

public:

  FrequencyExtractVector(const mxArray *ptr, double omega_an);

  inline double operator[] (int value) const { return vector[value]; }

  inline int size() const{ return n; }
};
