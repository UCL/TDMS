#pragma once
#include <string>
#include <stdexcept>
#include "matlabio.h"
#include "utils.h"


class XYZTensor3D {
public:
  double ***x = nullptr;
  double ***y = nullptr;
  double ***z = nullptr;

  double*** operator[] (char c) const{
    switch (c) {
      case 'x': return x;
      case 'y': return y;
      case 'z': return z;
      default: throw std::runtime_error("Have no element " + to_string(c));
    }
  }
};

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

template<typename T>
class Matrix{
protected:
  size_t n_rows = 0;
  size_t n_cols = 0;
  T** matrix = nullptr;

public:
  inline T* operator[] (int value) const { return matrix[value]; }

  bool has_elements(){ return matrix != nullptr; };
};

class GratingStructure: public Matrix<int>{

public:
  GratingStructure(const mxArray *ptr, int I_tot);

  ~GratingStructure();
};

template<typename T>
class Vector{
protected:
  size_t n = 0;        // Number of elements
  T* vector = nullptr; // Internal array

public:
  Vector() = default;

  explicit Vector(const mxArray *ptr);

  inline T operator[] (int value) const { return vector[value]; };

  inline size_t size() const{ return n; };
};

class FrequencyExtractVector: public Vector<double>{
public:
  FrequencyExtractVector(const mxArray *ptr, double omega_an);
};

class FrequencyVectors: public Vector<double>{
public:
  Vector x;
  Vector y;

  void initialise(const mxArray *ptr);
};
