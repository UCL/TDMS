#pragma once
#include <complex>
#include <string>
#include <stdexcept>
#include "matlabio.h"
#include "numeric.h"
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

// TODO: docstring
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

// TODO: docstring
class CCollection : public CCollectionBase {
private:
  void init_xyz_vectors(const mxArray *ptr, XYZVectors &arrays, const std::string &prefix);

public:
  bool is_multilayer = false;
  bool is_disp_ml = false;

  explicit CCollection(const mxArray *ptr);
};

// TODO: docstring
class CMaterial : public CCollectionBase, MaterialCollection {
public:
  explicit CMaterial(const mxArray *ptr);
};

class DCollectionBase {
public:
  XYZVectors a;
  XYZVectors b;
};

// TODO: docstring
class DCollection: public DCollectionBase{
private:
  static void init_xyz_vectors(const mxArray *ptr, XYZVectors &arrays, const std::string &prefix);

public:
  explicit DCollection(const mxArray *ptr);
};

// TODO: docstring
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
  int n_rows = 0;
  int n_cols = 0;
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
  int n = 0;        // Number of elements
  T* vector = nullptr; // Internal array

public:
  Vector() = default;

  explicit Vector(const mxArray *ptr);

  inline T operator[] (int value) const { return vector[value]; };

  inline int size() const{ return n; };
};

class FrequencyExtractVector: public Vector<double>{
public:
  FrequencyExtractVector(const mxArray *ptr, double omega_an);
};

class FrequencyVectors{
public:
  Vector<double> x;
  Vector<double> y;

  void initialise(const mxArray *ptr);
};

// TODO: docstring
class Pupil: public Matrix<double>{
public:
  Pupil() = default;

  void initialise(const mxArray *ptr, int n_rows, int n_cols);

  ~Pupil();
};

template<typename T>
class Tensor3D{
protected:
  int n_layers = 0;
  int n_cols = 0;
  int n_rows = 0;
  T*** tensor = nullptr;

public:
  bool is_matlab_initialised = false;

  Tensor3D() = default;

  Tensor3D(T*** tensor, int n_layers, int n_cols, int n_rows);

  inline T** operator[] (int value) const { return tensor[value]; };

  bool has_elements(){ return tensor != nullptr; };

  ~Tensor3D(){
    if (tensor != nullptr){
      if (is_matlab_initialised){
        free_cast_matlab_3D_array(tensor, n_layers);
      } else {
        destroy_3D_array(&tensor, n_cols, n_layers);
      }
    }
  };
};

class DTilde{
protected:
  int n_det_modes = 0;
  static Tensor3D<std::complex<double>> component_in(const mxArray *ptr, const std::string &name,
                                                     int n_rows, int n_cols);
public:
  inline int num_det_modes() const { return n_det_modes; };

  Tensor3D<std::complex<double>> x;
  Tensor3D<std::complex<double>> y;

  void initialise(const mxArray *ptr, int n_rows, int n_cols);
};

class IncidentField{
protected:
  static Tensor3D<double> component_in(const mxArray *ptr, const std::string &name);

public:
  Tensor3D<double> x;
  Tensor3D<double> y;

  explicit IncidentField(const mxArray *ptr);
};
