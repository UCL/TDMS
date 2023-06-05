/**
 * @file arrays.h
 * @brief Classes describing arrays, vertices etc.
 */
#pragma once

#include <complex>
#include <stdexcept>
#include <string>
#include <vector>

#include <fftw3.h>

#include "arrays/tensor3d.h"
#include "globals.h"
#include "matlabio.h"
#include "utils.h"

template<typename T>
struct xyz_vector {
  std::vector<T> x = {};
  std::vector<T> y = {};
  std::vector<T> z = {};
};

template<typename T>
class XYZTensor3D {
public:
  T ***x = nullptr;
  T ***y = nullptr;
  T ***z = nullptr;

  T ***operator[](char c) const {
    switch (c) {
      case 'x':
        return x;
      case 'y':
        return y;
      case 'z':
        return z;
      default:
        throw std::runtime_error("Have no element " + std::string(1, c));
    }
  }
  T ***operator[](AxialDirection d) const {
    switch (d) {
      case AxialDirection::X:
        return x;
      case AxialDirection::Y:
        return y;
      case AxialDirection::Z:
        return z;
      default:
        throw std::runtime_error("Have no element " + std::to_string(d));
    }
  }

  /**
   * @brief Allocates x, y, and z as (K_total+1) * (J_total+1) * (I_total+1)
   * arrays
   *
   * @param I_total,J_total,K_total Dimensions of the tensor size to set
   */
  void allocate(int I_total, int J_total, int K_total) {
    x = (T ***) malloc(K_total * sizeof(T **));
    y = (T ***) malloc(K_total * sizeof(T **));
    z = (T ***) malloc(K_total * sizeof(T **));
    for (int k = 0; k < K_total; k++) {
      x[k] = (T **) malloc(J_total * sizeof(T *));
      y[k] = (T **) malloc(J_total * sizeof(T *));
      z[k] = (T **) malloc(J_total * sizeof(T *));
      for (int j = 0; j < J_total; j++) {
        x[k][j] = (T *) malloc(I_total * sizeof(T));
        y[k][j] = (T *) malloc(I_total * sizeof(T));
        z[k][j] = (T *) malloc(I_total * sizeof(T));
      }
    }
  }
};

class XYZVectors {
public:
  double *x = nullptr;
  double *y = nullptr;
  double *z = nullptr;

  /**
   * Default constructor
   */
  XYZVectors() = default;

  /**
   * Set the pointer for one of the vectors in this collection with a name of c
   * @param c Character labeling the vector
   * @param ptr Pointer to assign
   */
  void set_ptr(char c, double *ptr);
  /**
   * Set the pointer for one of the vectors in this collection with a name of c
   * @param d AxialDirection labeling the vector
   * @param ptr Pointer to assign
   */
  void set_ptr(AxialDirection d, double *ptr);

  /**
   * @brief Determines whether all elements in the x, y, or z vector are less
   * than a given value.
   *
   * @param comparison_value Value to compare elements to
   * @param vector_length Number of elements to compare against
   * @param component Vector to compare elements against; x, y, or z
   * @param buffer_start Only compare elements between buffer_start (inclusive)
   * and buffer_start+vector_length-1 (inclusive)
   * @return true All elements are less than the comparison_value
   * @return false At least one element is not less than the comparison_value
   */
  bool all_elements_less_than(double comparison_value, int vector_length,
                              AxialDirection component,
                              int buffer_start = 0) const;
  /**
   * @brief Determines whether all elements in the x, y, AND z vectors are less
   * than a given value.
   *
   * @param comparison_value Value to compare elements to
   * @param nx,ny,nz Number of elements in the nx, ny, and nz vectors
   * respectively
   * @return true All elements are less than the comparison_value
   * @return false At least one element is not less than the comparison_value
   */
  bool all_elements_less_than(double comparison_value, int nx, int ny,
                              int nz) const;
};

// TODO: docstring
class MaterialCollection {
protected:
  static void init_xyz_vectors(const mxArray *ptr, XYZVectors &arrays,
                               const std::string &prefix);
};

/**
 * @brief A class to encapsulate collection of algebraic terms in the
 * discretized forms of Maxwells equations for E fields. The symbol chosen in
 * the original reference is \f$C\f$.
 *
 * @details Algebraic terms \f$C_{a,b,c}\f$ defined in Section 4.2 of Munro, P,.
 * "Application of numerical methods to high numerical aperture imaging", 2006,
 * PhD thesis, Imperial College London.
 *
 * The definitions are equations 4.13, 4.14 (pp 82-3). Part of Maxwell's E-field
 * equations in equations 4.7-9.
 */
class CCollectionBase {
public:
  XYZVectors a;
  XYZVectors b;
  XYZVectors c;
};

/*! @copydoc CCollectionBase */
class CCollection : public CCollectionBase {
private:
  void init_xyz_vectors(const mxArray *ptr, XYZVectors &arrays,
                        const std::string &prefix);

public:
  bool is_multilayer = false;
  bool is_disp_ml = false;

  explicit CCollection(const mxArray *ptr);
};

/*! @copydoc CCollectionBase */
class CMaterial : public CCollectionBase, MaterialCollection {
public:
  explicit CMaterial(const mxArray *ptr);
};

/**
 * @brief A class to encapsulate collection of algebraic terms in the
 * discretized forms of Maxwells equations for H fields. The symbol chosen in
 * the original reference is \f$D\f$.
 *
 * @details Algebraic terms \f$D_{a,b}\f$ defined in Section 4.2 of Munro, P,.
 * "Application of numerical methods to high numerical aperture imaging", 2006,
 * PhD thesis, Imperial College London.
 *
 * The definitions are equations 4.15, 4.16 (pp 82-3). Part of Maxwell's H-field
 * equations in equations 4.10-12.
 */
class DCollectionBase {
public:
  XYZVectors a;
  XYZVectors b;
};

/*! @copydoc DCollectionBase */
class DCollection : public DCollectionBase {
private:
  static void init_xyz_vectors(const mxArray *ptr, XYZVectors &arrays,
                               const std::string &prefix);

public:
  explicit DCollection(const mxArray *ptr);
};

/*! @copydoc DCollectionBase */
class DMaterial : public DCollectionBase, MaterialCollection {
public:
  explicit DMaterial(const mxArray *ptr);
};

struct DispersiveMultiLayer {
public:
  std::vector<double> alpha;
  std::vector<double> beta;
  std::vector<double> gamma;
  xyz_vector<double> kappa;
  xyz_vector<double> sigma;

  /**
   * @brief Determines whether the (background) medium is dispersive
   *
   * @param near_zero_tolerance Tolerance for non-zero gamma (attenuation)
   * values
   * @return true Background is dispersive
   * @return false Background is not dispersive
   */
  bool is_dispersive(double near_zero_tolerance = 1e-15) const;
};

template<typename T>
class Matrix {
protected:
  int n_rows = 0;
  int n_cols = 0;
  T **matrix = nullptr;

public:
  /**
   * @brief Construct a new Matrix object, without assigned elements
   *
   */
  Matrix() = default;
  /**
   * @brief Construct a new Matrix object, providing the dimensions
   *
   * @param n_rows,n_cols Number of rows and columns in the matrix
   * @param initial_value The initial value of the elements, defaults to 0 to
   * avoid initalised but unassigned values
   */
  Matrix(int n_rows, int n_cols) { allocate(n_rows, n_cols); }

  inline T *operator[](int value) const { return matrix[value]; }
  /**
   * @brief Check whether this matrix has elements assigned
   *
   * @return true If this matrix has assigned elements
   * @return false This matrix is currently unassigned
   */
  bool has_elements() { return matrix != nullptr; };

  /**
   * Allocate the memory for this matrix
   *
   * @param n_rows Number of rows
   * @param n_cols Number of columns
   */
  void allocate(int n_rows, int n_cols, T initial_value = 0) {
    this->n_rows = n_rows;
    this->n_cols = n_cols;

    matrix = (T **) malloc(sizeof(T *) * n_rows);
    for (int i = 0; i < n_rows; i++) {
      matrix[i] = (T *) malloc(sizeof(T) * n_cols);
    }
  };

  int get_n_cols() const { return n_cols; }
  int get_n_rows() const { return n_rows; }

  /**
   * Destructor. Must be defined in the header
   */
  ~Matrix() {
    if (has_elements()) {
      for (int i = 0; i < n_rows; i++) { free(matrix[i]); }
      free(matrix);
    }
  };
};

template<typename T>
class Vector {
protected:
  int n = 0;          // Number of elements
  T *vector = nullptr;// Internal array

public:
  Vector() = default;

  explicit Vector(const mxArray *ptr) {
    n = (int) mxGetNumberOfElements(ptr);
    vector = (T *) malloc((unsigned) (n * sizeof(T)));

    auto matlab_ptr = mxGetPr(ptr);
    for (int i = 0; i < n; i++) { vector[i] = (T) matlab_ptr[i]; }
  }

  bool has_elements() { return vector != nullptr; }

  inline T operator[](int value) const { return vector[value]; };

  inline int size() const { return n; };
};

class FrequencyExtractVector : public Vector<double> {
public:
  FrequencyExtractVector(const mxArray *ptr, double omega_an);

  double max();
};

struct FrequencyVectors {
  std::vector<double> x;
  std::vector<double> y;
};

/**
 * List of field components as integers
 */
class FieldComponentsVector : public Vector<int> {
public:
  FieldComponentsVector() = default;

  void initialise(const mxArray *ptr);

  /**
   * Get the index of a particular integer in this vector. If it does not exist
   * then return -1. Returns the first occurrence.
   * @param value value to find in the vector
   * @return index or -1
   */
  int index(int value);
};

class Vertices : public Matrix<int> {
public:
  Vertices() = default;

  void initialise(const mxArray *ptr);

  int n_vertices() { return n_rows; }

  ~Vertices() {
    if (has_elements()) { free_cast_matlab_2D_array(matrix); }
    matrix = nullptr;
  };
};

class DetectorSensitivityArrays {
public:
  fftw_complex *v = nullptr;          // Flat fftw vector
  fftw_plan plan = nullptr;           // fftw plan for the setup
  std::complex<double> **cm = nullptr;// Column major matrix

  void initialise(int n_rows, int n_cols);

  ~DetectorSensitivityArrays();
};

/**
 * Container for storing snapshots of the full-field
 */
class FullFieldSnapshot {
public:
  std::complex<double> Ex = 0.;//< x-component of the electric field
  std::complex<double> Ey = 0.;//< y-component of the electric field
  std::complex<double> Ez = 0.;//< z-component of the electric field
  std::complex<double> Hx = 0.;//< x-component of the magnetic field
  std::complex<double> Hy = 0.;//< y-component of the magnetic field
  std::complex<double> Hz = 0.;//< z-component of the magnetic field

  FullFieldSnapshot() = default;

  /**
   * @brief Return the component of the field corresponding to the index
   * provided.
   *
   * 0 = Ex, 1 = Ey, 2 = Ez, 3 = Hx, 4 = Hy, 5 = Hz.
   * This is the indexing order that other storage containers use.
   *
   * Throws an error if provided an index <0 or >5.
   *
   * @param index Field component index to fetch
   * @return std::complex<double> The field component
   */
  std::complex<double> operator[](int index) {
    switch (index) {
      case 0:
        return Ex;
        break;
      case 1:
        return Ey;
        break;
      case 2:
        return Ez;
        break;
      case 3:
        return Hx;
        break;
      case 4:
        return Hy;
        break;
      case 5:
        return Hz;
        break;
      default:
        throw std::runtime_error("Index " + std::to_string(index) +
                                 " does not correspond to a field component.");
        break;
    }
  }

  /**
   * @brief Multiplies the electric field components by `factor`.
   * @param factor to scale the field by
   */
  void multiply_E_by(std::complex<double> factor) {
    Ex *= factor;
    Ey *= factor;
    Ez *= factor;
  }
  /**
   * @brief Multiplies the magnetic field components by `factor`.
   * @param factor to scale the field by
   */
  void multiply_H_by(std::complex<double> factor) {
    Hx *= factor;
    Hy *= factor;
    Hz *= factor;
  }
};
