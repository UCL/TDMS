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
#include "arrays/xyz_vector.h"
#include "globals.h"
#include "matlabio.h"
#include "utils.h"

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

// TODO: docstring
class MaterialCollection {
protected:
  static void init_xyz_vectors(const mxArray *ptr, XYZVector &arrays,
                               const std::string &prefix);
};

/**
 * @brief A class to encapsulate collection of algebraic terms in the
 * discretized forms of Maxwell's equations for E fields. The symbol chosen in
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
  XYZVector a;
  XYZVector b;
  XYZVector c;
};

/*! @copydoc CCollectionBase */
class CCollection : public CCollectionBase {
private:
  void init_xyz_vectors(const mxArray *ptr, XYZVector &arrays,
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
 * discretized forms of Maxwell's equations for H fields. The symbol chosen in
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
  XYZVector a;
  XYZVector b;
};

/*! @copydoc DCollectionBase */
class DCollection : public DCollectionBase {
private:
  static void init_xyz_vectors(const mxArray *ptr, XYZVector &arrays,
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
  XYZVector kappa;
  XYZVector sigma;

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
