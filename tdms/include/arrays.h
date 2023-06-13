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

struct FrequencyVectors {
  std::vector<double> x;
  std::vector<double> y;
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
