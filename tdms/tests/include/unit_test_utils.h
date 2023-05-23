/**
 * @file unit_test_utils.h
 * @brief Functions common to multiple unit tests
 */
#pragma once

#include <complex>
#include <filesystem>
#include <random>
#include <string>

#include "globals.h"

using tdms_math_constants::DCPI;

namespace tdms_unit_test_data {

// Small test data files shipped with the test code.
#ifdef CMAKE_SOURCE_DIR
inline std::string
        tdms_object_data(std::string(CMAKE_SOURCE_DIR) +
                         "/tests/unit/hdf5_and_tdms_objects/class_data.mat");
#else
inline std::string
        tdms_object_data(std::filesystem::current_path() /
                         "../tests/unit/hdf5_and_tdms_objects/class_data.mat");
#endif

/* The struct_testdata is a .mat file containing a single structure array,
example_struct. This structure array has the following data.

MEMBER: value, MATLAB type
double_no_decimal: 1.0, double
double_half: 0.5, double
string: 'tdms', char array
boolean: 1, logical
uint_345: 3 * 4 * 5 array of 1s, uint8
double_22: [0.25, 0.5; 0.75, 1.], double
complex_22: [0., -i; i, 0], complex
*/
#ifdef CMAKE_SOURCE_DIR
inline std::string
        struct_testdata(std::string(CMAKE_SOURCE_DIR) +
                        "/tests/unit/hdf5_io_tests/structure_array.mat");
#else
inline std::string
        struct_testdata(std::filesystem::current_path() /
                        "../tests/unit/hdf5_io_tests/structure_array.mat");
#endif

}// namespace tdms_unit_test_data

namespace tdms_tests {

inline double TOLERANCE = 1e-16;//< Floating-point comparison tolerance

/**
 * @brief Determines whether a value is close to zero
 *
 * @param x Value to test
 * @param tol Tolerance for being "close" to zero
 */
inline bool near_zero(const double &x, double tol = TOLERANCE) {
  return std::abs(x) < tol;
}

/**
 * @brief Determines if two numerical values are close by relative comparison.
 *
 * Checks the truth value of the condition
 * abs(a - b) / max(abs(a), abs(b)) < tol.
 * If a and b are both close to zero, returns true.
 *
 * @param a,b Numerical values
 * @param tol Relative comparison tolerance
 * @param close_to_zero_tol Cutoff value for the "close to zero" criterion
 */
template<typename T>
inline bool is_close(T a, T b, double tol = 1E-10,
                     double close_to_zero_tol = 1E-30) {

  auto max_norm = std::max(std::abs(a), std::abs(b));

  if (near_zero(max_norm, close_to_zero_tol)) {// Prevent dividing by zero
    return true;
  }

  return std::abs(a - b) / max_norm < tol;
}

/**
 * @brief Determines whether an error value is better than a benchmark, or
 * sufficiently close to be insignificant.
 *
 * If the error to check is superior (IE, closer to zero in absolute value) than
 * the benchmark, return true. Otherwise, use relative comparison to determine
 * if the errors are sufficiently similar.
 *
 * @param to_check Numerical error to evaluate suitability of
 * @param to_beat Benchmark error
 * @param tol Relative comparison tolerance
 * @param close_to_zero_tol Cutoff value for the "close to zero" criterion
 * @return true to_check is a superior or equivalent error to to_beat
 * @return false to_check is an inferior error
 */
template<typename T>
inline bool is_close_or_better(T to_check, T to_beat, double tol = 1E-10,
                               double close_to_zero_tol = 1E-30) {
  if (std::abs(to_check) < std::abs(to_beat)) {
    // return true if the value to_check is better (closer to 0) than to_beat
    return true;
  } else {
    // determine if the numerical values are close by relative comparison
    return is_close(to_check, to_beat, tol, close_to_zero_tol);
  }
}

/**
 * @brief Compute the relative mean square difference of two arrays
 *
 * @param x,y Arrays to read from
 * @param n_elements Number of elements in the arrays
 * @param x_start,y_start Index to start reading the buffer from the x, y array
 * respectively (default 0)
 * @param close_to_zero_tol Tolerance for MSDs being zero (to avoid /0 errors)
 * @return double The relative mean square difference of x and y
 */
inline double
relative_mean_square_difference(double *x, double *y, int n_elements,
                                int x_start = 0, int y_start = 0,
                                double close_to_zero_tol = TOLERANCE) {
  double mean_sq_x = 0., mean_sq_y = 0., mean_sq_diff = 0.;
  for (int i = 0; i < n_elements; i++) {
    mean_sq_x += x[i + x_start] * x[i + x_start];
    mean_sq_y += y[i + y_start] * y[i + y_start];
    mean_sq_diff += (x[i + x_start] - y[i + y_start]) *
                    (x[i + x_start] - y[i + y_start]);
  }
  mean_sq_x = mean_sq_x / (double) n_elements;
  mean_sq_y = mean_sq_y / (double) n_elements;
  if (mean_sq_x < close_to_zero_tol && mean_sq_y < close_to_zero_tol) {
    return 0.;
  } else {
    mean_sq_diff = mean_sq_diff /
                   (((double) n_elements) * std::max(mean_sq_x, mean_sq_y));
    return mean_sq_diff;
  }
}

/**
 * @brief Computes the Euclidean norm of the vector provided
 *
 * @param v Vector or array
 * @param end (Inclusive) end of buffer to read vector from
 * @param start (Inclusive) start of buffer to read vector from
 * @return double Euclidean norm
 */
template<typename T>
inline double euclidean(T *v, int end, int start = 0) {
  double norm_val = 0.;
  for (int i = start; i < end; i++) { norm_val += std::norm(v[i]); }
  return std::sqrt(norm_val);
}

// returns the order of magnitude of the value x
inline int order_of_magnitude(double x) { return floor(log10(x)); }

/**
 * @brief Create a temporary directory for writing files.
 *
 * Creates a subdirectory (with randomised name) in the system tmp.  This can be
 * used for writing files when testing file I/O and field exporters and similar.
 *
 * @warning You should probably clean up after yourself: call
 *          `std::filesystem::remove_all` when you're done.
 *
 * @return std::filesystem::path Path to the temporary directory.
 */
inline std::filesystem::path create_tmp_dir() {
  // random number setup
  std::random_device device_seed;
  std::mt19937 random_number_generator(device_seed());

  // get system tmp directory (OS-dependent), add a uniquely named subdirectory
  auto tmp = std::filesystem::temp_directory_path();
  std::string subdir =
          "tdms_unit_tests_" + std::to_string(random_number_generator());
  auto path = tmp / subdir;

  // mkdir and return the path to the directory we've just created
  std::filesystem::create_directory(path);
  return path;
}

/** @brief Returns the char represented by a uint16 */
inline char uint16_to_char(const uint16_t &repr) { return *((char *) &repr); }

/**
 * @brief Returns the string composed of the characters represented by
 * subsequent uint16s.
 *
 * @param repr Buffer of uint16s that represent individual characters
 * @param buffer_length Length of the buffer
 * @return string The string composed of the converted characters
 */
inline std::string uint16s_to_string(uint16_t *repr, int buffer_length) {
  std::string output;
  for (int i = 0; i < buffer_length; i++) { output += uint16_to_char(repr[i]); }
  return output;
}

}// namespace tdms_tests
