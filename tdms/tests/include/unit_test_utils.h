/**
 * @file unit_test_utils.h
 * @brief Functions common to multiple unit tests
 */
#pragma once

#include <complex>
#include <filesystem>
#include <random>
#include <string>

namespace tdms_tests {

inline double TOLERANCE = 1e-16;//< Floating-point comparison tolerance

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
inline bool is_close(T a, T b, double tol = 1E-10, double close_to_zero_tol = 1E-30) {

  auto max_norm = std::max(std::abs(a), std::abs(b));

  if (max_norm < close_to_zero_tol) {// Prevent dividing by zero
    return true;
  }

  return std::abs(a - b) / std::max(std::abs(a), std::abs(b)) < tol;
}

/**
 * @brief Determines whether an error value is better than a benchmark, or sufficiently close to be insignificant.
 *
 * If the error to check is superior (IE, closer to zero in absolute value) than the benchmark, return true.
 * Otherwise, use relative comparison to determine if the errors are sufficiently similar.
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
  std::string subdir = "tdms_unit_tests_" + std::to_string(random_number_generator());
  auto path = tmp / subdir;

  // mkdir and return the path to the directory we've just created
  std::filesystem::create_directory(path);
  return path;
}

}// namespace tdms_tests
