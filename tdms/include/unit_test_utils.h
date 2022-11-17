/**
 * @file unit_test_utils.h
 * @brief Functions common to multiple unit tests
 */
#include <complex>

/**
 * @brief Determines if two numerical values are close by relative comparison.
 *
 * Checks the truth value of the condition
 * abs(a - b) / max(abs(a), abs(b)) < tol.
 * If a and b are both close to zero, returns true.
 *
 * @param a,b Numerical values
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
