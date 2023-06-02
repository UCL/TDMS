/**
 * @file utils.h
 * @brief Useful miscellaneous utility functions
 */
#pragma once

#include <algorithm>
#include <iterator>
#include <string>
#include <vector>

/**
 * @brief Throws a runtime error if a file is not found.
 *
 * @param filename The name of the file to check.
 * @param mode The mode to try and open with.
 */
void assert_can_open_file(const char *filename, const char *mode);

/**
 * @brief Check two strings are equal
 *
 * @param a The first string
 * @param b The second string
 * @return true if the strings are the same
 * @return false otherwise
 */
bool are_equal(const char *a, const char *b);

/**
 * @brief A collection of utility functions that we will be applying to
 * std::vectors multiple times throughout the codebase.
 */
namespace tdms_vector_utils {
/** @brief Return the maximum value stored in the vector. */
template<typename T>
T max(const std::vector<T> &v) {
  return *std::max_element(v.begin(), v.end());
}

/** @brief Return true if the vector has non-zero size. */
template<typename T>
bool has_elements(const std::vector<T> &v) {
  return v.size() != 0;
}

/**
 * Get the index of a particular integer in this vector. If it does not exist
 * then return -1. Returns the first occurrence.
 * @param v Vector to search through for the value
 * @param value value to find
 * @return index of the value in the vector, or -1 (if not found)
 */
template<typename T>
int index(const std::vector<T> &v, const T &value) {
  auto found_index = std::find(v.begin(), v.end(), value);
  // Return -1 if we didn't find the value
  if (found_index == v.end()) {
    return -1;
  } else {
    // Otherwise return the index of the occurrence of the value
    return std::distance(v.begin(), found_index);
  }
}

/**
 * @brief Static_casts an array of doubles to ints. ONLY INTENDED FOR CASES
 * WHERE WE KNOW THAT THE VALUES STORED AS DOUBLES ARE INTS.
 * @details Usage cases are restricted to when data is read from MATLAB .mat
 * files, which by default save all values as floats, however we are expecting
 * to receive integer values (cell indices, sizes, etc).
 * @param ints_that_are_doubles Vector of doubles that actually contain integer
 * values
 * @return std::vector<int> Converted values
 */
std::vector<int>
to_vector_int(const std::vector<double> &ints_that_are_doubles);

}// namespace tdms_vector_utils
