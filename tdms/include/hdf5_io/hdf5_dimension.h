#pragma once

#include <algorithm>
#include <vector>

#include <H5Cpp.h>

class H5Dimension : public std::vector<hsize_t> {
public:
  H5Dimension() = default;
  H5Dimension(const H5::DataSpace &data_space);
  H5Dimension(const H5::DataSet &data_set) : H5Dimension(data_set.getSpace()){};

  /**
   * @brief Whether these dimensions describe an array that is castable to a 1D
   * array.
   * @details In the event that these dimensions only have one entry, or at
   * most one of the entries is greater than 1, the shape described can be cast
   * to a 1D-array of length max_dim().
   *
   * @return true These dimensions describe (an object castable to) a 1D-array
   * @return false Otherwise
   */
  bool is_1D() const;

  /**
   * @brief Returns the dimension of the greatest extent.
   * @details For instances where is_1D() returns true, this conincides with the
   * number of elements in the array, and the length of a 1D array necessary to
   * hold all the elements.
   *
   * @return hsize_t
   */
  hsize_t max_dim() const { return *std::max_element(begin(), end()); };

  /** @brief The total number of elements (product of the dimensions) */
  int number_of_elements() const {
    int product = 1;
    for (hsize_t axis_size : *this) { product *= axis_size; }
    return product;
  }
};
