#include "arrays/xyz_vector.h"

using namespace std;

// void XYZVectors::set_ptr(const char c, double *ptr) {
//   switch (c) {
//     case 'x': {
//       x = ptr;
//       break;
//     }
//     case 'y': {
//       y = ptr;
//       break;
//     }
//     case 'z': {
//       z = ptr;
//       break;
//     }
//     default:
//       throw std::runtime_error("Have no element " + std::string(1, c));
//   }
// }
// void XYZVectors::set_ptr(AxialDirection d, double *ptr) {
//   switch (d) {
//     case AxialDirection::X: {
//       x = ptr;
//       break;
//     }
//     case AxialDirection::Y: {
//       y = ptr;
//       break;
//     }
//     case AxialDirection::Z: {
//       z = ptr;
//       break;
//     }
//     default:
//       throw std::runtime_error("Have no element " + to_string(d));
//   }
// }

bool XYZVector::all_elements_less_than(double comparison_value,
                                       AxialDirection component) const {
  // To save copying into a new vector, declare a pointer that will point to the
  // component we wish to compare to
  const vector<double> *component_pointer;
  switch (component) {
    case AxialDirection::X:
      component_pointer = &x;
      break;
    case AxialDirection::Y:
      component_pointer = &y;
      break;
    case AxialDirection::Z:
      component_pointer = &z;
      break;
    default:
      throw runtime_error("Error - component not recognised");
      break;
  }
  // Search through all elements of the vector and return false if any pass the
  // threshold
  for (const double &element : *component_pointer) {
    if (element > comparison_value) { return false; }
  }
  return true;
}

bool XYZVector::all_elements_less_than(double comparison_value) const {
  if (!all_elements_less_than(comparison_value, AxialDirection::X)) {
    return false;
  } else if (!all_elements_less_than(comparison_value, AxialDirection::Y)) {
    return false;
  } else if (!all_elements_less_than(comparison_value, AxialDirection::Z)) {
    return false;
  }
  return true;
}
