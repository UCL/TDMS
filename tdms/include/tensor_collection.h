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
