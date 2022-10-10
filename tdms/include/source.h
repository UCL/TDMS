#include <string>
#include <stdexcept>
#include "matlabio.h"

class Source{
public:
  double*** real = nullptr;
  double*** imag = nullptr;

  Source(const mxArray *ptr, int dim1, int dim2, const std::string &name);
};
