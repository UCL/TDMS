#pragma once

#include <complex.h>
#include <string>
#include <vector>

#include "arrays/tensor3d.h"
#include "matrix.h"

class DTilde {
protected:
  int n_det_modes = 0;
  static void set_component(Tensor3D<std::complex<double>> &tensor,
                            const mxArray *ptr, const std::string &name,
                            int n_rows, int n_cols);

public:
  DTilde() = default;

  inline int num_det_modes() const { return n_det_modes; };

  Tensor3D<std::complex<double>> x;
  Tensor3D<std::complex<double>> y;

  void initialise(const mxArray *ptr, int n_rows, int n_cols);
};
