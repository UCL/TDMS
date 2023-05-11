/**
 * @file incident_field.h
 * @brief Declaration for the IncidentField object.
 */
#pragma once

#include <string>

#include "arrays/tensor3d.h"
#include "matrix.h"

class IncidentField {
protected:
  void set_component(Tensor3D<double> &component, const mxArray *ptr,
                     const std::string &name);

public:
  Tensor3D<double> x;
  Tensor3D<double> y;

  explicit IncidentField(const mxArray *ptr);
};
