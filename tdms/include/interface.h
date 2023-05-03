/**
 * @file interface.h
 */
#pragma once

#include <string>

#include "mat_io.h"

/**
 * @brief Defines a plane over which a source/boundary condition is (or is not)
 * to be applied.
 *
 * There are 6 planes on which a source condition can be applied; I0, I1, J0,
 * J1, K0, and K1.
 * The {I,J,K} character indicates the axial direction to which the plane is
 * perpendicular, whilst the {0,1} character indicates whether this is the first
 * or second such plane perpendicular to that axial direction.
 *
 * The index member stores the value of the (constant) Yee cell index of all Yee
 * cells that lie in the plane defined. That is, index is the I-index of all Yee
 * cells in the I0 or I1 planes, the J-index for the J0 and J1 planes, and the
 * K-index of the K0 and K1 planes.
 *
 * The apply member flags whether an interface condition is to be applied across
 * that particular interface/plane.
 */
class InterfaceComponent {
public:
  /*! Whether or not a source or boundary condition is applied at this
   * interface */
  bool apply = false;
  /*! The value of the constant Yee-cell index for cells in this plane */
  int index = 0;

  InterfaceComponent() = default;
  InterfaceComponent(const mxArray *ptr, const std::string &name);
};
