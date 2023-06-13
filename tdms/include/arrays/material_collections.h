/**
 * @file material_collections.h
 * @author William Graham
 * @brief Storage containers for Yee-cell constants that appear in the update
 * equations.
 */
#pragma once

#include <string>

#include "arrays/xyz_vector.h"

/**
 * @brief A class to encapsulate collection of algebraic terms in the
 * discretized forms of Maxwell's equations for E fields. The symbol chosen in
 * the original reference is \f$C\f$.
 *
 * @details Algebraic terms \f$C_{a,b,c}\f$ defined in Section 4.2 of Munro, P,.
 * "Application of numerical methods to high numerical aperture imaging", 2006,
 * PhD thesis, Imperial College London.
 *
 * The definitions are equations 4.13, 4.14 (pp 82-3). Part of Maxwell's E-field
 * equations in equations 4.7-9.
 */
struct CMaterial {
  XYZVector a, b, c;
};

/** @copydoc CMaterial */
struct CCollection : public CMaterial {
  bool is_multilayer = false;
  bool is_disp_ml = false;
};

/**
 * @brief A class to encapsulate collection of algebraic terms in the
 * discretized forms of Maxwell's equations for H fields. The symbol chosen in
 * the original reference is \f$D\f$.
 *
 * @details Algebraic terms \f$D_{a,b}\f$ defined in Section 4.2 of Munro, P,.
 * "Application of numerical methods to high numerical aperture imaging", 2006,
 * PhD thesis, Imperial College London.
 *
 * The definitions are equations 4.15, 4.16 (pp 82-3). Part of Maxwell's H-field
 * equations in equations 4.10-12.
 */
template<bool is_material>
struct DBase {
  XYZVector a, b;
  const std::string input_field = is_material ? "Dmaterial" : "D";
};

/** @copydoc DBase */
typedef DBase<true> DMaterial;
/** @copydoc DBase */
typedef DBase<false> DCollection;
