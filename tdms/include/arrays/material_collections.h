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
  std::string input_field = is_material ? "Dmaterial" : "D";
};

// Try this for the C-stuff too? This one is easy since the construction methods
// _are_ the same barring the input_field we read from
// Also can we just do funky ? : assignment??

// struct DBase<true> {
//   XYZVector a, b;
//   std::string input_field = "Dmaterial";
// };
// struct DBase<false> {
//   XYZVector a, b;
//   std::string input_field = "D";
// };

typedef DBase<true> DMaterial;
typedef DBase<false> DCollection;
