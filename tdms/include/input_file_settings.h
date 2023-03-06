/**
 * @file input_file_settings.h
 * @author William Graham (ccaegra@ucl.ac.uk)
 * @brief Handles optional variables in the input file that are potentially not
 * present, setting them to their default values if they cannot be retrieved.
 * Also defines supporting enums for the settings these inputs control.
 */
#pragma once

#include "globals.h"
#include "mat_io.h"
#include "matlabio.h"

//! Solver methods implemented for time-propagation of the simulation
enum SolverMethod { PseudoSpectral = 0, FiniteDifference = 1 };

//! Interpolation methods that can be used to extract field values at Yee
//! cell centres
enum InterpolationMethod { Cubic = 1, BandLimited = 2 };

/**
 * @brief Class that handles retrieval of optional variables from the input
 * file. If a variable is not present, its default value (documented in the
 * corresponding member variable) is assigned instead.
 */
class InputFileSettings {
private:
  /*! Solver method for the time-propagation code. Reads from usecd in the input
   * file. If usecd is not present, or equal to 1, we use FDTD, otherwise PSTD
   */
  SolverMethod solver_method = SolverMethod::FiniteDifference;
  /*! The interpolation methods to use when extracting field values at Yee cell
   * centres. Reads from intmethod in the input file. If not present, or equal
   * to 1, we use Cubic interpolation, otherwise we prioritise band-limited
   * interpolation */
  InterpolationMethod pim = InterpolationMethod::Cubic;

public:
  InputFileSettings() = default;
  InputFileSettings(const char *input_filename) {
    set_from_input_file(input_filename);
  }

  /**
   * @brief Reads options (that may or may not be present) from the input file
   and sets the corresponding flags and member variables.
   *
   * @param input_filename Name of the input file to read options from
   */
  void set_from_input_file(const char *input_filename) {
    // Open input file for extraction of optional settings
    MATFile *mat_file = matOpen(input_filename, "r");

    // Set the preferred interpolation methods
    auto ptr_to_intmethod = matGetVariable(mat_file, "intmethod");
    if (ptr_to_intmethod != nullptr) {
      int intmethod = int_cast_from_double_in(ptr_to_intmethod, "intmethod");
      if (intmethod != 1) { pim = InterpolationMethod::BandLimited; }
    }

    // Set the solver method
    auto ptr_to_usecd = matGetVariable(mat_file, "usecd");
    if (ptr_to_usecd != nullptr) {
      int usecd = int_cast_from_double_in(ptr_to_usecd, "usecd");
      if (usecd == 0) { solver_method = SolverMethod::PseudoSpectral; }
    }
  }

  /*! @copydoc solver_method */
  SolverMethod solver() const { return solver_method; }
  /*! @copydoc pim */
  InterpolationMethod interpolation() const { return pim; }
};
