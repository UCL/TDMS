/**
 * @file objects_from_infile.h
 * @brief Classes that unpack variables from the input files the TDMS executable recieves and initalise fields.h and array.h datatypes to store this data.
 */
#pragma once

#include <stdexcept>

#include <spdlog/spdlog.h>

#include "arrays.h"
#include "cell_coordinate.h"
#include "field.h"
#include "fieldsample.h"
#include "globals.h"
#include "grid_labels.h"
#include "input_matrices.h"
#include "interface.h"
#include "matrix.h"
#include "shapes.h"
#include "simulation_parameters.h"
#include "source.h"
#include "vertex_phasors.h"

/**
 * @brief Class that constructs C++ arrays from an InputMatrices object and the SolverMethod enum.
 *
 * Arrays that are constructed by this class are independent of one another and other setup parameters - they can simply be read from the InputMatrices object.
 *
 * The ObjectsFromInfile class handles the additional matrices in InputMatrices that are inter-dependent.
 *
 */
class IndependentObjectsFromInfile {
public:
  /*! Either FDTD (default) or PSTD, the solver method */
  SolverMethod solver_method;
  int skip_tdf;//!< Either 1 if we are using PSTD, or 6 if using FDTD
  /*! Either Cubic (default) or BandLimited, the preferred interpolation methods */
  PreferredInterpolationMethods interpolation_methods;

  SimulationParameters params;//< The parameters for this simulation

  ElectricSplitField E_s;//!< The split electric-field values
  MagneticSplitField H_s;//!< The split magnetic-field values

  uint8_t ***materials;//!< TODO
  CMaterial Cmaterial; //!< TODO
  DMaterial Dmaterial; //!< TODO
  CCollection C;       //!< TODO
  DCollection D;       //!< TODO

  double *freespace_Cbx;       //!< freespace constants
  double *alpha, *beta, *gamma;//!< dispersion parameters

  /*! user-defined interface components */
  InterfaceComponent I0, I1, J0, J1, K0, K1;
  Cuboid cuboid;//!< user-defined surface to extract phasors over

  XYZVectors rho_cond;               //!< conductive aux
  DispersiveMultiLayer matched_layer;//!< dispersive aux

  IncidentField Ei;//< time-domain field

  FrequencyVectors f_vec;                //!< frequency vector
  Pupil pupil;                           //!< TODO
  DTilde D_tilde;                        //!< TODO
  TDFieldExporter2D ex_td_field_exporter;//!< two-dimensional field exporter

  GridLabels input_grid_labels;//!< cartesian labels of the Yee cells

  /* DERIVED VARIABLES FROM INDEPENDENT INPUTS */

  IJKDimensions IJK_tot;//!< total number of Yee cells in the x,y,z directions respectively
  int Nsteps;           //!< Number of dfts to perform before checking for phasor convergence

  IndependentObjectsFromInfile(InputMatrices matrices_from_input_file);

  /** @brief Set the solver method object using the flag from the input, and update dependent variables */
  void set_solver_method(int usecd_from_input_file) {
    switch (usecd_from_input_file) {
      case 0:
        solver_method = FiniteDifference;
        skip_tdf = 6;
        spdlog::info("Using finite-difference method (FDTD)");
        break;
      case 1:
        solver_method = PseudoSpectral;
        skip_tdf = 1;
        spdlog::info("Using pseudospectral method (PSTD)");
        break;
      default:
        throw std::runtime_error("Unknown solver method: " + to_string(usecd_from_input_file));
        break;
    }
  }

  /** @brief Set the preferred method of interpolation using the flag from the input, and update the fields about this change */
  void set_interpolation_method(const mxArray *interpolation_method_flag) {
    if (!mxIsEmpty(interpolation_method_flag)) {
      int flag_value = int_cast_from_double_in(interpolation_method_flag, "intmethod");
      switch (flag_value) {
        case 1:
          interpolation_methods = Cubic;
          spdlog::info("Restricting to cubic interpolation");
          break;
        case 2:
          interpolation_methods = BandLimited;
          spdlog::info("Using band-limited interpolation where possible");
          break;
        default:
          throw std::runtime_error("Interpolation methods not recognised:" + to_string(flag_value));
          break;
      }
    }
    // Enforce interpolation methods
    E_s.set_preferred_interpolation_methods(interpolation_methods);
    H_s.set_preferred_interpolation_methods(interpolation_methods);
  }

  ~IndependentObjectsFromInfile();
};

/**
 * @brief Class that handles the creation of C++ arrays from the matrices passed in an InputMatrices object, but also handling interdependencies between inputs from this file.
 *
 * The Sources, GratingStructure, and FrequencyExtractVector can only be initalised after the other arrays in the input file have been parsed.
 *
 * As such, this object inherits the setup of IndependentObjectsFromInfile, and then constructs the aforementioned arrays using a combination of information from the input file and previously mentioned setup.
 */
class ObjectsFromInfile : public IndependentObjectsFromInfile {
public:
  Source Isource, Jsource, Ksource;//< TODO
  GratingStructure structure;      //< TODO
  FrequencyExtractVector f_ex_vec; //< Vector of frequencies to extract field & phasors at

  ObjectsFromInfile(InputMatrices matrices_from_input_file);

  /** @brief Determine whether the {IJK}source terms are empty (true) or not (false) */
  bool all_sources_are_empty() {
    return Isource.is_empty() && Jsource.is_empty() && Ksource.is_empty();
  }

  ~ObjectsFromInfile();
};
