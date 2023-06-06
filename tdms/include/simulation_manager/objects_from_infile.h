/**
 * @file objects_from_infile.h
 * @brief Classes that unpack variables from the input files the TDMS executable
 * receives and initialise fields.h and array.h datatypes to store this data.
 */
#pragma once

#include <spdlog/spdlog.h>

#include "matrix.h"

#include "arrays.h"
#include "arrays/cuboid.h"
#include "arrays/dispersive_multilayer.h"
#include "arrays/dtilde.h"
#include "arrays/incident_field.h"
#include "arrays/material_collections.h"
#include "arrays/tdms_matrix.h"
#include "arrays/vector_typedefs.h"
#include "arrays/xyz_vector.h"
#include "cell_coordinate.h"
#include "field.h"
#include "fieldsample.h"
#include "globals.h"
#include "grid_labels.h"
#include "input_flags.h"
#include "input_matrices.h"
#include "interface.h"
#include "simulation_parameters.h"
#include "source.h"
#include "vertex_phasors.h"

/**
 * @brief Class that constructs C++ arrays from an InputMatrices object and the
 * SolverMethod enum.
 *
 * Arrays that are constructed by this class are independent of one another and
 * other setup parameters - they can simply be read from the InputMatrices
 * object.
 *
 * The ObjectsFromInfile class handles the additional matrices in InputMatrices
 * that are inter-dependent.
 *
 */
class IndependentObjectsFromInfile {
public:
  int skip_tdf;//< Either 1 if we are using PSTD, or 6 if using FDTD

  SimulationParameters params;//< The parameters for this simulation

  ElectricSplitField E_s;//< The split electric-field values
  MagneticSplitField H_s;//< The split magnetic-field values

  uint8_t ***materials;//< TODO
  CMaterial Cmaterial; //< TODO
  DMaterial Dmaterial; //< TODO
  CCollection C;       //< TODO
  DCollection D;       //< TODO

  double *freespace_Cbx;       //< freespace constants
  double *alpha, *beta, *gamma;//< dispersion parameters

  InterfaceComponent I0, I1, J0, J1, K0,
          K1;   //< user-defined interface components
  Cuboid cuboid;//< user-defined surface to extract phasors over

  XYZVector rho_cond;                //< conductive aux
  DispersiveMultiLayer matched_layer;//< dispersive aux

  IncidentField Ei;//< time-domain field

  FrequencyVectors f_vec;//!< Frequency vector
  Pupil pupil;           //!< Numerical aperture of the objective lens
  DTilde D_tilde;        //!< TODO
  TDFieldExporter2D ex_td_field_exporter;//!< two-dimensional field exporter

  GridLabels input_grid_labels;//!< cartesian labels of the Yee cells

  FrequencyExtractVector
          f_ex_vec;//< Vector of frequencies to extract field & phasors at

  /* DERIVED VARIABLES FROM INDEPENDENT INPUTS */

  IJKDimensions IJK_tot;//!< total number of Yee cells in the x,y,z directions
                        //!< respectively
  int Nsteps;//!< Number of dfts to perform before checking for phasor
             //!< convergence

  IndependentObjectsFromInfile(InputMatrices matrices_from_input_file,
                               const InputFlags &in_flags);

  ~IndependentObjectsFromInfile();
};

/**
 * @brief Class that handles the creation of C++ arrays from the matrices passed
 * in an InputMatrices object, but also handling interdependencies between
 * inputs from this file.
 *
 * The Sources and GratingStructure can only be initialised after the other
 * arrays in the input file have been parsed.
 *
 * As such, this object inherits the setup of IndependentObjectsFromInfile, and
 * then constructs the aforementioned arrays using a combination of information
 * from the input file and previously mentioned setup.
 */
class ObjectsFromInfile : public IndependentObjectsFromInfile {
public:
  Source Isource, Jsource, Ksource;//< TODO
  GratingStructure structure;      //< TODO

  ObjectsFromInfile(InputMatrices matrices_from_input_file,
                    const InputFlags &in_flags);

  /** @brief Determine whether the {IJK}source terms are empty (true) or not
   * (false) */
  bool all_sources_are_empty() {
    return Isource.is_empty() && Jsource.is_empty() && Ksource.is_empty();
  }

  ~ObjectsFromInfile();
};
