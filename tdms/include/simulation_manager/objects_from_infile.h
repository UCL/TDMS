/**
 * @file objects_from_infile.h
 * @brief Classes that unpack variables from the input files the TDMS executable recieves and initalise fields.h and array.h datatypes to store this data.
 */
#pragma once

#include "matrix.h"

#include "arrays.h"
#include "field.h"
#include "fieldsample.h"
#include "cell_coordinate.h"
#include "globals.h"
#include "grid_labels.h"
#include "input_matrices.h"
#include "interface.h"
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
  SolverMethod solver_method;               //< Either PSTD (default) or FDTD, the solver method
  int skip_tdf;                             //< Either 1 if we are using PSTD, or 6 if using FDTD
  SimulationParameters params;              //< The parameters for this simulation
  ElectricSplitField E_s;                   //< The split electric-field values
  MagneticSplitField H_s;                   //< The split magnetic-field values
  uint8_t ***materials;                     //< TODO
  CMaterial Cmaterial;                      //< TODO
  DMaterial Dmaterial;                      //< TODO
  CCollection C;                            //< TODO
  DCollection D;                            //< TODO
  double *freespace_Cbx;                    //< freespace constants
  double *alpha, *beta, *gamma;             //< dispersion parameters
  InterfaceComponent I0, I1, J0, J1, K0, K1;//< user-defined interface components
  Cuboid cuboid;                            //< user-defined surface to extract phasors over
  XYZVectors rho_cond;                      //< conductive aux
  DispersiveMultiLayer matched_layer;       //< dispersive aux
  IncidentField Ei;                         //< time-domain field
  FrequencyVectors f_vec;                   //< frequency vector
  Pupil pupil;                              //< TODO
  DTilde D_tilde;                           //< TODO
  TDFieldExporter2D ex_td_field_exporter;//< two-dimensional field exporter

  GridLabels input_grid_labels;//< cartesian labels of the Yee cells

  /* DERIVED VARIABLES FROM INDEPENDENT INPUTS */

  IJKDims IJK_tot;//< total number of Yee cells in the x,y,z directions respectively

  IndependentObjectsFromInfile(InputMatrices matrices_from_input_file,
                                        SolverMethod _solver_method);

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
    GratingStructure structure;//< TODO
    FrequencyExtractVector f_ex_vec;//< Vector of frequencies to extract field & phasors at

    ObjectsFromInfile(InputMatrices matrices_from_input_file, SolverMethod _solver_method);

    ~ObjectsFromInfile();
};
