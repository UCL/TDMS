#pragma once

#include <spdlog/spdlog.h>

#include "matrix.h"

#include "arrays.h"
#include "field.h"
#include "fieldsample.h"
#include "globals.h"
#include "grid_labels.h"
#include "input_matrices.h"
#include "interface.h"
#include "shapes.h"
#include "simulation_parameters.h"
#include "source.h"
#include "vertex_phasors.h"

class Iterator_IndependentObjectsFromInfile {
private:
  //Set the I_tot, J_tot, K_tot from the corresponding members of E_s
  void set_IJK_tot_from_E_s() {
    I_tot = E_s.I_tot;
    J_tot = E_s.J_tot;
    K_tot = E_s.K_tot;
  }
public:
    /* VARIABLES NEEDED IN THE MAIN LOOP,
        DERIVED FROM THE INPUT-FILE OR COMMAND-LINE */

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
  XYZVectors rho_cond;                      // < conductive aux
  DispersiveMultiLayer matched_layer;       //< dispersive aux
  IncidentField Ei;                         //< time-domain field
  FrequencyVectors f_vec;                   //< frequency vector
  Pupil pupil;                              //< TODO
  DTilde D_tilde;                           //< TODO
  FieldSample fieldsample;               //< handles field-value extraction at user-defined vertices
  VertexPhasors vertex_phasors;          //< handles phasor extraction at user-defined vertices
  TDFieldExporter2D ex_td_field_exporter;//< two-dimensional field exporter

  /* VARIABLES NOT NEEDED IN THE MAIN LOOP,
        DERIVED FROM THE INPUT-FILE / COMMAND-LINE */

  GridLabels input_grid_labels;//< cartesian labels of the Yee cells

  /* DERIVED VARIABLES FROM INDEPENDENT INPUTS */

  int I_tot, J_tot, K_tot;//< total number of Yee cells in the x,y,z directions respectively

  Iterator_IndependentObjectsFromInfile(InputMatrices matrices_from_input_file,
                                        SolverMethod solver_method);

  ~Iterator_IndependentObjectsFromInfile() {
    spdlog::info("Destroying Iterator_IndependentObjectsFromInfile");
  }
};

class Iterator_ObjectsFromInfile : public Iterator_IndependentObjectsFromInfile {
public:
    /* DECLARE VARIABLES NEEDED IN THE MAIN LOOP, DERIVED FROM INPUT FILE / COMMAND-LINE */

    Source Isource, Jsource, Ksource;//< TODO
    GratingStructure structure;//< TODO
    FrequencyExtractVector f_ex_vec;//< Vector of frequencies to extract field & phasors at

    /* DECLARE VARIABLES NOT NEEDED IN THE MAIN LOOP, DERIVED FROM INPUT FILE */

    /* DECLARE VARIABLES _DERIVED_ FROM VARIABLES OBTAINED FROM THE INPUT FILE */

    Iterator_ObjectsFromInfile(InputMatrices matrices_from_input_file, SolverMethod solver_method);

    ~Iterator_ObjectsFromInfile() { spdlog::info("Destroying Iterator_ObjectsFromInfile"); }
};
