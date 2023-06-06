/**
 * @file simulation_parameters.h
 * @brief Classes collecting parameters for the simulation.
 */
#pragma once

#include <string>

#include "arrays.h"
#include "arrays/vector_typedefs.h"
#include "input_matrices.h"

// Stores the thickness in each axial direction of the Perfectly Matched Layer
// (pml)
struct PerfectlyMatchedLayer {
  int Dxl = 0;//!< Thickness of lower pml in the x direction
  int Dxu = 0;//!< Thickness of upper pml in the x direction
  int Dyl = 0;//!< Thickness of lower pml in the y direction
  int Dyu = 0;//!< Thickness of upper pml in the y direction
  int Dzl = 0;//!< Thickness of lower pml in the z direction
  int Dzu = 0;//!< Thickness of upper pml in the z direction
};

/** @brief The x,y,z lengths of the cuboidal Yee cells */
struct YeeCellDimensions {
  double dx = 0.;//!< Yee cell width in the x direction
  double dy = 0.;//!< Yee cell width in the y direction
  double dz = 0.;//!< Yee cell width in the z direction
};

// Tracks whether the light source is pulsed or the simulation is at steady
// state
enum SourceMode { steadystate, pulsed };

enum RunMode { complete, analyse };

/** @brief Determines whether the simulation will compute all field components,
 * or only the TE or TM modal components */
enum Dimension {
  THREE,              //< Full dimensionality - compute all H and E components
  TRANSVERSE_ELECTRIC,//< Transverse electric - only compute Ex, Ey, and Hz
                      // components
  TRANSVERSE_MAGNETIC //< Transverse magnetic - only compute Hx, Hy, and Ez
                      // components
};

/**
 * @brief  A three-tuple of integers that contain the stride in each direction
 * to extract the phasors on.
 *
 * On the surface, a stride of EG x = 2 means extract from every 2nd Yee cell.
 */
struct SurfaceSpacingStride {
  int x = 1;
  int y = 1;
  int z = 1;
};

/**
 * @brief Class storing the various constants and behaviour flags for one
 * executation of the tdms executable.
 *
 * Stores physical constants like the angular frequency, and Yee cell
 * dimensions. Stores flags like whether the material is dispersive, or whether
 * we want to extract phasors on a predefined surface. Stores numerical method
 * constants like number of timesteps, length of timesteps.
 *
 */
class SimulationParameters {

public:
  SimulationParameters();

  double omega_an = 0.0;//< Angular ω
  unsigned int Nt = 0;  //< Number of simulation steps
  unsigned int Np = 0;  //< Number of times to extract the phasors
  unsigned int Npe = 0; //< Extract the phasors every this number of iterations
  unsigned int start_tind = 0;//< Starting iteration number for the time steps
  double dt = 0.0;            //< Time step
  bool has_tdfdir =
          false;//< Is the tdfdir (time domain field directory) defined?
  bool is_multilayer = false;//< Is this simulation of a multilayer?
  bool is_disp_ml = false;   //< Is the material dispersive?
  double to_l = 0.0;         //< Time delay of pulse
  double hwhm = 0.0;         //< Half width at half max of pulse
  PerfectlyMatchedLayer
          pml;//< Perfectly matched layer struct with size attributes
  bool exphasorsvolume =
          false;//< Should phasors be extracted in the whole volume?
  bool exphasorssurface = false;  //< Should phasors be extracted on a surface?
  bool intphasorssurface = false; //< Should phasors be extracted/interpolated?
  RunMode run_mode = complete;    //< Run mode
  SourceMode source_mode = pulsed;//< Source mode
  Dimension dimension = THREE;    //< Dimensions to calculate in
  bool is_structure = false;      //< Has a grating structure been defined?
  bool exdetintegral = false;     //< TODO: detector sensitivity evaluation ?
  int k_det_obs = 0;              //< TODO: no idea what this is?!
  double z_obs = 0.0;//< Value of the input grid_labels_z at k_det_obs
  bool air_interface_present = false;
  double air_interface = 0.0;//< TODO: what is this?!
  bool interp_mat_props =
          false;           //< Should the material properties be interpolated?
  bool exi_present = false;//< Is the time dependent x incident field present?
  bool eyi_present = false;//< Is the time dependent y incident field present?
  SurfaceSpacingStride spacing_stride;//< Surface stride for extracting phasors
                                      //(in matlab this is called 'phasorinc')
  YeeCellDimensions delta;            //< Yee cell dimensions (dx, dy, dz)

  void set_run_mode(std::string mode_string);

  void set_source_mode(std::string mode_string);

  void set_dimension(std::string mode_string);

  /**
   * @brief Set the surface spacing stride.
   * The x, y, z step size for extracting phasors (in matlab this is called
   * 'phasorinc')
   * @param vector the x, y, z steps (i.e. x = 2 means extract from every 2nd
   * Yee cell.)
   */
  void set_spacing_stride(const double *vector);

  /** @brief Calculate the number of times to extract the phasors (Np), and the
   * number of iterations between each extraction (Npe).
   */
  void set_Np_and_Npe(const FrequencyExtractVector &f_ex_vec);

  /**
   * @brief Unpacks all simulation parameters and flags from the matrix inputs
   * the executable received.
   *
   * @param in_matrices The input arrays from the input file(s) to tdms
   */
  void unpack_from_input_matrices(InputMatrices in_matrices);
};
