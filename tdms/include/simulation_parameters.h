#pragma once
#include <string>
#include "arrays.h"


struct PerfectlyMatchedLayer {
  int Dxl = 0;  // Thickness of lower pml in the x direction
  int Dxu = 0;  //              upper
  int Dyl = 0;
  int Dyu = 0;
  int Dzl = 0;
  int Dzu = 0;
};

enum SourceMode{
  steadystate,
  pulsed
};

enum RunMode{
  complete,
  analyse
};

enum Dimension{
  THREE,  // Full dimensionality - compute all H and E components
  TE,     // Transverse electric - only compute Ex, Ey, and Hz components
  TM      // Transverse magnetic - only compute Hx, Hy, and Ez components
};

enum InterpolationMethod{
  null,
  cubic,
  band_limited
};

/**
 * A three-tuple of integers that contain the stride in each direction
 * to extract the phasors on. the surface i.e. x = 2 means extract from
 * every 2nd Yee cell.
 */
struct PhasorInc{
  int x = 1;
  int y = 1;
  int z = 1;
};


/**
 * Enum defining a mapping to integers used in the MATLAB initialisation
 */
enum FieldComponents{
  Ex = 1,
  Ey,
  Ez,
  Hx,
  Hy,
  Hz
};

class SimulationParameters{

public:
    SimulationParameters();

    double       omega_an      = 0.0;       // Angular ω
    unsigned int Nt            = 0;         // Number of simulation steps
    unsigned int Np            = 0;         // Number of times to extract the phasors
    unsigned int Npe           = 0;         // Extract the phasors every this number of iterations
    unsigned int start_tind    = 0;         // Starting iteration number for the time steps
    double       dt            = 0.0;       // Time step
    bool         has_tdfdir    = false;     // Is the tdfdir (time domain field directory) defined?
    bool         is_multilayer = false;     // Is this simulation of a multilayer?
    bool         is_disp_ml    = false;     // Is the material dispersive?
    double       to_l          = 0.0;       // Time delay of pulse
    double       hwhm          = 0.0;       // Half width at half max of pulse
    PerfectlyMatchedLayer pml;              // Perfectly matched layer struct with size attributes
    bool         exphasorsvolume = false;   // Should phasors be extracted in the whole volume?
    bool         exphasorssurface = false;  // Should phasors be extracted on a surface?
    bool         intphasorssurface = false; // Should phasors be extracted/interpolated?
    RunMode      run_mode     = complete;   // Run mode
    SourceMode   source_mode  = pulsed;     // Source mode
    Dimension    dimension    = THREE;      // Dimensions to calculate in
    bool         is_structure = false;      // Has a grating structure been defined?
    bool         exdetintegral = false;     // TODO: detector sensitivity evaluation ?
    int          k_det_obs    = 0;          // TODO: no idea what this is?!
    double       z_obs        = 0.0;        // Value of the input grid_labels_z at k_det_obs
    bool         air_interface_present = false;
    double       air_interface = 0.0;       // TODO: what is this?!
    bool         interp_mat_props = false;  // Should the material properties be interpolated?
    InterpolationMethod interp_method = cubic; // Type of surface field interpolation to do
    bool         exi_present = false;       // Is the time dependent x incident field present?
    bool         eyi_present = false;       // Is the time dependent y incident field present?
    PhasorInc    phasorinc;                 // Surface stride for extracting phasors

    void set_run_mode(std::string mode_string);

    void set_source_mode(std::string mode_string);

    void set_dimension(std::string mode_string);

    void set_phasorinc(const double* vector);

    void set_Np(FrequencyExtractVector &f_ex_vec);
};
