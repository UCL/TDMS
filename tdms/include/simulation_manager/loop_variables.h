#pragma once

#include <complex>
#include <vector>

#include "arrays.h"
#include "field.h"
#include "matrix.h"
#include "objects_from_infile.h"

class LoopVariables {
private:
  // Pointers to the MATLAB structure that underlies E_copy, the convergence checker
  mxArray *E_copy_MATLAB_data[3] = {nullptr};

  /**
   * @brief Determine the dispersive properties and setup the corresponding arrays
   *
   * @param data The objects and data obtained from the input file
   */
  void setup_dispersive_properties(ObjectsFromInfile &data);
  /**
   * @brief Determine if we have a dispersive medium by searching for non-zero attenuation constants in gamma.
   *
   * @param materials
   * @param IJK_tot Number of Yee cells in each axial direction
   * @param attenuation_constants Material attenuation constants
   * @param dt Simulation timestep
   * @param non_zero_tol Tolerance for an attenuation constant being "non-zero"
   * @return true, we have a dispersive medium
   * @return false, we do not have a dispersive medium
   */
  bool is_dispersive(uint8_t ***materials, IJKDims IJK_tot, double *attenuation_constants,
                     double dt, double non_zero_tol = 1e-15);

  /**
   * @brief Performs an optimisation step in a 2D simulation (J_tot==0) when we have either TE or TM, but not both.
   *
   * In the J_tot==0 2D version; the 'TE' case involves components Ey, Hx and Hz, whilst the 'TM' case involves components Ex, Ez and Hy.
   *
   * The idea is to use an alternative upper limit to the loop over j when we have J_tot==0. As it stands, we are planning the loops listed below. We use the syntax (k k_min k_max : j j_min j_max : i i_min i_max) to represent the nested for loop:
   *    for( k = k_min; k < k_max; k++) {
   *        for( j = j_min; j < j_max; j++) {
   *            for( i = i_min; i < i_max; i++) {
   *            }
   *        }
   *    }
   *
   * Exy: Not involved in 2D, (k 0 K_tot+1 : j 1 J_tot : i 0 I_tot)
   * Exz: TM, (k 1 K_tot : j 0 J_tot+1 : i 0 I_tot)
   * Eyx: TE, (k 0 K_tot+1 : j 0 max(J_tot, 1) : i 1 I_tot)
   * Eyz: TE, (k 1 K_tot : j 0 max(J_tot, 1) : i 0 I_tot+1)
   * Ezx: TM, (k 0 K_tot : j 0 J_tot+1 : i 1 I_tot)
   * Ezy: Not involved in 2D, (k 0 K_tot : j 1 J_tot : i 0 I_tot+1)
   * Hxy: Not involved in 2D, (k 0 K_tot : j 0 J_tot : i 0 I_tot+1)
   * Hxz: TE, (k 0 K_tot : j 0 max(J_tot,1) : i 0 I_tot+1)
   * Hyx: TM, (k 0 K_tot : j 0 J_tot+1 : i 0 I_tot)
   * Hyz: TM, (k 0 K_tot : j 0 J_tot+1 : i 0 I_tot)
   * Hzx: TE, (k 0 K_tot+1 : j 0 max(J_tot,1) : i 0 I_tot)
   * Hzy: Not involved in 2D, (k 0 K_tot+1) : j 0 J_tot : i 0 I_tot)
   *
   * We see that in all cases, the TE update has the following loop on j:
   *    for(j=0;j<max(J_tot,1);j++)
   * whilst the TM case has:
   *    for(j=0;j<(J_tot+1);j++)
   *
   * So we can create variables
   * int J_tot_p1_bound, J_tot_bound
   * which would take the following values
   *    3D:
   *        J_tot_p1_bound = J_tot + 1;
   *        J_tot_bound = J_tot;
   *    2D:
   *        TE:
   *            J_tot_bound = 1;
   *        Not TE:
   *            J_tot_bound = 0;
   *        TM:
   *            J_tot_p1_bound = 1;
   *        Not TM:
   *            J_tot_p1_bound = 0;
   *
   * @param data Data and objects obtained from the input file
   * @param non_zero_tol Tolerance at which doubles are considered "non-zero"
   */
  void optimise_loop_J_range(ObjectsFromInfile &data, double non_zero_tol = 1e-15);

public:
  DetectorSensitivityArrays Ex_t, Ey_t;//< temporary storage for detector sensitivity evaluation

  ElectricField E_copy;//< Stores the phasors at the previous iteration, to check for convergence

  ElectricSplitField E_nm1;
  CurrentDensitySplitField
          J_c;//< The per-cell ( current density or conductivity ? ) of the material
  CurrentDensitySplitField J_s, J_nm1;

  EHVec eh_vec;

  std::vector<std::complex<double>>
          E_norm,//< Holds the E-field phasors norm at each extraction frequency
          H_norm;//< Holds the H-field phasors norm at each extraction frequency
  //std::complex<double> *E_norm = nullptr, *H_norm = nullptr;

  bool is_conductive, is_disp;//< Whether the materials are dispersive / conductive respectively

  double refind;//< refractive index of the first layer of the multilayer, or of the bulk of homogeneous

  int J_tot_p1_bound, J_tot_bound;
  int K;//< Number of non-pml cells in the K-direction (K_tot - Dxl - Dxu)

  LoopVariables(ObjectsFromInfile &data, IJKDims E_field_dims);

  ~LoopVariables();
};
