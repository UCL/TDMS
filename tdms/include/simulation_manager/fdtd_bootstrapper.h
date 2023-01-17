#pragma once

#include <complex>

#include "arrays.h"
#include "cell_coordinate.h"
#include "field.h"
#include "simulation_parameters.h"

/**
 * @brief Class that handles the FDTD bootstrapping variables.
 *
 * NOTE: Can be deleted from codebase without altering functionality. See https://github.com/UCL/TDMS/issues/214#issue-1532363539 - holding off doing so until meeting clarification.
 */
class FDTDBootstrapper {
    private:
      // Electric and magnetic source phasors, used in bootstrapping
      Matrix<std::complex<double>> Ex, Ey, Hx, Hy;
    public:
      FDTDBootstrapper() = default;
      FDTDBootstrapper(IJKDims IJK_tot) { allocate_memory(IJK_tot); }

      /**
       * @brief Allocate memory for the bootstrapping variables, given the number of Yee cells
       *
       * @param IJK_tot Number of Yee cells in each axial direction
       */
      void allocate_memory(IJKDims IJK_tot);

      /**
       * @brief Extract phasors on the user-defined plane
       *
       * @param E_s,H_s Split-field values
       * @param IJK_tot Number of Yee cells in each axial direction
       * @param K1
       * @param tind Current iteration number
       * @param params Simulation parameters for this run of TDMS
       */
      void extract_phasors_in_plane(ElectricSplitField &E_s, MagneticSplitField &H_s, IJKDims IJK_tot, int K1, int tind, SimulationParameters params);
};
