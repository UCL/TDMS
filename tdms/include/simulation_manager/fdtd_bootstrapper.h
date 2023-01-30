/**
 * @file fdtd_bootstrapper.h
 * @brief Handles the boostrapping variables used in the FDTD solver.
 *
 * NOTE: Can be deleted from codebase without altering functionality. See https://github.com/UCL/TDMS/issues/214#issue-1532363539 - holding off doing so until meeting clarification.
 */
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
      FDTDBootstrapper(const IJKDimensions& IJK_tot) { allocate_memory(IJK_tot); }

      /**
       * @brief Allocate memory for the bootstrapping variables, given the number of Yee cells
       *
       * @param IJK_tot Number of Yee cells in each axial direction
       */
      void allocate_memory(const IJKDimensions& IJK_tot);

      /**
       * @brief Extract phasors on the user-defined plane
       *
       * @param E_s,H_s Split-field values
       * @param IJK_tot Number of Yee cells in each axial direction
       * @param K1
       * @param tind Current iteration number
       * @param params Simulation parameters for this run of TDMS
       */
      void extract_phasors_in_plane(const ElectricSplitField &E_s, const MagneticSplitField &H_s,
                                    const IJKDimensions &IJK_tot, int K1, int tind,
                                    const SimulationParameters &params);
};
