#pragma once

#include <complex>

#include "arrays.h"
#include "cell_coordinate.h"
#include "field.h"
#include "simulation_parameters.h"

/**
 * @brief Class that handles the FDTD bootstrapping variables.
 */
class FDTDBootstrapper {
    private:
      // Electric and magnetic source phasors, used in bootstrapping
      Matrix<std::complex<double>> Ex, Ey, Hx, Hy;
    public:
      FDTDBootstrapper() = default;
      FDTDBootstrapper(IJKTotal IJK_tot) { allocate_memory(IJK_tot); }

      void allocate_memory(IJKTotal IJK_tot);

      void extract_phasors_in_plane(ElectricSplitField &E_s, MagneticSplitField &H_s, IJKTotal IJK_tot, int K1, int tind, SimulationParameters params);
};
