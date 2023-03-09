#include "simulation_manager/simulation_manager.h"

#include <omp.h>

using tdms_math_constants::DCPI, tdms_math_constants::IMAGINARY_UNIT;
using tdms_phys_constants::LIGHT_V;
using namespace std;

void SimulationManager::compute_detector_functions(unsigned int tind,
                                                   LoopVariables &lv) {
  /* This step is skipped if we are not running a pulsed simulation, in complete
  mode, with exdetintegral flagged, and are at an iteration number at which we
  want to extract phasors. Note the !negation at the start of this statement. */
  if (!(inputs.params.source_mode == SourceMode::pulsed &&// pulsed simulation
        inputs.params.run_mode == RunMode::complete &&    // complete run mode
        inputs.params.exdetintegral &&                    // exdetintegral flag
        (tind - inputs.params.start_tind) % inputs.params.Np ==
                0))// at an iteration number which we want the phasors
  {
    return;
  }

  // get the number of Yee cells in each axial direction
  IJKDimensions IJK_tot = n_Yee_cells();
  int I_tot = IJK_tot.i, J_tot = IJK_tot.j;

  double phaseTermE;
  complex<double> cphaseTermE;
  double lambda_an_t;//< Wavelength of light in free space at the current
                     // frequency
  complex<double> Idxt, Idyt, kprop;

  spdlog::debug("Setting Ex_t, Ey_t");

  /* First need to sum up the Ex and Ey values on a plane ready for FFT.
     Ex_t and Ey_t are in row-major format whilst split-field components are in
     column major format */
  for (int j = inputs.params.pml.Dyl; j < (J_tot - inputs.params.pml.Dyu); j++)
    for (int i = inputs.params.pml.Dxl; i < (I_tot - inputs.params.pml.Dxu);
         i++) {
      int m = j - inputs.params.pml.Dyl +
              (i - inputs.params.pml.Dxl) *
                      (J_tot - inputs.params.pml.Dyu - inputs.params.pml.Dyl);
      lv.Ex_t.v[m][0] = inputs.E_s.xy[inputs.params.k_det_obs][j][i] +
                        inputs.E_s.xz[inputs.params.k_det_obs][j][i];
      lv.Ex_t.v[m][1] = 0.;
      lv.Ey_t.v[m][0] = inputs.E_s.yx[inputs.params.k_det_obs][j][i] +
                        inputs.E_s.yz[inputs.params.k_det_obs][j][i];
      lv.Ey_t.v[m][1] = 0.;
    }

  // Fourier transform
  fftw_execute(lv.Ex_t.plan);
  fftw_execute(lv.Ey_t.plan);

  // Iterate over each mode
  for (int im = 0; im < inputs.D_tilde.num_det_modes(); im++) {
    // Convert back to column-major format
    for (int j = 0; j < (J_tot - inputs.params.pml.Dyu - inputs.params.pml.Dyl);
         j++)
      for (int i = 0;
           i < (I_tot - inputs.params.pml.Dxu - inputs.params.pml.Dxl); i++) {
        int m = j + i * (J_tot - inputs.params.pml.Dyu - inputs.params.pml.Dyl);
        lv.Ex_t.cm[j][i] = lv.Ex_t.v[m][0] + IMAGINARY_UNIT * lv.Ex_t.v[m][1];
        lv.Ey_t.cm[j][i] = lv.Ey_t.v[m][0] + IMAGINARY_UNIT * lv.Ey_t.v[m][1];
      }

    // Now multiply the pupil (NB: This is very sparse, typically non-zero at
    // only a few elements, speedup possible)
    for (int j = 0; j < (J_tot - inputs.params.pml.Dyu - inputs.params.pml.Dyl);
         j++)
      for (int i = 0;
           i < (I_tot - inputs.params.pml.Dxu - inputs.params.pml.Dxl); i++) {
        lv.Ex_t.cm[j][i] *= inputs.pupil[j][i] * inputs.D_tilde.x[j][i][im];
        lv.Ey_t.cm[j][i] *= inputs.pupil[j][i] * inputs.D_tilde.y[j][i][im];
      }

      /* Now iterate over each frequency we are extracting phasors at.
        We do this in parallel since, for each frequency, we only need to read
        values from the detector arrays. Ergo, we can handle computing phasors
        at different frequencies simultaneously. */
#pragma omp parallel default(shared) private(lambda_an_t, Idxt, Idyt, kprop,   \
                                             phaseTermE, cphaseTermE)
    {
#pragma omp for
      // For each frequency
      for (int ifx = 0; ifx < inputs.f_ex_vec.size(); ifx++) {
        // determine wavelength at this frequency
        lambda_an_t = LIGHT_V / inputs.f_ex_vec[ifx];
        Idxt = 0.;
        Idyt = 0.;

        // Loop over all angular frequencies
        for (int j = 0;
             j < (J_tot - inputs.params.pml.Dyu - inputs.params.pml.Dyl); j++)
          for (int i = 0;
               i < (I_tot - inputs.params.pml.Dxu - inputs.params.pml.Dxl);
               i++) {
            /* If the speed of light in the medium is less than that in free
              space:
              || lambda_an_t * f_vec ||^2 < 1  */
            if ((lambda_an_t * inputs.f_vec.x[i] * lambda_an_t *
                         inputs.f_vec.x[i] +
                 lambda_an_t * inputs.f_vec.y[j] * lambda_an_t *
                         inputs.f_vec.y[j]) < 1) {
              if (!inputs.params.air_interface_present) {
                // This had to be fixed since we must take into account the
                // refractive index of the medium.
                kprop = exp(
                        IMAGINARY_UNIT * inputs.params.z_obs * 2. * DCPI /
                        lambda_an_t * lv.refind *
                        sqrt(1. -
                             pow(lambda_an_t * inputs.f_vec.x[i] / lv.refind,
                                 2.) -
                             pow(lambda_an_t * inputs.f_vec.y[j] / lv.refind,
                                 2.)));
              } else {
                kprop = exp(IMAGINARY_UNIT *
                            (-inputs.params.air_interface +
                             inputs.params.z_obs) *
                            2. * DCPI / lambda_an_t * lv.refind *
                            sqrt(1. -
                                 pow(lambda_an_t * inputs.f_vec.x[i] /
                                             lv.refind,
                                     2.) -
                                 pow(lambda_an_t * inputs.f_vec.y[j] /
                                             lv.refind,
                                     2.))) *
                        exp(IMAGINARY_UNIT * inputs.params.air_interface * 2. *
                            DCPI / lambda_an_t *
                            sqrt(1. - pow(lambda_an_t * inputs.f_vec.x[i], 2.) -
                                 pow(lambda_an_t * inputs.f_vec.y[j], 2.)));
              }
            } else {
              kprop = 0.;
            }

            Idxt += lv.Ex_t.cm[j][i] * kprop;
            Idyt += lv.Ey_t.cm[j][i] * kprop;
          }
        phaseTermE = fmod(inputs.f_ex_vec[ifx] * 2. * DCPI * ((double) tind) *
                                  inputs.params.dt,
                          2 * DCPI);
        cphaseTermE = exp(phaseTermE * IMAGINARY_UNIT) * 1. /
                      ((double) inputs.params.Npe);

        outputs.ID.x[ifx][im] += Idxt * cphaseTermE;
        outputs.ID.y[ifx][im] += Idyt * cphaseTermE;

      }// end of loop on frequencies
    }  // end of pragma omp parallel
  }    // end of loop over each mode
}