/**
 * @file surface_phasors.h
 * @brief Contains a class that handles the complex amplitude extraction.
 */
#pragma once

#include "arrays.h"
#include "field.h"
#include "grid_labels.h"

class SurfacePhasors {
private:
  int **surface_vertices = nullptr;//< Pointer to the vertices on the surface
  int n_surface_vertices = 0;      //< Number of vertices on the surface

  mxArray *vertex_list = nullptr;         //< List of vertices
  double **vertex_list_data_ptr = nullptr;//< Buffer to place vertex data into MATLAB array

  mxArray *mx_surface_amplitudes = nullptr;//< Complex surface amplitudes
  int f_ex_vector_size = 0; //< Number of elements in the frequency extraction vector

  /* Storage for real and imag parts of mx_surface_amplitudes (these can be f_ex_vector_size * n_surface_vertices arrays of FullFieldSnapshots when MATLAB is removed!)

  Arrays are index by [frequency_index][field component][vertex_id/number]. Frequency index corresponds to the frequencies at which the user has requested we extract the amplitudes.
  */
  double ***surface_EHr = nullptr, ***surface_EHi = nullptr;

public:
  SurfacePhasors() = default;
  /**
   * @brief Construct a new Surface Phasors object, using the set_from_matlab_array() method
   */
  SurfacePhasors(mxArray *mx_surface_vertices, int _f_ex_vector_size);
  /**
   * @brief Sets the surface_vertices pointer and number of surface vertices tracker from the MATLAB array passed in
   *
   * @param mx_surface_vertices MATLAB array containing the vertex information
   * @param _f_ex_vector_size The FrequencyExtractionVector size, which we need to reserve sufficient memory
   */
  void set_from_matlab_array(mxArray *mx_surface_vertices, int f_ex_vector_size);

  /**
   * @brief Zeros the surface_EH{r,i} arrays
   */
  void zero_surface_EH() {
    for (int k = 0; k < f_ex_vector_size; k++) {
      for(int j = 0; j < 6; j++) { // 6 field components: Ex, y, z, and Hx, y, z
        for (int i = 0; i < n_surface_vertices; i++) {
          surface_EHr[k][j][i] = 0.;
          surface_EHi[k][j][i] = 0.;
        }
      }
    }
  }

  /**
   * @brief Get the number of surface vertices
   */
  int get_n_surface_vertices() { return n_surface_vertices; };

  /**
   * @brief Get the list of vertices
   */
  mxArray *get_vertex_list() { return vertex_list; };

  /**
   * @brief Get the array of complex surface amplitudes
   */
  mxArray *get_mx_surface_amplitudes() { return mx_surface_amplitudes; };

  /**
   * @brief Normalise the surface amplitudes at frequency_vector_index by the E- and H-norms provided.
   *
   * E-field components in surface_EH are divided by the (complex) Enorm.
   * H-field components in surface_EH are divided by the (complex) Hnorm.
   *
   * @param frequency_vector_index Frequency index, surface_EH{r,i}[frequency_vector_index] will be normalised
   * @param Enorm,Hnorm The {E,H}-norm to normalise the {E,H}-components by
   */
  void normalise_surface(int frequency_index, std::complex<double> Enorm, std::complex<double> Hnorm);

  /**
   * @brief Extract the phasor values at the vertices on the surface, for the given frequency index
   *
   * @param frequency_index The entries in surface_EH{r,i}[frequency_index] will be written to
   * @param E,H The electric,magnetic field
   * @param n Current timestep index
   * @param omega Angular frequency
   * @param Nt The number of timesteps in a sinusoidal period
   * @param J_tot Number of cells in the y-direction
   * @param params The parameters for this simulation
   * @param interpolate If true, perform interpolation on the fields when extracting phasors
   */
  void extractPhasorsSurface(int frequency_index, ElectricSplitField &E,
                             MagneticSplitField &H, int n, double omega, int Nt, int J_tot,
                             SimulationParameters &params, bool interpolate = true);

  /**
   * @brief Pulls the GridLabels information of vertices on the surface into vertex_list, a continuous block of memory.
   *
   * The vertex_list attribute will consist only of vertices that lie on the surface that the given class instance defines.
   *
   * @param input_grid_labels
   */
  void create_vertex_list(GridLabels input_grid_labels);

  /**
   * @brief Incriments surface_EH{r,i} at the given index by the field values provided.
   *
   * If we allow element-wise assignment, we are essentially performing the operations:
   * surface_EHr[frequency_vector_index][:][vertex_index] += F.real(),
   * surface_EHi[frequency_vector_index][:][vertex_index] += F.imag().
   *
   * @param frequency_index Frequency vector index (k) to assign to
   * @param vertex_index Vertex index (i) to assign to
   * @param F Field values to assign
   */
  void update_surface_EH(int frequency_index, int vertex_index, FullFieldSnapshot F);


  ~SurfacePhasors();
};
