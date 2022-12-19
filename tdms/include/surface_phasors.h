/**
 * @file surface_phasors.h
 * @brief
 */
#pragma once

#include "field.h"
#include "grid_labels.h"

class SurfacePhasors {
private:
  int **surface_vertices = nullptr;//< Pointer to the vertices on the surface
  int n_surface_vertices = 0;      //< Number of vertices on the surface

  mxArray *vertex_list = nullptr;         //< List of vertices
  double **vertex_list_data_ptr = nullptr;//< Buffer to place vertex data into MATLAB array

public:
  SurfacePhasors() = default;
  /**
   * @brief Construct a new Surface Phasors object, using the set_from_matlab_array() method
   */
  SurfacePhasors(mxArray *mx_surface_vertices);
  /**
   * @brief Sets the surface_vertices pointer and number of surface vertices tracker from the MATLAB array passed in
   *
   * @param mx_surface_vertices MATLAB array containing the vertex information
   */
  void set_from_matlab_array(mxArray *mx_surface_vertices);

  /**
   * @brief Get the n surface vertices object
   */
  int get_n_surface_vertices() { return n_surface_vertices; };

  /**
   * @brief Get the vertex list object
   */
  mxArray *get_vertex_list() { return vertex_list; };

  /**
   * @brief Extract the phasor values at the vertices on the surface
   *
   * @param surface_EHr,surface_EHi Location to write the real,imag parts of the surface amplitudes to
   * @param E,H The electric,magnetic field
   * @param n Current timestep index
   * @param omega Angular frequency
   * @param Nt The number of timesteps in a sinusoidal period
   * @param J_tot Number of cells in the y-direction
   * @param params The parameters for this simulation
   */
  void extractPhasorsSurface(double **surface_EHr, double **surface_EHi, ElectricSplitField &E,
                             MagneticSplitField &H, int n, double omega, int Nt, int J_tot,
                             SimulationParameters &params);

  void extractPhasorsSurfaceNoInterpolation(double **surface_EHr, double **surface_EHi,
                                            ElectricSplitField &E, MagneticSplitField &H, int n,
                                            double omega, int Nt, int J_tot,
                                            SimulationParameters &params);

  void create_vertex_list(GridLabels input_grid_labels);


  ~SurfacePhasors();
};
