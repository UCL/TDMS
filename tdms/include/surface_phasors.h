/**
 * @file surface_phasors.h
 * @brief
 */
#pragma once

#include "field.h"

class SurfacePhasors {
private:
  int **surface_vertices = nullptr; //< Pointer to the vertices on the surface
  int n_surface_vertices = 0;       //< Number of vertices on the surface

public:
  SurfacePhasors() = default;

  SurfacePhasors(mxArray *mx_surface_vertices);

  void set_from_matlab_array(mxArray *mx_surface_vertices);

  /**
   * @brief Get the n surface vertices object
   */
  int get_n_surface_vertices() { return n_surface_vertices; };

  void extractPhasorsSurface(double **surface_EHr, double **surface_EHi, ElectricSplitField &E,
                             MagneticSplitField &H, int n, double omega, int Nt, int J_tot,
                             SimulationParameters &params);

  void extractPhasorsSurfaceNoInterpolation(double **surface_EHr, double **surface_EHi,
                                            ElectricSplitField &E, MagneticSplitField &H, int n,
                                            double omega, int Nt, int J_tot,
                                            SimulationParameters &params);
};
