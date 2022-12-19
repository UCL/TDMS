#include "surface_phasors.h"

#include <complex>

#include "globals.h"

using namespace std;
using tdms_math_constants::DCPI;
using tdms_math_constants::IMAGINARY_UNIT;

SurfacePhasors::SurfacePhasors(mxArray *mx_surface_vertices) {
  set_from_matlab_array(mx_surface_vertices);
}

void SurfacePhasors::set_from_matlab_array(mxArray *mx_surface_vertices) {
  const mwSize *dimensions_pointer_out = mxGetDimensions(mx_surface_vertices);

  // set the number of vertices on the surface using the array dimensions
  n_surface_vertices = dimensions_pointer_out[0];

  //cast the vertex array as a 2-d integer array
  surface_vertices = cast_matlab_2D_array((int *) mxGetPr((mxArray *) mx_surface_vertices),
                                          dimensions_pointer_out[0], dimensions_pointer_out[1]);
}

void SurfacePhasors::extractPhasorsSurface(double **surface_EHr, double **surface_EHi,
                                           ElectricSplitField &E, MagneticSplitField &H, int n,
                                           double omega, int Nt, int J_tot,
                                           SimulationParameters &params, bool interpolate) {
  int vindex;
  FullFieldSnapshot F;
  complex<double> phaseTermE, phaseTermH, cphaseTermE, cphaseTermH;

  phaseTermE = fmod(omega * ((double) n) * params.dt, 2 * DCPI);
  phaseTermH = fmod(omega * ((double) n + 0.5) * params.dt, 2 * DCPI);

  cphaseTermH = exp(phaseTermH * IMAGINARY_UNIT) * 1. / ((double) Nt);
  cphaseTermE = exp(phaseTermE * IMAGINARY_UNIT) * 1. / ((double) Nt);

  //loop over every vertex in the list
#pragma omp parallel default(shared) private(F, phaseTermE, phaseTermH, vindex)
  {
    if (interpolate) {
#pragma omp for
      for (vindex = 0; vindex < n_surface_vertices; vindex++) {
        CellCoordinate current_cell(surface_vertices[0][vindex], surface_vertices[1][vindex],
                                    surface_vertices[2][vindex]);
        switch (params.dimension) {
          case Dimension::THREE:
            F.Ex = E.interpolate_to_centre_of(AxialDirection::X, current_cell);
            F.Ey = E.interpolate_to_centre_of(AxialDirection::Y, current_cell);
            F.Ez = E.interpolate_to_centre_of(AxialDirection::Z, current_cell);
            F.Hx = H.interpolate_to_centre_of(AxialDirection::X, current_cell);
            F.Hy = H.interpolate_to_centre_of(AxialDirection::Y, current_cell);
            F.Hz = H.interpolate_to_centre_of(AxialDirection::Z, current_cell);
            break;
          case Dimension::TRANSVERSE_ELECTRIC:
            F.Ex = E.interpolate_to_centre_of(AxialDirection::X, current_cell);
            F.Ey = E.interpolate_to_centre_of(AxialDirection::Y, current_cell);
            F.Hz = H.interpolate_to_centre_of(AxialDirection::Z, current_cell);
            break;
          case Dimension::TRANSVERSE_MAGNETIC:
            F.Ez = E.interpolate_to_centre_of(AxialDirection::Z, current_cell);
            F.Hx = H.interpolate_to_centre_of(AxialDirection::X, current_cell);
            F.Hy = H.interpolate_to_centre_of(AxialDirection::Y, current_cell);
            break;
          default:
            throw runtime_error("Dimension was not recognised!");
            break;
        }
        // multiply by phase factors
        F.multiply_E_by(cphaseTermE);
        F.multiply_H_by(cphaseTermH);

        // update the master arrays
        update_surface_EH(surface_EHr, surface_EHi, vindex, F);
      }
    } else {
#pragma omp for
      for (vindex = 0; vindex < n_surface_vertices; vindex++) {
        CellCoordinate current_cell(surface_vertices[0][vindex], surface_vertices[1][vindex],
                                    surface_vertices[2][vindex]);

        F.Ex = E.xy[current_cell] + E.xz[current_cell];
        F.Ey = E.yx[current_cell] + E.yz[current_cell];
        F.Ez = E.zx[current_cell] + E.zy[current_cell];
        F.Hx = H.xy[current_cell] + H.xz[current_cell];
        F.Hy = H.yx[current_cell] + H.yz[current_cell];
        F.Hz = H.zx[current_cell] + H.zy[current_cell];
        // multiply by phase factors
        F.multiply_E_by(cphaseTermE);
        F.multiply_H_by(cphaseTermH);

        // update the master arrays
        update_surface_EH(surface_EHr, surface_EHi, vindex, F);
      }
    }//end parallel region
  }
}

void SurfacePhasors::create_vertex_list(GridLabels input_grid_labels) {
  // setup vertex list dimensions
  mwSize dims[2];
  dims[0] = n_surface_vertices;
  dims[1] = 3;

  // create vertex list and infer MATLAB data storage location
  vertex_list = mxCreateNumericArray(2, (const mwSize *) dims, mxDOUBLE_CLASS, mxREAL);
  vertex_list_data_ptr = cast_matlab_2D_array(mxGetPr((mxArray *) vertex_list), dims[0], dims[1]);

  //now populate the vertex list
  for (int i = 0; i < n_surface_vertices; i++) {
    vertex_list_data_ptr[0][i] = input_grid_labels.x[surface_vertices[0][i]];
    vertex_list_data_ptr[1][i] = input_grid_labels.y[surface_vertices[1][i]];
    vertex_list_data_ptr[2][i] = input_grid_labels.z[surface_vertices[2][i]];
  }
}

SurfacePhasors::~SurfacePhasors() {
  // free the surface_vertices array if it has been cast
  if (surface_vertices) { free_cast_matlab_2D_array(surface_vertices); }

  // free the vertices_list if it was created and allocated
  if (vertex_list_data_ptr) { free_cast_matlab_2D_array(vertex_list_data_ptr); }
}

void update_surface_EH(double **surface_EHr, double **surface_EHi, int vindex, FullFieldSnapshot F) {
  // Ex
  surface_EHr[0][vindex] += real(F.Ex);
  surface_EHi[0][vindex] += imag(F.Ex);
  // Hx
  surface_EHr[3][vindex] += real(F.Hx);
  surface_EHi[3][vindex] += imag(F.Hx);
  // Ey
  surface_EHr[1][vindex] += real(F.Ey);
  surface_EHi[1][vindex] += imag(F.Ey);
  //Hy
  surface_EHr[4][vindex] += real(F.Hy);
  surface_EHi[4][vindex] += imag(F.Hy);
  //Ez
  surface_EHr[2][vindex] += real(F.Ez);
  surface_EHi[2][vindex] += imag(F.Ez);
  // Hz
  surface_EHr[5][vindex] += real(F.Hz);
  surface_EHi[5][vindex] += imag(F.Hz);
}
