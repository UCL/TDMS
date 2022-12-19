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
                                           SimulationParameters &params) {
  int vindex;
  double Ex, Ey, Ez, Hx, Hy, Hz;
  complex<double> phaseTermE, phaseTermH, subResultE, subResultH, cphaseTermE, cphaseTermH;

  phaseTermE = fmod(omega * ((double) n) * params.dt, 2 * DCPI);
  phaseTermH = fmod(omega * ((double) n + 0.5) * params.dt, 2 * DCPI);

  cphaseTermH = exp(phaseTermH * IMAGINARY_UNIT) * 1. / ((double) Nt);
  cphaseTermE = exp(phaseTermE * IMAGINARY_UNIT) * 1. / ((double) Nt);

  //loop over every vertex in the list
#pragma omp parallel default(shared) private(Ex, Ey, Ez, Hx, Hy, Hz, phaseTermE, phaseTermH,       \
                                             subResultE, subResultH, vindex)
  {
#pragma omp for
    for (vindex = 0; vindex < n_surface_vertices; vindex++) {
      CellCoordinate current_cell(surface_vertices[0][vindex], surface_vertices[1][vindex],
                                  surface_vertices[2][vindex]);
      if (params.dimension == Dimension::THREE) {
        // these should adapt to use 2/1D interpolation depending on whether the y-direction is available (hence no J_tot check)
        Hx = H.interpolate_to_centre_of(AxialDirection::X, current_cell);
        Hy = H.interpolate_to_centre_of(AxialDirection::Y, current_cell);
        Hz = H.interpolate_to_centre_of(AxialDirection::Z, current_cell);
        if (J_tot != 0) {
          Ex = E.interpolate_to_centre_of(AxialDirection::X, current_cell);
          Ey = E.interpolate_to_centre_of(AxialDirection::Y, current_cell);
          Ez = E.interpolate_to_centre_of(AxialDirection::Z, current_cell);
        } else {
          Ex = E.interpolate_to_centre_of(AxialDirection::X, current_cell);
          Ey = E.yx[current_cell] + E.yz[current_cell];
          Ez = E.interpolate_to_centre_of(AxialDirection::Z, current_cell);
        }
      } else if (params.dimension == Dimension::TRANSVERSE_ELECTRIC) {
        Ex = E.interpolate_to_centre_of(AxialDirection::X, current_cell);
        Ey = E.interpolate_to_centre_of(AxialDirection::Y, current_cell);
        Ez = 0.;
        Hx = 0.;
        Hy = 0.;
        Hz = H.interpolate_to_centre_of(AxialDirection::Z, current_cell);
      } else {
        Ex = 0.;
        Ey = 0.;
        Ez = E.interpolate_to_centre_of(AxialDirection::Z, current_cell);
        Hx = H.interpolate_to_centre_of(AxialDirection::X, current_cell);
        Hy = H.interpolate_to_centre_of(AxialDirection::Y, current_cell);
        Hz = 0.;
      }
      //    fprintf(stderr,"2nd interp donezn");

      /*Ex and Hx*/
      subResultH = Hx * cphaseTermH;//exp(phaseTermH * IMAGINARY_UNIT) * 1./((double) Nt);
      subResultE = Ex * cphaseTermE;//exp(phaseTermE * IMAGINARY_UNIT) * 1./((double) Nt);

      //now update the master array
      surface_EHr[0][vindex] = surface_EHr[0][vindex] + real(subResultE);
      surface_EHi[0][vindex] = surface_EHi[0][vindex] + imag(subResultE);

      surface_EHr[3][vindex] = surface_EHr[3][vindex] + real(subResultH);
      surface_EHi[3][vindex] = surface_EHi[3][vindex] + imag(subResultH);

      /*Ey and Hy*/
      subResultH = Hy * cphaseTermH;//exp(phaseTermH * IMAGINARY_UNIT) * 1./((double) Nt);
      subResultE = Ey * cphaseTermE;//exp(phaseTermE * IMAGINARY_UNIT) * 1./((double) Nt);

      //now update the master array
      surface_EHr[1][vindex] = surface_EHr[1][vindex] + real(subResultE);
      surface_EHi[1][vindex] = surface_EHi[1][vindex] + imag(subResultE);

      surface_EHr[4][vindex] = surface_EHr[4][vindex] + real(subResultH);
      surface_EHi[4][vindex] = surface_EHi[4][vindex] + imag(subResultH);


      /*Ez and Hz*/
      subResultH = Hz * cphaseTermH;//exp(phaseTermH * IMAGINARY_UNIT) * 1./((double) Nt);
      subResultE = Ez * cphaseTermE;//exp(phaseTermE * IMAGINARY_UNIT) * 1./((double) Nt);

      //now update the master array
      surface_EHr[2][vindex] = surface_EHr[2][vindex] + real(subResultE);
      surface_EHi[2][vindex] = surface_EHi[2][vindex] + imag(subResultE);

      surface_EHr[5][vindex] = surface_EHr[5][vindex] + real(subResultH);
      surface_EHi[5][vindex] = surface_EHi[5][vindex] + imag(subResultH);
    }
  }//end parallel region
}

void SurfacePhasors::extractPhasorsSurfaceNoInterpolation(double **surface_EHr, double **surface_EHi,
                                          ElectricSplitField &E, MagneticSplitField &H, int n,
                                          double omega, int Nt, int J_tot,
                                          SimulationParameters &params) {
  int vindex;
  double Ex, Ey, Ez, Hx, Hy, Hz;
  complex<double> phaseTermE, phaseTermH, subResultE, subResultH, cphaseTermE, cphaseTermH;

  phaseTermE = fmod(omega * ((double) n) * params.dt, 2 * DCPI);
  phaseTermH = fmod(omega * ((double) n + 0.5) * params.dt, 2 * DCPI);

  cphaseTermH = exp(phaseTermH * IMAGINARY_UNIT) * 1. / ((double) Nt);
  cphaseTermE = exp(phaseTermE * IMAGINARY_UNIT) * 1. / ((double) Nt);

  //loop over every vertex in the list
#pragma omp parallel default(shared) private(Ex, Ey, Ez, Hx, Hy, Hz, phaseTermE, phaseTermH,       \
                                             subResultE, subResultH, vindex)
  {
#pragma omp for
    for (vindex = 0; vindex < n_surface_vertices; vindex++) {
      //    fprintf(stderr,"vindex: %d: (%d %d %d)\n",vindex,surface_vertices[0][vindex],surface_vertices[1][vindex],surface_vertices[2][vindex]);
      Ex = E.xy[surface_vertices[2][vindex]][surface_vertices[1][vindex]]
               [surface_vertices[0][vindex]] +
           E.xz[surface_vertices[2][vindex]][surface_vertices[1][vindex]]
               [surface_vertices[0][vindex]];
      Ey = E.yx[surface_vertices[2][vindex]][surface_vertices[1][vindex]]
               [surface_vertices[0][vindex]] +
           E.yz[surface_vertices[2][vindex]][surface_vertices[1][vindex]]
               [surface_vertices[0][vindex]];
      Ez = E.zx[surface_vertices[2][vindex]][surface_vertices[1][vindex]]
               [surface_vertices[0][vindex]] +
           E.zy[surface_vertices[2][vindex]][surface_vertices[1][vindex]]
               [surface_vertices[0][vindex]];

      Hx = H.xy[surface_vertices[2][vindex]][surface_vertices[1][vindex]]
               [surface_vertices[0][vindex]] +
           H.xz[surface_vertices[2][vindex]][surface_vertices[1][vindex]]
               [surface_vertices[0][vindex]];
      Hy = H.yx[surface_vertices[2][vindex]][surface_vertices[1][vindex]]
               [surface_vertices[0][vindex]] +
           H.yz[surface_vertices[2][vindex]][surface_vertices[1][vindex]]
               [surface_vertices[0][vindex]];
      Hz = H.zx[surface_vertices[2][vindex]][surface_vertices[1][vindex]]
               [surface_vertices[0][vindex]] +
           H.zy[surface_vertices[2][vindex]][surface_vertices[1][vindex]]
               [surface_vertices[0][vindex]];

      /*Ex and Hx*/
      subResultH = Hx * cphaseTermH;//exp(phaseTermH * IMAGINARY_UNIT) * 1./((double) Nt);
      subResultE = Ex * cphaseTermE;//exp(phaseTermE * IMAGINARY_UNIT) * 1./((double) Nt);

      //now update the master array
      surface_EHr[0][vindex] = surface_EHr[0][vindex] + real(subResultE);
      surface_EHi[0][vindex] = surface_EHi[0][vindex] + imag(subResultE);

      surface_EHr[3][vindex] = surface_EHr[3][vindex] + real(subResultH);
      surface_EHi[3][vindex] = surface_EHi[3][vindex] + imag(subResultH);

      /*Ey and Hy*/
      subResultH = Hy * cphaseTermH;//exp(phaseTermH * IMAGINARY_UNIT) * 1./((double) Nt);
      subResultE = Ey * cphaseTermE;//exp(phaseTermE * IMAGINARY_UNIT) * 1./((double) Nt);

      //now update the master array
      surface_EHr[1][vindex] = surface_EHr[1][vindex] + real(subResultE);
      surface_EHi[1][vindex] = surface_EHi[1][vindex] + imag(subResultE);

      surface_EHr[4][vindex] = surface_EHr[4][vindex] + real(subResultH);
      surface_EHi[4][vindex] = surface_EHi[4][vindex] + imag(subResultH);


      /*Ez and Hz*/
      subResultH = Hz * cphaseTermH;//exp(phaseTermH * IMAGINARY_UNIT) * 1./((double) Nt);
      subResultE = Ez * cphaseTermE;//exp(phaseTermE * IMAGINARY_UNIT) * 1./((double) Nt);

      //now update the master array
      surface_EHr[2][vindex] = surface_EHr[2][vindex] + real(subResultE);
      surface_EHi[2][vindex] = surface_EHi[2][vindex] + imag(subResultE);

      surface_EHr[5][vindex] = surface_EHr[5][vindex] + real(subResultH);
      surface_EHi[5][vindex] = surface_EHi[5][vindex] + imag(subResultH);
    }
  }//end parallel region
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
