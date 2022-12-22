#include "surface_phasors.h"

#include <complex>

#include "globals.h"

using namespace std;
using tdms_math_constants::DCPI;
using tdms_math_constants::IMAGINARY_UNIT;

SurfacePhasors::SurfacePhasors(mxArray *mx_surface_vertices, int _f_ex_vector_size) {
  set_from_matlab_array(mx_surface_vertices, _f_ex_vector_size);
}

void SurfacePhasors::set_from_matlab_array(mxArray *mx_surface_vertices, int _f_ex_vector_size) {
  f_ex_vector_size = _f_ex_vector_size;
  const mwSize *dimensions_pointer_out = mxGetDimensions(mx_surface_vertices);

  // set the number of vertices on the surface using the array dimensions
  n_surface_vertices = dimensions_pointer_out[0];

  //cast the vertex array as a 2-d integer array
  surface_vertices = cast_matlab_2D_array((int *) mxGetPr((mxArray *) mx_surface_vertices),
                                          dimensions_pointer_out[0], dimensions_pointer_out[1]);

  // Create space for the complex amplitudes E and H around the surface. These will be in a large complex
  // array with each line being of the form Re(Ex) Im(Ex) Re(Ey) ... Im(Hz). Each line corresponds to the
  // the vertex with the same line as in surface_phasors.surface_vertices.
  mwSize dims[3] = {n_surface_vertices, 6, f_ex_vector_size};
  mx_surface_amplitudes = mxCreateNumericArray(3, (const mwSize *) dims, mxDOUBLE_CLASS, mxCOMPLEX);
  surface_EHr = cast_matlab_3D_array(mxGetPr((mxArray *) mx_surface_amplitudes), dims[0], dims[1],
                                     dims[2]);
  surface_EHi = cast_matlab_3D_array(mxGetPi((mxArray *) mx_surface_amplitudes), dims[0], dims[1],
                                     dims[2]);
}

void SurfacePhasors::normalise_surface(int frequency_index, complex<double> Enorm,
                                       complex<double> Hnorm) {
  complex<double> normalised_amplitude;

  // E-field components
  for (int vindex = 0; vindex < n_surface_vertices; vindex++) {
    for (int i = 0; i < 3; i++) {
      // again, if surface_EH was just one array of FullFieldSnapshots, we could use the multiply_by methods of that class rather than doing this awkward cast-then-recast. Alas, MATLAB woes
      normalised_amplitude = complex<double>(surface_EHr[frequency_index][i][vindex],
                                             surface_EHi[frequency_index][i][vindex]) /
                             Enorm;
      surface_EHr[frequency_index][i][vindex] = normalised_amplitude.real();
      surface_EHi[frequency_index][i][vindex] = normalised_amplitude.imag();
    }
  }
  // H-field components
  for (int vindex = 0; vindex < n_surface_vertices; vindex++) {
    for (int i = 3; i < 6; i++) {
      normalised_amplitude = complex<double>(surface_EHr[frequency_index][i][vindex],
                                             surface_EHi[frequency_index][i][vindex]) /
                             Hnorm;
      surface_EHr[frequency_index][i][vindex] = normalised_amplitude.real();
      surface_EHi[frequency_index][i][vindex] = normalised_amplitude.imag();
    }
  }
}

void SurfacePhasors::extractPhasorsSurface(int frequency_index,
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
  /* Loop over every vertex in the surface.

  Since the value of the phasors at each vertex is entirely determined from the previously calculated fields,
  these computations can be done in parallel as the computation of the phasors is independent of one another.
  Ergo, we use a parallel loop.
  */
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
        update_surface_EH(frequency_index, vindex, F);
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
        update_surface_EH(frequency_index, vindex, F);
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

void SurfacePhasors::update_surface_EH(int frequency_index, int vertex_index,
                       FullFieldSnapshot F) {
  // FullFieldSnapshot uses the same indexing scheme as surface_EH{r,i} do, so we can just loop
  for(int component = 0; component < 6; component++) {
    surface_EHr[frequency_index][component][vertex_index] += real(F[component]);
    surface_EHi[frequency_index][component][vertex_index] += imag(F[component]);
  }
}

SurfacePhasors::~SurfacePhasors() {
  // free the surface_vertices array if it has been cast
  if (surface_vertices) { free_cast_matlab_2D_array(surface_vertices); }

  // free the vertices_list if it was created and allocated
  if (vertex_list_data_ptr) { free_cast_matlab_2D_array(vertex_list_data_ptr); }

  // free the surfaceEH arrays if they were allocated
  if (surface_EHr) { free_cast_matlab_3D_array(surface_EHr, f_ex_vector_size); }
  if (surface_EHi) { free_cast_matlab_3D_array(surface_EHi, f_ex_vector_size); }
}
