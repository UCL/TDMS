/**
 * @file field.h
 * @brief Classes for the electric and magnetic (split) fields on a grid.
 */
#pragma once

#include <complex>

#include "arrays.h"
#include "cell_coordinate.h"
#include "dimensions.h"
#include "globals.h"
#include "mat_io.h"
#include "simulation_parameters.h"

/**
 * A generic grid entity. For example:
 *
 *          ______________________ (2, 1, 0)
 *         /           /         /
 *        /           /         /
 *       /___________/_________/
 *  (0, 0, 0)   (1, 0, 0)     (2, 0, 0)
 *
 *  has I_tot = 2, J_tot = 1, K_tot = 0.
 *
 * NOTE: For storage purposes, this means that field values associated to cells
 * are stored _to the left_. That is, Grid(0,0,0) is associated to the cell
 * (-1,-1,-1). This is contrary to the way values are associated to cells, where
 * cell (0,0,0) is associated to the field values (0,0,0).
 */
class Grid {
protected:
  // the preferred interpolation methods (pim) for interpolating between the
  // grid values, default is BandLimited
  PreferredInterpolationMethods pim =
          PreferredInterpolationMethods::BandLimited;

public:
  // The {IJK}_tot values of this grid
  IJKDimensions tot = {0, 0, 0};

  /**
   * Maximum value out of I_tot, J_tot and K_tot
   * @return value
   */
  int max_IJK_tot() const { return tot.max(); };

  /**
   * @brief Set the preferred interpolation methods
   */
  void set_preferred_interpolation_methods(PreferredInterpolationMethods _pim) {
    pim = _pim;
  };
};

class SplitFieldComponent : public Tensor3D<double> {
public:
  int n_threads = 1;// Number of threads this component was chunked with
  fftw_plan *plan_f = nullptr;// Forward fftw plan
  fftw_plan *plan_b = nullptr;// Backward fftw plan

  double **operator[](int value) const { return tensor[value]; };
  double &operator[](ijk cell) const { return tensor[cell.k][cell.j][cell.i]; }

  void initialise_from_matlab(double ***tensor, Dimensions &dims);

  /**
   * Initialise a vector of 1d discrete Fourier transform plans
   * @param n_threads Number of threads that will be used
   * @param size Length of the vector
   * @param eh_vec TODO: what is this?
   */
  void initialise_fftw_plan(int n_threads, int size, EHVec &eh_vec);

  ~SplitFieldComponent();
};

/**
 * @brief A split field defined over a grid.
 *
 * To reconstruct the components we have e.g.: Ex = Exy + Exz multiplied by a
 * phase factor
 */
class SplitField : public Grid {
protected:
  virtual int delta_n() = 0;// TODO: no idea what this is or why it's needed

public:
  /*! Magnitude of the xy component at each grid point (i,j,k) */
  SplitFieldComponent xy;
  /*! Magnitude of the xz component at each grid point (i,j,k) */
  SplitFieldComponent xz;
  /*! Magnitude of the yx component at each grid point (i,j,k) */
  SplitFieldComponent yx;
  /*! Magnitude of the yz component at each grid point (i,j,k) */
  SplitFieldComponent yz;
  /*! Magnitude of the zx component at each grid point (i,j,k) */
  SplitFieldComponent zx;
  /*! Magnitude of the zy component at each grid point (i,j,k) */
  SplitFieldComponent zy;

  /**
   * Default no arguments constructor
   */
  SplitField() = default;

  /**
   * Constructor of the field with a defined size in the x, y, z Cartesian
   * dimensions
   */
  SplitField(int I_total, int J_total, int K_total);

  /**
   * Allocate the memory appropriate for all the 3D tensors associated with
   * this split field
   */
  void allocate();

  /**
   * Set all the values of all components of the field to zero
   */
  void zero();

  /**
   * Allocate and set to zero all components of the field
   */
  void allocate_and_zero() {
    allocate();
    zero();
  }

  /**
   * Initialise the fftw plans for all components
   * @param n_threads Number of threads to split over
   * @param eh_vec // TODO
   */
  void initialise_fftw_plan(int n_threads, EHVec &eh_vec);

  /**
   * @brief Fetches the largest absolute value of the field.
   *
   * Split field values are sums of the corresponding components, so this
   * function returns the largest absolute value of the entries in (xy + xz),
   * (yx + yz), (zx + zy)
   *
   * @return double Largest (by absolute value) field value
   */
  double largest_field_value();

  /**
   * @brief Interpolates a SplitField component to the centre of a Yee cell
   *
   * @param d SplitField component to interpolate
   * @param cell Index (i,j,k) of the Yee cell to interpolate to the centre of
   * @return double The interpolated field value
   */
  virtual double interpolate_to_centre_of(AxialDirection d,
                                          CellCoordinate cell) = 0;
};

class ElectricSplitField : public SplitField {
protected:
  int delta_n() override {
    return -1;
  };// TODO: no idea what this is or why it's needed

public:
  ElectricSplitField() = default;

  /**
   * Constructor of the field with a defined size in the x, y, z Cartesian
   * dimensions
   */
  ElectricSplitField(int I_total, int J_total, int K_total)
      : SplitField(I_total, J_total, K_total){};

  /**
   * @brief Interpolates a split E-field component to the centre of a Yee cell
   *
   * @param d Field component to interpolate
   * @param cell Index (i,j,k) of the Yee cell to interpolate to the centre of
   * @return double The interpolated component value
   */
  double interpolate_to_centre_of(AxialDirection d,
                                  CellCoordinate cell) override;
};

class MagneticSplitField : public SplitField {
protected:
  int delta_n() override {
    return 0;
  };// TODO: no idea what this is or why it's needed

public:
  MagneticSplitField() = default;

  /**
   * Constructor of the field with a defined size in the x, y, z Cartesian
   * dimensions
   */
  MagneticSplitField(int I_total, int J_total, int K_total)
      : SplitField(I_total, J_total, K_total){};

  /**
   * @brief Interpolates a split E-field component to the centre of a Yee cell
   *
   * @param d Field component to interpolate
   * @param cell Index (i,j,k) of the Yee cell to interpolate to the centre of
   * @return double The interpolated component value
   */
  double interpolate_to_centre_of(AxialDirection d,
                                  CellCoordinate cell) override;
};

class CurrentDensitySplitField : public SplitField {
protected:
  int delta_n() override { return 0; }

public:
  CurrentDensitySplitField() = default;

  /**
   * Constructor of the field with a defined size in the x, y, z Cartesian
   * dimensions
   */
  CurrentDensitySplitField(int I_total, int J_total, int K_total)
      : SplitField(I_total, J_total, K_total){};

  double interpolate_to_centre_of(AxialDirection d,
                                  CellCoordinate cell) override {
    return 0.;
  };
};

/**
 * A complex field defined over a grid. Has real and imaginary components
 * at each (x, y, z) grid point
 */
class Field : public Grid {
public:
  double ft = 0.;// TODO: an explanation of what this is

  std::complex<double> angular_norm = 0.;

  // TODO: this is likely better as a set of complex arrays - use
  // XYZTensor3D<std::complex<double>> This also makes implimenting
  // normalise_volume, and the interpolation schemes, much easier...
  XYZTensor3D<double> real;
  XYZTensor3D<double> imag;

  /**
   * Upper (u) and lower (l) indices in the x,y,z directions. e.g.
   * il is the first non-pml cell in the i direction and iu the last in the
   * corresponding split field grid
   */
  int il = 0, iu = 0, jl = 0, ju = 0, kl = 0, ku = 0;

  /**
   * Default no arguments constructor
   */
  Field() = default;

  /**
   * Constructor of the field with a defined size in the x, y, z Cartesian
   * dimensions
   */
  Field(int I_total, int J_total, int K_total);

  /**
   * Allocate the memory appropriate for all the 3D tensors associated with
   * this split field
   */
  void allocate();

  /**
   * Set all the values of all components of the field to zero
   */
  void zero();

  /**
   * Allocate and set to zero all components of the field
   */
  void allocate_and_zero() {
    allocate();
    zero();
  }

  /**
   * @brief Normalises the field entries by dividing by the angular norm.
   *
   * Specifically,
   * real[c][k][j][i] + i imag[c][k][j][i] = ( real[c][k][j][i] + i
   * imag[c][k][j][i] ) / angular_norm
   *
   */
  void normalise_volume();

  /**
   * Set the phasors for this field, given a split field. Result gives field
   * according to the exp(-iwt) convention
   * @param F The split-field to read values from
   * @param n The current timestep
   * @param omega Angular frequency
   * @param dt Timestep
   * @param Nt Number of timesteps in a sinusoidal period
   */
  void set_phasors(SplitField &F, int n, double omega, double dt, int Nt);

  // TODO: Docstring
  /**
   * @brief Compute the phasor_norm of the current field and add it to the
   * current norm-value.
   *
   * Note that the calls:
   * add_to_angular_norm(n, Nt, params), and
   * add_to_angular_norm(n % Nt, Nt, params),
   *
   * are equivalent up to numerical precision when using n = tind, and Nt =
   * Nsteps as set in the main simulation.
   *
   * @param n The number of timesteps INTO the current sinusoidal period
   * @param Nt The number of timesteps in a sinusoidal period
   * @param params The simulation parameters
   */
  void add_to_angular_norm(int n, int Nt, SimulationParameters &params);

  // TODO: Docstring
  std::complex<double> phasor_norm(double f, int n, double omega, double dt,
                                   int Nt);

  virtual double phase(int n, double omega, double dt) = 0;
  /**
   * @brief Interpolates a Field component to the centre of a Yee cell
   *
   * @param d Field component to interpolate
   * @param cell Index (i,j,k) of the Yee cell to interpolate to the centre of
   * @return std::complex<double> The interpolated field value
   */
  virtual std::complex<double>
  interpolate_to_centre_of(AxialDirection d, CellCoordinate cell) = 0;

  /**
   * @brief Interpolates the Field over the range provided.
   *
   * Default range is to interpolate to the midpoint of all consecutive points.
   *
   * @param[out] x_out,y_out,z_out Output arrays for interpolated values
   * @param i_lower,j_lower,k_lower Lower index for interpolation in the i,j,k
   * directions, respectively
   * @param i_upper,j_upper,k_upper Upper index for interpolation in the i,j,k
   * directions, respectively
   * @param mode Determines which field components to compute, based on the
   * simulation Dimension
   */
  void interpolate_over_range(mxArray *x_out, mxArray *y_out, mxArray *z_out,
                              int i_lower, int i_upper, int j_lower,
                              int j_upper, int k_lower, int k_upper,
                              Dimension mode = Dimension::THREE);
  /*! @copydoc interpolate_over_range */
  void interpolate_over_range(mxArray *x_out, mxArray *y_out, mxArray *z_out,
                              Dimension mode = Dimension::THREE);

  /**
   * @brief Interpolates the Field's transverse electric components to the
   * centre of Yee cell i,j,k
   *
   * @param[in] cell Yee cell index
   * @param[out] x_at_centre,y_at_centre,z_at_centre Addresses to write
   * interpolated values for the x,y,z components (respectively)
   */
  virtual void interpolate_transverse_electric_components(
          CellCoordinate cell, std::complex<double> *x_at_centre,
          std::complex<double> *y_at_centre,
          std::complex<double> *z_at_centre) = 0;
  /**
   * @brief Interpolates the Field's transverse magnetic components to the
   * centre of Yee cell i,j,k
   *
   * @param[in] cell Yee cell index
   * @param[out] x_at_centre,y_at_centre,z_at_centre Addresses to write
   * interpolated values for the x,y,z components (respectively)
   */
  virtual void interpolate_transverse_magnetic_components(
          CellCoordinate cell, std::complex<double> *x_at_centre,
          std::complex<double> *y_at_centre,
          std::complex<double> *z_at_centre) = 0;

  /**
   * Set the values of all components in this field from another, equally sized
   * field
   */
  void set_values_from(Field &other);

  /**
   * @brief Computes the maximum pointwise absolute difference of the other
   * field to this one, divided by the largest absolute value of this field's
   * components.
   *
   * Specifically, we compute and return the quantity
   * \f$ \frac{ \max_{i,j,k}(this[k][j][i] - other[k][j][i]) }{
   * \max_{i,j,k}this[k][j][i] } \f$ This quantity is used to determine
   * convergence of phasors.
   *
   * The other field must have the same dimensions as this field.
   *
   * @param other The other field
   * @return double Maximum relative pointwise absolute difference
   */
  double normalised_difference(Field &other);

  ~Field();
};

class ElectricField : public Field {

private:
  double phase(int n, double omega, double dt) override;

public:
  ElectricField() = default;
  ElectricField(int I_total, int J_total, int K_total)
      : Field(I_total, J_total, K_total){};

  /**
   * @brief Interpolates an E-field component to the centre of a Yee cell
   *
   * @param d Field component to interpolate
   * @param i,j,k Index (i,j,k) of the Yee cell to interpolate to the centre of
   * @return std::complex<double> The interpolated component value
   */
  std::complex<double> interpolate_to_centre_of(AxialDirection d,
                                                CellCoordinate cell) override;

  /**
   * @brief Interpolates the transverse electric components to the centre of Yee
   * cell i,j,k.
   *
   * Ex and Ey are interpolated. Ez is set to a placeholder (default) value.
   *
   * @param[in] cell Yee cell index
   * @param[out] x_at_centre,y_at_centre,z_at_centre Addresses to write
   * interpolated values for the x,y,z components (respectively)
   */
  void interpolate_transverse_electric_components(
          CellCoordinate cell, std::complex<double> *x_at_centre,
          std::complex<double> *y_at_centre,
          std::complex<double> *z_at_centre) override;
  /**
   * @brief Interpolates the transverse magnetic components to the centre of Yee
   * cell i,j,k.
   *
   * Ez is interpolated. Ex and Ey are set to a placeholder (default) values.
   *
   * @param[in] cell Yee cell index
   * @param[out] x_at_centre,y_at_centre,z_at_centre Addresses to write
   * interpolated values for the x,y,z components (respectively)
   */
  void interpolate_transverse_magnetic_components(
          CellCoordinate cell, std::complex<double> *x_at_centre,
          std::complex<double> *y_at_centre,
          std::complex<double> *z_at_centre) override;
};

class MagneticField : public Field {

private:
  double phase(int n, double omega, double dt) override;

public:
  MagneticField() = default;
  MagneticField(int I_total, int J_total, int K_total)
      : Field(I_total, J_total, K_total){};

  /**
   * @brief Interpolates an H-field component to the centre of a Yee cell
   *
   * @param d Field component to interpolate
   * @param cell Index (i,j,k) of the Yee cell to interpolate to the centre of
   * @return std::complex<double> The interpolated component value
   */
  std::complex<double> interpolate_to_centre_of(AxialDirection d,
                                                CellCoordinate cell) override;

  /**
   * @brief Interpolates the transverse electric components to the centre of Yee
   * cell i,j,k.
   *
   * Hz is interpolated. Hx and Hy are set to a placeholder (default) values.
   *
   * @param[in] cell Yee cell index
   * @param[out] x_at_centre,y_at_centre,z_at_centre Addresses to write
   * interpolated values for the x,y,z components (respectively)
   */
  void interpolate_transverse_electric_components(
          CellCoordinate cell, std::complex<double> *x_at_centre,
          std::complex<double> *y_at_centre,
          std::complex<double> *z_at_centre) override;
  /**
   * @brief Interpolates the transverse magnetic components to the centre of Yee
   * cell i,j,k.
   *
   * Hx and Hy are interpolated. Hz is set to a placeholder (default) value.
   *
   * @param[in] cell Yee cell index
   * @param[out] x_at_centre,y_at_centre,z_at_centre Addresses to write
   * interpolated values for the x,y,z components (respectively)
   */
  void interpolate_transverse_magnetic_components(
          CellCoordinate cell, std::complex<double> *x_at_centre,
          std::complex<double> *y_at_centre,
          std::complex<double> *z_at_centre) override;
};

/**
 * Structure to hold a field and allow saving it to a file
 */
class TDFieldExporter2D {
private:
  int nI = 0;//< Array size in the i direction
  int nK = 0;//< Array size in the k direction
public:
  mxArray *matlab_array = nullptr;
  double **array = nullptr;
  std::string folder_name = "";

  /**
   * Allocate the arrays to hold the field
   */
  void allocate(int _nI, int _nK);

  /**
   * Export/save a field
   *
   * @param F Field to save
   * @param stride Interval to compute the field component at
   * @param iteration Iteration number of the main loop. Used in the filename
   */
  void export_field(SplitField &F, int stride, int iteration) const;

  ~TDFieldExporter2D();
};
