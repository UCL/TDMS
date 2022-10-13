/**
 * @file field.h
 * @brief Classes for the electric and magnetic (split) fields on a grid.
 */
#pragma once
#include <complex>
#include <stdexcept>
#include "arrays.h"
#include "dimensions.h"
#include "mat_io.h"
#include "simulation_parameters.h"
#include "utils.h"


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
 */
class Grid{

public:
    int I_tot = 0;
    int J_tot = 0;
    int K_tot = 0;

    /**
     * Maximum value out of I_tot, J_tot and K_tot
     * @return value
     */
    int max_IJK_tot() const {return max(I_tot, J_tot, K_tot); };
};

class SplitFieldComponent: public Tensor3D<double>{
public:
  int n_threads = 1;             // Number of threads this component was chunked with
  fftw_plan* plan_f = nullptr;  // Forward fftw plan
  fftw_plan* plan_b = nullptr;  // Backward fftw plan

  void initialise_from_matlab(double*** tensor, Dimensions &dims);

  /**
   * Initialise a vector of 1d discrete Fourier transform plans
   * @param n_threads Number of threads that will be used
   * @param size Length of the vector
   * @param eh_vec // TODO: what is this?
   */
  void initialise_fftw_plan(int n_threads, int size, EHVec &eh_vec);

  ~SplitFieldComponent();
};

/**
 * A split field defined over a grid.
 * To reconstruct the components we have e.g.: Ex = Exy + Exz multiplied by
 * a phase factor
 */
class SplitField : public Grid{
protected:
  virtual int delta_n() = 0;  // TODO: no idea what this is or why it's needed

public:
    // Pointers (3D arrays) which hold the magnitude of the split field
    // component at each grid point (i, j, k)
    SplitFieldComponent xy;
    SplitFieldComponent xz;
    SplitFieldComponent yx;
    SplitFieldComponent yz;
    SplitFieldComponent zx;
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
    void allocate_and_zero(){
      allocate();
      zero();
    }

    /**
     * Initialise the fftw plans for all components
     * @param n_threads Number of threads to split over
     * @param eh_vec // TODO
     */
    void initialise_fftw_plan(int n_threads, EHVec &eh_vec);
};

class ElectricSplitField: public SplitField{
protected:
  int delta_n() override { return -1; };  // TODO: no idea what this is or why it's needed

public:
    ElectricSplitField() = default;

    /**
     * Constructor of the field with a defined size in the x, y, z Cartesian
     * dimensions
     */
    ElectricSplitField(int I_total, int J_total, int K_total) :
            SplitField(I_total, J_total, K_total){};

    /**
     * @brief Interpolate the x-component to the centre of cell (i,j,k)
     * 
     * @param i,j,k Index of the cell to interpolate to the centre of 
     * @return double Value of the Ex component at the cell centre
     */
    double interpolate_x_to_centre(int i, int j, int k);
    /**
     * @brief Interpolate the y-component to the centre of cell (i,j,k)
     * 
     * @param i,j,k Index of the cell to interpolate to the centre of 
     * @return double Value of the Ey component at the cell centre
     */
    double interpolate_y_to_centre(int i, int j, int k);
    /**
     * @brief Interpolate the z-component to the centre of cell (i,j,k)
     * 
     * @param i,j,k Index of the cell to interpolate to the centre of 
     * @return double Value of the Ez component at the cell centre
     */
    double interpolate_z_to_centre(int i, int j, int k);
};

class MagneticSplitField: public SplitField{
protected:
  int delta_n() override { return 0; };  // TODO: no idea what this is or why it's needed

public:
    MagneticSplitField() = default;

    /**
     * Constructor of the field with a defined size in the x, y, z Cartesian
     * dimensions
     */
    MagneticSplitField(int I_total, int J_total, int K_total) :
            SplitField(I_total, J_total, K_total){};

    /**
     * @brief Interpolate the x-component to the centre of cell (i,j,k)
     * 
     * @param i,j,k Index of the cell to interpolate to the centre of 
     * @return double Value of the Hx component at the cell centre
     */
    double interpolate_x_to_centre(int i, int j, int k);
    /**
     * @brief Interpolate the y-component to the centre of cell (i,j,k)
     * 
     * @param i,j,k Index of the cell to interpolate to the centre of 
     * @return double Value of the Hy component at the cell centre
     */
    double interpolate_y_to_centre(int i, int j, int k);
    /**
     * @brief Interpolate the z-component to the centre of cell (i,j,k)
     * 
     * @param i,j,k Index of the cell to interpolate to the centre of 
     * @return double Value of the Hz component at the cell centre
     */
    double interpolate_z_to_centre(int i, int j, int k);
};

class CurrentDensitySplitField: public SplitField{
protected:
  int delta_n() override { return 0; }

public:
    CurrentDensitySplitField() = default;

    /**
     * Constructor of the field with a defined size in the x, y, z Cartesian
     * dimensions
     */
    CurrentDensitySplitField(int I_total, int J_total, int K_total) :
            SplitField(I_total, J_total, K_total){};
};

/**
 * A complex field defined over a grid. Has real and imaginary components
 * at each (x, y, z) grid point
 */
class Field : public Grid{

public:
  double ft = 0.;  // TODO: an explanation of what this is

  std::complex<double> angular_norm = 0.;

  // TODO: this is likely better as a set of complex arrays - use Tensor3d<std::complex<double>>
  XYZTensor3D real;
  XYZTensor3D imag;

  /**
     * Upper (u) and lower (l) indices in the x,y,z directions. e.g.
     * il is the first non-pml cell in the i direction and iu the last in the corresponding split
     * field grid
     */
  int il = 0, iu = 0, jl = 0, ju = 0, kl = 0, ku = 0;

  /**
   * @brief Normalise all field values by the angular_norm
   * 
   * If the field values are stored as complex numbers, this amounts to dividing each stored field value by angular_norm.
   * 
   * If the field values are stored in their real and imaginary parts separately, the expressions are more complicated.
   */
  void normalise_volume();

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
   * Set the phasors for this field, given a split field. Result gives field according to the
   * exp(-iwt) convention
   * @param F
   * @param n
   * @param omega
   * @param dt
   * @param Nt
   */
  void set_phasors(SplitField &F, int n, double omega, double dt, int Nt);

  // TODO: Docstring
  void add_to_angular_norm(int n, int Nt, SimulationParameters &params);

  // TODO: Docstring
  std::complex<double> phasor_norm(double f, int n, double omega, double dt, int Nt);

  virtual double phase(int n, double omega, double dt) = 0;
  virtual std::complex<double> interpolate_x_to_centre(int i, int j, int k) = 0;
  virtual std::complex<double> interpolate_y_to_centre(int i, int j, int k) = 0;
  virtual std::complex<double> interpolate_z_to_centre(int i, int j, int k) = 0;

  /**
   * @brief Interpolate the Field to the centre of all Yee cells within the provided index range.
   * 
   * Interpolation is performed across all cells if no cutoffs are specified.
   * 
   * @param x_out,y_out,z_out Arrays in which to store the interpolated x, y, z components respectively. 
   * @param i_lower_cutoff,j_lower_cutoff,k_lower_cutoff First cell index in the x,y,z direction (respectively) to interpolate to the centre of.
   * @param i_upper_cutoff,j_upper_cutoff,k_upper_cutoff Last cell index in the x,y,z direction (respectively) to interpolate to the centre of.
   */
  void interpolate_across_range(mxArray **x_out, mxArray **y_out, mxArray **z_out,
                                int i_lower_cutoff, int i_upper_cutoff, int j_lower_cutoff,
                                int j_upper_cutoff, int k_lower_cutoff, int k_upper_cutoff);
  void interpolate_across_range(mxArray **x_out, mxArray **y_out, mxArray **z_out);

  /**
   * Set the values of all components in this field from another, equally sized field
   */
   void set_values_from(Field &other);

  ~Field();
};

class ElectricField: public Field{

private:
  double phase(int n, double omega, double dt) override;

public:
  ElectricField() = default;
  ElectricField(int I_total, int J_total, int K_total) : Field(I_total, J_total, K_total){};

  /**
   * @brief Interpolate the Ex component to the centre of cell (i,j,k)
   * 
   * @param i,j,k Cell to interpolate to the centre of
   * @return double Interpolated value
   */
  std::complex<double> interpolate_x_to_centre(int i, int j, int k);
  /**
   * @brief Interpolate the Ey component to the centre of cell (i,j,k)
   * 
   * @param i,j,k Cell to interpolate to the centre of 
   * @return double Interpolated value
   */
  std::complex<double> interpolate_y_to_centre(int i, int j, int k);
  /**
   * @brief Interpolate the Ez component to the centre of cell (i,j,k)
   * 
   * @param i,j,k Cell to interpolate to the centre of 
   * @return double Interpolated value
   */
  std::complex<double> interpolate_z_to_centre(int i, int j, int k);
};

class MagneticField: public Field{

private:
  double phase(int n, double omega, double dt) override;

public:
  MagneticField() = default;
  MagneticField(int I_total, int J_total, int K_total) : Field(I_total, J_total, K_total){};

  /**
   * @brief Interpolate the Hx component to the centre of cell (i,j,k)
   * 
   * @param i,j,k Cell to interpolate to the centre of
   * @return double Interpolated value
   */
  std::complex<double> interpolate_x_to_centre(int i, int j, int k);
  /**
   * @brief Interpolate the Hy component to the centre of cell (i,j,k)
   * 
   * @param i,j,k Cell to interpolate to the centre of 
   * @return double Interpolated value
   */
  std::complex<double> interpolate_y_to_centre(int i, int j, int k);
  /**
   * @brief Interpolate the Hz component to the centre of cell (i,j,k)
   * 
   * @param i,j,k Cell to interpolate to the centre of 
   * @return double Interpolated value
   */
  std::complex<double> interpolate_z_to_centre(int i, int j, int k);
};

/**
 * Structure to hold a field and allow saving it to a file
 */
class TDFieldExporter2D{

public:
  mxArray* matlab_array = nullptr;
  double** array = nullptr;
  const char* folder_name = nullptr;

  /**
   * Allocate the arrays to hold the field
   */
  void allocate(int nI, int nJ);

  /**
   * Export/save a field
   *
   * @param F Field to save
   * @param stride Interval to compute the field component at
   * @param iteration Iteration number of the main loop. Used in the filename
   */
  void export_field(SplitField& F, int stride, int iteration) const;

  ~TDFieldExporter2D();
};
