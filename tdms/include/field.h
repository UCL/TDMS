#include <complex>
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
 */
class Grid{

protected:
    int I_tot = 0;
    int J_tot = 0;
    int K_tot = 0;

    /**
     * Is this grid 0 x 0 x 0?
     * @return True if there are no components
     */
    bool has_no_elements() const {return I_tot == 0 && J_tot == 0 && K_tot == 0;};
};


struct xyz_arrays{
    double ***x;
    double ***y;
    double ***z;
};


/**
 * A complex vector field defined over a grid. Has three real and imaginary vector components
 * at each (x, y, z) grid point
 */
class Field : public Grid{

public:
    double ft = 0.;  // TODO: an explanation of what this is

    std::complex<double> angular_norm = 0.;

    // TODO: this is likely better as a set of complex arrays
    xyz_arrays real = xyz_arrays{nullptr, nullptr, nullptr};
    xyz_arrays imag = xyz_arrays{nullptr, nullptr, nullptr};

    void add_to_angular_norm(int n, int Nt, SimulationParameters &params);

    virtual std::complex<double> phasor_norm(double f, int n, double omega, double dt, int Nt) = 0;
};


/**
 * A split field defined over a grid.
 * To reconstruct the components we have e.g.: Ex = Exy + Exz multiplied by
 * a phase factor
 */
class SplitField : public Grid{

public:
    // Pointers (3D arrays) which hold the magnitude of the split field
    // component at each grid point (i, j, k)
    double ***xy = nullptr;
    double ***xz = nullptr;
    double ***yx = nullptr;
    double ***yz = nullptr;
    double ***zx = nullptr;
    double ***zy = nullptr;

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
     * Destructor. Frees all the allocated memory
     */
    ~SplitField();
};


class ElectricSplitField: public SplitField{

public:
    ElectricSplitField() = default;

    /**
     * Constructor of the field with a defined size in the x, y, z Cartesian
     * dimensions
     */
    ElectricSplitField(int I_total, int J_total, int K_total) :
            SplitField(I_total, J_total, K_total){};
};

class MagneticSplitField: public SplitField{

public:
    MagneticSplitField() = default;

    /**
     * Constructor of the field with a defined size in the x, y, z Cartesian
     * dimensions
     */
    MagneticSplitField(int I_total, int J_total, int K_total) :
            SplitField(I_total, J_total, K_total){};
};

class CurrentDensitySplitField: public SplitField{

public:
    CurrentDensitySplitField() = default;

    /**
     * Constructor of the field with a defined size in the x, y, z Cartesian
     * dimensions
     */
    CurrentDensitySplitField(int I_total, int J_total, int K_total) :
            SplitField(I_total, J_total, K_total){};
};

class ElectricField: public Field{

private:
    std::complex<double> phasor_norm(double f, int n, double omega, double dt, int Nt) override;
};

class MagneticField: public Field{

private:
    std::complex<double> phasor_norm(double f, int n, double omega, double dt, int Nt) override;
};
