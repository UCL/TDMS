#include "complex"
#include "simulation_parameters.h"


class Field{

public:
    std::complex<double> angular_norm = 0.;

    std::complex<double> ***x = nullptr;
    std::complex<double> ***y = nullptr;
    std::complex<double> ***z = nullptr;
};

class SplitField{

public:
    // Pointers (3D arrays) which hold the magnitude of the split field
    // component at each grid point (i, j, k)
    double ***xy = nullptr;
    double ***xz = nullptr;
    double ***yx = nullptr;
    double ***yz = nullptr;
    double ***zx = nullptr;
    double ***zy = nullptr;

    // To reconstruct the components we have e.g.: Ex = Exy + Exz multiplied by
    // a phase factor
};


class ElectricSplitField: public SplitField{};

class MagneticSplitField: public SplitField{};

class CurrentDensitySplitField: public SplitField{};

class ElectricField: public Field{

private:
    static std::complex<double> phasor_norm(double f, int n, double omega, double dt, int Nt);

public:
    double ft = 0.;  // TODO: an explanation of what this is

    void add_to_angular_norm(int n, int Nt, SimulationParameters &params);

};

class MagneticField: public Field{

private:
    static std::complex<double> phasor_norm(double f, int n, double omega, double dt, int Nt);

public:
    double ft = 0.;  // TODO: an explanation of what this is

    void add_to_angular_norm(int n, int Nt, SimulationParameters &params);

};
