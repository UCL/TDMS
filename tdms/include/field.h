#include "complex"
#include "simulation_parameters.h"


class SplitField{

public:
    std::complex<double> angular_norm = 0.;

    double ***xy = nullptr;  // xy component of the split field (3D array)
    double ***xz = nullptr;  // xz
    double ***yx = nullptr;  // yx
    double ***yz = nullptr;  // yz
    double ***zx = nullptr;  // zx
    double ***zy = nullptr;  // zy

    // To reconstruct the components we have e.g.: Ex = Exy + Exz
};


class ElectricField: public SplitField{

private:
    static std::complex<double> phasor_norm(double f, int n, double omega, double dt, int Nt);

public:
    double ft = 0.;  // TODO: an explanation of what this is

    void add_to_angular_norm(int n, int Nt, SimulationParameters &params);

};


class MagneticField: public SplitField{

private:
    static std::complex<double> phasor_norm(double f, int n, double omega, double dt, int Nt);

public:
    double ft = 0.;  // TODO: an explanation of what this is

    void add_to_angular_norm(int n, int Nt, SimulationParameters &params);

};


class CurrentDensityField: public SplitField{

};
