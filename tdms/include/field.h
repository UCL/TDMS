#include <complex>
#include "simulation_parameters.h"


class Field{

public:
    std::complex<double> angular_norm = 0.;
};


class ElectricField: public Field{

private:
    static std::complex<double> phasor_norm(double f, int n, double omega, double dt, int Nt);

public:
    void add_to_angular_norm(double f, int n, int Nt, SimulationParameters &params);

};


class MagneticField: public Field{

private:
    static std::complex<double> phasor_norm(double f, int n, double omega, double dt, int Nt);

public:
    void add_to_angular_norm(double f, int n, int Nt, SimulationParameters &params);

};

