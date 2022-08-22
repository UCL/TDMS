
class SimulationParameters{

public:
    SimulationParameters();

    double       omega_an = 0.0;          // Angular ω
    unsigned int Nsteps   = 0;            // Number of simulation steps
    double       dt       = 0.0;          // Time step
    bool         has_tdfdir = false;      // Is the tdfdir (time domain field directory) defined?
};
