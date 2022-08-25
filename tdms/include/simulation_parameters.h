
class SimulationParameters{

public:
    SimulationParameters();

    double       omega_an      = 0.0;      // Angular ω
    unsigned int Nsteps        = 0;        // Number of simulation steps
    double       dt            = 0.0;      // Time step
    bool         has_tdfdir    = false;    // Is the tdfdir (time domain field directory) defined?
    bool         is_multilayer = false;    // Is this simulation of a multilayer?
    bool         is_disp_ml    = false;    // Is the material dispersive?
    double       to_l          = 0.0;      // time delay of pulse
    double       hwhm          = 0.0;      // hwhm of pulse
};
