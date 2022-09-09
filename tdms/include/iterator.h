#include "arrays.h"
#include "complex"
#include "field.h"
#include "grid_labels.h"


void extractPhasorsPlane( double **iwave_lEx_Rbs, double **iwave_lEx_Ibs, double **iwave_lEy_Rbs, double **iwave_lEy_Ibs, 
			  double **iwave_lHx_Rbs, double **iwave_lHx_Ibs, double **iwave_lHy_Rbs, double **iwave_lHy_Ibs, 
        ElectricSplitField &E, MagneticSplitField &H,
			  int I_tot, int J_tot, int K1, int n, double omega, double dt, int Nt);

void initialiseDouble3DArray(double ***inArray, int i_lim, int j_lim, int k_lim);

void initialiseDouble2DArray(double **inArray, int i_lim, int j_lim);

double linearRamp(double t, double period, double rampwidth);

double checkPhasorConvergence(ElectricField &A, ElectricField &B);

void copyPhasors(ElectricField &from, ElectricField &to, int nelements);

void setGridLabels(GridLabels &input_labels, GridLabels &output_labels,
                   int i_l, int i_u, int j_l, int j_u, int k_l, int k_u);

void extractPhasorsSurface(double **surface_EHr, double **surface_EHi, ElectricSplitField &E, MagneticSplitField &H,
			    int **surface_vertices, int n_surface_vertices, int n, double omega, int Nt, int J_tot, SimulationParameters &params);

void extractPhasorsSurfaceNoInterpolation(double **surface_EHr, double **surface_EHi, ElectricSplitField &E, MagneticSplitField &H,
					   int **surface_vertices, int n_surface_vertices, int n, double omega, int Nt, int J_tot, SimulationParameters &params);

bool is_dispersive(unsigned char ***materials,double *gamma, double dt, int I_tot, int J_tot, int K_tot);

void extractPhasorENorm(std::complex<double> *Enorm, double ft, int n, double omega, double dt, int Nt);

void extractPhasorHNorm(std::complex<double> *Hnorm, double ft, int n, double omega, double dt, int Nt);

void normaliseSurface( double **surface_EHr, double **surface_EHi ,
		       int **surface_vertices, int n_surface_vertices,  std::complex<double> Enorm , std::complex<double> Hnorm );

void normaliseVertices( double **EHr, double **EHi, ComplexAmplitudeSample &campssample, std::complex<double> Enorm , std::complex<double> Hnorm );

void update_EH(double **EHr, double **EHi, int vindex, int idx, std::complex<double> &phase_term, double &value);

bool is_conductive(const XYZVectors &rho, int I_tot, int J_tot, int K_tot);

bool is_dispersive_ml(const DispersiveMultiLayer &ml, int K_tot);

void extractPhasorsVertices(double **EHr, double **EHi, ElectricSplitField &E, MagneticSplitField &H,
                            ComplexAmplitudeSample &campssample, int n, double omega, double dt, int Nt,
                            int dimension,int J_tot,int intmethod );
