#include <complex>
using namespace std;


void extractPhasorsPlane( double **iwave_lEx_Rbs, double **iwave_lEx_Ibs, double **iwave_lEy_Rbs, double **iwave_lEy_Ibs, 
			  double **iwave_lHx_Rbs, double **iwave_lHx_Ibs, double **iwave_lHy_Rbs, double **iwave_lHy_Ibs, 
			  double ***Exz, double ***Eyz, double ***Hxz, double ***Hyz,
			  double ***Exy, double ***Eyx, double ***Hxy, double ***Hyx,
			  int I_tot, int J_tot, int K1, int n, double omega, double dt, int Nt);

void extractPhasorsVolume(double ***ExR, double ***ExI, double ***EyR, double ***EyI, double ***EzR, double ***EzI,
			  double ***Exy, double ***Exz, double ***Eyx, double ***Eyz, double ***Ezx, double ***Ezy,
			  int i_l, int i_u, int j_l, int j_u, int k_l, int k_u, int n, double omega, double dt, int Nt);

void extractPhasorsVolumeH(double ***HxR, double ***HxI, double ***HyR, double ***HyI, double ***HzR, double ***HzI,
			  double ***Hxy, double ***Hxz, double ***Hyx, double ***Hyz, double ***Hzx, double ***Hzy,
			   int i_l, int i_u, int j_l, int j_u, int k_l, int k_u, int n, double omega, double dt, int Nt);

void initialiseDouble3DArray(double ***inArray, int i_lim, int j_lim, int k_lim);
void initialiseDouble2DArray(double **inArray, int i_lim, int j_lim);
double linearRamp(double t, double period, double rampwidth);
double complexAbs(complex<double> z);
double checkPhasorConvergence(double ***ExR, double ***ExI, double ***EyR, double ***EyI, double ***EzR, double ***EzI,
			    double ***ExR2, double ***ExI2, double ***EyR2, double ***EyI2, double ***EzR2, double ***EzI2,
			    double ***Exy, double ***Exz, double ***Eyx, double ***Eyz, double ***Ezx, double ***Ezy,
			      int i_l, int i_u, int j_l, int j_u, int k_l, int k_u, int n, double omega, double dt, int Nt);

void copyPhasors(double *ExR, double *ExI, double *EyR, double *EyI, double *EzR, double *EzI,
		 double *ExR2, double *ExI2, double *EyR2, double *EyI2, double *EzR2, double *EzI2,
		 int nelements);
void setGridLabels(double *x_labels_in , double * y_labels_in , double *z_labels_in ,
		   double *x_labels_out, double * y_labels_out, double *z_labels_out,
		   int i_l, int i_u, int j_l, int j_u, int k_l, int k_u);
//void mexFunction(int nlhs, mxArray *plhs[], int nrhs,const mxArray *prhs[]);
void extractPhasorsSurface( double **surface_EHr, double **surface_EHi,  
			    double ***Hxy, double ***Hxz, double ***Hyx, double ***Hyz, double ***Hzx, double ***Hzy,
			    double ***Exy, double ***Exz, double ***Eyx, double ***Eyz, double ***Ezx, double ***Ezy,
			    int **surface_vertices, int n_surface_vertices, int n, double omega, double dt, int Nt, int dimension );
void extractPhasorsSurface( double **surface_EHr, double **surface_EHi,  
			    double ***Hxy, double ***Hxz, double ***Hyx, double ***Hyz, double ***Hzx, double ***Hzy,
			    double ***Exy, double ***Exz, double ***Eyx, double ***Eyz, double ***Ezx, double ***Ezy,
			    int **surface_vertices, int n_surface_vertices, int n, double omega, double dt, int Nt, int dimension,int J_tot,int intmethod );
void extractPhasorsSurfaceNoInterpolation( double **surface_EHr, double **surface_EHi,  
			    double ***Hxy, double ***Hxz, double ***Hyx, double ***Hyz, double ***Hzx, double ***Hzy,
			    double ***Exy, double ***Exz, double ***Eyx, double ***Eyz, double ***Ezx, double ***Ezy,
					   int **surface_vertices, int n_surface_vertices, int n, double omega, double dt, int Nt, int dimension,int J_tot );
int is_dispersive(unsigned char ***materials,double *gamma, double dt, int I_tot, int J_tot, int K_tot);
void allocate_auxilliary_mem(int I_tot, int J_tot, int K_tot,
			     double ****Exy, double ****Exz, 
			     double ****Eyx, double ****Eyz,  
			     double ****Ezx, double ****Ezy,
			     double ****Jxy, double ****Jxz, 
			     double ****Jyx, double ****Jyz,  
			     double ****Jzx, double ****Jzy,
			     double ****Jxy2, double ****Jxz2, 
			     double ****Jyx2, double ****Jyz2,  
			     double ****Jzx2, double ****Jzy2);

void destroy_auxilliary_mem(int I_tot, int J_tot, int K_tot,
			     double ****Exy, double ****Exz, 
			     double ****Eyx, double ****Eyz,  
			     double ****Ezx, double ****Ezy,
			     double ****Jxy, double ****Jxz, 
			     double ****Jyx, double ****Jyz,  
			    double ****Jzx, double ****Jzy,
			     double ****Jxy2, double ****Jxz2, 
			     double ****Jyx2, double ****Jyz2,  
			     double ****Jzx2, double ****Jzy2);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,const mxArray *prhs[]);
void extractPhasorENorm(complex<double> *Enorm, double ft, int n, double omega, double dt, int Nt);
void extractPhasorHNorm(complex<double> *Hnorm, double ft, int n, double omega, double dt, int Nt);
void normaliseVolume(double ***ExR, double ***ExI, double ***EyR, double ***EyI, double ***EzR, double ***EzI,
		     int i_l, int i_u, int j_l, int j_u, int k_l, int k_u,  complex<double> norm);
void normaliseSurface( double **surface_EHr, double **surface_EHi ,
		       int **surface_vertices, int n_surface_vertices,  complex<double> Enorm , complex<double> Hnorm );
void normaliseVertices( double **EHr, double **EHi ,
			int **vertices, int nvertices,
			int *components, int ncomponents,
			complex<double> Enorm , complex<double> Hnorm );
int is_conductive(double *rho_x, double *rho_y, double *rho_z, int I_tot, int J_tot, int K_tot);
int is_dispersive_ml(double *ml_gamma, int K_tot);
void allocate_auxilliary_mem_conductive(int I_tot, int J_tot, int K_tot,
					double ****Jxy, double ****Jxz, 
					double ****Jyx, double ****Jyz,  
					double ****Jzx, double ****Jzy);
void destroy_auxilliary_mem_conductive(int I_tot, int J_tot, int K_tot,
				       double ****Jxy, double ****Jxz, 
				       double ****Jyx, double ****Jyz,  
				       double ****Jzx, double ****Jzy);
void extractPhasorsVertices( double **EHr, double **EHi,  
			     double ***Hxy, double ***Hxz, double ***Hyx, double ***Hyz, double ***Hzx, double ***Hzy,
			     double ***Exy, double ***Exz, double ***Eyx, double ***Eyz, double ***Ezx, double ***Ezy,
			     int **vertices, int nvertices, int *components, int ncomponents,
			     int n, double omega, double dt, int Nt, int dimension,int J_tot,int intmethod );
int find(int *a, int na, int b);
