#include "matrix.h"
#include "field.h"

void interpolateFieldCentralE( double ***Ex_yee, double ***Ey_yee, double ***Ez_yee,
			      double ***Ex    , double ***Ey    , double ***Ez    ,
                              int       I     , int       J     , int       K     ,
			      int i_l, int i_u, int j_l, int j_u, int k_l, int k_u);
void interpolateFieldCentralE_TE( double ***Ex_yee, double ***Ey_yee, double ***Ez_yee,
			      double ***Ex    , double ***Ey    , double ***Ez    ,
                              int       I     , int       J     , int       K     ,
			      int i_l, int i_u, int j_l, int j_u, int k_l, int k_u);
void mxInterpolateFieldCentralE( mxArray *Ex_yee , mxArray *Ey_yee , mxArray *Ez_yee,
			        mxArray **Ex    , mxArray **Ey    , mxArray **Ez    ,
				int i_l, int i_u, int j_l, int j_u, int k_l, int k_u);
void mxInterpolateFieldCentralE_TE( mxArray *Ex_yee , mxArray *Ey_yee , mxArray *Ez_yee,
			        mxArray **Ex    , mxArray **Ey    , mxArray **Ez    ,
				int i_l, int i_u, int j_l, int j_u, int k_l, int k_u);
void mxInterpolateFieldCentralE_TM( mxArray *Ex_yee , mxArray *Ey_yee , mxArray *Ez_yee,
			        mxArray **Ex    , mxArray **Ey    , mxArray **Ez    ,
				int i_l, int i_u, int j_l, int j_u, int k_l, int k_u);
void interpolateFieldCentralH( double ***Hx_yee, double ***Hy_yee, double ***Hz_yee,
			      double ***Hx    , double ***Hy    , double ***Hz    ,
                              int       I     , int       J     , int       K     ,
			      int i_l, int i_u, int j_l, int j_u, int k_l, int k_u);
void interpolateFieldCentralH_TE( double ***Hx_yee, double ***Hy_yee, double ***Hz_yee,
			      double ***Hx    , double ***Hy    , double ***Hz    ,
                              int       I     , int       J     , int       K     ,
			      int i_l, int i_u, int j_l, int j_u, int k_l, int k_u);
void mxInterpolateFieldCentralH( mxArray *Hx_yee , mxArray *Hy_yee , mxArray *Hz_yee,
			        mxArray **Hx    , mxArray **Hy    , mxArray **Hz    ,
				int i_l, int i_u, int j_l, int j_u, int k_l, int k_u);
void mxInterpolateFieldCentralH_TE( mxArray *Hx_yee , mxArray *Hy_yee , mxArray *Hz_yee,
			        mxArray **Hx    , mxArray **Hy    , mxArray **Hz    ,
				int i_l, int i_u, int j_l, int j_u, int k_l, int k_u);
void mxInterpolateFieldCentralH_TM( mxArray *Hx_yee , mxArray *Hy_yee , mxArray *Hz_yee,
			        mxArray **Hx    , mxArray **Hy    , mxArray **Hz    ,
				int i_l, int i_u, int j_l, int j_u, int k_l, int k_u);
void interpolateTimeDomainFieldCentralE(  SplitFieldComponent& Exy, SplitFieldComponent& Exz, SplitFieldComponent& Eyx, SplitFieldComponent& Eyz, SplitFieldComponent& Ezx, SplitFieldComponent& Ezy,
			       int i, int j, int k,
					 double *Ex, double *Ey, double *Ez);
void interpolateTimeDomainFieldCentralEBandLimited(  SplitFieldComponent& Exy, SplitFieldComponent& Exz, SplitFieldComponent& Eyx, SplitFieldComponent& Eyz, SplitFieldComponent& Ezx, SplitFieldComponent& Ezy,
					 int i, int j, int k,
						    double *Ex, double *Ey, double *Ez);
void interpolateTimeDomainFieldCentralE_2Dy(  SplitFieldComponent& Exy, SplitFieldComponent& Exz, SplitFieldComponent& Eyx, SplitFieldComponent& Eyz, SplitFieldComponent& Ezx, SplitFieldComponent& Ezy,
			       int i, int j, int k,
			       double *Ex, double *Ey, double *Ez);
void interpolateTimeDomainFieldCentralE_TE(  SplitFieldComponent& Exy, SplitFieldComponent& Exz, SplitFieldComponent& Eyx, SplitFieldComponent& Eyz, SplitFieldComponent& Ezx, SplitFieldComponent& Ezy,
			       int i, int j, int k,
			       double *Ex, double *Ey, double *Ez);
void interpolateTimeDomainFieldCentralE_TM(  SplitFieldComponent& Exy, SplitFieldComponent& Exz, SplitFieldComponent& Eyx, SplitFieldComponent& Eyz, SplitFieldComponent& Ezx, SplitFieldComponent& Ezy,
			       int i, int j, int k,
			       double *Ex, double *Ey, double *Ez);
void interpolateTimeDomainFieldCentralH( SplitFieldComponent& Hxy, SplitFieldComponent& Hxz, SplitFieldComponent& Hyx, SplitFieldComponent& Hyz, SplitFieldComponent& Hzx, SplitFieldComponent& Hzy,
					int i, int j, int k,
					double *Hx, double *Hy, double *Hz);
void interpolateTimeDomainFieldCentralH_2Dy( SplitFieldComponent& Hxy, SplitFieldComponent& Hxz, SplitFieldComponent& Hyx, SplitFieldComponent& Hyz, SplitFieldComponent& Hzx, SplitFieldComponent& Hzy,
					int i, int j, int k,
					double *Hx, double *Hy, double *Hz);
void interpolateTimeDomainFieldCentralH_TE( SplitFieldComponent& Hxy, SplitFieldComponent& Hxz, SplitFieldComponent& Hyx, SplitFieldComponent& Hyz, SplitFieldComponent& Hzx, SplitFieldComponent& Hzy,
					int i, int j, int k,
					double *Hx, double *Hy, double *Hz);
void interpolateTimeDomainFieldCentralH_TM( SplitFieldComponent& Hxy, SplitFieldComponent& Hxz, SplitFieldComponent& Hyx, SplitFieldComponent& Hyz, SplitFieldComponent& Hzx, SplitFieldComponent& Hzy,
					int i, int j, int k,
					double *Hx, double *Hy, double *Hz);
void interpolateTimeDomainFieldCentralHBandLimited( SplitFieldComponent& Hxy, SplitFieldComponent& Hxz, SplitFieldComponent& Hyx, SplitFieldComponent& Hyz, SplitFieldComponent& Hzx, SplitFieldComponent& Hzy,
					int i, int j, int k,
						   double *Hx, double *Hy, double *Hz);