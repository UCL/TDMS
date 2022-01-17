double interp1(double v1, double v2, double v3, double v4);
double interp2(double v1, double v2, double v3, double v4);
double interp3(double v1, double v2, double v3, double v4);
int interpolateFieldCentralE( double ***Ex_yee, double ***Ey_yee, double ***Ez_yee,
			      double ***Ex    , double ***Ey    , double ***Ez    ,
                              int       I     , int       J     , int       K     ,
			      int i_l, int i_u, int j_l, int j_u, int k_l, int k_u);
int interpolateFieldCentralE_TE( double ***Ex_yee, double ***Ey_yee, double ***Ez_yee,
			      double ***Ex    , double ***Ey    , double ***Ez    ,
                              int       I     , int       J     , int       K     ,
			      int i_l, int i_u, int j_l, int j_u, int k_l, int k_u);
int interpolateFieldCentralE_TH( double ***Ex_yee, double ***Ey_yee, double ***Ez_yee,
			      double ***Ex    , double ***Ey    , double ***Ez    ,
                              int       I     , int       J     , int       K     ,
			      int i_l, int i_u, int j_l, int j_u, int k_l, int k_u);
int mxInterpolateFieldCentralE( mxArray *Ex_yee , mxArray *Ey_yee , mxArray *Ez_yee,
			        mxArray **Ex    , mxArray **Ey    , mxArray **Ez    ,
				int i_l, int i_u, int j_l, int j_u, int k_l, int k_u);
int mxInterpolateFieldCentralE_TE( mxArray *Ex_yee , mxArray *Ey_yee , mxArray *Ez_yee,
			        mxArray **Ex    , mxArray **Ey    , mxArray **Ez    ,
				int i_l, int i_u, int j_l, int j_u, int k_l, int k_u);
int mxInterpolateFieldCentralE_TM( mxArray *Ex_yee , mxArray *Ey_yee , mxArray *Ez_yee,
			        mxArray **Ex    , mxArray **Ey    , mxArray **Ez    ,
				int i_l, int i_u, int j_l, int j_u, int k_l, int k_u);
int interpolateFieldCentralH( double ***Hx_yee, double ***Hy_yee, double ***Hz_yee,
			      double ***Hx    , double ***Hy    , double ***Hz    ,
                              int       I     , int       J     , int       K     ,
			      int i_l, int i_u, int j_l, int j_u, int k_l, int k_u);
int interpolateFieldCentralH_TE( double ***Hx_yee, double ***Hy_yee, double ***Hz_yee,
			      double ***Hx    , double ***Hy    , double ***Hz    ,
                              int       I     , int       J     , int       K     ,
			      int i_l, int i_u, int j_l, int j_u, int k_l, int k_u);
int interpolateFieldCentralH_TH( double ***Hx_yee, double ***Hy_yee, double ***Hz_yee,
			      double ***Hx    , double ***Hy    , double ***Hz    ,
                              int       I     , int       J     , int       K     ,
			      int i_l, int i_u, int j_l, int j_u, int k_l, int k_u);
int mxInterpolateFieldCentralH( mxArray *Hx_yee , mxArray *Hy_yee , mxArray *Hz_yee,
			        mxArray **Hx    , mxArray **Hy    , mxArray **Hz    ,
				int i_l, int i_u, int j_l, int j_u, int k_l, int k_u);
int mxInterpolateFieldCentralH_TE( mxArray *Hx_yee , mxArray *Hy_yee , mxArray *Hz_yee,
			        mxArray **Hx    , mxArray **Hy    , mxArray **Hz    ,
				int i_l, int i_u, int j_l, int j_u, int k_l, int k_u);
int mxInterpolateFieldCentralH_TM( mxArray *Hx_yee , mxArray *Hy_yee , mxArray *Hz_yee,
			        mxArray **Hx    , mxArray **Hy    , mxArray **Hz    ,
				int i_l, int i_u, int j_l, int j_u, int k_l, int k_u);
int interpolateTimeDomainFieldCentralE(  double ***Exy, double ***Exz, double ***Eyx, double ***Eyz, double ***Ezx, double ***Ezy,
			       int i, int j, int k,
					 double *Ex, double *Ey, double *Ez);
int interpolateTimeDomainFieldCentralEBandLimited(  double ***Exy, double ***Exz, double ***Eyx, double ***Eyz, double ***Ezx, double ***Ezy,
					 int i, int j, int k,
						    double *Ex, double *Ey, double *Ez);
int interpolateTimeDomainFieldCentralE_2Dy(  double ***Exy, double ***Exz, double ***Eyx, double ***Eyz, double ***Ezx, double ***Ezy,
			       int i, int j, int k,
			       double *Ex, double *Ey, double *Ez);
int interpolateTimeDomainFieldCentralE_TE(  double ***Exy, double ***Exz, double ***Eyx, double ***Eyz, double ***Ezx, double ***Ezy,
			       int i, int j, int k,
			       double *Ex, double *Ey, double *Ez);
int interpolateTimeDomainFieldCentralE_TM(  double ***Exy, double ***Exz, double ***Eyx, double ***Eyz, double ***Ezx, double ***Ezy,
			       int i, int j, int k,
			       double *Ex, double *Ey, double *Ez);
int interpolateTimeDomainFieldCentralH( double ***Hxy, double ***Hxz, double ***Hyx, double ***Hyz, double ***Hzx, double ***Hzy,
					int i, int j, int k,
					double *Hx, double *Hy, double *Hz);
int interpolateTimeDomainFieldCentralH_2Dy( double ***Hxy, double ***Hxz, double ***Hyx, double ***Hyz, double ***Hzx, double ***Hzy,
					int i, int j, int k,
					double *Hx, double *Hy, double *Hz);
int interpolateTimeDomainFieldCentralH_TE( double ***Hxy, double ***Hxz, double ***Hyx, double ***Hyz, double ***Hzx, double ***Hzy,
					int i, int j, int k,
					double *Hx, double *Hy, double *Hz);
int interpolateTimeDomainFieldCentralH_TM( double ***Hxy, double ***Hxz, double ***Hyx, double ***Hyz, double ***Hzx, double ***Hzy,
					int i, int j, int k,
					double *Hx, double *Hy, double *Hz);
int interpolateTimeDomainFieldCentralHBandLimited( double ***Hxy, double ***Hxz, double ***Hyx, double ***Hyz, double ***Hzx, double ***Hzy,
					int i, int j, int k,
						   double *Hx, double *Hy, double *Hz);
