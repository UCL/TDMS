/**
 * @file interpolate.h
 * @brief Interpolation of field values within FDTD grid.
 * 
 * The Yee cell specifies the 6 field components at different points in space.
 * Interpolation is required to calculate the 6 field components at the same
 * spatial position.
 */
#pragma once
#include "matrix.h"
#include "field.h"

/**
 * @brief Interpolate the electric field to the origin of the Yee cell
 *
 * @param[in] Ex_yee Steady state x component of electric field calculated at points in the Yee cell
 * @param[in] Ey_yee Steady state y component of electric field calculated at points in the Yee cell
 * @param[in] Ez_yee Steady state z component of electric field calculated at points in the Yee cell
 * @param[out] Ex Steady state x component of electric field interpolated to Yee cell origin
 * @param[out] Ey Steady state y component of electric field interpolated to Yee cell origin
 * @param[out] Ez Steady state z component of electric field interpolated to Yee cell origin
 * @param[in] I Number of elements in the i direction of the FDTD grid
 * @param[in] J Number of elements in the j direction of the FDTD grid
 * @param[in] K Number of elements in the k direction of the FDTD grid
 * @param[in] i_l Least i index into the FDTD grid to evaluate the field at. Should be >= 2
 * @param[in] i_u Greatest i index into the FDTD grid to evaluate the field at. Should be <= I-2
 * @param[in] j_l Least j index into the FDTD grid to evaluate the field at. Should be >= 2
 * @param[in] j_u Greatest j index into the FDTD grid to evaluate the field at. Should be <= J-2
 * @param[in] k_l Least k index into the FDTD grid to evaluate the field at. Should be >= 2
 * @param[in] k_u Greatest k index into the FDTD grid to evaluate the field at. Should be <= K-2
 */
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