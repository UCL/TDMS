/******************************************************************************
 *  Application.:  interpolation of field values within FDTD grid
 *
 *  Description.:  The Yee cell specifies the 6 field components at
 *                 different points in space. Interpolation is required
 *                 to calculate the 6 field components at the same 
 *                 spatial position.
 *****************************************************************************/
#include <stdexcept>
#include "interpolate.h"
#include "interpolation_methods.h"
#include "matlabio.h"

using namespace std;

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
			      int i_l, int i_u, int j_l, int j_u, int k_l, int k_u){
  
  int i,j,k;
 
  //first check that the limits of field extraction of within range
  checkInterpolationPoints(i_l, i_u, j_l, j_u, k_l, k_u, I, J, K);

  //fprintf(stderr, "interpolateFieldCentralE 01\n");
  if(j_u<j_l) {
    //fprintf(stderr, "interpolateFieldCentralE 02\n");
    j=0;
    for(k=k_l;k<=k_u;k++)
      for(i=i_l;i<=i_u;i++) {
	      Ex[k-k_l][0][i-i_l] = interp1(Ex_yee[k][j][i-2], Ex_yee[k][j][i-1], Ex_yee[k][j][i], Ex_yee[k][j][i+1]);
	      Ey[k-k_l][0][i-i_l] = Ey_yee[k][0][i];
	      Ez[k-k_l][0][i-i_l] = interp1(Ez_yee[k-2][j][i], Ez_yee[k-1][j][i], Ez_yee[k][j][i], Ez_yee[k+1][j][i]);
	    }
  //fprintf(stderr, "interpolateFieldCentralE 03\n");
  }
  else
    for(k=k_l;k<=k_u;k++)
      for(j=j_l;j<=j_u;j++)
	      for(i=i_l;i<=i_u;i++){
	        Ex[k-k_l][j-j_l][i-i_l] = interp1(Ex_yee[k][j][i-2], Ex_yee[k][j][i-1], Ex_yee[k][j][i], Ex_yee[k][j][i+1]);
	        Ey[k-k_l][j-j_l][i-i_l] = interp1(Ey_yee[k][j-2][i], Ey_yee[k][j-1][i], Ey_yee[k][j][i], Ey_yee[k][j+1][i]);
	        Ez[k-k_l][j-j_l][i-i_l] = interp1(Ez_yee[k-2][j][i], Ez_yee[k-1][j][i], Ez_yee[k][j][i], Ez_yee[k+1][j][i]);
	      }
}

/**
 * @brief Interpolate the TE electric field (x,y components) to the origin of the Yee cell
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
void interpolateFieldCentralE_TE( double ***Ex_yee, double ***Ey_yee, double ***Ez_yee,
				 double ***Ex    , double ***Ey    , double ***Ez    ,
				 int       I     , int       J     , int       K     ,
				 int i_l, int i_u, int j_l, int j_u, int k_l, int k_u){
  
  int i,j,k;
 
  //first check that the limits of field extraction of within range
  checkInterpolationPoints(i_l, i_u, j_l, j_u, k_l, k_u, I, J, K);

  for(k=k_l;k<=k_u;k++)
    for(j=j_l;j<=j_u;j++)
      for(i=i_l;i<=i_u;i++){
	      Ex[k-k_l][j-j_l][i-i_l] = interp1(Ex_yee[k][j][i-2], Ex_yee[k][j][i-1], Ex_yee[k][j][i], Ex_yee[k][j][i+1]);
	      Ey[k-k_l][j-j_l][i-i_l] = interp1(Ey_yee[k][j-2][i], Ey_yee[k][j-1][i], Ey_yee[k][j][i], Ey_yee[k][j+1][i]);
	      Ez[k-k_l][j-j_l][i-i_l] = 0.;
      }
}

/**
 * @brief Interpolate the TM electric field (z component) to the origin of the Yee cell
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
void interpolateFieldCentralE_TM( double ***Ex_yee, double ***Ey_yee, double ***Ez_yee,
				 double ***Ex    , double ***Ey    , double ***Ez    ,
				 int       I     , int       J     , int       K     ,
				 int i_l, int i_u, int j_l, int j_u, int k_l, int k_u){
  
  int i,j,k;
 
  //first check that the limits of field extraction of within range
  checkInterpolationPoints(i_l, i_u, j_l, j_u, k_l, k_u, I, J, K);

  for(k=k_l;k<=k_u;k++)
    for(j=j_l;j<=j_u;j++)
      for(i=i_l;i<=i_u;i++){
	      Ex[k-k_l][j-j_l][i-i_l] = 0.;
	      Ey[k-k_l][j-j_l][i-i_l] = 0.;
	      Ez[k-k_l][j-j_l][i-i_l] = Ez_yee[k][j][i];
      }
}

/**
 * @brief Interpolate the electric field to the origin of the Yee cell - MATLAB interface
 *
 * @param[in] Ex_yee Input Yee cell Ex field matrix
 * @param[in] Ey_yee Input Yee cell Ey field matrix
 * @param[in] Ez_yee Input Yee cell Ez field matrix
 * @param[out] Ex Output Yee cell Ex field matrix
 * @param[out] Ey Output Yee cell Ey field matrix
 * @param[out] Ez Output Yee cell Ez field matrix
 * @param[in] i_l Least i index into the FDTD grid to evaluate the field at. Should be >= 2
 * @param[in] i_u Greatest i index into the FDTD grid to evaluate the field at. Should be <= I-2
 * @param[in] j_l Least j index into the FDTD grid to evaluate the field at. Should be >= 2
 * @param[in] j_u Greatest i index into the FDTD grid to evaluate the field at. Should be <= J-2
 * @param[in] k_l Least k index into the FDTD grid to evaluate the field at. Should be >= 2
 * @param[in] k_u Greatest i index into the FDTD grid to evaluate the field at. Should be <= K-2
 * 
 * I, J, K are the number of elements in the i, j, k directions of the FDTD grid, respectively.
 */
void mxInterpolateFieldCentralE( mxArray *Ex_yee , mxArray *Ey_yee , mxArray *Ez_yee,
				mxArray **Ex    , mxArray **Ey    , mxArray **Ez    ,
				int i_l, int i_u, int j_l, int j_u, int k_l, int k_u){

  double ***ExR, ***ExI, ***EyR, ***EyI, ***EzR, ***EzI,***Ex_yee_R, ***Ex_yee_I, ***Ey_yee_R, ***Ey_yee_I, ***Ez_yee_R, ***Ez_yee_I;
  const int *indims;
  int outdims[3];

  //fprintf(stderr, "mxInterpolateFieldCentralE Pos 00\n");
  if( (int)mxGetNumberOfDimensions( (const mxArray *)Ex_yee) < 3){
    throw runtime_error("Error in mxInterpolateFieldCentralE, Ex_yee does not have 3 dimensions\n");   
  }
  //fprintf(stderr, "mxInterpolateFieldCentralE Pos 01\n");

  indims = (int *)mxGetDimensions( (mxArray *)Ex_yee);
  //fprintf(stderr, "mxInterpolateFieldCentralE(indims): (%d,%d,%d)\n",indims[0],indims[1],indims[2]);
  //fprintf(stderr, "mxInterpolateFieldCentralE: j_u: %d, j_l: %d\n",j_u,j_l);
  //assume that all matrices have the same dimensions - DANGEROUS, SHOULDN'T WE CHECK THIS??
  if( !mxIsComplex(Ex_yee) ){
    mexErrMsgTxt("Ex_yee is not complex");
  }
  if( !mxIsComplex(Ey_yee) ){
    mexErrMsgTxt("Ey_yee is not complex");
  }
  if( !mxIsComplex(Ez_yee) ){
    mexErrMsgTxt("Ez_yee is not complex");
  }

  // cast data in the Yee cells to complex MATLAB arrays
  //fprintf(stderr, "mxInterpolateFieldCentralE Pos 02\n");
  Ex_yee_R = castMatlab3DArray(mxGetPr(Ex_yee), indims[0], indims[1], indims[2]);
  Ex_yee_I = castMatlab3DArray(mxGetPi(Ex_yee), indims[0], indims[1], indims[2]);

  Ey_yee_R = castMatlab3DArray(mxGetPr(Ey_yee), indims[0], indims[1], indims[2]);
  Ey_yee_I = castMatlab3DArray(mxGetPi(Ey_yee), indims[0], indims[1], indims[2]);

  Ez_yee_R = castMatlab3DArray(mxGetPr(Ez_yee), indims[0], indims[1], indims[2]);
  Ez_yee_I = castMatlab3DArray(mxGetPi(Ez_yee), indims[0], indims[1], indims[2]);
  //fprintf(stderr, "mxInterpolateFieldCentralE Pos 03\n");

  //now construct the output matrices
  int ndims = 3;
  
  outdims[0] = i_u - i_l + 1;
  outdims[1] = j_u - j_l + 1;
  if(outdims[1]<1)
    outdims[1]=1;
  outdims[2] = k_u - k_l + 1;
  //fprintf(stderr, "mxInterpolateFieldCentralE(outdims): (%d,%d,%d)\n",outdims[0],outdims[1],outdims[2]);
    //fprintf(stderr, "mxInterpolateFieldCentralE Pos 04\n");
  *Ex = mxCreateNumericArray( ndims, (const mwSize *)outdims, mxDOUBLE_CLASS, mxCOMPLEX);
  *Ey = mxCreateNumericArray( ndims, (const mwSize *)outdims, mxDOUBLE_CLASS, mxCOMPLEX);
  *Ez = mxCreateNumericArray( ndims, (const mwSize *)outdims, mxDOUBLE_CLASS, mxCOMPLEX);
//fprintf(stderr, "mxInterpolateFieldCentralE Pos 04a\n");
  ExR = castMatlab3DArray(mxGetPr(*Ex), outdims[0], outdims[1], outdims[2]);
  ExI = castMatlab3DArray(mxGetPi(*Ex), outdims[0], outdims[1], outdims[2]);
//fprintf(stderr, "mxInterpolateFieldCentralE Pos 04b\n");
  EyR = castMatlab3DArray(mxGetPr(*Ey), outdims[0], outdims[1], outdims[2]);
  EyI = castMatlab3DArray(mxGetPi(*Ey), outdims[0], outdims[1], outdims[2]);
//fprintf(stderr, "mxInterpolateFieldCentralE Pos 04c\n");
  EzR = castMatlab3DArray(mxGetPr(*Ez), outdims[0], outdims[1], outdims[2]);
  EzI = castMatlab3DArray(mxGetPi(*Ez), outdims[0], outdims[1], outdims[2]);
  //fprintf(stderr, "mxInterpolateFieldCentralE Pos 05\n");

  //now interpolate fields
  // imaginary part
  interpolateFieldCentralE( Ex_yee_I, Ey_yee_I, Ez_yee_I, ExI, EyI, EzI, indims[0], indims[1], indims[2], i_l, i_u, j_l, j_u, k_l, k_u);
  // real part
  interpolateFieldCentralE( Ex_yee_R, Ey_yee_R, Ez_yee_R, ExR, EyR, EzR, indims[0], indims[1], indims[2], i_l, i_u, j_l, j_u, k_l, k_u);
  
  //free the extra memory used by casting array
  freeCastMatlab3DArray(Ex_yee_R,indims[2]);
  freeCastMatlab3DArray(Ex_yee_I,indims[2]);
  freeCastMatlab3DArray(Ey_yee_R,indims[2]);
  freeCastMatlab3DArray(Ey_yee_I,indims[2]);
  freeCastMatlab3DArray(Ez_yee_R,indims[2]);
  freeCastMatlab3DArray(Ez_yee_I,indims[2]);

  freeCastMatlab3DArray(ExR,outdims[2]);
  freeCastMatlab3DArray(ExI,outdims[2]);
  freeCastMatlab3DArray(EyR,outdims[2]);
  freeCastMatlab3DArray(EyI,outdims[2]);
  freeCastMatlab3DArray(EzR,outdims[2]);
  freeCastMatlab3DArray(EzI,outdims[2]);
}

/**
 * @brief Interpolate the TE electric field (x,y components) to the origin of the Yee cell - MATLAB interface
 *
 * @param[in] Ex_yee Input Yee cell Ex field matrix
 * @param[in] Ey_yee Input Yee cell Ey field matrix
 * @param[in] Ez_yee Input Yee cell Ez field matrix
 * @param[out] Ex Output Yee cell Ex field matrix
 * @param[out] Ey Output Yee cell Ey field matrix
 * @param[out] Ez Output Yee cell Ez field matrix
 * @param[in] i_l Least i index into the FDTD grid to evaluate the field at. Should be >= 2
 * @param[in] i_u Greatest i index into the FDTD grid to evaluate the field at. Should be <= I-2
 * @param[in] j_l Least j index into the FDTD grid to evaluate the field at. Should be >= 2
 * @param[in] j_u Greatest i index into the FDTD grid to evaluate the field at. Should be <= J-2
 * @param[in] k_l Least k index into the FDTD grid to evaluate the field at. Should be >= 2
 * @param[in] k_u Greatest i index into the FDTD grid to evaluate the field at. Should be <= K-2
 *
 * I, J, K are the number of elements in the i, j, k directions of the FDTD grid, respectively.
 */
void mxInterpolateFieldCentralE_TE( mxArray *Ex_yee , mxArray *Ey_yee , mxArray *Ez_yee,
				   mxArray **Ex    , mxArray **Ey    , mxArray **Ez    ,
				   int i_l, int i_u, int j_l, int j_u, int k_l, int k_u){

  double ***ExR, ***ExI, ***EyR, ***EyI, ***EzR, ***EzI,***Ex_yee_R, ***Ex_yee_I, ***Ey_yee_R, ***Ey_yee_I, ***Ez_yee_R, ***Ez_yee_I;
  const int *indims;
  int outdims[3], ndims;
  
  if( mxGetNumberOfDimensions(Ex_yee) !=2 ){
    throw runtime_error("Error in mxInterpolateFieldCentralE_TE, Ex_yee does not have 2 dimensions\n");
    
    
  }
  
  indims = (int *)mxGetDimensions( (mxArray *)Ex_yee);
  //assume that all matrices have the same dimensions
  if( !mxIsComplex(Ex_yee) ){
    mexErrMsgTxt("Ex_yee is not complex");
  }
  if( !mxIsComplex(Ey_yee) ){
    mexErrMsgTxt("Ey_yee is not complex");
  }
  if( !mxIsComplex(Ez_yee) ){
    mexErrMsgTxt("Ez_yee is not complex");
  }

  Ex_yee_R = castMatlab3DArray(mxGetPr(Ex_yee), indims[0], indims[1],0);
  Ex_yee_I = castMatlab3DArray(mxGetPi(Ex_yee), indims[0], indims[1],0);

  Ey_yee_R = castMatlab3DArray(mxGetPr(Ey_yee), indims[0], indims[1],0);
  Ey_yee_I = castMatlab3DArray(mxGetPi(Ey_yee), indims[0], indims[1],0);

  Ez_yee_R = castMatlab3DArray(mxGetPr(Ez_yee), indims[0], indims[1],0);
  Ez_yee_I = castMatlab3DArray(mxGetPi(Ez_yee), indims[0], indims[1],0);

  //now construct the output matrices
  ndims = 3;
  
  outdims[0] = i_u - i_l + 1;
  outdims[1] = j_u - j_l + 1;
  outdims[2] = 1;
   
  *Ex = mxCreateNumericArray( ndims, (const mwSize *)outdims, mxDOUBLE_CLASS, mxCOMPLEX);
  *Ey = mxCreateNumericArray( ndims, (const mwSize *)outdims, mxDOUBLE_CLASS, mxCOMPLEX);
  *Ez = mxCreateNumericArray( ndims, (const mwSize *)outdims, mxDOUBLE_CLASS, mxCOMPLEX);

  ExR = castMatlab3DArray(mxGetPr(*Ex), outdims[0], outdims[1], outdims[2]);
  ExI = castMatlab3DArray(mxGetPi(*Ex), outdims[0], outdims[1], outdims[2]);

  EyR = castMatlab3DArray(mxGetPr(*Ey), outdims[0], outdims[1], outdims[2]);
  EyI = castMatlab3DArray(mxGetPi(*Ey), outdims[0], outdims[1], outdims[2]);

  EzR = castMatlab3DArray(mxGetPr(*Ez), outdims[0], outdims[1], outdims[2]);
  EzI = castMatlab3DArray(mxGetPi(*Ez), outdims[0], outdims[1], outdims[2]);
  //now finally ready for interpolation
  interpolateFieldCentralE_TE( Ex_yee_I, Ey_yee_I, Ez_yee_I,
			       ExI     , EyI     , EzI    ,
			       indims[0]     , indims[1]     , indims[2]     ,
			       i_l, i_u, j_l, j_u, k_l, k_u);

  interpolateFieldCentralE_TE( Ex_yee_R, Ey_yee_R, Ez_yee_R,
			       ExR     , EyR     , EzR    ,
			       indims[0]     , indims[1]     , indims[2]     ,
			       i_l, i_u, j_l, j_u, k_l, k_u);
  //free the extra memory used by casting array
  freeCastMatlab3DArray(Ex_yee_R,0);
  freeCastMatlab3DArray(Ex_yee_I,0);
  freeCastMatlab3DArray(Ey_yee_R,0);
  freeCastMatlab3DArray(Ey_yee_I,0);
  freeCastMatlab3DArray(Ez_yee_R,0);
  freeCastMatlab3DArray(Ez_yee_I,0);

  freeCastMatlab3DArray(ExR,outdims[2]);
  freeCastMatlab3DArray(ExI,outdims[2]);
  freeCastMatlab3DArray(EyR,outdims[2]);
  freeCastMatlab3DArray(EyI,outdims[2]);
  freeCastMatlab3DArray(EzR,outdims[2]);
  freeCastMatlab3DArray(EzI,outdims[2]);
}

/**
 * @brief Interpolate the TM electric field (z component) to the origin of the Yee cell - MATLAB interface
 *
 * @param[in] Ex_yee Input Yee cell Ex field matrix
 * @param[in] Ey_yee Input Yee cell Ey field matrix
 * @param[in] Ez_yee Input Yee cell Ez field matrix
 * @param[out] Ex Output Yee cell Ex field matrix
 * @param[out] Ey Output Yee cell Ey field matrix
 * @param[out] Ez Output Yee cell Ez field matrix
 * @param[in] i_l Least i index into the FDTD grid to evaluate the field at. Should be >= 2
 * @param[in] i_u Greatest i index into the FDTD grid to evaluate the field at. Should be <= I-2
 * @param[in] j_l Least j index into the FDTD grid to evaluate the field at. Should be >= 2
 * @param[in] j_u Greatest i index into the FDTD grid to evaluate the field at. Should be <= J-2
 * @param[in] k_l Least k index into the FDTD grid to evaluate the field at. Should be >= 2
 * @param[in] k_u Greatest i index into the FDTD grid to evaluate the field at. Should be <= K-2
 *
 * I, J, K are the number of elements in the i, j, k directions of the FDTD grid, respectively.
 */
void mxInterpolateFieldCentralE_TM( mxArray *Ex_yee , mxArray *Ey_yee , mxArray *Ez_yee,
				   mxArray **Ex    , mxArray **Ey    , mxArray **Ez    ,
				   int i_l, int i_u, int j_l, int j_u, int k_l, int k_u){

  double ***ExR, ***ExI, ***EyR, ***EyI, ***EzR, ***EzI,***Ex_yee_R, ***Ex_yee_I, ***Ey_yee_R, ***Ey_yee_I, ***Ez_yee_R, ***Ez_yee_I;
  const int *indims;
  int outdims[3], ndims;
  
  if( mxGetNumberOfDimensions(Ex_yee) != 2){
    throw runtime_error("Error in mxInterpolateFieldCentralE_TM, Ex_yee does not have 2 dimensions\n");
    
    
  }
  
  indims = (int *)mxGetDimensions(Ex_yee);
  //assume that all matrices have the same dimensions
  if( !mxIsComplex(Ex_yee) ){
    mexErrMsgTxt("Ex_yee is not complex");
  }
  if( !mxIsComplex(Ey_yee) ){
    mexErrMsgTxt("Ey_yee is not complex");
  }
  if( !mxIsComplex(Ez_yee) ){
    mexErrMsgTxt("Ez_yee is not complex");
  }


  Ex_yee_R = castMatlab3DArray(mxGetPr(Ex_yee), indims[0], indims[1],0);
  Ex_yee_I = castMatlab3DArray(mxGetPi(Ex_yee), indims[0], indims[1],0);

  Ey_yee_R = castMatlab3DArray(mxGetPr(Ey_yee), indims[0], indims[1],0);
  Ey_yee_I = castMatlab3DArray(mxGetPi(Ey_yee), indims[0], indims[1],0);

  Ez_yee_R = castMatlab3DArray(mxGetPr(Ez_yee), indims[0], indims[1],0);
  Ez_yee_I = castMatlab3DArray(mxGetPi(Ez_yee), indims[0], indims[1],0);

  //now construct the output matrices
  ndims = 3;
  
  outdims[0] = i_u - i_l + 1;
  outdims[1] = j_u - j_l + 1;
  outdims[2] = 1;
   
  *Ex = mxCreateNumericArray( ndims, (const mwSize *)outdims, mxDOUBLE_CLASS, mxCOMPLEX);
  *Ey = mxCreateNumericArray( ndims, (const mwSize *)outdims, mxDOUBLE_CLASS, mxCOMPLEX);
  *Ez = mxCreateNumericArray( ndims, (const mwSize *)outdims, mxDOUBLE_CLASS, mxCOMPLEX);

  ExR = castMatlab3DArray(mxGetPr(*Ex), outdims[0], outdims[1], outdims[2]);
  ExI = castMatlab3DArray(mxGetPi(*Ex), outdims[0], outdims[1], outdims[2]);

  EyR = castMatlab3DArray(mxGetPr(*Ey), outdims[0], outdims[1], outdims[2]);
  EyI = castMatlab3DArray(mxGetPi(*Ey), outdims[0], outdims[1], outdims[2]);

  EzR = castMatlab3DArray(mxGetPr(*Ez), outdims[0], outdims[1], outdims[2]);
  EzI = castMatlab3DArray(mxGetPi(*Ez), outdims[0], outdims[1], outdims[2]);

  //now finally ready for interpolation
  interpolateFieldCentralE_TM( Ex_yee_I, Ey_yee_I, Ez_yee_I,
			       ExI     , EyI     , EzI    ,
			       indims[0]     , indims[1]     , indims[2]     ,
			       i_l, i_u, j_l, j_u, k_l, k_u);


  interpolateFieldCentralE_TM( Ex_yee_R, Ey_yee_R, Ez_yee_R,
			       ExR     , EyR     , EzR    ,
			       indims[0]     , indims[1]     , indims[2]     ,
			       i_l, i_u, j_l, j_u, k_l, k_u);
  
  //free the extra memory used by casting array
  freeCastMatlab3DArray(Ex_yee_R,0);
  freeCastMatlab3DArray(Ex_yee_I,0);
  freeCastMatlab3DArray(Ey_yee_R,0);
  freeCastMatlab3DArray(Ey_yee_I,0);
  freeCastMatlab3DArray(Ez_yee_R,0);
  freeCastMatlab3DArray(Ez_yee_I,0);

  freeCastMatlab3DArray(ExR,outdims[2]);
  freeCastMatlab3DArray(ExI,outdims[2]);
  freeCastMatlab3DArray(EyR,outdims[2]);
  freeCastMatlab3DArray(EyI,outdims[2]);
  freeCastMatlab3DArray(EzR,outdims[2]);
  freeCastMatlab3DArray(EzI,outdims[2]);
}

/**
 * @brief Interpolate the magnetic field to the origin of the Yee cell
 *
 * @param[in] Hx_yee Steady state x component of magnetic field calculated at points in the Yee cell
 * @param[in] Hy_yee Steady state y component of magnetic field calculated at points in the Yee cell
 * @param[in] Hz_yee Steady state z component of magnetic field calculated at points in the Yee cell
 * @param[out] Hx Steady state x component of magnetic field interpolated to Yee cell origin
 * @param[out] Hy Steady state y component of magnetic field interpolated to Yee cell origin
 * @param[out] Hz Steady state z component of magnetic field interpolated to Yee cell origin
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
void interpolateFieldCentralH( double ***Hx_yee, double ***Hy_yee, double ***Hz_yee,
			      double ***Hx    , double ***Hy    , double ***Hz    ,
                              int       I     , int       J     , int       K     ,
			      int i_l, int i_u, int j_l, int j_u, int k_l, int k_u){
  
  int i,j,k;
  double res1, res2, res3, res4;
 
  //  throw runtime_error("Entering interpolateFieldCentralH\n");
  //first check that the limits of field extraction of within range
  checkInterpolationPoints(i_l, i_u, j_l, j_u, k_l, k_u, I, J, K);

  if(j_u<j_l){
    //fprintf(stderr, "interpolateFieldCentralH 02\n");
    j=0;
    for(k=k_l;k<=k_u;k++)
      for(i=i_l;i<=i_u;i++){
	      res1 = Hx_yee[k+1][j][i];
	      res2 = Hx_yee[k  ][j][i];
	      res3 = Hx_yee[k-1][j][i];
	      res4 = Hx_yee[k-2][j][i];
	      Hx[k-k_l][j][i-i_l] = interp1(res1, res2, res3, res4);
	
	      res1 = interp1(Hy_yee[k+1][j][i+1],Hy_yee[k+1][j][i],Hy_yee[k+1][j][i-1],Hy_yee[k+1][j][i-2]);
	      res2 = interp1(Hy_yee[k  ][j][i+1],Hy_yee[k  ][j][i],Hy_yee[k  ][j][i-1],Hy_yee[k  ][j][i-2]);
	      res3 = interp1(Hy_yee[k-1][j][i+1],Hy_yee[k-1][j][i],Hy_yee[k-1][j][i-1],Hy_yee[k-1][j][i-2]);
	      res4 = interp1(Hy_yee[k-2][j][i+1],Hy_yee[k-2][j][i],Hy_yee[k-2][j][i-1],Hy_yee[k-2][j][i-2]);
	      Hy[k-k_l][j][i-i_l] = interp1(res1, res2, res3, res4);

	      Hz[k-k_l][j][i-i_l] = interp1(Hz_yee[k][j  ][i+1],Hz_yee[k][j  ][i],Hz_yee[k][j  ][i-1],Hz_yee[k][j  ][i-2]);
      }

  }
  else
    for(k=k_l;k<=k_u;k++)
      for(j=j_l;j<=j_u;j++)
	      for(i=i_l;i<=i_u;i++){
        /*
          res1 = interp1(Hx_yee[k+1][j][i+1],Hx_yee[k+1][j][i],Hx_yee[k+1][j][i-1],Hx_yee[k+1][j][i-2]);
          res2 = interp1(Hx_yee[k  ][j][i+1],Hx_yee[k  ][j][i],Hx_yee[k  ][j][i-1],Hx_yee[k  ][j][i-2]);
          res3 = interp1(Hx_yee[k-1][j][i+1],Hx_yee[k-1][j][i],Hx_yee[k-1][j][i-1],Hx_yee[k-1][j][i-2]);
          res4 = interp1(Hx_yee[k-2][j][i+1],Hx_yee[k-2][j][i],Hx_yee[k-2][j][i-1],Hx_yee[k-2][j][i-2]);
          Hx[k-k_l][j-j_l][i-i_l] = interp1(res1, res2, res3, res4);
          
          res1 = interp1(Hy_yee[k+1][j+1][i],Hy_yee[k+1][j][i],Hy_yee[k+1][j-1][i],Hy_yee[k+1][j-2][i]);
          res2 = interp1(Hy_yee[k  ][j+1][i],Hy_yee[k  ][j][i],Hy_yee[k  ][j-1][i],Hy_yee[k  ][j-2][i]);
          res3 = interp1(Hy_yee[k-1][j+1][i],Hy_yee[k-1][j][i],Hy_yee[k-1][j-1][i],Hy_yee[k-1][j-2][i]);
          res4 = interp1(Hy_yee[k-2][j+1][i],Hy_yee[k-2][j][i],Hy_yee[k-2][j-1][i],Hy_yee[k-2][j-2][i]);
          Hy[k-k_l][j-j_l][i-i_l] = interp1(res1, res2, res3, res4);
          
          res1 = interp1(Hz_yee[k][j+1][i+1],Hz_yee[k][j+1][i],Hz_yee[k][j+1][i-1],Hz_yee[k][j+1][i-2]);
          res2 = interp1(Hz_yee[k][j  ][i+1],Hz_yee[k][j  ][i],Hz_yee[k][j  ][i-1],Hz_yee[k][j  ][i-2]);
          res3 = interp1(Hz_yee[k][j-1][i+1],Hz_yee[k][j-1][i],Hz_yee[k][j-1][i-1],Hz_yee[k][j-1][i-2]);
          res4 = interp1(Hz_yee[k][j-2][i+1],Hz_yee[k][j-2][i],Hz_yee[k][j-2][i-1],Hz_yee[k][j-2][i-2]);
          Hz[k-k_l][j-j_l][i-i_l] = interp1(res1, res2, res3, res4);
        */
        res1 = interp1(Hx_yee[k+1][j+1][i],Hx_yee[k+1][j][i],Hx_yee[k+1][j-1][i],Hx_yee[k+1][j-2][i]);
        res2 = interp1(Hx_yee[k  ][j+1][i],Hx_yee[k  ][j][i],Hx_yee[k  ][j-1][i],Hx_yee[k  ][j-2][i]);
        res3 = interp1(Hx_yee[k-1][j+1][i],Hx_yee[k-1][j][i],Hx_yee[k-1][j-1][i],Hx_yee[k-1][j-2][i]);
        res4 = interp1(Hx_yee[k-2][j+1][i],Hx_yee[k-2][j][i],Hx_yee[k-2][j-1][i],Hx_yee[k-2][j-2][i]);
        Hx[k-k_l][j-j_l][i-i_l] = interp1(res1, res2, res3, res4);
        
        res1 = interp1(Hy_yee[k+1][j][i+1],Hy_yee[k+1][j][i],Hy_yee[k+1][j][i-1],Hy_yee[k+1][j][i-2]);
        res2 = interp1(Hy_yee[k  ][j][i+1],Hy_yee[k  ][j][i],Hy_yee[k  ][j][i-1],Hy_yee[k  ][j][i-2]);
        res3 = interp1(Hy_yee[k-1][j][i+1],Hy_yee[k-1][j][i],Hy_yee[k-1][j][i-1],Hy_yee[k-1][j][i-2]);
        res4 = interp1(Hy_yee[k-2][j][i+1],Hy_yee[k-2][j][i],Hy_yee[k-2][j][i-1],Hy_yee[k-2][j][i-2]);
        Hy[k-k_l][j-j_l][i-i_l] = interp1(res1, res2, res3, res4);
        
        res1 = interp1(Hz_yee[k][j+1][i+1],Hz_yee[k][j+1][i],Hz_yee[k][j+1][i-1],Hz_yee[k][j+1][i-2]);
        res2 = interp1(Hz_yee[k][j  ][i+1],Hz_yee[k][j  ][i],Hz_yee[k][j  ][i-1],Hz_yee[k][j  ][i-2]);
        res3 = interp1(Hz_yee[k][j-1][i+1],Hz_yee[k][j-1][i],Hz_yee[k][j-1][i-1],Hz_yee[k][j-1][i-2]);
        res4 = interp1(Hz_yee[k][j-2][i+1],Hz_yee[k][j-2][i],Hz_yee[k][j-2][i-1],Hz_yee[k][j-2][i-2]);
        Hz[k-k_l][j-j_l][i-i_l] = interp1(res1, res2, res3, res4);
	    }
}

/**
 * @brief Interpolate the TE magnetic field (z component) to the origin of the Yee cell
 *
 * @param[in] Hx_yee Steady state x component of magnetic field calculated at points in the Yee cell
 * @param[in] Hy_yee Steady state y component of magnetic field calculated at points in the Yee cell
 * @param[in] Hz_yee Steady state z component of magnetic field calculated at points in the Yee cell
 * @param[out] Hx Steady state x component of magnetic field interpolated to Yee cell origin
 * @param[out] Hy Steady state y component of magnetic field interpolated to Yee cell origin
 * @param[out] Hz Steady state z component of magnetic field interpolated to Yee cell origin
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
void interpolateFieldCentralH_TE( double ***Hx_yee, double ***Hy_yee, double ***Hz_yee,
				 double ***Hx    , double ***Hy    , double ***Hz    ,
				 int       I     , int       J     , int       K     ,
				 int i_l, int i_u, int j_l, int j_u, int k_l, int k_u){
  
  int i,j,k;
  double res1, res2, res3, res4;
 
  //first check that the limits of field extraction of within range
  checkInterpolationPoints(i_l, i_u, j_l, j_u, k_l, k_u, I, J, K);

  for(k=k_l;k<=k_u;k++)
    for(j=j_l;j<=j_u;j++)
      for(i=i_l;i<=i_u;i++){
      /*
        res1 = interp1(Hx_yee[k+1][j][i+1],Hx_yee[k+1][j][i],Hx_yee[k+1][j][i-1],Hx_yee[k+1][j][i-2]);
        res2 = interp1(Hx_yee[k  ][j][i+1],Hx_yee[k  ][j][i],Hx_yee[k  ][j][i-1],Hx_yee[k  ][j][i-2]);
        res3 = interp1(Hx_yee[k-1][j][i+1],Hx_yee[k-1][j][i],Hx_yee[k-1][j][i-1],Hx_yee[k-1][j][i-2]);
        res4 = interp1(Hx_yee[k-2][j][i+1],Hx_yee[k-2][j][i],Hx_yee[k-2][j][i-1],Hx_yee[k-2][j][i-2]);
        Hx[k-k_l][j-j_l][i-i_l] = interp1(res1, res2, res3, res4);

        res1 = interp1(Hy_yee[k+1][j+1][i],Hy_yee[k+1][j][i],Hy_yee[k+1][j-1][i],Hy_yee[k+1][j-2][i]);
        res2 = interp1(Hy_yee[k  ][j+1][i],Hy_yee[k  ][j][i],Hy_yee[k  ][j-1][i],Hy_yee[k  ][j-2][i]);
        res3 = interp1(Hy_yee[k-1][j+1][i],Hy_yee[k-1][j][i],Hy_yee[k-1][j-1][i],Hy_yee[k-1][j-2][i]);
        res4 = interp1(Hy_yee[k-2][j+1][i],Hy_yee[k-2][j][i],Hy_yee[k-2][j-1][i],Hy_yee[k-2][j-2][i]);
        Hy[k-k_l][j-j_l][i-i_l] = interp1(res1, res2, res3, res4);

        res1 = interp1(Hz_yee[k][j+1][i+1],Hz_yee[k][j+1][i],Hz_yee[k][j+1][i-1],Hz_yee[k][j+1][i-2]);
        res2 = interp1(Hz_yee[k][j  ][i+1],Hz_yee[k][j  ][i],Hz_yee[k][j  ][i-1],Hz_yee[k][j  ][i-2]);
        res3 = interp1(Hz_yee[k][j-1][i+1],Hz_yee[k][j-1][i],Hz_yee[k][j-1][i-1],Hz_yee[k][j-1][i-2]);
        res4 = interp1(Hz_yee[k][j-2][i+1],Hz_yee[k][j-2][i],Hz_yee[k][j-2][i-1],Hz_yee[k][j-2][i-2]);
        Hz[k-k_l][j-j_l][i-i_l] = interp1(res1, res2, res3, res4);
      */

      Hx[k-k_l][j-j_l][i-i_l] = 0.;

      Hy[k-k_l][j-j_l][i-i_l] = 0.;

      res1 = interp1(Hz_yee[k][j+1][i+1],Hz_yee[k][j+1][i],Hz_yee[k][j+1][i-1],Hz_yee[k][j+1][i-2]);
      res2 = interp1(Hz_yee[k][j  ][i+1],Hz_yee[k][j  ][i],Hz_yee[k][j  ][i-1],Hz_yee[k][j  ][i-2]);
      res3 = interp1(Hz_yee[k][j-1][i+1],Hz_yee[k][j-1][i],Hz_yee[k][j-1][i-1],Hz_yee[k][j-1][i-2]);
      res4 = interp1(Hz_yee[k][j-2][i+1],Hz_yee[k][j-2][i],Hz_yee[k][j-2][i-1],Hz_yee[k][j-2][i-2]);
      Hz[k-k_l][j-j_l][i-i_l] = interp1(res1, res2, res3, res4);
      }
}

/**
 * @brief Interpolate the TE magnetic field (x,y components) to the origin of the Yee cell
 *
 * @param[in] Hx_yee Steady state x component of magnetic field calculated at points in the Yee cell
 * @param[in] Hy_yee Steady state y component of magnetic field calculated at points in the Yee cell
 * @param[in] Hz_yee Steady state z component of magnetic field calculated at points in the Yee cell
 * @param[out] Hx Steady state x component of magnetic field interpolated to Yee cell origin
 * @param[out] Hy Steady state y component of magnetic field interpolated to Yee cell origin
 * @param[out] Hz Steady state z component of magnetic field interpolated to Yee cell origin
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
void interpolateFieldCentralH_TM( double ***Hx_yee, double ***Hy_yee, double ***Hz_yee,
				 double ***Hx    , double ***Hy    , double ***Hz    ,
				 int       I     , int       J     , int       K     ,
				 int i_l, int i_u, int j_l, int j_u, int k_l, int k_u){
  
  int i,j,k;

  //first check that the limits of field extraction of within range
  checkInterpolationPoints(i_l, i_u, j_l, j_u, k_l, k_u, I, J, K);

  for(k=k_l;k<=k_u;k++)
    for(j=j_l;j<=j_u;j++)
      for(i=i_l;i<=i_u;i++){
      /*
        res1 = interp1(Hx_yee[k+1][j][i+1],Hx_yee[k+1][j][i],Hx_yee[k+1][j][i-1],Hx_yee[k+1][j][i-2]);
        res2 = interp1(Hx_yee[k  ][j][i+1],Hx_yee[k  ][j][i],Hx_yee[k  ][j][i-1],Hx_yee[k  ][j][i-2]);
        res3 = interp1(Hx_yee[k-1][j][i+1],Hx_yee[k-1][j][i],Hx_yee[k-1][j][i-1],Hx_yee[k-1][j][i-2]);
        res4 = interp1(Hx_yee[k-2][j][i+1],Hx_yee[k-2][j][i],Hx_yee[k-2][j][i-1],Hx_yee[k-2][j][i-2]);
        Hx[k-k_l][j-j_l][i-i_l] = interp1(res1, res2, res3, res4);

        res1 = interp1(Hy_yee[k+1][j+1][i],Hy_yee[k+1][j][i],Hy_yee[k+1][j-1][i],Hy_yee[k+1][j-2][i]);
        res2 = interp1(Hy_yee[k  ][j+1][i],Hy_yee[k  ][j][i],Hy_yee[k  ][j-1][i],Hy_yee[k  ][j-2][i]);
        res3 = interp1(Hy_yee[k-1][j+1][i],Hy_yee[k-1][j][i],Hy_yee[k-1][j-1][i],Hy_yee[k-1][j-2][i]);
        res4 = interp1(Hy_yee[k-2][j+1][i],Hy_yee[k-2][j][i],Hy_yee[k-2][j-1][i],Hy_yee[k-2][j-2][i]);
        Hy[k-k_l][j-j_l][i-i_l] = interp1(res1, res2, res3, res4);

        res1 = interp1(Hz_yee[k][j+1][i+1],Hz_yee[k][j+1][i],Hz_yee[k][j+1][i-1],Hz_yee[k][j+1][i-2]);
        res2 = interp1(Hz_yee[k][j  ][i+1],Hz_yee[k][j  ][i],Hz_yee[k][j  ][i-1],Hz_yee[k][j  ][i-2]);
        res3 = interp1(Hz_yee[k][j-1][i+1],Hz_yee[k][j-1][i],Hz_yee[k][j-1][i-1],Hz_yee[k][j-1][i-2]);
        res4 = interp1(Hz_yee[k][j-2][i+1],Hz_yee[k][j-2][i],Hz_yee[k][j-2][i-1],Hz_yee[k][j-2][i-2]);
        Hz[k-k_l][j-j_l][i-i_l] = interp1(res1, res2, res3, res4);
      */

      Hx[k-k_l][j-j_l][i-i_l] = interp1(Hx_yee[k][j-2][i], Hx_yee[k][j-1][i], Hx_yee[k][j][i], Hx_yee[k][j+1][i]);

      Hy[k-k_l][j-j_l][i-i_l] = interp1(Hy_yee[k][j][i-2], Hy_yee[k][j][i-1], Hy_yee[k][j][i], Hy_yee[k][j][i+1]);

      Hz[k-k_l][j-j_l][i-i_l] = 0.;
      }
}

/**
 * @brief Interpolate the magnetic field to the origin of the Yee cell - MATLAB interface
 *
 * @param[in] Hx_yee Input Yee cell Hx field matrix
 * @param[in] Hy_yee Input Yee cell Hy field matrix
 * @param[in] Hz_yee Input Yee cell Hz field matrix
 * @param[out] Hx Output Yee cell Hx field matrix
 * @param[out] Hy Output Yee cell Hy field matrix
 * @param[out] Hz Output Yee cell Hz field matrix
 * @param[in] i_l Least i index into the FDTD grid to evaluate the field at. Should be >= 2
 * @param[in] i_u Greatest i index into the FDTD grid to evaluate the field at. Should be <= I-2
 * @param[in] j_l Least j index into the FDTD grid to evaluate the field at. Should be >= 2
 * @param[in] j_u Greatest i index into the FDTD grid to evaluate the field at. Should be <= J-2
 * @param[in] k_l Least k index into the FDTD grid to evaluate the field at. Should be >= 2
 * @param[in] k_u Greatest i index into the FDTD grid to evaluate the field at. Should be <= K-2
 *
 * I, J, K are the number of elements in the i, j, k directions of the FDTD grid, respectively.
 */
void mxInterpolateFieldCentralH( mxArray *Hx_yee , mxArray *Hy_yee , mxArray *Hz_yee,
			        mxArray **Hx    , mxArray **Hy    , mxArray **Hz    ,
				int i_l, int i_u, int j_l, int j_u, int k_l, int k_u){

  double ***HxR, ***HxI, ***HyR, ***HyI, ***HzR, ***HzI,***Hx_yee_R, ***Hx_yee_I, ***Hy_yee_R, ***Hy_yee_I, ***Hz_yee_R, ***Hz_yee_I;
  const int *indims;
  int outdims[3], ndims;
  
  if( mxGetNumberOfDimensions(Hx_yee) < 3){
    throw runtime_error("Error in mxInterpolateFieldCentralH, Ex_yee does not have 3 dimensions\n");
  }
  
  indims = (int *)mxGetDimensions(Hx_yee);
  //assume that all matrices have the same dimensions
  if( !mxIsComplex(Hx_yee) ){
    mexErrMsgTxt("Ex_yee is not complex");
  }
  if( !mxIsComplex(Hy_yee) ){
    mexErrMsgTxt("Ey_yee is not complex");
  }
  if( !mxIsComplex(Hz_yee) ){
    mexErrMsgTxt("Ez_yee is not complex");
  }


  Hx_yee_R = castMatlab3DArray(mxGetPr(Hx_yee), indims[0], indims[1], indims[2]);
  Hx_yee_I = castMatlab3DArray(mxGetPi(Hx_yee), indims[0], indims[1], indims[2]);

  Hy_yee_R = castMatlab3DArray(mxGetPr(Hy_yee), indims[0], indims[1], indims[2]);
  Hy_yee_I = castMatlab3DArray(mxGetPi(Hy_yee), indims[0], indims[1], indims[2]);

  Hz_yee_R = castMatlab3DArray(mxGetPr(Hz_yee), indims[0], indims[1], indims[2]);
  Hz_yee_I = castMatlab3DArray(mxGetPi(Hz_yee), indims[0], indims[1], indims[2]);

  //now construct the output matrices
  ndims = 3;
  
  outdims[0] = i_u - i_l + 1;
  outdims[1] = j_u - j_l + 1;
  if(outdims[1]<1)
    outdims[1]=1;
  outdims[2] = k_u - k_l + 1;
   
  *Hx = mxCreateNumericArray( ndims, (const mwSize *)outdims, mxDOUBLE_CLASS, mxCOMPLEX);
  *Hy = mxCreateNumericArray( ndims, (const mwSize *)outdims, mxDOUBLE_CLASS, mxCOMPLEX);
  *Hz = mxCreateNumericArray( ndims, (const mwSize *)outdims, mxDOUBLE_CLASS, mxCOMPLEX);

  HxR = castMatlab3DArray(mxGetPr(*Hx), outdims[0], outdims[1], outdims[2]);
  HxI = castMatlab3DArray(mxGetPi(*Hx), outdims[0], outdims[1], outdims[2]);

  HyR = castMatlab3DArray(mxGetPr(*Hy), outdims[0], outdims[1], outdims[2]);
  HyI = castMatlab3DArray(mxGetPi(*Hy), outdims[0], outdims[1], outdims[2]);

  HzR = castMatlab3DArray(mxGetPr(*Hz), outdims[0], outdims[1], outdims[2]);
  HzI = castMatlab3DArray(mxGetPi(*Hz), outdims[0], outdims[1], outdims[2]);

  //now finally ready for interpolation
  interpolateFieldCentralH( Hx_yee_I, Hy_yee_I, Hz_yee_I,
			    HxI     , HyI     , HzI    ,
                            indims[0]     , indims[1]     , indims[2]     ,
			    i_l, i_u, j_l, j_u, k_l, k_u);


  interpolateFieldCentralH( Hx_yee_R, Hy_yee_R, Hz_yee_R,
			    HxR     , HyR     , HzR    ,
			    indims[0]     , indims[1]     , indims[2]     ,
			    i_l, i_u, j_l, j_u, k_l, k_u);
  
  //free the extra memory used by casting array
  freeCastMatlab3DArray(Hx_yee_R,indims[2]);
  freeCastMatlab3DArray(Hx_yee_I,indims[2]);
  freeCastMatlab3DArray(Hy_yee_R,indims[2]);
  freeCastMatlab3DArray(Hy_yee_I,indims[2]);
  freeCastMatlab3DArray(Hz_yee_R,indims[2]);
  freeCastMatlab3DArray(Hz_yee_I,indims[2]);

  freeCastMatlab3DArray(HxR,outdims[2]);
  freeCastMatlab3DArray(HxI,outdims[2]);
  freeCastMatlab3DArray(HyR,outdims[2]);
  freeCastMatlab3DArray(HyI,outdims[2]);
  freeCastMatlab3DArray(HzR,outdims[2]);
  freeCastMatlab3DArray(HzI,outdims[2]);
}

/**
 * @brief Interpolate the TE magnetic field (z component) to the origin of the Yee cell - MATLAB interface
 *
 * @param[in] Hx_yee Input Yee cell Hx field matrix
 * @param[in] Hy_yee Input Yee cell Hy field matrix
 * @param[in] Hz_yee Input Yee cell Hz field matrix
 * @param[out] Hx Output Yee cell Hx field matrix
 * @param[out] Hy Output Yee cell Hy field matrix
 * @param[out] Hz Output Yee cell Hz field matrix
 * @param[in] i_l Least i index into the FDTD grid to evaluate the field at. Should be >= 2
 * @param[in] i_u Greatest i index into the FDTD grid to evaluate the field at. Should be <= I-2
 * @param[in] j_l Least j index into the FDTD grid to evaluate the field at. Should be >= 2
 * @param[in] j_u Greatest i index into the FDTD grid to evaluate the field at. Should be <= J-2
 * @param[in] k_l Least k index into the FDTD grid to evaluate the field at. Should be >= 2
 * @param[in] k_u Greatest i index into the FDTD grid to evaluate the field at. Should be <= K-2
 *
 * I, J, K are the number of elements in the i, j, k directions of the FDTD grid, respectively.
 */
void mxInterpolateFieldCentralH_TE( mxArray *Hx_yee , mxArray *Hy_yee , mxArray *Hz_yee,
				   mxArray **Hx    , mxArray **Hy    , mxArray **Hz    ,
				   int i_l, int i_u, int j_l, int j_u, int k_l, int k_u){

  double ***HxR, ***HxI, ***HyR, ***HyI, ***HzR, ***HzI,***Hx_yee_R, ***Hx_yee_I, ***Hy_yee_R, ***Hy_yee_I, ***Hz_yee_R, ***Hz_yee_I;
  const int *indims;
  int outdims[3], ndims;
  
  if( mxGetNumberOfDimensions(Hx_yee) != 2){
    throw runtime_error("Error in mxInterpolateFieldCentralH_TE, Ex_yee does not have 3 dimensions\n");
    
    
  }
  
  indims = (int *)mxGetDimensions(Hx_yee);
  //assume that all matrices have the same dimensions
  if( !mxIsComplex(Hx_yee) ){
    mexErrMsgTxt("Ex_yee is not complex");
  }
  if( !mxIsComplex(Hy_yee) ){
    mexErrMsgTxt("Ey_yee is not complex");
  }
  if( !mxIsComplex(Hz_yee) ){
    mexErrMsgTxt("Ez_yee is not complex");
  }


  Hx_yee_R = castMatlab3DArray(mxGetPr(Hx_yee), indims[0], indims[1], 0);
  Hx_yee_I = castMatlab3DArray(mxGetPi(Hx_yee), indims[0], indims[1], 0);

  Hy_yee_R = castMatlab3DArray(mxGetPr(Hy_yee), indims[0], indims[1], 0);
  Hy_yee_I = castMatlab3DArray(mxGetPi(Hy_yee), indims[0], indims[1], 0);

  Hz_yee_R = castMatlab3DArray(mxGetPr(Hz_yee), indims[0], indims[1], 0);
  Hz_yee_I = castMatlab3DArray(mxGetPi(Hz_yee), indims[0], indims[1], 0);

  //now construct the output matrices
  ndims = 3;
  
  outdims[0] = i_u - i_l + 1;
  outdims[1] = j_u - j_l + 1;
  outdims[2] = 1;
   
  *Hx = mxCreateNumericArray( ndims, (const mwSize *)outdims, mxDOUBLE_CLASS, mxCOMPLEX);
  *Hy = mxCreateNumericArray( ndims, (const mwSize *)outdims, mxDOUBLE_CLASS, mxCOMPLEX);
  *Hz = mxCreateNumericArray( ndims, (const mwSize *)outdims, mxDOUBLE_CLASS, mxCOMPLEX);

  HxR = castMatlab3DArray(mxGetPr(*Hx), outdims[0], outdims[1], outdims[2]);
  HxI = castMatlab3DArray(mxGetPi(*Hx), outdims[0], outdims[1], outdims[2]);

  HyR = castMatlab3DArray(mxGetPr(*Hy), outdims[0], outdims[1], outdims[2]);
  HyI = castMatlab3DArray(mxGetPi(*Hy), outdims[0], outdims[1], outdims[2]);

  HzR = castMatlab3DArray(mxGetPr(*Hz), outdims[0], outdims[1], outdims[2]);
  HzI = castMatlab3DArray(mxGetPi(*Hz), outdims[0], outdims[1], outdims[2]);

  //now finally ready for interpolation
  interpolateFieldCentralH_TE( Hx_yee_I, Hy_yee_I, Hz_yee_I,
			       HxI     , HyI     , HzI    ,
			       indims[0]     , indims[1]     , indims[2]     ,
			       i_l, i_u, j_l, j_u, k_l, k_u);


  interpolateFieldCentralH_TE( Hx_yee_R, Hy_yee_R, Hz_yee_R,
			       HxR     , HyR     , HzR    ,
			       indims[0]     , indims[1]     , indims[2]     ,
			       i_l, i_u, j_l, j_u, k_l, k_u);
  
  //free the extra memory used by casting array
  freeCastMatlab3DArray(Hx_yee_R,0);
  freeCastMatlab3DArray(Hx_yee_I,0);
  freeCastMatlab3DArray(Hy_yee_R,0);
  freeCastMatlab3DArray(Hy_yee_I,0);
  freeCastMatlab3DArray(Hz_yee_R,0);
  freeCastMatlab3DArray(Hz_yee_I,0);

  freeCastMatlab3DArray(HxR,outdims[2]);
  freeCastMatlab3DArray(HxI,outdims[2]);
  freeCastMatlab3DArray(HyR,outdims[2]);
  freeCastMatlab3DArray(HyI,outdims[2]);
  freeCastMatlab3DArray(HzR,outdims[2]);
  freeCastMatlab3DArray(HzI,outdims[2]);
}

/**
 * @brief Interpolate the TM magnetic field (x,y components) to the origin of the Yee cell - MATLAB interface
 *
 * @param[in] Hx_yee Input Yee cell Hx field matrix
 * @param[in] Hy_yee Input Yee cell Hy field matrix
 * @param[in] Hz_yee Input Yee cell Hz field matrix
 * @param[out] Hx Output Yee cell Hx field matrix
 * @param[out] Hy Output Yee cell Hy field matrix
 * @param[out] Hz Output Yee cell Hz field matrix
 * @param[in] i_l Least i index into the FDTD grid to evaluate the field at. Should be >= 2
 * @param[in] i_u Greatest i index into the FDTD grid to evaluate the field at. Should be <= I-2
 * @param[in] j_l Least j index into the FDTD grid to evaluate the field at. Should be >= 2
 * @param[in] j_u Greatest i index into the FDTD grid to evaluate the field at. Should be <= J-2
 * @param[in] k_l Least k index into the FDTD grid to evaluate the field at. Should be >= 2
 * @param[in] k_u Greatest i index into the FDTD grid to evaluate the field at. Should be <= K-2
 *
 * I, J, K are the number of elements in the i, j, k directions of the FDTD grid, respectively.
 */
void mxInterpolateFieldCentralH_TM( mxArray *Hx_yee , mxArray *Hy_yee , mxArray *Hz_yee,
				   mxArray **Hx    , mxArray **Hy    , mxArray **Hz    ,
				   int i_l, int i_u, int j_l, int j_u, int k_l, int k_u){

  double ***HxR, ***HxI, ***HyR, ***HyI, ***HzR, ***HzI,***Hx_yee_R, ***Hx_yee_I, ***Hy_yee_R, ***Hy_yee_I, ***Hz_yee_R, ***Hz_yee_I;
  const int *indims;
  int outdims[3], ndims;
  
  if( mxGetNumberOfDimensions(Hx_yee) != 2){
    throw runtime_error("Error in mxInterpolateFieldCentralH_TM, Ex_yee does not have 3 dimensions\n");
  }
  
  indims = (int *)mxGetDimensions(Hx_yee);
  //assume that all matrices have the same dimensions
  if( !mxIsComplex(Hx_yee) ){
    mexErrMsgTxt("Ex_yee is not complex");
  }
  if( !mxIsComplex(Hy_yee) ){
    mexErrMsgTxt("Ey_yee is not complex");
  }
  if( !mxIsComplex(Hz_yee) ){
    mexErrMsgTxt("Ez_yee is not complex");
  }


  Hx_yee_R = castMatlab3DArray(mxGetPr(Hx_yee), indims[0], indims[1], 0);
  Hx_yee_I = castMatlab3DArray(mxGetPi(Hx_yee), indims[0], indims[1], 0);

  Hy_yee_R = castMatlab3DArray(mxGetPr(Hy_yee), indims[0], indims[1], 0);
  Hy_yee_I = castMatlab3DArray(mxGetPi(Hy_yee), indims[0], indims[1], 0);

  Hz_yee_R = castMatlab3DArray(mxGetPr(Hz_yee), indims[0], indims[1], 0);
  Hz_yee_I = castMatlab3DArray(mxGetPi(Hz_yee), indims[0], indims[1], 0);

  //now construct the output matrices
  ndims = 3;
  
  outdims[0] = i_u - i_l + 1;
  outdims[1] = j_u - j_l + 1;
  outdims[2] = 1;
   
  *Hx = mxCreateNumericArray( ndims, (const mwSize *)outdims, mxDOUBLE_CLASS, mxCOMPLEX);
  *Hy = mxCreateNumericArray( ndims, (const mwSize *)outdims, mxDOUBLE_CLASS, mxCOMPLEX);
  *Hz = mxCreateNumericArray( ndims, (const mwSize *)outdims, mxDOUBLE_CLASS, mxCOMPLEX);

  HxR = castMatlab3DArray(mxGetPr(*Hx), outdims[0], outdims[1], outdims[2]);
  HxI = castMatlab3DArray(mxGetPi(*Hx), outdims[0], outdims[1], outdims[2]);

  HyR = castMatlab3DArray(mxGetPr(*Hy), outdims[0], outdims[1], outdims[2]);
  HyI = castMatlab3DArray(mxGetPi(*Hy), outdims[0], outdims[1], outdims[2]);

  HzR = castMatlab3DArray(mxGetPr(*Hz), outdims[0], outdims[1], outdims[2]);
  HzI = castMatlab3DArray(mxGetPi(*Hz), outdims[0], outdims[1], outdims[2]);

  //now finally ready for interpolation
  interpolateFieldCentralH_TM( Hx_yee_I, Hy_yee_I, Hz_yee_I,
			       HxI     , HyI     , HzI    ,
			       indims[0]     , indims[1]     , indims[2]     ,
			       i_l, i_u, j_l, j_u, k_l, k_u);


  interpolateFieldCentralH_TM( Hx_yee_R, Hy_yee_R, Hz_yee_R,
			       HxR     , HyR     , HzR    ,
			       indims[0]     , indims[1]     , indims[2]     ,
			       i_l, i_u, j_l, j_u, k_l, k_u);
  
  //free the extra memory used by casting array
  freeCastMatlab3DArray(Hx_yee_R,0);
  freeCastMatlab3DArray(Hx_yee_I,0);
  freeCastMatlab3DArray(Hy_yee_R,0);
  freeCastMatlab3DArray(Hy_yee_I,0);
  freeCastMatlab3DArray(Hz_yee_R,0);
  freeCastMatlab3DArray(Hz_yee_I,0);

  freeCastMatlab3DArray(HxR,outdims[2]);
  freeCastMatlab3DArray(HxI,outdims[2]);
  freeCastMatlab3DArray(HyR,outdims[2]);
  freeCastMatlab3DArray(HyI,outdims[2]);
  freeCastMatlab3DArray(HzR,outdims[2]);
  freeCastMatlab3DArray(HzI,outdims[2]);
}

/*Interpolate the electric field to the origin of the Yee cell in the time domain
 *
 *[in]Hxy to Hzy are the split components of the Yee cell
 *[in](i,j,k) denotes the particular Yee cell we are interested in
 *
 *[out]Hx to Hz are the field components interpolated to the origin of the Yee cell
 */
void interpolateTimeDomainFieldCentralH( double ***Hxy, double ***Hxz, double ***Hyx, double ***Hyz, double ***Hzx, double ***Hzy,
					int i, int j, int k,
					double *Hx, double *Hy, double *Hz){
  

  double res1, res2, res3, res4;

  res1 = interp1(Hxy[k+1][j+1][i]+Hxz[k+1][j+1][i],Hxy[k+1][j][i]+Hxz[k+1][j][i],Hxy[k+1][j-1][i]+Hxz[k+1][j-1][i],Hxy[k+1][j-2][i]+Hxz[k+1][j-2][i]);
  res2 = interp1(Hxy[k  ][j+1][i]+Hxz[k  ][j+1][i],Hxy[k  ][j][i]+Hxz[k  ][j][i],Hxy[k  ][j-1][i]+Hxz[k  ][j-1][i],Hxy[k  ][j-2][i]+Hxz[k  ][j-2][i]);
  res3 = interp1(Hxy[k-1][j+1][i]+Hxz[k-1][j+1][i],Hxy[k-1][j][i]+Hxz[k-1][j][i],Hxy[k-1][j-1][i]+Hxz[k-1][j-1][i],Hxy[k-1][j-2][i]+Hxz[k-1][j-2][i]);
  res4 = interp1(Hxy[k-2][j+1][i]+Hxz[k-2][j+1][i],Hxy[k-2][j][i]+Hxz[k-2][j][i],Hxy[k-2][j-1][i]+Hxz[k-2][j-1][i],Hxy[k-2][j-2][i]+Hxz[k-2][j-2][i]);
  *Hx  = interp1(res1, res2, res3, res4);
  
  res1 = interp1(Hyx[k+1][j][i+1]+Hyz[k+1][j][i+1],Hyx[k+1][j][i]+Hyz[k+1][j][i],Hyx[k+1][j][i-1]+Hyz[k+1][j][i-1],Hyx[k+1][j][i-2]+Hyz[k+1][j][i-2]);
  res2 = interp1(Hyx[k  ][j][i+1]+Hyz[k  ][j][i+1],Hyx[k  ][j][i]+Hyz[k  ][j][i],Hyx[k  ][j][i-1]+Hyz[k  ][j][i-1],Hyx[k  ][j][i-2]+Hyz[k  ][j][i-2]);
  res3 = interp1(Hyx[k-1][j][i+1]+Hyz[k-1][j][i+1],Hyx[k-1][j][i]+Hyz[k-1][j][i],Hyx[k-1][j][i-1]+Hyz[k-1][j][i-1],Hyx[k-1][j][i-2]+Hyz[k-1][j][i-2]);
  res4 = interp1(Hyx[k-2][j][i+1]+Hyz[k-2][j][i+1],Hyx[k-2][j][i]+Hyz[k-2][j][i],Hyx[k-2][j][i-1]+Hyz[k-2][j][i-1],Hyx[k-2][j][i-2]+Hyz[k-2][j][i-2]);
  *Hy  = interp1(res1, res2, res3, res4);
  
  res1 = interp1(Hzx[k][j+1][i+1]+Hzy[k][j+1][i+1],Hzx[k][j+1][i]+Hzy[k][j+1][i],Hzx[k][j+1][i-1]+Hzy[k][j+1][i-1],Hzx[k][j+1][i-2]+Hzy[k][j+1][i-2]);
  res2 = interp1(Hzx[k][j  ][i+1]+Hzy[k][j  ][i+1],Hzx[k][j  ][i]+Hzy[k][j  ][i],Hzx[k][j  ][i-1]+Hzy[k][j  ][i-1],Hzx[k][j  ][i-2]+Hzy[k][j  ][i-2]);
  res3 = interp1(Hzx[k][j-1][i+1]+Hzy[k][j-1][i+1],Hzx[k][j-1][i]+Hzy[k][j-1][i],Hzx[k][j-1][i-1]+Hzy[k][j-1][i-1],Hzx[k][j-1][i-2]+Hzy[k][j-1][i-2]);
  res4 = interp1(Hzx[k][j-2][i+1]+Hzy[k][j-2][i+1],Hzx[k][j-2][i]+Hzy[k][j-2][i],Hzx[k][j-2][i-1]+Hzy[k][j-2][i-1],Hzx[k][j-2][i-2]+Hzy[k][j-2][i-2]);
  *Hz  = interp1(res1, res2, res3, res4);
}

/*Interpolate the electric field to the origin of the Yee cell in the time domain using band limited interpolation
 *
 *[in]Hxy to Hzy are the split components of the Yee cell
 *[in](i,j,k) denotes the particular Yee cell we are interested in
 *
 *[out]Hx to Hz are the field components interpolated to the origin of the Yee cell
 */
void interpolateTimeDomainFieldCentralHBandLimited( double ***Hxy, double ***Hxz, double ***Hyx, double ***Hyz, double ***Hzx, double ***Hzy,
					int i, int j, int k,
					double *Hx, double *Hy, double *Hz){

  /*Array for performing bandwidth limited interpolation obtained using Matlab's interp function*/
  const int Nbvec = 8;
  const double bvec[Nbvec] = {-0.006777513830606,0.039457774230186,-0.142658093428622,0.609836360661632,0.609836360661632,-0.142658093428622,0.039457774230186,-0.006777513830606};
  *Hx = 0.;
  *Hy = 0.;
  *Hz = 0.;

  double hx = 0;
  double hy = 0;
  double hz = 0;
  

  for(int ind1=0;ind1<Nbvec;ind1++){
    hx = 0.;
    hy = 0.;
    hz = 0.;
    for(int ind2=0;ind2<Nbvec;ind2++){
      hx += bvec[ind2]*(Hxy[k-Nbvec/2+ind2][j-Nbvec/2+ind1][i]+Hxz[k-Nbvec/2+ind2][j-Nbvec/2+ind1][i]);
      hy += bvec[ind2]*(Hyx[k-Nbvec/2+ind2][j][i-Nbvec/2+ind1]+Hyz[k-Nbvec/2+ind2][j][i-Nbvec/2+ind1]);
      hz += bvec[ind2]*(Hzx[k][j-Nbvec/2+ind2][i-Nbvec/2+ind1]+Hzy[k][j-Nbvec/2+ind2][i-Nbvec/2+ind1]);
    }
    *Hx += bvec[ind1]*hx;
    *Hy += bvec[ind1]*hy;
    *Hz += bvec[ind1]*hz;
  }
}

/*Interpolate the electric field to the origin of the Yee cell in the time domain
 *
 *[in]Hxy to Hzy are the split components of the Yee cell
 *[in](i,j,k) denotes the particular Yee cell we are interested in
 *
 *[out]Hx to Hz are the field components interpolated to the origin of the Yee cell
 */
void interpolateTimeDomainFieldCentralH_2Dy( double ***Hxy, double ***Hxz, double ***Hyx, double ***Hyz, double ***Hzx, double ***Hzy,
					int i, int j, int k,
					double *Hx, double *Hy, double *Hz){
  

  double res1, res2, res3, res4;
  
  *Hx  = interp1(Hxy[k+1][j][i]+Hxz[k+1][j][i],Hxy[k][j][i]+Hxz[k][j][i],Hxy[k-1][j][i]+Hxz[k-1][j][i],Hxy[k-2][j][i]+Hxz[k-2][j][i]);
  
  res1 = interp1(Hyx[k+1][j][i+1]+Hyz[k+1][j][i+1],Hyx[k+1][j][i]+Hyz[k+1][j][i],Hyx[k+1][j][i-1]+Hyz[k+1][j][i-1],Hyx[k+1][j][i-2]+Hyz[k+1][j][i-2]);
  res2 = interp1(Hyx[k  ][j][i+1]+Hyz[k  ][j][i+1],Hyx[k  ][j][i]+Hyz[k  ][j][i],Hyx[k  ][j][i-1]+Hyz[k  ][j][i-1],Hyx[k  ][j][i-2]+Hyz[k  ][j][i-2]);
  res3 = interp1(Hyx[k-1][j][i+1]+Hyz[k-1][j][i+1],Hyx[k-1][j][i]+Hyz[k-1][j][i],Hyx[k-1][j][i-1]+Hyz[k-1][j][i-1],Hyx[k-1][j][i-2]+Hyz[k-1][j][i-2]);
  res4 = interp1(Hyx[k-2][j][i+1]+Hyz[k-2][j][i+1],Hyx[k-2][j][i]+Hyz[k-2][j][i],Hyx[k-2][j][i-1]+Hyz[k-2][j][i-1],Hyx[k-2][j][i-2]+Hyz[k-2][j][i-2]);
  *Hy  = interp1(res1, res2, res3, res4);
  *Hz  = interp1(Hzx[k][j][i+1]+Hzy[k][j][i+1],Hzx[k][j][i]+Hzy[k][j][i],Hzx[k][j][i-1]+Hzy[k][j][i-1],Hzx[k][j][i-1]+Hzy[k][j][i-1]);
}


/*Interpolate the electric field to the origin of the Yee cell in the time domain
 *
 *[in]Hxy to Hzy are the split components of the Yee cell
 *[in](i,j,k) denotes the particular Yee cell we are interested in
 *
 *[out]Hx to Hz are the field components interpolated to the origin of the Yee cell
 */
void interpolateTimeDomainFieldCentralH_TE( double ***Hxy, double ***Hxz, double ***Hyx, double ***Hyz, double ***Hzx, double ***Hzy,
					   int i, int j, int k,
					   double *Hx, double *Hy, double *Hz){
  

  double res1, res2, res3, res4;

  *Hx  = 0.;
  *Hy  = 0.;
  
  res1 = interp1(Hzx[k][j+1][i+1]+Hzy[k][j+1][i+1],Hzx[k][j+1][i]+Hzy[k][j+1][i],Hzx[k][j+1][i-1]+Hzy[k][j+1][i-1],Hzx[k][j+1][i-2]+Hzy[k][j+1][i-2]);
  res2 = interp1(Hzx[k][j  ][i+1]+Hzy[k][j  ][i+1],Hzx[k][j  ][i]+Hzy[k][j  ][i],Hzx[k][j  ][i-1]+Hzy[k][j  ][i-1],Hzx[k][j  ][i-2]+Hzy[k][j  ][i-2]);
  res3 = interp1(Hzx[k][j-1][i+1]+Hzy[k][j-1][i+1],Hzx[k][j-1][i]+Hzy[k][j-1][i],Hzx[k][j-1][i-1]+Hzy[k][j-1][i-1],Hzx[k][j-1][i-2]+Hzy[k][j-1][i-2]);
  res4 = interp1(Hzx[k][j-2][i+1]+Hzy[k][j-2][i+1],Hzx[k][j-2][i]+Hzy[k][j-2][i],Hzx[k][j-2][i-1]+Hzy[k][j-2][i-1],Hzx[k][j-2][i-2]+Hzy[k][j-2][i-2]);
  *Hz  = interp1(res1, res2, res3, res4);
}

/*Interpolate the electric field to the origin of the Yee cell in the time domain
 *
 *[in]Hxy to Hzy are the split components of the Yee cell
 *[in](i,j,k) denotes the particular Yee cell we are interested in
 *
 *[out]Hx to Hz are the field components interpolated to the origin of the Yee cell
 */
void interpolateTimeDomainFieldCentralH_TM( double ***Hxy, double ***Hxz, double ***Hyx, double ***Hyz, double ***Hzx, double ***Hzy,
					   int i, int j, int k,
					   double *Hx, double *Hy, double *Hz){
  

  *Hx  = interp1(Hxy[k][j-2][i]+Hxz[k][j-2][i], Hxy[k][j-1][i]+Hxz[k][j-1][i], Hxy[k][j][i]+Hxz[k][j][i], Hxy[k][j+1][i]+Hxz[k][j+1][i]); 
  *Hy  = interp1(Hyx[k][j][i-2]+Hyz[k][j][i-2], Hyx[k][j][i-1]+Hyz[k][j][i-1], Hyx[k][j][i]+Hyz[k][j][i], Hyx[k][j][i+1]+Hyz[k][j][i+1]);
  *Hz  = 0.;
}

/*Interpolate the electric field to the origin of the Yee cell in the time domain
 * 
 *[in]Exy to Ezy are the split components of the Yee cell
 *[in](i,j,k) denotes the particular Yee cell we are interested in
 *
 *[out]Ex to Ez are the field components interpolated to the origin of the Yee cell
 *
 */
void interpolateTimeDomainFieldCentralE(  double ***Exy, double ***Exz, double ***Eyx, double ***Eyz, double ***Ezx, double ***Ezy,
					 int i, int j, int k,
					 double *Ex, double *Ey, double *Ez){
 
  *Ex = interp1(Exy[k][j][i-2]+Exz[k][j][i-2], Exy[k][j][i-1]+Exz[k][j][i-1], Exy[k][j][i]+Exz[k][j][i], Exy[k][j][i+1]+Exz[k][j][i+1]);
  *Ey = interp1(Eyx[k][j-2][i]+Eyz[k][j-2][i], Eyx[k][j-1][i]+Eyz[k][j-1][i], Eyx[k][j][i]+Eyz[k][j][i], Eyx[k][j+1][i]+Eyz[k][j+1][i]);
  *Ez = interp1(Ezy[k-2][j][i]+Ezx[k-2][j][i], Ezy[k-1][j][i]+Ezx[k-1][j][i], Ezy[k][j][i]+Ezx[k][j][i], Ezy[k+1][j][i]+Ezx[k+1][j][i]);
}


/*Interpolate the electric field to the origin of the Yee cell in the time domain using band limited interpolation
 * 
 *[in]Exy to Ezy are the split components of the Yee cell
 *[in](i,j,k) denotes the particular Yee cell we are interested in
 *
 *[out]Ex to Ez are the field components interpolated to the origin of the Yee cell
 *
 */
void interpolateTimeDomainFieldCentralEBandLimited(  double ***Exy, double ***Exz, double ***Eyx, double ***Eyz, double ***Ezx, double ***Ezy,
					 int i, int j, int k,
					 double *Ex, double *Ey, double *Ez){
  
  /*Array for performing bandwidth limited interpolation obtained using Matlab's interp function*/
  const int Nbvec = 8;
  const double bvec[Nbvec] = {-0.006777513830606,0.039457774230186,-0.142658093428622,0.609836360661632,0.609836360661632,-0.142658093428622,0.039457774230186,-0.006777513830606};
  
  *Ex = 0.;
  *Ey = 0.;
  *Ez = 0.;
  if(j != 0)
    for(int ind=0;ind<Nbvec;ind++){
      *Ex += (Exy[k][j][i-Nbvec/2+ind]+Exz[k][j][i-Nbvec/2+ind])*bvec[ind];
      *Ey += (Eyx[k][j-Nbvec/2+ind][i]+Eyz[k][j-Nbvec/2+ind][i])*bvec[ind];
      *Ez += (Ezy[k-Nbvec/2+ind][j][i]+Ezx[k-Nbvec/2+ind][j][i])*bvec[ind];
    }
  else{
    for(int ind=0;ind<Nbvec;ind++){
      *Ex += (Exy[k][j][i-Nbvec/2+ind]+Exz[k][j][i-Nbvec/2+ind])*bvec[ind];
      //*Ey += (Eyx[k][j-Nbvec/2+ind][i]+Eyz[k][j-Nbvec/2+ind][i])*bvec[ind];
      *Ez += (Ezy[k-Nbvec/2+ind][j][i]+Ezx[k-Nbvec/2+ind][j][i])*bvec[ind];
    }
    *Ey = (Eyx[k][j][i]+Eyz[k][j][i]);
  }

}

/*Interpolate the electric field to the origin of the Yee cell in the time domain
 * 
 *[in]Exy to Ezy are the split components of the Yee cell
 *[in](i,j,k) denotes the particular Yee cell we are interested in
 *
 *[out]Ex to Ez are the field components interpolated to the origin of the Yee cell
 *
 */
void interpolateTimeDomainFieldCentralE_2Dy(  double ***Exy, double ***Exz, double ***Eyx, double ***Eyz, double ***Ezx, double ***Ezy,
					 int i, int j, int k,
					 double *Ex, double *Ey, double *Ez){
  *Ex = interp1(Exy[k][j][i-2]+Exz[k][j][i-2], Exy[k][j][i-1]+Exz[k][j][i-1], Exy[k][j][i]+Exz[k][j][i], Exy[k][j][i+1]+Exz[k][j][i+1]);
  *Ey = Eyx[k][j][i]+Eyz[k][j][i];
  *Ez = interp1(Ezy[k-2][j][i]+Ezx[k-2][j][i], Ezy[k-1][j][i]+Ezx[k-1][j][i], Ezy[k][j][i]+Ezx[k][j][i], Ezy[k+1][j][i]+Ezx[k+1][j][i]);
}


/*Interpolate the electric field to the origin of the Yee cell in the time domain
 * 
 *[in]Exy to Ezy are the split components of the Yee cell
 *[in](i,j,k) denotes the particular Yee cell we are interested in
 *
 *[out]Ex to Ez are the field components interpolated to the origin of the Yee cell
 *
 */
void interpolateTimeDomainFieldCentralE_TE(  double ***Exy, double ***Exz, double ***Eyx, double ***Eyz, double ***Ezx, double ***Ezy,
					    int i, int j, int k,
					    double *Ex, double *Ey, double *Ez){
 
  *Ex = interp1(Exy[k][j][i-2]+Exz[k][j][i-2], Exy[k][j][i-1]+Exz[k][j][i-1], Exy[k][j][i]+Exz[k][j][i], Exy[k][j][i+1]+Exz[k][j][i+1]);
  *Ey = interp1(Eyx[k][j-2][i]+Eyz[k][j-2][i], Eyx[k][j-1][i]+Eyz[k][j-1][i], Eyx[k][j][i]+Eyz[k][j][i], Eyx[k][j+1][i]+Eyz[k][j+1][i]);
  *Ez = 0.;
}

/*Interpolate the electric field to the origin of the Yee cell in the time domain
 * 
 *[in]Exy to Ezy are the split components of the Yee cell
 *[in](i,j,k) denotes the particular Yee cell we are interested in
 *
 *[out]Ex to Ez are the field components interpolated to the origin of the Yee cell
 *
 */
void interpolateTimeDomainFieldCentralE_TM(  double ***Exy, double ***Exz, double ***Eyx, double ***Eyz, double ***Ezx, double ***Ezy,
					    int i, int j, int k,
					    double *Ex, double *Ey, double *Ez){
 
  *Ex = 0.;
  *Ey = 0.;
  *Ez = Ezy[k][j][i]+Ezx[k][j][i];
}