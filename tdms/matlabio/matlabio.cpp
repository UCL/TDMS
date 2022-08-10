/*****************************************************************
 *
 *  Project.....:  isotropic FDTD code
 *  Application.:  matlab data structure manipulation
 *  Module......:  matlabio.cpp
 *  Description.:  code for processing the matlab data structures
 *  Compiler....:  g++
 *  Written by..:  Peter Munro, Imperial College London, 2002-2008
 *  Environment.:  Linux
 *  Modified....:  Numerous times
 *
 ******************************************************************/

/*---------------------------------------------------------------*/
//                        INCLUDE section
/*---------------------------------------------------------------*/

#include "mat_io.h"
#include <complex>
using namespace std;
#include "matlabio.h"

/* Update the matlab array. This ensures no mistakes are made when indexing the matlab array.

  mArray - pointer to the memory used to store the matlab array
  arrayVal - the value which should be inserted into mArray[i,j]
  i - row number of value to update
  j - column number of value to update
  nrows - the number of rows in the matrix. This is equivalent to the number of
          x elements under my standard of storing the xy data.

  */

void updateMatlabArray(double *mArray, double arrayVal, int i, int j, int nrows){
	*(mArray + j*nrows + i) = arrayVal;
}

/* Get the matlab array element. This ensures no mistakes are made when indexing the matlab array.

  mArray - pointer to the memory used to store the matlab array
  arrayVal - the value which should be inserted into mArray[i,j]
  i - row number of value to update
  j - column number of value to update
  nrows - the number of rows in the matrix. This is equivalent to the number of
          x elements under my standard of storing the xy data.

  */

double getMatlabArray(double *mArray, int i, int j, int nrows){
	return *(mArray + j*nrows + i);
}

/* Update the 3-dimensional matlab array. This ensures no mistakes are made when indexing the matlab array.

  mArray - pointer to the memory used to store the matlab array
  arrayVal - the value which should be inserted into mArray[i,j]
  i - row number of value to update
  j - column number of value to update
  k - index of third dimension
  nrows - the number of rows in the matrix. This is equivalent to the number of
          x elements under my standard of storing the xy data.
  ncols - the number of coumns in the matrix. This is equivalent to the number of 
          elements in the y direction under my standard of storing xy data

  */

void updateMatlab3DArray(double *mArray, double arrayVal, int i, int j, int k, int nrows, int ncols){
	*(mArray +k*nrows*ncols+ j*nrows + i) = arrayVal;
}

/* Gets an element from the 3-dimensional matlab array. This ensures no mistakes are made when indexing the matlab array.

  mArray - pointer to the memory used to store the matlab array
  i - row number of value to update
  j - column number of value to update
  k - index of third dimension
  nrows - the number of rows in the matrix. This is equivalent to the number of
          x elements under my standard of storing the xy data.
  ncols - the number of coumns in the matrix. This is equivalent to the number of 
          elements in the y direction under my standard of storing xy data

  */

double getMatlab3DArray(double *mArray, int i, int j, int k, int nrows, int ncols){
	return *(mArray +k*nrows*ncols+ j*nrows + i);
}

#ifdef MIOFLAG
void MatlabDisplayDouble(char *desc,double x){
	mxArray *prhs[2];
	char buf[50];
	sprintf(buf, "%s: %.5f\n",desc,x);
	prhs[0] = mxCreateDoubleMatrix(1,1, mxREAL);
	*mxGetPr(prhs[0]) = 1;
	prhs[1] = mxCreateString(buf);
	mexCallMATLAB(0, NULL, 2, prhs ,"fprintf");

}
#endif
#ifdef MIOFLAG
void MatlabDisplayComplex(char *desc,complex<double> x){
	mxArray *prhs[2];
	char buf[50];
	sprintf(buf, "%s: %f+i%f\n",desc,real(x),imag(x));
	prhs[0] = mxCreateDoubleMatrix(1,1, mxREAL);
	*mxGetPr(prhs[0]) = 1;
	prhs[1] = mxCreateString(buf);
	mexCallMATLAB(0, NULL, 2, prhs ,"fprintf");
}
#endif
#ifdef MIOFLAG
void MatlabDisplayString(char *desc,char *outstring){
	mxArray *prhs[2];
	char buf[50];
	sprintf(buf, "%s: %s\n",desc,outstring);
	prhs[0] = mxCreateDoubleMatrix(1,1, mxREAL);
	*mxGetPr(prhs[0]) = 1;
	prhs[1] = mxCreateString(buf);
	mexCallMATLAB(0, NULL, 2, prhs ,"fprintf");
}
#endif
#ifdef MIOFLAG
void MatlabDisplayString(char *desc,char outchar){
	mxArray *prhs[2];
	char buf[50];
	sprintf(buf, "%s: %c\n",desc,outchar);
	prhs[0] = mxCreateDoubleMatrix(1,1, mxREAL);
	*mxGetPr(prhs[0]) = 1;
	prhs[1] = mxCreateString(buf);
	mexCallMATLAB(0, NULL, 2, prhs ,"fprintf");
}
#endif
/*
  Checks whether or not to be verbose in reporting
  
  Returns 1 if it is necessary to report verbosely
          0 if no verbose reporting
  Just checks to see if a variable valled verbose is defined, if it is it returns 
  1 if verbose is non zero or 0 if it is zero. If it is not defined, verbose output
  is the default
*/
#ifdef MIOFLAG
int ReportVerbosely(){
	const mxArray *prhs;
	prhs = mexGetVariablePtr("base","verbose");
	if(prhs == NULL)
	  prhs = mexGetVariablePtr( "caller","verbose");
	if(prhs == NULL)
	  prhs = mexGetVariablePtr( "global","verbose");
	if(prhs == NULL){
	  return 0;
	}
	else{
		double *verb = mxGetPr(prhs);
		if( *verb > 0)
			return 1;
		else
			return 0;
	}
}
#endif

/*Casts a 4-dimensional array such that it may be indexed according to the 
  usual array indexing scheme array[l,k,j,i].

  array is a point to a matlab 4 dimensional array
  nrows the number of rows in the array
  ncols the number of columns in the array
  nlayers the number of layers, each of dimension nrows*ncols 
  nblocks the number of blocks, each of dimension nrows*ncols*nlayers 

*/

double ****castMatlab4DArray(double *array, int nrows, int ncols, int nlayers, int nblocks){
  double ****p;
  int j,k,l;

  if( nlayers==0 )
	  nlayers++;
  
  if( nblocks==0 )
	  nblocks++;
  
  p = (double ****)malloc((unsigned) (nblocks*sizeof(double ***)));
  for(l=0; l<nblocks; l++)
    p[l] = (double ***)malloc((unsigned) (nlayers*sizeof(double **)));
  
  for(l=0; l<nblocks; l++)
    for(k =0; k<nlayers;k++)
      p[l][k] = (double **)malloc((unsigned) (ncols*sizeof(double *)));

  for(l=0; l<nblocks; l++)
    for(k =0; k<nlayers;k++)
      for(j =0; j<ncols;j++)
	p[l][k][j] = (array + l*nrows*ncols*nlayers + k*nrows*ncols + j*nrows);
  return p;

}

/*Frees the axilliary memory used by the castMatlab3DArray
 */
void freeCastMatlab4DArray(double ****castArray, int nlayers, int nblocks){
  for(int l =0; l<nblocks; l++)
    for(int k =0; k<nlayers;k++)
      free(castArray[l][k]);

  for(int l =0; l<nblocks; l++)
    free(castArray[l]);
  free(castArray);
}


/*Casts a 3-dimensional array such that it may be indexed according to the 
  usual array indexing scheme array[k,j,i].

  array is a point to a matlab 3 dimensional array
  nrows the number of rows in the array
  ncols the number of columns in the array
  nlayers the number of layers, each of dimension nrows*ncols 

*/

double ***castMatlab3DArray(double *array, int nrows, int ncols, int nlayers){
  double ***p;
  int j,k;

  if( nlayers==0 )
	  nlayers++;

  p = (double ***)malloc((unsigned) (nlayers*sizeof(double **)));
  for(k =0; k<nlayers;k++)
    p[k] = (double **)malloc((unsigned) (ncols*sizeof(double *)));
  
  for(k =0; k<nlayers;k++)
    for(j =0; j<ncols;j++)
      p[k][j] = (array + k*nrows*ncols+ j*nrows);
  return p;

}

/*Frees the axilliary memory used by the castMatlab3DArray
 */
void freeCastMatlab3DArray(double ***castArray, int nlayers){
  for(int k =0; k<nlayers;k++)
    free(castArray[k]);
  free(castArray);
}

/*Casts a 2-dimensional array such that it may be indexed according to the 
  usual array indexing scheme array[j,i].

  array is a point to a matlab 2 dimensional array
  nrows the number of rows in the array
  ncols the number of columns in the array

*/

double **castMatlab2DArray(double *array, int nrows, int ncols){
  double **p;
  int j;

  p = (double **)malloc((unsigned) (ncols*sizeof(double *)));
  for(j =0; j<ncols;j++)
      p[j] = (array + j*nrows);
  return p;

}

/*Frees the axilliary memory used by the castMatlab2DArray
 */
void freeCastMatlab2DArray(double **castArray){
  free(castArray);
}
/*Casts a 3-dimensional array such that it may be indexed according to the 
  usual array indexing scheme array[k,j,i].

  array is a point to a matlab 3 dimensional array
  nrows the number of rows in the array
  ncols the number of columns in the array
  nlayers the number of layers, each of dimension nrows*ncols 

*/

unsigned char ***castMatlab3DArrayUint8(unsigned char *array, int nrows, int ncols, int nlayers){
  unsigned char ***p;
  int j,k;
  
  if( nlayers==0 )
    nlayers++;
  
  p = (unsigned char ***)malloc((unsigned) (nlayers*sizeof(unsigned char **)));
  for(k =0; k<nlayers;k++)
    p[k] = (unsigned char **)malloc((unsigned) (ncols*sizeof(unsigned char *)));
  
  for(k =0; k<nlayers;k++)
    for(j =0; j<ncols;j++)
      p[k][j] = (array + k*nrows*ncols+ j*nrows);
  return p;

}

/*Frees the axilliary memory used by the castMatlab3DArray
 */
void freeCastMatlab3DArrayUint8(unsigned char ***castArray, int nlayers){
  for(int k =0; k<nlayers;k++)
    free(castArray[k]);
  free(castArray);
}

/*Casts a 3-dimensional array such that it may be indexed according to the 
  usual array indexing scheme array[k,j,i].

  array is a point to a matlab 3 dimensional array
  nrows the number of rows in the array
  ncols the number of columns in the array
  nlayers the number of layers, each of dimension nrows*ncols 

*/

int ***castMatlab3DArrayInt(int *array, int nrows, int ncols, int nlayers){
  int ***p;
  int j,k;
  
  if( nlayers==0 )
    nlayers++;

  p = (int ***)malloc((unsigned) (nlayers*sizeof(int **)));
  for(k =0; k<nlayers;k++)
    p[k] = (int **)malloc((unsigned) (ncols*sizeof(int *)));
  
  for(k =0; k<nlayers;k++)
    for(j =0; j<ncols;j++)
      p[k][j] = (array + k*nrows*ncols+ j*nrows);
  return p;

}

/*Frees the axilliary memory used by the castMatlab3DArray
 */
void freeCastMatlab3DArrayInt(int ***castArray, int nlayers){
  for(int k =0; k<nlayers;k++)
    free(castArray[k]);
  free(castArray);
}

/*Casts a 2-dimensional array such that it may be indexed according to the 
  usual array indexing scheme array[j,i].

  array is a point to a matlab 2 dimensional array
  nrows the number of rows in the array
  ncols the number of columns in the array

*/

int **castMatlab2DArrayInt(int *array, int nrows, int ncols){
  int **p;
  int j;

  p = (int **)malloc((unsigned) (ncols*sizeof(int *)));
  for(j =0; j<ncols;j++)
      p[j] = (array + j*nrows);
  return p;

}

/*Frees the axilliary memory used by the castMatlab2DArray
 */
void freeCastMatlab2DArrayInt(int **castArray){
  free(castArray);
}
