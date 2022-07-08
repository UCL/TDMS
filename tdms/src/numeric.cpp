/*****************************************************************
 *
 *  Project.....:  isotropic FDTD code
 *  Application.:  very basic data structure initialisation and
 *                 destruction
 *  Module......:  numeric.cpp
 *  Compiler....:  g++
 *  Written by..:  Peter Munro, Imperial College London, 2002-2008
 *  Environment.:  Linux
 *  Modified....:  Numerous times
 *
 ******************************************************************/

/*---------------------------------------------------------------*/
//                        INCLUDE section
/*---------------------------------------------------------------*/

#include "math.h"
#include "stdlib.h"

void cons3dArray(double ****E, int I, int J, int K){
  *E = (double ***)malloc(K*sizeof(double *));
  for(int k=0;k<K;k++){
    *(*E+k) = (double **)malloc(J*sizeof(double *));
  }

  for(int k=0;k<K;k++){
    for(int j=0;j<J;j++){
      *(*(*E+k)+j) = (double *)malloc(I*sizeof(double));
    }
  }
}

void destroy3DArray(double ****E, int J, int K){
  for(int k=0;k<K;k++){
    for(int j=0;j<J;j++){
      free(*(*(*E+k)+j));
    }
  }

  for(int k=0;k<K;k++){
    free(*(*E+k));
  }

  

}


