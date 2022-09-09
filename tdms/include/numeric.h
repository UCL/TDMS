#pragma once
#include <cstdlib>

template<typename T, typename S>
void construct_3d_array(T ****E, S I, S J, S K){

  *E = (T ***)malloc(K*sizeof(T *));

  for(S k=0;k<K;k++){
    *(*E+k) = (T **)malloc(J*sizeof(T *));
  }

  for(S k=0;k<K;k++){
    for(S j=0;j<J;j++){
      *(*(*E+k)+j) = (T *)malloc(I*sizeof(T));
    }
  }
}

template<typename T, typename S>
void destroy_3D_array(T ****E, S J, S K){

  for(S k=0;k<K;k++){
    for(S j=0;j<J;j++){
      free(*(*(*E+k)+j));
    }
    free(*(*E+k));
  }
  free(*E);
}
