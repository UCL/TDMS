#pragma once
#include "matlabio.h"

class Dimensions{
private:
  int i = 0;
  int j = 0;
  int k = 0;

  bool are_nd(int n) const{return (bool(i) + bool(j) + bool(k)) == n;}


public:
  int operator[] (int value) const;

  explicit Dimensions(const mxArray* ptr);

  bool are_1d() const{return are_nd(1);}
  bool are_2d() const{return are_nd(2);}
};
