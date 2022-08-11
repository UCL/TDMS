#include "numeric.h"
#include "field.h"

using namespace std;


void SplitField::allocate() {

  for (auto component : {xy, xz, yx, yz, zx, zy}){
    construct3dArray(reinterpret_cast<double ****>(component), I_tot + 1, J_tot + 1, K_tot + 1);
  }
}

SplitField::SplitField(int I_total, int J_total, int K_total) {

  I_tot = I_total;
  J_tot = J_total;
  K_tot = K_total;
}

void SplitField::zero() {

  for (auto component : {xy, xz, yx, yz, zx, zy})
    for (int k = 0; k < (K_tot + 1); k++)
      for (int j = 0; j < (J_tot + 1); j++)
        for (int i = 0; i < (I_tot + 1); i++) {
          (*reinterpret_cast<double ****>(component))[k][j][i] = 0.;
        }
}

SplitField::~SplitField() {

  for (auto component : {xy, xz, yx, yz, zx, zy}){
    if (component != nullptr && !has_no_elements()) {
      destroy3DArray(reinterpret_cast<double ****>(component), J_tot+1, K_tot+1);
    }
  }
}
