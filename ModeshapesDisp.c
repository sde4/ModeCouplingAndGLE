/*
 * This function finds the mode combinations which lead to non-zero coupling constant
 */

#include "struct_def.h"

void ModeshapesDisp(gsl_vector* A, int mord , int nord, disc_const dc){
  int i, j;
  double dat, x, y;

  for (j = 0; j < dc.Nfy; j++){
    for (i = 0; i < dc.Nfx; i++){
        y = 0.0+(j)*dc.hy;
        x = 0.0+(i)*dc.hx;
        dat= sin(mord*PI*x/dc.Lx)*sin(nord*PI*y/dc.Ly);
        gsl_vector_set(A, i + (j)*dc.Nfx, dat);
    }
  }
  
  return;
}
