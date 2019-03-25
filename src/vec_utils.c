#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "cprini.h"
#include "vec_utils.h"

void vec_diff(int n, double *x, double *y, double *diff) {

  int i;
  for (i=0; i<n; i++) {
    diff[i] = x[i] - y[i];
  }
  
  return;
}




void vec_dist(int n, double *x, double *y, double *dist) {

  double dx, dtot;
  dtot = 0;
  
  int i;
  for (i=0; i<n; i++) {
    dx = x[i]-y[i];
    dtot = dtot + dx*dx;
  }
  
  *dist = sqrt( dtot);

  return;
}




void vec_cross(double *a, double *b, double *c) {
  // cross product: c = a x b
  c[0] = a[1]*b[2]-a[2]*b[1];
  c[1] = a[2]*b[0]-a[0]*b[2];
  c[2] = a[0]*b[1]-a[1]*b[0];
}


