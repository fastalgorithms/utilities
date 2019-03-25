#ifndef VEC_UTILS_H
#define VEC_UTILS_H


// This file contains the basic math utils to be used across ytri
// and yquad



// MACROS and math utils (todo: move where?)
// 2-norm in R3 only:
#define VEC_LEN(v) sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])

void vec_cross(double *a, double *b, double *c);

void vec_diff(int n, double *x, double *y, double *diff);

void vec_dist(int n, double *x, double *y, double *dist);


#endif
