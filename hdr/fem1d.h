#ifndef FEM1D_H
#define FEM1D_H 1

//defining tolerence for residuals etc.
extern const double TOL;

int gsl_linalg_gauss_seidel(gsl_matrix *A,gsl_vector *x,gsl_vector *b, int matrixDim, double e);


#endif
