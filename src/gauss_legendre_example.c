#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gauss_legendre.h"

#ifndef PI
	#define PI 3.1415926535897932384626433832795028841971693993751
#endif
#ifndef FABS
	#define FABS(a) ((a)>=0?(a):-(a))
#endif


double f(double x, double y,void* data)
{
	double *test, *coeff;
	double I;
	coeff = (double *)data;

	test = (double *)malloc(sizeof(double)*3);

	test[0] = coeff[0]; test[1] = coeff[1]; test[2] = coeff[2]; 

	I = x*y + test[0] + test[1];
	free(test);
	return I;
}

int main(int argc, char* argv[])
{

	/* numerical approximation of integral */
	double approx;		

	/* true value of int(sin(x), x=0..Pi) = 2.0*/
	double exact = (2.0/3)+20.0+(2.0/5); 

	int i = 20;
	double *coeff;

	coeff = (double *)malloc(sizeof(double)*3);
	coeff[0] = 1; coeff[1] = 1; coeff[2] = 1;

	approx = gauss_legendre_2D_cube(5, f, coeff,0, 1, 0, 1);
	
	double error = FABS(approx-exact);
	
	printf("error = %.15g\n",approx);

}

