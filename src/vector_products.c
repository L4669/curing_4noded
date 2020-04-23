#include <math.h>

void vector_cross_product(const double *u, const double *v, double *product)
{
    product[0] = u[1]*v[2]- u[2]*v[1];
    product[1] = u[2]*v[0]- u[0]*v[2];
    product[2] = u[0]*v[1]- u[1]*v[0];

}

double vector_dot_product(const double *u, const double *v)
{
    double product = u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
	return product;
}

double vector_mod(double *vector)
{
	double value = sqrt(vector[0]*vector[0] + vector[1]*vector[1]
			+ vector[2]*vector[2]);
	return value;
}

