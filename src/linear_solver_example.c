#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_permutation.h>
#include "fem2d.h"

int main()
{
	int i, j, s;

	gsl_matrix *A;
	gsl_vector *x, *B;
	A = gsl_matrix_alloc(4, 4);
	B = gsl_vector_alloc(4);
	x = gsl_vector_calloc(4);

	double a_data[] = { 3, 0.60, 0.57, 0.96,
                      0.41, 3, 0.99, 0.58,
                      0.14, 0.30, 3, 0.66,
                      0.51, 0.13, 0.19, 3 };

	double b_data[] = { 1.0, 2.0, 3.0, 4.0 };

	for( i = 0; i < 4; i++)
	{
		for(j = 0; j < 4; j++)
		{
			gsl_matrix_set(A, i, j, a_data[4*i+j]);
		}
	}

	for(i = 0;i < 4; i++)
	{
		gsl_vector_set(B, i, b_data[i]);
		gsl_vector_set(x, i, 1);
	}


	/*	
	printf("Exact solution\n");
	gsl_permutation *p = gsl_permutation_alloc(4);

	gsl_linalg_LU_decomp(A, p, &s); //modifies matrix 'A'
	gsl_linalg_LU_solve(A, p, B, x);

	gsl_vector_fprintf(stdout, x, "%g");
	gsl_permutation_free(p);
	*/
	/*	
	printf("Gauss Seidel Result\n");
	gsl_linalg_gauss_seidel(A, x, B, 4, 10e-7);
	gsl_vector_fprintf(stdout, x, "%g");
	*/

	int itr;
	double exact_solution[] = {-0.266783, 0.216177, 0.699326, 1.32503};
	itr = gsl_linalg_BiCGStab(A, x, B, 4, 10e-15);

	printf("Exact solution\n");
	for(i = 0; i < 4; i++)
	{
		printf("%g\n", exact_solution[i]);
	}
	printf("BiCGStab results\n");
	printf("Iteration converged in %d\n", itr);
	gsl_vector_fprintf(stdout, x, "%g");

	gsl_matrix_free(A);
	gsl_vector_free(B);
	gsl_vector_free(x);

	return 0;
}





