#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include "fem1d.h"

const double TOL = 10e-7;

int main()
{
	int no_of_elements, element, i, j;
	double x_l, x_u, h, M_ij, M_local_ij;
	gsl_matrix *M, *M_local; // M is the global mass matrix, M_local is per element

	no_of_elements = 10;
	x_l = 0; x_u = 0.1;
	h = (x_u - x_l)/no_of_elements;

	//Mass Matrix
	
	//Matrix memory allocation and initialization
	M = gsl_matrix_calloc(no_of_elements+1, no_of_elements+1);
	M_local = gsl_matrix_alloc(2,2);
	
	//See Eq.3.71 for values of local element mass matrix 
	gsl_matrix_set(M_local, 0, 0, 2);
	gsl_matrix_set(M_local, 0, 1, 1);
	gsl_matrix_set(M_local, 1, 0, 1);
	gsl_matrix_set(M_local, 1, 1, 2);

	gsl_matrix_scale(M_local, (const double)h/6);


	for(element = 0; element < no_of_elements; element++)
	{
		for(i = 0; i < 2; i++)
		{
			for(j = 0;j <2; j++)
			{
				M_local_ij = gsl_matrix_get(M_local, i, j);
				M_ij = gsl_matrix_get(M, i+element, j+element);
				gsl_matrix_set(M, i+element, j+element, M_local_ij+M_ij) ;
			}
		}
	}

	//gsl_matrix_scale(M, (const double)120);
	/*
	printf("Mass Matrix\n");
	gsl_matrix_scale(M, (const double)120);
	gsl_matrix_fprintf(stdout, M, "%g");
	*/

	//Conductance Matrix
	
	double alpha, K_ij, K_local_ij;
	/*
	alpha = (double *)malloc(sizeof(double)*no_of_elements);
	alpha[0] = 1.125*pow(10,-5);
	alpha[1] = 1.375*pow(10,-5);
	*/
	alpha = 1.125*10e-5;

	//Matrix declaration and memory allocation
	gsl_matrix *K, *K_local;

	K = gsl_matrix_calloc(no_of_elements+1, no_of_elements+1);
	K_local = gsl_matrix_alloc(2,2);
	
	//See Eq.3.72 for values of local element conductance matrix
	gsl_matrix_set(K_local, 0, 0, 1);
    gsl_matrix_set(K_local, 0, 1, -1);
    gsl_matrix_set(K_local, 1, 0, -1);
	gsl_matrix_set(K_local, 1, 1, 1); 
	
	gsl_matrix_scale(K_local, (const double)1/h);
	
	for(element = 0; element < no_of_elements; element++)
	{
		for(i = 0; i < 2; i++)
		{
			for(j = 0;j < 2; j++)
			{
				K_local_ij = alpha*gsl_matrix_get(K_local, i, j);
				K_ij = gsl_matrix_get(K, i+element, j+element);
				gsl_matrix_set(K, i+element, j+element, K_local_ij+K_ij) ;
			}
		}
	}
	
	double h_bar = 2.5*10e-5;
	double T_inf = 400; //deg. celcius
	double T_L = 39.18;
	//Implementing flux bc at x=0
	K_ij = gsl_matrix_get(K, 0, 0) + h_bar;
	gsl_matrix_set(K, 0, 0, K_ij);
	
	//gsl_matrix_scale(K, (const double)120);
	/*
	printf("Conductance Matrix\n");
	gsl_matrix_scale(K, (const double)120);
	gsl_matrix_fprintf(stdout, K, "%g");
	*/

	//load vector F
	gsl_vector *F;
	F = gsl_vector_calloc(no_of_elements+1);
	gsl_vector_set(F, 0, h_bar*T_inf);
	gsl_vector_set(F, no_of_elements, T_L);
	
	//making super global matrix see Eq.3.109
	
	double theta = 1;
	double dt = 100;
	gsl_matrix *A, *B;
	gsl_vector *C;
	A = gsl_matrix_calloc(no_of_elements+1, no_of_elements+1);
	B = gsl_matrix_calloc(no_of_elements+1, no_of_elements+1);
	C = gsl_vector_alloc(no_of_elements+1);

	gsl_matrix_memcpy(A, K);
	gsl_matrix_scale(A, (const double)theta*dt);
	gsl_matrix_add(A, M);

	//gsl_matrix_fprintf(stdout, A, "%g");
	
	gsl_matrix_add(B, K);
	gsl_matrix_scale(B, (const double)(theta-1)*dt);
	gsl_matrix_add(B, M);

	//gsl_matrix_fprintf(stdout, B, "%g");

	gsl_vector *T_n, *T_o;
	T_n = gsl_vector_alloc(no_of_elements+1);
	T_o = gsl_vector_alloc(no_of_elements+1);

	//initial condition
	for(i = 0; i <= no_of_elements; i++)
	{
		gsl_vector_set(T_o, i, 39.18);
	}

	gsl_vector_memcpy(C, F);
	gsl_blas_dgemv(CblasNoTrans, 1.0, B, T_o, dt, C);

	gsl_vector_fprintf(stdout, C, "%g");

	//implementing dirichlet BC at x=L
	gsl_vector_set(C, no_of_elements, T_L);
	for(i = 0; i < no_of_elements; i++)
	{
		gsl_matrix_set(A, no_of_elements, i, 0);
	}
	gsl_matrix_set(A, no_of_elements, no_of_elements, 1);

	//solving the system of equations with gauss_seidel routine
	gsl_vector *x;
	x = gsl_vector_calloc(no_of_elements+1);

	for(i = 0;i < no_of_elements; i++)
	{
		gsl_vector_set(x, i, 1);
	}

	gsl_linalg_gauss_seidel(A, x, C, no_of_elements+1, TOL);

	printf("\n");
	gsl_vector_fprintf(stdout, x, "%g");
	

	gsl_matrix_free(M);
	gsl_matrix_free(M_local);
	gsl_matrix_free(K);
	gsl_matrix_free(K_local);
	gsl_matrix_free(A);
	gsl_matrix_free(B);
	gsl_vector_free(C);
	gsl_vector_free(F);
	gsl_vector_free(T_o);
	gsl_vector_free(T_n);
	gsl_vector_free(x);

	return 0;
}

