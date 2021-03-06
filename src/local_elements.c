#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "gauss_legendre.h"
#include "fem2d.h"



/*2D Bilinear Quadrilateral Element*/

double stiffness_matrix_nm_bqe(double zeta, double ita, void* data)
{
	double *alpha;
	alpha = (double *)malloc(sizeof(double)*3);
	double K_nm,
		   A_Nx, B_Nx_1, B_Nx_2, A_Ny, B_Ny_1, B_Ny_2,
		   A_Mx, B_Mx_1, B_Mx_2, A_My, B_My_1, B_My_2;
	double detJ, C_Nx, C_Mx, C_Ny, C_My;
	double *coeff;

	coeff = (double *)data;

	alpha[0] = coeff[0]; alpha[1] = coeff[1]; alpha[2] = coeff[2];

	A_Nx = coeff[3]; B_Nx_1 = coeff[4]; B_Nx_2 = coeff[5];
	A_Mx = coeff[6]; B_Mx_1 = coeff[7]; B_Mx_2 = coeff[8];

	A_Ny = coeff[9]; B_Ny_1 = coeff[10]; B_Ny_2 = coeff[11];
	A_My = coeff[12]; B_My_1 = coeff[13]; B_My_2 = coeff[14];

	//Coordinate Tranformation Jacobian Determinant
	detJ = (alpha[0] + alpha[1]*zeta + alpha[2]*ita)/8;
	//detJ = 1;

	C_Nx = (A_Nx + B_Nx_1*zeta + B_Nx_2*ita)/(8*detJ);
	C_Mx = (A_Mx + B_Mx_1*zeta + B_Mx_2*ita)/(8*detJ);

	C_Ny = (A_Ny + B_Ny_1*zeta + B_Ny_2*ita)/(8*detJ);
	C_My = (A_My + B_My_1*zeta + B_My_2*ita)/(8*detJ);

	//[n][m] element of stiffness matrix
	K_nm = (C_Nx*C_Mx*K_xx + C_Ny*C_My*K_zz)*detJ;


	free(alpha);
	return K_nm;
}

void stiffness_matrix_bqe(double *cord_x, double *cord_y,
		gsl_matrix *K_local_NM)
{
	double *alpha, *A_Nx, *A_Ny, *B_Nx_1, *B_Ny_1, *B_Nx_2, *B_Ny_2;
	double *transfer, K_nm;

	int i, n, m;

	alpha = (double *)malloc(sizeof(double)*3);
	A_Nx = (double *)malloc(sizeof(double)*4);
	A_Ny = (double *)malloc(sizeof(double)*4);
	B_Nx_1 = (double *)malloc(sizeof(double)*4);
	B_Ny_1 = (double *)malloc(sizeof(double)*4);
	B_Nx_2 = (double *)malloc(sizeof(double)*4);
	B_Ny_2 = (double *)malloc(sizeof(double)*4);
	transfer = (double *)malloc(sizeof(double)*15);


	alpha[0] = (cord_x[3] - cord_x[1])*(cord_y[0] - cord_y[2]) -
		(cord_x[0] - cord_x[2])*(cord_y[3] - cord_y[1]);
	alpha[1] = (cord_x[2] - cord_x[3])*(cord_y[0] - cord_y[1]) -
		(cord_x[0] - cord_x[1])*(cord_y[2] - cord_y[3]);
	alpha[2] = (cord_x[3] - cord_x[0])*(cord_y[1] - cord_y[2]) -
		(cord_x[1] - cord_x[2])*(cord_y[3] - cord_y[0]);

	for(i = 0; i < 4; i++)
	{
		A_Nx[i] = cord_y[(i+1)%4] - cord_y[(i+3)%4];
		A_Ny[i] = cord_x[(i+3)%4] - cord_x[(i+1)%4];
	}

	B_Nx_1[0] = (cord_y[3] - cord_y[2]);
	B_Nx_1[2] = (cord_y[0] - cord_y[1]);
	B_Nx_1[1] = -B_Nx_1[0]; B_Nx_1[3] = -B_Nx_1[2];

	B_Ny_1[0] = (cord_x[2] - cord_x[3]);
	B_Ny_1[2] = (cord_x[1] - cord_x[0]);
	B_Ny_1[1] = -B_Ny_1[0]; B_Ny_1[3] = -B_Ny_1[2];

	B_Nx_2[0] = (cord_y[2] - cord_y[1]);
	B_Nx_2[2] = (cord_y[3] - cord_y[0]);
	B_Nx_2[1] = -B_Nx_2[2]; B_Nx_2[3] = -B_Nx_2[0];

	B_Ny_2[0] = (cord_x[1] - cord_x[2]);
	B_Ny_2[2] = (cord_x[0] - cord_x[3]);
	B_Ny_2[1] = -B_Ny_2[2]; B_Ny_2[3] = -B_Ny_2[0];

	for(i = 0; i <= 2; i++)
	{
		transfer[i] = alpha[i];
	}

	//for different n,m values of stiffness matrix
	for(n = 0; n < 4; n++)
	{
		for(m = 0; m< 4; m++)
		{

			transfer[3] = A_Nx[n]; transfer[4] = B_Nx_1[n];
			transfer[5] = B_Nx_2[n];

			transfer[6] = A_Nx[m]; transfer[7] = B_Nx_1[m];
			transfer[8] = B_Nx_2[m];

			transfer[9] = A_Ny[n]; transfer[10] = B_Ny_1[n];
			transfer[11] = B_Ny_2[n];

			transfer[12] = A_Ny[m]; transfer[13] = B_Ny_1[m];
			transfer[14] = B_Ny_2[m];


		K_nm = gauss_legendre_2D_cube
			(5, stiffness_matrix_nm_bqe, transfer, -1, 1, -1, 1);

		gsl_matrix_set(K_local_NM, n, m, K_nm);
		}
	}

	free(alpha);
	free(A_Nx); free(A_Ny); free(B_Nx_1);
	free(B_Nx_2); free(B_Ny_1); free(B_Ny_2);
	free(transfer);
}



double mass_matrix_nm_bqe(double zeta, double ita, void *data)
{
	int i;
	double N_n, N_m, M_nm, zeta_n, ita_n,
		   zeta_m, ita_m, detJ;
	double *coeff, *alpha;

	coeff = (double *)data;
	alpha = (double *)malloc(sizeof(double)*3);

	alpha[0] = coeff[0]; alpha[1] = coeff[1]; alpha[2] = coeff[2];
	zeta_n = coeff[3]; ita_n = coeff[4];
	zeta_m = coeff[5]; ita_m = coeff[6];

	N_n = (1 + zeta_n*zeta)*(1 + ita_n*ita)/4;
	N_m = (1 + zeta_m*zeta)*(1 + ita_m*ita)/4;

	detJ = (alpha[0] + alpha[1]*zeta + alpha[2]*ita)/8;

	M_nm = N_n*N_m*detJ;

	free(alpha);

	return M_nm;
}

void mass_matrix_bqe(double *cord_x, double *cord_y,
		gsl_matrix *M_local_NM)
{

	int n, m;
	double M_nm;
	double *alpha, *transfer, *zeta_n, *ita_n;

	alpha = (double *)malloc(sizeof(double)*3);
	zeta_n = (double *)malloc(sizeof(double)*4);
	ita_n = (double *)malloc(sizeof(double)*4);
	transfer = (double *)malloc(sizeof(double)*7);

	alpha[0] = (cord_x[3] - cord_x[1])*(cord_y[0] - cord_y[2]) -
		(cord_x[0] - cord_x[2])*(cord_y[3] - cord_y[1]);
	alpha[1] = (cord_x[2] - cord_x[3])*(cord_y[0] - cord_y[1]) -
		(cord_x[0] - cord_x[1])*(cord_y[2] - cord_y[3]);
	alpha[2] = (cord_x[3] - cord_x[0])*(cord_y[1] - cord_y[2]) -
		(cord_x[1] - cord_x[2])*(cord_y[3] - cord_y[0]);

	zeta_n[0] = -1; zeta_n[1] = 1; zeta_n[2] = 1; zeta_n[3] = -1;
	ita_n[0] = -1; ita_n[1] = -1; ita_n[2] = 1; ita_n[3] = 1;

	for(n = 0; n <= 2; n++)
	{
		transfer[n] = alpha[n];
	}

	for(n = 0; n < 4; n++)
	{
		for(m = 0; m < 4; m++)
		{
			transfer[3] = zeta_n[n]; transfer[4] = ita_n[n];
			transfer[5] = zeta_n[m]; transfer[6] = ita_n[m];

			M_nm = gauss_legendre_2D_cube(5, mass_matrix_nm_bqe,
					transfer, -1, 1, -1, 1);
			gsl_matrix_set(M_local_NM, n, m, M_nm);
		}
	}

	free(alpha); free(transfer); free(zeta_n); free(ita_n);
}

//Generating BC Matrix of Type I
double boundary_matrix_n_bqe_type1(double xi, void *data)
{
	double F_n, xi_n;
	double *coeff, line_jacobian;

	coeff = (double *)data;
	line_jacobian = coeff[0];
	xi_n = coeff[1];

	F_n = (1 + xi_n*xi)*line_jacobian/2;

	return F_n;
}


void boundary_matrix_bqe_type1(int boundary, double *x_cord, double *y_cord,
		gsl_vector *F_local_N)

{
	double F_n, *transfer, line_jacobian, b_x, b_y, c_x, c_y;
	transfer = (double *)malloc(sizeof(double)*2);

	if(boundary == 12)
	{
		b_x = -x_cord[0] + x_cord[1] + x_cord[2] - x_cord[3];
		b_y = -y_cord[0] + y_cord[1] + y_cord[2] - y_cord[3];
		line_jacobian = sqrt(b_x*b_x + b_y*b_y)/4;
		transfer[0] = line_jacobian;

		transfer[1] = -1;
		F_n = gauss_legendre(5, boundary_matrix_n_bqe_type1,
				transfer, -1, 1);
		gsl_vector_set(F_local_N, 0, F_n);

		transfer[1] = 1;
		F_n = gauss_legendre(5, boundary_matrix_n_bqe_type1,
				transfer, -1, 1);
		gsl_vector_set(F_local_N, 1, F_n);

		gsl_vector_set(F_local_N, 2, 0);

		gsl_vector_set(F_local_N, 3, 0);

	}
	else if(boundary == 23)
	{
		c_x = -x_cord[0] - x_cord[1] + x_cord[2] + x_cord[3];
		c_y = -y_cord[0] - y_cord[1] + y_cord[2] + y_cord[3];
		line_jacobian = sqrt(c_x*c_x + c_y*c_y)/4;
		transfer[0] = line_jacobian;

		gsl_vector_set(F_local_N, 0, 0);

		transfer[1] = -1;
		F_n = gauss_legendre(5, boundary_matrix_n_bqe_type1,
				transfer, -1, 1);
		gsl_vector_set(F_local_N, 1, F_n);

		transfer[1] = 1;
		F_n = gauss_legendre(5, boundary_matrix_n_bqe_type1,
				transfer, -1, 1);
		gsl_vector_set(F_local_N, 2, F_n);

		gsl_vector_set(F_local_N, 3, 0);
	}
	else if(boundary == 34)
	{
		b_x = -x_cord[0] + x_cord[1] + x_cord[2] - x_cord[3];
		b_y = -y_cord[0] + y_cord[1] + y_cord[2] - y_cord[3];
		line_jacobian = sqrt(b_x*b_x + b_y*b_y)/4;
		transfer[0] = line_jacobian;

		gsl_vector_set(F_local_N, 0, 0);

		gsl_vector_set(F_local_N, 1, 0);

		transfer[1] = 1;
		F_n = gauss_legendre(5, boundary_matrix_n_bqe_type1,
				transfer, -1, 1);
		gsl_vector_set(F_local_N, 2, F_n);

		transfer[1] = -1;
		F_n = gauss_legendre(5, boundary_matrix_n_bqe_type1,
				transfer, -1, 1);
		gsl_vector_set(F_local_N, 3, F_n);
	}
	else
	{
		c_x = -x_cord[0] - x_cord[1] + x_cord[2] + x_cord[3];
		c_y = -y_cord[0] - y_cord[1] + y_cord[2] + y_cord[3];
		line_jacobian = sqrt(c_x*c_x + c_y*c_y)/4;
		transfer[0] = line_jacobian;

		transfer[1] = -1;
		F_n = gauss_legendre(5, boundary_matrix_n_bqe_type1,
				transfer, -1, 1);
		gsl_vector_set(F_local_N, 0, F_n);

		gsl_vector_set(F_local_N, 1, 0);

		gsl_vector_set(F_local_N, 2, 0);

		transfer[1] = 1;
		F_n = gauss_legendre(5, boundary_matrix_n_bqe_type1,
				transfer, -1, 1);
		gsl_vector_set(F_local_N, 3, F_n);
	}

	free(transfer);
}

//Generating BC Matrix of Type II
double boundary_matrix_n_bqe_type2(double zeta, void *data)
{
	double Kb_nm, ita, ita_n, ita_m, zeta_n, zeta_m;
	double *coeff, line_jacobian;

	coeff = (double *)data;
	line_jacobian = coeff[0];
	zeta_n = coeff[1]; zeta_m = coeff[2];
	ita_n = coeff[3]; ita_m = coeff[4];
	ita = coeff[5];

	Kb_nm = (1 + zeta_n*zeta)*(1 + zeta_m*zeta)*
		(1 + ita_n*ita)*(1 + ita_m*ita)*line_jacobian/16;

	return Kb_nm;
}

double boundary_matrix_m_bqe_type2(double ita, void *data)
{
	double Kb_nm, zeta, zeta_n, zeta_m, ita_n, ita_m;
	double *coeff, line_jacobian;

	coeff = (double *)data;
	line_jacobian = coeff[0];
	zeta_n = coeff[1]; zeta_m = coeff[2];
	ita_n = coeff[3]; ita_m = coeff[4];
	zeta = coeff[5];

	Kb_nm = (1 + ita_n*ita)*(1 + ita_m*ita)*(1 + zeta_n*zeta)
		*(1 + zeta_m*zeta)*line_jacobian/16;


	return Kb_nm;
}

void boundary_matrix_bqe_type2(int boundary, double *x_cord, double *y_cord,
		gsl_matrix *Kb_local_NM)
{
	int n, m;
	double Kb_nm, *transfer, ita_n[4], zeta_n[4], ita, zeta,
		   line_jacobian, b_x, b_y, c_x, c_y;
	transfer = (double *)malloc(sizeof(double)*6);

	zeta_n[0] = -1; zeta_n[1] = 1; zeta_n[2] = 1; zeta_n[3] = -1;
	ita_n[0] = -1; ita_n[1] = -1; ita_n[2] = 1; ita_n[3] = 1;

	if(boundary == 12)
	{
		ita = -1;

		b_x = -x_cord[0] + x_cord[1] + x_cord[2] - x_cord[3];
		b_y = -y_cord[0] + y_cord[1] + y_cord[2] - y_cord[3];
		line_jacobian = sqrt(b_x*b_x + b_y*b_y)/4;
		transfer[0] = line_jacobian;

		for(n = 0; n < 4; n++)
		{
			for(m = 0; m < 4; m++)
			{
				transfer[1] = zeta_n[n]; transfer[2] = zeta_n[m];
				transfer[3] = ita_n[n]; transfer[4] = ita_n[m];
				transfer[5] = ita;
				Kb_nm = gauss_legendre(5, boundary_matrix_n_bqe_type2,
				transfer, -1, 1);
				gsl_matrix_set(Kb_local_NM, n, m, Kb_nm);
			}
		}
	}
	else if(boundary == 23)
	{
		zeta = 1;

		c_x = -x_cord[0] - x_cord[1] + x_cord[2] + x_cord[3];
		c_y = -y_cord[0] - y_cord[1] + y_cord[2] + y_cord[3];
		line_jacobian = sqrt(c_x*c_x + c_y*c_y)/4;
		transfer[0] = line_jacobian;

		for(n = 0; n < 4; n++)
		{
			for(m = 0; m < 4; m++)
			{
				transfer[1] = zeta_n[n]; transfer[2] = zeta_n[m];
				transfer[3] = ita_n[n]; transfer[4] = ita_n[m];
				transfer[5] = zeta;
				Kb_nm = gauss_legendre(5, boundary_matrix_m_bqe_type2,
				transfer, -1, 1);
				gsl_matrix_set(Kb_local_NM, n, m, Kb_nm);
			}
		}
	}
	else if(boundary == 34)
	{
		ita = 1;

		b_x = -x_cord[0] + x_cord[1] + x_cord[2] - x_cord[3];
		b_y = -y_cord[0] + y_cord[1] + y_cord[2] - y_cord[3];
		line_jacobian = sqrt(b_x*b_x + b_y*b_y)/4;
		transfer[0] = line_jacobian;

		for(n = 0; n < 4; n++)
		{
			for(m = 0; m < 4; m++)
			{
				transfer[1] = zeta_n[n]; transfer[2] = zeta_n[m];
				transfer[3] = ita_n[n]; transfer[4] = ita_n[m];
				transfer[5] = ita;
				Kb_nm = gauss_legendre(5, boundary_matrix_n_bqe_type2,
				transfer, -1, 1);
				gsl_matrix_set(Kb_local_NM, n, m, Kb_nm);
			}
		}
	}
	else
	{
		zeta = -1;

		c_x = -x_cord[0] - x_cord[1] + x_cord[2] + x_cord[3];
		c_y = -y_cord[0] - y_cord[1] + y_cord[2] + y_cord[3];
		line_jacobian = sqrt(c_x*c_x + c_y*c_y)/4;
		transfer[0] = line_jacobian;

		for(n = 0; n < 4; n++)
		{
			for(m = 0; m < 4; m++)
			{
				transfer[1] = zeta_n[n]; transfer[2] = zeta_n[m];
				transfer[3] = ita_n[n]; transfer[4] = ita_n[m];
				transfer[5] = zeta;
				Kb_nm = gauss_legendre(5, boundary_matrix_m_bqe_type2,
				transfer, -1, 1);
				gsl_matrix_set(Kb_local_NM, n, m, Kb_nm);
			}
		}
	}

	free(transfer);
}

//Constant heat source vector

double heat_source_matrix_n_bqe_1(double zeta, double ita, void *data)
{
	double N_n, F_n, zeta_n, ita_n, detJ;
	double *coeff, *alpha;

	coeff = (double *)data;
	alpha = (double *)malloc(sizeof(double)*3);

	alpha[0] = coeff[0]; alpha[1] = coeff[1]; alpha[2] = coeff[2];
	zeta_n = coeff[3]; ita_n = coeff[4];

	N_n = (1 + zeta_n*zeta)*(1 + ita_n*ita)/4;

	detJ = (alpha[0] + alpha[1]*zeta + alpha[2]*ita)/8;

  	int i;
    double dalpha_dt_temp=0, alpha_gauss = 0;
    double b,dalpha_dt_gauss = 0, *shape_fun, *gauss_temp;
    shape_fun = (double *)malloc(sizeof(double)*4);
    gauss_temp = (double *)malloc(sizeof(double)*4);
	
	gauss_temp[0] = coeff[5]; gauss_temp[1] = coeff[6]; 
	gauss_temp[2] = coeff[7]; gauss_temp[3] = coeff[8];
	
 
    shape_fun[0] = (1-zeta)*(1-ita)*0.25;
    shape_fun[1] = (1+zeta)*(1-ita)*0.25;
    shape_fun[2] = (1+zeta)*(1+ita)*0.25;
    shape_fun[3] = (1-zeta)*(1+ita)*0.25;

    for(i=0; i<4; i++)
    {
        dalpha_dt_temp += shape_fun[i]*gauss_temp[i];
		alpha_gauss += shape_fun[i]*coeff[9+i];
    }
	b = A_c*exp(-E_c/(GSL_CONST_MKSA_MOLAR_GAS*dalpha_dt_temp));
	dalpha_dt_gauss = b*pow(alpha_gauss, m_c)*pow(1-alpha_gauss, n_c);

	    

	F_n = N_n*detJ*dalpha_dt_gauss;

	free(alpha);
	free(shape_fun);
	free(gauss_temp);

	return F_n;
}

void heat_source_matrix_bqe_1(double *cord_x, double *cord_y,
		gsl_vector *F_local_N, double *gauss_temp, double *alpha_gauss)
{

	int n;
	double F_n;
	double *alpha, *transfer, *zeta_n, *ita_n;

	alpha = (double *)malloc(sizeof(double)*3);
	zeta_n = (double *)malloc(sizeof(double)*4);
	ita_n = (double *)malloc(sizeof(double)*4);
	transfer = (double *)malloc(sizeof(double)*13);

	alpha[0] = (cord_x[3] - cord_x[1])*(cord_y[0] - cord_y[2]) -
		(cord_x[0] - cord_x[2])*(cord_y[3] - cord_y[1]);
	alpha[1] = (cord_x[2] - cord_x[3])*(cord_y[0] - cord_y[1]) -
		(cord_x[0] - cord_x[1])*(cord_y[2] - cord_y[3]);
	alpha[2] = (cord_x[3] - cord_x[0])*(cord_y[1] - cord_y[2]) -
		(cord_x[1] - cord_x[2])*(cord_y[3] - cord_y[0]);

	zeta_n[0] = -1; zeta_n[1] = 1; zeta_n[2] = 1; zeta_n[3] = -1;
	ita_n[0] = -1; ita_n[1] = -1; ita_n[2] = 1; ita_n[3] = 1;

	for(n = 0; n <= 2; n++)
	{
		transfer[n] = alpha[n];
	}
	
	transfer[5] = gauss_temp[0]; 
	transfer[6] = gauss_temp[1]; 
	transfer[7] = gauss_temp[2]; 
	transfer[8] = gauss_temp[3];
	transfer[9] = alpha_gauss[0];  
	transfer[10] = alpha_gauss[1];  
	transfer[11] = alpha_gauss[2];  
	transfer[12] = alpha_gauss[3];  

	for(n = 0; n < 4; n++)
	{
		transfer[3] = zeta_n[n]; transfer[4] = ita_n[n];

		F_n = gauss_legendre_2D_cube(5, heat_source_matrix_n_bqe_1,
					transfer, -1, 1, -1, 1);
		gsl_vector_set(F_local_N, n, F_n);
	}
	free(alpha); free(transfer); free(zeta_n); free(ita_n);
}
/*
double heat_source_matrix_n_bqe_2(double zeta, double ita, void *data)
{
	double N_n, F_n, zeta_n, ita_n, detJ;
	double *coeff, *alpha;

	coeff = (double *)data;
	alpha = (double *)malloc(sizeof(double)*3);

	alpha[0] = coeff[0]; alpha[1] = coeff[1]; alpha[2] = coeff[2];
	zeta_n = coeff[3]; ita_n = coeff[4];

	N_n = (1 + zeta_n*zeta)*(1 + ita_n*ita)/4;

	detJ = (alpha[0] + alpha[1]*zeta + alpha[2]*ita)/8;

  	int i;
    double dalpha_dt_temp=0, alpha_gauss = 0;
    double b,dalpha_dt_gauss = 0, *shape_fun, *gauss_temp;
    shape_fun = (double *)malloc(sizeof(double)*4);
    gauss_temp = (double *)malloc(sizeof(double)*4);
	
	gauss_temp[0] = coeff[5]; gauss_temp[1] = coeff[6]; 
	gauss_temp[2] = coeff[7]; gauss_temp[3] = coeff[8];
	
 
    shape_fun[0] = (1-zeta)*(1-ita)*0.25;
    shape_fun[1] = (1+zeta)*(1-ita)*0.25;
    shape_fun[2] = (1+zeta)*(1+ita)*0.25;
    shape_fun[3] = (1-zeta)*(1+ita)*0.25;

    for(i=0; i<4; i++)
    {
        dalpha_dt_temp += shape_fun[i]*gauss_temp[i];
		alpha_gauss += shape_fun[i]*coeff[9+i];
    }
	b = A_c*exp(-E_c/(GSL_CONST_MKSA_MOLAR_GAS*dalpha_dt_temp))*
		(E_c/(dalpha_dt_temp*dalpha_dt_temp*GSL_CONST_MKSA_MOLAR_GAS));

	dalpha_dt_gauss = b*pow(alpha_gauss, m_c)*pow(1-alpha_gauss, n_c);

	
	F_n = N_n*detJ*dalpha_dt_gauss;

	free(alpha);
	free(shape_fun);
	free(gauss_temp);

	return F_n;
}
*/
double heat_source_matrix_n_bqe_2(double zeta, double ita, void *data)
{
	double N_n, F_n, zeta_n, ita_n, detJ;
	double *coeff, *alpha;

	coeff = (double *)data;
	alpha = (double *)malloc(sizeof(double)*3);

	alpha[0] = coeff[0]; alpha[1] = coeff[1]; alpha[2] = coeff[2];
	zeta_n = coeff[3]; ita_n = coeff[4];

	N_n = (1 + zeta_n*zeta)*(1 + ita_n*ita)/4;

	detJ = (alpha[0] + alpha[1]*zeta + alpha[2]*ita)/8;

  	int i;
    double dalpha_dt_temp=0, alpha_gauss = 0, b1, b2, b3, dt_dT_gauss = 0;
    double b,dalpha_dt_gauss = 0, *shape_fun, *gauss_temp, d_dalpha_dt_gauss = 0;
    shape_fun = (double *)malloc(sizeof(double)*4);
    gauss_temp = (double *)malloc(sizeof(double)*4);
	
	gauss_temp[0] = coeff[5]; gauss_temp[1] = coeff[6]; 
	gauss_temp[2] = coeff[7]; gauss_temp[3] = coeff[8];
	
 
    shape_fun[0] = (1-zeta)*(1-ita)*0.25;
    shape_fun[1] = (1+zeta)*(1-ita)*0.25;
    shape_fun[2] = (1+zeta)*(1+ita)*0.25;
    shape_fun[3] = (1-zeta)*(1+ita)*0.25;

    for(i=0; i<4; i++)
    {
        dalpha_dt_temp += shape_fun[i]*gauss_temp[i];
		alpha_gauss += shape_fun[i]*coeff[9+i];
		dt_dT_gauss += shape_fun[i]*coeff[13+i];
    }
	b = A_c*exp(-E_c/(GSL_CONST_MKSA_MOLAR_GAS*dalpha_dt_temp));
	dalpha_dt_gauss = b*pow(alpha_gauss, m_c)*pow(1-alpha_gauss, n_c);

	b1 = E_c/(dalpha_dt_temp*dalpha_dt_temp*GSL_CONST_MKSA_MOLAR_GAS);
	b2 = m_c/alpha_gauss;
	b3 = n_c/(alpha_gauss - 1);

	d_dalpha_dt_gauss = b1*dalpha_dt_gauss + (b2 + b3)*pow(dalpha_dt_gauss, 2)*
		dt_dT_gauss;
   
	//d_dalpha_dt_gauss = b1*dalpha_dt_gauss;

	F_n = N_n*detJ*d_dalpha_dt_gauss;

	free(alpha);
	free(shape_fun);
	free(gauss_temp);

	return F_n;
}

void heat_source_matrix_bqe_2(double *cord_x, double *cord_y,
		gsl_vector *F_local_N, double *gauss_temp, double *alpha_gauss, double *dt_dT)
{

	int n;
	double F_n;
	double *alpha, *transfer, *zeta_n, *ita_n;

	alpha = (double *)malloc(sizeof(double)*3);
	zeta_n = (double *)malloc(sizeof(double)*4);
	ita_n = (double *)malloc(sizeof(double)*4);
	transfer = (double *)malloc(sizeof(double)*17);

	alpha[0] = (cord_x[3] - cord_x[1])*(cord_y[0] - cord_y[2]) -
		(cord_x[0] - cord_x[2])*(cord_y[3] - cord_y[1]);
	alpha[1] = (cord_x[2] - cord_x[3])*(cord_y[0] - cord_y[1]) -
		(cord_x[0] - cord_x[1])*(cord_y[2] - cord_y[3]);
	alpha[2] = (cord_x[3] - cord_x[0])*(cord_y[1] - cord_y[2]) -
		(cord_x[1] - cord_x[2])*(cord_y[3] - cord_y[0]);

	zeta_n[0] = -1; zeta_n[1] = 1; zeta_n[2] = 1; zeta_n[3] = -1;
	ita_n[0] = -1; ita_n[1] = -1; ita_n[2] = 1; ita_n[3] = 1;

	for(n = 0; n <= 2; n++)
	{
		transfer[n] = alpha[n];
	}
	
	transfer[5] = gauss_temp[0]; 
	transfer[6] = gauss_temp[1]; 
	transfer[7] = gauss_temp[2]; 
	transfer[8] = gauss_temp[3];
	transfer[9] = alpha_gauss[0];  
	transfer[10] = alpha_gauss[1];  
	transfer[11] = alpha_gauss[2];  
	transfer[12] = alpha_gauss[3]; 
	transfer[13] = dt_dT[0]; transfer[14] = dt_dT[1];
	transfer[15] = dt_dT[2]; transfer[16] = dt_dT[4];

	for(n = 0; n < 4; n++)
	{
		transfer[3] = zeta_n[n]; transfer[4] = ita_n[n];

		F_n = gauss_legendre_2D_cube(5, heat_source_matrix_n_bqe_2,
					transfer, -1, 1, -1, 1);
		gsl_vector_set(F_local_N, n, F_n);
	}
	free(alpha); free(transfer); free(zeta_n); free(ita_n);
}

/*
void heat_source_matrix_bqe_2(double *cord_x, double *cord_y,
		gsl_vector *F_local_N, double *gauss_temp, double *alpha_gauss)
{

	int n;
	double F_n;
	double *alpha, *transfer, *zeta_n, *ita_n;

	alpha = (double *)malloc(sizeof(double)*3);
	zeta_n = (double *)malloc(sizeof(double)*4);
	ita_n = (double *)malloc(sizeof(double)*4);
	transfer = (double *)malloc(sizeof(double)*13);

	alpha[0] = (cord_x[3] - cord_x[1])*(cord_y[0] - cord_y[2]) -
		(cord_x[0] - cord_x[2])*(cord_y[3] - cord_y[1]);
	alpha[1] = (cord_x[2] - cord_x[3])*(cord_y[0] - cord_y[1]) -
		(cord_x[0] - cord_x[1])*(cord_y[2] - cord_y[3]);
	alpha[2] = (cord_x[3] - cord_x[0])*(cord_y[1] - cord_y[2]) -
		(cord_x[1] - cord_x[2])*(cord_y[3] - cord_y[0]);

	zeta_n[0] = -1; zeta_n[1] = 1; zeta_n[2] = 1; zeta_n[3] = -1;
	ita_n[0] = -1; ita_n[1] = -1; ita_n[2] = 1; ita_n[3] = 1;

	for(n = 0; n <= 2; n++)
	{
		transfer[n] = alpha[n];
	}
	
	transfer[5] = gauss_temp[0]; 
	transfer[6] = gauss_temp[1]; 
	transfer[7] = gauss_temp[2]; 
	transfer[8] = gauss_temp[3];
	transfer[9] = alpha_gauss[0];  
	transfer[10] = alpha_gauss[1];  
	transfer[11] = alpha_gauss[2];  
	transfer[12] = alpha_gauss[3];  

	for(n = 0; n < 4; n++)
	{
		transfer[3] = zeta_n[n]; transfer[4] = ita_n[n];

		F_n = gauss_legendre_2D_cube(5, heat_source_matrix_n_bqe_2,
					transfer, -1, 1, -1, 1);
		gsl_vector_set(F_local_N, n, F_n);
	}
	free(alpha); free(transfer); free(zeta_n); free(ita_n);
}
*/

