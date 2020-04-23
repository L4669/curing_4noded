#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>

#include "fem2d.h"
#include "config.h"

//def of global variables 
const int MAX_ELEMENT = 200, MAX_NODE = 231;
const int neu_edges = 20, dirich_nodes = 0;
const double theta = 0.5, delta_t = 0.5, maxTime = 195*60;

//printing gsl_matrix array in m x n format at stdout
void matrix_printf_hr(int order, gsl_matrix *mat)
{
	int n, m;

	for(n = 0; n < order; n++)
	{
		for(m = 0; m < order; m++)
		{
			printf("%.4g ", gsl_matrix_get(mat, n , m));
		}
		printf("\n");
	}
}


int main()
{
	
	int status, n;
	double T_inf, maxResidual = 1;

	char *fname = "./ext/ACTUAL_PROBLEM_20x10.hmascii";
	
	struct node_list *node_arr = (struct node_list *)
		malloc(sizeof(struct node_list) * MAX_NODE);
    	struct ele_list *ele_arr = (struct ele_list *)
		malloc(sizeof(struct ele_list) * MAX_ELEMENT);
	struct neumann_boundary_edgelist *nb_array;
		nb_array = (struct neumann_boundary_edgelist *)
			malloc(sizeof(struct neumann_boundary_edgelist)*neu_edges);

	
	gsl_matrix *K_global_NM, *M_global_NM, *MK_LHS, *MK_RHS, *dR_dT;
	gsl_vector *F_global_N, *F_global_N_i, *F_global_N_im1, *F_global_N_const,
			   *T_global_N, *F_RHS, *T_global_N_i, *T_global_N_im1, *R;
	gsl_vector *CD_global_E, *CDR_global_E,
			   *CD_global_E_i, *CDR_global_E_i, *dT;

	
	T_global_N = gsl_vector_calloc(MAX_NODE);
	T_global_N_i = gsl_vector_calloc(MAX_NODE);
	T_global_N_im1 = gsl_vector_calloc(MAX_NODE);
	dT = gsl_vector_calloc(MAX_NODE);
	CD_global_E = gsl_vector_calloc(MAX_NODE);
	CD_global_E_i = gsl_vector_calloc(MAX_NODE);
	CDR_global_E = gsl_vector_calloc(MAX_NODE);
	CDR_global_E_i = gsl_vector_calloc(MAX_NODE);
	
	node_connectivity(node_arr, ele_arr, fname);



/* TRANSIENT SOLVER BEGIN */
#ifdef TRANSIENT
	double time = delta_t;
	MK_LHS = gsl_matrix_calloc(MAX_NODE, MAX_NODE);
	MK_RHS = gsl_matrix_calloc(MAX_NODE, MAX_NODE);
	F_RHS = gsl_vector_calloc(MAX_NODE);
	F_global_N = gsl_vector_calloc(MAX_NODE);
	F_global_N_const = gsl_vector_calloc(MAX_NODE);
	F_global_N_i = gsl_vector_calloc(MAX_NODE);
	F_global_N_im1 = gsl_vector_calloc(MAX_NODE);
	R = gsl_vector_calloc(MAX_NODE);
	dR_dT = gsl_matrix_calloc(MAX_NODE, MAX_NODE);
	

	//assuming uniform temperature throughout the domain at t = 0
	double T_initial = 25;
	for(n = 0; n < MAX_NODE; n++)
	{
		gsl_vector_set(T_global_N, n, T_initial);
		gsl_vector_set(T_global_N_i, n, T_initial);
		gsl_vector_set(dT, n, 1);
	}
	

#ifdef INDIRECT_METHOD
	
	//stiffness matrix formation (local to global)
	K_global_NM = gsl_matrix_calloc(MAX_NODE, MAX_NODE);
	//printf("Siffness Matrix Formation Over\n");

	//global mass matrix formation
	M_global_NM = gsl_matrix_calloc(MAX_NODE, MAX_NODE);
	//printf("Mass Matrix formation over\n");

	while(time <= maxTime)
	{

		gsl_matrix_set_zero(K_global_NM);
		gsl_matrix_set_zero(M_global_NM);
		global_stiffness_matrix(node_arr, ele_arr, K_global_NM);
		global_mass_matrix(node_arr, ele_arr, M_global_NM);
		gsl_matrix_scale(M_global_NM, rho*cp);

		//printf("asd\n");
		//Imposing Neumann BC
		gsl_vector_set_zero(F_global_N);
		
		T_inf = autoclave_input_cycle(time);
		//printf("%g %g\n", time, T_inf);

		char *tag2 = "TOP ROBIN BC"; //top face
		neumann_boundary_edges(nb_array, ele_arr, node_arr, fname, tag2, neu_edges);
		double h = nb_array[0].value*K_zz;
		global_neumann_bc(nb_array, K_global_NM, F_global_N, node_arr, ele_arr, 
				T_inf, h, neu_edges);
		
		char *tag3 = "BOTTOM ROBIN BC"; //bottom face
		neumann_boundary_edges(nb_array, ele_arr, node_arr, fname, tag3, 
				neu_edges);
		h = nb_array[0].value*K_zz;
		global_neumann_bc(nb_array, K_global_NM, F_global_N, node_arr, ele_arr,
				T_inf, h, neu_edges);
		
		gsl_matrix_memcpy(MK_LHS, K_global_NM);
		gsl_matrix_memcpy(MK_RHS, K_global_NM);
		//free(K_global_NM);
		//LHS matrix formation
		gsl_matrix_scale(MK_LHS, theta*delta_t);
		gsl_matrix_add(MK_LHS, M_global_NM);

		//RHS matrix formation
		gsl_matrix_scale(MK_RHS, (theta-1)*delta_t);
		gsl_matrix_add(MK_RHS, M_global_NM);

		//free(M_global_NM);

		//printf("%g %g\n", time, gsl_vector_max(CD_global_E));
		gsl_vector_memcpy(F_global_N_const, F_global_N);
		global_cure_rate_vector(T_global_N, CD_global_E, CD_global_E_i,
				CDR_global_E_i, ele_arr);

		global_heat_source_matrix(node_arr, ele_arr, F_global_N, CDR_global_E_i);
		
		gsl_vector_memcpy(F_RHS, F_global_N);

		gsl_blas_dgemv(CblasNoTrans, 1.0, MK_RHS, T_global_N, 
				delta_t, F_RHS);

		//first iteration with explicit method
		gsl_linalg_BiCGStab(MK_LHS, T_global_N_i, F_RHS, 10e-7);
		//printf("%g %g %g\n", time, gsl_vector_max(T_global_N_i), 
		//		gsl_vector_min(T_global_N_i));

	//	return 0;

		
		global_cure_rate_vector(T_global_N_i, CD_global_E, CD_global_E_i,
				CDR_global_E_i, ele_arr);
	
		gsl_vector_memcpy(F_global_N_i, F_global_N_const);
		global_heat_source_matrix(node_arr, ele_arr, F_global_N_i, CDR_global_E_i);
		
		gsl_vector_memcpy(F_RHS, F_global_N_i);
		gsl_blas_dgemv(CblasNoTrans, 1.0, MK_RHS, T_global_N, 
				delta_t*(1-theta), F_RHS);
		
		//gsl_vector_memcpy(CD_global_E, CD_global_E_i);
		gsl_vector_memcpy(T_global_N_im1, T_global_N);
		gsl_vector_memcpy(F_global_N_i, F_global_N);
		
		while(1)
		{

			gsl_vector_set_all(dT, 1);
			//gsl_vector_memcpy(CD_global_E_im1, CD_global_E_i);
			gsl_vector_memcpy(F_global_N_im1, F_global_N_i);
			gsl_vector_memcpy(F_global_N_i, F_global_N_const);

			global_cure_rate_vector(T_global_N_i, CD_global_E, CD_global_E_i, 
					CDR_global_E_i, ele_arr);
			//Imposing heat source
			global_heat_source_matrix(node_arr, ele_arr, F_global_N_i, CDR_global_E_i);
	
			//gsl_vector_memcpy(F_RHS, F_global_N_i);
			//gsl_blas_dgemv(CblasNoTrans, 1.0, MK_RHS, T_global_N_i, 
			//		delta_t*(1- theta), F_RHS);
			gsl_vector_memcpy(R, F_global_N_i);
			gsl_blas_dgemv(CblasNoTrans, 1.0, MK_LHS, T_global_N_i, -delta_t*theta, R);
			gsl_vector_sub(R, F_RHS);

			//printf("%g %g\n", time, gsl_vector_get(T_global_N_i, 134));
			double maxR = fabs(gsl_vector_max(R));
			double minR = fabs(gsl_vector_min(R));

			if(maxR > minR)
			{
				maxResidual = maxR;
			}
			else
			{
				maxResidual = minR;
			}
			
			if(maxResidual < 1e-5)
			{
				//printf("%g %g\n", time, maxResidual);
				break;
			}
			//printf("%g %g\n", time, maxResidual);
			
			
			gsl_vector_sub(F_global_N_im1, F_global_N_i);
			gsl_vector_sub(T_global_N_im1, T_global_N_i);
			status = gsl_vector_div(F_global_N_im1, T_global_N_im1);
			if(status)
			{
				printf("asd\n");
			}
			gsl_vector_scale(F_global_N_im1, -theta*delta_t);
			gsl_vector_memcpy(T_global_N_im1, T_global_N_i);

			residual_derivative_matrix(dR_dT, MK_LHS, F_global_N_im1);
			gsl_vector_scale(R, -1);
			gsl_linalg_BiCGStab(dR_dT, dT, R, 10e-7);
			
			//printf("%g %g\n", time, gsl_vector_get(dT, 134));

			gsl_vector_add(T_global_N_i, dT);	
		}

		gsl_vector_memcpy(CD_global_E, CD_global_E_i);
		gsl_vector_memcpy(CDR_global_E, CDR_global_E_i);
		gsl_vector_memcpy(T_global_N, T_global_N_i);
		
		printf("%g %g\n", time, gsl_vector_get(T_global_N, 134));
		time += delta_t;
		maxResidual = 1;

	}	
	//printf("Tmax = %g\n", gsl_vector_max(T_global_N));
	//printf("Tmin = %g\n", gsl_vector_min(T_global_N));
#endif

	gsl_matrix_free(MK_LHS);
	gsl_matrix_free(MK_RHS);
	gsl_vector_free(F_RHS);
	gsl_matrix_free(M_global_NM);
#endif
/* TRANSIENT SOLVER END */


	

#ifdef FILE_OUTPUT
	printf("Writing to file...\n");
	FILE *fp;
	//wrting gnuplot file 
	fp = fopen("./output/ACTUAL_PROBLEM_20x10.dat", "w+");
	for(n = 0; n < MAX_NODE; n++)
	{
		fprintf(fp, "%g %g %g\n", node_arr[n].x, node_arr[n].y, 
				gsl_vector_get(T_global_N, n));
	}

	fclose(fp);

	//writing VTK file output
	fp  = fopen("./output/ACTUAL_PROBLEM_20x10.vtk", "w+");
	//writing VTK headers
	fprintf(fp, "%s\n", "# vtk DataFile Version 4.2");
	fprintf(fp, "%s\n", "ACTUAL_PROBLEM_20x10 Data File");
	fprintf(fp, "%s\n", "ASCII");
	fprintf(fp, "%s\n", "DATASET UNSTRUCTURED_GRID");
	//fprintf(fp, "%s\n", "POINTS MAX_NODE double");
	fprintf(fp, "%s\n", "POINTS 231 double");
	for(n = 0; n < MAX_NODE; n++)
	{
		fprintf(fp, "%g %g %g\n", node_arr[n].x, node_arr[n].y, 0);
	}

	//fprintf(fp, "%s\n", "CELLS MAX_ELEMENT MAX_ELE*5");
	fprintf(fp, "%s\n", "CELLS 200 1000");
	for(n = 0; n < MAX_ELEMENT; n++)
	{
		fprintf(fp, "%d %d %d %d %d\n", 4, ele_arr[n].ele_node_no[0]-1,
				ele_arr[n].ele_node_no[1]-1,
				ele_arr[n].ele_node_no[2]-1,
				ele_arr[n].ele_node_no[3]-1);
	}
	fprintf(fp, "%s\n", "CELL_TYPES 200");
	for(n = 0; n < MAX_ELEMENT; n++)
	{
		fprintf(fp, "%d\n", 9);
	}

	fprintf(fp, "%s\n", "POINT_DATA 231");
	fprintf(fp, "%s\n", "SCALARS temperature double");
	fprintf(fp, "%s\n", "LOOKUP_TABLE default");
	for(n = 0; n < MAX_NODE; n++)
	{
		fprintf(fp, "%g\n", gsl_vector_get(T_global_N, n));
	}
	fclose(fp);

#endif

	gsl_vector_free(T_global_N);
	gsl_vector_free(T_global_N_i);
	gsl_vector_free(T_global_N_im1);
	gsl_vector_free(CD_global_E); gsl_vector_free(CDR_global_E);
	free(nb_array);
	free(node_arr); free(ele_arr);

	return 0;
}
