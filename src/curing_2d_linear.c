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
const double theta = 0.5, maxTime = 275*60;
double delta_t = 0.5; 

int main()
{
	
	int status, n;
	double T_inf, maxResidual = 1;

	char *fname = "./ext/ACTUAL_PROBLEM_20x10_HM_GEOM.hmascii";
	
	struct node_list *node_arr = (struct node_list *)
		malloc(sizeof(struct node_list) * MAX_NODE);
    struct ele_list *ele_arr = (struct ele_list *)
		malloc(sizeof(struct ele_list) * MAX_ELEMENT);
	struct neumann_boundary_edgelist *nb_array_1, *nb_array_2;
		nb_array_1 = (struct neumann_boundary_edgelist *)
			malloc(sizeof(struct neumann_boundary_edgelist)*neu_edges);

		nb_array_2 = (struct neumann_boundary_edgelist *)
			malloc(sizeof(struct neumann_boundary_edgelist)*neu_edges);
	
	gsl_matrix *K_global_NM, *M_global_NM, *MK_LHS, *MK_RHS;
	gsl_vector *F_global_N, *T_global_N, *F_RHS;
	gsl_vector *CD_global_E, *CDR_global_E, *CD_global_E_i; 
	gsl_vector *F_global_N_const, *F_global_N_i, *F_global_N_im1, 
			   *T_global_N_i, *T_global_N_im1, *nulll;
	
	K_global_NM = gsl_matrix_calloc(MAX_NODE, MAX_NODE);
	M_global_NM = gsl_matrix_calloc(MAX_NODE, MAX_NODE);
	T_global_N = gsl_vector_calloc(MAX_NODE);
	T_global_N_i = gsl_vector_calloc(MAX_NODE);
	CD_global_E = gsl_vector_calloc(MAX_NODE);
	CD_global_E_i = gsl_vector_calloc(MAX_NODE);
	CDR_global_E = gsl_vector_calloc(MAX_NODE);
	F_global_N = gsl_vector_calloc(MAX_NODE);
	nulll = gsl_vector_calloc(MAX_NODE);
	
	node_connectivity(node_arr, ele_arr, fname);


	//assuming uniform temperature throughout the domain at t = 0
	//absolute temperature giving erroneous results ??
	double T_initial = 20 + 273, h1 = 0, h2 = 0;

	gsl_vector_set_all(T_global_N, T_initial);
	gsl_vector_set_all(CD_global_E, 1e-5);

	char *tag2 = "TOP ROBIN BC"; //top face
	neumann_boundary_edges(nb_array_1, ele_arr, node_arr, fname, tag2, neu_edges);
	h1 = nb_array_1[0].value*K_zz;

	char *tag3 = "BOTTOM ROBIN BC"; //bottom face
	neumann_boundary_edges(nb_array_2, ele_arr, node_arr, fname, tag3, neu_edges);
	h2 = nb_array_2[0].value*K_zz;

/* TRANSIENT SOLVER BEGIN */
#ifdef TRANSIENT
	double time = 0;
	MK_LHS = gsl_matrix_calloc(MAX_NODE, MAX_NODE);
	MK_RHS = gsl_matrix_calloc(MAX_NODE, MAX_NODE);
	F_RHS = gsl_vector_calloc(MAX_NODE);
	
/*
	//LHS matrix formation
	gsl_matrix_memcpy(MK_LHS, K_global_NM);
	gsl_matrix_scale(MK_LHS, theta*delta_t);
	gsl_matrix_add(MK_LHS, M_global_NM);

	//RHS matrix formation
	gsl_matrix_memcpy(MK_RHS, K_global_NM);
	gsl_matrix_scale(MK_RHS, (theta-1)*delta_t);
	gsl_matrix_add(MK_RHS, M_global_NM);
*/
/*
#ifdef DIRECT_METHOD
    int s;
    gsl_permutation *p = gsl_permutation_alloc(MAX_NODE);
	gsl_linalg_LU_decomp(MK_LHS, p, &s);
    while(time < 0.5)
    {
        gsl_vector_memcpy(F_RHS, F_global_N);
        gsl_blas_dgemv(CblasNoTrans, 1.0, MK_RHS, T_global_N, delta_t, F_RHS);

        gsl_linalg_LU_solve(MK_LHS, p, F_RHS, T_global_N_i);
		gsl_vector_memcpy(T_global_N, T_global_N_i);
		gsl_vector_set_zero(T_global_N_i);
        time += delta_t;
    }
    gsl_permutation_free(p);
    printf("T at (0.6,0.2) = %g\n", gsl_vector_get(T_global_N, 134));
    printf("Tmax = %g\n", gsl_vector_max(T_global_N));
    printf("Tmin = %g\n", gsl_vector_min(T_global_N));
#endif
*/

#ifdef INDIRECT_METHOD
	FILE *fp1;
	fp1 = fopen("data_steady_10.dat", "w+");
	while(time <= maxTime)
	{
	/*	
		if(time >= 100*60)
		{
			delta_t = 0.7;
		}
	*/	
		
		T_inf = autoclave_input_cycle_2(time) + 273;

		gsl_matrix_set_zero(K_global_NM);
		gsl_matrix_set_zero(M_global_NM);
		gsl_vector_set_zero(F_global_N);

		global_stiffness_matrix(node_arr, ele_arr, K_global_NM);
		global_mass_matrix(node_arr, ele_arr, M_global_NM);
		gsl_matrix_scale(M_global_NM, rho*cp);

		global_neumann_bc(nb_array_1, K_global_NM, F_global_N, node_arr, ele_arr, 
				T_inf, h1, neu_edges);
			
		global_neumann_bc(nb_array_2, K_global_NM, F_global_N, node_arr, ele_arr,
				T_inf, h2, neu_edges);
		

		//LHS matrix formation
		gsl_matrix_memcpy(MK_LHS, K_global_NM);
		gsl_matrix_scale(MK_LHS, theta*delta_t);
		gsl_matrix_add(MK_LHS, M_global_NM);

		//RHS matrix formation
		gsl_matrix_memcpy(MK_RHS, K_global_NM);
		gsl_matrix_scale(MK_RHS, (theta-1)*delta_t);
		gsl_matrix_add(MK_RHS, M_global_NM);

		global_cure_rate_vector(T_global_N, CD_global_E, CD_global_E_i,
                CDR_global_E, ele_arr, nulll, time);
        global_heat_source_matrix(node_arr, ele_arr, F_global_N, CDR_global_E,
				nulll, nulll);

		gsl_vector_memcpy(F_RHS, F_global_N);

		// for linear solution approach uncomment the following two steps
		gsl_blas_dgemv(CblasNoTrans, 1.0, MK_RHS, T_global_N, 
				delta_t, F_RHS);
		status = gsl_linalg_BiCGStab(MK_LHS, T_global_N, F_RHS, 10e-7);

		
		
		gsl_vector_memcpy(CD_global_E, CD_global_E_i);
		printf("%g %g %g\n", time/60, gsl_vector_get(T_global_N, 134)-273,
				gsl_vector_get(CD_global_E, 51));
		fprintf(fp1, "%g %g %g %g\n", time/60, gsl_vector_get(T_global_N, 134)-273,
				gsl_vector_get(CD_global_E, 51), gsl_vector_get(CDR_global_E, 51));

		fflush(fp1);
		time += delta_t;

	}	
	//printf("Tmax = %g\n", gsl_vector_max(T_global_N));
	//printf("Tmin = %g\n", gsl_vector_min(T_global_N));
#endif
	fclose(fp1);
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
		fprintf(fp, "%g %g %d\n", node_arr[n].x, node_arr[n].y, 0);
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
	gsl_vector_free(CD_global_E); gsl_vector_free(CDR_global_E);
	//free(nb_array);
	free(node_arr); free(ele_arr);

	return 0;
}
