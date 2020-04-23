#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>
#include "fem2d.h"
#include "config.h"

//def of global variables 
const double K = 0.5, rho = 1000, cp = 2100;
const int MAX_ELEMENT = 1152, MAX_NODE = 1248;
const double T_inf = 25, h = 100, q_cf = 100, Q = 1e6;
const int dirich_nodes = 121, neu_edges = 24, neu_edges_cf = 24;
const double theta = 0.5, delta_t = 1, maxTime = 1000;

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

	char *fname = "./ext/MAPPED_MESH_ANSYS_PROB.hmascii";
	
	struct node_list *node_arr = (struct node_list *)
		malloc(sizeof(struct node_list) * MAX_NODE);
    struct ele_list *ele_arr = (struct ele_list *)
		malloc(sizeof(struct ele_list) * MAX_ELEMENT);
	
	gsl_matrix *K_global_NM, *M_global_NM, *MK_LHS, *MK_RHS;
	gsl_vector *F_global_N, *T_global_N, *F_RHS, *T_global_old_N;
	
	K_global_NM = gsl_matrix_calloc(MAX_NODE, MAX_NODE);
	F_global_N = gsl_vector_calloc(MAX_NODE);
	T_global_N = gsl_vector_calloc(MAX_NODE);
	T_global_old_N = gsl_vector_calloc(MAX_NODE);
	
	node_connectivity(node_arr, ele_arr, fname);

	//stiffness matrix formation (local to global)
	global_stiffness_matrix(node_arr, ele_arr, K_global_NM);

	printf("Siffness Matrix Formation Over\n");


	//Imposing Neumann BC
	char *tag2 = "CONVECTIVE BC"; //bottom face
	struct neumann_boundary_edgelist *nb_array;

	nb_array = (struct neumann_boundary_edgelist *)
		malloc(sizeof(struct neumann_boundary_edgelist)*neu_edges);

	neumann_boundary_edges(nb_array, ele_arr, node_arr, fname, tag2, neu_edges);
	global_neumann_bc(nb_array, K_global_NM, F_global_N, node_arr, ele_arr, 
			T_inf, h, neu_edges);
	/*
	char *tag3 = "CONVECTIVE BC-2"; //circular face
	neumann_boundary_edges(nb_array, ele_arr, node_arr, fname, tag3, 
			neu_edges_circular);
	global_neumann_bc(nb_array, K_global_NM, F_global_N, node_arr, ele_arr,
			T_inf_circular, h_circular, neu_edges_circular);
	*/

	//constant flux neumann bc	
	char *tag3 = "CONSTANT HEAT FLUX";
	struct neumann_boundary_edgelist *nb_array_cf;

	nb_array_cf = (struct neumann_boundary_edgelist *)
		malloc(sizeof(struct neumann_boundary_edgelist)*neu_edges_cf);

	neumann_boundary_edges(nb_array_cf, ele_arr, node_arr, fname, tag3, neu_edges_cf);
	global_neumann_bc_cf(nb_array_cf, F_global_N, node_arr, ele_arr, neu_edges_cf);
	
	printf("Neumann Boundary Formation Over\n");


	//Imposing constant & uniform heat source
	global_heat_source_matrix(node_arr, ele_arr, F_global_N);

	//Extracting dirichlet BC
	char *tag1 = "FIXED TEMP";
	struct dirichlet_boundary_nodelist *db_array;

	db_array = (struct dirichlet_boundary_nodelist *)
		malloc(sizeof(struct dirichlet_boundary_nodelist)*dirich_nodes);

	dirichlet_boundary_nodes(db_array, fname, tag1);
	

/* TRANSIENT SOLVER BEGIN */
#ifdef TRANSIENT
	double time = 0;
	M_global_NM = gsl_matrix_calloc(MAX_NODE, MAX_NODE);
	MK_LHS = gsl_matrix_calloc(MAX_NODE, MAX_NODE);
	MK_RHS = gsl_matrix_calloc(MAX_NODE, MAX_NODE);
	F_RHS = gsl_vector_calloc(MAX_NODE);
	
	//global mass matrix formation
	global_mass_matrix(node_arr, ele_arr, M_global_NM);
	gsl_matrix_scale(M_global_NM, rho*cp);
	printf("Mass Matrix formation over\n");
	
	//assuming uniform temperature throughout the domain at t = 0
	double T_initial = 37;
	for(n = 0; n < MAX_NODE; n++)
	{
		gsl_vector_set(T_global_N, n, T_initial);
	}
	

	//LHS matrix formation
	gsl_matrix_memcpy(MK_LHS, K_global_NM);
	gsl_matrix_scale(MK_LHS, theta*delta_t);
	gsl_matrix_add(MK_LHS, M_global_NM);

	//RHS matrix formation
	gsl_matrix_memcpy(MK_RHS, K_global_NM);
	gsl_matrix_scale(MK_RHS, (theta-1)*delta_t);
	gsl_matrix_add(MK_RHS, M_global_NM);

	//memory not used in further computation
	gsl_matrix_free(K_global_NM);
	gsl_matrix_free(M_global_NM);

#ifdef DIRECT_METHOD
	int s;
	gsl_permutation *p = gsl_permutation_alloc(MAX_NODE);
	while(time < maxTime)
	{
		gsl_vector_memcpy(F_RHS, F_global_N);
		gsl_blas_dgemv(CblasNoTrans, 1.0, MK_RHS, T_global_N, delta_t, F_RHS);
		//imposing dirichlet BC
		global_dirichlet_bc(db_array, MK_LHS, F_RHS);

		gsl_linalg_LU_decomp(MK_LHS, p, &s);
		gsl_linalg_LU_solve(MK_LHS, p, F_RHS, T_global_N);
		printf("T at (0.6,0.2) @ %g sec = %g\n", time, 
				gsl_vector_get(T_global_N, 26));
		time += delta_t;
	}
	gsl_permutation_free(p);
	printf("T at (0.6,0.2) = %g\n", gsl_vector_get(T_global_N, 26));
	printf("Tmax = %g\n", gsl_vector_max(T_global_N));
	printf("Tmin = %g\n", gsl_vector_min(T_global_N));
#endif

#ifdef INDIRECT_METHOD
	while(time <= maxTime)
	{
		gsl_vector_memcpy(F_RHS, F_global_N);
		gsl_blas_dgemv(CblasNoTrans, 1.0, MK_RHS, T_global_N, delta_t, F_RHS);
		//imposing dirichlet BC
		global_dirichlet_bc(db_array, MK_LHS, F_RHS);
	//	printf("asd\n");

		status = gsl_linalg_BiCGStab(MK_LHS, T_global_N, F_RHS, 10e-7);
		//printf("T at (0.1,0.3) @ %g sec = %g\n", time, 
		//		gsl_vector_get(T_global_N, 26));
		time += delta_t;
	}	
	//printf("T at left boundary middle = %g\n", gsl_vector_get(T_global_N, 176));
	printf("Tmax = %g\n", gsl_vector_max(T_global_N));
	printf("Tmin = %g\n", gsl_vector_min(T_global_N));
#endif

	gsl_matrix_free(MK_LHS);
	gsl_matrix_free(MK_RHS);
	gsl_vector_free(F_RHS);
#endif
/* TRANSIENT SOLVER END */


/* STEADY STATE SOLVER BEGIN */
#ifdef STEADY_STATE
	
	global_dirichlet_bc(db_array, K_global_NM, F_global_N);

#ifdef DIRECT_METHOD	
	// Direct method to solve system of linear equations
	int s;
	gsl_permutation *p = gsl_permutation_alloc(MAX_NODE);
	gsl_linalg_LU_decomp(K_global_NM, p, &s);
	gsl_linalg_LU_solve(K_global_NM, p, F_global_N, T_global_N);

	printf("T at (0.6,0.2) = %g\n", gsl_vector_get(T_global_N, 20));
	printf("Tmax = %g\n", gsl_vector_max(T_global_N));
	printf("Tmin = %g\n", gsl_vector_min(T_global_N));
	gsl_permutation_free(p);

#endif

#ifdef INDIRECT_METHOD 	
	//Indirect method to solve systems of linear equation
	double T_initial = 100;
	for(n = 0; n < MAX_NODE; n++)
	{
		gsl_vector_set(T_global_N, n, T_initial);
	}
	status = gsl_linalg_BiCGStab(K_global_NM, T_global_N, F_global_N, 10e-7);

	if(status > 0)
	{
		printf("Iteration converged in %d\n", status);
		//printf("T at (0.6,0.2) = %g\n", gsl_vector_get(T_global_N, 20));
		printf("Tmax = %g\n", gsl_vector_max(T_global_N));
		printf("Tmin = %g\n", gsl_vector_min(T_global_N));
	}
	else
	{
		printf("Error in iterative solver module\n");
	}
	
	gsl_matrix_free(K_global_NM);
	gsl_vector_free(F_global_N);

#endif

#endif
/* STEADY STATE SOLVER END */
	

#ifdef FILE_OUTPUT
	printf("Writing to file...\n");
	FILE *fp;
	//wrting gnuplot file 
	fp = fopen("./output/test_case_4.dat", "w+");
	for(n = 0; n < MAX_NODE; n++)
	{
		fprintf(fp, "%g %g %g\n", node_arr[n].x, node_arr[n].y, 
				gsl_vector_get(T_global_N, n));
	}

	fclose(fp);

	//writing VTK file output
	fp  = fopen("./output/MAPPED_MESH_ANSYS_PROB.vtk", "w+");
	//writing VTK headers
	fprintf(fp, "%s\n", "# vtk DataFile Version 4.2");
	fprintf(fp, "%s\n", "MAPPED_MESH_ANSYS_PROB Data File");
	fprintf(fp, "%s\n", "ASCII");
	fprintf(fp, "%s\n", "DATASET UNSTRUCTURED_GRID");
	//fprintf(fp, "%s\n", "POINTS MAX_NODE double");
	fprintf(fp, "%s\n", "POINTS 1248 double");
	for(n = 0; n < MAX_NODE; n++)
	{
		fprintf(fp, "%g %g %g\n", node_arr[n].x, node_arr[n].y, 0);
	}

	//fprintf(fp, "%s\n", "CELLS MAX_ELEMENT MAX_ELE*5");
	fprintf(fp, "%s\n", "CELLS 1152 5760");
	for(n = 0; n < MAX_ELEMENT; n++)
	{
		fprintf(fp, "%d %d %d %d %d\n", 4, ele_arr[n].ele_node_no[0]-1,
				ele_arr[n].ele_node_no[1]-1,
				ele_arr[n].ele_node_no[2]-1,
				ele_arr[n].ele_node_no[3]-1);
	}
	fprintf(fp, "%s\n", "CELL_TYPES 1152");
	for(n = 0; n < MAX_ELEMENT; n++)
	{
		fprintf(fp, "%d\n", 9);
	}

	fprintf(fp, "%s\n", "POINT_DATA 1248");
	fprintf(fp, "%s\n", "SCALARS temperature double");
	fprintf(fp, "%s\n", "LOOKUP_TABLE default");
	for(n = 0; n < MAX_NODE; n++)
	{
		fprintf(fp, "%g\n", gsl_vector_get(T_global_N, n));
	}
	fclose(fp);

#endif

	gsl_vector_free(T_global_N);
	free(db_array);
	//free(nb_array_cf);
	free(nb_array);
	free(node_arr); free(ele_arr);

	return 0;
}
