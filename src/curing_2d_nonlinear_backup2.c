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
const double theta = 0.5, maxTime = 195*60;
double delta_t = 5; 

int main()
{
	
	int status, n;
	double T_inf, maxResidual = 1, scale_const;

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
	
	gsl_matrix *K_global_NM, *M_global_NM, *MK_LHS, *MK_RHS, *dR_dT;
	gsl_vector *F_global_N, *T_global_N, *F_RHS;
	gsl_vector *CD_global_E, *CDR_global_E, *CD_global_E_i; 
	gsl_vector *F_global_N_const, *F_global_N_i, *F_global_N_im1, 
			   *T_global_N_i, *T_global_N_im1, *dT, *R, *CDR_global_E_i,
			   *d_CDR;
	gsl_vector *Null;
	
	K_global_NM = gsl_matrix_calloc(MAX_NODE, MAX_NODE);
	M_global_NM = gsl_matrix_calloc(MAX_NODE, MAX_NODE);
	T_global_N = gsl_vector_calloc(MAX_NODE);
	T_global_N_i = gsl_vector_calloc(MAX_NODE);
	T_global_N_im1 = gsl_vector_calloc(MAX_NODE);
	CD_global_E = gsl_vector_calloc(MAX_NODE);
	CD_global_E_i = gsl_vector_calloc(MAX_NODE);
	CDR_global_E = gsl_vector_calloc(MAX_NODE);
	CDR_global_E_i = gsl_vector_calloc(MAX_NODE);
	F_global_N = gsl_vector_calloc(MAX_NODE);
	F_global_N_i = gsl_vector_calloc(MAX_NODE);
	F_global_N_im1 = gsl_vector_calloc(MAX_NODE);
	F_global_N_const = gsl_vector_calloc(MAX_NODE);
	dT = gsl_vector_calloc(MAX_NODE);
	R = gsl_vector_calloc(MAX_NODE);
	dR_dT = gsl_matrix_calloc(MAX_NODE, MAX_NODE);
	d_CDR = gsl_vector_calloc(MAX_NODE);
	Null = gsl_vector_calloc(MAX_NODE);
	
	node_connectivity(node_arr, ele_arr, fname);


	//assuming uniform temperature throughout the domain at t = 0
	double T_initial = 25+273, h1 = 0, h2 = 0;

	gsl_vector_set_all(T_global_N, T_initial);
	gsl_vector_set_all(CD_global_E, 1e-15);

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
	

#ifdef INDIRECT_METHOD
	double w = 1;
	FILE *fp1;
	fp1 = fopen("./output/data_nl_1.dat", "w+");
	int s;
    gsl_permutation *p = gsl_permutation_alloc(MAX_NODE);
	while(time <= maxTime)
	{
		
		if(time >= 115*60)
		{
			delta_t = 0.5;
		}
		
		T_inf = autoclave_input_cycle(time)+273;
		
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
	
		//printf("%g %g\n", gsl_matrix_max(MK_LHS), gsl_matrix_max(MK_RHS));
		//return -1;

		gsl_vector_memcpy(F_global_N_const, F_global_N);

        global_cure_rate_vector(T_global_N, CD_global_E, CD_global_E_i,
                CDR_global_E_i, ele_arr, Null, time);

        global_heat_source_matrix(node_arr, ele_arr, F_global_N, CDR_global_E_i,
				Null, Null);
        
        gsl_vector_memcpy(F_RHS, F_global_N);

        gsl_blas_dgemv(CblasNoTrans, 1.0, MK_RHS, T_global_N, 
                delta_t*(1-theta), F_RHS);


		int itr = 0;
        gsl_vector_memcpy(T_global_N_i, T_global_N);
		//gsl_vector_set_all(T_global_N_i, T_initial);
		while(1)
        {

			itr++;
            gsl_vector_memcpy(F_global_N_i, F_global_N_const);
			gsl_vector_set_all(F_global_N_im1, 0);

            global_cure_rate_vector_2(T_global_N_i, CD_global_E, CD_global_E_i,
                    CDR_global_E_i, ele_arr, d_CDR, time);
            //Imposing heat source
            global_heat_source_matrix(node_arr, ele_arr, F_global_N_i, 
					CDR_global_E_i, d_CDR, F_global_N_im1);

            gsl_vector_memcpy(R, F_global_N_i);
            gsl_blas_dgemv(CblasNoTrans, 1.0, MK_LHS, T_global_N_i, -delta_t*theta, R);
            gsl_vector_sub(R, F_RHS);

			maxResidual = gsl_blas_dnrm2(R);
            if(maxResidual < 1e-7)
            {
                break;
            }

			gsl_vector_scale(F_global_N_im1, -theta*delta_t);
            residual_derivative_matrix(dR_dT, MK_LHS, F_global_N_im1, ele_arr);
            gsl_vector_scale(R, -1);
			
	/*		
            status = gsl_linalg_BiCGStab(dR_dT, dT, R, 10e-7);
			if(status < 0)
			{
				printf("Error in BICGStab Solver: %d\n", status);
				return -1;
			}
	*/		
			
			
			gsl_linalg_LU_decomp(dR_dT, p, &s);
			gsl_linalg_LU_solve(dR_dT, p, R, dT);
			
			gsl_vector_add(T_global_N_i, dT);
			//printf("%g\n", gsl_vector_max(dT));
        }
        //printf("%g %g %d %g %g\n", time/60, gsl_vector_get(T_global_N_i, 134)-273, 
		//		itr, maxResidual, gsl_vector_get(CD_global_E_i, 51));
        fprintf(fp1, "%g %g %d %g %g\n", time/60, gsl_vector_get(T_global_N_i, 134)-273					,itr, maxResidual, gsl_vector_get(CD_global_E_i, 51));
		fflush(fp1);
        gsl_vector_memcpy(CD_global_E, CD_global_E_i);
        gsl_vector_memcpy(CDR_global_E, CDR_global_E_i);
        gsl_vector_memcpy(T_global_N, T_global_N_i);

        //printf("%g %g\n", time, gsl_vector_get(T_global_N, 134));
        time += delta_t;
        maxResidual = 1;
		itr = 1;
	}
	gsl_permutation_free(p);
	fclose(fp1);
#endif

	gsl_matrix_free(MK_LHS);
	gsl_matrix_free(MK_RHS);
	gsl_vector_free(F_RHS);
	gsl_matrix_free(M_global_NM);
#endif
/* TRANSIENT SOLVER END */


	gsl_vector_free(T_global_N);
	gsl_vector_free(CD_global_E); gsl_vector_free(CDR_global_E);
	//free(nb_array);
	free(node_arr); free(ele_arr);

	return 0;
}
