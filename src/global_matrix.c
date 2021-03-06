#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "fem2d.h"
#include "config.h"

void global_stiffness_matrix(struct node_list *node_arr, struct ele_list *ele_arr,
		gsl_matrix *K_global_NM)
{
	int ele, node_local, node_global, node_global_adj, eq;
	double *x, *y, temp;
	gsl_matrix *K_local_NM;

	x = (double *)malloc(sizeof(double)*4);
    y = (double *)malloc(sizeof(double)*4);
	K_local_NM = gsl_matrix_calloc(4, 4);

	for(ele = 0; ele < MAX_ELEMENT; ele++)
    {
        for(node_local = 0; node_local < 4; node_local++)
        {

            node_global = ele_arr[ele].ele_node_no[node_local] - 1;
            x[node_local] = node_arr[node_global].x;
            y[node_local] = node_arr[node_global].y;
        }

        stiffness_matrix_bqe(x, y, K_local_NM);
		//included K in the local matrix formation
		//due to directionl nature of K
        //gsl_matrix_scale(K_local_NM, K);

        for(eq = 0; eq < 4; eq++)
        {
            node_global = ele_arr[ele].ele_node_no[eq] - 1;
            for(node_local = 0; node_local < 4; node_local++)
            {
                node_global_adj = ele_arr[ele].ele_node_no[node_local]
                    - 1;

                temp = gsl_matrix_get(K_global_NM, node_global,
                        node_global_adj)
                    + gsl_matrix_get(K_local_NM, eq, node_local);

                gsl_matrix_set(K_global_NM, node_global,
                        node_global_adj, temp);
            }
        }
    }

	free(x); free(y);
	gsl_matrix_free(K_local_NM);

}

void global_neumann_bc(struct neumann_boundary_edgelist *nb_array, 
		gsl_matrix *K_global_NM, gsl_vector *F_global_N, struct node_list *node_arr,
		struct ele_list *ele_arr, double T_inf, double h, int neu_edges)
{
	int n, m, ele, node_local, node_global, node_global_adj, eq;
	double *x, *y ,temp;
	gsl_matrix *Kb_local_NM;
	gsl_vector *F_local_N;
	
	x = (double *)malloc(sizeof(double)*4);
    y = (double *)malloc(sizeof(double)*4);
	Kb_local_NM = gsl_matrix_calloc(4, 4);
	F_local_N = gsl_vector_calloc(4);
	
    //type1 and type2   matrix formation    
    for(n = 0; n < neu_edges; n++)
    {
        ele = nb_array[n].element - 1;
        for(m = 0; m < 4; m++)
        {
            node_global = ele_arr[ele].ele_node_no[m] - 1;
            x[m] = node_arr[node_global].x;
            y[m] = node_arr[node_global].y;
        }
        boundary_matrix_bqe_type2(nb_array[n].edge_local, x, y, Kb_local_NM);
        boundary_matrix_bqe_type1(nb_array[n].edge_local, x, y, F_local_N);

        gsl_matrix_scale(Kb_local_NM, h);
        gsl_vector_scale(F_local_N, h*T_inf);

        for(eq = 0; eq < 4; eq++)
        {
            node_global = ele_arr[ele].ele_node_no[eq] - 1;
            for(node_local = 0; node_local < 4; node_local++)
            {
                node_global_adj = ele_arr[ele].ele_node_no[node_local] - 1;
                temp = gsl_matrix_get(K_global_NM, node_global, node_global_adj)
                    + gsl_matrix_get(Kb_local_NM, eq, node_local);
                gsl_matrix_set(K_global_NM, node_global, node_global_adj, temp);
            }

            temp = gsl_vector_get(F_global_N, node_global) +
                gsl_vector_get(F_local_N, eq);
            gsl_vector_set(F_global_N, node_global, temp);
        }
    }

	free(x); free(y);
	gsl_matrix_free(Kb_local_NM);
	gsl_vector_free(F_local_N);

}

/*
void global_neumann_bc_cf(struct neumann_boundary_edgelist *nb_array, 
		gsl_vector *F_global_N, struct node_list *node_arr, struct ele_list *ele_arr,
		int neu_edges_cf)
{
	int n, m, ele, node_local, node_global, node_global_adj, eq;
	double *x, *y ,temp;
	gsl_vector *F_local_N;
	
	x = (double *)malloc(sizeof(double)*4);
    y = (double *)malloc(sizeof(double)*4);
	F_local_N = gsl_vector_calloc(4);
	
    //type1 matrix formation    
    for(n = 0; n < neu_edges_cf; n++)
    {
        ele = nb_array[n].element - 1;
        for(m = 0; m < 4; m++)
        {
            node_global = ele_arr[ele].ele_node_no[m] - 1;
            x[m] = node_arr[node_global].x;
            y[m] = node_arr[node_global].y;
        }

        boundary_matrix_bqe_type1(nb_array[n].edge_local, x, y, F_local_N);

        gsl_vector_scale(F_local_N, -q_cf);

        for(eq = 0; eq < 4; eq++)
        {
            node_global = ele_arr[ele].ele_node_no[eq] - 1;

            temp = gsl_vector_get(F_global_N, node_global) +
                gsl_vector_get(F_local_N, eq);
            gsl_vector_set(F_global_N, node_global, temp);
        }
    }

	free(x); free(y);
	gsl_vector_free(F_local_N);

}
*/


void global_dirichlet_bc(struct dirichlet_boundary_nodelist *db_array, 
		gsl_matrix *K_global_NM, gsl_vector *F_global_N)
{
	int n, m, node_global;

	for(n = 0; n < dirich_nodes; n++)
    {
        node_global = db_array[n].node - 1;
        gsl_matrix_set(K_global_NM, node_global, node_global, 1);
        gsl_vector_set(F_global_N, node_global, db_array[n].value);
        for(m = 0; m < MAX_NODE; m++)
        {
            if(m == node_global)
            {
                continue;
            }
            else
            {
                gsl_matrix_set(K_global_NM, node_global, m, 0);
            }
        }
    }

}


void global_mass_matrix(struct node_list *node_arr, struct ele_list *ele_arr,
		gsl_matrix *M_global_NM)
{
	int ele, node_local, node_global, node_global_adj, eq;
	double *x, *y, temp;
	gsl_matrix *M_local_NM;

	x = (double *)malloc(sizeof(double)*4);
    y = (double *)malloc(sizeof(double)*4);
	M_local_NM = gsl_matrix_calloc(4, 4);
	
	for(ele = 0; ele < MAX_ELEMENT; ele++)
    {
        for(node_local = 0; node_local < 4; node_local++)
        {

            node_global = ele_arr[ele].ele_node_no[node_local] - 1;
            x[node_local] = node_arr[node_global].x;
            y[node_local] = node_arr[node_global].y;
        }

        mass_matrix_bqe(x, y, M_local_NM);

        for(eq = 0; eq < 4; eq++)
        {
            node_global = ele_arr[ele].ele_node_no[eq] - 1;
            for(node_local = 0; node_local < 4; node_local++)
            {
                node_global_adj = ele_arr[ele].ele_node_no[node_local]
                    - 1;

                temp = gsl_matrix_get(M_global_NM, node_global,
                        node_global_adj)
                    + gsl_matrix_get(M_local_NM, eq, node_local);

                gsl_matrix_set(M_global_NM, node_global,
                        node_global_adj, temp);
            }
        }
    }
	free(x); free(y); 
	gsl_matrix_free(M_local_NM);

}
/*
void global_heat_source_matrix(struct node_list *node_arr, struct ele_list *ele_arr,
		gsl_vector *F_global_N, gsl_vector *CDR_global_E, gsl_vector *d_CDR,
		gsl_vector *dF)
{
	int ele, node_local, node_global, node_global_adj, eq;
	double *x, *y ,temp, Q;
	gsl_vector *F_local_N, *F_local_N_1;
	
	x = (double *)malloc(sizeof(double)*4);
    y = (double *)malloc(sizeof(double)*4);
	F_local_N = gsl_vector_calloc(4);
	F_local_N_1 = gsl_vector_calloc(4);
	
    for(ele = 0; ele < MAX_ELEMENT; ele++)
    {
        for(node_local = 0; node_local < 4; node_local++)
        {
            node_global = ele_arr[ele].ele_node_no[node_local] - 1;
            x[node_local] = node_arr[node_global].x;
            y[node_local] = node_arr[node_global].y;
        }

        heat_source_matrix_bqe(x, y, F_local_N);
		gsl_vector_memcpy(F_local_N_1, F_local_N);
		Q = gsl_vector_get(CDR_global_E, ele)*rho*H_r;
        gsl_vector_scale(F_local_N, Q);

		Q = gsl_vector_get(d_CDR, ele)*rho*H_r;
        gsl_vector_scale(F_local_N_1, Q);

        for(eq = 0; eq < 4; eq++)
        {
            node_global = ele_arr[ele].ele_node_no[eq] - 1;

            temp = gsl_vector_get(F_global_N, node_global) +
                gsl_vector_get(F_local_N, eq);
            gsl_vector_set(F_global_N, node_global, temp);

            temp = gsl_vector_get(dF, node_global) + gsl_vector_get(F_local_N_1, eq);
            gsl_vector_set(dF, node_global, temp);
        }
    }

	free(x); free(y);
	gsl_vector_free(F_local_N);

}


void global_cure_rate_vector(gsl_vector *T_global_N, gsl_vector *CD_global_E,
		gsl_vector *CD_global_E_i, gsl_vector *CDR_global_E_i, 
		struct ele_list *ele_arr, gsl_vector *d_CDR, double time)
{
	int node_local = 0, n = 0, node_global = 0;
	double alpha, dalpha_dt, alpha_old, T_avg = 0, b;

	for(n = 0; n < MAX_ELEMENT; n++)
	{
		T_avg = 0;
        for(node_local = 0; node_local < 4; node_local++)
        {
            node_global = ele_arr[n].ele_node_no[node_local] - 1;
			T_avg += gsl_vector_get(T_global_N, node_global);
        }
		T_avg = T_avg/4;

		b = A_c*exp(-E_c/(GSL_CONST_MKSA_MOLAR_GAS*(0+T_avg)));

		alpha_old = gsl_vector_get(CD_global_E, n);
		//alpha_old = 1e-15;
		alpha = degree_of_cure_2(alpha_old, time, b);
		
		dalpha_dt = b*pow(alpha, m_c)*pow(1-alpha, n_c);

		gsl_vector_set(d_CDR, n, 0.25*dalpha_dt*E_c/((0+T_avg)*(0+T_avg)*
						GSL_CONST_MKSA_MOLAR_GAS));

		gsl_vector_set(CD_global_E_i, n, alpha);
		gsl_vector_set(CDR_global_E_i, n, dalpha_dt);
	}

}
void global_cure_rate_vector_2(gsl_vector *T_global_N, gsl_vector *CD_global_E,
		gsl_vector *CD_global_E_i, gsl_vector *CDR_global_E_i, 
		struct ele_list *ele_arr, gsl_vector *d_CDR, double time)
{
	int node_local = 0, n = 0, node_global = 0;
	double alpha, dalpha_dt, alpha_old, T_avg = 0, b;

	for(n = 0; n < MAX_ELEMENT; n++)
	{
		T_avg = 0;
        for(node_local = 0; node_local < 4; node_local++)
        {
            node_global = ele_arr[n].ele_node_no[node_local] - 1;
			T_avg += gsl_vector_get(T_global_N, node_global);
        }
		T_avg = T_avg/4;

		b = A_c*exp(-E_c/(GSL_CONST_MKSA_MOLAR_GAS*(0+T_avg)));

		alpha = gsl_vector_get(CD_global_E_i, n);
		
		dalpha_dt = b*pow(alpha, m_c)*pow(1-alpha, n_c);

		gsl_vector_set(d_CDR, n, 0.25*dalpha_dt*E_c/((0+T_avg)*(0+T_avg)*
						GSL_CONST_MKSA_MOLAR_GAS));

		//gsl_vector_set(CD_global_E_i, n, alpha);
		gsl_vector_set(CDR_global_E_i, n, dalpha_dt);
	}

}
*/
void residual_derivative_matrix(gsl_matrix *dR_dT, gsl_matrix *MK_LHS, gsl_vector *dF,
		struct ele_list *ele_arr)
{
	int n, m, p, node_local, node_global;
	double temp;
	gsl_vector *nodes;
	nodes = gsl_vector_alloc(4);

	gsl_matrix_memcpy(dR_dT, MK_LHS);
	
	for(n = 0; n < MAX_NODE; n++)
	{
		temp = gsl_matrix_get(dR_dT, n, n) + gsl_vector_get(dF, n);
		//temp = gsl_matrix_get(dR_dT, n, n);
		gsl_matrix_set(dR_dT, n, n, temp);
	}
	

/*	
	for(n = 0; n < MAX_ELEMENT; n++)
	{
		for(node_local = 0; node_local < 4; node_local++)
		{
			node_global = ele_arr[n].ele_node_no[node_local]-1;
			gsl_vector_set(nodes, node_local, node_global);
		}

		for(p = 0; p < 4; p++)
		{
			for(m = 0; m < 4; m++)
			{
				temp = gsl_vector_get(dF, gsl_vector_get(nodes, p)) + 
					gsl_matrix_get(dR_dT, gsl_vector_get(nodes, p), 
							gsl_vector_get(nodes,m));
				gsl_matrix_set(dR_dT, gsl_vector_get(nodes, p), 
							gsl_vector_get(nodes,m), temp);

			}
		}
	}
*/
	gsl_vector_free(nodes);
}


void global_cure_rate_vector(gsl_vector *T_global_N, gsl_vector *CD_global_E,
		gsl_vector *CD_global_E_i,struct ele_list *ele_arr, double time)
{
	int node_local = 0, n = 0, node_global = 0;
	double alpha, dalpha_dt, alpha_old, T_avg = 0, b;

	for(n = 0; n < MAX_NODE; n++)
	{
		
		T_avg = gsl_vector_get(T_global_N, n);

		//b = A_c*exp(-E_c/(GSL_CONST_MKSA_MOLAR_GAS*(0+T_avg)));
		
		alpha_old = gsl_vector_get(CD_global_E, n);
		alpha = degree_of_cure_3(alpha_old, time, T_avg);
		
		//alpha = alpha_old + b*pow(alpha_old, m_c)*pow(1-alpha_old, n_c)*delta_t;
		//dalpha_dt = b*pow(alpha, m_c)*pow(1-alpha, n_c);

		//gsl_vector_set(d_CDR, n, 0.25*dalpha_dt*E_c/((0+T_avg)*(0+T_avg)*
		//				GSL_CONST_MKSA_MOLAR_GAS));

		gsl_vector_set(CD_global_E_i, n, alpha);
		//gsl_vector_set(CDR_global_E_i, n, dalpha_dt);
	}

}

void global_cure_rate_vector_2(gsl_vector *T_global_N, gsl_vector *CD_global_E,
		struct ele_list *ele_arr, gsl_vector *CDR_global_E_i, double time)
{
	int node_local = 0, n = 0, node_global = 0;
	long double alpha, dalpha_dt, alpha_old, T_avg = 0, b;

	for(n = 0; n < MAX_NODE; n++)
	{
		
		T_avg = gsl_vector_get(T_global_N, n);

		b = A_c*exp(-E_c/(GSL_CONST_MKSA_MOLAR_GAS*(0+T_avg)));
		
		alpha = gsl_vector_get(CD_global_E, n);
		
		dalpha_dt = b*pow(alpha, m_c)*pow(1-alpha, n_c);

		gsl_vector_set(CDR_global_E_i, n, dalpha_dt);
	}

}

void global_heat_source_matrix(struct node_list *node_arr, struct ele_list *ele_arr,
		gsl_vector *F_global_N, gsl_vector *CD_global_E, gsl_vector *d_CDR,
		gsl_vector *dF, gsl_vector *T_global_N_i, gsl_vector *T_global_N)
{
	int ele, node_local, node_global, node_global_adj, eq;
	double *x, *y, *gauss_temp, *alpha_gauss, *dt_dT;
	double temp;
	gsl_vector *F_local_N, *F_local_N_1, *Q;
	
	x = (double *)malloc(sizeof(double)*4);
    y = (double *)malloc(sizeof(double)*4);
    gauss_temp = (double *)malloc(sizeof(double)*4);
    alpha_gauss = (double *)malloc(sizeof(double)*4);
    dt_dT = (double *)malloc(sizeof(double)*4);
    
	//Q = gsl_vector_calloc(4);
	F_local_N = gsl_vector_calloc(4);
	F_local_N_1 = gsl_vector_calloc(4);
	
    for(ele = 0; ele < MAX_ELEMENT; ele++)
    {
        for(node_local = 0; node_local < 4; node_local++)
        {
            node_global = ele_arr[ele].ele_node_no[node_local] - 1;
            x[node_local] = node_arr[node_global].x;
            y[node_local] = node_arr[node_global].y;
			gauss_temp[node_local] = gsl_vector_get(T_global_N_i, node_global);
			dt_dT[node_local] = delta_t/(gsl_vector_get(T_global_N_i, node_global) - 
				0.99*gsl_vector_get(T_global_N, node_global));
			alpha_gauss[node_local] = gsl_vector_get(CD_global_E, node_global);
        }

        heat_source_matrix_bqe_1(x, y, F_local_N, gauss_temp, alpha_gauss);
        heat_source_matrix_bqe_2(x, y, F_local_N_1, gauss_temp, alpha_gauss, dt_dT);

		/*
		for(eq = 0; eq < 4; eq++)
		{
			node_global = ele_arr[ele].ele_node_no[eq] - 1;
			temp = gsl_vector_get(CDR_global_E, node_global)*rho*H_r;
			gsl_vector_set(Q, eq, temp);
		}
		*/
		gsl_vector_scale(F_local_N, rho*H_r);
		
		/*
		for(eq = 0; eq < 4; eq++)
		{
			node_global = ele_arr[ele].ele_node_no[eq] - 1;
			temp = gsl_vector_get(d_CDR, node_global)*rho*H_r;
			gsl_vector_set(Q, eq, temp);
		}
		*/
        gsl_vector_scale(F_local_N_1, rho*H_r);


        for(eq = 0; eq < 4; eq++)
        {
            node_global = ele_arr[ele].ele_node_no[eq] - 1;

            temp = gsl_vector_get(F_global_N, node_global) +
                gsl_vector_get(F_local_N, eq);
            gsl_vector_set(F_global_N, node_global, temp);
            
			temp = gsl_vector_get(dF, node_global) + gsl_vector_get(F_local_N_1, eq);
            gsl_vector_set(dF, node_global, temp);

        }
    }

	free(x); free(y); free(dt_dT);
	free(alpha_gauss); free(gauss_temp);
	gsl_vector_free(F_local_N);
	gsl_vector_free(F_local_N_1);

}

void global_neumann_bc_2(struct neumann_boundary_edgelist *nb_array, 
		gsl_vector *dF, struct node_list *node_arr, struct ele_list *ele_arr, 
		double h, int neu_edges, double dT_inf, gsl_vector *T_global_N_i, 
		gsl_vector *T_global_N)
{
	int n, m, ele, node_local, node_global, node_global_adj, eq;
	double *x, *y ,temp;
	gsl_vector *F_local_N, *dt_dT;
	
	x = (double *)malloc(sizeof(double)*4);
    y = (double *)malloc(sizeof(double)*4);
	F_local_N = gsl_vector_calloc(4);
	dt_dT = gsl_vector_calloc(4);
	
    //type1 and type2   matrix formation    
    for(n = 0; n < neu_edges; n++)
    {
        ele = nb_array[n].element - 1;
        for(m = 0; m < 4; m++)
        {
            node_global = ele_arr[ele].ele_node_no[m] - 1;
            x[m] = node_arr[node_global].x;
            y[m] = node_arr[node_global].y;
			temp = delta_t/(gsl_vector_get(T_global_N_i, node_global) - 
				0.99*gsl_vector_get(T_global_N, node_global));
			gsl_vector_set(dt_dT, node_local, temp);

        }
        boundary_matrix_bqe_type1(nb_array[n].edge_local, x, y, F_local_N);

        gsl_vector_scale(F_local_N, h*dT_inf);
		gsl_vector_mul(F_local_N, dt_dT);

        for(eq = 0; eq < 4; eq++)
        {
            node_global = ele_arr[ele].ele_node_no[eq] - 1;
            temp = gsl_vector_get(dF, node_global) +
                gsl_vector_get(F_local_N, eq);
            gsl_vector_set(dF, node_global, temp);
        }
    }

	free(x); free(y);
	gsl_vector_free(dt_dT);
	gsl_vector_free(F_local_N);

}

