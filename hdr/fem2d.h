#ifndef FEM2D_H  //prevents the inclusion of the header file more than once in a given compilation
#define FEM2D_H 1

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_const_mksa.h>
/*
// Material Properties
extern const double rho_r;  // Resin density (kg/m^3)
extern const double cp_r;   // Specific heat of resin (J/kg/K)
extern const double k_r;    // Thermal conductivity of resin (W/m/K)
extern const double rho_f;  // Fibre density (kg/m^3)
extern const double cp_f;   // Specific heat of fibre (J/kg/K)
extern const double k_f;    // Thermal conductivity of fibre (W/m/K)
extern const double K_xx, K_zz; // thermal conductivities of composite in x and y direction

// Chemical Properties
extern const double A1;     // Pre-exponential factor (min^-1)
extern const double A2;     // Pre-exponential factor (min^-1)
extern const double A3;     // Pre-exponential factor (min^-1)
extern const double del_E1; // Activation energy (J/mol)
extern const double del_E2; // Activation energy (J/mol)
extern const double del_E3; // Activation energy (J/mol)
extern const double Hr;     // Heat of reaction (J/kg)
extern const double U;      // Activation energy for viscosity (J/mol)
extern const double mu_inff;// Viscosity constant (Pa-s)
extern const double d;      // Viscosity constant
extern const double r_f;    // Fibre radius (m)
extern const double Kz;     // Kozeny constant
extern const double V0;     // Initial fibre volume fraction
extern const double b_bar;     //constant in rate of degree of cure

//properties for linear benchmark problems 
extern const double K;
extern double T_inf;
extern double h;
extern const double rho;
extern const double cp;
extern const double q_cf;
extern const double Q;
*/

//Glass-Polyster [ref. 2D Cure Simulation of Thick Thermosetting Composites]
extern const double rho;
extern const double cp;
extern const double m_c;
extern const double n_c;
extern const double A_c;
extern const double E_c;
extern const double H_r;
extern double K_xx;
extern double K_zz;


//mesh file parameters
extern const int MAX_ELEMENT, MAX_NODE, dirich_nodes, neu_edges, neu_edges_cf;

//transient solver parameters
extern const double theta, maxTime;
extern double delta_t; 
//structure declaration to handle different parameters
struct node_list
{
    int node;
    double x;
    double y;
    double z;
};

struct ele_list
{
    int ele_no;
	int ele_node_no[4];
};

struct dirichlet_boundary_nodelist
{
	int node;
	double value;
};

struct neumann_boundary_edgelist
{
	int edge_local;
	int element;
	double length;
	double value;
	double dcosine[3];
};

/*function prototypes */

//def in nodes_connectivity.c
void node_connectivity(struct node_list *, struct ele_list *, char *);
void dirichlet_boundary_nodes(struct dirichlet_boundary_nodelist *db_array,
		        char *fname, char *tag);
void neumann_boundary_edges(struct neumann_boundary_edgelist *nb_array,
		        struct ele_list *ele_arr, struct node_list *node_arr, 
				char *fname, char *tag, int neu_edges);

//def in gauss_seidel.c
int gsl_linalg_gauss_seidel(gsl_matrix *A,gsl_vector *x,gsl_vector *b, double e);
int gsl_linalg_BiCGStab(const gsl_matrix *A, gsl_vector *x, const gsl_vector *b, 
		double e);

//def in material_properties.c
double rho_c(double Vf);
double k_x(double Vf);
double k_y(double Vf);
double k_z(double Vf, int model);
int F(double t, const double alpha[], double dalpha_dt[], void *params);
double degree_of_cure(double alpha_old, double t_final, double b,
        int NODES_PER_ELEMENT);
double degree_of_cure_2(double alpha_old, double t_final, double b);
double degree_of_cure_3(double alpha_old, double t_final, double b);


//def in local_elements.c
double stiffness_matrix_nm_bqe(double zeta, double ita, void *data);
void stiffness_matrix_bqe(double *cord_x, double *cord_y, 
		gsl_matrix *K_local_NM);

double mass_matrix_nm_bqe(double zeta, double ita, void *data);
void mass_matrix_bqe(double *cord_x, double *cord_y,
		gsl_matrix *M_local_NM);

double boundary_matrix_n_bqe_type1(double xi, void *data);
void boundary_matrix_bqe_type1(int boundary, double *x_cord, double *y_cord,
		gsl_vector *F_local_N);

double boundary_matrix_n_bqe_type2(double zeta, void *data);
double boundary_matrix_m_bqe_type2(double zeta, void *data);
void boundary_matrix_bqe_type2(int boundary, double *x_cord, double *y_cord,
		        gsl_matrix *Kb_local_NM);

double heat_source_matrix_n_bqe_1(double zeta, double ita, void *data);
void heat_source_matrix_bqe_1(double *cord_x, double *cord_y,
		        gsl_vector *F_local_N, double *gauss_temp,
				double *alpha_gauss);
double heat_source_matrix_n_bqe_2(double zeta, double ita, void *data);
void heat_source_matrix_bqe_2(double *cord_x, double *cord_y,
		        gsl_vector *F_local_N, double *gauss_temp,
				double *alpha_gauss, double *dt_dT);

//def in vector_products.c
void vector_cross_product(const double *u, const double *v, double *product);
double vector_dot_product(const double *u, const double *v);
double vector_mod(double *vector);

//def in global_matrix.c
void global_stiffness_matrix(struct node_list *node_arr, struct ele_list *ele_arr,
		gsl_matrix *K_global_NM);
void global_neumann_bc(struct neumann_boundary_edgelist *nb_array,
        gsl_matrix *K_global_NM, gsl_vector *F_global_N, struct node_list *node_arr,
        struct ele_list *ele_arr, double T_inf, double h, int neu_edges);
void global_dirichlet_bc(struct dirichlet_boundary_nodelist *db_array, 
        gsl_matrix *K_global_NM, gsl_vector *F_global_N);
void global_mass_matrix(struct node_list *node_arr, struct ele_list *ele_arr,
        gsl_matrix *M_global_NM);
void global_neumann_bc_cf(struct neumann_boundary_edgelist *nb_array,
        gsl_vector *F_global_N, struct node_list *node_arr, struct ele_list *ele_arr,
		int neu_edges_cf);
void global_heat_source_matrix(struct node_list *node_arr, struct ele_list *ele_arr,
		        gsl_vector *F_global_N, gsl_vector *CDR_global_E,
				gsl_vector *d_CDR, gsl_vector *dF, gsl_vector *T_global_N_i,
				gsl_vector *T_global_N);

void global_cure_rate_vector(gsl_vector *T_global_N, gsl_vector *CD_global_E,
		        gsl_vector *CD_global_E_i,struct ele_list *ele_arr, double time);
void residual_derivative_matrix(gsl_matrix *dR_dT, gsl_matrix *MK_LHS, gsl_vector *dF,
		struct ele_list *ele_arr);
void global_cure_rate_vector_2(gsl_vector *T_global_N, gsl_vector *CD_global_E,
		        struct ele_list *ele_arr, gsl_vector *d_CDR, double time);

void global_neumann_bc_2(struct neumann_boundary_edgelist *nb_array,
        gsl_vector *dF, struct node_list *node_arr, struct ele_list *ele_arr, 
        double h, int neu_edges, double dT_inf, gsl_vector *T_global_N_i, 
        gsl_vector *T_global_N);

//def in curing_cycle.c
double autoclave_input_cycle(double time);
double autoclave_input_cycle_2(double time);

#endif
