#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include "fem2d.h"


/*
// Hercules AS4/3501-6

// Material Properties
const double rho_r = 1.26e3;
const double cp_r = 1.26e3;
const double k_r = 0.167;
const double rho_f = 1.79e3;
const double cp_f = 7.12e2;
const double k_f = 26;

// Chemical Properties
const double A1 = 2.102e9/60;
const double A2 = -2.014e9/60;
const double A3 = 1.96e5/60;
const double del_E1 = 8.07e4;
const double del_E2 = 7.78e4;
const double del_E3 = 5.66e4;
const double Hr = 198.9e3;
const double U = 9.08e4;
const double mu_inff = 7.93e-4;
const double d = 14.1;
const double r_f = 3.5e-6;
const double Kz = 6;
const double V0 = 0.53;
const double b_bar = 0.47/60;
*/

// Glass-Polyster [ref. 2D Cure Simulation of Thick Thermosetting Composites]
// All in SI units
const double rho = 1.89e3;
const double cp = 1.26e3;
const double m_c = 0.524;
const double n_c = 1.476;
const double A_c = 3.7e22/60;
const double E_c = 1.674e5;
const double H_r = 77.5e3; 
double K_xx = 2.163e-1, K_zz = 2*2.163e-1;
/*
double rho_c(double Vf)
{
    double rho_composite = rho_f*Vf + rho_r*(1-Vf);
    return rho_composite;
}

double k_x(double Vf)
{
    double k_x = k_f*Vf + k_r*(1-Vf);
    return k_x;
}


double k_y(double Vf)
{
    double k_y = k_f*Vf + k_r*(1-Vf);
    return k_y;
}

double k_z(double Vf, int model)
{
	double k_z;
    if(model == 1)
    {
        // Springer-Tsai Model
        double D,temp1,temp2,temp3;
        D = 2*(k_r/k_f-1);
        temp1 = (sqrt(1-(D*D*Vf/M_PI)))/(1+D*sqrt(Vf/M_PI));
        temp2 = 4/sqrt(1-D*D*Vf/M_PI)*atan(temp1);
        temp3 = (1-2*sqrt(Vf/M_PI));
        k_z = k_r*(temp3+(M_PI-temp2)/D);
        return k_z;
    }
    else if(model == 2)
    {
        // Tasi-Halpin Model
        k_z = k_r*((k_f+k_r+(k_f-k_r)*Vf)/(k_f+k_r-(k_f-k_r)*Vf));
        return k_z;
    }
    else
        return -1;

}
*/
int F(double t, const double alpha[], double dalpha_dt[], void *params)
{
	double *coeff = (double *)params;
	dalpha_dt[0] = coeff[0]*pow(alpha[0], m_c)*pow(1 - alpha[0], n_c);

	return GSL_SUCCESS;
}

double degree_of_cure(double alpha_old, double t_final, double b, 
		int NODES_PER_ELEMENT)
{
	double tau = t_final/(1000.0), t = 0;
	
	double params[1] = {b};
	gsl_odeiv2_system system = {F, NULL, 1, params};
	
	const gsl_odeiv2_step_type *st = gsl_odeiv2_step_rk4;


	gsl_odeiv2_step *stepping = gsl_odeiv2_step_alloc(st, 1);

	double alpha[1] = {alpha_old};
	double err[1] = {0};

	while(t <= t_final)
	{
		gsl_odeiv2_step_apply(stepping, t, tau, alpha, err, 
				NULL, NULL, &system);
		t += tau;
	}

	gsl_odeiv2_step_free(stepping);

	return alpha[0];
}


double degree_of_cure_2(double alpha_old, double t_final, double b)
{
	int i;
	double t = t_final, ti;
	double params[1] = {b};
	gsl_odeiv2_system system = {F, NULL, 1, params};

	gsl_odeiv2_driver *d = 
		gsl_odeiv2_driver_alloc_y_new(&system, gsl_odeiv2_step_rk4,
				1e-15, 1e-15, 0.0);

	double alpha[1] = {alpha_old};
	
	for(i = 0; i <= 200; i++ )
	{
		ti = t_final + i*delta_t/200.0;
		int status = gsl_odeiv2_driver_apply(d, &t, ti, alpha);
		if( status != GSL_SUCCESS)
		{
			printf("Error in gsl_driver_apply, return value= %d\n", status);
			return -1;
		}
	}

	gsl_odeiv2_driver_free(d);

	return alpha[0];
}

double f(double alpha, double b)
{
	double alpha_h = b*pow(alpha, m_c)*pow(1 - alpha, n_c);
	return alpha_h;
}

double degree_of_cure_3(double alpha_old, double time, double T_old)
{
	double t = time, k1, k2, k3, k4, h, b, dalpha_dt; 
	double alpha = alpha_old, T_M = T_old, T = T_old;
	int m, M = 500;

	h = delta_t/M;
	
	for(m = 1; m <= M; m++)
	{
		dalpha_dt = A_c*exp(-E_c/(GSL_CONST_MKSA_MOLAR_GAS*(0+T_M)))*
			pow(alpha,m_c)*pow(1-alpha, n_c);
		alpha += dalpha_dt*h;
		T_M +=  (H_r/cp)*dalpha_dt*h;
	}
	
	alpha = alpha_old; 
	while(t <= time+delta_t)
	{
		T += (T_M - T_old)*(t- time)/delta_t;
		b = A_c*exp(-E_c/(GSL_CONST_MKSA_MOLAR_GAS*T));
		k1 = f(alpha, b);
		k2 = f(alpha + 0.5*h*k1, b);
		k3 = f(alpha + 0.5*h*k2, b);
		k4 = f(alpha + h*k3, b);

		alpha = alpha + h*(k1 + 2*k2 + 2*k3 + k4)/6;
		t += h;
	}
	return alpha;
}
