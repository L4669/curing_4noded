#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "fem2d.h"



int gsl_linalg_gauss_seidel(gsl_matrix *A,gsl_vector *x,gsl_vector *b, 
		double e)
{
    int i, j, row, col;
    gsl_vector *residual;
	double var, maxE, sum;
    int iterations = 1000, itr = 0, matrixDim = MAX_NODE;

    residual = gsl_vector_alloc(matrixDim);

    while(itr < iterations)
    {
        for(row = 0; row < matrixDim; row++)
        {
            sum = 0;
            for(col = 0; col < row; col++)
				sum += gsl_matrix_get(A, row, col)*gsl_vector_get(x, col);
			
            for(col = row + 1; col < matrixDim; col++)
				sum += gsl_matrix_get(A, row, col)*gsl_vector_get(x, col);
			
			
			var = (gsl_vector_get(b, row) - sum)/gsl_matrix_get(A, row, row);
			gsl_vector_set(residual, row, 
					fabs(gsl_vector_get(x, row)- var)); 
			gsl_vector_set(x, row, var);
        }

		if(((fabs(gsl_vector_max(residual)) < e)) 
				&& (fabs(gsl_vector_min(residual)) < e))
		{
			break;
		}
		itr++;
    }

	gsl_vector_free(residual);

	if(itr < iterations)
	{
		return itr; //residuals are below tolerance supplied
	}
	else
	{
		return -1; //iterations didn't converge
	}

}


int gsl_linalg_BiCGStab(const gsl_matrix *A, gsl_vector *x, const gsl_vector *b, 
		double e)
{
	gsl_vector *r, *r0, *v, *p, *s, *t, *residual;
	int itr = 0, maxItr = 100000, matrixDim = MAX_NODE;
	double rho, rho_old = 1, alpha = 1, omega = 1, beta, temp, resNorm;

	r = gsl_vector_calloc(matrixDim);
	r0 = gsl_vector_calloc(matrixDim);
	v = gsl_vector_calloc(matrixDim);
	p = gsl_vector_calloc(matrixDim);
	s = gsl_vector_calloc(matrixDim);
	t = gsl_vector_calloc(matrixDim);
	residual = gsl_vector_calloc(matrixDim);


	gsl_vector_memcpy(r, b);
	gsl_blas_dgemv(CblasNoTrans, -1.0, A, x, 1.0, r);
	gsl_vector_memcpy(r0, r);

	while(itr < maxItr)
	{
		//gsl_vector_memcpy(residual, x);
		gsl_blas_ddot(r0, r, &rho);
		beta = rho*alpha/(rho_old*omega);
		if(beta == 0)
		{
			return -1;
		}
		gsl_blas_daxpy(-omega, v, p);
		gsl_vector_scale(p, beta);
		gsl_blas_daxpy(1.0, r, p);
		gsl_blas_dgemv(CblasNoTrans, 1.0, A, p, 0, v);
		gsl_blas_ddot(r0, v, &temp);
		if(gsl_isnan(temp) == 1)
		{
			return -2;
		}
		alpha = rho/temp;
		gsl_vector_memcpy(s, r);
		gsl_blas_daxpy(-alpha, v, s);
		gsl_blas_dgemv(CblasNoTrans, 1.0, A, s, 0, t);
		gsl_blas_ddot(t, s, &omega);
		gsl_blas_ddot(t, t, &temp);
		omega = omega/temp;
		if(gsl_isnan(temp) == 1)
		{
			return -2;
		}
		gsl_blas_daxpy(alpha, p, x);
		gsl_blas_daxpy(omega, s, x);
		gsl_vector_memcpy(r, s);
		gsl_blas_daxpy(-omega, t, r);
		
		//gsl_vector_sub(residual, x);
		gsl_vector_memcpy(residual, b);
		gsl_blas_dgemv(CblasNoTrans, -1.0, A, x, 1.0, residual);
		//printf("%g %g\n", fabs(gsl_vector_max(residual)), fabs(gsl_vector_min(residual)));
		resNorm = gsl_blas_dnrm2(residual);
		if(resNorm < e)
		{
			break;
		}
		rho_old = rho;

		itr++;
	}

	gsl_vector_free(r);
	gsl_vector_free(r0);
	gsl_vector_free(v);
	gsl_vector_free(p);
	gsl_vector_free(s);
	gsl_vector_free(t);
	gsl_vector_free(residual);

	if(itr < maxItr)
	{
		return itr;
	}
	else
	{
		return -1;
	}

}

