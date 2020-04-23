#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>

/*
    //check for the positive definitness of the stiffness matrix
    gsl_error_handler_t *default_handler;
    default_handler = gsl_set_error_handler_off();

    if(gsl_linalg_cholesky_decomp(K_global_NM) == GSL_EDOM)
    {
        printf("Matrix is not positive definite\n");
    }
    else
    {
        printf("Matrix is positive definite\n");
    }
    */

    /*
    //check for symmetricity
    gsl_matrix *transpose_matrix;
    transpose_matrix = gsl_matrix_alloc(MAX_NODE, MAX_NODE);
    gsl_matrix_transpose_memcpy(transpose_matrix, K_global_NM);

    if(gsl_matrix_equal(transpose_matrix, K_global_NM) == 1)
    {
        printf("matrix is symmetric\n");
    }
    else
    {
        printf("matrix is not symmetric\n");
    }
    */

