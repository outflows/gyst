#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "definitions.h"
#include "declarations.h"
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>


void get_1stderivative(int i, int j, int k, double ***J_cs)
{
	double dJdx_ip, dJdx_im, dJdy_jp, dJdy_jm, dJdz_kp, dJdz_km;

	//1st derivative -- must be outside of this function
	if (i+1 > N1-1) dJdx_ip = J_cs[i][j][k];
	else dJdx_ip = J_cs[i+1][j][k];
	if (i-1 < 0) dJdx_im = J_cs[i][j][k];
	else dJdx_im = J_cs[i-1][j][k];
	dJdx[i][j][k] = (dJdx_ip - dJdx_im)/2.;

	if (j+1 > N2-1) dJdy_jp = J_cs[i][j][k];
	else dJdy_jp = J_cs[i][j+1][k];
	if (j-1 < 0) dJdy_jm = J_cs[i][j][k];
	else dJdy_jm = J_cs[i][j-1][k];
	dJdy[i][j][k] = (dJdy_jp - dJdy_jm)/2.;

	if (k+1 > N3-1) dJdz_kp = J_cs[i][j][k];
	else dJdz_kp = J_cs[i][j][k+1];
	if (k-1 < 0) dJdz_km = J_cs[i][j][k];
	else dJdz_km = J_cs[i][j][k-1];
	dJdz[i][j][k] = (dJdz_kp - dJdz_km)/2.;
}

void get_hessian(int i, int j, int k, double **Hess3D)
{
    // IMPORTANT WHEN CALLING THIS FUNCTION! (not anymore, we changed to void)
    // https://www.tutorialspoint.com/cprogramming/c_return_arrays_from_function.htm

	int m, n;
	double d2Jdx2_ip, d2Jdx2_im, d2Jdy2_jp, d2Jdy2_jm, d2Jdz2_kp, d2Jdz2_km;
	double d2Jdxdy_jp, d2Jdxdy_jm, d2Jdydx_ip, d2Jdydx_im;
	double d2Jdxdz_kp, d2Jdxdz_km, d2Jdydz_kp, d2Jdydz_km;
	double d2Jdzdx_ip, d2Jdzdx_im, d2Jdzdy_jp, d2Jdzdy_jm;
	double d2Jdx2, d2Jdxdy, d2Jdxdz;
	double d2Jdydx, d2Jdy2, d2Jdydz;
	double d2Jdzdx, d2Jdzdy, d2Jdz2;

	for (m = 0; m < 3; m++) {
		for (n = 0; n < 3; n++) {
			Hess3D[m][n] = 0.;
		}
	}

	//2nd derivative (Hessian)
	if (i+1 > N1-1) d2Jdx2_ip = dJdx[i][j][k];
	else d2Jdx2_ip = dJdx[i+1][j][k];
	if (i-1 < 0) d2Jdx2_im = dJdx[i][j][k];
	else d2Jdx2_im = dJdx[i-1][j][k];
	d2Jdx2 = (d2Jdx2_ip - d2Jdx2_im)/2.;

	if (j+1 > N2-1) d2Jdy2_jp = dJdy[i][j][k];
	else d2Jdy2_jp = dJdy[i][j+1][k];
	if (j-1 < 0) d2Jdy2_jm = dJdy[i][j][k];
	else d2Jdy2_jm = dJdy[i][j-1][k];
	d2Jdy2 = (d2Jdy2_jp - d2Jdy2_jm)/2.;

	if (k+1 > N3-1) d2Jdz2_kp = dJdz[i][j][k];
	else d2Jdz2_kp = dJdz[i][j][k+1];
	if (k-1 < 0) d2Jdz2_km = dJdz[i][j][k];
	else d2Jdz2_km = dJdz[i][j][k-1];
	d2Jdz2 = (d2Jdz2_kp - d2Jdz2_km)/2.;

	if (j+1 > N2-1) d2Jdxdy_jp = dJdx[i][j][k];
	else d2Jdxdy_jp = dJdx[i][j+1][k];
	if (j-1 < 0) d2Jdxdy_jm = dJdx[i][j][k];
	else d2Jdxdy_jm = dJdx[i][j-1][k];
	d2Jdxdy = (d2Jdxdy_jp - d2Jdxdy_jm)/2.;

	if (k+1 > N3-1) d2Jdxdz_kp = dJdx[i][j][k];
	else d2Jdxdz_kp = dJdx[i][j][k+1];
	if (k-1 < 0) d2Jdxdz_km = dJdx[i][j][k];
	else d2Jdxdz_km = dJdx[i][j][k-1];
	d2Jdxdz = (d2Jdxdz_kp - d2Jdxdz_km)/2.;

	if (i+1 > N1-1) d2Jdydx_ip = dJdy[i][j][k];
	else d2Jdydx_ip = dJdy[i+1][j][k];
	if (i-1 < 0) d2Jdydx_im = dJdy[i][j][k];
	else d2Jdydx_im = dJdy[i-1][j][k];
	d2Jdydx = (d2Jdydx_ip - d2Jdydx_im)/2.;

	if (k+1 > N3-1) d2Jdydz_kp = dJdy[i][j][k];
	else d2Jdydz_kp = dJdy[i][j][k+1];
	if (k-1 < 0) d2Jdydz_km = dJdy[i][j][k];
	else d2Jdydz_km = dJdy[i][j][k-1];
	d2Jdydz = (d2Jdydz_kp - d2Jdydz_km)/2.;

	if (i+1 > N1-1) d2Jdzdx_ip = dJdz[i][j][k];
	else d2Jdzdx_ip = dJdz[i+1][j][k];
	if (i-1 < 0) d2Jdzdx_im = dJdz[i][j][k];
	else d2Jdzdx_im = dJdz[i-1][j][k];
	d2Jdzdx = (d2Jdzdx_ip - d2Jdzdx_im)/2.;

	if (j+1 > N2-1) d2Jdzdy_jp = dJdz[i][j][k];
	else d2Jdzdy_jp = dJdz[i][j+1][k];
	if (j-1 < 0) d2Jdzdy_jm = dJdz[i][j][k];
	else d2Jdzdy_jm = dJdz[i][j-1][k];
	d2Jdzdy = (d2Jdzdy_jp - d2Jdzdy_jm)/2.;

	Hess3D[0][0] = d2Jdx2 - grid_conn[i][j][k][1][1][1]*dJdx[i][j][k] -
							grid_conn[i][j][k][2][1][1]*dJdy[i][j][k] -
							grid_conn[i][j][k][3][1][1]*dJdz[i][j][k];

	Hess3D[0][1] = d2Jdxdy - grid_conn[i][j][k][1][1][2]*dJdx[i][j][k] -
							 grid_conn[i][j][k][2][1][2]*dJdy[i][j][k] -
							 grid_conn[i][j][k][3][1][2]*dJdz[i][j][k];

	Hess3D[0][2] = d2Jdxdz - grid_conn[i][j][k][1][1][3]*dJdx[i][j][k] -
							 grid_conn[i][j][k][2][1][3]*dJdy[i][j][k] -
							 grid_conn[i][j][k][3][1][3]*dJdz[i][j][k];

	Hess3D[1][0] = d2Jdydx - grid_conn[i][j][k][1][2][1]*dJdx[i][j][k] -
						     grid_conn[i][j][k][2][2][1]*dJdy[i][j][k] -
						     grid_conn[i][j][k][3][2][1]*dJdz[i][j][k];

	Hess3D[1][1] = d2Jdy2 - grid_conn[i][j][k][1][2][2]*dJdx[i][j][k] -
						    grid_conn[i][j][k][2][2][2]*dJdy[i][j][k] -
						    grid_conn[i][j][k][3][2][2]*dJdz[i][j][k];

	Hess3D[1][2] = d2Jdydz - grid_conn[i][j][k][1][2][3]*dJdx[i][j][k] -
							 grid_conn[i][j][k][2][2][3]*dJdy[i][j][k] -
							 grid_conn[i][j][k][3][2][3]*dJdz[i][j][k];

	Hess3D[2][0] = d2Jdzdx - grid_conn[i][j][k][1][3][1]*dJdx[i][j][k] -
						     grid_conn[i][j][k][2][3][1]*dJdy[i][j][k] -
						     grid_conn[i][j][k][3][3][1]*dJdz[i][j][k];

	Hess3D[2][1] = d2Jdzdy - grid_conn[i][j][k][1][3][2]*dJdx[i][j][k] -
						     grid_conn[i][j][k][2][3][2]*dJdy[i][j][k] -
						     grid_conn[i][j][k][3][3][2]*dJdz[i][j][k];

	Hess3D[2][2] = d2Jdz2 - grid_conn[i][j][k][1][3][3]*dJdx[i][j][k] -
							grid_conn[i][j][k][2][3][3]*dJdy[i][j][k] -
						    grid_conn[i][j][k][3][3][3]*dJdz[i][j][k];
}


double *get_evec(double **Hess)
{
    // returns eigenvector corresponding to highest eigenvalue
	// based on https://www.gnu.org/software/gsl/doc/html/eigen.html
	// see also https://www.tutorialspoint.com/cprogramming/c_return_arrays_from_function.htm

	int i, j, k, ind_max, ind_mid, ind_min, n;
	double *hessian;
	double max_eigval = -SMALL;
	double min_eigval = 1e20;
	static double *max_eigvec, *mid_eigvec, *min_eigvec;

	if (N3 > 1) n = 3;
	else n = 2;
	hessian = (double *) malloc(n*n*sizeof(double*));

	k = 0;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			hessian[k] = Hess[i][j];
			k++;
		}
	}

    max_eigvec = (double *) malloc(n*sizeof(double*));
	mid_eigvec = (double *) malloc(n*sizeof(double*));
	min_eigvec = (double *) malloc(n*sizeof(double*));
	eigvec =  = (double *) malloc(n*n*sizeof(double*));

	gsl_matrix_view m = gsl_matrix_view_array(hessian, n, n);
	gsl_vector *eval = gsl_vector_alloc(n);
	gsl_matrix *evec = gsl_matrix_alloc(n, n);
	gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(n);

	gsl_eigen_symmv(&m.matrix, eval, evec, w);
	gsl_eigen_symmv_free(w);
	gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);

	for (i = 0; i < n; i++) {
    	double eval_i = gsl_vector_get(eval, i);
    	//gsl_vector_view evec_i = gsl_matrix_column(evec, i);

        //printf ("eigenvalue = %g\n", eval_i);
        //printf ("eigenvector = \n");
        //gsl_vector_fprintf(stdout, &evec_i.vector, "%g");

        ind_max = gsl_vector_max_index(eval);
		if (eval_i > max_eigval) {
			max_eigval = eval_i;
            for (j = 0; j < n; j++) {
                max_eigvec[j] = gsl_matrix_get(evec,j,ind_max);
            }
		}

		ind_min = gsl_vector_min_index(eval);
		if (eval_i < min_eigval) {
			min_eigval = eval_i;
            for (j = 0; j < n; j++) {
                min_eigvec[j] = gsl_matrix_get(evec,j,ind_min);
            }
		}

		if (ind_max+ind_min == 1) {
			ind_mid = 2;
		} else if (ind_max+ind_min == 2) {
			ind_mid = 1;
		} else {
			ind_mid = 0;
		}
		for (j = 0; j < n; j++) {
			mid_eigvec[j] = gsl_matrix_get(evec,j,ind_mid);
		}
	}

	if (max_eigval < SMALL) {
		for (i = 0; i < n; i++) {
			max_eigvec[i] = 0.;
		}
	}
	// don't forget to add in main code a condition to skip analysis
	// if max_eigval < SMALL

    gsl_vector_free(eval);
    gsl_matrix_free(evec);
	free(hessian);

	for (i = 0; i < n; i++) {
		eigvec[i] = max_eigvec[i];
		eigvec[i+3] = mid_eigvec[i];
		eigvec[i+6] = min_eigvec[i];
	}

    return eigvec;
}



void get_hessian2D(int i, int j, double Hess2D[2][2])
{
    // IMPORTANT WHEN CALLING THIS FUNCTION! (not anymore, we changed to void)
    // https://www.tutorialspoint.com/cprogramming/c_return_arrays_from_function.htm

	int m, n;
	double dJdx_ip, dJdx_im, dJdy_jp, dJdy_jm;
	double d2Jdx2_ip, d2Jdx2_im, d2Jdy2_jp, d2Jdy2_jm;
	double d2Jdxdy_jp, d2Jdxdy_jm, d2Jdydx_ip, d2Jdydx_im;
	double ***dJdx, ***dJdy, ***d2Jdx2, ***d2Jdy2, ***d2Jdxdy, ***d2Jdydx;

	dJdx = (double ***) malloc_rank3_cont(N1, N2, N3);
	dJdy = (double ***) malloc_rank3_cont(N1, N2, N3);
	d2Jdx2 = (double ***) malloc_rank3_cont(N1, N2, N3);
	d2Jdy2 = (double ***) malloc_rank3_cont(N1, N2, N3);
	d2Jdxdy = (double ***) malloc_rank3_cont(N1, N2, N3);
	d2Jdydx = (double ***) malloc_rank3_cont(N1, N2, N3);

	for (m = 0; m < 2; m++) {
		for (n = 0; n < 2; n++) {
			Hess2D[m][n] = 0.;
		}
	}

	//1st derivative
	if (i+1 > N1-1) dJdx_ip = J_cs[i][j][0];
	else dJdx_ip = J_cs[i+1][j][0];
	if (i-1 < 0) dJdx_im = J_cs[i][j][0];
	else dJdx_im = J_cs[i-1][j][0];
	dJdx[i][j][0] = (dJdx_ip - dJdx_im)/2.;

	if (j+1 > N2-1) dJdy_jp = J_cs[i][j][0];
	else dJdy_jp = J_cs[i][j+1][0];
	if (j-1 < 0) dJdy_jm = J_cs[i][j][0];
	else dJdy_jm = J_cs[i][j-1][0];
	dJdy[i][j][0] = (dJdy_jp - dJdy_jm)/2.;

	//2nd derivative (Hessian)
	if (i+1 > N1-1) d2Jdx2_ip = dJdx[i][j][0];
	else d2Jdx2_ip = dJdx[i+1][j][0];
	if (i-1 < 0) d2Jdx2_im = dJdx[i][j][0];
	else d2Jdx2_im = dJdx[i-1][j][0];
	d2Jdx2[i][j][0] = (d2Jdx2_ip - d2Jdx2_im)/2.;

	if (j+1 > N2-1) d2Jdy2_jp = dJdy[i][j][0];
	else d2Jdy2_jp = dJdy[i][j+1][0];
	if (j-1 < 0) d2Jdy2_jm = dJdy[i][j][0];
	else d2Jdy2_jm = dJdy[i][j-1][0];
	d2Jdy2[i][j][0] = (d2Jdy2_jp - d2Jdy2_jm)/2.;

	if (j+1 > N2-1) d2Jdxdy_jp = dJdx[i][j][0];
	else d2Jdxdy_jp = dJdx[i][j+1][0];
	if (j-1 < 0) d2Jdxdy_jm = dJdx[i][j][0];
	else d2Jdxdy_jm = dJdx[i][j-1][0];
	d2Jdxdy[i][j][0] = (d2Jdxdy_jp - d2Jdxdy_jm)/2.;

	if (i+1 > N1-1) d2Jdydx_ip = dJdy[i][j][0];
	else d2Jdydx_ip = dJdy[i+1][j][0];
	if (i-1 < 0) d2Jdydx_im = dJdy[i][j][0];
	else d2Jdydx_im = dJdy[i-1][j][0];
	d2Jdydx[i][j][0] = (d2Jdydx_ip - d2Jdydx_im)/2.;

	Hess2D[0][0] = d2Jdx2[i][j][0] - grid_conn[i][j][0][1][1][1]*dJdx[i][j][0] -
								     grid_conn[i][j][0][2][1][1]*dJdy[i][j][0];

	Hess2D[0][1] = d2Jdxdy[i][j][0] - grid_conn[i][j][0][1][1][2]*dJdx[i][j][0] -
								      grid_conn[i][j][0][2][1][2]*dJdy[i][j][0];

	Hess2D[1][0] = d2Jdy2[i][j][0] - grid_conn[i][j][0][1][2][1]*dJdx[i][j][0] -
								     grid_conn[i][j][0][2][2][1]*dJdy[i][j][0];

	Hess2D[1][1] = d2Jdydx[i][j][0] - grid_conn[i][j][0][1][2][2]*dJdx[i][j][0] -
								      grid_conn[i][j][0][2][2][2]*dJdy[i][j][0];

    free(dJdx);
    free(dJdy);
	free(d2Jdx2);
    free(d2Jdy2);
	free(d2Jdxdy);
    free(d2Jdydx);
}

double *get_evec2D(double hessian[4])
{
    // returns eigenvector corresponding to highest eigenvalue
	// based on https://www.gnu.org/software/gsl/doc/html/eigen.html

	int i, j, ind_max;
	double max_eigval = SMALL;
	static double *max_eigvec;


    max_eigvec = (double *) malloc(2*sizeof(double*));

	gsl_matrix_view m = gsl_matrix_view_array(hessian, 2, 2);
	gsl_vector *eval = gsl_vector_alloc(2);
	gsl_matrix *evec = gsl_matrix_alloc(2, 2);
	gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(2);

	gsl_eigen_symmv(&m.matrix, eval, evec, w);
	gsl_eigen_symmv_free(w);
	gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);

	for (i = 0; i < 2; i++) {
    	double eval_i = gsl_vector_get(eval, i);
    	//gsl_vector_view evec_i = gsl_matrix_column(evec, i);

        //printf ("eigenvalue = %g\n", eval_i);
        //printf ("eigenvector = \n");
        //gsl_vector_fprintf(stdout, &evec_i.vector, "%g");

        ind_max = gsl_vector_max_index(eval);
		if (eval_i > max_eigval) {
			max_eigval = eval_i;
            for (j = 0; j < 2; j++) {
                max_eigvec[j] = gsl_matrix_get(evec,j,ind_max);
            }
		}
	}

	if (max_eigval < SMALL) {
		max_eigvec[0] = 0.;
		max_eigvec[1] = 0.;
		// don't forget to add in main code a condition to skip analysis
		// if max_eigval < SMALL
	}

    gsl_vector_free(eval);
    gsl_matrix_free(evec);
	free(hessian);

    return max_eigvec;
}


double *get_evec3D(double hessian[9])
{
    // returns eigenvector corresponding to highest eigenvalue
	// based on https://www.gnu.org/software/gsl/doc/html/eigen.html

	int i, j, ind_max;
	double max_eigval = SMALL;
    static double *max_eigvec;

	max_eigvec = (double *) malloc(3*sizeof(double*));

	gsl_matrix_view m = gsl_matrix_view_array(hessian, 3, 3);
	gsl_vector *eval = gsl_vector_alloc(3);
	gsl_matrix *evec = gsl_matrix_alloc(3, 3);
	gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(3);

	gsl_eigen_symmv(&m.matrix, eval, evec, w);
	gsl_eigen_symmv_free(w);
	gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);

	for (i = 0; i < 3; i++) {
    	double eval_i = gsl_vector_get(eval, i);
    	//gsl_vector_view evec_i = gsl_matrix_column(evec, i);

        //printf ("eigenvalue = %g\n", eval_i);
        //printf ("eigenvector = \n");
        //gsl_vector_fprintf(stdout, &evec_i.vector, "%g");

        ind_max = gsl_vector_max_index(eval);
		if (eval_i > max_eigval) {
			max_eigval = eval_i;
            for (j = 0; j < 3; j++) {
                max_eigvec[j] = gsl_matrix_get(evec,j,ind_max);
            }
		}
	}
	if (max_eigval < SMALL) {
		max_eigvec[0] = 0.;
		max_eigvec[1] = 0.;
		max_eigvec[2] = 0.;
		// don't forget to add in main code a condition to skip analysis
		// if max_eigval < SMALL
	}

    gsl_vector_free(eval);
    gsl_matrix_free(evec);
	free(hessian);

    return max_eigvec;
}
