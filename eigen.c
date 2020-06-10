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

double *get_hessian2D(int i, int j)
{
    // IMPORTANT WHEN CALLING THIS FUNCTION!
    // https://www.tutorialspoint.com/cprogramming/c_return_arrays_from_function.htm

	double ip, im, jp, jm;
	double ***Dx, ***Dy;
	static double H[4];

	Dx = (double ***) malloc_rank3_cont(N1, N2, N3);
	Dy = (double ***) malloc_rank3_cont(N1, N2, N3);

	//1st derivative
	if (i+1 > N1) ip = J_cs[i][j][0];
	else ip = J_cs[i+1][j][0];
	if (i-1 < N1) im = J_cs[i][j][0];
	else im = J_cs[i-1][j][0];
	if (j+1 > N2) jp = J_cs[i][j][0];
	else jp = J_cs[i][j+1][0];
	if (j-1 < N2) jm = J_cs[i][j][0];
	else jm = J_cs[i][j-1][0];

	Dx[i][j][0] = (ip - im)/2.;
	Dy[i][j][0] = (jp - jm)/2.;

	//2nd derivative (Hessian)

	// second derivative of Dx: Hxx, Hxy
	if (i+1 > N1) ip = Dx[i][j][0];
	else ip = Dx[i+1][j][0];
	if (i-1 < N1) im = Dx[i][j][0];
	else im = Dx[i-1][j][0];
	if (j+1 > N2) jp = Dx[i][j][0];
	else jp = Dx[i][j+1][0];
	if (j-1 < N2) jm = Dx[i][j][0];
	else jm = Dx[i][j-1][0];

	H[0] = (ip - im)/2.;
	H[1] = (jp - jm)/2.;

	// second derivative of Dy: Hyx, Hyy
	if (i+1 > N1) ip = Dy[i][j][0];
	else ip = Dy[i+1][j][0];
	if (i-1 < N1) im = Dy[i][j][0];
	else im = Dy[i-1][j][0];
	if (j+1 > N2) jp = Dy[i][j][0];
	else jp = Dy[i][j+1][0];
	if (j-1 < N2) jm = Dy[i][j][0];
	else jm = Dy[i][j-1][0];

	H[2] = (ip - im)/2.;
	H[3] = (jp - jm)/2.;

    free(Dx);
    free(Dy);

	return H;
}

double *get_hessian3D(int i, int j, int k)
{
    // IMPORTANT WHEN CALLING THIS FUNCTION!
    // https://www.tutorialspoint.com/cprogramming/c_return_arrays_from_function.htm

	double ip, im, jp, jm, kp, km;
	double ***Dx, ***Dy, ***Dz;
	static double H[9];

	Dx = (double ***) malloc_rank3_cont(N1, N2, N3);
	Dy = (double ***) malloc_rank3_cont(N1, N2, N3);
	Dz = (double ***) malloc_rank3_cont(N1, N2, N3);

	//1st derivative
	if (i+1 > N1) ip = J_cs[i][j][k];
	else ip = J_cs[i+1][j][k];
	if (i-1 < N1) im = J_cs[i][j][k];
	else im = J_cs[i-1][j][k];
	if (j+1 > N2) jp = J_cs[i][j][k];
	else jp = J_cs[i][j+1][k];
	if (j-1 < N2) jm = J_cs[i][j][k];
	else jm = J_cs[i][j-1][k];
	if (k+1 > N3) kp = J_cs[i][j][k];
	else kp = J_cs[i][j][k+1];
	if (k-1 < N3) km = J_cs[i][j][k];
	else km = J_cs[i][j][k-1];

	Dx[i][j][k] = (ip - im)/2.;
	Dy[i][j][k] = (jp - jm)/2.;
	Dz[i][j][k] = (kp - km)/2.;

	//2nd derivative (Hessian)

	// second derivative of Dx: Hxx, Hxy, Hxz
	if (i+1 > N1) ip = Dx[i][j][k];
	else ip = Dx[i+1][j][k];
	if (i-1 < N1) im = Dx[i][j][k];
	else im = Dx[i-1][j][k];
	if (j+1 > N2) jp = Dx[i][j][k];
	else jp = Dx[i][j+1][k];
	if (j-1 < N2) jm = Dx[i][j][k];
	else jm = Dx[i][j-1][k];
	if (k+1 > N3) kp = Dx[i][j][k];
	else kp = Dx[i][j][k+1];
	if (k-1 < N3) km = Dx[i][j][k];
	else km = Dx[i][j][k-1];

	H[0] = (ip - im)/2.;
	H[1] = (jp - jm)/2.;
	H[2] = (kp - km)/2.;

	// second derivative of Dy: Hyx, Hyy, Hyz
	if (i+1 > N1) ip = Dy[i][j][k];
	else ip = Dy[i+1][j][k];
	if (i-1 < N1) im = Dy[i][j][k];
	else im = Dy[i-1][j][k];
	if (j+1 > N2) jp = Dy[i][j][k];
	else jp = Dy[i][j+1][k];
	if (j-1 < N2) jm = Dy[i][j][k];
	else jm = Dy[i][j-1][k];
	if (k+1 > N3) kp = Dy[i][j][k];
	else kp = Dy[i][j][k+1];
	if (k-1 < N3) km = Dy[i][j][k];
	else km = Dy[i][j][k-1];

	H[3] = (ip - im)/2.;
	H[4] = (jp - jm)/2.;
	H[5] = (kp - km)/2.;

	// second derivative of Dz: Hzx, Hzy, Hzz
	if (i+1 > N1) ip = Dz[i][j][k];
	else ip = Dz[i+1][j][k];
	if (i-1 < N1) im = Dz[i][j][k];
	else im = Dz[i-1][j][k];
	if (j+1 > N2) jp = Dz[i][j][k];
	else jp = Dz[i][j+1][k];
	if (j-1 < N2) jm = Dz[i][j][k];
	else jm = Dz[i][j-1][k];
	if (k+1 > N3) kp = Dz[i][j][k];
	else kp = Dz[i][j][k+1];
	if (k-1 < N3) km = Dz[i][j][k];
	else km = Dz[i][j][k-1];

	H[6] = (ip - im)/2.;
	H[7] = (jp - jm)/2.;
	H[8] = (kp - km)/2.;

    free(Dx);
    free(Dy);
    free(Dz);

	return H;
}

double *get_evec(double *hessian) {
    // returns eigenvector corresponding to highest eigenvalue
	// based on https://www.gnu.org/software/gsl/doc/html/eigen.html
	// see also https://www.tutorialspoint.com/cprogramming/c_return_arrays_from_function.htm

	int i, j, ind_max, n;
    static double *max_eigvec;
	double max_eigval = SMALL;

	if (N3 > 1) n = 3;
	else n = 2;

    max_eigvec = (double *) malloc(n*sizeof(double*));

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
	}

    gsl_vector_free(eval);
    gsl_matrix_free(evec);

    return max_eigvec;
}



double *get_evec2D(double hessian[4]) {
    // returns eigenvector corresponding to highest eigenvalue
	// based on https://www.gnu.org/software/gsl/doc/html/eigen.html

	int i, j, ind_max;
    static double *max_eigvec;
	double max_eigval = SMALL;

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

    gsl_vector_free(eval);
    gsl_matrix_free(evec);

    return max_eigvec;
}

double *get_evec3D(double hessian[9]) {
    // returns eigenvector corresponding to highest eigenvalue
	// based on https://www.gnu.org/software/gsl/doc/html/eigen.html

	int i, j, ind_max;
    static double *max_eigvec;
	double max_eigval = SMALL;

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

    gsl_vector_free(eval);
    gsl_matrix_free(evec);

    return max_eigvec;
}
