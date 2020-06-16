#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "definitions.h"
#include "declarations.h"

double *malloc_rank1(int n1, int size)
{
    // Allocates a 1D array

	void *A;

	if ((A = (void *) malloc(n1 * size)) == NULL) {
		fprintf(stderr, "malloc failure in malloc_rank1\n");
		exit(123);
	}

	return A;
}


double **malloc_rank2_cont(int n1, int n2)
{
    // Allocates a 2D array

	double **A;
	double *space;
	int i;

	space = malloc_rank1(n1 * n2, sizeof(double));

	A = malloc_rank1(n1, sizeof(double *));

	for (i = 0; i < n1; i++)
		A[i] = &(space[i * n2]);

	return A;
}


double ***malloc_rank3_cont(int n1, int n2, int n3)
{
    // Allocates a 3D array

	double*** arr;
	int i,j;

	arr = (double ***) malloc(n1*sizeof(double**));
	for (i = 0; i < n1; i++) {
		arr[i] = (double **) malloc(n2*sizeof(double*));
        for(j = 0; j < n2; j++) {
        	arr[i][j] = (double *) malloc(n3 * sizeof(double));
        }
    }

	return arr;
}


double ****malloc_rank4_cont(int n1, int n2, int n3, int n4)
{
    // Allocates a 4D array

	double***** arr;
	int i,j,k;

	arr = (double *****) malloc(n1*sizeof(double**));
	for (i = 0; i < n1; i++) {
		arr[i] = (double ****) malloc(n2*sizeof(double*));
        for(j = 0; j < n2; j++) {
        	arr[i][j] = (double ***) malloc(n3 * sizeof(double));
            for(k = 0; k < n3; k++) {
            	arr[i][j][k] = (double **) malloc(n4 * sizeof(double));
            }
        }
    }

	return arr;
}


double *****malloc_rank5_cont(int n1, int n2, int n3, int n4, int n5)
{
    // Allocates a 5D array

	double***** arr;
	int i,j,k,l;

	arr = (double *****) malloc(n1*sizeof(double**));
	for (i = 0; i < n1; i++) {
		arr[i] = (double ****) malloc(n2*sizeof(double*));
        for(j = 0; j < n2; j++) {
        	arr[i][j] = (double ***) malloc(n3 * sizeof(double));
            for(k = 0; k < n3; k++) {
            	arr[i][j][k] = (double **) malloc(n4 * sizeof(double));
                for(l = 0; l < n4; l++) {
                	arr[i][j][k][l] = (double *) malloc(n5 * sizeof(double));
                }
            }
        }
    }

	return arr;
}


double ******malloc_rank6_cont(int n1, int n2, int n3, int n4, int n5, int n6)
{
    // Allocates a 6D array

	double****** arr;
	int i,j,k,l,m;

    arr = (double *****) malloc(n1*sizeof(double**));
	for (i = 0; i < n1; i++) {
		arr[i] = (double ****) malloc(n2*sizeof(double*));
        for(j = 0; j < n2; j++) {
        	arr[i][j] = (double ***) malloc(n3 * sizeof(double));
            for(k = 0; k < n3; k++) {
            	arr[i][j][k] = (double **) malloc(n4 * sizeof(double));
                for(l = 0; l < n4; l++) {
                	arr[i][j][k][l] = (double *) malloc(n5 * sizeof(double));
                    for(m = 0; m < n5; m++) {
                    	arr[i][j][k][l][m] = (double *) malloc(n6 * sizeof(double));
                    }
                }
            }
        }
    }

	return arr;
}


void init_storage_data()
{
    printf("Allocating memory for data...\n");

    flag = (int ***) malloc_rank3_cont(N1, N2, N3);
    flag_buffer = (int ***) malloc_rank3_cont(N1, N2, N3);
    a_r = (double ***) malloc_rank3_cont(N1, N2, N3);
    a_th = (double ***) malloc_rank3_cont(N1, N2, N3);
    a_phi = (double ***) malloc_rank3_cont(N1, N2, N3);
    B = (double ****) malloc_rank4_cont(NDIM, N1, N2, N3);
    ucov = (double ****) malloc_rank4_cont(NDIM, N1, N2, N3);
    ucon = (double ****) malloc_rank4_cont(NDIM, N1, N2, N3);
    bcov = (double ****) malloc_rank4_cont(NDIM, N1, N2, N3);
    bcon = (double ****) malloc_rank4_cont(NDIM, N1, N2, N3);
    jcov = (double ****) malloc_rank4_cont(NDIM, N1, N2, N3);
    jcon = (double ****) malloc_rank4_cont(NDIM, N1, N2, N3);
    J = (double ***) malloc_rank3_cont(N1, N2, N3);
    J_cs = (double ***) malloc_rank3_cont(N1, N2, N3);
	J_cs_peak = (double ***) malloc_rank3_cont(N1, N2, N3);
	J_cs_char = (double ***) malloc_rank3_cont(N1, N2, N3);
    betapl = (double ***) malloc_rank3_cont(N1, N2, N3);
    sigmaphi = (double ***) malloc_rank3_cont(N1, N2, N3);
    Sigmaphi = (double ***) malloc_rank3_cont(N1, N2, N3);
	dJdx = (double ***) malloc_rank3_cont(N1, N2, N3); // must be global
	dJdy = (double ***) malloc_rank3_cont(N1, N2, N3); // must be global
	dJdz = (double ***) malloc_rank3_cont(N1, N2, N3); // must be global


}


void init_storage_metric()
{
    printf("Allocating memory for grid...\n");

    grid_gcov = (double *****) malloc_rank5_cont(N1, N2, N3, NDIM, NDIM);
    grid_gcon = (double *****) malloc_rank5_cont(N1, N2, N3, NDIM, NDIM);
    grid_gdet = (double ***) malloc_rank3_cont(N1, N2, N3);
    grid_conn = (double ******) malloc_rank6_cont(N1, N2, N3, NDIM, NDIM, NDIM);
}
