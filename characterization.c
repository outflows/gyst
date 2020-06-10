/***********************************************

NOTE!!! THIS IS UNFINISHED AND PRONE TO ERRORS!

DO NOT USE THIS!

***********************************************/

// METHOD 1: NORMAL
// METHOD 2: HESSIAN

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "definitions.h"
#include "declarations.h"


void find_normal()
{
    int sheet_box_lower_i, sheet_box_lower_j, sheet_box_lower_k;
    int sheet_box_upper_i, sheet_box_upper_j, sheet_box_upper_k;
    int count, i, j, k, ii, jj, kk;
    int halfboxsize;
    double z_c, r_c, z_p, r_p;
    double sum, at, th_mean, m, ***m_perp;

    m_perp = (double ***) malloc_rank3_cont(N1, N2, N3);

    printf("Finding slopes and normals...\n");

    // initialize m_perp
    for (i = 0; i < N1; i++) {
        for (j = 0; j < N2; j++) {
            for (k = 0; k < N3; k++) {
                m_perp[i][j][k] = 0;
            }
        }
    }

    for (i = 0; i < N1; i++) {
        for (j = 0; j < N2; j++) {
            for (k = 0; k < N3; k++) {
                if (J_cs[i][j][k] != 0) {
                    // first, we must check the limits to avoid issues if the
                    // point is too close to the edges

                    halfboxsize = SHEET_HALF_BOX_SIZE;
                    check_box_limits(i, j, k, sheet_box_lower_i, sheet_box_lower_j, sheet_box_lower_k, sheet_box_upper_i, sheet_box_upper_j, sheet_box_upper_k, halfboxsize);

                    sum = 0.0;
                    count = 0;
                    z_c = r[i][j][k]*cos(phi[i][j][k]);
                    r_c = r[i][j][k]*sin(phi[i][j][k]);
                    // now, we check the points within the box
                    for (ii = sheet_box_lower_i; ii < sheet_box_upper_i + 1; ii++) {
                        for (jj = sheet_box_lower_j; jj < sheet_box_upper_j + 1; jj++) {
                            for (kk = sheet_box_lower_k; kk < sheet_box_upper_k + 1; kk++) {
                                if (J_cs[ii][jj][kk] != 0) {
                                    z_p = r[ii][jj][kk]*cos(phi[ii][jj][kk]);
                                    r_p = r[ii][jj][kk]*sin(phi[ii][jj][kk]);
                                    at = atan((z_p - z_c)/(r_p - r_c));
                                    if (isnan(at)) {
                                        at = 0;
                                    }
                                    sum = sum + at;
                                    count++;
                                }
                            }
                        }
                    }
                    th_mean = sum/count;
                    m = tan(th_mean); // slope at central point
                    m_perp[i][j][k] = - 1./m; // normal to central point
                }
            }
        }
    }
    //return m_perp;


    //FILE *fp;
    //fp = fopen(fname, "w");
    //printf("Writing normals into file '%s'...\n", fname);
    //for (i = 0; i < N1; i++) {
    //    for (j = 0; j < N2; j++) {
    //        for (k = 0; k < N3; k++) {
    //            fprintf(fp, "%lf\n", m_perp[i][j][k]);
    //        }
    //    }
    //}
    printf("Done!\n\n");

    fclose(fp);

}

void check_box_limits(int i, int j, int k, int box_lower_i, int box_lower_j, int box_lower_k, int box_upper_i, int box_upper_j, int box_upper_k, int halfboxsize) {

    // i limits
    if (i - halfboxsize < 0)
        box_lower_i = 0;
    else
        box_lower_i = i - halfboxsize;
    if (i + halfboxsize >= N1)
        box_upper_i = N1 - 1;
    else
        box_upper_i = i + halfboxsize;

    // j limits
    if (j - halfboxsize < 0)
        box_lower_j = 0;
    else
        box_lower_j = j - halfboxsize;
    if (j + halfboxsize >= N2)
        box_upper_j = N2 - 1;
    else
        box_upper_j = j + halfboxsize;

    if (N3 > 1) {
        // k limits
        if (i - halfboxsize < 0)
            box_lower_k = 0;
        else
            box_lower_k = k - halfboxsize;
        if (k + halfboxsize >= N3)
            box_upper_k = N3 - 1;
        else
            box_upper_k = k + halfboxsize;
    }

}


/*
// FIND points at same normal

// get b at cell (r,z), (interpolate)
// increase distances
// get b at new cells (interpolate)

r_plus = r_c + d_plus/(sqrt(1. + m_perp*m_perp));
r_minus = r_c - d_minus/(sqrt(1. + m_perp*m_perp));
z_plus = z_c + m_perp*(r_plus - r_c);
z_minus = z_c + m_perp*(r_minus - r_c);
*/
