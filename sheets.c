#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "definitions.h"
#include "declarations.h"


void get_current_sheets()
{
    // Calculates current sheets according to an algorithm similar to
    // Ball et al. 2016 and Zhdankin et al. 2013.

    int i, j, k;
    double J_peak, J_peak_min, J_max, percent;
    double *J_mean_inner, *J_mean_outer;
    int N1_inside, N1_outside;

    J_mean_inner = (double *) malloc(N1*sizeof(double*));
    J_mean_outer = (double *) malloc(N1*sizeof(double*));

    fprintf(stdout, "Initializing current sheets...\n");

    for (i = 0; i < N1; i++) {
        J_mean_inner[i] = 0.;
        J_mean_outer[i] = 0.;
    }

    // Initialize cell flags, J_cs and calculate J_mean in both patches.
    // J_cs is a 3D array whose elements form the current sheets.
    // "cs" : current sheet
    // Flag notation:
    // flag[i][j][k] = 0: cell has not been checked yet
    // flag[i][j][k] = 1: cell has been checked and may be skipped
    N1_inside = N1/5;
    N1_outside = N1 - N1_inside;
    for (i = 0; i < N1; i++) {
        for (j = 0; j < N2; j++) {
            for (k = 0; k < N3; k++) {
                flag[i][j][k] = 0;
                J_cs[i][j][k] = 0.;
                if (i <= N1_inside) {
                    J_mean_inner[i] += J[i][j][k];
                } else {
                    J_mean_outer[i] += J[i][j][k];
                }
            }
        }
    }
    for (i = 0; i < N1; i++) {
        J_mean_inner[i] = J_mean_inner[i]/(N1_inside*N2*N3);
        J_mean_outer[i] = J_mean_outer[i]/(N1_outside*N2*N3);
    }

    // Now, loop over the grid a second time, this time to add to J_cs the
    // "peak" values of the current.
    // Also check for the smallest and largest of these "peak" values.
    // The reason we do this will be clear soon.
    J_peak_min = pow(N1*N2*N3, 8); // just a high enough number
    J_max = 0;
    peak_count = 0;

    for (i = 0; i < N1; i++) {
        for (j = 0; j < N2; j++) {
            for (k = 0; k < N3; k++) {
                if (i <= N1_inside) {
                    //if ((J[i][j][k] > J_FAC*J_mean_inner[i]) && (betapl[i][j][k] > BETAPL_THR)) {
                    if ((J[i][j][k] > J_FAC*J_mean_inner[i]) && (sigmaphi[i][j][k] < SIGMAPHI_THR)) {
                        peak_count++;
                        J_cs[i][j][k] = J[i][j][k];
                        if (J_cs[i][j][k] < J_peak_min) {
                            J_peak_min = J_cs[i][j][k];
                        }
                        if (J_cs[i][j][k] > J_max) {
                            J_max = J_cs[i][j][k];
                        }
                    }
                } else {
                    //if ((J[i][j][k] > J_FAC*J_mean_outer[i]) && (betapl[i][j][k] > BETAPL_THR)) {
                    if ((J[i][j][k] > J_FAC*J_mean_outer[i]) && (sigmaphi[i][j][k] < SIGMAPHI_THR)) {
                        peak_count++;
                        J_cs[i][j][k] = J[i][j][k];
                        if (J_cs[i][j][k] < J_peak_min) {
                            J_peak_min = J_cs[i][j][k];
                        }
                        if (J_cs[i][j][k] > J_max) {
                            J_max = J_cs[i][j][k];
                        }
                    }
                }
            }
        }
    }

    fprintf(stdout, "Checking initial candidates and getting local maxima...\n");
    // Loop over all J_cs points. For each point, look within a box
    // (see HALF_BOX_SIZE in definitions.h for box size) surrounding it.
    // If any point inside this box has a larger current, then the initial
    // point is not a local maximum and is excluded from the sample.
    get_locmax();

    // Calculate amount and print percentage of points which are above the
    // current sheet threshold.
    percent = (double)(peak_count)/(N1*N2*N3)*100.;
    fprintf(stdout, "%d (%.2lf%%) points are above the threshold.\n", peak_count, percent);

    // TO DO: adapt jpeak to GRMONTY grid (0, 1 thing, like flag)
    // TO DO: print to another file the 0-1 grid for the new jpeak

    // Loop to get J_cs_peak and flag cells with low current
    // After get_locmax(), J_cs contains only the peak values; it will later
    // be filled with the other points found to belong to current sheets.
    for (i = 0; i < N1; i++) {
        for (j = 0; j < N2; j++) {
            for (k = 0; k < N3; k++) {
                J_cs_peak[i][j][k] = J_cs[i][j][k];
                if (J[i][j][k] < JPEAK_FAC*J_peak_min || (sigmaphi[i][j][k] > SIGMAPHI_THR)) {
                    flag[i][j][k] = 1;
                }
            }
        }
    }
    // Print the maxima to a file
    //write_current_sheets(jpeak_output, J_cs_peak);
    fprintf(stdout, "Finished current sheet initialization.\n\n");

    // Main loop
    fprintf(stdout, "Getting current sheets...\n");
    fprintf(stdout, "Entering main loop...\n");
    for (i = 0; i < N1; i++) {
        for (j = 0; j < N2; j++) {
            for (k = 0; k < N3; k++) {
                if ((i <= N1_inside) &&
                    (J_cs_peak[i][j][k] > J_FAC*J_mean_inner[i]) &&
                    (sigmaphi[i][j][k] < SIGMAPHI_THR) &&
                    (flag[i][j][k] == 0)) {
                    J_peak = J_cs_peak[i][j][k];
                    fill(i, j, k, J_peak);
                } else if ((i > N1_inside) &&
                           (J_cs_peak[i][j][k] > J_FAC*J_mean_outer[i]) &&
                           (sigmaphi[i][j][k] < SIGMAPHI_THR) &&
                           (flag[i][j][k] == 0)) {
                    J_peak = J_cs_peak[i][j][k];
                    fill(i, j, k, J_peak);
                }
            }
        }
    }
    fprintf(stdout, "NOTE: normalization of current left to plotting script!\n");
    write_current_sheets(jcs_output, J_cs);
    fprintf(stdout, "Finished getting current sheets.\n");
}


void fill(int i, int j, int k, double threshold)
{
    // Flood fill algorithm to get cells belonging to current sheets

    if (BETWEEN(i, 0, N1) && BETWEEN(j, 0, N2) && BETWEEN(k, 0, N3) && (flag[i][j][k] == 0)) {
        flag[i][j][k] = 1;
        if ((J[i][j][k] > JPEAK_FAC*threshold) && (sigmaphi[i][j][k] < SIGMAPHI_THR)) {
            J_cs[i][j][k] = J[i][j][k];
            fill(i+1, j, k, threshold);
            fill(i-1, j, k, threshold);
            fill(i, j+1, k, threshold);
            fill(i, j-1, k, threshold);
            if (N3 > 1) {
                fill(i, j, k+1, threshold);
                fill(i, j, k-1, threshold);
            }
        }
    }
}


void get_locmax()
{
    int i, j, k, ii, jj, kk, endloop;
    int box_lower_i, box_lower_j, box_lower_k;
    int box_upper_i, box_upper_j, box_upper_k;

    for (i = 0; i < N1; i++) {
        for (j = 0; j < N2; j++) {
            for (k = 0; k < N3; k++) {
                if (J_cs[i][j][k] != 0) {
                    // first, we must check the limits to avoid issues if the
                    // point is too close to the edges

                    // i limits
                    if (i - HALF_BOX_SIZE < 0) box_lower_i = 0;
                    else box_lower_i = i - HALF_BOX_SIZE;
                    if (i + HALF_BOX_SIZE >= N1) box_upper_i = N1 - 1;
                    else box_upper_i = i + HALF_BOX_SIZE;

                    // j limits
                    if (j - HALF_BOX_SIZE < 0) box_lower_j = 0;
                    else box_lower_j = j - HALF_BOX_SIZE;
                    if (j + HALF_BOX_SIZE >= N2) box_upper_j = N2 - 1;
                    else box_upper_j = j + HALF_BOX_SIZE;

                    // k limits
                    if (k - HALF_BOX_SIZE < 0) box_lower_k = 0;
                    else box_lower_k = k - HALF_BOX_SIZE;
                    if (k + HALF_BOX_SIZE >= N3) box_upper_k = N3 - 1;
                    else box_upper_k = k + HALF_BOX_SIZE;

                    // now, we check the points within the box
                    for (ii = box_lower_i; ii < box_upper_i + 1; ii++) {
                        for (jj = box_lower_j; jj < box_upper_j + 1; jj++) {
                            for (kk = box_lower_k; kk < box_upper_k + 1; kk++) {
                                if (J[ii][jj][kk] > J_cs[i][j][k]) {
                                    J_cs[i][j][k] = 0.;
                                    peak_count--;
                                    endloop = 1;
                                }
                                // thou shalt not goto
                                if (endloop) break;
                            }
                            // thou shalt not goto
                            if (endloop) break;
                        }
                        // thou shalt not goto
                        if (endloop) break;
                    }
                }
            }
        }
    }
}
