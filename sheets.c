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

    int i, j, k, m, ii, jj, kk, mm;
    int count;
    double percent;
    int cells, dead_cells, max_cells;
    int size_current, i_ctr, j_ctr, k_ctr, i_adj, j_adj, k_adj;
    double *J_mean, J_peak, J_peak_min, J_max;
    double *J_mean_inner, *J_mean_outer;

    J_mean = (double *) malloc(N1*sizeof(double*));
    J_mean_inner = (double *) malloc(N1*sizeof(double*));
    J_mean_outer = (double *) malloc(N1*sizeof(double*));

    printf("Initializing current sheets...\n");
    // Initialize cell flags, J_cs and calculate J_mean.
    // J_cs is a 3D array whose elements form the current sheets.
    //
    // Flag notation:
    // flag = 0: cell has not been checked yet
    // flag = 1: cell has been checked and may be skipped

    int N1_inside = N1/5;
    int N1_outside = N1 - N1_inside;
    //J_mean = 0;
    for (i = 0; i < N1; i++) {
        J_mean[i] = 0.;
        J_mean_inner[i] = 0.;
        J_mean_outer[i] = 0.;
    }
    //J_nonzero = 0;
    for (i = 0; i < N1; i++) {
        for (j = 0; j < N2; j++) {
            for (k = 0; k < N3; k++) {
                flag[i][j][k] = 0;
                flag_buffer[i][j][k] = 0;
                J_cs[i][j][k] = 0.;
                J_mean[i] += J[i][j][k];
                if (i <= N1_inside) {
                    J_mean_inner[i] += J[i][j][k];
                } else {
                    J_mean_outer[i] += J[i][j][k];
                }

                //if (J[i][j][k] > 0.) {
                //    J_mean += J[i][j][k];
                //    J_nonzero++;
                //}
            }
        }
    }
    //J_mean = J_mean/(N1*N2*N3);
    for (i = 0; i < N1; i++) {
        J_mean[i] = J_mean[i]/(N2*N3);
        J_mean_inner[i] = J_mean_inner[i]/(N1_inside*N2*N3);
        J_mean_outer[i] = J_mean_outer[i]/(N1_outside*N2*N3);
    }

    // Now, loop over the grid a second time, this time to add to J_cs the
    // "peak" values of the current. Note that "cs" means current sheet.
    //
    // We also check for the smallest and largest of these "peak" values.
    // The reason we do this will be clear soon.
    J_peak_min = pow(N1*N2*N3, 8); // just a high enough number
    J_max = 0;
    count = 0;

    //double aux;
    for (i = 0; i < N1; i++) {
        for (j = 0; j < N2; j++) {
            for (k = 0; k < N3; k++) {
                if (i <= N1_inside) {
                    //if ((J[i][j][k] > J_FAC*J_mean_inner[i]) && (betapl[i][j][k] > BETAPL_THR)) {
                    if ((J[i][j][k] > J_FAC*J_mean_inner[i]) && (sigmaphi[i][j][k] < SIGMAPHI_THR)) {
                        count++;
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
                        count++;
                        J_cs[i][j][k] = J[i][j][k];
                        if (J_cs[i][j][k] < J_peak_min) {
                            J_peak_min = J_cs[i][j][k];
                        }
                        if (J_cs[i][j][k] > J_max) {
                            J_max = J_cs[i][j][k];
                        }
                    }
                }

                //if (J[i][j][k] > J_FAC*J_mean[i]) {
/*              if ((J[i][j][k] > J_FAC*J_mean[i]) && (betapl[i][j][k] > BETAPL_THR)) {
                //if ((J[i][j][k] > J_FAC*J_mean[i]) && (sigmaphi[i][j][k] < SIGMAPHI_THR)) {
                //if ((J[i][j][k] > J_FAC*J_mean[i]) && (Sigmaphi[i][j][k] < SIGMAPHI_THR)) {
                //aux = bu3[i][j][k];
                //if ((J[i][j][k] > J_FAC*J_mean[i]) && (BETWEEN(aux, -BU3_THR, BU3_THR))) {
                    count++;
                    J_cs[i][j][k] = J[i][j][k];
                    if (J_cs[i][j][k] < J_peak_min) {
                        J_peak_min = J_cs[i][j][k];
                    }
                    if (J_cs[i][j][k] > J_max) {
                        J_max = J_cs[i][j][k];
                    }
                }
*/
            }
        }
    }

    printf("Checking initial candidates and getting local maxima...\n");
    // Now, we loop over all J_cs points. For each point, we look within a 7x7x7
    // box surrounding it. If any point inside this box has a larger current,
    // then our initial point is not a local maximum and we exclude it from our
    // sample.
    get_locmax();

    // Calculate amount and print percentage of points which are above the
    // current sheet threshold
    percent = (double)(count)/(N1*N2*N3)*100.;
    printf("%d (%.2lf%%) points are above the threshold.\n", count, percent);

    // TO DO: adapt jpeak to GRMONTY grid (0, 1 thing, like flag)
    // TO DO: print to another file the 0-1 grid for the new jpeak

    for (i = 0; i < N1; i++) {
        for (j = 0; j < N2; j++) {
            for (k = 0; k < N3; k++) {
                J_cs_peak[i][j][k] = J_cs[i][j][k];
            }
        }
    }
    // We now print the maxima to a file
    write_current_sheets(jpeak_output, J_cs_peak);

    // Now, we loop over the entire grid again, but this time we check if the
    // current at cell (i,j, k) is smaller than half the smallest J_peak, which
    // we found above.
    dead_cells = 0; // "dead cells" are the cells that we won't have to check
    for (i = 0; i < N1; i++) {
        for (j = 0; j < N2; j++) {
            for (k = 0; k < N3; k++) {
                if (J[i][j][k] < JPEAK_FAC*J_peak_min) {
                    flag_buffer[i][j][k] = 1;
                    dead_cells++;
                }
            }
        }
    }
    max_cells = N1*N2*N3 - dead_cells; // max number of cells we have to check
    // Since the flag at these cells is 1, we won't have to bother about them
    // when we look for the current sheets, as the minimum value of J at a
    // current sheet will ALWAYS be equal to at least half the smallest J_peak.

    // NOTE: THIS PART IS ONLY FOR THE 2D VERSION, NOT USED HERE!!!
    // Now, we create two arrays, each of which will have (later, upon malloc)
    // a number of elements equal to N1*N2*N3-dead_cells. In other words, the
    // number of elements in these arrays is equal to the number of grid cells
    // that may be part of current sheets. It is also the size of the largest
    // possible current sheet, so we can guarantee that no current sheets will
    // have more elements than N1*N2*N3-dead_cells.
    //int *icoords;// = (int *)malloc(max_cells*sizeof(int));
    //int *jcoords;// = (int *)malloc(max_cells*sizeof(int));
    //int *kcoords;// = (int *)malloc(max_cells*sizeof(int));
    // The idea is that these arrays will contain the coordinates (i,j, k) of
    // the grid cells belonging to the current sheets, so that we will loop over
    // this array instead of looping over the entire grid a number of times
    // equal to (N1*N2*N3 - dead_cells). The advantage of this method is that
    // we won't have to check each grid cell every time we are trying to find
    // the elements of our current. Instead, we will go straight to the cells
    // that must be checked. This will be useful if we apply this algorithm to
    // higher resolution 3d grids.

    printf("Finished current sheet initialization.\n\n");

    // Main loop
    printf("Entering main loop...\n");
    clock_t begin_loop = clock();
    clock_t now_loop;
    float time_spent_loop;
    for (i = 0; i < N1; i++) {
        if (i % 20 == 0) {
            now_loop = clock();
            time_spent_loop = (double)(now_loop - begin_loop) / CLOCKS_PER_SEC;
            printf("i = %d, %.2f%% of N1, t (main loop) = %.2f seconds\n", i, (float)i/N1*100, time_spent_loop);
        }
        for (j = 0; j < N2; j++) {
            for (k = 0; k < N3; k++) {

                if (i <= N1_inside) {

                    //if ((J_cs[i][j][k] > J_FAC*J_mean_inner[i]) && (betapl[i][j][k] > BETAPL_THR)) {
                    if ((J_cs[i][j][k] > J_FAC*J_mean_inner[i]) && (sigmaphi[i][j][k] < SIGMAPHI_THR)) {
                        // reset flags
                        for (ii = 0; ii < N1; ii++) {
                            for (jj = 0; jj < N2; jj++) {
                                for (kk = 0; kk < N3; kk++) {
                                    flag[ii][jj][kk] = flag_buffer[ii][jj][kk];
                                }
                            }
                        }
                        J_peak = J_cs[i][j][k];

                        // initialize icoords, jcoords, kcoords
                        mm = 0;
                        icoords = (int *)malloc(max_cells*sizeof(int));
                        jcoords = (int *)malloc(max_cells*sizeof(int));
                        if (N3 > 1) {
                            kcoords = (int *)malloc(max_cells*sizeof(int));
                        } else {
                            kcoords = (int *)malloc(1*sizeof(int));
                        }
                        icoords[mm] = i;
                        jcoords[mm] = j;
                        if (N3 > 1) {
                            kcoords[mm] = k;
                        } else {
                            kcoords[mm] = 0;
                        }
                        size_current = 1;

                        cells = 1; // number of cells left to be checked in sheet
                        // While there are cells to be checked in current sheet...
                        while(cells) {
                            cells--;

                            // ...loop over (i,j,k) coords of cells to be checked
                            for (m = 0; m < size_current; m++) {
                                i_ctr = icoords[m];
                                if (i_ctr > N1_inside || i_ctr+1 > N1_inside) break;
                                j_ctr = jcoords[m];
                                if (N3 > 1) {
                                    k_ctr = kcoords[m];
                                } else {
                                    k_ctr = 0;
                                }

                                //if (flag[i_ctr][j_ctr][k_ctr] == 0) {
                                if (flag_buffer[i_ctr][j_ctr][k_ctr] == 0) {
                                    // Check adjacency:
                                    //    _______________
                                    //   /|             /|
                                    //  / |     UU     / |
                                    // /__|___________/  |
                                    // |  |           |  |
                                    // |  |      NN   |  |
                                    // |WW|    CTR    |EE|
                                    // |  |___SS______|__|
                                    // |  /           |  /
                                    // | /     DD     | /
                                    // |/_____________|/
                                    //
                                    // straightforward to add diagonals if needed

                                    // NN
                                    i_adj = i_ctr;
                                    j_adj = j_ctr - 1;
                                    k_adj = k_ctr;
                                    if (i_adj > N1_inside) break;
                                    check_adjacency(i_adj, j_adj, k_adj,
                                                    J_peak, mm, cells, size_current,
                                                    icoords, jcoords, kcoords);

                                    // EE
                                    i_adj = i_ctr + 1;
                                    j_adj = j_ctr;
                                    k_adj = k_ctr;
                                    if (i_adj > N1_inside) break;
                                    check_adjacency(i_adj, j_adj, k_adj,
                                                    J_peak, mm, cells, size_current,
                                                    icoords, jcoords, kcoords);

                                    // SS
                                    i_adj = i_ctr;
                                    j_adj = j_ctr + 1;
                                    k_adj = k_ctr;
                                    if (i_adj > N1_inside) break;
                                    check_adjacency(i_adj, j_adj, k_adj,
                                                    J_peak, mm, cells, size_current,
                                                    icoords, jcoords, kcoords);

                                    // WW
                                    i_adj = i_ctr - 1;
                                    j_adj = j_ctr;
                                    k_adj = k_ctr;
                                    if (i_adj > N1_inside) break;
                                    check_adjacency(i_adj, j_adj, k_adj,
                                                    J_peak, mm, cells, size_current,
                                                    icoords, jcoords, kcoords);

                                    if (N3 > 1) {
                                        // UU
                                        i_adj = i_ctr;
                                        j_adj = j_ctr;
                                        k_adj = k_ctr + 1;
                                        if (i_adj > N1_inside) break;
                                        check_adjacency(i_adj, j_adj, k_adj,
                                                    J_peak, mm, cells, size_current,
                                                    icoords, jcoords, kcoords);

                                        // DD
                                        i_adj = i_ctr;
                                        j_adj = j_ctr;
                                        k_adj = k_ctr - 1;
                                        if (i_adj > N1_inside) break;
                                        check_adjacency(i_adj, j_adj, k_adj,
                                                    J_peak, mm, cells, size_current,
                                                    icoords, jcoords, kcoords);
                                    }

                                    // after checking the central cell's adjacency,
                                    // flag central cell so we don't return to it
                                    // when checking the adjacency of its own
                                    // adjacent cells
                                    flag[i_ctr][j_ctr][k_ctr] = 1;
                                    flag_buffer[i_ctr][j_ctr][k_ctr] = 1;
                                }
                            }
                        }

                        // empty icoords, jcoords and kcoords to prepare for the
                        // next current sheet
                        free(icoords);
                        free(jcoords);
                        free(kcoords);
                    }

                } else {

                    //if ((J_cs[i][j][k] > J_FAC*J_mean_outer[i]) && (betapl[i][j][k] > BETAPL_THR)) {
                    if ((J_cs[i][j][k] > J_FAC*J_mean_outer[i]) && (sigmaphi[i][j][k] < SIGMAPHI_THR)) {
                        // reset flags
                        for (ii = 0; ii < N1; ii++) {
                            for (jj = 0; jj < N2; jj++) {
                                for (kk = 0; kk < N3; kk++) {
                                    flag[ii][jj][kk] = flag_buffer[ii][jj][kk];
                                }
                            }
                        }
                        J_peak = J_cs[i][j][k];

                        // initialize icoords, jcoords, kcoords
                        mm = 0;
                        icoords = (int *)malloc(max_cells*sizeof(int));
                        jcoords = (int *)malloc(max_cells*sizeof(int));
                        if (N3 > 1) {
                            kcoords = (int *)malloc(max_cells*sizeof(int));
                        } else {
                            kcoords = (int *)malloc(1*sizeof(int));
                        }
                        icoords[mm] = i;
                        jcoords[mm] = j;
                        if (N3 > 1) {
                            kcoords[mm] = k;
                        } else {
                            kcoords[mm] = 0;
                        }
                        size_current = 1;

                        cells = 1; // number of cells left to be checked in sheet
                        // While there are cells to be checked in current sheet...
                        while(cells) {
                            cells--;

                            // ...loop over (i,j,k) coords of cells to be checked
                            for (m = 0; m < size_current; m++) {
                                i_ctr = icoords[m];
                                if (i_ctr <= N1_inside || i_ctr+1 <= N1_inside) break;
                                j_ctr = jcoords[m];
                                if (N3 > 1) {
                                    k_ctr = kcoords[m];
                                } else {
                                    k_ctr = 0;
                                }

                                //if (flag[i_ctr][j_ctr][k_ctr] == 0) {
                                if (flag_buffer[i_ctr][j_ctr][k_ctr] == 0) {
                                    // Check adjacency:
                                    //    _______________
                                    //   /|             /|
                                    //  / |     UU     / |
                                    // /__|___________/  |
                                    // |  |           |  |
                                    // |  |      NN   |  |
                                    // |WW|    CTR    |EE|
                                    // |  |___SS______|__|
                                    // |  /           |  /
                                    // | /     DD     | /
                                    // |/_____________|/
                                    //
                                    // straightforward to add diagonals if needed

                                    // NN
                                    i_adj = i_ctr;
                                    j_adj = j_ctr - 1;
                                    k_adj = k_ctr;
                                    if (i_adj <= N1_inside) break;
                                    check_adjacency(i_adj, j_adj, k_adj,
                                                    J_peak, mm, cells, size_current,
                                                    icoords, jcoords, kcoords);

                                    // EE
                                    i_adj = i_ctr + 1;
                                    j_adj = j_ctr;
                                    k_adj = k_ctr;
                                    if (i_adj <= N1_inside) break;
                                    check_adjacency(i_adj, j_adj, k_adj,
                                                    J_peak, mm, cells, size_current,
                                                    icoords, jcoords, kcoords);

                                    // SS
                                    i_adj = i_ctr;
                                    j_adj = j_ctr + 1;
                                    k_adj = k_ctr;
                                    if (i_adj <= N1_inside) break;
                                    check_adjacency(i_adj, j_adj, k_adj,
                                                    J_peak, mm, cells, size_current,
                                                    icoords, jcoords, kcoords);

                                    // WW
                                    i_adj = i_ctr - 1;
                                    j_adj = j_ctr;
                                    k_adj = k_ctr;
                                    if (i_adj <= N1_inside) break;
                                    check_adjacency(i_adj, j_adj, k_adj,
                                                    J_peak, mm, cells, size_current,
                                                    icoords, jcoords, kcoords);

                                    if (N3 > 1) {
                                        // UU
                                        i_adj = i_ctr;
                                        j_adj = j_ctr;
                                        k_adj = k_ctr + 1;
                                        if (i_adj <= N1_inside) break;
                                        check_adjacency(i_adj, j_adj, k_adj,
                                                    J_peak, mm, cells, size_current,
                                                    icoords, jcoords, kcoords);

                                        // DD
                                        i_adj = i_ctr;
                                        j_adj = j_ctr;
                                        k_adj = k_ctr - 1;
                                        if (i_adj <= N1_inside) break;
                                        check_adjacency(i_adj, j_adj, k_adj,
                                                    J_peak, mm, cells, size_current,
                                                    icoords, jcoords, kcoords);
                                    }
                                    // after checking the central cell's adjacency,
                                    // flag central cell so we don't return to it
                                    // when checking the adjacency of its own
                                    // adjacent cells
                                    flag[i_ctr][j_ctr][k_ctr] = 1;
                                    flag_buffer[i_ctr][j_ctr][k_ctr] = 1;
                                }
                            }
                        }

                        // empty icoords, jcoords and kcoords to prepare for the
                        // next current sheet
                        free(icoords);
                        free(jcoords);
                        free(kcoords);
                    }

                }

            }
        }
    }

    printf("Finished getting current sheets.\n");
    printf("NOTE: normalization of current left to plotting script!\n");
    write_current_sheets(jcs_output, J_cs);
}

void get_locmax() {
    int i, j, k, ii, jj, kk, count, endloop;
    int box_lower_i, box_lower_j, box_lower_k;
    int box_upper_i, box_upper_j, box_upper_k;
    //double local_maximum;

    for (i = 0; i < N1; i++) {
        for (j = 0; j < N2; j++) {
            for (k = 0; k < N3; k++) {
                if (J_cs[i][j][k] != 0) {
                    // first, we must check the limits to avoid issues if the
                    // point is too close to the edges

                    // i limits
                    if (i - HALF_BOX_SIZE < 0)
                        box_lower_i = 0;
                    else
                        box_lower_i = i - HALF_BOX_SIZE;
                    if (i + HALF_BOX_SIZE >= N1)
                        box_upper_i = N1 - 1;
                    else
                        box_upper_i = i + HALF_BOX_SIZE;

                    // j limits
                    if (j - HALF_BOX_SIZE < 0)
                        box_lower_j = 0;
                    else
                        box_lower_j = j - HALF_BOX_SIZE;
                    if (j + HALF_BOX_SIZE >= N2)
                        box_upper_j = N2 - 1;
                    else
                        box_upper_j = j + HALF_BOX_SIZE;

                    // k limits
                    if (k - HALF_BOX_SIZE < 0)
                        box_lower_k = 0;
                    else
                        box_lower_k = k - HALF_BOX_SIZE;
                    if (k + HALF_BOX_SIZE >= N3)
                        box_upper_k = N3 - 1;
                    else
                        box_upper_k = k + HALF_BOX_SIZE;

                    //local_maximum = J_cs[i][j][k];
                    // now, we check the points within the box
                    for (ii = box_lower_i; ii < box_upper_i + 1; ii++) {
                        for (jj = box_lower_j; jj < box_upper_j + 1; jj++) {
                            for (kk = box_lower_k; kk < box_upper_k + 1; kk++) {
                                if (J_cs[ii][jj][kk] > J_cs[i][j][k]) {
                                    J_cs[i][j][k] = 0;
                                    count--;
                                    endloop = 1;
                                }
                                if (endloop) break;
                            }
                            if (endloop) break;
                        }
                        if (endloop) break;
                    }
                }
            }
        }
    }
}

void check_adjacency (int i_adj, int j_adj, int k_adj,
                      double J_peak, int mm, int cells, int size_current,
                      int *icoords, int *jcoords, int *kcoords) {
    // Checks that the adjacent cell satisfies 2 conditions:
    // 1) is within the grid
    // 2) has not been checked yet
    // If 1, 2 satisfied, adds adjacent cell to current sheet if current at
    // adjacent cell is large enough, or flags it if it's not large enough.

    //double aux;

    if (BETWEEN(i_adj, 0, N1) && BETWEEN(j_adj, 0, N2) && BETWEEN(k_adj, 0, N3) && (flag[i_adj][j_adj][k_adj] == 0)) {
    // check if current in adjacent cell is large enough
    // for the cell to be added to the current sheet
        if (J[i_adj][j_adj][k_adj] > JPEAK_FAC*J_peak) {
    	//if ((J[i_adj][j_adj][k_adj] > JPEAK_FAC*J_peak) && (betapl[i_adj][j_adj][k_adj] > BETAPL_THR)) {
        //if ((J[i_adj][j_adj][k_adj] > JPEAK_FAC*J_peak) && (sigmaphi[i_adj][j_adj][k_adj] < SIGMAPHI_THR)) {
        //if ((J[i_adj][j_adj][k_adj] > JPEAK_FAC*J_peak) && (Sigmaphi[i_adj][j_adj][k_adj] < SIGMAPHI_THR)) {
        //aux = bu3[i_adj][j_adj][k_adj];
        //if ((J[i_adj][j_adj][k_adj] > 0.2*J_peak) && BETWEEN(aux, -BU3_THR, BU3_THR)) {
        	J_cs[i_adj][j_adj][k_adj] = J[i_adj][j_adj][k_adj];
			mm++;
			icoords[mm] = i_adj;
			jcoords[mm] = j_adj;
            if (N3 > 1) {
                kcoords[mm] = k_adj;
            }
			cells++;
			size_current++;
        }
        // if it's not large enough, flag adjacent cell
        else {
        	flag[i_adj][j_adj][k_adj] = 1;
            flag_buffer[i_adj][j_adj][k_adj] = 1;
        }
    }
}



void write_current_sheets(char *fname, double ***sheets) {

    // Writes J_cs on a file to be plotted by a Python script. J_cs will
    // be written as a 1D array of size N1*N2*N3. This means that we must use
    // numpy.reshape(J_cs, (N1, N2, N3)) to plot it with our Python script.
    //
    // I mean, I could've made it write the current as a matrix, but... lazy :)

    FILE *fp;
    int i, j, k;

    fp = fopen(fname, "wb");

    printf("Writing into file '%s'...\n", fname);
    for (i = 0; i < N1; i++) {
        for (j = 0; j < N2; j++) {
            for (k = 0; k < N3; k++) {
                fwrite(&sheets[i][j][k], sizeof(double), 1, fp);
            }
        }
    }

    fclose(fp);
}

/*
void read_current_sheets(char *fname)
{
    // Reads J_cs from a file generated by write_current_sheets.

    FILE *fp;
    int i, j, k;

    fp = fopen(fname, "r");

    printf("Reading current sheets from file '%s'...\n", fname);
    for (i = 0; i < N1; i++) {
        for (j = 0; j < N2; j++) {
            for (k = 0; k < N3; k++) {
                scanf(fp, "%lf\n", J_cs[i][j][k]);
            }
        }
    }
    printf("Done!\n");

    fclose(fp);
}
*/
