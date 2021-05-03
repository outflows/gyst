#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "definitions.h"
#include "declarations.h"


void reverse_array(double *array, int n)
{
    int i;
    int end = n - 1;
    double temp;

    for (i = 0; i < n/2; i++) {
        temp = array[i];
        array[i] = array[end];
        array[end] = temp;
        end--;
    }
}

void reverse_array_int(int *array, int n)
{
    int i;
    int end = n - 1;
    int temp;

    for (i = 0; i < n/2; i++) {
        temp = array[i];
        array[i] = array[end];
        array[end] = temp;
        end--;
    }
}

void merge_arrays(double *a1, double *a2, double *newa, int size_a1, int size_a2)
{
    int i, j;

    j = 0;
    for (i = 0; i < size_a1; i++) {
        newa[j] = a1[i];
        j++;
    }
    for (i = 0; i < size_a2; i++) {
        newa[j] = a2[i];
        j++;
    }
}

void merge_arrays_int(int *a1, int *a2, int *newa, int size_a1, int size_a2)
{
    int i, j;

    j = 0;
    for (i = 0; i < size_a1; i++) {
        newa[j] = a1[i];
        j++;
    }
    for (i = 0; i < size_a2; i++) {
        newa[j] = a2[i];
        j++;
    }
}

void characterize()
{
    int i, j, k, ii, jj, kk, m, iter;
    double rr, tth, pphi, p_r, p_th, p_phi, x, y, p_x, p_y;
    double Hess[2][2], e1[2], e1_cart[2], e2[2], evec[4];
    double ap_x[NDIM], p_cart[NDIM], Vprim[NDIM], Bprim[NDIM];
    double del[NDIM];
    double new_position[3];
    double V_proj_e1, B_proj_e1, B_proj_e2;
    double step, stepsize, B_upper, B_lower, V_upper, V_lower, V_A, V_rec;
    double Br_center, Bth_center, Bphi_center, beta_center, sigma_center;

    int is_good_count, is_good, is_good_lower, is_good_upper, size_lower, size_upper, size_B;
    char i_str[6], j_str[6], k_str[6], sheet_file[200], count_file[200];
    FILE *fp;

    double Br_upper_arr[100], Br_lower_arr[100], Bth_upper_arr[100], Bth_lower_arr[100], Bphi_upper_arr[100], Bphi_lower_arr[100];
    double J_upper_arr[100], J_lower_arr[100], beta_upper_arr[100], beta_lower_arr[100];
    double sigma_upper_arr[100], sigma_lower_arr[100], sigmaphi_upper_arr[100], sigmaphi_lower_arr[100];
    int i_upper_arr[100], i_lower_arr[100], j_upper_arr[100], j_lower_arr[100], k_upper_arr[100], k_lower_arr[100];
    double Br_arr[200], Bth_arr[200], Bphi_arr[200], J_arr[200], sigma_arr[200], sigmaphi_arr[200], beta_arr[200];
    int i_arr[200], j_arr[200], k_arr[200];

    int i_am, j_am, k_am;
    double B_test, B_min, B_max;
    int more_to_go_lower, more_to_go_upper;
    double iaux, jaux, kaux;

    if (SHEETS) fprintf(stdout, "\n");
    fprintf(stdout, "Beginning characterization...\n");

    fprintf(stdout, "Getting 1st derivatives...\n");
    for (i = 0; i < N1; i++) {
        for (j = 0; j < N2; j++) {
            for (k = 0; k < N3; k++) {
                get_1stderivative2D(i, j, k, J_cs);
                J_cs_char[i][j][k] = 0.;
            }
        }
    }

    is_good_count = 0;
    fprintf(stdout, "Entering main loop...\n");
    for (i = 0; i < ILIM; i++) {
        for (j = 0; j < N2; j++) {
            for (k = 0; k < N3; k++) {
                // disc: > BE_THR
                // jet,outflows: < BE_THR
                if ( (J_cs[i][j][k] > 0.) && (bernoulli[i][j][k] < BE_THR) && (va[i][j][k] > VA_THR) ) {
                    rr = a_r[i][j][k];
                    tth = a_th[i][j][k];
                    pphi = a_phi[i][j][k];

                    Br_center = B[1][i][j][k];
                    Bth_center = B[2][i][j][k];
                    Bphi_center = B[3][i][j][k];
                    beta_center = betapl[i][j][k];
                    sigma_center = sigma[i][j][k];

                    //1. get hessian (2nd derivative)
                    for (ii = 0; ii < 2; ii++) {
                        for (jj = 0; jj < 2; jj++) {
                            Hess[ii][jj] = 0.;
                        }
                    }
                    for (ii = 0; ii < 4; ii++) evec[ii] = 0.;
                    get_hessian2D(i, j, k, Hess);

                    //2. get eigenvectors (max: e1, min: e3)
                    get_evec2D(Hess, evec);
                    for (ii = 0; ii < 2; ii++) {
                        e1[ii] = evec[ii];
                        e2[ii] = evec[ii+2];
                    }

                    //3. move along eigenvector (upper)
                    i_am = i;
                    j_am = j;
                    k_am = k;

                    step = STEP1;
                    stepsize = step;
                    m = 0;
                    iter = 0;
                    is_good_upper = 0;
                    more_to_go_upper = 0;
                    B_upper = 0.;
                    V_upper = 0.;
                    size_upper = 0;
                    B_test = Bphi_center;

                    B_min = 0.;
                    B_max = 0.;

                    while(iter < 50) {
                        p_cart[0] = rr;
                        p_cart[1] = tth;

                        //5.1 find cartesian coordinates of point (find p_x, p_y, p_z)
                        for (ii = 0; ii < 2; ii++) new_position[ii] = p_cart[ii] + e1[ii]*step;
                        ap_x[1] = new_position[0];
                        ap_x[2] = new_position[1];
                        ap_x[3] = 0.;
                        step += stepsize;

                        //7. find closest cell corresponding to spherical
                        Xtoijk_new(ap_x, &ii, &jj, &kk, del);

                        //7.1 test value of J_cs
                        if ( (J[ii][jj][kk] > 0.5*J_cs_peak[i][j][k]) ) {
                            if ((ii != i_am) || (jj != j_am)) {
                                i_am = ii;
                                j_am = jj;

                                Vprim[0] = 0.;
                                Vprim[1] = V[1][ii][jj][kk];
                                Vprim[2] = V[2][ii][jj][kk];
                                Vprim[3] = V[3][ii][jj][kk];
                                Bprim[0] = 0.;
                                Bprim[1] = B[1][ii][jj][kk];
                                Bprim[2] = B[2][ii][jj][kk];
                                Bprim[3] = B[3][ii][jj][kk];

                                B_upper = Bprim[3];
                                Br_upper_arr[m] = Bprim[1];
                                Bth_upper_arr[m] = Bprim[2];
                                Bphi_upper_arr[m] = Bprim[3];
                                J_upper_arr[m] = J[ii][jj][kk];
                                beta_upper_arr[m] = betapl[ii][jj][kk];
                                sigma_upper_arr[m] = sigma[ii][jj][kk];
                                sigmaphi_upper_arr[m] = sigmaphi[ii][jj][kk];
                                i_upper_arr[m] = ii;
                                j_upper_arr[m] = jj;
                                k_upper_arr[m] = 0;

                                if (Bphi_upper_arr[m] < B_min) {
                                    B_min = Bphi_upper_arr[m];
                                }
                                if (Bphi_upper_arr[m] > B_max) {
                                    B_max = Bphi_upper_arr[m];
                                }

                                m++;
                                size_upper = m;
                                is_good_upper = 1;
                            }
                            iter++;
                        }
                        else break;
                    }

                    //5. move along eigenvector (lower)
                    i_am = i;
                    j_am = j;
                    k_am = k;

                    step = STEP1;
                    stepsize = step;
                    m = 0;
                    iter = 0;
                    is_good_lower = 0;
                    more_to_go_lower = 0;
                    B_lower = 0.;
                    V_lower = 0.;
                    size_lower = 0;
                    B_test = Bphi_center;
                    while(iter < 50) {
                        p_cart[0] = rr;
                        p_cart[1] = tth;

                        //5.1 find cartesian coordinates of point (find p_x, p_y)
                        for (ii = 0; ii < 2; ii++) new_position[ii] = p_cart[ii] - e1[ii]*step;
                        ap_x[1] = new_position[0];
                        ap_x[2] = new_position[1];
                        ap_x[3] = 0.;
                        step += stepsize;

                        //7. find closest cell corresponding to spherical using Xtoijk
                        Xtoijk_new(ap_x, &ii, &jj, &kk, del);

                        //7.1 test value of J_cs
                        if (J[ii][jj][kk] > 0.5*J_cs_peak[i][j][k]) {
                            if ((ii != i_am) || (jj != j_am)) {
                                i_am = ii;
                                j_am = jj;

                                Vprim[0] = 0.;
                                Vprim[1] = V[1][ii][jj][kk];
                                Vprim[2] = V[2][ii][jj][kk];
                                Vprim[3] = V[3][ii][jj][kk];
                                Bprim[0] = 0.;
                                Bprim[1] = B[1][ii][jj][kk];
                                Bprim[2] = B[2][ii][jj][kk];
                                Bprim[3] = B[3][ii][jj][kk];

                                B_lower = Bprim[3];
                                Br_lower_arr[m] = Bprim[1];
                                Bth_lower_arr[m] = Bprim[2];
                                Bphi_lower_arr[m] = Bprim[3];
                                J_lower_arr[m] = J[ii][jj][kk];
                                beta_lower_arr[m] = betapl[ii][jj][kk];
                                sigma_lower_arr[m] = sigma[ii][jj][kk];
                                sigmaphi_lower_arr[m] = sigmaphi[ii][jj][kk];
                                i_lower_arr[m] = ii;
                                j_lower_arr[m] = jj;
                                k_lower_arr[m] = 0;

                                if (Bphi_lower_arr[m] < B_min) {
                                    B_min = Bphi_lower_arr[m];
                                }
                                if (Bphi_lower_arr[m] > B_max) {
                                    B_max = Bphi_lower_arr[m];
                                }

                                m++;
                                size_lower = m;

                                is_good_lower = 1;
                                more_to_go_lower = 1;
                            }
                            iter++;
                        }
                        else break;
                    }

                    if (is_good_lower && is_good_upper) {
                        is_good = 1;

                        // don't consider point if too asymetrical
                        if (fabs(B_upper/B_lower) < 0.8) {
                            is_good = 0;
                        }
                        if (fabs(B_lower/B_upper) < 0.8) {
                            is_good = 0;
                        }

                        if (B_upper * B_lower > 0) {
                            is_good = 0;
                        }

                        size_B = size_lower + size_upper + 1;
                        if (size_B <= 3) {
                            is_good = 0;
                        }

                        if (is_good) {
                            is_good_count++;

                            for (m = 0; m < 200; m++) {
                                Br_arr[m] = 0.;
                                Bth_arr[m] = 0.;
                                Bphi_arr[m] = 0.;
                                J_arr[m] = 0.;
                                beta_arr[m] = 0.;
                                sigma_arr[m] = 0.;
                                sigmaphi_arr[m] = 0.;
                                i_arr[m] = 0;
                                j_arr[m] = 0;
                                k_arr[m] = 0;
                            }

                            reverse_array(Br_lower_arr, size_lower);
                            reverse_array(Bth_lower_arr, size_lower);
                            reverse_array(Bphi_lower_arr, size_lower);
                            reverse_array(J_lower_arr, size_lower);
                            reverse_array(beta_lower_arr, size_lower);
                            reverse_array(sigma_lower_arr, size_lower);
                            reverse_array(sigmaphi_lower_arr, size_lower);
                            reverse_array_int(i_lower_arr, size_lower);
                            reverse_array_int(j_lower_arr, size_lower);
                            reverse_array_int(k_lower_arr, size_lower);

                            Br_lower_arr[size_lower] = B[1][i][j][k];
                            Bth_lower_arr[size_lower] = B[2][i][j][k];
                            Bphi_lower_arr[size_lower] = B[3][i][j][k];
                            J_lower_arr[size_lower] = J[i][j][k];
                            beta_lower_arr[size_lower] = betapl[i][j][k];
                            sigma_lower_arr[size_lower] = sigma[i][j][k];
                            sigmaphi_lower_arr[size_lower] = sigmaphi[i][j][k];
                            i_lower_arr[size_lower] = i;
                            j_lower_arr[size_lower] = j;
                            k_lower_arr[size_lower] = 0;

                            merge_arrays(Br_lower_arr, Br_upper_arr, Br_arr, size_lower+1, size_upper);
                            merge_arrays(Bth_lower_arr, Bth_upper_arr, Bth_arr, size_lower+1, size_upper);
                            merge_arrays(Bphi_lower_arr, Bphi_upper_arr, Bphi_arr, size_lower+1, size_upper);
                            merge_arrays(J_lower_arr, J_upper_arr, J_arr, size_lower+1, size_upper);
                            merge_arrays(beta_lower_arr, beta_upper_arr, beta_arr, size_lower+1, size_upper);
                            merge_arrays(sigma_lower_arr, sigma_upper_arr, sigma_arr, size_lower+1, size_upper);
                            merge_arrays(sigmaphi_lower_arr, sigmaphi_upper_arr, sigmaphi_arr, size_lower+1, size_upper);
                            merge_arrays_int(i_lower_arr, i_upper_arr, i_arr, size_lower+1, size_upper);
                            merge_arrays_int(j_lower_arr, j_upper_arr, j_arr, size_lower+1, size_upper);
                            merge_arrays_int(k_lower_arr, k_upper_arr, k_arr, size_lower+1, size_upper);

                            // check absolute max and min
                            B_max = 0.;
                            B_min = 0.;
                            for (m = 0; m < size_B; m++) {
                                if (Bphi_arr[m] < B_min) {
                                    B_min = Bphi_arr[m];
                                }
                                if (Bphi_arr[m] > B_max) {
                                    B_max = Bphi_arr[m];
                                }
                            }

                            if ( (fabs(B_max/B_min) < 0.8) || (fabs(B_min/B_max) < 0.8) ) {
                                is_good = 0;
                                is_good_count--;
                            }

                            if ( (RUN_JET) && (is_good != 0) ) {
                                if ( (sigma_arr[0] < 1.) && (sigma_arr[size_lower+size_upper] < 1.) ) {
                                    is_good = 0;
                                    is_good_count--;
                                }
                            }

                            if (is_good) {
                                J_cs_char[i][j][k] = J[i][j][k];
                                sprintf(i_str, "%d", i);
                                sprintf(j_str, "%d", j);
                                sprintf(k_str, "%d", k);
                                #if (RUN_DISC)
                                {
                                    //strcpy(sheet_file, "/work/gustavo/gyst/sheets_disc/"); // DESKTOP IAG
                                    strcpy(sheet_file, SHEETS_DISC_DIR);
                                }
                                #elif (RUN_OUTFLOWS)
                                {
                                    //strcpy(sheet_file, "/work/gustavo/gyst/sheets_outflows/"); // DESKTOP IAG
                                    strcpy(sheet_file, SHEETS_OUTFLOWS_DIR);
                                }
                                #elif (RUN_JET)
                                {
                                    //strcpy(sheet_file, "/work/gustavo/gyst/sheets_jet/"); // DESKTOP IAG
                                    strcpy(sheet_file, SHEETS_JET_DIR);
                                }
                                #endif
                                strcat(sheet_file, dump_file);
                                strcat(sheet_file, "_");
                                strcat(sheet_file, i_str);
                                strcat(sheet_file, "_");
                                strcat(sheet_file, j_str);
                                strcat(sheet_file, "_");
                                strcat(sheet_file, k_str);
                                strcat(sheet_file, "_s");

                                fp = fopen(sheet_file, "wb");
                                for (m = 0; m < size_B; m++) {
                                    fwrite(&Br_arr[m], sizeof(double), 1, fp);
                                    fwrite(&Bth_arr[m], sizeof(double), 1, fp);
                                    fwrite(&Bphi_arr[m], sizeof(double), 1, fp);
                                    fwrite(&J_arr[m], sizeof(double), 1, fp);
                                    fwrite(&beta_arr[m], sizeof(double), 1, fp);
                                    fwrite(&sigma_arr[m], sizeof(double), 1, fp);
                                    fwrite(&sigmaphi_arr[m], sizeof(double), 1, fp);
                                    iaux = i_arr[m];
                                    jaux = j_arr[m];
                                    kaux = 0.;
                                    fwrite(&iaux, sizeof(double), 1, fp);
                                    fwrite(&jaux, sizeof(double), 1, fp);
                                    fwrite(&kaux, sizeof(double), 1, fp);
                                    //fwrite(&V_rec, sizeof(double), 1, fp);
                                }
                                fclose(fp);
                            }
                        }
                    }
                }
                //else (J_cs_char[i][j][k] = 0.);
            }
        }
    }
    fprintf(stdout, "Identified %d possible reconnection sites.\n", is_good_count);
    FILE *sitecount;
    #if (RUN_DISC)
    {
        strcpy(count_file, RUN_DIR);
        strcat(count_file, "sitecount_2D_disc.dat");
        sitecount = fopen(count_file, "a"); //// DONT FORGET TO CHANGE!!!!!!
    }
    #elif (RUN_OUTFLOWS)
    {
        strcpy(count_file, RUN_DIR);
        strcat(count_file, "sitecount_2D_outflows.dat");
        sitecount = fopen(count_file, "a"); //// DONT FORGET TO CHANGE!!!!!!
    }
    #elif (RUN_JET)
    {
        strcpy(count_file, RUN_DIR);
        strcat(count_file, "sitecount_2D_jet.dat");
        sitecount = fopen(count_file, "a"); //// DONT FORGET TO CHANGE!!!!!!
    }
    #endif
    fprintf(sitecount, "%s,%d\n", dump_file, is_good_count);
    fclose(sitecount);
    write_current_sheets(jchar_output, J_cs_char);
    fprintf(stdout, "Finished characterization of sheets.\n");
}
