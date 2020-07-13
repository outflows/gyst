#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "definitions.h"
#include "declarations.h"

void characterize()
{
    int i, j, k, m;
    int ii, jj, kk, newk;
    double rr, tth, pphi, p_r, p_th, p_phi;
    double x, y, z, p_x, p_y, p_z;
    double Hess[3][3];
    double e1[3], e1_cart[3], e2[3], e3[3], evec[9];
    double ap_x[NDIM], p_cart[NDIM], Vprim[NDIM], Bprim[NDIM];
    double del[NDIM];
    double new_position[3];
    double V_proj_e1, B_proj_e1, B_proj_e2, B_proj_e3;
    double step, stepsize, B_upper, B_lower, V_upper, V_lower, V_A, V_rec;
    double B_upper_arr[100], B_lower_arr[100];
    double *B_arr;
    int is_good, is_good_lower, is_good_upper, size_lower, size_upper, size_B;
    char sheet_file[50], i_str[6], j_str[6], k_str[6];
    FILE *fp;

    if (SHEETS) printf("\n");
    printf("Beginning characterization...\n");

    printf("Getting 1st derivatives...\n");
    for (i = 0; i < N1; i++) {
        for (j = 0; j < N2; j++) {
            for (k = 0; k < N3; k++) {
                get_1stderivative(i, j, k, J_cs);
                J_cs_char[i][j][k] = J_cs_peak[i][j][k];
            }
        }
    }

    printf("Entering main loop...\n");
    for (i = 0; i < N1; i++) {
        for (j = 0; j < N2; j++) {
            for (k = 0; k < N3; k++) {

                if (J_cs_peak[i][j][k] > 0.) {
                    rr = a_r[i][j][k];
                    tth = a_th[i][j][k];
                    pphi = a_phi[i][j][k];

                    //1. get hessian
                    for (ii = 0; ii < 3; ii++) {
                        for (jj = 0; jj < 3; jj++) {
                            Hess[ii][jj] = 0.;
                        }
                    }
                    for (ii = 0; ii < 9; ii++) evec[ii] = 0.;
                    get_hessian(i, j, k, Hess);

                    //2. get eigenvectors (max: e1, min: e3)
                    get_evec(Hess, evec);
                    for (ii = 0; ii < 3; ii++) {
                        e1[ii] = evec[ii];
                        e2[ii] = evec[ii+3];
                        e3[ii] = evec[ii+6];
                    }

                    //3. convert r, theta, phi to cartesian
                    rthphi_to_xyz(rr, tth, pphi, &x, &y, &z);

                    //4. convert max_eigvec to cartesian
                    vec_sph_to_cart(e1, e1_cart, tth, pphi);

                    //5. move along eigenvector (upper)
                    step = STEP;
                    stepsize = step;
                    m = 0;
                    size_upper = 1;
                    //B_upper_arr = (double*) malloc(size_upper*sizeof(double));
                    for (int index = 0; index < 100; index++) B_upper_arr[index] = 0.;
                    is_good_upper = 0;
                    while(m < 100) {
                        p_cart[0] = x;
                        p_cart[1] = y;
                        p_cart[2] = z;

                        //5.1 find cartesian coordinates of point (find p_x, p_y, p_z)
                        for (ii = 0; ii < 3; ii++) new_position[ii] = p_cart[ii] + e1_cart[ii]*step;
                        printf("%lf %lf %lf\n", e1_cart[0], e1_cart[1], e1_cart[2]);
                        printf("%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", e1[0], e1[1], e1[2], e1[3], e1[4], e1[5], e1[6], e1[7], e1[8]);
                        p_x = new_position[0];
                        p_y = new_position[1];
                        p_z = new_position[2];
                        step += stepsize;

                        //5.2 convert point to spherical
                        xyz_to_rthphi(p_x, p_y, p_z, &p_r, &p_th, &p_phi);

                        //6. find x1, x2, x3 corresponding to them
                        ap_x[1] = zbrent(find_x1_cyl, p_r, 0, 8., 1E-6);
                        x1in = ap_x[1];
                        ap_x[2] = zbrent(find_x2_cyl, p_th, -1.0, 1.0, 1E-6); // check this
                        ap_x[3] = p_phi;

                        //7. find closest cell corresponding to spherical
                        Xtoijk(ap_x, &ii, &jj, &kk, del);

                        //printf("upper %g %g %g %d %d %d %g\n", ap_x[1], ap_x[2], ap_x[3], ii, jj, kk, J_cs[ii][jj][kk]);

                        //7.1 test value of J_cs
                        if (J_cs[ii][jj][kk] > 0.5*J_cs_peak[i][j][k]) {
                            printf("size upper %d, %d %d %d %d %g %lf\n", size_upper, i, j, size_upper*sizeof(double), m, B_upper_arr[m], p_x);
                            //realloc(B_upper_arr, size_upper*sizeof(double));
                            size_upper++;

                            //8. project V along max_evec to find V_in at this cell
                            // NOTE: may have to change, use 4vec and recalculate V, B
                            Vprim[0] = 0.;
                            Vprim[1] = V[1][ii][jj][kk];
                            Vprim[2] = V[2][ii][jj][kk];
                            Vprim[3] = V[3][ii][jj][kk];
                            Bprim[0] = 0.;
                            Bprim[1] = B[1][ii][jj][kk];
                            Bprim[2] = B[2][ii][jj][kk];
                            Bprim[3] = B[3][ii][jj][kk];
                            V_proj_e1 = dot3(Vprim, e1) / (sqrt(dot3(e1,e1)) + SMALL);

                            //9. project B along e1, e2, e3 to find B along these eigenvectors
                            B_proj_e1 = dot3(Bprim, e1) / (sqrt(dot3(e1,e1)) + SMALL);
                            B_proj_e2 = dot3(Bprim, e2) / (sqrt(dot3(e2,e2)) + SMALL);
                            B_proj_e3 = dot3(Bprim, e3) / (sqrt(dot3(e3,e3)) + SMALL);
                            B_upper = abs(B_proj_e1);
                            B_upper_arr[m] = B_proj_e1;
                            //printf("%d %lf\n", m, B_upper_arr[m]);
                            m++;

                            //10. use B_proj to find v_alfven
                            V_A = sqrt((B_proj_e1*B_proj_e1 +
                                        B_proj_e2*B_proj_e2 +
                                        B_proj_e3*B_proj_e3))/(rho[ii][jj][kk] + SMALL);
                            V_upper = V_proj_e1/V_A;
                            is_good_upper = 1;
                        }
                        else break;
                    }

                    //5. move along eigenvector (lower)
                    step = STEP;
                    stepsize = step;
                    m = 0;
                    size_lower = 1;
                    //B_lower_arr = (double*) malloc(size_lower*sizeof(double));
                    for (int index = 0; index < 100; index++) B_lower_arr[index] = 0.;
                    is_good_lower = 0;
                    while(m < 100) {
                        p_cart[0] = x;
                        p_cart[1] = y;
                        p_cart[2] = z;

                        //5.1 find cartesian coordinates of point (find p_x, p_y, p_z)
                        for (ii = 0; ii < 3; ii++) new_position[ii] = p_cart[ii] - e1_cart[ii]*step;
                        p_x = new_position[0];
                        p_y = new_position[1];
                        p_z = new_position[2];
                        step += stepsize;

                        //5.2 convert point to spherical
                        xyz_to_rthphi(p_x, p_y, p_z, &p_r, &p_th, &p_phi);

                        //6. find x1, x2, x3 corresponding to them
                        ap_x[1] = zbrent(find_x1_cyl, p_r, 0, 8., 1E-6);
                        x1in = ap_x[1];
                        ap_x[2] = zbrent(find_x2_cyl, p_th, -1.0, 1.0, 1E-6);
                        ap_x[3] = p_phi;

                        //7. find closest cell corresponding to spherical using Xtoijk
                        Xtoijk(ap_x, &ii, &jj, &kk, del);

                        //printf("lower %g %g %g %d %d %d %g\n", ap_x[1], ap_x[2], ap_x[3], ii, jj, kk, J_cs[ii][jj][kk]);

                        //7.1 test value of J_cs
                        if (J_cs[ii][jj][kk] > 0.5*J_cs_peak[i][j][k]) {
                            //printf("size lower %d, %d %d %d\n", size_lower, i, j, size_lower*sizeof(double));
                            //realloc(B_lower_arr, size_lower*sizeof(double));
                            size_lower++;

                            //8. project V along max_evec to find V_in at this cell
                            // NOTE: may have to change, use 4vec and recalculate V, B
                            Vprim[0] = 0.;
                            Vprim[1] = V[1][ii][jj][kk];
                            Vprim[2] = V[2][ii][jj][kk];
                            Vprim[3] = V[3][ii][jj][kk];
                            Bprim[0] = 0.;
                            Bprim[1] = B[1][ii][jj][kk];
                            Bprim[2] = B[2][ii][jj][kk];
                            Bprim[3] = B[3][ii][jj][kk];
                            V_proj_e1 = dot3(Vprim, e1) / (sqrt(dot3(e1,e1)) + SMALL);

                            //9. project B along e1, e2, e3 to find B along these eigenvectors
                            B_proj_e1 = dot3(Bprim, e1) / (sqrt(dot3(e1,e1)) + SMALL);
                            B_proj_e2 = dot3(Bprim, e2) / (sqrt(dot3(e2,e2)) + SMALL);
                            B_proj_e3 = dot3(Bprim, e3) / (sqrt(dot3(e3,e3)) + SMALL);
                            B_lower = abs(B_proj_e1);
                            B_lower_arr[m] = B_proj_e1;
                            m++;

                            //10. use B_proj to find v_alfven
                            V_A = sqrt((B_proj_e1*B_proj_e1 +
                                        B_proj_e2*B_proj_e2 +
                                        B_proj_e3*B_proj_e3))/(rho[ii][jj][kk] + SMALL);
                            V_lower = V_proj_e1/V_A;
                            is_good_lower = 1;
                        }
                        else break;
                    }

                    //printf("%d %d\n", is_good_lower, is_good_upper);
                    if (is_good_lower && is_good_upper) {
                        //printf("both good\n");
                        //12. get V_in/V_A
                        V_rec = 0.5*(V_lower - V_upper);

                        is_good = 1;
                        // don't consider point if smallish reconnection velocity
                        if (V_rec < 0.7) {
                            J_cs_char[i][j][k] = 0.;
                            is_good = 0;
                        }
                        // don't consider point if too asymetrical
                        if (B_upper < 0.75*B_lower || B_lower < 0.75*B_upper) {
                            J_cs_char[i][j][k] = 0.;
                            is_good = 0;
                        }
                        if (is_good) {
                            //printf("IT'S GOOD\n");
                            size_B = size_lower + size_upper - 1;
                            B_arr = (double*) malloc(size_B*sizeof(double));
                            reverse_array(B_lower_arr, size_lower);
                            merge_arrays(B_lower_arr, B_upper_arr, B_arr, size_lower - 1, size_upper);
                            sprintf(i_str, "%d", i);
                            sprintf(j_str, "%d", j);
                            sprintf(k_str, "%d", k);

                            strcpy(sheet_file, dump_file);
                            strcat(sheet_file, "_");
                            strcat(sheet_file, i_str);
                            strcat(sheet_file, "_");
                            strcat(sheet_file, j_str);
                            strcat(sheet_file, "_");
                            strcat(sheet_file, k_str);
                            strcat(sheet_file, "_s");

                            fp = fopen(sheet_file, "w");
                            for (m = 0; m < size_B; m++) {
                                fprintf(fp, "%g\n", B_arr[m]);
                            }
                            fclose(fp);
                        }
                        //free(B_upper_arr);
                        //free(B_lower_arr);
                        free(B_arr);
                    }
                }
            }
        }
    }
    write_current_sheets(jchar_output, J_cs_char);
    printf("Finished characterization of sheets.\n");

}


void xyz_to_rthphi(double x, double y, double z, double *r, double *th, double *phi)
{
    // change to BL?
    *r = sqrt(x*x + y*y + z*z);
    *th = atan(sqrt(x*x + y*y)/z);
    *phi = atan(y/x);
    //printf("during %g %g %g %g %g %g\n", x, y, z, *r, *th, *phi);
}

void rthphi_to_xyz(double r, double th, double phi, double *x, double *y, double *z)
{
    // change to BL?
    *x = r*sin(th)*cos(phi);
    *y = r*sin(th)*sin(phi);
    *z = r*cos(th);
    //if (N3 == 1) y = 0.;
    //printf("during %g %g %g %g %g %g\n", r, th, phi, *x, *y, *z);
}

void vec_sph_to_cart(double eig_sph[3], double eig_cart[3], double theta, double phi)
{
    int i;
    double sph_to_cart[3][3];

    sph_to_cart[0][0] = sin(theta)*cos(phi);
    sph_to_cart[0][1] = cos(theta)*cos(phi);
    sph_to_cart[0][2] = -sin(phi);
    sph_to_cart[1][0] = sin(theta)*sin(phi);
    sph_to_cart[1][1] = cos(theta)*sin(phi);
    sph_to_cart[1][2] = cos(phi);
    sph_to_cart[2][0] = cos(theta);
    sph_to_cart[2][1] = -sin(phi);
    sph_to_cart[2][2] = 0.;

    for (i = 0; i < 3; i++) {
        eig_cart[i] = 0.;
    }

    for (i = 0; i < 3; i++) {
        eig_cart[0] += eig_sph[i]*sph_to_cart[i][0];
        eig_cart[1] += eig_sph[i]*sph_to_cart[i][1];
        eig_cart[2] += eig_sph[i]*sph_to_cart[i][2];
    }
}


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

/*

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
                    check_box_limits(i, j, k,
                        sheet_box_lower_i, sheet_box_lower_j, sheet_box_lower_k,
                        sheet_box_upper_i, sheet_box_upper_j, sheet_box_upper_k,
                        halfboxsize);

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

void check_box_limits(int i, int j, int k,
                      int box_lower_i, int box_lower_j, int box_lower_k,
                      int box_upper_i, int box_upper_j, int box_upper_k,
                      int halfboxsize) {

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



// FIND points at same normal

// get b at cell (r,z), (interpolate)
// increase distances
// get b at new cells (interpolate)

r_plus = r_c + d_plus/(sqrt(1. + m_perp*m_perp));
r_minus = r_c - d_minus/(sqrt(1. + m_perp*m_perp));
z_plus = z_c + m_perp*(r_plus - r_c);
z_minus = z_c + m_perp*(r_minus - r_c);

*/












/*



theta = th[i][j][k];
phi = phi[i][j][k];
void sph_to_cart2D(double eig_sph[2], double eig_cart[2], double theta)
{
    int i;
    double sph_to_cart[2][2];

    sph_to_cart[0][0] = sin(theta)*cos(phi);
    sph_to_cart[0][1] = cos(theta)*cos(phi);
    sph_to_cart[1][0] = sin(theta)*sin(phi);
    sph_to_cart[1][1] = cos(theta)*sin(phi);


    for (i = 0; i < 2; i++) {
        eig_cart[i] = 0.;
    }

    for (i = 0; i < 2; i++) {
        eig_cart[0] += eig_sph[i]*sph_to_cart[i][0];
        eig_cart[1] += eig_sph[i]*sph_to_cart[i][1];
    }
}


theta = th[i][j][k];
phi = phi[i][j][k];
void cart_to_sph2D(double eig_cart[2], double eig_sph[2], double theta)
{
    int i;
    double cart_to_sph[2][2];

    cart_to_sph[0][0] = sin(theta)*cos(phi);
    cart_to_sph[0][1] = sin(theta)*sin(phi);
    cart_to_sph[1][0] = cos(theta)*cos(phi);
    cart_to_sph[1][1] = cos(theta)*sin(phi);

    // transform from cart to sph
    for (i = 0; i < 2; i++) {
        eig_sph[i] = 0;
    }

    for (i = 0; i < 2; i++) {
        eig_sph[0] += eig_cart[i]*cart_to_sph[i][0];
        eig_sph[1] += eig_cart[i]*cart_to_sph[i][1];
    }
}


theta = th[i][j][k];
phi = phi[i][j][k];
void sph_to_cart3D(double eig_sph[3], double eig_cart[3], double theta, double phi)
{
    int i;
    double sph_to_cart[3][3];

    sph_to_cart[0][0] = sin(theta)*cos(phi);
    sph_to_cart[0][1] = cos(theta)*cos(phi);
    sph_to_cart[0][2] = -sin(phi);
    sph_to_cart[1][0] = sin(theta)*sin(phi);
    sph_to_cart[1][1] = cos(theta)*sin(phi);
    sph_to_cart[1][2] = cos(phi);
    sph_to_cart[2][0] = cos(theta);
    sph_to_cart[2][1] = -sin(phi);
    sph_to_cart[2][2] = 0.;

    for (i = 0; i < 3; i++) {
        eig_cart[i] = 0.;
    }

    for (i = 0; i < 3; i++) {
        eig_cart[0] += eig_sph[i]*sph_to_cart[i][0];
        eig_cart[1] += eig_sph[i]*sph_to_cart[i][1];
        eig_cart[2] += eig_sph[i]*sph_to_cart[i][2];
    }
}


theta = th[i][j][k];
phi = phi[i][j][k];
void cart_to_sph3D(double eig_cart[3], double eig_sph[3], double theta, double phi)
{
    int i;
    double cart_to_sph[3][3];

    cart_to_sph[0][0] = sin(theta)*cos(phi);
    cart_to_sph[0][1] = sin(theta)*sin(phi);
    cart_to_sph[0][2] = cos(phi);
    cart_to_sph[1][0] = cos(theta)*cos(phi);
    cart_to_sph[1][1] = cos(theta)*sin(phi);
    cart_to_sph[1][2] = -sin(phi);
    cart_to_sph[2][0] = -sin(phi);
    cart_to_sph[2][1] = cos(phi);
    cart_to_sph[2][2] = 0.;

    // transform from cart to sph
    for (i = 0; i < 3; i++) {
        eig_sph[i] = 0;
    }

    for (i = 0; i < 3; i++) {
        eig_sph[0] += eig_cart[i]*cart_to_sph[i][0];
        eig_sph[1] += eig_cart[i]*cart_to_sph[i][1];
        eig_sph[2] += eig_cart[i]*cart_to_sph[i][2];
    }
}


theta = th[i][j][k];
phi = phi[i][j][k];
void cart_to_sph(double eig_cart[3], double eig_sph[3], double theta, double phi)
{
    int i;
    double cart_to_sph[3][3];

    cart_to_sph[0][0] = sin(theta)*cos(phi);
    cart_to_sph[0][1] = sin(theta)*sin(phi);
    cart_to_sph[0][2] = cos(phi);
    cart_to_sph[1][0] = cos(theta)*cos(phi);
    cart_to_sph[1][1] = cos(theta)*sin(phi);
    cart_to_sph[1][2] = -sin(phi);
    cart_to_sph[2][0] = -sin(phi);
    cart_to_sph[2][1] = cos(phi);
    cart_to_sph[2][2] = 0.;

    // transform from cart to sph
    for (i = 0; i < 3; i++) {
        eig_sph[i] = 0;
    }

    for (i = 0; i < 3; i++) {
        eig_sph[0] += eig_cart[i]*cart_to_sph[i][0];
        eig_sph[1] += eig_cart[i]*cart_to_sph[i][1];
        eig_sph[2] += eig_cart[i]*cart_to_sph[i][2];
    }
}

*/












void characterize2D()
{
    int i, j, k, ii, jj, kk, m;
    double rr, tth, pphi, p_r, p_th, p_phi, x, y, p_x, p_y;
    double Hess[2][2], e1[2], e1_cart[2], e2[2], evec[4];
    double ap_x[NDIM], p_cart[NDIM], Vprim[NDIM], Bprim[NDIM];
    double del[NDIM];
    double new_position[3];
    double V_proj_e1, B_proj_e1, B_proj_e2;
    double step, stepsize, B_upper, B_lower, V_upper, V_lower, V_A, V_rec;
    double Br_center, Bth_center, Bphi_center, beta_center, sigma_center;

    int is_good_count, is_good, is_good_lower, is_good_upper, size_lower, size_upper, size_B;
    char sheet_file[100], i_str[6], j_str[6], k_str[6];
    FILE *fp;

//    double B_upper_arr[100], B_lower_arr[100], J_upper_arr[100], J_lower_arr[100];
//    double *B_arr, *J_arr;

    double Br_upper_arr[100], Br_lower_arr[100], Bth_upper_arr[100], Bth_lower_arr[100], Bphi_upper_arr[100], Bphi_lower_arr[100];
    double J_upper_arr[100], J_lower_arr[100], beta_upper_arr[100], beta_lower_arr[100];
    double sigma_upper_arr[100], sigma_lower_arr[100], sigmaphi_upper_arr[100], sigmaphi_lower_arr[100];
    int i_upper_arr[100], i_lower_arr[100], j_upper_arr[100], j_lower_arr[100];

    double Br_arr[200], Bth_arr[200], Bphi_arr[200], J_arr[200], sigma_arr[200], sigmaphi_arr[200], beta_arr[200];
    int i_arr[200], j_arr[200];

    if (SHEETS) printf("\n");
    printf("Beginning characterization...\n");

    printf("Getting 1st derivatives...\n");
    for (i = 0; i < N1; i++) {
        for (j = 0; j < N2; j++) {
            for (k = 0; k < N3; k++) {
                get_1stderivative2D(i, j, k, J_cs);
                J_cs_char[i][j][k] = J_cs_peak[i][j][k];
            }
        }
    }

    is_good_count = 0;
    printf("Entering main loop...\n");
    for (i = 0; i < N1; i++) {
        for (j = 0; j < N2; j++) {
            for (k = 0; k < N3; k++) {

                //if(J_cs_peak[i][j][k] <= 0) printf("jcs %d %d %d %lf\n", i, j, k, J_cs[i][j][k]);
                if (J_cs_peak[i][j][k] > 0.) {
                    //printf("%d %d\n", i, j);
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

                    //3. convert r, theta, to cartesian
                    rth_to_xy(rr, tth, &x, &y);

                    //4. convert max_eigvec to cartesian
                    vec_pol_to_cart(e1, e1_cart, tth);

                    //printf("%d %d %d %lf %lf\n", i, j, k, B[3][i][j][k], J[i][j][k]);
                    //5. move along eigenvector (upper)
                    step = STEP;
                    stepsize = step;
                    m = 0;
                    //for (int index = 0; index < 100; index++) B_upper_arr[index] = 0.;
                    is_good_upper = 0;
                    while(m < 80) {
                        p_cart[0] = x;
                        p_cart[1] = y;

                        //5.1 find cartesian coordinates of point (find p_x, p_y, p_z)
                        for (ii = 0; ii < 2; ii++) new_position[ii] = p_cart[ii] + e1_cart[ii]*step;
                        p_x = new_position[0];
                        p_y = new_position[1];
                        step += stepsize;

                        //printf("p_cart %lf %lf\n", p_cart[0], p_cart[1]);
                        //printf("e1 %lf %lf\n", e1[0], e1[1]);
                        //printf("e1_cart %lf %lf\n", e1_cart[0], e1_cart[1]);
                        //printf("new_p %lf %lf\n", p_x, p_y);

                        //5.2 convert point to polar
                        xy_to_rth(p_x, p_y, &p_r, &p_th);

                        //6. find x1, x2, x3 corresponding to them
                        ap_x[1] = zbrent(find_x1_cyl, p_r, 0, 8., 1E-6);
                        x1in = ap_x[1];
                        ap_x[2] = zbrent(find_x2_cyl, p_th, -1.0, 1.0, 1E-6); // check this
                        ap_x[3] = p_phi;

                        //7. find closest cell corresponding to spherical
                        Xtoijk(ap_x, &ii, &jj, &kk, del);

                        //7.1 test value of J_cs
                        if (J[ii][jj][kk] > 0.5*J_cs_peak[i][j][k]) {

                            Vprim[0] = 0.;
                            Vprim[1] = V[1][ii][jj][kk];
                            Vprim[2] = V[2][ii][jj][kk];
                            Vprim[3] = V[3][ii][jj][kk];
                            Bprim[0] = 0.;
                            Bprim[1] = B[1][ii][jj][kk];
                            Bprim[2] = B[2][ii][jj][kk];
                            Bprim[3] = B[3][ii][jj][kk];

                            //9. project V, B along e1, e2, e3 to find V, B along these eigenvectors
                            // NOTE: may have to change, use 4vec and recalculate V, B
                            /*
                            V_proj_e1 = dot3(Vprim, e1) / (sqrt(dot3(e1,e1)) + SMALL);
                            B_proj_e1 = dot3(Bprim, e1) / (sqrt(dot3(e1,e1)) + SMALL);
                            B_proj_e2 = dot3(Bprim, e2) / (sqrt(dot3(e2,e2)) + SMALL);
                            B_upper = abs(B_proj_e1);
                            B_upper_arr[m] = B_proj_e1;
                            */

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

                            m++;
                            size_upper = m;

                            //10. use B_proj to find v_alfven
                            /*
                            V_A = sqrt((B_proj_e1*B_proj_e1 +
                                        B_proj_e2*B_proj_e2)/(rho[ii][jj][kk] + SMALL));
                            V_upper = V_proj_e1/V_A;
                            */

                            V_A = sqrt((Bprim[1]*Bprim[1] +
                                        Bprim[2]*Bprim[2] +
                                        Bprim[3]*Bprim[3])/(rho[ii][jj][kk] + SMALL));
                            V_upper = Vprim[1]/V_A;

                            is_good_upper = 1;
                        }
                        else break;
                    }

                    //5. move along eigenvector (lower)
                    step = STEP;
                    stepsize = step;
                    m = 0;
                    //for (int index = 0; index < 100; index++) B_lower_arr[index] = 0.;
                    is_good_lower = 0;
                    while(m < 80) {
                        p_cart[0] = x;
                        p_cart[1] = y;

                        //5.1 find cartesian coordinates of point (find p_x, p_y)
                        for (ii = 0; ii < 2; ii++) new_position[ii] = p_cart[ii] - e1_cart[ii]*step;
                        p_x = new_position[0];
                        p_y = new_position[1];
                        step += stepsize;

                        //5.2 convert point to spherical
                        xy_to_rth(p_x, p_y, &p_r, &p_th);

                        //6. find x1, x2, x3 corresponding to them
                        ap_x[1] = zbrent(find_x1_cyl, p_r, 0, 8., 1E-6);
                        x1in = ap_x[1];
                        ap_x[2] = zbrent(find_x2_cyl, p_th, -1.0, 1.0, 1E-6);
                        ap_x[3] = p_phi;

                        //7. find closest cell corresponding to spherical using Xtoijk
                        Xtoijk(ap_x, &ii, &jj, &kk, del);

                        //7.1 test value of J_cs
                        if (J[ii][jj][kk] > 0.5*J_cs_peak[i][j][k]) {

                            //8. project V along max_evec to find V_in at this cell
                            // NOTE: may have to change, use 4vec and recalculate V, B
                            Vprim[0] = 0.;
                            Vprim[1] = V[1][ii][jj][kk];
                            Vprim[2] = V[2][ii][jj][kk];
                            Vprim[3] = V[3][ii][jj][kk];
                            Bprim[0] = 0.;
                            Bprim[1] = B[1][ii][jj][kk];
                            Bprim[2] = B[2][ii][jj][kk];
                            Bprim[3] = B[3][ii][jj][kk];

                            //9. project B along e1, e2, e3 to find B along these eigenvectors
                            /*
                            V_proj_e1 = dot3(Vprim, e1) / (sqrt(dot3(e1,e1)) + SMALL);
                            B_proj_e1 = dot3(Bprim, e1) / (sqrt(dot3(e1,e1)) + SMALL);
                            B_proj_e2 = dot3(Bprim, e2) / (sqrt(dot3(e2,e2)) + SMALL);
                            B_lower = abs(B_proj_e1);
                            B_lower_arr[m] = B_proj_e1;
                            */

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

                            m++;
                            size_lower = m;

                            //10. use B_proj to find v_alfven
                            /*
                            V_A = sqrt((B_proj_e1*B_proj_e1 +
                                        B_proj_e2*B_proj_e2)/(rho[ii][jj][kk] + SMALL));
                            V_lower = V_proj_e1/V_A;
                            */

                            V_A = sqrt((Bprim[1]*Bprim[1] +
                                        Bprim[2]*Bprim[2] +
                                        Bprim[3]*Bprim[3])/(rho[ii][jj][kk] + SMALL));
                            V_lower = Vprim[1]/V_A;

                            is_good_lower = 1;
                        }
                        else break;
                    }

                    if (is_good_lower && is_good_upper) {
                        //12. get V_in/V_A
                        V_rec = 0.5*(V_lower - V_upper);

                        is_good = 1;
                        // don't consider point if smallish reconnection velocity
                        if (abs(V_rec) < 0.6) {
                            J_cs_char[i][j][k] = 0.;
                            is_good = 0;
                        }

                        // don't consider point if too asymetrical
                        if (abs(B_upper) < 0.9*abs(B_lower) || abs(B_lower) < 0.9*abs(B_upper)) {
                            J_cs_char[i][j][k] = 0.;
                            is_good = 0;
                        }

                        if (B_upper * B_lower > 0) {
                            J_cs_char[i][j][k] = 0.;
                            is_good = 0;
                        }

                        //if (is_good) printf("is good %d %d %d\n", i, j, k);

                        if (is_good) {
                            is_good_count++;
                            size_B = size_lower + size_upper - 2;

                            //B_arr = (double*) malloc(size_B*sizeof(double));
                            //J_arr = (double*) malloc(size_B*sizeof(double));
                            //reverse_array(B_lower_arr, size_lower);
                            //merge_arrays(B_lower_arr, B_upper_arr, B_arr, size_lower, size_upper);
                            //reverse_array(J_lower_arr, size_lower);
                            //merge_arrays(J_lower_arr, J_upper_arr, J_arr, size_lower, size_upper);

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
                            }

                            reverse_array(Br_lower_arr, size_lower);
                            merge_arrays(Br_lower_arr, Br_upper_arr, Br_arr, size_lower, size_upper);
                            reverse_array(Bth_lower_arr, size_lower);
                            merge_arrays(Bth_lower_arr, Bth_upper_arr, Bth_arr, size_lower, size_upper);
                            reverse_array(Bphi_lower_arr, size_lower);
                            merge_arrays(Bphi_lower_arr, Bphi_upper_arr, Bphi_arr, size_lower, size_upper);
                            reverse_array(J_lower_arr, size_lower);
                            merge_arrays(J_lower_arr, J_upper_arr, J_arr, size_lower, size_upper);
                            reverse_array(beta_lower_arr, size_lower);
                            merge_arrays(beta_lower_arr, beta_upper_arr, beta_arr, size_lower, size_upper);
                            reverse_array(sigma_lower_arr, size_lower);
                            merge_arrays(sigma_lower_arr, sigma_upper_arr, sigma_arr, size_lower, size_upper);
                            reverse_array(sigmaphi_lower_arr, size_lower);
                            merge_arrays(sigmaphi_lower_arr, sigmaphi_lower_arr, sigmaphi_arr, size_lower, size_upper);
                            reverse_array_int(i_lower_arr, size_lower);
                            merge_arrays_int(i_lower_arr, i_lower_arr, i_arr, size_lower, size_upper);
                            reverse_array_int(j_lower_arr, size_lower);
                            merge_arrays_int(j_lower_arr, j_lower_arr, j_arr, size_lower, size_upper);

                            sprintf(i_str, "%d", i);
                            sprintf(j_str, "%d", j);
                            sprintf(k_str, "%d", k);

                            strcpy(sheet_file, "/home/gustavo/work/gyst/sheets/");
                            strcat(sheet_file, dump_file);
                            strcat(sheet_file, "_");
                            strcat(sheet_file, i_str);
                            strcat(sheet_file, "_");
                            strcat(sheet_file, j_str);
                            strcat(sheet_file, "_");
                            strcat(sheet_file, k_str);
                            strcat(sheet_file, "_s");

                            fp = fopen(sheet_file, "w");
                            //fprintf(fp, "%lf\t%lf\t%lf\t%lf\t%lf\n",
                            //    Br_center, Bth_center, Bphi_center, beta_center, sigma_center);
                            for (m = 0; m < size_B; m++) {
                                //fprintf(fp, "%lf\t%lf\n", B_arr[m], J_arr[m]);
                                fprintf(fp, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\t%d\t0\n",
                                Br_arr[m], Bth_arr[m], Bphi_arr[m], J_arr[m],
                                beta_arr[m], sigma_arr[m], sigmaphi_arr[m], i_arr[m], j_arr[m]);
                            }
                            fclose(fp);
                        }
                    }
                }
            }
        }
    }
    printf("Identified %d possible reconnection sites.\n", is_good_count);
    write_current_sheets(jchar_output, J_cs_char);
    printf("Finished characterization of sheets.\n");
}


void xy_to_rth(double x, double y, double *r, double *th)
{
    // change to BL?
    *r = sqrt(x*x + y*y);
    *th = atan(y/x);
}

void rth_to_xy(double r, double th, double *x, double *y)
{
    // change to BL?
    *x = r*cos(th);
    *y = r*sin(th);
}

void vec_pol_to_cart(double eig_pol[2], double eig_cart[2], double theta)
{
    int i;
    double pol_to_cart[2][2];

    pol_to_cart[0][0] = cos(theta);
    pol_to_cart[0][1] = -sin(theta);
    pol_to_cart[1][0] = sin(theta);
    pol_to_cart[1][1] = cos(theta);

    for (i = 0; i < 2; i++) {
        eig_cart[i] = 0.;
    }

    for (i = 0; i < 2; i++) {
        eig_cart[0] += eig_pol[i]*pol_to_cart[i][0];
        eig_cart[1] += eig_pol[i]*pol_to_cart[i][1];
    }
}




/******************************************************************************/



void characterize2D_single()
{
    int i, j, k, ii, jj, kk, m;
    double rr, tth, pphi, p_r, p_th, p_phi, x, y, p_x, p_y;
    double Hess[2][2], e1[2], e1_cart[2], e2[2], evec[4];
    double ap_x[NDIM], p_cart[NDIM], Vprim[NDIM], Bprim[NDIM];
    double del[NDIM];
    double new_position[3];
    double V_proj_e1, B_proj_e1, B_proj_e2;
    double step, stepsize, B_upper, B_lower, V_upper, V_lower, V_A, V_rec;
    double Br_center, Bth_center, Bphi_center, beta_center, sigma_center;

    int is_good_count, is_good, is_good_lower, is_good_upper, size_lower, size_upper, size_B;
    char sheet_file[100], i_str[6], j_str[6], k_str[6];
    FILE *fp;

//    double B_upper_arr[100], B_lower_arr[100], J_upper_arr[100], J_lower_arr[100];
//    double *B_arr, *J_arr;

    double Br_upper_arr[100], Br_lower_arr[100], Bth_upper_arr[100], Bth_lower_arr[100], Bphi_upper_arr[100], Bphi_lower_arr[100];
    double J_upper_arr[100], J_lower_arr[100], beta_upper_arr[100], beta_lower_arr[100];
    double sigma_upper_arr[100], sigma_lower_arr[100], sigmaphi_upper_arr[100], sigmaphi_lower_arr[100];
    int i_upper_arr[100], i_lower_arr[100], j_upper_arr[100], j_lower_arr[100];

    double Br_arr[200], Bth_arr[200], Bphi_arr[200], J_arr[200], sigma_arr[200], sigmaphi_arr[200], beta_arr[200];
    int i_arr[200], j_arr[200];

    if (SHEETS) printf("\n");
    printf("Beginning characterization...\n");

    printf("Getting 1st derivatives...\n");
    for (i = 0; i < N1; i++) {
        for (j = 0; j < N2; j++) {
            for (k = 0; k < N3; k++) {
                get_1stderivative2D(i, j, k, J_cs);
                J_cs_char[i][j][k] = J_cs_peak[i][j][k];
            }
        }
    }

    i = 650;
    j = 327;
    //i = 534;
    //j = 72;
    k = 0;


    is_good_count = 0;
    printf("Entering main loop...\n");

    //if(J_cs_peak[i][j][k] <= 0) printf("jcs %d %d %d %lf\n", i, j, k, J_cs[i][j][k]);
    if (J_cs_peak[i][j][k] > 0.) {
        printf("%d %d\n", i, j);
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

        for (ii = 0; ii < 2; ii++) {
            for (jj = 0; jj < 2; jj++) {
                printf("%lf\n", Hess[ii][jj]);
            }
        }

        //2. get eigenvectors (max: e1, min: e3)
        get_evec2D(Hess, evec);
        for (ii = 0; ii < 2; ii++) {
            e1[ii] = evec[ii];
            e2[ii] = evec[ii+2];
        }

        //3. convert r, theta, to cartesian
        rth_to_xy(rr, tth, &x, &y);

        //4. convert max_eigvec to cartesian
        vec_pol_to_cart(e1, e1_cart, tth);
        for (ii = 0; ii < 2; ii++) {
            printf("%lf %lf\n", e1[ii], e1_cart[ii]);
        }

        //printf("%d %d %d %lf %lf\n", i, j, k, B[3][i][j][k], J[i][j][k]);
        //5. move along eigenvector (upper)
        step = STEP;
        stepsize = step;
        m = 0;
        //for (int index = 0; index < 100; index++) B_upper_arr[index] = 0.;
        is_good_upper = 0;
        while(m < 80) {
            p_cart[0] = x;
            p_cart[1] = y;

            //5.1 find cartesian coordinates of point (find p_x, p_y, p_z)
            for (ii = 0; ii < 2; ii++) new_position[ii] = p_cart[ii] + e1_cart[ii]*step;
            p_x = new_position[0];
            p_y = new_position[1];
            step += stepsize;

            //printf("p_cart %lf %lf\n", p_cart[0], p_cart[1]);
            //printf("e1 %lf %lf\n", e1[0], e1[1]);
            //printf("e1_cart %lf %lf\n", e1_cart[0], e1_cart[1]);
            //printf("new_p %lf %lf\n", p_x, p_y);

            //5.2 convert point to polar
            xy_to_rth(p_x, p_y, &p_r, &p_th);

            //6. find x1, x2, x3 corresponding to them
            ap_x[1] = zbrent(find_x1_cyl, p_r, 0, 8., 1E-6);
            x1in = ap_x[1];
            ap_x[2] = zbrent(find_x2_cyl, p_th, -1.0, 1.0, 1E-6); // check this
            ap_x[3] = p_phi;

            //7. find closest cell corresponding to spherical
            Xtoijk(ap_x, &ii, &jj, &kk, del);
            printf("%lf %lf %lf %lf %lf %lf %lf %lf\n", x, y, p_x, p_y, p_r, p_th, ap_x[1], ap_x[2]);
            printf("%d %d %lf %lf\n", ii, jj, ap_x[1], new_position[2]);

            //7.1 test value of J_cs
            if (J[ii][jj][kk] > 0.5*J_cs_peak[i][j][k]) {

                Vprim[0] = 0.;
                Vprim[1] = V[1][ii][jj][kk];
                Vprim[2] = V[2][ii][jj][kk];
                Vprim[3] = V[3][ii][jj][kk];
                Bprim[0] = 0.;
                Bprim[1] = B[1][ii][jj][kk];
                Bprim[2] = B[2][ii][jj][kk];
                Bprim[3] = B[3][ii][jj][kk];

                //9. project V, B along e1, e2, e3 to find V, B along these eigenvectors
                // NOTE: may have to change, use 4vec and recalculate V, B
                /*
                V_proj_e1 = dot3(Vprim, e1) / (sqrt(dot3(e1,e1)) + SMALL);
                B_proj_e1 = dot3(Bprim, e1) / (sqrt(dot3(e1,e1)) + SMALL);
                B_proj_e2 = dot3(Bprim, e2) / (sqrt(dot3(e2,e2)) + SMALL);
                B_upper = abs(B_proj_e1);
                B_upper_arr[m] = B_proj_e1;
                */

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

                printf("upper %lf %lf %lf %lf %lf %lf %lf %d %d\n",
                    Br_upper_arr[m], Bth_upper_arr[m], Bphi_upper_arr[m], J_upper_arr[m],
                    beta_upper_arr[m], sigma_upper_arr[m], sigmaphi_upper_arr[m],
                    i_upper_arr[m], j_upper_arr[m]);

                //printf("upper %d %d %d %lf %lf\n", ii, jj, kk, Bprim[3], J[ii][jj][kk]);

                m++;
                size_upper = m;

                //10. use B_proj to find v_alfven
                /*
                V_A = sqrt((B_proj_e1*B_proj_e1 +
                            B_proj_e2*B_proj_e2)/(rho[ii][jj][kk] + SMALL));
                V_upper = V_proj_e1/V_A;
                */

                V_A = sqrt((Bprim[1]*Bprim[1] +
                            Bprim[2]*Bprim[2] +
                            Bprim[3]*Bprim[3])/(rho[ii][jj][kk] + SMALL));
                V_upper = Vprim[1]/V_A;

                is_good_upper = 1;
            }
            else break;
        }

        //5. move along eigenvector (lower)
        step = STEP;
        stepsize = step;
        m = 0;
        //for (int index = 0; index < 100; index++) B_lower_arr[index] = 0.;
        is_good_lower = 0;
        while(m < 80) {
            p_cart[0] = x;
            p_cart[1] = y;

            //5.1 find cartesian coordinates of point (find p_x, p_y)
            for (ii = 0; ii < 2; ii++) new_position[ii] = p_cart[ii] - e1_cart[ii]*step;
            p_x = new_position[0];
            p_y = new_position[1];
            step += stepsize;

            //5.2 convert point to spherical
            xy_to_rth(p_x, p_y, &p_r, &p_th);

            //6. find x1, x2, x3 corresponding to them
            ap_x[1] = zbrent(find_x1_cyl, p_r, 0, 8., 1E-6);
            x1in = ap_x[1];
            ap_x[2] = zbrent(find_x2_cyl, p_th, -1.0, 1.0, 1E-6);
            ap_x[3] = p_phi;

            //7. find closest cell corresponding to spherical using Xtoijk
            Xtoijk(ap_x, &ii, &jj, &kk, del);

            //7.1 test value of J_cs
            if (J[ii][jj][kk] > 0.5*J_cs_peak[i][j][k]) {

                //8. project V along max_evec to find V_in at this cell
                // NOTE: may have to change, use 4vec and recalculate V, B
                Vprim[0] = 0.;
                Vprim[1] = V[1][ii][jj][kk];
                Vprim[2] = V[2][ii][jj][kk];
                Vprim[3] = V[3][ii][jj][kk];
                Bprim[0] = 0.;
                Bprim[1] = B[1][ii][jj][kk];
                Bprim[2] = B[2][ii][jj][kk];
                Bprim[3] = B[3][ii][jj][kk];

                //9. project B along e1, e2, e3 to find B along these eigenvectors
                /*
                V_proj_e1 = dot3(Vprim, e1) / (sqrt(dot3(e1,e1)) + SMALL);
                B_proj_e1 = dot3(Bprim, e1) / (sqrt(dot3(e1,e1)) + SMALL);
                B_proj_e2 = dot3(Bprim, e2) / (sqrt(dot3(e2,e2)) + SMALL);
                B_lower = abs(B_proj_e1);
                B_lower_arr[m] = B_proj_e1;
                */

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

                printf("lower %lf %lf %lf %lf %lf %lf %lf %d %d\n",
                    Br_lower_arr[m], Bth_lower_arr[m], Bphi_lower_arr[m], J_lower_arr[m],
                    beta_lower_arr[m], sigma_lower_arr[m], sigmaphi_lower_arr[m],
                    i_lower_arr[m], j_lower_arr[m]);

                //printf("size lower %d, m %d, B %lf\n", size_lower, m, B_lower_arr[m]);
                m++;
                size_lower = m;

                //10. use B_proj to find v_alfven
                /*
                V_A = sqrt((B_proj_e1*B_proj_e1 +
                            B_proj_e2*B_proj_e2)/(rho[ii][jj][kk] + SMALL));
                V_lower = V_proj_e1/V_A;
                */

                V_A = sqrt((Bprim[1]*Bprim[1] +
                            Bprim[2]*Bprim[2] +
                            Bprim[3]*Bprim[3])/(rho[ii][jj][kk] + SMALL));
                V_lower = Vprim[1]/V_A;

                is_good_lower = 1;
            }
            else break;
        }

        if (is_good_lower && is_good_upper) {
            //12. get V_in/V_A
            V_rec = 0.5*(V_lower - V_upper);

            is_good = 1;
            // don't consider point if smallish reconnection velocity
            if (abs(V_rec) < 0.6) {
                J_cs_char[i][j][k] = 0.;
                is_good = 0;
            }

            // don't consider point if too asymetrical
            if (abs(B_upper) < 0.9*abs(B_lower) || abs(B_lower) < 0.9*abs(B_upper)) {
                J_cs_char[i][j][k] = 0.;
                is_good = 0;
            }

            if (B_upper * B_lower > 0) {
                J_cs_char[i][j][k] = 0.;
                is_good = 0;
            }

            //if (is_good) printf("is good %d %d %d\n", i, j, k);

            if (is_good) {
                is_good_count++;
                size_B = size_lower + size_upper - 2;

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
                }

                reverse_array(Br_lower_arr, size_lower);
                merge_arrays(Br_lower_arr, Br_upper_arr, Br_arr, size_lower, size_upper);
                reverse_array(Bth_lower_arr, size_lower);
                merge_arrays(Bth_lower_arr, Bth_upper_arr, Bth_arr, size_lower, size_upper);
                reverse_array(Bphi_lower_arr, size_lower);
                merge_arrays(Bphi_lower_arr, Bphi_upper_arr, Bphi_arr, size_lower, size_upper);
                reverse_array(J_lower_arr, size_lower);
                merge_arrays(J_lower_arr, J_upper_arr, J_arr, size_lower, size_upper);
                reverse_array(beta_lower_arr, size_lower);
                merge_arrays(beta_lower_arr, beta_upper_arr, beta_arr, size_lower, size_upper);
                reverse_array(sigma_lower_arr, size_lower);
                merge_arrays(sigma_lower_arr, sigma_upper_arr, sigma_arr, size_lower, size_upper);
                reverse_array(sigmaphi_lower_arr, size_lower);
                merge_arrays(sigmaphi_lower_arr, sigmaphi_lower_arr, sigmaphi_arr, size_lower, size_upper);
                reverse_array_int(i_lower_arr, size_lower);
                merge_arrays_int(i_lower_arr, i_lower_arr, i_arr, size_lower, size_upper);
                reverse_array_int(j_lower_arr, size_lower);
                merge_arrays_int(j_lower_arr, j_lower_arr, j_arr, size_lower, size_upper);

                //for (m = 0; m < size_B; m++) printf("total %d\n", i_arr[m]);

                sprintf(i_str, "%d", i);
                sprintf(j_str, "%d", j);
                sprintf(k_str, "%d", k);
                strcpy(sheet_file, "/home/gustavo/work/gyst/sheets/");
                strcat(sheet_file, dump_file);
                strcat(sheet_file, "_");
                strcat(sheet_file, i_str);
                strcat(sheet_file, "_");
                strcat(sheet_file, j_str);
                strcat(sheet_file, "_");
                strcat(sheet_file, k_str);
                strcat(sheet_file, "_s");

                fp = fopen(sheet_file, "w");
                //fprintf(fp, "%lf\t%lf\t%lf\t%lf\t%lf\n",
                //    Br_center, Bth_center, Bphi_center, beta_center, sigma_center);
                for (m = 0; m < size_B; m++) {
                    fprintf(fp, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\t%d\t0\n",
                    Br_arr[m], Bth_arr[m], Bphi_arr[m], J_arr[m],
                    beta_arr[m], sigma_arr[m], sigmaphi_arr[m], i_arr[m], j_arr[m]);

                    printf("%lf,%lf,%lf,%lf,%lf,%lf,%lf,%d,%d,0\n",
                    Br_arr[m], Bth_arr[m], Bphi_arr[m], J_arr[m],
                    beta_arr[m], sigma_arr[m], sigmaphi_arr[m], i_arr[m], j_arr[m]);
                }
                fclose(fp);
            }
        }
    }

    printf("Identified %d possible reconnection sites.\n", is_good_count);
    write_current_sheets(jchar_output, J_cs_char);
    printf("Finished characterization of sheets.\n");
}




/******************************************************************************/



void characterize2D_single2()
{
    int i, j, k, ii, jj, kk, m;
    double rr, tth, pphi, p_r, p_th, p_phi, x, y, p_x, p_y;
    double Hess[2][2], e1[2], e1_cart[2], e2[2], evec[4];
    double ap_x[NDIM], p_cart[NDIM], Vprim[NDIM], Bprim[NDIM];
    double del[NDIM];
    double new_position[3];
    double V_proj_e1, B_proj_e1, B_proj_e2;
    double step, stepsize, B_upper, B_lower, V_upper, V_lower, V_A, V_rec;
    double Br_center, Bth_center, Bphi_center, beta_center, sigma_center;

    int is_good_count, is_good, is_good_lower, is_good_upper, size_lower, size_upper, size_B;
    char sheet_file[100], i_str[6], j_str[6], k_str[6];
    FILE *fp;

//    double B_upper_arr[100], B_lower_arr[100], J_upper_arr[100], J_lower_arr[100];
//    double *B_arr, *J_arr;

    double Br_upper_arr[100], Br_lower_arr[100], Bth_upper_arr[100], Bth_lower_arr[100], Bphi_upper_arr[100], Bphi_lower_arr[100];
    double J_upper_arr[100], J_lower_arr[100], beta_upper_arr[100], beta_lower_arr[100];
    double sigma_upper_arr[100], sigma_lower_arr[100], sigmaphi_upper_arr[100], sigmaphi_lower_arr[100];
    int i_upper_arr[100], i_lower_arr[100], j_upper_arr[100], j_lower_arr[100];

    double Br_arr[200], Bth_arr[200], Bphi_arr[200], J_arr[200], sigma_arr[200], sigmaphi_arr[200], beta_arr[200];
    int i_arr[200], j_arr[200];

    if (SHEETS) printf("\n");
    printf("Beginning characterization...\n");

    printf("Getting 1st derivatives...\n");
    for (i = 0; i < N1; i++) {
        for (j = 0; j < N2; j++) {
            for (k = 0; k < N3; k++) {
                get_1stderivative2D(i, j, k, J_cs);
                J_cs_char[i][j][k] = J_cs_peak[i][j][k];
            }
        }
    }

    i = 650;
    j = 327;
    //i = 534;
    //j = 72;
    k = 0;

    is_good_count = 0;
    printf("Entering main loop...\n");

    //if(J_cs_peak[i][j][k] <= 0) printf("jcs %d %d %d %lf\n", i, j, k, J_cs[i][j][k]);
    if (J_cs_peak[i][j][k] > 0.) {
        printf("%d %d\n", i, j);
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

        for (ii = 0; ii < 2; ii++) {
            for (jj = 0; jj < 2; jj++) {
                printf("%lf\n", Hess[ii][jj]);
            }
        }

        //2. get eigenvectors (max: e1, min: e3)
        get_evec2D(Hess, evec);
        for (ii = 0; ii < 2; ii++) {
            e1[ii] = evec[ii];
            printf("%lf\n", e1[ii]);
            e2[ii] = evec[ii+2];
        }

        printf("%d %d %lf %lf\n", i, j, rr, tth);

        //3. convert r, theta, to cartesian
        //rth_to_xy(rr, tth, &x, &y);

        //4. convert max_eigvec to cartesian
        //vec_pol_to_cart(e1, e1_cart, tth);

        //printf("%d %d %d %lf %lf\n", i, j, k, B[3][i][j][k], J[i][j][k]);
        //5. move along eigenvector (upper)
        step = STEP2;
        stepsize = step;
        m = 0;
        is_good_upper = 0;
        while(m < 80) {
            p_cart[0] = rr;
            p_cart[1] = tth;

            //5.1 find cartesian coordinates of point (find p_x, p_y, p_z)
            for (ii = 0; ii < 2; ii++) new_position[ii] = p_cart[ii] + e1[ii]*step;
            ap_x[1] = new_position[0];
            ap_x[2] = new_position[1];
            ap_x[3] = 0.;
            step += stepsize;

            //printf("p_cart %lf %lf\n", p_cart[0], p_cart[1]);
            //printf("e1 %lf %lf\n", e1[0], e1[1]);
            //printf("e1_cart %lf %lf\n", e1_cart[0], e1_cart[1]);
            //printf("new_p %lf %lf\n", p_x, p_y);

            //5.2 convert point to polar
            //xy_to_rth(p_x, p_y, &p_r, &p_th);

            //6. find x1, x2, x3 corresponding to them
            //ap_x[1] = zbrent(find_x1_cyl, p_r, 0, 8., 1E-6);
            //x1in = ap_x[1];
            //ap_x[2] = zbrent(find_x2_cyl, p_th, -1.0, 1.0, 1E-6); // check this
            //ap_x[3] = p_phi;

            //7. find closest cell corresponding to spherical
            Xtoijk_new(ap_x, &ii, &jj, &kk, del);
            printf("%d %d %lf %lf\n", ii, jj, new_position[0], new_position[1]);


            //7.1 test value of J_cs
            if (J[ii][jj][kk] > 0.5*J_cs_peak[i][j][k]) {

                Vprim[0] = 0.;
                Vprim[1] = V[1][ii][jj][kk];
                Vprim[2] = V[2][ii][jj][kk];
                Vprim[3] = V[3][ii][jj][kk];
                Bprim[0] = 0.;
                Bprim[1] = B[1][ii][jj][kk];
                Bprim[2] = B[2][ii][jj][kk];
                Bprim[3] = B[3][ii][jj][kk];

                //9. project V, B along e1, e2, e3 to find V, B along these eigenvectors
                // NOTE: may have to change, use 4vec and recalculate V, B
                /*
                V_proj_e1 = dot3(Vprim, e1) / (sqrt(dot3(e1,e1)) + SMALL);
                B_proj_e1 = dot3(Bprim, e1) / (sqrt(dot3(e1,e1)) + SMALL);
                B_proj_e2 = dot3(Bprim, e2) / (sqrt(dot3(e2,e2)) + SMALL);
                B_upper = abs(B_proj_e1);
                B_upper_arr[m] = B_proj_e1;
                */

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

                printf("upper %lf %lf %lf %lf %lf %lf %lf %d %d\n",
                    Br_upper_arr[m], Bth_upper_arr[m], Bphi_upper_arr[m], J_upper_arr[m],
                    beta_upper_arr[m], sigma_upper_arr[m], sigmaphi_upper_arr[m],
                    i_upper_arr[m], j_upper_arr[m]);

                //printf("upper %d %d %d %lf %lf\n", ii, jj, kk, Bprim[3], J[ii][jj][kk]);

                m++;
                size_upper = m;

                //10. use B_proj to find v_alfven
                /*
                V_A = sqrt((B_proj_e1*B_proj_e1 +
                            B_proj_e2*B_proj_e2)/(rho[ii][jj][kk] + SMALL));
                V_upper = V_proj_e1/V_A;
                */

                V_A = sqrt((Bprim[1]*Bprim[1] +
                            Bprim[2]*Bprim[2] +
                            Bprim[3]*Bprim[3])/(rho[ii][jj][kk] + SMALL));
                V_upper = Vprim[1]/V_A;

                is_good_upper = 1;
            }
            else break;
        }

        //5. move along eigenvector (lower)
        step = STEP2;
        stepsize = step;
        m = 0;
        //for (int index = 0; index < 100; index++) B_lower_arr[index] = 0.;
        is_good_lower = 0;
        while(m < 80) {
            p_cart[0] = rr;
            p_cart[1] = tth;

            //5.1 find cartesian coordinates of point (find p_x, p_y)
            for (ii = 0; ii < 2; ii++) new_position[ii] = p_cart[ii] - e1[ii]*step;
            ap_x[1] = new_position[0];
            ap_x[2] = new_position[1];
            ap_x[3] = 0.;
            step += stepsize;

            //5.2 convert point to spherical
            //xy_to_rth(p_x, p_y, &p_r, &p_th);

            //6. find x1, x2, x3 corresponding to them
            //ap_x[1] = zbrent(find_x1_cyl, p_r, 0, 8., 1E-6);
            //x1in = ap_x[1];
            //ap_x[2] = zbrent(find_x2_cyl, p_th, -1.0, 1.0, 1E-6);
            //ap_x[3] = p_phi;

            //7. find closest cell corresponding to spherical using Xtoijk
            Xtoijk_new(ap_x, &ii, &jj, &kk, del);

            //7.1 test value of J_cs
            if (J[ii][jj][kk] > 0.5*J_cs_peak[i][j][k]) {

                //8. project V along max_evec to find V_in at this cell
                // NOTE: may have to change, use 4vec and recalculate V, B
                Vprim[0] = 0.;
                Vprim[1] = V[1][ii][jj][kk];
                Vprim[2] = V[2][ii][jj][kk];
                Vprim[3] = V[3][ii][jj][kk];
                Bprim[0] = 0.;
                Bprim[1] = B[1][ii][jj][kk];
                Bprim[2] = B[2][ii][jj][kk];
                Bprim[3] = B[3][ii][jj][kk];

                //9. project B along e1, e2, e3 to find B along these eigenvectors
                /*
                V_proj_e1 = dot3(Vprim, e1) / (sqrt(dot3(e1,e1)) + SMALL);
                B_proj_e1 = dot3(Bprim, e1) / (sqrt(dot3(e1,e1)) + SMALL);
                B_proj_e2 = dot3(Bprim, e2) / (sqrt(dot3(e2,e2)) + SMALL);
                B_lower = abs(B_proj_e1);
                B_lower_arr[m] = B_proj_e1;
                */

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

                printf("lower %lf %lf %lf %lf %lf %lf %lf %d %d\n",
                    Br_lower_arr[m], Bth_lower_arr[m], Bphi_lower_arr[m], J_lower_arr[m],
                    beta_lower_arr[m], sigma_lower_arr[m], sigmaphi_lower_arr[m],
                    i_lower_arr[m], j_lower_arr[m]);

                //printf("size lower %d, m %d, B %lf\n", size_lower, m, B_lower_arr[m]);
                m++;
                size_lower = m;

                //10. use B_proj to find v_alfven
                /*
                V_A = sqrt((B_proj_e1*B_proj_e1 +
                            B_proj_e2*B_proj_e2)/(rho[ii][jj][kk] + SMALL));
                V_lower = V_proj_e1/V_A;
                */

                V_A = sqrt((Bprim[1]*Bprim[1] +
                            Bprim[2]*Bprim[2] +
                            Bprim[3]*Bprim[3])/(rho[ii][jj][kk] + SMALL));
                V_lower = Vprim[1]/V_A;

                is_good_lower = 1;
            }
            else break;
        }

        if (is_good_lower && is_good_upper) {
            //12. get V_in/V_A
            V_rec = 0.5*(V_lower - V_upper);

            is_good = 1;
            // don't consider point if smallish reconnection velocity
            if (abs(V_rec) < 0.6) {
                J_cs_char[i][j][k] = 0.;
                is_good = 0;
            }

            // don't consider point if too asymetrical
            if (abs(B_upper) < 0.9*abs(B_lower) || abs(B_lower) < 0.9*abs(B_upper)) {
                J_cs_char[i][j][k] = 0.;
                is_good = 0;
            }

            if (B_upper * B_lower > 0) {
                J_cs_char[i][j][k] = 0.;
                is_good = 0;
            }

            //if (is_good) printf("is good %d %d %d\n", i, j, k);

            if (is_good) {
                is_good_count++;
                size_B = size_lower + size_upper - 2;

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
                }

                reverse_array(Br_lower_arr, size_lower);
                merge_arrays(Br_lower_arr, Br_upper_arr, Br_arr, size_lower, size_upper);
                reverse_array(Bth_lower_arr, size_lower);
                merge_arrays(Bth_lower_arr, Bth_upper_arr, Bth_arr, size_lower, size_upper);
                reverse_array(Bphi_lower_arr, size_lower);
                merge_arrays(Bphi_lower_arr, Bphi_upper_arr, Bphi_arr, size_lower, size_upper);
                reverse_array(J_lower_arr, size_lower);
                merge_arrays(J_lower_arr, J_upper_arr, J_arr, size_lower, size_upper);
                reverse_array(beta_lower_arr, size_lower);
                merge_arrays(beta_lower_arr, beta_upper_arr, beta_arr, size_lower, size_upper);
                reverse_array(sigma_lower_arr, size_lower);
                merge_arrays(sigma_lower_arr, sigma_upper_arr, sigma_arr, size_lower, size_upper);
                reverse_array(sigmaphi_lower_arr, size_lower);
                merge_arrays(sigmaphi_lower_arr, sigmaphi_lower_arr, sigmaphi_arr, size_lower, size_upper);
                reverse_array_int(i_lower_arr, size_lower);
                merge_arrays_int(i_lower_arr, i_lower_arr, i_arr, size_lower, size_upper);
                reverse_array_int(j_lower_arr, size_lower);
                merge_arrays_int(j_lower_arr, j_lower_arr, j_arr, size_lower, size_upper);

                //for (m = 0; m < size_B; m++) printf("total %d\n", i_arr[m]);

                sprintf(i_str, "%d", i);
                sprintf(j_str, "%d", j);
                sprintf(k_str, "%d", k);
                strcpy(sheet_file, "/home/gustavo/work/gyst/sheets/");
                strcat(sheet_file, dump_file);
                strcat(sheet_file, "_");
                strcat(sheet_file, i_str);
                strcat(sheet_file, "_");
                strcat(sheet_file, j_str);
                strcat(sheet_file, "_");
                strcat(sheet_file, k_str);
                strcat(sheet_file, "_s");

                fp = fopen(sheet_file, "w");
                //fprintf(fp, "%lf\t%lf\t%lf\t%lf\t%lf\n",
                //    Br_center, Bth_center, Bphi_center, beta_center, sigma_center);
                for (m = 0; m < size_B; m++) {
                    fprintf(fp, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\t%d\t0\n",
                    Br_arr[m], Bth_arr[m], Bphi_arr[m], J_arr[m],
                    beta_arr[m], sigma_arr[m], sigmaphi_arr[m], i_arr[m], j_arr[m]);

                    printf("%lf,%lf,%lf,%lf,%lf,%lf,%lf,%d,%d,0\n",
                    Br_arr[m], Bth_arr[m], Bphi_arr[m], J_arr[m],
                    beta_arr[m], sigma_arr[m], sigmaphi_arr[m], i_arr[m], j_arr[m]);
                }
                fclose(fp);
            }
        }
    }

    printf("Identified %d possible reconnection sites.\n", is_good_count);
    write_current_sheets(jchar_output, J_cs_char);
    printf("Finished characterization of sheets.\n");
}



/******************************************************************************/


void characterize2D2()
{
    int i, j, k, ii, jj, kk, m;
    double rr, tth, pphi, p_r, p_th, p_phi, x, y, p_x, p_y;
    double Hess[2][2], e1[2], e1_cart[2], e2[2], evec[4];
    double ap_x[NDIM], p_cart[NDIM], Vprim[NDIM], Bprim[NDIM];
    double del[NDIM];
    double new_position[3];
    double V_proj_e1, B_proj_e1, B_proj_e2;
    double step, stepsize, B_upper, B_lower, V_upper, V_lower, V_A, V_rec;
    double Br_center, Bth_center, Bphi_center, beta_center, sigma_center;

    int is_good_count, is_good, is_good_lower, is_good_upper, size_lower, size_upper, size_B;
    char sheet_file[100], i_str[6], j_str[6], k_str[6];
    FILE *fp;

//    double B_upper_arr[100], B_lower_arr[100], J_upper_arr[100], J_lower_arr[100];
//    double *B_arr, *J_arr;

    double Br_upper_arr[100], Br_lower_arr[100], Bth_upper_arr[100], Bth_lower_arr[100], Bphi_upper_arr[100], Bphi_lower_arr[100];
    double J_upper_arr[100], J_lower_arr[100], beta_upper_arr[100], beta_lower_arr[100];
    double sigma_upper_arr[100], sigma_lower_arr[100], sigmaphi_upper_arr[100], sigmaphi_lower_arr[100];
    int i_upper_arr[100], i_lower_arr[100], j_upper_arr[100], j_lower_arr[100];

    double Br_arr[200], Bth_arr[200], Bphi_arr[200], J_arr[200], sigma_arr[200], sigmaphi_arr[200], beta_arr[200];
    int i_arr[200], j_arr[200];

    if (SHEETS) printf("\n");
    printf("Beginning characterization...\n");

    printf("Getting 1st derivatives...\n");
    for (i = 0; i < N1; i++) {
        for (j = 0; j < N2; j++) {
            for (k = 0; k < N3; k++) {
                get_1stderivative2D(i, j, k, J_cs);
                J_cs_char[i][j][k] = J_cs_peak[i][j][k];
            }
        }
    }

    is_good_count = 0;
    printf("Entering main loop...\n");
    for (i = 0; i < N1; i++) {
        for (j = 0; j < N2; j++) {
            for (k = 0; k < N3; k++) {

                //if(J_cs_peak[i][j][k] <= 0) printf("jcs %d %d %d %lf\n", i, j, k, J_cs[i][j][k]);
                if (J_cs_peak[i][j][k] > 0.) {
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

                    //3. convert r, theta, to cartesian
                    //rth_to_xy(rr, tth, &x, &y);

                    //4. convert max_eigvec to cartesian
                    //vec_pol_to_cart(e1, e1_cart, tth);

                    //printf("%d %d %d %lf %lf\n", i, j, k, B[3][i][j][k], J[i][j][k]);
                    //5. move along eigenvector (upper)
                    //step = STEP1;

                    if (BETWEEN(i, 0, 180)) step = STEP1;
                    else if (BETWEEN(i, 181, 510)) step = STEP2;
                    else if (BETWEEN(i, 511, 850)) step = STEP3;
                    else if (BETWEEN(i, 851, N1)) step = STEP4;
                    step = STEP1;
                    stepsize = step;
                    m = 0;
                    is_good_upper = 0;
                    while(m < 80) {
                        p_cart[0] = rr;
                        p_cart[1] = tth;

                        //5.1 find cartesian coordinates of point (find p_x, p_y, p_z)
                        for (ii = 0; ii < 2; ii++) new_position[ii] = p_cart[ii] + e1[ii]*step;
                        ap_x[1] = new_position[0];
                        ap_x[2] = new_position[1];
                        ap_x[3] = 0.;
                        step += stepsize;

                        //printf("p_cart %lf %lf\n", p_cart[0], p_cart[1]);
                        //printf("e1 %lf %lf\n", e1[0], e1[1]);
                        //printf("e1_cart %lf %lf\n", e1_cart[0], e1_cart[1]);
                        //printf("new_p %lf %lf\n", p_x, p_y);

                        //5.2 convert point to polar
                        //xy_to_rth(p_x, p_y, &p_r, &p_th);

                        //6. find x1, x2, x3 corresponding to them
                        //ap_x[1] = zbrent(find_x1_cyl, p_r, 0, 8., 1E-6);
                        //x1in = ap_x[1];
                        //ap_x[2] = zbrent(find_x2_cyl, p_th, -1.0, 1.0, 1E-6); // check this
                        //ap_x[3] = p_phi;

                        //7. find closest cell corresponding to spherical
                        Xtoijk_new(ap_x, &ii, &jj, &kk, del);

                        //7.1 test value of J_cs
                        if (J[ii][jj][kk] > 0.5*J_cs_peak[i][j][k]) {

                            Vprim[0] = 0.;
                            Vprim[1] = V[1][ii][jj][kk];
                            Vprim[2] = V[2][ii][jj][kk];
                            Vprim[3] = V[3][ii][jj][kk];
                            Bprim[0] = 0.;
                            Bprim[1] = B[1][ii][jj][kk];
                            Bprim[2] = B[2][ii][jj][kk];
                            Bprim[3] = B[3][ii][jj][kk];

                            //9. project V, B along e1, e2, e3 to find V, B along these eigenvectors
                            // NOTE: may have to change, use 4vec and recalculate V, B
                            /*
                            V_proj_e1 = dot3(Vprim, e1) / (sqrt(dot3(e1,e1)) + SMALL);
                            B_proj_e1 = dot3(Bprim, e1) / (sqrt(dot3(e1,e1)) + SMALL);
                            B_proj_e2 = dot3(Bprim, e2) / (sqrt(dot3(e2,e2)) + SMALL);
                            B_upper = abs(B_proj_e1);
                            B_upper_arr[m] = B_proj_e1;
                            */

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
                            //upper_arr[m] = "upper";


                            //printf("upper %d %d %d %lf %lf\n", ii, jj, kk, Bprim[3], J[ii][jj][kk]);

                            m++;
                            size_upper = m;

                            //10. use B_proj to find v_alfven
                            /*
                            V_A = sqrt((B_proj_e1*B_proj_e1 +
                                        B_proj_e2*B_proj_e2)/(rho[ii][jj][kk] + SMALL));
                            V_upper = V_proj_e1/V_A;
                            */

                            V_A = sqrt((Bprim[1]*Bprim[1] +
                                        Bprim[2]*Bprim[2] +
                                        Bprim[3]*Bprim[3])/(rho[ii][jj][kk] + SMALL));
                            V_upper = Vprim[1]/V_A;

                            is_good_upper = 1;
                        }
                        else break;
                    }

                    //5. move along eigenvector (lower)
                    if (BETWEEN(i, 0, 180)) step = STEP1;
                    else if (BETWEEN(i, 181, 510)) step = STEP2;
                    else if (BETWEEN(i, 511, 850)) step = STEP3;
                    else if (BETWEEN(i, 851, N1)) step = STEP4;
                    step = STEP1;
                    stepsize = step;
                    m = 0;
                    //for (int index = 0; index < 100; index++) B_lower_arr[index] = 0.;
                    is_good_lower = 0;
                    while(m < 80) {
                        p_cart[0] = rr;
                        p_cart[1] = tth;

                        //5.1 find cartesian coordinates of point (find p_x, p_y)
                        for (ii = 0; ii < 2; ii++) new_position[ii] = p_cart[ii] - e1[ii]*step;
                        ap_x[1] = new_position[0];
                        ap_x[2] = new_position[1];
                        ap_x[3] = 0.;
                        step += stepsize;

                        //5.2 convert point to spherical
                        //xy_to_rth(p_x, p_y, &p_r, &p_th);

                        //6. find x1, x2, x3 corresponding to them
                        //ap_x[1] = zbrent(find_x1_cyl, p_r, 0, 8., 1E-6);
                        //x1in = ap_x[1];
                        //ap_x[2] = zbrent(find_x2_cyl, p_th, -1.0, 1.0, 1E-6);
                        //ap_x[3] = p_phi;

                        //7. find closest cell corresponding to spherical using Xtoijk
                        Xtoijk_new(ap_x, &ii, &jj, &kk, del);

                        //7.1 test value of J_cs
                        if (J[ii][jj][kk] > 0.5*J_cs_peak[i][j][k]) {

                            //8. project V along max_evec to find V_in at this cell
                            // NOTE: may have to change, use 4vec and recalculate V, B
                            Vprim[0] = 0.;
                            Vprim[1] = V[1][ii][jj][kk];
                            Vprim[2] = V[2][ii][jj][kk];
                            Vprim[3] = V[3][ii][jj][kk];
                            Bprim[0] = 0.;
                            Bprim[1] = B[1][ii][jj][kk];
                            Bprim[2] = B[2][ii][jj][kk];
                            Bprim[3] = B[3][ii][jj][kk];

                            //9. project B along e1, e2, e3 to find B along these eigenvectors
                            /*
                            V_proj_e1 = dot3(Vprim, e1) / (sqrt(dot3(e1,e1)) + SMALL);
                            B_proj_e1 = dot3(Bprim, e1) / (sqrt(dot3(e1,e1)) + SMALL);
                            B_proj_e2 = dot3(Bprim, e2) / (sqrt(dot3(e2,e2)) + SMALL);
                            B_lower = abs(B_proj_e1);
                            B_lower_arr[m] = B_proj_e1;
                            */

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
                            //lower_arr[m] = "lower";


                            //printf("size lower %d, m %d, B %lf\n", size_lower, m, B_lower_arr[m]);
                            m++;
                            size_lower = m;

                            //10. use B_proj to find v_alfven
                            /*
                            V_A = sqrt((B_proj_e1*B_proj_e1 +
                                        B_proj_e2*B_proj_e2)/(rho[ii][jj][kk] + SMALL));
                            V_lower = V_proj_e1/V_A;
                            */

                            V_A = sqrt((Bprim[1]*Bprim[1] +
                                        Bprim[2]*Bprim[2] +
                                        Bprim[3]*Bprim[3])/(rho[ii][jj][kk] + SMALL));
                            V_lower = Vprim[1]/V_A;

                            is_good_lower = 1;
                        }
                        else break;
                    }

                    if (is_good_lower && is_good_upper) {
                        //12. get V_in/V_A
                        V_rec = 0.5*(V_lower - V_upper);

                        is_good = 1;
                        // don't consider point if smallish reconnection velocity
                        if (abs(V_rec) < 0.6) {
                            J_cs_char[i][j][k] = 0.;
                            is_good = 0;
                        }

                        // don't consider point if too asymetrical
                        if (abs(B_upper) < 0.9*abs(B_lower) || abs(B_lower) < 0.9*abs(B_upper)) {
                            J_cs_char[i][j][k] = 0.;
                            is_good = 0;
                        }

                        if (B_upper * B_lower > 0) {
                            J_cs_char[i][j][k] = 0.;
                            is_good = 0;
                        }

                        //if (is_good) printf("is good %d %d %d\n", i, j, k);

                        if (is_good) {
                            is_good_count++;
                            size_B = size_lower + size_upper - 2;

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
                            }

                            reverse_array(Br_lower_arr, size_lower);
                            merge_arrays(Br_lower_arr, Br_upper_arr, Br_arr, size_lower, size_upper);
                            reverse_array(Bth_lower_arr, size_lower);
                            merge_arrays(Bth_lower_arr, Bth_upper_arr, Bth_arr, size_lower, size_upper);
                            reverse_array(Bphi_lower_arr, size_lower);
                            merge_arrays(Bphi_lower_arr, Bphi_upper_arr, Bphi_arr, size_lower, size_upper);
                            reverse_array(J_lower_arr, size_lower);
                            merge_arrays(J_lower_arr, J_upper_arr, J_arr, size_lower, size_upper);
                            reverse_array(beta_lower_arr, size_lower);
                            merge_arrays(beta_lower_arr, beta_upper_arr, beta_arr, size_lower, size_upper);
                            reverse_array(sigma_lower_arr, size_lower);
                            merge_arrays(sigma_lower_arr, sigma_upper_arr, sigma_arr, size_lower, size_upper);
                            reverse_array(sigmaphi_lower_arr, size_lower);
                            merge_arrays(sigmaphi_lower_arr, sigmaphi_upper_arr, sigmaphi_arr, size_lower, size_upper);
                            reverse_array_int(i_lower_arr, size_lower);
                            merge_arrays_int(i_lower_arr, i_upper_arr, i_arr, size_lower, size_upper);
                            reverse_array_int(j_lower_arr, size_lower);
                            merge_arrays_int(j_lower_arr, j_upper_arr, j_arr, size_lower, size_upper);

                            //reverse_array_str(lower_arr, size_lower);
                            //merge_arrays_str(lower_arr, j_lower_arr, j_arr, size_lower, size_upper);

                            //for (m = 0; m < size_B; m++) printf("total %d\n", i_arr[m]);

                            sprintf(i_str, "%d", i);
                            sprintf(j_str, "%d", j);
                            sprintf(k_str, "%d", k);
                            strcpy(sheet_file, "/home/gustavo/work/gyst/sheets/");
                            strcat(sheet_file, dump_file);
                            strcat(sheet_file, "_");
                            strcat(sheet_file, i_str);
                            strcat(sheet_file, "_");
                            strcat(sheet_file, j_str);
                            strcat(sheet_file, "_");
                            strcat(sheet_file, k_str);
                            strcat(sheet_file, "_s");

                            fp = fopen(sheet_file, "w");
                            fprintf(fp, "Br\tBth\tBph\tJ\tbeta\tsigma\tsigmaphi\ti\tj\tk\n");
                            for (m = 0; m < size_B; m++) {
                                fprintf(fp, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\t%d\t0\n",
                                Br_arr[m], Bth_arr[m], Bphi_arr[m], J_arr[m],
                                beta_arr[m], sigma_arr[m], sigmaphi_arr[m], i_arr[m], j_arr[m]);

                            }
                            fclose(fp);
                        }
                    }
                }
            }
        }
    }
    printf("Identified %d possible reconnection sites.\n", is_good_count);
    write_current_sheets(jchar_output, J_cs_char);
    printf("Finished characterization of sheets.\n");
}
