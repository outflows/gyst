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

void characterize(int i, int j, int k)
{
    int ii, jj, kk;
    double rr, tth, pphi, p_r, p_th, p_phi;
    double x, y, z, p_x, p_y, p_z;
    double Hess[3][3];
    double e1[3], e1_cart[3], e2[3], e3[3];
    double my_x[NDIM];
    double del[NDIM];
    double V_proj_e1, B_proj_e1, B_proj_e2, B_proj_e3;
    double step;

    for (ii = 0; ii < N1; ii++) {
        for (jj = 0; jj < N2; jj++) {
            for (kk = 0; kk < N3; kk++) {
                J_cs_char[ii][jj][kk] = J_cs_peak[ii][jj][kk];
                get_1stderivative(ii, jj, kk, J_cs);
            }
        }
    }

    if (J_cs_peak[i][j][k] > 0.) {
        rr = r[i][j][k];
        tth = th[i][j][k];
        pphi = phi[i][j][k];
        //1. get hessian
        get_hessian(i, j, k, Hess);
        //2. get eigenvectors (max: e1, min: e3)
        evec = get_evec(Hess);
        for (ii = 0; ii < 3; ii++) {
            e1[ii] = evec[ii];
            e2[ii] = evec[ii+3];
            e3[ii] = evec[ii+6];
        }
        //3. convert r, theta, phi to cartesian
        rthphi_to_xyz(rr, tth, pphi, x, y, z);
        //4. convert max_eigvec to cartesian
        vec_sph_to_cart(e1, e1_cart, tth, pphi);

        //5. move along eigenvector (upper)
        step = 0.1;
        while(1) {
            p_cart[0] = x;
            p_cart[1] = y;
            p_cart[2] = z;
            //5.1 find cartesian coordinates of point (find p_x, p_y, p_z)
            for (ii = 0; ii < 3; ii++) new_position[ii] = p_cart[ii] + e1_cart[ii]*step;
            p_x = new_position[0];
            p_y = new_position[1];
            p_z = new_position[2];
            step += step;

            //5.2 convert point to spherical
            xyz_to_rthphi(p_x, p_y, p_z, p_r, p_th, p_phi);

            //6. find x1, x2, x3 corresponding to them
            p_x[1] = zbrent(find_x1_cyl, p_r, 0, 8., 1E-6);
            p_x[2] = zbrent(find_x2_cyl, p_th, -1.0, 1.0, 1E-6);
            p_x[3] = p_phi;

            //7. find closest cell corresponding to spherical
            Xtoijk(p_x, &ii, &jj, &kk, del);

            //7.1 test value of J_cs
            if (J_cs[ii][jj][kk] < 0.5*J_cs_peak[i][j][k]) break;

            //8. project V along max_evec to find V_in at this cell
            // NOTE: may have to change, use 4vec and recalculate V, B
            V_proj_e1 = dot3(V, e1) / sqrt(dot3(e1));
            //9. project B along e1, e2, e3 to find B along these eigenvectors
            B_proj_e1 = dot3(B, e1) / sqrt(dot3(e1));
            B_proj_e2 = dot3(B, e2) / sqrt(dot3(e2));
            B_proj_e3 = dot3(B, e3) / sqrt(dot3(e3));
            B_upper = abs(B_proj_e3);

            //10. use B_proj to find v_alfven
            V_A = sqrt((B_proj_e1*B_proj_e1 +
                        B_proj_e2*B_proj_e2 +
                        B_proj_e3*B_proj_e3)/rho[ii][jj][kk]);
            V_upper = V_proj_e1/V_A;
        }

        //5. move along eigenvector (lower)
        step = 0.1;
        while(1) {
            p_cart[0] = x;
            p_cart[1] = y;
            p_cart[2] = z;

            //5.1 find cartesian coordinates of point (find p_x, p_y, p_z)
            for (ii = 0; ii < 3; ii++) new_position[ii] = p_cart[ii] - e1_cart[ii]*step;
            p_x = new_position[0];
            p_y = new_position[1];
            p_z = new_position[2];
            step += step;

            //5.2 convert point to spherical
            xyz_to_rthphi(p_x, p_y, p_z, p_r, p_th, p_phi);
            //6. find x1, x2, x3 corresponding to them
            p_x[1] = zbrent(find_x1_cyl, p_r, 0, 8., 1E-6);
            p_x[2] = zbrent(find_x2_cyl, p_th, -1.0, 1.0, 1E-6);
            p_x[3] = p_phi;
            //7. find closest cell corresponding to spherical using Xtoijk
            Xtoijk(my_x, &ii, &jj, &kk, del);

            //7.1 test value of J_cs
            if (J_cs[ii][jj][kk] < 0.5*J_cs_peak[i][j][k]) break;

            //8. project V along max_evec to find V_in at this cell
            // NOTE: may have to change, use 4vec and recalculate V, B
            V_proj_e1 = dot3(V, e1) / sqrt(dot3(e1));
            //9. project B along e1, e2, e3 to find B along these eigenvectors
            B_proj_e1 = dot3(B, e1) / sqrt(dot3(e1));
            B_proj_e2 = dot3(B, e2) / sqrt(dot3(e2));
            B_proj_e3 = dot3(B, e3) / sqrt(dot3(e3));
            B_lower = abs(B_proj_e3);
            //10. use B_proj to find v_alfven
            V_A = sqrt((B_proj_e1*B_proj_e1 +
                        B_proj_e2*B_proj_e2 +
                        B_proj_e3*B_proj_e3)/rho[ii][jj][kk]);
            V_lower = V_proj_e1/V_A;
        }
        //12. get V_in/V_A
        V_rec = 0.5*(V_lower - V_upper);

        // don't consider point if smallish reconnection velocity
        if (V_rec < 0.7) {
            J_cs_char[i][j][k] = 0.;
        }
        // don't consider point if too asymetrical
        if (B_upper < 0.75*B_lower && B_lower < 0.75*B_upper) {
            J_cs_char[i][j][k] = 0.;
        }

    }

}


void xyz_to_rthphi(double x, double y, double z, double r, double th, double phi)
{
    // change to BL?
    r = sqrt(x*x + y*y + z*z);
    th = atan(sqrt(x*x + y*y)/z);
    phi = atan(y/x);
}

void rthphi_to_xyz(double r, double th, double phi, double x, double y, double z)
{
    // change to BL?
    x = r*sin(th)*cos(phi);
    y = r*sin(th)*sin(phi);
    z = r*cos(th);
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
