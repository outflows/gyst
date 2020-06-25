#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "definitions.h"
#include "declarations.h"

// grmonty units
void assign_units() {
    // Calculates the units for radiative transfer calculations.
    // Not needed here.

    printf("Assigning units...\n");

    L_unit = GNEWT * MBH / (CL * CL);
    T_unit = L_unit / CL;
    RHO_unit = MUNIT / (L_unit*L_unit*L_unit);
    U_unit = RHO_unit * CL * CL;
    B_unit = CL * sqrt(4. * M_PI * RHO_unit);
    Ne_unit = RHO_unit / (MP + ME);

    //printf("MBH: %g\n", MBH);
    //printf("M_unit: %g\n", MUNIT);
    //printf("L_unit: %g\n", L_unit);
    //printf("RHO_unit: %g\n", RHO_unit);
    //printf("U_unit: %g\n", U_unit);
    //printf("B_unit: %g\n", B_unit);
    //printf("Ne_unit: %g\n\n", Ne_unit);
}

void read_data(char *fname) {
    // Reads the data from a HARMPI dump file, stores these data in the
    // appropriate variables and uses them to calculate a few derived
    // quantities.
    //
    // This function is largely based on init_harm2dv3_data.c from grmonty,
    // with a few ideas taken from Gammie's file to read iharm2d data when
    // plotting.

    FILE *fp;
    int i, j, k;
    double tf, cour, DTd, DTl, DTi, DTr, DTr01, lim, failed;
    int N1tot, N2tot, N3tot, N1G, N2G, N3G;
    int nstep, dump_cnt, image_cnt, rdump_cnt, rdump01_cnt, NPR, DOKTOT;
    double fractheta, fracphi;
    //double ti, tj, tk, x1, x2, x3;//, r, th, phi;
    //double ktot, v1min, v1max, v2min, v2max, v3min, v3max;

    float var[50]; // must be float because HARMPI dump data is saved as float!

    printf("Opening dump file...\n");
    fp = fopen(fname, "r");
    if (fp == NULL) {
		fprintf(stderr, "dump file not found.\n");
		exit(1);
	} else {
		fprintf(stderr, "Successfully opened dump file.\n");
	}

    // get standard HARMPI header
    printf("Reading header...\n");

    fscanf(fp, "%lf ", &t);
    fscanf(fp, "%d ", &N1);
    fscanf(fp, "%d ", &N2);
    fscanf(fp, "%d ", &N3);
    fscanf(fp, "%d ", &N1tot);
    fscanf(fp, "%d ", &N2tot);
    fscanf(fp, "%d ", &N3tot);
    fscanf(fp, "%d ", &N1G);
    fscanf(fp, "%d ", &N2G);
    fscanf(fp, "%d ", &N3G);
    fscanf(fp, "%lf ", &startx[1]);
    fscanf(fp, "%lf ", &startx[2]);
    fscanf(fp, "%lf ", &startx[3]);
    fscanf(fp, "%lf ", &dx[1]);
    fscanf(fp, "%lf ", &dx[2]);
    fscanf(fp, "%lf ", &dx[3]);
    fscanf(fp, "%lf ", &tf);
    fscanf(fp, "%d ", &nstep);
    fscanf(fp, "%lf ", &a);
    fscanf(fp, "%lf ", &gam);
    fscanf(fp, "%lf ", &cour);
    fscanf(fp, "%lf ", &DTd);
    fscanf(fp, "%lf ", &DTl);
    fscanf(fp, "%lf ", &DTi);
    fscanf(fp, "%lf ", &DTr);
    fscanf(fp, "%d ", &DTr01);
    fscanf(fp, "%d ", &dump_cnt);
    fscanf(fp, "%d ", &image_cnt);
    fscanf(fp, "%d ", &rdump_cnt);
    fscanf(fp, "%d ", &rdump01_cnt);
    fscanf(fp, "%lf ", &dt);
    fscanf(fp, "%d ", &lim);
    fscanf(fp, "%d ", &failed);
    fscanf(fp, "%lf ", &Rin);
    fscanf(fp, "%lf ", &Rout);
    fscanf(fp, "%lf ", &hslope);
    fscanf(fp, "%lf ", &R0);
    fscanf(fp, "%d ", &NPR);
    fscanf(fp, "%d ", &DOKTOT);
    if(CYL) fscanf(fp, "%d ", &DOCYLINDRIFYCOORDS);
    fscanf(fp, "%lf ", &fractheta);
    fscanf(fp, "%lf ", &fracphi);
    fscanf(fp, "%lf ", &rbr);
    fscanf(fp, "%lf ", &npow2);
    fscanf(fp, "%lf ", &cpow2);
    if(CYL) fscanf(fp, "%lf ", &global_x10);
    if(CYL) fscanf(fp, "%lf ", &global_x20);
    if(CYL) fscanf(fp, "%lf ", &fracdisk);
    if(CYL) fscanf(fp, "%lf ", &fracjet);
    if(CYL) fscanf(fp, "%lf ", &r0disk);
    if(CYL) fscanf(fp, "%lf ", &rdiskend);
    if(CYL) fscanf(fp, "%lf ", &r0jet);
    if(CYL) fscanf(fp, "%lf ", &rjetend);
    if(CYL) fscanf(fp, "%lf ", &jetnu);
    if(CYL) fscanf(fp, "%lf ", &rsjet);
    if(CYL) fscanf(fp, "%lf ", &r0grid);
    fscanf(fp, "%d ", &BL);

    // finish reading out the line
    fseek(fp, 0, SEEK_SET);
    while ((i=fgetc(fp)) != '\n');

    stopx[0] = 1.;
	stopx[1] = startx[1] + N1 * dx[1];
	stopx[2] = startx[2] + N2 * dx[2];
	stopx[3] = startx[3] + N3 * dx[3];
    x1br = log(rbr - R0);

    // initialize malloc for variables
    init_storage_data();

    printf("Simulation grid: %d x %d x %d\n", N1, N2, N3);

    // get data
    printf("Reading body...\n");
    for (i = 0; i < N1; i++) {
        for (j = 0; j < N2; j++) {
            for (k = 0; k < N3; k++) {

                fread(var, sizeof(float), 50, fp);
                //ti  = var[0];
                //tj  = var[1];
                //tk  = var[2];
                //X1  = var[3];
                //X2  = var[4];
                //X3  = var[5];
                a_r[i][j][k]   = var[6];
                a_th[i][j][k]  = var[7];
                a_phi[i][j][k] = var[8];
                if (N3 == 1) a_phi[i][j][k] = 2.*M_PI;

                // primitive variables
                rho[i][j][k] = var[9];
                ug  = var[10];
                V[1][i][j][k] = var[11];
                V[2][i][j][k] = var[12];
                V[3][i][j][k] = var[13];
                B[1][i][j][k] = var[14];
                B[2][i][j][k] = var[15];
                B[3][i][j][k] = var[16];

                //ktot = var[17];
                //divb = var[18];

                // 4-vec U and B (cov and con)
                ucov[0][i][j][k] = var[19];
                ucov[1][i][j][k] = var[20];
                ucov[2][i][j][k] = var[21];
                ucov[3][i][j][k] = var[22];
                ucon[0][i][j][k] = var[23];
                ucon[1][i][j][k] = var[24];
                ucon[2][i][j][k] = var[25];
                ucon[3][i][j][k] = var[26];
                bcov[0][i][j][k] = var[27];
                bcov[1][i][j][k] = var[28];
                bcov[2][i][j][k] = var[29];
                bcov[3][i][j][k] = var[30];
                bcon[0][i][j][k] = var[31];
                bcon[1][i][j][k] = var[32];
                bcon[2][i][j][k] = var[33];
                bcon[3][i][j][k] = var[34];

                //v1min = var[35];
                //v1max = var[36];
                //v2min = var[37];
                //v3min = var[38];
                //v2max = var[39];
                //v3max = var[40];

                // metric determinant
                gdet = var[41];

                // current
                jcov[0][i][j][k] = var[42];
                jcov[1][i][j][k] = var[43];
                jcov[2][i][j][k] = var[44];
                jcov[3][i][j][k] = var[45];
                jcon[0][i][j][k] = var[46];
                jcon[1][i][j][k] = var[47];
                jcon[2][i][j][k] = var[48];
                jcon[3][i][j][k] = var[49];

                // get some derived quantities
                jsq = jcov[0][i][j][k]*jcon[0][i][j][k] +
                      jcov[1][i][j][k]*jcon[1][i][j][k] +
                      jcov[2][i][j][k]*jcon[2][i][j][k] +
                      jcov[3][i][j][k]*jcon[3][i][j][k];
                jdotu = jcov[0][i][j][k]*ucon[0][i][j][k] +
                        jcov[1][i][j][k]*ucon[1][i][j][k] +
                        jcov[2][i][j][k]*ucon[2][i][j][k] +
                        jcov[3][i][j][k]*ucon[3][i][j][k];
                Jsq   = jsq + jdotu*jdotu;
                gJsq  = gdet*Jsq;
                J[i][j][k] = sqrt(Jsq);
                if (isnan(J[i][j][k])) J[i][j][k] = 1e-20;
                pg  = (gam - 1.)*ug;
                bsq = bcov[0][i][j][k]*bcon[0][i][j][k] +
                      bcov[1][i][j][k]*bcon[1][i][j][k] +
                      bcov[2][i][j][k]*bcon[2][i][j][k] +
                      bcov[3][i][j][k]*bcon[3][i][j][k];
                betapl[i][j][k] = 2.*pg/bsq;
                sigmaphi[i][j][k] = bcov[3][i][j][k]*bcon[3][i][j][k]/(rho[i][j][k]);
                Sigmaphi[i][j][k] = B[3][i][j][k]*B[3][i][j][k]/(rho[i][j][k]);
            }
        }
    }
    printf("Finished reading fluid data.\n\n");
    fclose(fp);
/*
    // OLD WAY
    printf("Getting derived quantities...\n");
    // get derived quantities
    for (i = 0; i < N1; i++) {
        for (j = 0; j < N2; j++) {
            for (k = 0; k < N3; k++) {
                jsq[i][j][k]   = ju0[i][j][k]*jd0[i][j][k] +
                                 ju1[i][j][k]*jd1[i][j][k] +
                                 ju2[i][j][k]*jd2[i][j][k] +
                                 ju3[i][j][k]*jd3[i][j][k];
                jdotu[i][j][k] = ju0[i][j][k]*ud0[i][j][k] +
                                 ju1[i][j][k]*ud1[i][j][k] +
                                 ju2[i][j][k]*ud2[i][j][k] +
                                 ju3[i][j][k]*ud3[i][j][k];
                Jsq[i][j][k]   = jsq[i][j][k] + jdotu[i][j][k]*jdotu[i][j][k];
                gJsq[i][j][k]  = gdet[i][j][k]*Jsq[i][j][k];
                J[i][j][k]     = sqrt(Jsq[i][j][k]);
                if (isnan(J[i][j][k])) J[i][j][k] = 1e-20;

                pg[i][j][k]    = (gam - 1.)*ug[i][j][k];
                //K[i][j][k]     = pg[i][j][k]*pow(rho[i][j][k], -gam);
                //EF[i][j][k]    = rho[i][j][k] + gam*ug[i][j][k];
                bsq[i][j][k]   = bu0[i][j][k]*bd0[i][j][k] +
                                 bu1[i][j][k]*bd1[i][j][k] +
                                 bu2[i][j][k]*bd2[i][j][k] +
                                 bu3[i][j][k]*bd3[i][j][k];
                //EE[i][j][k]    = bsq[i][j][k] + EF[i][j][k];
                //va2[i][j][k]   = bsq[i][j][k]/EE[i][j][k];
                //cs2[i][j][k]   = gam*(gam - 1.)*ug[i][j][k]/EF[i][j][k];
                //cms2[i][j][k]  = cs2[i][j][k] + va2[i][j][k] -
                //                 cs2[i][j][k]*va2[i][j][k];
                //T[i][j][k]     = (pg[i][j][k]/rho[i][j][k])*918.059;
                //ptot[i][j][k]  = pg[i][j][k] + 0.5*bsq[i][j][k];
                betapl[i][j][k] = 2.*pg[i][j][k]/bsq[i][j][k];
                //ibetapl[i][j][k] = 0.5*bsq[i][j][k]/pg[i][j][k];
                //sigma[i][j][k] = 0.5*bsq[i][j][k]/rho[i][j][k];
                //sigmaphi[i][j][k] = bu3[i][j][k]*bd3[i][j][k]/(rho[i][j][k]);

                //eflem[i][j][k] = bsq[i][j][k]*uu1[i][j][k]*ud0[i][j][k]
                //                - bu1[i][j][k]*bd0[i][j][k];
                //eflma[i][j][k] = (rho[i][j][k]+pg[i][j][k]+ug[i][j][k])
                //                *uu1[i][j][k]*ud0[i][j][k];

                //lflem[i][j][k] = bsq[i][j][k]*uu1[i][j][k]*ud3[i][j][k]
                //                - bu1[i][j][k]*bd3[i][j][k];
                //lflma[i][j][k] = (rho[i][j][k]+pg[i][j][k]+ug[i][j][k])
                //                *uu1[i][j][k]*ud3[i][j][k];
            }
        }
    }
*/
}


void read_gdump(char *fname)
{
    FILE *fp;
    int i, j, k, l, m, n, index;
    double tf, cour, DTd, DTl, DTi, DTr, DTr01, lim, failed;
    int N1tot, N2tot, N3tot, N1G, N2G, N3G;
    int nstep, dump_cnt, image_cnt, rdump_cnt, rdump01_cnt, NPR, DOKTOT;
    double fractheta, fracphi;
    //double mpi_startn[NDIM];
    //double rr, tth, pphi, X1, X2, X3;
    double dxdxp[NDIM][NDIM];

    float var[58];

    printf("Opening gdump file...\n");
    fp = fopen(fname, "rb");
	if (fp == NULL) {
		fprintf(stderr, "gdump not found.\n");
		exit(1);
	} else {
		fprintf(stderr, "Successfully opened gdump file.\n");
	}

    // get standard HARMPI header
    printf("Reading header...\n");

    fscanf(fp, "%lf ", &t);
    fscanf(fp, "%d ", &N1);
    fscanf(fp, "%d ", &N2);
    fscanf(fp, "%d ", &N3);
    fscanf(fp, "%d ", &N1tot);
    fscanf(fp, "%d ", &N2tot);
    fscanf(fp, "%d ", &N3tot);
    fscanf(fp, "%d ", &N1G);
    fscanf(fp, "%d ", &N2G);
    fscanf(fp, "%d ", &N3G);
    fscanf(fp, "%lf ", &startx[1]);
    fscanf(fp, "%lf ", &startx[2]);
    fscanf(fp, "%lf ", &startx[3]);
    fscanf(fp, "%lf ", &dx[1]);
    fscanf(fp, "%lf ", &dx[2]);
    fscanf(fp, "%lf ", &dx[3]);
    fscanf(fp, "%lf ", &tf);
    fscanf(fp, "%d ", &nstep);
    fscanf(fp, "%lf ", &a);
    fscanf(fp, "%lf ", &gam);
    fscanf(fp, "%lf ", &cour);
    fscanf(fp, "%lf ", &DTd);
    fscanf(fp, "%lf ", &DTl);
    fscanf(fp, "%lf ", &DTi);
    fscanf(fp, "%lf ", &DTr);
    fscanf(fp, "%d ", &DTr01);
    fscanf(fp, "%d ", &dump_cnt);
    fscanf(fp, "%d ", &image_cnt);
    fscanf(fp, "%d ", &rdump_cnt);
    fscanf(fp, "%d ", &rdump01_cnt);
    fscanf(fp, "%lf ", &dt);
    fscanf(fp, "%d ", &lim);
    fscanf(fp, "%d ", &failed);
    fscanf(fp, "%lf ", &Rin);
    fscanf(fp, "%lf ", &Rout);
    fscanf(fp, "%lf ", &hslope);
    fscanf(fp, "%lf ", &R0);
    fscanf(fp, "%d ", &NPR);
    fscanf(fp, "%d ", &DOKTOT);
    if(CYL) fscanf(fp, "%d ", &DOCYLINDRIFYCOORDS);
    fscanf(fp, "%lf ", &fractheta);
    fscanf(fp, "%lf ", &fracphi);
    fscanf(fp, "%lf ", &rbr);
    fscanf(fp, "%lf ", &npow2);
    fscanf(fp, "%lf ", &cpow2);
    if(CYL) fscanf(fp, "%lf ", &global_x10);
    if(CYL) fscanf(fp, "%lf ", &global_x20);
    if(CYL) fscanf(fp, "%lf ", &fracdisk);
    if(CYL) fscanf(fp, "%lf ", &fracjet);
    if(CYL) fscanf(fp, "%lf ", &r0disk);
    if(CYL) fscanf(fp, "%lf ", &rdiskend);
    if(CYL) fscanf(fp, "%lf ", &r0jet);
    if(CYL) fscanf(fp, "%lf ", &rjetend);
    if(CYL) fscanf(fp, "%lf ", &jetnu);
    if(CYL) fscanf(fp, "%lf ", &rsjet);
    if(CYL) fscanf(fp, "%lf ", &r0grid);
    fscanf(fp, "%d ", &BL);

    // finish reading out the line
    fseek(fp, 0, SEEK_SET);
    while ((i=fgetc(fp)) != '\n');

    stopx[0] = 1.;
	stopx[1] = startx[1] + N1 * dx[1];
	stopx[2] = startx[2] + N2 * dx[2];
	stopx[3] = startx[3] + N3 * dx[3];
    x1br = log(rbr - R0);

    // initialize malloc for metric and connection coefficients
    init_storage_metric();

    // Read metric
    printf("Reading body...\n");
    for (i = 0; i < N1; i++) {
        for (j = 0; j < N2; j++) {
            for (k = 0; k < N3; k++) {

                fread(var, sizeof(float), 58, fp);

                //mpi_startn[1] = var[0];
                //mpi_startn[2] = var[1];
                //mpi_startn[3] = var[2];
                //X1 = var[3];
                //X2 = var[4];
                //X3 = var[5];
                //rr = var[6];
                //tth = var[7];
                //pphi = var[8];

                // gcov
                index = 9;
                for (m = 0; m < NDIM; m++) {
                    for(n = 0; n < NDIM; n++) {
                        grid_gcov[i][j][k][m][n] = var[index];
                        index++;
                    }
                }

                // gcon
                for (m = 0; m < NDIM; m++) {
                    for(n = 0; n < NDIM; n++) {
                        grid_gcon[i][j][k][m][n] = var[index];
                        index++;
                    }
                }

                // g
                grid_gdet[i][j][k] = var[index];
                index++;

                // dxdxp
                int index = 9;
                for (m = 0; m < NDIM; m++) {
                    for(n = 0; n < NDIM; n++) {
                        dxdxp[m][n] = var[index];
                        index++;
                    }
                }
            }
        }
    }

    printf("Finished reading gdump.\n\n");
    fclose(fp);
}


void read_metric(char *fname)
{
    int i, j, k;
    int l, m, n;
    int index;
    double var[96];
    FILE *fp;

    fp = fopen(fname, "rb");

    init_storage_metric();

    if(fp==NULL) {
      fprintf(stderr,"error opening metricfile\n") ;
      exit(1) ;
    }

    for (i = 0; i < N1; i++) {
        for (j = 0; j < N2; j++) {
            for (k = 0; k < N3; k++) {

                fread(var, sizeof(double), 96, fp);
                index = 0;

                for (m = 0; m < NDIM; m++) {
                    for (n = 0; n < NDIM; n++) {
                        grid_gcov[i][j][k][m][n] = var[index];
                        index++;
                    }
                }

                for (m = 0; m < NDIM; m++) {
                    for (n = 0; n < NDIM; n++) {
                        grid_gcon[i][j][k][m][n] = var[index];
                        index++;
                    }
                }

                for (l = 0; l < NDIM; l++) {
                    for (m = 0; m < NDIM; m++) {
                        for (n = 0; n < NDIM; n++) {
                            grid_conn[i][j][k][l][m][n] = var[index];
                            index++;
                        }
                    }
                }

            }
        }
    }
    printf("Finished reading metric.\n\n");
    fclose(fp);
}


void write_metric(char *fname)
{
    int i, j, k;
    int l, m, n;
    FILE *fp;

    fp = fopen(fname, "wb");

    if(fp==NULL) {
      fprintf(stderr,"error opening metricfile\n") ;
      exit(1) ;
    }

    for (i = 0; i < N1; i++) {
        for (j = 0; j < N2; j++) {
            for (k = 0; k < N3; k++) {

                for (m = 0; m < NDIM; m++) {
                    for (n = 0; n < NDIM; n++) {
                        fwrite(&grid_gcov[i][j][k][m][n], sizeof(double), 1, fp);
                    }
                }

                for (m = 0; m < NDIM; m++) {
                    for (n = 0; n < NDIM; n++) {
                        fwrite(&grid_gcon[i][j][k][m][n], sizeof(double), 1, fp);
                    }
                }

                for (l = 0; l < NDIM; l++) {
                    for (m = 0; m < NDIM; m++) {
                        for (n = 0; n < NDIM; n++) {
                            fwrite(&grid_conn[i][j][k][l][m][n], sizeof(double), 1, fp);
                        }
                    }
                }

            }
        }
    }
    fclose(fp);
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


void read_current_sheets(char *fname, double ***sheets)
{
    // Reads current from a file generated by write_current_sheets.

    FILE *fp;
    int i, j, k;
    double var[1];

    fp = fopen(fname, "rb");

    printf("Reading current sheets from file '%s'...\n", fname);
    for (i = 0; i < N1; i++) {
        for (j = 0; j < N2; j++) {
            for (k = 0; k < N3; k++) {
                fread(var, sizeof(double), 1, fp);
                sheets[i][j][k] = var[0];
            }
        }
    }
    fclose(fp);
}
