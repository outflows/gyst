#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "definitions.h"
#include "declarations.h"

// grmonty units
void assign_units() {
    /*
     * Calculates the units for radiative transfer calculations.
     * Not needed here.
    */

    printf("Assigning units...\n");

    L_unit = GNEWT * MBH / (CL * CL);
    T_unit = L_unit / CL;
    RHO_unit = MUNIT / (L_unit*L_unit*L_unit);
    U_unit = RHO_unit * CL * CL;
    B_unit = CL * sqrt(4. * M_PI * RHO_unit);
    Ne_unit = RHO_unit / (MP + ME);

    printf("MBH: %g\n", MBH);
    printf("M_unit: %g\n", MUNIT);
    printf("L_unit: %g\n", L_unit);
    printf("RHO_unit: %g\n", RHO_unit);
    printf("U_unit: %g\n", U_unit);
    printf("B_unit: %g\n", B_unit);
    printf("Ne_unit: %g\n\n", Ne_unit);
}

void read_data(char *fname) {
    /*
     * Reads the data from an iharm2d dump file, stores these data in the
     * appropriate variables and uses them to calculate a few derived
     * quantities.
     *
     * This function is largely based on init_harm2dv3_data.c from grmonty,
     * with a few ideas taken from Gammie's file to read iharm2d data when
     * plotting.
    */

    FILE *fp;
    int i, j, k;
    //double startx[NDIM];//, dx[NDIM];
    double tf, cour, DTd, DTl, DTi, DTr, DTr01, lim, failed;
    int N1tot, N2tot, N3tot, N1G, N2G, N3G;
    int nstep, dump_cnt, image_cnt, rdump_cnt, rdump01_cnt, NPR, DOKTOT;
    double fractheta, fracphi;
    // cylindrification stuff
    //double x10, x20, fracdisk, fracjet, r0disk, rdiskend, r0jet, rjetend;
    //double jetnu, rsjet, r0grid;

    //double ti, tj, tk, x1, x2, x3;//, r, th, phi;
    //double ktot, v1min, v1max, v2min, v2max, v3min, v3max;

    float var[50]; // must be float because HARMPI data is saved as float!

    printf("Opening dump file '%s'...\n", fname);
    fp = fopen(fname, "r");
    if (fp == NULL) {
		fprintf(stderr, "dump file not found.\n");
		exit(1);
	} else {
		fprintf(stderr, "Successfully opened '%s'.\n", fname);
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

    //printf("Done!\n");
    init_storage();

    printf("Simulation grid: %d x %d x %d\n", N1, N2, N3);

    // get data
    printf("Reading body...\n");
    for (i = 0; i < N1; i++) {
        for (j = 0; j < N2; j++) {
            for (k = 0; k < N3; k++) {

                // float, not double, because HARMPI data is saved as float
                fread(var, sizeof(float), 50, fp);
                //ti  = var[0];
                //tj  = var[1];
                //tk  = var[2];
                //X1  = var[3];
                //X2  = var[4];
                //X3  = var[5];
                r[i][j][k]   = var[6];
                th[i][j][k]  = var[7];
                phi[i][j][k] = var[8];

                // primitive variables
                rho = var[9];
                ug  = var[10];
                v1  = var[11];
                v2  = var[12];
                v3  = var[13];
                B1  = var[14];
                B2  = var[15];
                B3[i][j][k] = var[16];

                //ktot = var[17];
                //divb = var[18];

                // 4-vec U and B (con and cov)
                //uu0 = var[19];
                //uu1 = var[20];
                //uu2 = var[21];
                //uu3 = var[22];
                ud0 = var[23];
                ud1 = var[24];
                ud2 = var[25];
                ud3 = var[26];
                bu0 = var[27];
                bu1 = var[28];
                bu2 = var[29];
                bu3[i][j][k] = var[30];
                bd0 = var[31];
                bd1 = var[32];
                bd2 = var[33];
                bd3 = var[34];

                //v1min = var[35];
                //v1max = var[36];
                //v2min = var[37];
                //v3min = var[38];
                //v2max = var[39];
                //v3max = var[40];

                // metric determinant
                gdet = var[41];

                // current; added for my version of HARMPI which calculates currents
                ju0 = var[42];
                ju1 = var[43];
                ju2 = var[44];
                ju3 = var[45];
                jd0 = var[46];
                jd1 = var[47];
                jd2 = var[48];
                jd3 = var[49];

                // get some derived quantities
                jsq   = ju0*jd0 + ju1*jd1 + ju2*jd2 + ju3*jd3;
                jdotu = ju0*ud0 + ju1*ud1 + ju2*ud2 + ju3*ud3;
                Jsq   = jsq + jdotu*jdotu;
                gJsq  = gdet*Jsq;
                J[i][j][k] = sqrt(Jsq);
                if (isnan(J[i][j][k])) J[i][j][k] = 1e-20;
                pg  = (gam - 1.)*ug;
                bsq = bu0*bd0 + bu1*bd1 + bu2*bd2 + bu3[i][j][k]*bd3;
                betapl[i][j][k] = 2.*pg/bsq;
                sigmaphi[i][j][k] = bu3[i][j][k]*bd3/(rho);
                Sigmaphi[i][j][k] = B3[i][j][k]*B3[i][j][k]/(rho);
            }
        }
    }

    //printf("Done!\n");
    printf("Finished reading data.\n\n");
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
    double mpi_startn[NDIM];
    double rr, tth, pphi, X1, X2, X3;
    double dxdxp[NDIM][NDIM];

    float var[58];

    printf("Opening gdump file '%s'...\n", fname);
    fp = fopen(fname, "rb");
	if (fp == NULL) {
		fprintf(stderr, "gdump not found.\n");
		exit(1);
	} else {
		fprintf(stderr, "Successfully opened '%s'.\n", fname);
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

    // Read metric
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
                        gcov[m][n] = var[index];
                        index++;
                    }
                }

                // gcon
                for (m = 0; m < NDIM; m++) {
                    for(n = 0; n < NDIM; n++) {
                        gcon[m][n] = var[index];
                        index++;
                    }
                }

                // g
                gdet = var[index];
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


double *malloc_rank1(int n1, int size)
{
    /*
        Allocates a 1D array
    */

	void *A;

	if ((A = (void *) malloc(n1 * size)) == NULL) {
		fprintf(stderr, "malloc failure in malloc_rank1\n");
		exit(123);
	}

	return A;
}

double **malloc_rank2_cont(int n1, int n2)
{
    // Allocates a multidimensional array

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
    // Allocates a multidimensional array

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


void init_storage(void)
{
    printf("Allocating memory...\n\n");

    flag = (int ***) malloc_rank3_cont(N1, N2, N3);
    flag_buffer = (int ***) malloc_rank3_cont(N1, N2, N3);
    r = (double ***) malloc_rank3_cont(N1, N2, N3);
    th = (double ***) malloc_rank3_cont(N1, N2, N3);
    phi = (double ***) malloc_rank3_cont(N1, N2, N3);
    B3 = (double ***) malloc_rank3_cont(N1, N2, N3);
    bu3 = (double ***) malloc_rank3_cont(N1, N2, N3);
    J = (double ***) malloc_rank3_cont(N1, N2, N3);
    J_cs = (double ***) malloc_rank3_cont(N1, N2, N3);
    betapl = (double ***) malloc_rank3_cont(N1, N2, N3);
    sigmaphi = (double ***) malloc_rank3_cont(N1, N2, N3);
    Sigmaphi = (double ***) malloc_rank3_cont(N1, N2, N3);

}

void init_storage_metric()
{
    printf("Allocating memory for grid...\n\n");
}
