// Welcome to GYST (pronounced with a hard g, as in Get Your Sheets Together)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "definitions.h"
#include "declarations.h"

int main(int argc, char *argv[])
{
    int i, j, k;
    double X[NDIM];
    struct of_geom geom;
    //FILE *gdumpfile;
    char dump[100], gdump[100];
    char path[100] = "/home/gustavo/Dropbox/trabalho/programs/grmonty_harmpi/";

    clock_t begin = clock();

    if (argc < 2) {
		fprintf(stderr, "ERROR: dumpfile must be specified. To run this program properly, you must type:\n\n");
        fprintf(stderr, "./gyst dumpfile_name\n\n");
        //fprintf(stderr, "where dumpfile_name is your dumpfile.\n");
		exit(0);
    }

    init_storage();
    assign_units();

    strcpy(dump, path);
    strcpy(gdump, path);

    sscanf(argv[1], "%s", dump_file);
    strcat(dump, dump_file);
    printf("Path to dump file: '%s'\n", dump);

    strcpy(jcs_output, dump_file);
    strcat(jcs_output, "_jcs");
    printf("Current sheets output file: '%s'\n", jcs_output);

    strcpy(jpeak_output, dump_file);
    strcat(jpeak_output, "_jcs_peak");
    printf("Peak current output file: '%s'\n\n", jpeak_output);

    strcpy(normal_output, dump_file);
    strcat(normal_output, "_normal");

    if (!CHARACTERIZE) {
        printf("CHARACTERIZE set to 0 by user: this run WILL NOT characterize the current sheets\n\n");
    }
    else {
        printf("CHARACTERIZE set to 1 by user: this run WILL characterize the current sheets\n\n");
    }

    clock_t begin_read = clock();

    read_data(dump);

    strcpy(gdump, path);
    if (N3 > 1) {
        strcat(gdump, "gdump3D");
    }
    else {
        strcat(gdump, "gdump2D");
    }
    FILE *metricfile, *gdumpfile;
    char metricfilename[30];
    gdumpfile = fopen(gdump,"rb");
    if (N3 > 1) {
        metricfile = fopen("metricfile3D","rb");
    }
    else {
        metricfile = fopen("metricfile2D","rb");
    }

    if (gdumpfile == NULL) {
        printf("gdump not found.\n");
        if (metricfile == NULL) {
            printf("metricfile not found.\n");
            exit(1);
        }
        else {
            fclose(metricfile); // we'll open it in function
            printf("metricfile found, reading metric...\n");
            if (N3 > 1) {
                strcpy(metricfilename, "metricfile3D");
            }
            else {
                strcpy(metricfilename, "metricfile2D");
            }
            read_metric(metricfilename);
        }
    } else if (metricfile == NULL) {
        printf("gdump file found.\n");
        printf("metricfile not found.\n");

        fclose(gdumpfile); // we'll open it in function
        printf("Reading gdump...\n");
        read_gdump(gdump);

        printf("Getting Christoffel symbols...\n");
        for (i = 0; i < N1; i++) {
            for (j = 0; j < N2; j++) {
                for (k = 0; k < N3; k++) {
                    coord(i, j, k, X);
                    //get_connection(X, grid_gcon, grid_conn);
                }
            }
        }

        if (N3 > 1) {
            strcpy(metricfilename, "metricfile3D");
        }
        else {
            strcpy(metricfilename, "metricfile2D");
        }
        printf("Writing metric and Christoffel symbols into %s...\n", metricfilename);
        write_metric(metricfilename); // write metric and christoffels to file so we won't have to do it again
    }
    else {
        fclose(gdumpfile);
        fclose(metricfile); // we'll open it in function
        printf("gdump file found, but not used.\n");
        printf("metricfile found.\n");

        if (N3 > 1) {
            strcpy(metricfilename, "metricfile3D");
        }
        else {
            strcpy(metricfilename, "metricfile2D");
        }
        printf("Reading metricfile...\n");
        read_metric(metricfilename);
    }
    printf("Finished reading metric file(s).\n\n");



/*
    strcat(gdump, "gdump2D");
    printf("gdump file: '%s'\n", gdump);
    gdumpfile = fopen(gdump, "rb");
    if (gdumpfile == NULL) {
        printf("gdump not found, preparing grid...\n");

        for (i = 0; i < N1; i++) {
            for (j = 0; j < N2; j++) {
                for (k = 0; k < N3; k++) {
                    coord(i, j, k, X);
                    gcov_func(X, grid_gcov[i][j][k]);
                    grid_gdet[i][j][k] = gdet_func(grid_gcov[i][j][k]);
                    gcon_func(grid_gcov[i][j][k], grid_gcon[i][j][k]);
                }
            }
        }
    }
    else {
        fclose(gdumpfile);
        printf("gdump found, reading gdump...\n");
        read_gdump(gdump); // will give me grid_gcov, grid_gcon, grid_gdet
    }

    printf("Calculating connection coefficients...\n");
    for (i = 0; i < N1; i++) {
        for (j = 0; j < N2; j++) {
            for (k = 0; k < N3; k++) {
                get_geometry(i, j, k, &geom);
                get_connection(X, &geom, grid_conn[i][j][k]);
            }
        }
    }
*/

    clock_t end_read = clock();
    double time_spent_read = (double)(end_read - begin_read) / CLOCKS_PER_SEC;

    clock_t begin_current = clock();
    get_current_sheets();
    write_current_sheets(jcs_output);
    //if (CHARACTERIZE) find_normal(normal_output);
    if (CHARACTERIZE) find_normal();
    clock_t end_current = clock();
    double time_spent_current = (double)(end_current - begin_current) / CLOCKS_PER_SEC;

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

    printf("Running time (initialization): %lf seconds.\n", time_spent_read);
    printf("Running time (get current sheets): %lf seconds.\n", time_spent_current);
    printf("Total running time: %lf seconds.\n", time_spent);
    printf("Total running time: %.2lf minutes.\n", time_spent/60.);

    return(0);
}
