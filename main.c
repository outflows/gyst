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
    FILE *gdumpfile, *metricfile;
    char dump[100], gdump[100];
    char path[100] = "/home/gustavo/Dropbox/trabalho/programs/grmonty_harmpi/";

    clock_t begin = clock();

    if (argc < 2) {
        fprintf(stderr, "ERROR: dumpfile must be specified. To run this program properly, you must type:\n\n");
        fprintf(stderr, "./gyst dumpfile_name\n\n");
        exit(0);
    }

    strcpy(dump, path);
    strcpy(gdump, path);
    sscanf(argv[1], "%s", dump_file);
    strcat(dump, dump_file);
    strcpy(jcs_output, dump_file);
    strcat(jcs_output, "_jcs");
    strcpy(jpeak_output, dump_file);
    strcat(jpeak_output, "_jcs_peak");
    strcpy(jchar_output, dump_file);
    strcat(jchar_output, "_jcs_char");

    if (!CHARACTERIZE) {
        printf("CHARACTERIZE set to 0 by user: this run WILL NOT characterize the current sheets\n\n");
    }
    else {
        printf("CHARACTERIZE set to 1 by user: this run WILL characterize the current sheets\n\n");
    }

    clock_t begin_read = clock();
    assign_units();
    printf("dump file: '%s'\n", dump);
    read_data(dump);

    metricfile = fopen("metric", "rb");
    if (metricfile != NULL) {
        printf("metric file found, reading metric...");
        read_metric("metric");
        fclose(metricfile);
    }
    else {
        strcat(gdump, "gdump2D");
        gdumpfile = fopen(gdump, "rb");
        if (gdumpfile == NULL) {
            printf("gdump not found, setting up grid...\n");

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
            printf("gdump file: '%s'\n", gdump);
            printf("gdump found, reading gdump...\n");
            read_gdump(gdump);
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
        printf("writing metric information into 'metric'\n");
        write_metric("metric");
    }

    printf("Finished data initialization.\n\n");
    clock_t end_read = clock();
    double time_spent_read = (double)(end_read - begin_read) / CLOCKS_PER_SEC;

    clock_t begin_current = clock();
    get_current_sheets();
    if (CHARACTERIZE) {
        characterize();
        write_current_sheets(jchar_output, J_cs_char);
    }
    clock_t end_current = clock();
    double time_spent_current = (double)(end_current - begin_current) / CLOCKS_PER_SEC;

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

    printf("\n");
    printf("Running time (initialization): %lf seconds.\n", time_spent_read);
    printf("Running time (get current sheets): %lf seconds.\n", time_spent_current);
    printf("Total running time: %.2lf seconds.\n", time_spent);
    printf("Total running time: %.2lf minutes.\n", time_spent/60.);

    return(0);
}
