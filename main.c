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
    char dump[200], gdump[200], metric[200], jcs_dir[200];

    #if (RUN_DISC)
    {
        strcpy(jcs_dir, JCS_DISC_DIR);
    }
    #elif (RUN_OUTFLOWS)
    {
        strcpy(jcs_dir, JCS_OUTFLOWS_DIR);
    }
    #elif (RUN_JET)
    {
        strcpy(jcs_dir, JCS_JET_DIR);
    }
    #endif

    clock_t begin = clock();

    if (argc < 2) {
        fprintf(stderr, "ERROR: dumpfile must be specified. To run this program properly, you must type:\n\n");
        fprintf(stderr, "./gyst dumpfile_name\n\n");
        exit(0);
    }
    strcpy(dump, DUMP_DIR);
    strcpy(gdump, DUMP_DIR);

    sscanf(argv[1], "%s", dump_file);
    strcat(dump, dump_file);
    strcat(jcs_dir, dump_file);
    strcpy(jcs_output, jcs_dir);
    strcat(jcs_output, "_jcs");
    //strcpy(jcs_output, dump_file);
    //strcpy(jpeak_output, dump_file);
    strcpy(jpeak_output, jcs_dir);
    strcat(jpeak_output, "_jcs_peak");
    //strcpy(jchar_output, dump_file);
    strcpy(jchar_output, jcs_dir);
    strcat(jchar_output, "_jcs_char");

    strcpy(metric, DUMP_DIR);
    strcat(metric, "metric");
    //strcpy(metric, dump_file);
    //strcat(metric, "_metric");

    if (SHEETS) {
        fprintf(stdout, "SHEETS set to 1 by user: this run WILL find the current sheets\n");
    }
    else {
        fprintf(stdout, "SHEETS set to 0 by user: this run WILL NOT find the current sheets\n");
    }

    if (CHARACTERIZE) {
        fprintf(stdout, "CHARACTERIZE set to 1 by user: this run WILL characterize the current sheets\n\n");
        if (RUN_DISC) {
            fprintf(stdout, "RUN_DISC set to 1\n\n");
        }
        else if (RUN_OUTFLOWS) {
            fprintf(stdout, "RUN_OUTFLOWS set to 1\n\n");
        }
        else if (RUN_JET) {
            fprintf(stdout, "RUN_JET set to 1\n\n");
        }

    }
    else {
        printf(stdout, "CHARACTERIZE set to 0 by user: this run WILL NOT characterize the current sheets\n\n");
    }

    clock_t begin_read = clock();
    assign_units();
    fprintf(stdout, "dump file: '%s'\n", dump);
    read_data(dump);

    metricfile = fopen(metric, "rb");
    if (metricfile != NULL) {
        fprintf(stdout, "metric file: %s\n", metric);
        fprintf(stdout, "metric file found, reading metric...\n");
        fclose(metricfile);
        read_metric(metric);
    }
    else {
        strcat(gdump, "gdump");
        gdumpfile = fopen(gdump, "rb");
        if (gdumpfile == NULL) {
            fprintf(stdout, "gdump not found, setting up grid...\n");

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
            fprintf(stdout, "gdump file: '%s'\n", gdump);
            fprintf(stdout, "gdump found, reading gdump...\n");
            read_gdump(gdump);
        }

        fprintf(stdout, "Calculating connection coefficients...\n");
        for (i = 0; i < N1; i++) {
            for (j = 0; j < N2; j++) {
                for (k = 0; k < N3; k++) {
                    coord(i, j, k, X);
                    get_geometry(i, j, k, &geom);
                    get_connection(X, &geom, grid_conn, i, j, k);
                }
            }
        }
        fprintf(stdout, "writing metric information into 'metric'\n");
        write_metric(metric);
    }

    fprintf(stdout, "Finished data initialization.\n\n");
    clock_t end_read = clock();
    double time_spent_read = (double)(end_read - begin_read) / CLOCKS_PER_SEC;

    clock_t begin_current = clock();
    if (SHEETS) {
        get_current_sheets();
    }
    clock_t end_current = clock();
    double time_spent_current = (double)(end_current - begin_current) / CLOCKS_PER_SEC;

    clock_t begin_char = clock();
    if (CHARACTERIZE) {
        if (!SHEETS) {
            read_current_sheets(jcs_output, J_cs);
            read_current_sheets(jpeak_output, J_cs_peak);
        }
        if (N3 > 1) characterize(); // not doing 3D here, this is just legacy code!
        else characterize();
        //else characterize2D_single();
    }
    clock_t end_char = clock();
    double time_spent_char = (double)(end_char - begin_char) / CLOCKS_PER_SEC;

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

    fprintf(stdout, "\n");
    fprintf(stdout, "Running time (initialization): %.2lf seconds.\n", time_spent_read);
    if (SHEETS) fprintf(stdout, "Running time (get current sheets): %.2lf seconds.\n", time_spent_current);
    if (CHARACTERIZE) fprintf(stdout, "Running time (characterization): %.2lf seconds.\n", time_spent_char);
    fprintf(stdout, "Total running time: %.2lf seconds.\n", time_spent);
    fprintf(stdout, "Total running time: %.2lf minutes.\n", time_spent/60.);

    return(0);
}
