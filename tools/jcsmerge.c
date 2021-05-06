//gcc jcsmerge.c -o jcsmerge
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int main(int argc, char *argv[])
{
	FILE *fp1, *fp2, *fp3;
	int i, j, k;
	double var1[50], var2[50];
	double J_cs_char_disc, J_cs_char_jet;
	double J_cs_char[1024][512][1];
	char DISCNAME[100] = "/work/gustavo/gyst/jcs_files_disc/";
	char JETNAME[100]  = "/work/gustavo/gyst/jcs_files/";
	//char JCSNAME[100]  = "/work/gustavo/grmonty_harmpi/jcs_files/";
	char JCSNAME[100]  = "/work/gustavo/jcs_files_merged/";
	char dumpname[100];

	sscanf(argv[1], "%s", &dumpname); //dumpXXX

    strcat(JETNAME, dumpname);
    strcat(JETNAME, "_jcs_char");

    strcat(DISCNAME, dumpname);
    strcat(DISCNAME, "_jcs_char");

    strcat(JCSNAME, dumpname);
    strcat(JCSNAME, "_jcs_char");

	fp1 = fopen(DISCNAME, "rb");
	fp2 = fopen(JETNAME, "rb");
	fp3 = fopen(JCSNAME, "wb");

	for (i = 0; i < 1024; i++) {
		for (j = 0; j < 512; j++) {
			for (k = 0; k < 1; k++) {
				fread(var1, sizeof(double), 1, fp1);
				J_cs_char_disc = var1[0];

				fread(var2, sizeof(double), 1, fp2);
				J_cs_char_jet = var2[0];

				J_cs_char[i][j][k] = J_cs_char_disc + J_cs_char_jet;
				if (J_cs_char[i][j][k] > 1.) {
					J_cs_char[i][j][k] = 1.;
				}

                fwrite(&J_cs_char[i][j][k], sizeof(double), 1, fp3);
			}
		}
	}

	return 0;
}

