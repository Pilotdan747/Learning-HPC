#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "helperFuncs.h"
#include "lambert.h"

#define sec2day 1/24.0/3600.0
#define day2sec 24*3600

#define DIM1 1000
#define DIM2 1001

int main() {
    double progStart = omp_get_wtime();

    double muSun, aEarth, aMars, TE, TM, Tsynodic;

    muSun = 1.327124e11;
    aEarth = 1.49598e8;
    aMars = 2.27939e8;

    double launchT[DIM1];
    double transT[DIM2];

    double epoch = julianDate(2023, 0, 0, 0, 0, 0);

    TE = 2*M_PI/sqrt(muSun)*pow(aEarth, 3.0/2.0);
    TM = 2*M_PI/sqrt(muSun)*pow(aMars, 3.0/2.0);
    Tsynodic = 1/fabs(1/TE - 1/TM);

    for (int i = 0; i < DIM1; i++) {
        launchT[i] = (double)i/(double)(DIM1-1)*Tsynodic*sec2day;
    }

    for (int i = 0; i < DIM2; i++) {
        transT[i] = (double)i/(double)(DIM2-1)*(200 - 90) + 90;
    }

    struct vector REarth, VEarth, RMars, VMars;
    struct vector Vs[2];

    double **VinfE = (double**) malloc(sizeof(double)*DIM1);
    double **VinfM = (double**) malloc(sizeof(double)*DIM1);

    for (int i = 0; i <DIM1; i++) {
        VinfE[i] = (double*) malloc(sizeof(double)*DIM2);
        VinfM[i] = (double*) malloc(sizeof(double)*DIM2);
    }

    struct vector earthState[2];
    struct vector marsState[2];
    
    double start = omp_get_wtime();

    int printFlag = 1;
    for (int i = 0; i < DIM1; i++) {
        planets_SV_JD(3, launchT[i] + epoch, earthState);
        REarth = earthState[0];
        VEarth = earthState[1];

        for (int j = 0; j < DIM2; j++) {
            planets_SV_JD(4, transT[j] + launchT[i] + epoch, marsState);
            RMars = marsState[0];
            VMars = marsState[1];

            lambert_battin(REarth, RMars, transT[j]*day2sec, muSun, 0, Vs);

            VinfE[i][j] = norm(vinf(VEarth, Vs[0]));
            VinfM[i][j] = norm(vinf(VMars, Vs[1]));

            if (i == 0 && j == 0) {
                printf("VinfE: %f\n", VinfE[i][j]);
                printf("REarth X: %f Y: %f Z: %f\n", REarth.x, REarth.y, REarth.z);
                printf("VEarth X: %f Y: %f Z: %f\n", VEarth.x, VEarth.y, VEarth.z);
                printf("V1 X: %f Y: %f Z: %f\n", Vs[0].x, Vs[0].y, Vs[0].z);
            }
        }
    }

    double end = omp_get_wtime();

    printf("Time to Run: %f seconds\n", end - start);

    FILE *out;
    out = fopen("SerialOut.txt", "w");

    for (int i = 0; i < DIM1; i++) {
        for (int j = 0; j < DIM2; j++) {
            if (j != DIM2 - 1)
                fprintf(out, "%f,", VinfE[i][j]);
            else
                fprintf(out, "%f\n", VinfE[i][j]);
        }
    }

    fclose(out);

    double progEnd = omp_get_wtime();
    printf("Time to Run Prog: %f seconds\n", progEnd - progStart);

    return 0;
}