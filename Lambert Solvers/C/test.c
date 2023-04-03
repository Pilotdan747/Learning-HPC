#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "helperFuncs.h"
#include "lambert.h"

#define sec2day 1/24.0/3600.0
#define day2sec 24*3600

#define DIM1 1000
#define DIM2 1001

int main() {
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

    struct vector earthState[2];
    struct vector marsState[2];

    int i = 519; int j = 1000;

    planets_SV_JD(3, launchT[i] + epoch, earthState);
    REarth = earthState[0];
    VEarth = earthState[1];

    planets_SV_JD(4, transT[j] + launchT[i] + epoch, marsState);
    RMars = marsState[0];
    VMars = marsState[1];

    lambert_battin(REarth, RMars, transT[j]*day2sec, muSun, 0, Vs);

    double VinfE = norm(vinf(VEarth, Vs[0]));
    double VinfM = norm(vinf(VMars, Vs[1]));

    printf("VinfE: %f\n", VinfE);
    printf("VinfM: %f\n", VinfM);

    return 0;
}