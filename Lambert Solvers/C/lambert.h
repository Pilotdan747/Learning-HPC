#include <math.h>
#include <stdio.h>
#include <complex.h>

#include "helperFuncs.h"

double complex battin_xi(double complex x);
double complex battin_K(double complex u);

void lambert_battin(struct vector R1, struct vector R2, double dt, double mu, int dir, struct vector Vs[2]);