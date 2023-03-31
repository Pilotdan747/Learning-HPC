#include <math.h>

#include "helperFuncs.h"

double battin_xi(double x);
double battin_K(double u);

void lambert_battin(struct vector R1, struct vector R2, double dt, double mu, int dir, struct vector Vs[2]);