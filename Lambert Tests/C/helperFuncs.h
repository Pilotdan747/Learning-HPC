//
// Created by Daniel Owen on 2019-05-15.
//

#ifndef C___HELPERFUNCS_H
#define C___HELPERFUNCS_H

struct vector {
    double x;
    double y;
    double z;
};

#ifndef printFlag
int printFlag = 1;
#endif

double norm(struct vector v);

struct vector cross(struct vector a, struct vector b);
double dot(struct vector a, struct vector b);

struct vector vinf(struct vector v, struct vector vPlanet);

void lambert(struct vector R1, struct vector R2, double dT, double mu, int k, struct vector V[2]);

double julianDate(int y, int m, int d, int hr, int min, int s);

void planets_SV_JD(int id, double JD, struct vector stateOut[2]);

void perifocal2Inertial(double RAAN, double inc, double omega, struct vector Rin, struct vector Vin, struct vector stateOut[2]);

#endif //C___HELPERFUNCS_H
