//
// Created by Daniel Owen on 2019-05-15.
//

#include "helperFuncs.h"
#include <math.h>
#include <stdio.h>

#define pi 4*atan(1)

double J2000_elements[9][6] = {{0.38709927,0.20563593,7.00497902,48.33076593,77.45779628,252.25032350},  
{0.72333566,0.00677672,3.39467605,76.67984255,131.60246718,181.97909950}, 
{1.00000261,0.01671123,-0.00001531,0.0,102.93768193,100.46457166},
{1.52371034,0.09339410,1.84969142,49.55953891,-23.94362959,-4.55343205},
{5.20288700,0.04838624,1.30439695,100.47390909,14.72847983,34.39644501},
{9.53667594,0.05386179,2.48599187,113.66242448,92.59887831,49.95424423},
{19.18916464,0.04725744,0.77263783,74.01692503,170.95427630,313.23810451},
{30.06992276,0.00859048,1.77004347,131.78422574,44.96476227,-55.12002969},
{39.48211675,0.24882730,17.14001206,110.30393684,224.06891629,238.92903833}};

double cent_rates[9][6] = {{0.00000037,0.00001906,-0.00594749,-0.12534081,0.16047689,149472.67411175},  
{0.00000390,-0.00004107,-0.00078890,-0.27769418,0.00268329,58517.81538729},
{0.00000562,-0.00004392,-0.01294668,0.0,0.32327364,35999.37244981},
{0.0001847,0.00007882,-0.00813131,-0.29257343,0.44441088,19140.30268499},
{-0.00011607,-0.00013253,-0.00183714,0.20469106,0.21252668,3034.74612775},
{-0.00125060,-0.00050991,0.00193609,-0.28867794,-0.41897216,1222.49362201},
{-0.00196176,-0.00004397,-0.00242939,0.04240589,0.40805281,428.48202785},
{0.00026291,0.00005105,0.00035372,-0.00508664,-0.32241464,218.45945325},
{-0.00031596,0.00005170,0.00004818,-0.01183482,-0.04062942,145.20780515}}; 

double norm(struct vector v) {
    return sqrt(pow(v.x, 2) + pow(v.y, 2) + pow(v.z, 2));
}

struct vector vinf(struct vector v, struct vector vPlanet) {
    struct vector vInf;

    vInf.x = v.x - vPlanet.x;
    vInf.y = v.y - vPlanet.y;
    vInf.z = v.z - vPlanet.z;

    return vInf;
}

struct vector cross(struct vector a, struct vector b) {
    struct vector c;

    c.x = a.y * b.z - a.z * b.y;
    c.y = a.z * b.x - a.x * b.z;
    c.z = a.x * b.y - a.y * b.x;

    return c;
}

double dot(struct vector a, struct vector b) {
    double c = 0;

    c += a.x*b.x;
    c += a.y*b.y;
    c += a.z*b.z;

    return c;
}

void lambert(struct vector R1, struct vector R2, double dT, double mu, int k, struct vector V[2]) {

    double r1 = norm(R1);
    double r2 = norm(R2);

    struct vector c = cross(R1, R2);

    double dTheta;

    if (k > 0) {
        if (c.z > 0) {
            dTheta = acos(dot(R1, R2) / (r1 * r2));
        } else {
            dTheta = 2 * pi - acos(dot(R1, R2) / (r1 * r2));
        }
    } else {
        if (c.z > 0) {
            dTheta = 2 * pi - acos(dot(R1, R2) / (r1 * r2));
        } else {
            dTheta = acos(dot(R1, R2) / (r1 * r2));
        }
    }
    double A = sin(dTheta)*sqrt(r1*r2/(1 - cos(dTheta)));

    double Znew, Zold, dZ, C, S, y, F, FPri;
    Zold = pi;
    int count = 0;
    int biSect = 0;
    do {
        C = 0.5 - Zold/24 + pow(Zold, 2)/720;
        S = 1.0/6.0 - Zold/120 + pow(Zold, 2)/5040;
        y = r1 + r2 + A*(Zold*S - 1)/sqrt(C);

        F = pow(y/C, 3.0/2.0)*S + A*sqrt(y) - sqrt(mu)*dT;
        FPri = pow(y/C, 3.0/2.0)*(1/(2*Zold)*(C - 3*S/(2*C)) + 3*pow(S, 2)/(4*C)) + A/8*(3*S/C*sqrt(y) + A*sqrt(C/y));

        Znew = Zold - F/FPri;
        dZ = fabs(Znew - Zold);
        Zold = Znew;
        count++;
        if (count > 1000) {
            biSect = 1;
        }
    } while (dZ > 0.00001 && count < 1001);

    double Zl, Zu, Z;

    if (biSect) {
        Zl = -4*pi;
        Zu = 4*pow(pi, 2);
        Z = 0;

        while (fabs(F) >= 0.001) {
            C = 0.5 - Z/24 + pow(Z, 2)/720 - pow(Z, 3)/40320;
            S = 1.0/6.0 - Z/120 + pow(Z, 2)/5040 - pow(Z, 3)/362880;
            y = r1 + r2 + A*(Z*S - 1)/sqrt(C);

            F = pow(y/C, 3.0/2.0)*S + A*sqrt(y) - sqrt(mu)*dT;

            if (F > 0) {
                Zu = Z;
            } else {
                Zl = Z;
            }

            Z = (Zl + Zu)/2;
        }
    }

    y = r1 + r2 + A*(Zold*S - 1)/sqrt(C);

    double f, g, gDot;

    f = 1 - y/r1;
    g = A*sqrt(y/mu);
    gDot = 1 - y/r2;

    struct vector v1;
    v1.x = 1/g*(R2.x - f*R1.x); v1.y = 1/g*(R2.y - f*R1.y); v1.z = 1/g*(R2.z - f*R1.z);

    struct vector v2;
    v2.x = 1/g*(gDot*R2.x - R1.x); v2.y = 1/g*(gDot*R2.y - R1.y); v2.z = 1/g*(gDot*R2.z - R1.z);

    V[0] = v1;
    V[1] = v2;
}

double julianDate(int y, int m, int d, int hr, int min, int s) {
    double yy, mm, dd, hrhr, minmin, ss;
    yy = (double) y; mm = (double) m; dd = (double) d; hrhr = (double) hr; minmin = (double) min; ss = (double) s;

    double j0 = 367*yy-floor(7*(yy+floor((mm+9)/12))/4)+floor(275*mm/9)+dd + 1721013.5;
    double ut = hrhr + minmin/60 + ss/3600;
    return j0 + ut/24;
}

void planets_SV_JD(int id, double JD, struct vector stateOut[2]) {
    double muSun = 1.327124e11;
    
    double T = (JD - 2451545)/(double)36525;

    double elements[6];

    for (int i = 0; i < 6; i++) {
        elements[i] = J2000_elements[id-1][i] + cent_rates[id-1][i]*T;
    }

    double AU = 149597870.7;
    elements[0] = elements[0]*AU;
    elements[2] = elements[2]*M_PI/180; elements[3] = elements[3]*M_PI/180; elements[4] = elements[4]*M_PI/180; elements[5] = elements[5]*M_PI/180; 

    double a = elements[0];
    double e = elements[1];
    double h = sqrt(muSun*a*(1-pow(e,2)));
    double inc = elements[2];
    double RAAN = fmod(elements[3], 2*M_PI);
    double wHat = fmod(elements[4], 2*M_PI);
    double L = fmod(elements[5], 2*M_PI);
    double omega = fmod(wHat - RAAN, 2*M_PI);
    double M = fmod(L - wHat, 2*M_PI); 

    double tol = 1e-12;
    int maxIter = 1000;
    
    double E = M;

    double f, df, E_New;
    for (int i = 0; i < maxIter; i++) {
        f = E - e*sin(E) - M;
        df = 1 - e*cos(E);

        E_New = E - f/df;

        if (fabs(E - E_New) < tol) {
            E = E_New;
            break;
        }

        E = E_New;
    }

    double theta = fmod(2*atan(sqrt((1 + e)/(1 - e))*tan(E/2)), 2*M_PI);
    double r = (pow(h, 2)/muSun)/(1 + e*cos(theta));

    struct vector rp, vp;
    rp.x = r*cos(theta); rp.y = r*sin(theta); rp.z = 0;
    vp.x = -muSun/h*sin(theta); vp.y = muSun/h*(e + cos(theta)); vp.z = 0;

    perifocal2Inertial(RAAN, inc, omega, rp, vp, stateOut);
}

void perifocal2Inertial(double RAAN, double inc, double omega, struct vector Rin, struct vector Vin, struct vector stateOut[2]) {
    struct vector R1, R2, R3, V1, V2, V3;

    R1.x = cos(omega)*Rin.x + -sin(omega)*Rin.y; R1.y = sin(omega)*Rin.x + cos(omega)*Rin.y; R1.z = Rin.z;
    V1.x = cos(omega)*Vin.x + -sin(omega)*Vin.y; V1.y = sin(omega)*Vin.x + cos(omega)*Vin.y; V1.z = Vin.z;

    R2.x = R1.x; R2.y = cos(inc)*R1.y + -sin(inc)*R1.z; R2.z = sin(inc)*R1.y + cos(inc)*R1.z;
    V2.x = V1.x; V2.y = cos(inc)*V1.y + -sin(inc)*V1.z; V2.z = sin(inc)*V1.y + cos(inc)*V1.z;

    R3.x = cos(RAAN)*R2.x + -sin(RAAN)*R2.y; R3.y = sin(RAAN)*R2.x + cos(RAAN)*R2.y; R3.z = R2.z;
    V3.x = cos(RAAN)*V2.x + -sin(RAAN)*V2.y; V3.y = sin(RAAN)*V2.x + cos(RAAN)*V2.y; V3.z = V2.z;

    stateOut[0] = R3;
    stateOut[1] = V3;
}
