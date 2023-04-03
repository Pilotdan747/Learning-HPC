#include "lambert.h"

void lambert_battin(struct vector R1, struct vector R2, double dt, double mu, int dir, struct vector Vs[2]) {
    double r1 = norm(R1); 
    double r2 = norm(R2);
    double kk = cross(R1,R2).z;
    double th = acos(dot(R1,R2)/r1/r2);
    if ((dir == 0 && kk < 0) || (dir > 0 && kk >= 0))
        th = 2*M_PI - th;

    double c = sqrt(pow(r1, 2.0) + pow(r2, 2.0) - 2*r1*r2*cos(th));
    double s = 0.5*(r1 + r2 + c);

    double L = sqrt((s-c)/s);
    if (th > M_PI)
        L = -L; 

    double T = sqrt(8*mu/pow(s, 3.0))*dt;

    double r0p = 0.25*s*pow(1+L, 2.0);
    double l = pow((1-L)/(1+L), 2.0);
    double m = pow(T, 2.0)/pow(1+L, 6.0);

    double Tp = 4.0/3.0*pow(1-L, 3.0);
    double x0 = 0;
    if (T <= Tp)
        x0 = 0;
    else 
        x0 = l;
    double x = x0;

    double complex z, den, h1, h2, B, u, K, y;
    for (int i = 0; i < 100; i++) {
        z = battin_xi(x);
        den = (1 + 2*x + l)*(4*x + z*(3+x));
        h1 = cpow(l+x, 2.0)*(1+3*x+z)/den;
        h2 = m*(x - l + z)/den;
        B = 0.25*27*h2/cpow(1+h1, 3.0);
        u = 0.5*B/(1 + csqrt(1+B));
        K = battin_K(u);
        y = (1+h1)/3.0*(2 + csqrt(1+B)/(1+2*u*cpow(K, 2.0)));
        x = csqrt(0.25*cpow(1-l, 2.0)+m/cpow(y, 2.0)) - 0.5*(1+l);
        if (cabs(x-x0) < 1e-14)
            break;
        else
            x0 = x;
    }

    double a = mu*pow(dt, 2.0)/16.0/pow(r0p, 2.0)/creal(x)/pow(creal(y), 2.0);

    double b, amin, tmin, ae, dE, f, g, gdot, ah, bh, dH;
    if (a > 0) {
        b = 2*asin(sqrt(0.5*(s-c)/a));
        if (th > M_PI)
            b = -b;

        amin = 0.5*s;
        tmin = sqrt(pow(amin, 3.0)/mu)*(M_PI - b + sin(b));
        ae = 2*asin(sqrt(0.5*s/a));
        if (dt > tmin)
            ae = 2*M_PI - ae;

        dE = ae - b;
        f = 1 - a/r1*(1 - cos(dE));
        g = dt - sqrt(pow(a, 3.0)/mu)*(dE - sin(dE));
        gdot = 1 - a/r2*(1 - cos(dE));
    } else {
        ah = 2*asinh(sqrt(-0.5*s/a));
        bh = 2*asinh(sqrt(-0.5*(s-c)/a));
        dH = ah - bh;
        f = 1 - a/r1*(1 - cosh(dH));
        g = dt - sqrt(-pow(a, 3.0)/mu)*(sinh(dH) - dH);
        gdot = 1 - a/r2*(1 - cosh(dH));
    }

    struct vector V1, V2;
    V1.x = (R2.x - f*R1.x)/g; V1.y = (R2.y - f*R1.y)/g; V1.z = (R2.z - f*R1.z)/g;
    V2.x = (gdot*R2.x - R1.x)/g; V2.y = (gdot*R2.y - R1.y)/g;  V2.z = (gdot*R2.z - R1.z)/g; 

    Vs[0] = V1;
    Vs[1] = V2;
}


double complex battin_xi(double complex x) {
    double tiny = 1e-30;
    double complex d = sqrt(1+x)+1;
    double complex n = x/(d*d);

    double complex f0 = tiny; double complex C0 = f0; double complex D0 = 0;

    //stage 1
    double complex D = 3 + 8*d*D0;
    if (cabs(D) < tiny) 
        D = tiny;
    double complex C = 3 + 8*d/C0;
    if (cabs(C) < tiny) 
        C = tiny; 
    D = 1/D;  double complex Del = C*D;
    double complex f = f0*Del;
    f0 = f; C0 = C; D0 = D;

    // stage 2
    D = 5+n + 1*D0;
    if (cabs(D) < tiny) 
        D = tiny;
    C = 5+n + 1/C0;
    if (cabs(C) < tiny)
        C = tiny;
    D = 1/D;  Del = C*D;
    f = f0*Del;
    f0 = f; C0 = C; D0 = D;

    // stage 3
    D = 1 + 9.0/7.0*n*D0;
    if (cabs(D) < tiny)
        D = tiny;
    C = 1 + 9.0/7.0*n/C0;
    if (cabs(C) < tiny)
        C = tiny;
    D = 1/D;  Del = C*D;
    f = f0*Del;
    f0 = f; C0 = C; D0 = D;

    double complex c;
    for (int i = 1; i <= 100; i++) {
        c = cpow(i+3, 2.0)/(cpow(2*(i+3), 2.0)-1);
        D = 1 + c*n*D0;
        if (cabs(D) < tiny)
            D = tiny;
        C = 1 + c*n/C0;
        if (cabs(C) < tiny)
            C = tiny;
        D = 1/D;  Del = C*D;
        f = f0*Del;
        if (cabs(Del-1) < 1e-14)
            break;
        else
            f0 = f; C0 = C; D0 = D;
    }

    return f;
}

double complex battin_K(double complex u) {
    double tiny = 1e-30;

    double complex f0 = tiny; double complex C0 = f0; double complex D0 = 0;

    // stage 1
    double complex D = 1 + 1.0/3.0*D0;
    if (cabs(D) < tiny)
        D = tiny;
    double complex C = 1 + 1.0/3.0/C0;
    if (cabs(C) < tiny)
        C = tiny;
    D = 1/D; double complex Del = C*D;
    double complex f = f0*Del;
    f0 = f; C0 = C; D0 = D;

    // stage 2
    D = 1 + 4.0/27.0*u*D0;
    if (cabs(D) < tiny)
        D = tiny;
    C = 1 + 4.0/27.0*u/C0;
    if (cabs(C) < tiny)
        C = tiny;
    D = 1/D;  Del = C*D;
    f = f0*Del;
    f0 = f; C0 = C; D0 = D;

    double complex c1, c2;
    for (int i = 1; i <= 100; i++) {
        c1 = 2.0*(3.0*i+1.0)*(6.0*i-1.0)/9.0/(4.0*i-1.0)/(4.0*i+1.0);
        c2 = 2.0*(3.0*i+2.0)*(6.0*i+1.0)/9.0/(4.0*i+1.0)/(4.0*i+3.0);
        D = 1 + c1*u*D0;
        if (cabs(D) < tiny)
            D = tiny;
        C = 1 + c1*u/C0;
        if (cabs(C) < tiny)
            C = tiny;
        D = 1/D;  Del = C*D;
        f = f0*Del;
        
        f0 = f; C0 = C; D0 = D;
        D = 1 + c2*u*D0;
        if (cabs(D) < tiny)
            D = tiny;
        C = 1 + c2*u/C0;
        if (cabs(C) < tiny)
            C = tiny;
        D = 1/D;  Del = C*D;
        f = f0*Del;

        if (cabs(Del-1) < 1e-14)
            break;
        else
            f0 = f; C0 = C; D0 = D;
    }

    return f;
}
   

