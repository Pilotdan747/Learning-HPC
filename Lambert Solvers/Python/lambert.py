from numpy import linalg as LA
from numpy import cross, dot
from math import sqrt, cos, sin, acos, asin, cosh, sinh, asinh, pi

def battin_xi(x):
    pass

def battin_K(u):
    pass

def lambert_battin(R1, R2, dt, mu, dir):
    r1 = LA.norm(R1)
    r2 = LA.norm(R2)
    kk = cross(R1, R2)[3]
    th = acos(dot(R1, R2)/r1/r2)

    if ((dir == 0 and kk < 0) or (dir > 0 and kk >= 0)):
        th = 2*pi - th

    c = sqrt(r1**2 + r2**2 - 2*r1*r2*cos(th))
    s = 0.5*(r1 + r2 + c)

    L = sqrt((s - c)/s)

    if (th > pi):
        L = -L

    T = sqrt(8*mu/s**3)*dt

    r0p = 0.25*s*(1 + L)**2
    l = ((1 - L)/(1 + L))**2
    m = T**2/(1 + L)**6

    Tp = 4/3*(1 - L)**3
    if (T <= Tp):
        x0 = 0
    else:
        x0 = l
    x = x0

    for i in range(0, 100):
        z = battin_xi(x)
        den = (1 + 2*x + l)*(4*x + x*(3+x))
        h1 = (l + x)**2*(1 + 3*x + z)/den
        h2 = m*(x - l + z)/den
        B = 0.25*27*h2/(1 + h1)**3
        u = 0.5*B/(1 + sqrt(1 + B))
        K = battin_K(u)
        y = (1 + h1)/3*(2 + sqrt(1 + B)/(1 + 2*u*K**2))
        x = sqrt(0.25*(1 - l)**2 + m/y**2) - 0.5*(1 + l)
        if (abs(x - x0) < 1e-14):
            break
        else:
            x0 = x

    a = mu*dt**2/16/r0p**2/x/y**2

    if (a > 0):
        b = 2*asin(sqrt(0.5*(s-c)/a))
        if (th > pi):
            b = -b

        amin = 0.5*s
        tmin = sqrt(amin**3/mu)*(pi - b + sin(b))
        ae = 2*asin(sqrt(0.5*s/a))
        if (dt > tmin):
            ae = 2*pi - ae
        
        dE = ae - B
        f = 1 - a/r1*(1 - cos(dE))
        g = dt - sqrt(a**3/mu)*(dE - sin(dE))
        gdot = 1 - a/r2*(1 - cos(dE))
    else:
        ah = 2*asinh(sqrt(-0.5*s/a))
        bh = 2*asinh(sqrt(-0.5*(s - c)/a))
        dH = ah - bh
        f = 1 - a/r1*(1 - cosh(dH))
        g = dt - sqrt(-1*a**3/mu)*(sinh(dH) - dH)
        gdot = 1 - a/r2*(1 - cosh(dH))
    
    V1 = (R2 - f*R1)/g
    V2 = (gdot*R2 - R1)/g

    return V1, V2

