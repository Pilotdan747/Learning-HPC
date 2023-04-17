from helperFuncs import julianDate, planets_SV_JD
from lambert import lambert_battin
from math import sqrt, pi
from numpy import linspace
from numpy import matlib as MAT
from numpy import linalg as LA

muSun = 1.327124e11
aEarth = 1.49598e8
aMars = 2.27939e8

dim1 = 1000
dim2 = 1001

epoch = julianDate(2023, 0, 0, 0, 0, 0)

TE = 2*pi/sqrt(muSun)*aEarth**(3/2)
TM = 2*pi/sqrt(muSun)*aMars**(3/2)
Tsynodic = 1/abs(1/TE - 1/TM)

launchT = linspace(0, Tsynodic, dim1)/24/3600
transT = linspace(90, 200, dim2)

VinfE = MAT.empty((dim1, dim2))
VinfM = MAT.empty((dim1, dim2))
with open('VinfEPython.txt', 'w') as out:
    for i, lT in enumerate(launchT):
        for j, tT in enumerate(transT):
            REarth, VEarth = planets_SV_JD(3, lT + epoch)
            RMars, VMars = planets_SV_JD(4, tT + lT + epoch)

            if (i == 2 and j == 8):
                check = 1

            V1, V2 = lambert_battin(REarth, RMars, tT*24*3600, muSun, 0)

            VinfE[i, j] = LA.norm(VEarth - V1)
            VinfM[i, j] = LA.norm(VMars - V2)

            out.write('{EarthVinf},'.format(EarthVinf = VinfE[i, j]))
        out.write('\n')
    out.close()
