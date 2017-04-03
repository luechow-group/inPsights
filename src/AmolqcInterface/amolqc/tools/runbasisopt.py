#!/usr/bin/env python
import numpy as np
import scipy.optimize
from Basis import basisoptg09

# this changes the first basis function for the given element
# to an STO with given zeta

zeta = 9.4
elem = 'Ne'
g09comFile = 'optump2'     # base name of com file. Must contain basis as 'opt.gbs'
absFile = '6-311Gdp.abs'   # original basis set in amolqc basis set format
key = 'OPT14'

x0 = np.array([zeta])
bfList = [[elem,0]]
scf = basisoptg09.scfEnergy(absFile,bfList,g09comFile,'sto2gto.dat',key)
e1 = scf(x0)

print "optimizing STO orbital coefficient for ",elem
print " expanding the STO with ",key

print " Start energy with zeta = ",zeta,": E = ",e1
print " optimizing (running g09) ..."
xopt = scipy.optimize.fmin(scf,x0,maxiter=150,ftol=1.e-5)
e1 = scf(xopt)
print " Optimized zeta = ",xopt[0]," with E = ",e1
print " "

