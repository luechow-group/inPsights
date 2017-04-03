#!/usr/bin/env python
import sys,os 
import subprocess
import numpy as np
import scipy.optimize
import basisset
import basis

def runGaussianAndGetEnergy(baseName):
    """ runGaussianAndGetEnergy: 
    run g09 with baseName.com and return first SCF energy
    This file may contain a local basis set (@file.gbs)""" 
    cmd = ['g09',baseName+'.com']
    subprocess.check_call(cmd)
    logfile = open(baseName+'.log',"r")
    while True:
        line = logfile.readline()
        if 'SCF Done' in line:
            words = line.split()
            energy = float(words[4])
            break
    return energy

class scfEnergy:
    """construct GTO basis set mixed basis set, expand STO to GTOs
    currently only 1S STO expansion, run gaussian and return SCF energy"""

    def __init__(self,amolqcBasisFile,bfList,gaussianBaseName,cvtFile,cvtkey=0):
        comments,bs = basisset.readAmolqcBasisFile(amolqcBasisFile)
        self.bs = bs
        self.bfList = bfList # list containing element and index of basis function for element
        self.baseName = gaussianBaseName
        self.cvt = basis.STOConverter(cvtFile)
        self.cvtkey = cvtkey
        keys = []
        for bf in self.bs:
            keys.append(bf[0])
        # generate list of ptrs to those basis functions that
        # are to be modified
        self.ptrList = []
        for bf in bfList:
            elem,k = bf
            idx = keys.index(elem)
            atombs = self.bs[idx][1]
            self.ptrList.append([atombs,k])
            
    def __call__(self,x):
        """call Gaussian and return SCF energy"""
        assert len(x)==len(self.ptrList),"scfEnergy: vector length must match basis list"
        for p,zeta in zip(self.ptrList,x):
            atombs,idx = p
            if self.cvtkey==0:
                sto = basis.STO('1S',zeta,self.cvt)
            else:
                sto = basis.STO('1S',zeta,self.cvt,self.cvtkey)
            atombs.changeBasisFunction(idx,sto)
        basisset.writeGBSBasisFile(self.bs,'opt.gbs')       
        energy = runGaussianAndGetEnergy(self.baseName)
        return energy

if __name__ == "__main__":
    print "optimizing STO orbital coefficients"
    x0 = np.array([8.0])
    bfList = [['O',0]]

    scf = scfEnergy('6-311Gdp.abs',bfList,'optump2','sto2gto.dat','HUZ10')
    e1 = scf(x0)
    print e1

    xopt = scipy.optimize.fmin(scf,x0,maxiter=150,ftol=1.e-5)
    e1 = scf(xopt)
    print xopt


