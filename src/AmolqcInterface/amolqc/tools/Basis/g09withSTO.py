#!/usr/bin/env python
import sys,os 
import subprocess
from numpy import *
from basisset import *


def writeGaussianFile(handle,R,gbsFile):
    handle.write("# gfinput HF/gen int=nobasistransform units=au scf=tight\n")
    handle.write("\n")
    handle.write(" minimal STO basis\n")
    handle.write("\n")
    handle.write("0 1\n")
    handle.write("H  0.0   0.0   "+str(-R/2.0)+"\n")
    handle.write("H  0.0   0.0   "+str(R/2.0)+"\n")
    handle.write("\n")
    handle.write("@/home/hx016lu/Projekte/Eqmc/Basis/"+gbsFile)
    handle.write("\n")
    handle.write("\n")
    handle.write("\n")

def runGaussianAndGetEnergy(baseName):
    """ run g09 with baseName.com and return first SCF energy
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
    """construct basis set from vector, run gaussian and
    return SCF energy"""

    def __init__(self,amolqcBasisFile,bfList,gaussianBaseName,cvtFile,cvtkey=0):
        comments,bs = readAmolqcBasisFile(amolqcBasisFile)
        self.bs = bs
        self.bfList = bfList # list containing element and index of basis function for element
        self.baseName = gaussianBaseName
        self.cvt = STOConverter(cvtFile)
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
                sto = STO('1S',zeta,self.cvt)
            else:
                sto = STO('1S',zeta,self.cvt,self.cvtkey)
            atombs.changeBasisFunction(idx,sto)
        writeGBSBasisFile(self.bs,self.baseName+'.gbs')       
        energy = runGaussianAndGetEnergy(self.baseName)
        print "scf:",x,energy
        return energy

    def value(self,x):
        """call Gaussian and return SCF energy"""
        assert len(x)==len(self.ptrList),"scfEnergy: vector length must match basis list"
        for p,zeta in zip(self.ptrList,x):
            atombs,idx = p
            if self.cvtkey==0:
                sto = STO('1S',zeta,self.cvt)
            else:
                sto = STO('1S',zeta,self.cvt,self.cvtkey)
            atombs.changeBasisFunction(idx,sto)
        writeGBSBasisFile(self.bs,self.baseName+'.gbs')       
        energy = runGaussianAndGetEnergy(self.baseName)
        return energy


def EH2minimalbasis(zeta,R_AB):
    comFile = 'h2sto.com'
    gbsFile = 'h2sto.gbs'
    basename = 'h2sto'
    gfile = open('h2sto.com','w')
    writeGaussianFile(gfile,R_AB,gbsFile)
    gfile.close() 
    x0 = array([zeta])
    bfList = [['H',0]]
    scf = scfEnergy('STO-SZ.abs',bfList,basename,'sto2gto.dat','HUZ10')
    e = scf(x0)
    return e

E = EH2minimalbasis(1.0,1.4)
print "H2 mit zeta=1.0, R=1.4:",E
print "BE = ",(-1.0-E)*2625.5

E = EH2minimalbasis(1.2,1.4)
print "H2 mit zeta=1.2, R=1.4:",E
print "BE = ",(-1.0-E)*2625.5

