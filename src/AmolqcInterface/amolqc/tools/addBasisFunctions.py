#!/usr/bin/env python
#
# reads and writes amolqc wf files.
# allow to add STO basis functions with zero coeffs in MOs
#
from numpy import *
import sys,getopt
from Basis import basisset
from Basis import basis
from Basis import pse
#

def readBlock(wfname,keyword):
    block = []
    try: 
        wffile = open(wfname,"r")
        found = False
        while True:
            line = wffile.readline()
            if not line: break
            if keyword in line: 
                found = True
                break
        while found:
            line = wffile.readline()
            assert line,'block not terminated with $end'
            if '$end' in line:
                break
            block.append(line)
        wffile.close()
    except IOError,EOFError:
        print "error reading block"
        raise
    return block

def writeBlock(handle,block,blockName):
    handle.write(blockName+"\n")
    for line in block:
        handle.write(line)
    handle.write("$end\n")

def readBasis(wfname):
    """read (abs) basis from .wf file and return basis set data structure
    and the number of individual basis functions
    """
    thebasis=[]
    try:
        bfidx = 0
        wffile = open(wfname,"r")
        while True:
            line = wffile.readline()
            assert line, '$basis block missing'
            if '$basis' in line:
                break
        while True:
            line = wffile.readline()
            if not line:
                break
            if '$end' in line:
                break
            words = line.split()
            elem = words[0]
            Z = pse.Z(elem)
            atom = basisset.AtomBasisSet(Z)
            if verbose: print "now read basis functions for ",elem,Z
            atom.read(wffile)
            atombfidx=[]
            for bf in atom.basis:
                atombfidx.append(bfidx)
                if "S" in bf.typ:
                    bfidx += 1
                elif "P" in bf.typ:
                    bfidx += 3
                elif "D" in bf.typ:
                    bfidx +=6
                elif "F" in bf.typ:
                    bfidx += 10
            thebasis.append([elem,atom,atombfidx])
        wffile.close()
    except IOError,EOFError:
        print "error reading basis"
        raise
    return thebasis,bfidx


def countNBF(thebasis):
    nbf1 = 0
    for atomEntry in thebasis:
        atom = atomEntry[1]
        for bf in atom.basis:
            if "S" in bf.typ:
                nbf1 += 1
            elif "P" in bf.typ:
                nbf1 += 3
            elif "D" in bf.typ:
                nbf1 +=6
            elif "F" in bf.typ:
                nbf1 += 10
    aidx = thebasis[-1][2]
    nbf2 = aidx[-1]
    return nbf1, nbf2


def addBasis(thebasis,newbf,nuc,baspos):
    """add (=insert) new basis function newbf at position baspos 
    in basis set of atom nuc as given by full basis set thebasis
    """
    assert nuc<len(thebasis),"addBasis: atom no available in basis"
    elem,atom,aidx = thebasis[nuc]
    atom.basis.insert(baspos,newbf)
    ii = aidx[baspos]
    aidx.insert(baspos,ii)
    for i in range(baspos+1,len(aidx)):
        if "S" in newbf.typ:
            aidx[i] += 1
        elif "P" in newbf.typ:
            aidx[i] += 3
        elif "D" in newbf.typ:
            aidx[i] += 6
        elif "F" in newbf.typ:
            aidx[i] += 10
    for elem,atom,aidx in thebasis[nuc+1:]:
        for i in range(len(aidx)):
            if "S" in newbf.typ:
                aidx[i] += 1
            elif "P" in newbf.typ:
                aidx[i] += 3
            elif "D" in newbf.typ:
                aidx[i] += 6
            elif "F" in newbf.typ:
                aidx[i] += 10
    return thebasis

def writeBasis(handle,thebasis):
    handle.write("$basis\n")
    for elem,atom,aidx in thebasis:
        atom.writeAmolqc(newwf)
    newwf.write("$end\n")

def readMOs(wfname,nbf):
    """read MOs from wf file, for nbf individual basis functions
    """
    mos = []
    try: 
        nmos = 0
        wffile = open(wfname,"r")
        while True:
            line = wffile.readline()
            assert line, '$mos block missing'
            if '$mos' in line:
                break
        while True:
            line = wffile.readline()
            if not line:
                break
            if '$end' in line:
                break
            words = line.split()
            nmos = int(words[0])
            line = wffile.readline()
            for i in range(1,nmos+1):
                line = wffile.readline()
                words = line.split()
                ii = int(words[0])
                assert i==ii,"readCMO: form error in mos in wf file"
                cmo = []
                nlines = nbf/5
                if (nbf%5)>0: nlines += 1
                for k in range(nlines):
                    line = wffile.readline()
                    for l in range(min(5,nbf-k*5)):
                        s = line[l*15:(l+1)*15]
                        s1 = s.replace("D","E")
                        cmo.append(float(s1))
                mos.append(cmo)
        wffile.close()
    except IOError,EOFError:
        print "error reading basis"
        raise
    return mos

def addBFinMOs(mos,bfidx,n):
    """add n basis functions entries into basis at idx bfidx.
    Add zeros in mos.
    """
    newmos = []
    for moi in mos:
        for i in range(n):
            moi.insert(bfidx,0.0)
        newmos.append(moi)
    return newmos

def writeMOs(handle,mos):
    handle.write("$mos\n")
    handle.write(str(len(mos))+"\n")
    handle.write(" \n")
    nbf = len(mos[1])
    for i in range(len(mos)):
        handle.write("   "+str(i+1)+"\n")
        nlines = nbf/5
        if (nbf%5)>0: nlines += 1
        for k in range(nlines):
            line = ""
            for l in range(min(5,nbf-k*5)):
                s = '{0:15.8e}'.format(mos[i][5*k+l])
                s1 = s.replace("e","D")
                line += s1
            handle.write(line+"\n")
    handle.write("$end\n")


###########################################

wfin = 't1.wf'
wfout = 'tnn.wf'

verbose = False

genBlock = readBlock(wfin,'$general')
geomBlock = readBlock(wfin,'$geom')
thebasis,nbf = readBasis(wfin)
jasBlock = readBlock(wfin,'$jastrow')
mos = readMOs(wfin,nbf)
detBlock = readBlock(wfin,'$csfs')


print "done reading wf"
print "# of individ. basis functions: ",nbf
nbf1,nbf2 = countNBF(thebasis)
print "nbf1/2=",nbf1,nbf2

newbf = basis.STO('2P',7.777)
atomnr = 1
atombf = 7
newbasis = addBasis(thebasis,newbf,atomnr,atombf)
bfidx = thebasis[atomnr][2][atombf]
newmos = addBFinMOs(mos,bfidx,3)
thebasis = newbasis
mos = newmos

nbf1,nbf2 = countNBF(thebasis)
print "nbf1/2=",nbf1,nbf2

atomnr = 0
atombf = 7
newbasis = addBasis(thebasis,newbf,atomnr,atombf)
bfidx = thebasis[atomnr][2][atombf]
newmos = addBFinMOs(mos,bfidx,3)
thebasis = newbasis
mos = newmos

nbf1,nbf2 = countNBF(thebasis)
print "nbf1/2=",nbf1,nbf2

newwf = open(wfout,"w")
writeBlock(newwf,genBlock,"$general")
writeBlock(newwf,geomBlock,"$geom")
writeBasis(newwf,thebasis)
writeBlock(newwf,jasBlock,"$jastrow")
writeMOs(newwf,mos)
writeBlock(newwf,detBlock,"$csfs")
newwf.close()





