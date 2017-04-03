#
import numpy as np
import sys
#
#
#
pse = ['X','H','He',
       'Li','Be','B','C','N','O','F','Ne',
       'Na','Mg','Al','Si','P','S','Cl','Ar',
       'K','Ca',
           'Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn',
                'Ga','Ge','As','Se','Br','Kr']

def elem(idx):
    assert idx>0 and idx<len(pse),"elem: internal pse too small"
    return pse[idx]

def Z(elem):
    assert elem in pse,"Z: element is not in internal pse"
    return pse.index(elem)

def norm2(a):
    return(np.sqrt(np.dot(a,a)))

def distance(a,b):
    return(np.sqrt(np.dot(a-b,a-b)))

def angle(b,a,c):
    return 180.0/np.pi*np.arccos(np.dot(b-a,c-a)/(norm2(b-a)*norm2(c-a)))

        
    
class Walkers:
    """
    class working on nuclei und reference coordinates
    allowing several analyses such as distances and
    a characterization of the reference for use in
    energy partitioning
    """
    def __init__(self,refFileName=None,charge=0,mult=1):
        assert refFileName != None
        self.bohr2angs = 0.529177249
        self.nElec = charge

        refFile = open(refFileName,"r")
        # read geometry
        line = refFile.readline()
        words = line.split()
        self.nnuc = int(words[0])
        self.elemList = []
        self.nucCoord = []
        for i in range(self.nnuc):
            line = refFile.readline()
            words = line.split()
            self.nElec += Z(words[1])
            self.elemList.append(words[1])
            x = float(words[2])
            y = float(words[3])
            z = float(words[4])
            self.nucCoord.append(np.array([x,y,z]))
        # read next line
        line = refFile.readline()
        self.fh = refFile
        self.nAlpha = (self.nElec + mult-1) / 2

    def getFileHandle(self):
        return self.fh

    def readNextWalker(self):
        coords = []
        n = 0
        counter = 0
        line = self.fh.readline()
        if line:
            words = line.split()
            n = int(words[1])
            counter = int(words[2])
            line = self.fh.readline()
            words = line.split()
            nelec = int(words[0])
            coords = []
            for j in range(nelec):
                line = self.fh.readline()
                words = line.split()
                coords.append(np.array([float(words[0]),float(words[1]),float(words[2])]))
            if self.nElec==0: 
                self.nElec = nelec
        return n,counter,coords            

    def readNextWalkerWithSpin(self,nalpha):
        """ ref format with electron coordinates and index (e.g. xyz files) as
            x y z ii   where ii is the original electron index, i.e. alpha when ii <= nalpha
            return 1 for alpha and -1 for beta spin 
        """
        coords = []
        spin = []
        n = 0
        counter = 0
        line = self.fh.readline()
        if line:
            words = line.split()
            n = int(words[1])
            counter = int(words[2])
            line = self.fh.readline()
            words = line.split()
            nelec = int(words[0])
            coords = []
            for j in range(nelec):
                line = self.fh.readline()
                words = line.split()
                coords.append(np.array([float(words[0]),float(words[1]),float(words[2])]))
                if int(words[3]) <= self.nAlpha:
                    spin.append(1)
                else:
                    spin.append(-1)
            if self.nElec == 0:
                self.nElec = nelec
        return n,counter,coords,spin            

    def getNMax(self):
        return self.nMax

    def getNElec(self):
        return self.nElec

    def getNAlpha(self):
        return self.nAlpha

    def getNNuc(self):
        return self.nnuc

    def getBohr2angs(self):
        return self.bohr2angs

    def getElecNucDistances(self,coords):
        Dist = np.zeros((self.nElec,self.nnuc))
        for i in range(self.nElec):
            for j in range(self.nnuc):
                Dist[i,j] = distance(coords[i],self.nucCoord[j])
        return Dist

    def getElecNucSortedDistances(self,coords):
        Dist = []
        for i in range(self.nElec):
            l = []
            for j in range(self.nnuc):
                d = distance(coords[i],self.nucCoord[j])
                l.append((j,d))
            Dist.append(sorted(l,key=lambda pair: pair[1]))
        return Dist

    def getElecDistances(self,coords):
        Dist = np.zeros((self.nElec,self.nElec))
        for i in range(self.nElec):
            for j in range(i+1,self.nElec):
                Dist[i,j] = distance(coords[i],coords[j])
                Dist[j,i] = Dist[i,j]
        return Dist

    def getNucDistances(self):
        Dist = np.zeros((self.nnuc,self.nnuc))
        for i in range(self.nnuc):
            for j in range(i+1,self.nnuc):
                Dist[i,j] = distance(self.nucCoord[i],self.nucCoord[j])
                Dist[j,i] = Dist[i,j]
        return Dist

    def getDistance(self,a,b,coords):
        """ a,b strings of the form "[e|n]123". Initial e refers to electron
        initial n to nucleus. Indices are standard (1,2..) not Python (from 0)"""
        assert a[0:1]=="e" or a[0:1]=="n" and b[0:1]=="e" or b[0:1]=="n"
        if a[0:1]=="e":
            ra = coords[int(a[1:])-1]
        else:
            ra = self.nucCoord[int(a[1:])-1]
        if b[0:1]=="e":
            rb = coords[int(b[1:])-1]
        else:
            rb = self.nucCoord[int(b[1:])-1]
        return distance(ra,rb)
                    

    def getAngle(self,b,a,c,coords):
        """ b,a,c strings of the form "[e|n]123". Initial e refers to electron
        initial n to nucleus. Return angle b,a,c in degree and the two distances
        b-a and c-a"""
        assert a[0:1]=="e" or a[0:1]=="n" and b[0:1]=="e" or b[0:1]=="n" and c[0:1]=="e" or c[0:1]=="n"
        if a[0:1]=="e":
            ra = coords[int(a[1:])-1]
        else:
            ra = self.nucCoord[int(a[1:])-1]
        if b[0:1]=="e":
            rb = coords[int(b[1:])-1]
        else:
            rb = self.nucCoord[int(b[1:])-1]
        if c[0:1]=="e":
            rc = coords[int(c[1:])-1]
        else:
            rc = self.nucCoord[int(c[1:])-1]
        dab = distance(ra,rb)
        dac = distance(ra,rc)
        alp = angle(rb,ra,rc)
        return (alp,dab,dac)

        
