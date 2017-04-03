#
import numpy as np
import sys
#
#
#
def norm2(a):
    return(np.sqrt(np.dot(a,a)))

def distance(a,b):
    return(np.sqrt(np.dot(a-b,a-b)))

def angle(b,a,c):
    return 180.0/np.pi*np.arccos(np.dot(b-a,c-a)/(norm2(b-a)*norm2(c-a)))

def generateRotMatrix(a,b):
    R = np.zeros((3,3))
    a0 = a / norm2(a)
    baproj = a0 * np.dot(a0,b)
    b1 = b - baproj
    b10 = b1 / norm2(b1)
    c0 = np.cross(a0,b10) 
    R[:,0] = a0
    R[:,1] = b10
    R[:,2] = c0
    return R
    
class Molecule:
    """
    class for operations on molecules
    """
    def __init__(self,wfFileName=None):
        assert wfFileName != None
        self.bohr2angs = 0.529177249

        wfFile = open(wfFileName,"r")
        line = ''
        while True:
            line = wfFile.readline()
            if not line: break
            if '$geom' in line: break
        if not line:
            print 'no $geom block found'
            sys.exit(1)

        # read geometry
        line = wfFile.readline()
        words = line.split()
        self.nnuc = int(words[0])
        self.elemList = []
        self.nucCoord = []
        for i in range(self.nnuc):
            line = wfFile.readline()
            words = line.split()
            self.elemList.append(words[0])
            x = float(words[1])
            y = float(words[2])
            z = float(words[3])
            self.nucCoord.append(np.array([x,y,z]))

    def printMolecule(self):
        print self.nnuc
        for i in range(self.nnuc):
            print self.elemList[i],self.nucCoord[i][0],self.nucCoord[i][1],self.nucCoord[i][2]

    def getNNuc(self):
        return self.nnuc

    def getBohr2angs(self):
        return self.bohr2angs

    def translate(self,vec):
        for i in range(self.nnuc):
            self.nucCoord[i] = self.nucCoord[i] + vec

    def rotate(self,R):
        for i in range(self.nnuc):
            atomvec = self.nucCoord[i]
            print "DBG rotate: ",i
            print atomvec
            atomvec = np.dot(R,atomvec)
            print atomvec
            self.nucCoord[i] = atomvec

