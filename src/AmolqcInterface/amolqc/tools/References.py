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

def dict2plist(pdict):
    plist = []
    for key in sorted(pdict):
        plist.append([pdict[key],[],key])
    return plist

def printplist(plist):
    print "["
    for e in plist:
        print e,","
    print "]"

def printcgentry(cglist):
    print "chemgroups=",len(cglist)
    for i,e in enumerate(cglist):
        line = " "+str(i+1)+" "
        line +=  e[0][0:4]+" "
        line +=  str(len(e[2]))+" "+str(len(e[1][0]))+" "+str(len(e[1][1]))+"  "
        for t in e[2]:
            line += str(t)+" "
        line += " "
        for t in e[1][0]:
            line += str(t)+" "
        line += " "
        for t in e[1][1]:
            line += str(t)+" "
        print line
    print ")"

def cglist2enlist(cglist):
    """
    converts chemical group list (from characterizeRef1) into
    the enlist format 
    """
    enlist = []
    for i in range(len(cglist)):
        entry = cglist[i]
        typ = entry[0]
        if typ[:4]=="bdtH":
            elecs = entry[2]
            nucs = []
            nucs.append(entry[1][0][0])
            adjcore = entry[1][1][0] - 1
            nucs.append(cglist[adjcore][1][0][0])
            enlist.append([elecs,nucs,typ])
        elif typ[:4]=="bond":
            elecs = entry[2]
            nucs = []
            adjcore1 = entry[1][1][0] - 1
            adjcore2 = entry[1][1][1] - 1
            nucs.append(cglist[adjcore1][1][0][0])
            nucs.append(cglist[adjcore2][1][0][0])
            enlist.append([elecs,nucs,typ])
        elif typ[:4]=="lone":
            elecs = entry[2]
            nucs = []
            adjcore1 = entry[1][1][0] - 1
            nucs.append(cglist[adjcore1][1][0][0])
            enlist.append([elecs,nucs,typ])
    return enlist
 
def analyseRefSignature(nl,VERBOSE=True):
    nBonds = [0,0,0,0,0,0,0,0,0]
    nCore = 0
    nlp = 0
    for key,elist in nl.iteritems():
        nb = len(elist)
        if 'core' in key:
            if VERBOSE: print "    {0:2d}e {1:10s} :".format(nb,key),elist 
            nCore += nb
    for key,elist in nl.iteritems():
        nb = len(elist)
        if 'bond' in key:
            if nb==1:
                if VERBOSE: print "   half {1:10s} :".format(nb,key),elist 
                nBonds[1] +=  1
            if nb==2:
                if VERBOSE: print " single {1:10s} :".format(nb,key),elist 
                nBonds[2] += 1
            if nb==3:
                if VERBOSE: print "oneHalf {1:10s} :".format(nb,key),elist 
                nBonds[3] += 1
            if nb==4:
                if VERBOSE: print " double {1:10s} :".format(nb,key),elist 
                nBonds[4] += 1
            if nb==5:
                if VERBOSE: print "twoHalf {1:10s} :".format(nb,key),elist 
                nBonds[5] += 1
            if nb==6:
                if VERBOSE: print " triple {1:10s} :".format(nb,key),elist 
                nBonds[6] += 1
    for key,elist in nl.iteritems():
        nb = len(elist)
        if 'lp' in key:
            if VERBOSE: print "    {0:2d}e {1:10s} :".format(nb,key),elist 
            nlp += nb
    if VERBOSE: print "\nlist summary:"
    if nCore > 0:
        if VERBOSE: print " {0:4d} core electrons".format(nCore)
    if nBonds[1] > 0:
        if VERBOSE: print " {0:4d} half bonds".format(nBonds[1])
    if nBonds[2] > 0:
        if VERBOSE: print " {0:4d} single bonds".format(nBonds[2])
    if nBonds[3] > 0:
        if VERBOSE: print " {0:4d} oneHalf bonds".format(nBonds[3])
    if nBonds[4] > 0:
        if VERBOSE: print " {0:4d} double bonds".format(nBonds[4])
    if nBonds[5] > 0:
        if VERBOSE: print " {0:4d} twoHalf bonds".format(nBonds[5])
    if nBonds[6] > 0:
        if VERBOSE: print " {0:4d} triple bonds".format(nBonds[6])
    if nlp > 0:
        if VERBOSE: print " {0:4d} lone pair-like electrons".format(nlp)
    return nCore,nBonds,nlp

    
class References:
    """
    class working on nuclei und reference coordinates
    allowing several analyses such as distances and
    a characterization of the reference for use in
    energy partitioning
    """
    def __init__(self,refFileName=None,strno=1):
        assert refFileName != None
        self.bohr2angs = 0.529177249
        self.NUC_THRESH = 0.01
        self.CORE_THRESH = 2.0
        self.BOND_DIFF = 0.4

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
            self.elemList.append(words[1])
            x = float(words[2])
            y = float(words[3])
            z = float(words[4])
            self.nucCoord.append(np.array([x,y,z]))
        # read # of references
        line = refFile.readline()
        words = line.split()
        self.nMax = int(words[0])
        self.mMax = int(words[3])
        self.refList = (self.nMax+1)*[None]
        for i in range(self.nMax+1):
            self.refList[i] = (self.mMax+1)*[None]
        # read ref list
        while True:
            line = refFile.readline()
            if not line: break
            words = line.split()
            n = int(words[1])
            m = int(words[2])
            if (len(words) >= 5):
                value = float(words[4])
            else:
                value = 0.0
            if (len(words) >= 7):
                found = int(words[6])
                self.nAlpha = int(words[8])
            else:
                found = 0
            line = refFile.readline()
            words = line.split()
            nelec = int(words[0])
            coords = []
            for j in range(nelec):
                line = refFile.readline()
                words = line.split()
                coords.append(np.array([float(words[0]),float(words[1]),float(words[2])]))
            self.refList[n][m]=[value,found,coords]
        self.nElec = nelec
        refFile.close()

    def getMultiplicity(self):
        return 2*self.nAlpha-self.nElec+1

    def getNAlpha(self):
        return self.nAlpha
        
    def getNMax(self):
        return self.nMax

    def getMMax(self):
        return self.mMax

    def DoesExist(self,n,m):
        return self.refList[n][m] != None

    def getNElec(self):
        return self.nElec

    def getNNuc(self):
        return self.nnuc

    def getBohr2angs(self):
        return self.bohr2angs

    def getRefCoords(self,n=1,m=1):
        coords = self.refList[n][m][2]
        return coords
    
    def getNucCoords(self,n=1,m=1):
        coords = self.nucCoord
        return coords

    def getValueAndFreq(self,n=1,m=1):
        value = self.refList[n][m][0]
        freq = self.refList[n][m][1]
        return (value,freq)

    def getElecNucDistances(self,n=1,m=1):
        Dist = np.zeros((self.nElec,self.nnuc))
        coords = self.refList[n][m][2]
        for i in range(self.nElec):
            for j in range(self.nnuc):
                Dist[i,j] = distance(coords[i],self.nucCoord[j])
        return Dist

    def getElecNucSortedDistances(self,n=1,m=1):
        Dist = []
        coords = self.refList[n][m][2]
        for i in range(self.nElec):
            l = []
            for j in range(self.nnuc):
                d = distance(coords[i],self.nucCoord[j])
                l.append((j,d))
            Dist.append(sorted(l,key=lambda pair: pair[1]))
        return Dist

    def getNearestElecSortedDistances(self,n=1,m=1):
        eDist = np.zeros((self.nElec))
        coords = self.refList[n][m][2]
        nElecsDists = [[] for _ in range(self.nnuc)]
        for i in range(self.nElec):
            minidx = 0
            mind = 10000
            for j in range(self.nnuc):
                d = distance(coords[i],self.nucCoord[j])
                if d < mind:
                    minidx = j
                    mind = d
            nElecsDists[minidx].append((i,mind))
        for j in range(self.nnuc):
            nElecsDists[j] = sorted(nElecsDists[j],key=lambda pair: pair[1])
        return nElecsDists

    def getElecDistances(self,n=1,m=1):
        Dist = np.zeros((self.nElec,self.nElec))
        coords = self.refList[n][m][2]
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

    def getDistance(self,a,b,n=1,m=1):
        """ a,b strings of the form "[e|n]123". Initial e refers to electron
        initial n to nucleus. Indices are standard (1,2..) not Python (from 0)"""
        assert a[0:1]=="e" or a[0:1]=="n" and b[0:1]=="e" or b[0:1]=="n"
        coords = self.refList[n][m][2]
        if a[0:1]=="e":
            ra = coords[int(a[1:])-1]
        else:
            ra = self.nucCoord[int(a[1:])-1]
        if b[0:1]=="e":
            rb = coords[int(b[1:])-1]
        else:
            rb = self.nucCoord[int(b[1:])-1]
        return distance(ra,rb)
                    

    def getAngle(self,b,a,c,n=1,m=1):
        """ b,a,c strings of the form "[e|n]123". Initial e refers to electron
        initial n to nucleus. Return angle b,a,c in degree and the two distances
        b-a and c-a"""
        assert a[0:1]=="e" or a[0:1]=="n" and b[0:1]=="e" or b[0:1]=="n" and c[0:1]=="e" or c[0:1]=="n"
        coords = self.refList[n][m][2]
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


    def getNucElecVector(self,j,i,n=1,m=1):
        coords = self.refList[n][m][2]
        relec = coords[i]
        rnuc = self.nucCoord[j]
        return relec - rnuc

        
    def getMidPointRef(self,pairList,n=1,m=1):
        """ pairList with format [ [ [a,b], [nuc] ], ... ]
            as for energies """
        coords = self.refList[n][m][2]
        midPoints = np.zeros((self.nElec,3))
        for i in range(len(coords)):
            midPoints[i,:] = coords[i]
        for p in pairList:
            assert len(p[0])==1 or len(p[0])==2
            if len(p[0])==2:
                i = p[0][0] - 1; j = p[0][1] - 1
                midvec = 0.5*(midPoints[i,0:3] + midPoints[j,0:3])
                midPoints[i,0:3] = midvec
                midPoints[j,0:3] = midvec
        return midPoints

    def writeMidPointRef(self,pairList,n=1,m=1):
        """ pairList with format [ [ [a,b], [nuc] ], ... ]
            as for energies """
        mid = self.getMidPointRef(pairList,n,m)
        value,freq = self.getValueAndFreq(n,m)
        print "MidpointRef: {0:3d} {1:3d} F(ref): {2:10.5f} Found: {3:10d}".format(n,m,value,freq) 
        print self.nElec
        for i in range(self.nElec):
            print " {0:10.5f} {1:10.5f} {2:10.5f}".format(mid[i,0],mid[i,1],mid[i,2])

    def writeNucAndBondDistances(self,enList,n=1,m=1,unit='angs'):
        """ enlist with format [ [e,n,type], ... ], where e and
            n are electron and nuclei lists resp., and type is a string
            nuc list must contain one element (lone pair) or two (bond)
            in the latter case the bond distance, the ref distance to nuclei
            and the distance of the ref position to the bond line is printed
            All distances are printed per default in angstrom or bohr 
        """
        b2a = self.bohr2angs
        if (unit=='bohr'):
            b2a = 1.0
        coords = self.refList[n][m][2]
        enDist = self.getElecNucDistances(n,m)
        eeDist = self.getElecDistances(n,m)
        for en in enList:
            elist = en[0]; nlist = en[1]
            assert len(nlist)==1 or len(nlist)==2
            print
            print "******** ",en[2]
            print "ee dists:"
            if len(nlist)==1:
                for i in range(len(elist)):
                    ei = elist[i]
                    for j in range(i+1,len(elist)):
                        ej = elist[j]
                        print ' {0:3d} {1:3d} {2:10.4f}'.format(ei, ej, b2a*eeDist[ei-1,ej-1])
                print "en dists:"
                n = nlist[0]
                for e in elist:
                    print ' {0:3d} {1:3d} {2:10.4f}'.format(e,n,b2a*enDist[e-1,n-1])
            elif len(nlist)==2:
                for i in range(len(elist)):
                    ei = elist[i]
                    for j in range(i+1,len(elist)):
                        ej = elist[j]
                        print ' {0:3d} {1:3d} {2:10.4f}'.format(ei, ej, b2a*eeDist[ei-1,ej-1])
                print "en dists:"
                for e in elist:
                    for n in nlist:
                        print ' {0:3d} {1:3d} {2:10.4f}'.format(e,n,b2a*enDist[e-1,n-1])
                n1 = nlist[0]; n2 = nlist[1]
                bdvec = self.nucCoord[n2-1] - self.nucCoord[n1-1]
                bdlen = distance(self.nucCoord[n1-1],self.nucCoord[n2-1])
                print 'bond len = {0:8.4f}'.format(b2a*bdlen)
                # old code with vector to bond
                #print " e dist to bond (with direction) and projection onto bond:"
                #for e in elist:
                #    evec = coords[e-1] - self.nucCoord[n1-1]
                #    # vector from elec to projection of elec on bond vector
                #    bddistvec = evec-np.dot(bdvec0,evec)*bdvec0
                #    bddist = norm2(evec-np.dot(bdvec0,evec)*bdvec0)
                #    bdv = bddistvec / bddist
                #    evecproj = self.nucCoord[n1-1] + np.dot(bdvec0,evec)*bdvec0
                #    print ' {0:3d} {1:10.4f}  ({2:5.2f},{3:5.2f},{4:5.2f}) '.format(e,bddist,bdv[0],bdv[1],bdv[2])
                #    print '       {0:11.5f},{1:11.5f},{2:11.5f} '.format(evecproj[0],evecproj[1],evecproj[2])
                print 'dist to bond proj and dists to nucs:'
                bdvec0 = bdvec/bdlen
                for e in elist:
                    evec = coords[e-1] - self.nucCoord[n1-1]
                    # vector from elec to projection of elec on bond vector
                    bddistvec = evec-np.dot(bdvec0,evec)*bdvec0
                    bddist = norm2(evec-np.dot(bdvec0,evec)*bdvec0)
                    evecproj = self.nucCoord[n1-1] + np.dot(bdvec0,evec)*bdvec0
                    print ' {0:3d} {1:10.4f}'.format(e,b2a*bddist)
                    for n in nlist:
                        print ' {0:3d} {1:3d} {2:10.4f}'.format(e,n,b2a*distance(evecproj,self.nucCoord[n-1]))

    def characterizeRef(self,n=1,m=1,VERBOSE=False):
        nnDist = self.getNucDistances()
        sDist = self.getElecNucSortedDistances(n,m)
        pairdict = {}
    
        for i,d in enumerate(sDist):
            if d[0][1] < self.NUC_THRESH and self.elemList[d[0][0]] != "H":
                elem = self.elemList[d[0][0]]
                f = '1d' if d[0][0] < 9 else '2d'
                s = "core {0:1s}{1:"+f+"}"
                s1 = " {0:3d} core {1:1s}{2:"+f+"}  {3:8.4f}"
                key = s.format(elem,d[0][0]+1)
                if key in pairdict:
                    pairdict[key].append(i+1)
                else:
                    pairdict[key] = [i+1]
                if VERBOSE: print s1.format(i+1,elem,d[0][0]+1,d[0][1])
            elif len(d)>1 and d[0][1]+d[1][1] < nnDist[d[0][0],d[1][0]] + self.BOND_DIFF:
                # sort bond nucs according to geom index
                d1 = d[0] if d[0][0] < d[1][0] else d[1]
                d2 = d[1] if d[0][0] < d[1][0] else d[0]
                elem1 = self.elemList[d1[0]]
                elem2 = self.elemList[d2[0]]
                
                f1 = '1d' if d1[0] < 9 else '2d'
                f2 = '1d' if d2[0] < 9 else '2d'
                s = "bond {0:1s}{1:"+f1+"}-{2:1s}{3:"+f2+"}"
                s1 = " {0:3d} bond {1:1s}{2:"+f1+"}-{3:1s}{4:"+f2+   \
                     "}  {5:8.4f} {6:8.4f} {7:8.4f}"
                key = s.format(elem1,d1[0]+1,elem2,d2[0]+1)
                if key in pairdict:
                    pairdict[key].append(i+1)
                else:
                    pairdict[key] = [i+1]
                if VERBOSE:
                    print s1.format(i+1,elem1,d1[0]+1,elem2,d2[0]+1,
                                   d1[1],d2[1],nnDist[d1[0],d2[0]])
            else:
                elem = self.elemList[d[0][0]]
                f = '1d' if d[0][0] < 9 else '2d'
                s = "lp {0:1s}{1:"+f+"}"
                key = s.format(elem,d[0][0]+1)
                if key in pairdict:
                    pairdict[key].append(i+1)
                else:
                    pairdict[key] = [i+1]
                s1 = " {0:3d} lp {1:1s}{2:"+f+"} {3:8.4f}"
                if VERBOSE: print s1.format(i+1,elem,d[0][0]+1,d[0][1])
        return pairdict

    
    def characterizeRef1(self,n=1,m=1,VERBOSE=False):
        nnDist = self.getNucDistances()
        sDist = self.getElecNucSortedDistances(n,m)
        cgl = []
    
        for i,d in enumerate(sDist):
            if d[0][1] < self.NUC_THRESH and self.elemList[d[0][0]] != "H":
                elem = self.elemList[d[0][0]]
                f = '1d' if d[0][0] < 9 else '2d'
                s = "core {0:1s}{1:"+f+"}"
                s1 = " {0:3d} core {1:1s}{2:"+f+"}  {3:8.4f}"
                key = s.format(elem,d[0][0]+1)
                ltmp = [e[0] for e in cgl]
                if key in ltmp:
                    k = ltmp.index(key)
                    cgl[k][2].append(i+1)
                else:
                    nuclist = [[d[0][0]+1],[]]
                    cgl.append([key,nuclist,[i+1]])
                if VERBOSE: print s1.format(i+1,elem,d[0][0]+1,d[0][1])
            elif len(d)>1 and d[0][1]+d[1][1] < nnDist[d[0][0],d[1][0]] + self.BOND_DIFF:
                # sort bond nucs according to geom index
                d1 = d[0] if d[0][0] < d[1][0] else d[1]
                d2 = d[1] if d[0][0] < d[1][0] else d[0]
                elem1 = self.elemList[d1[0]]
                elem2 = self.elemList[d2[0]]
                
                f1 = '1d' if d1[0] < 9 else '2d'
                f2 = '1d' if d2[0] < 9 else '2d'
                if elem2=='H':
                    stype='bdtH'
                else:
                    stype='bond'
                s = stype+" {0:1s}{1:"+f1+"}-{2:1s}{3:"+f2+"}"
                s1 = " {0:3d} "+stype+" {1:1s}{2:"+f1+"}-{3:1s}{4:"+f2+   \
                     "}  {5:8.4f} {6:8.4f} {7:8.4f}"
                key = s.format(elem1,d1[0]+1,elem2,d2[0]+1)
                ltmp = [e[0] for e in cgl]
                if key in ltmp:
                    k = ltmp.index(key)
                    cgl[k][2].append(i+1)
                else:
                    if elem2=='H':
                        nuclist = [[d2[0]+1],[d1[0]+1]]
                    else:
                        nuclist = [[],[d1[0]+1,d2[0]+1]]
                    cgl.append([key,nuclist,[i+1]])
                if VERBOSE:
                    print s1.format(i+1,elem1,d1[0]+1,elem2,d2[0]+1,
                                   d1[1],d2[1],nnDist[d1[0],d2[0]])
            else:
                elem = self.elemList[d[0][0]]
                f = '1d' if d[0][0] < 9 else '2d'
                s = "lone {0:1s}{1:"+f+"}"
                key = s.format(elem,d[0][0]+1)
                ltmp = [e[0] for e in cgl]
                if key in ltmp:
                    k = ltmp.index(key)
                    cgl[k][2].append(i+1)
                else:
                    nuclist = [[],[d[0][0]+1]]
                    cgl.append([key,nuclist,[i+1]])
                s1 = " {0:3d} lp {1:1s}{2:"+f+"} {3:8.4f}"
                if VERBOSE: print s1.format(i+1,elem,d[0][0]+1,d[0][1])
        return cgl


    def characterizeAndPrintRef(self,n=1,m=1,VERBOSE=False):
        pdict = self.characterizeRef(n,m,VERBOSE)
        plist = dict2plist(pdict)
        printplist(plist)

                
    def getRefSignature(self,n=1,m=1,VERBOSE=False):
        # characterize all atoms
        # bond:sb,mb,lp,core
        # bonds aaaa,aaab etc.: closest spin, d for double elec
        # signature would be number code
        # call characterizeRef and analyze dict for atoms
        nnDist = self.getNucDistances()
        sDist = self.getElecNucSortedDistances(n,m)
        pairdict = {}
    
        for i,d in enumerate(sDist):
            if d[0][1] < self.NUC_THRESH and self.elemList[d[0][0]] != "H":
                elem = self.elemList[d[0][0]]
                f = '1d' if d[0][0] < 9 else '2d'
                s = "core {0:1s}{1:"+f+"}"
                s1 = " {0:3d} core {1:1s}{2:"+f+"}  {3:8.4f}"
                key = s.format(elem,d[0][0]+1)
                if key in pairdict:
                    pairdict[key].append(i+1)
                else:
                    pairdict[key] = [i+1]
                if VERBOSE: print s1.format(i+1,elem,d[0][0]+1,d[0][1])
            #elif d[0][1]+d[1][1] < nnDist[d[0][0],d[1][0]] + self.BOND_DIFF:
            elif d[0][1]**2 + d[1][1]**2 < nnDist[d[0][0],d[1][0]]**2 + self.NUC_THRESH:
                # sort bond nucs according to geom index
                d1 = d[0] if d[0][0] < d[1][0] else d[1]
                d2 = d[1] if d[0][0] < d[1][0] else d[0]
                elem1 = self.elemList[d1[0]]
                elem2 = self.elemList[d2[0]]
                
                f1 = '1d' if d1[0] < 9 else '2d'
                f2 = '1d' if d2[0] < 9 else '2d'
                s = "bond {0:1s}{1:"+f1+"}-{2:1s}{3:"+f2+"}"
                s1 = " {0:3d} bond {1:1s}{2:"+f1+"}-{3:1s}{4:"+f2+   \
                     "}  {5:8.4f} {6:8.4f} {7:8.4f}"
                key = s.format(elem1,d1[0]+1,elem2,d2[0]+1)
                if key in pairdict:
                    pairdict[key].append(i+1)
                else:
                    pairdict[key] = [i+1]
                if VERBOSE:
                    print s1.format(i+1,elem1,d1[0]+1,elem2,d2[0]+1,
                                   d1[1],d2[1],nnDist[d1[0],d2[0]])
            else:
                elem = self.elemList[d[0][0]]
                f = '1d' if d[0][0] < 9 else '2d'
                s = "lp {0:1s}{1:"+f+"}"
                key = s.format(elem,d[0][0]+1)
                if key in pairdict:
                    pairdict[key].append(i+1)
                else:
                    pairdict[key] = [i+1]
                s1 = " {0:3d} lp {1:1s}{2:"+f+"} {3:8.4f}"
                if VERBOSE: print s1.format(i+1,elem,d[0][0]+1,d[0][1])
        return pairdict


        
    def getNucSignature(self,n=1,m=1,VERBOSE=False):
        # return list of counts of electrons at the nucleus
        enDist = self.getElecNucDistances(n,m)
        nlist = []
    
        for j in range(self.nnuc):
            s = 0
            for i in range(self.nElec):
                if enDist[i,j] < self.NUC_THRESH:
                    s += 1
            nlist.append(s)
        return nlist    
