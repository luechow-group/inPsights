#!/usr/bin/env python
from numpy import *
import sys,getopt
import References
#
#
if len(sys.argv) < 2:
    print 'arguments missing. use --help for help'
    sys.exit(0)

if '-h' in sys.argv or '--help' in sys.argv:
    print '''
    read and analyse a .ref file
    usage:
    analyseref.py [-n refn] [-m refm] [-cg] [-cgdist] [-d d1 d2] [-a a1 a2 a3] reffile
    -n refn: use ref # 'refn' (default 1)
    -cg determine the "chemical groups" as required for energy partitioning
    -cgdist show distances to nearest nuclei and to linear bond path
    -d d1 d2 calculate distance between d1 and d2. Format "ex" or "nx" where
       x is the number of the electron (ex) or the nucleus (nx)
    -a a1 a2 a3 calculate angle a1-a2-a3. Same format as -d
    '''
    sys.exit(0)



refn = 1
refm = 1
showCG = False
showCGDist = False
showRefSig = False
showBondSig = False
showDist = False
showAngle = False
nearestElecs = False
args = sys.argv[1:]
while (len(args)>1):
    if args[0] == '-n':
        refn = int(args[1])
        args = args[2:]
    if args[0] == '-m':
        refm = int(args[1])
        args = args[2:]
    if args[0] == '-ne':
        nearestElecs = True
        args = args[1:]
    if args[0] == '-cg':
        showCG = True
        args = args[1:]
    if args[0] == '-cgdist':
        showCGDist = True
        args = args[1:]
    if args[0] == '-refsig':
        showRefSig = True
        args = args[1:]
    if args[0] == '-bondsigs':
        showBondSig = True
        args = args[1:]
    if args[0] == '-d':
        showDist = True
        d1 = args[1]
        d2 = args[2]
        args = args[3:]
    if args[0] == '-a':
        showAngle = True
        a1 = args[1]
        a2 = args[2]
        a3 = args[3]
        args = args[4:]
        
fname = args[0]

r = References.References(fname)
ne = r.getNElec()
nnuc = r.getNNuc()
print fname,"ref",refn,refm," with",ne," electrons",nnuc," nuclei"

if showCG or showCGDist:
    cgl = r.characterizeRef1(refn,refm,VERBOSE=False)

if showCG:
    References.printcgentry(cgl)

if showCGDist:
    enlist = References.cglist2enlist(cgl)
    r.writeNucAndBondDistances(enlist,refn,refm)

if showRefSig:
    if r.DoesExist(refn,refm):
        nl = r.getRefSignature(refn,refm,VERBOSE=False)
        ncore,nBonds,nlp = References.analyseRefSignature(nl,VERBOSE=True)
    else:
        print "reference ",refn," ",refm," does not exist"

if showBondSig:
    nmax = r.getNMax()
    mmax = r.getMMax()
    for n in range(1,nmax+1):
        for m in range(1,mmax+1):
            if r.DoesExist(n,m):
                nl = r.getRefSignature(n,m,VERBOSE=False)
                ncore,nBonds,nlp = References.analyseRefSignature(nl,VERBOSE=False)
                if nBonds[4] >= 6:
                    print n,m,nBonds[2],nBonds[4]

if showDist:
    d = r.getDistance(d1,d2,refn,refm)
    print "distance "+d1+"-"+d2+" :",d*r.getBohr2angs()

if showAngle:
    a,ad1,ad2 = r.getAngle(a1,a2,a3,refn,refm)
    print "angle ("+a1+","+a2+","+a3+") :",a,ad1,ad2

if nearestElecs:
    l = r.getNearestElecSortedDistances(refn,refm)
    for i in range(len(l)):
        print "nuc: ",i
        print l[i]





