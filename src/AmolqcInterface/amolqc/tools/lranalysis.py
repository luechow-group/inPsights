#!/usr/bin/env python
from numpy import *
import sys,getopt
import Walkers

def SimpleLeftRightAnalysis():
   assert n1 != 0 and n2 != 0,"-ep argument missing"

   w = Walkers.Walkers(fname)
   nnuc = w.getNNuc()
   print fname,"xyz file with",nnuc," nuclei"

   n,counter,coords = w.readNextWalker()
   theCount = 0
   llCount = 0
   lrCount = 0
   rlCount = 0
   rrCount = 0
   z1sum = 0
   z12sum = 0
   z2sum = 0
   z22sum = 0
   z1z2sum = 0
   while (n != 0):
      if (n==nRef):
         z1sum = z1sum + coords[n1-1][2]
         z2sum = z2sum + coords[n2-1][2]
         z12sum = z12sum + coords[n1-1][2]**2
         z22sum = z22sum + coords[n2-1][2]**2
         z1z2sum = z1z2sum + coords[n1-1][2]*coords[n2-1][2]
         theCount = theCount+1
         if coords[n1-1][2] > 0:
            if coords[n2-1][2] > 0:
               rrCount = rrCount + 1
            else:
               rlCount = rlCount + 1
         else:
            if coords[n2-1][2] > 0:
               lrCount = lrCount + 1
            else:
               llCount = llCount + 1
      n,counter,coords = w.readNextWalker()

   print theCount,llCount,lrCount,rlCount,rrCount
   print float(llCount)/theCount,float(lrCount)/theCount,float(rlCount)/theCount,float(rrCount)/theCount
   print "ll+rr: ",float(llCount+rrCount)/theCount," lr+rl: ",float(rlCount+lrCount)/theCount
   z1 = z1sum/theCount
   z2 = z2sum/theCount
   varz1 = z12sum/theCount
   varz2 = z22sum/theCount
   z1z2 = z1z2sum/theCount
   cov = z1z2-z1*z2
   rho = cov / sqrt(varz1*varz2)
   print "<z1>=",z1," sigma(z1)=",sqrt(varz1), "<z2>=",z2," sigma(varz2)=",sqrt(varz2)
   print "cov(z1,z2)=",cov," rho(z1,z2)=",rho

def LeftMiddleRightAnalysis(halfbondlen,ratio):
   assert n1 != 0 and n2 != 0,"-ep argument missing"
   assert halfbondlen > 0,"-hbl argument missing"

   w = Walkers.Walkers(fname)
   nnuc = w.getNNuc()
   print fname,"xyz file with",nnuc," nuclei"

   n,counter,coords = w.readNextWalker()
   theCount = 0
   lCount = 0
   mCount = 0
   rCount = 0
   while (n != 0):
      if (n==nRef):
         theCount = theCount+1
         isum = 0
         if coords[n1-1][2] > ratio*halfbondlen:
            isum = isum + 1
         elif coords[n1-1][2] < -ratio*halfbondlen:
            isum = isum - 1
         if coords[n2-1][2] > ratio*halfbondlen:
            isum = isum + 1
         elif coords[n2-1][2] < -ratio*halfbondlen:
            isum = isum - 1
         if isum == 2:
            rCount = rCount + 1
         elif isum == -2:
            lCount = lCount + 1
         else:
            mCount = mCount + 1
      n,counter,coords = w.readNextWalker()

   print "left-middle-right analysis (bond assumed symmetric on z axis)"
   print "ratio (middle/bond length)= {0:8.3f}".format(ratio)
   print "bond length (bohr)= {0:10.4f}".format(2*halfbondlen)
   print "total count:",theCount,"left:",lCount,"middle:",mCount,"right:",rCount
   print "ratios left/middle/right: {0:8.3f}  {1:8.3f}  {2:8.3f}".format(float(lCount)/theCount,float(mCount)/theCount,float(rCount)/theCount)
   print "ratios:  {0:8.3f} ionic and {1:8.3f} covalent".format(float(lCount+rCount)/theCount,float(mCount)/theCount)


if len(sys.argv) < 2:
   print 'arguments missing. use --help for help'
   sys.exit(0)

if '-h' in sys.argv or '--help' in sys.argv:
   print '''
   read and analyse a .xyz file
   usage:
   lranalysis.py [-n nRef] [-ep n1 n2] xyzfile
   -ep n1 n2: left right analysis for electron pair n1 n2
   -n nRef: left right analysis for reference nRef (default: nRef=1)
   -hbl ddd: half bond length (distance from origin along z axis in angstrom(!))
   -ratio ddd: ratio middle region/bond length, default 1/3
   NOTE: electron and reference counting start with 1!
   '''
   sys.exit(0)


n1 = 0
n2 = 0
nRef = 1
mode = ''
ratio = 1.0/3.0
halfbondlen = 0.0
BOHR2ANGS = 0.529177249
args = sys.argv[1:]
while (len(args)>1):
   if args[0] == '-ep':
      n1 = int(args[1]) 
      n2 = int(args[2])
      args = args[3:]
   if args[0] == '-n':
      nRef = int(args[1])
      args = args[2:]
   if args[0] == '-lr':
      mode = 'lr'
      args = args[1:]
   if args[0] == '-lmr':
      mode ='lmr'
      args = args[1:]
   if args[0] == '-hbl':
      halfbondlen = float(args[1])
      args = args[2:]
   if args[0] == '-ratio':
      ratio = float(args[1])
      args = args[2:]

fname = args[0]


if (mode=='lr'):
   SimpleLeftRightAnalysis()
elif (mode=='lmr'):
   LeftMiddleRightAnalysis(halfbondlen,ratio)
else:
   print 'mode argument missing'

