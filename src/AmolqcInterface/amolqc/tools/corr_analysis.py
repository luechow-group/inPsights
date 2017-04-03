#!/usr/bin/env python
import numpy as np
import numpy.linalg
import sys,getopt
import Walkers

def correlationAnalysis(xyzfile,n1,n2,nRef,charge,mult):
   """ do a correlation analysis for pair n1,n2 of reference nRef
       using electron configuration of file xyzfile
       nalpha is the number of alpha electrons.
       This version analyses only spin pairs, ignoring other configurations 
   """
   assert n1 > 0 and n2 > 0 and nRef > 0," illegal argument"

   w = Walkers.Walkers(xyzfile,charge,mult)
   nnuc = w.getNNuc()
   print xyzfile," xyz file with",nnuc," nuclei"
   print " # electrons      : ",w.getNElec()
   print " # alpha electrons: ",w.getNAlpha()
   nalpha = w.getNAlpha()
   print " correlation analysis for maximum positions ",n1," ",n2
   print " in ref = ",nRef

   covmat1 = np.zeros((3,3))
   covmat2 = np.zeros((3,3))

   mean1 = np.zeros(3)
   mean2 = np.zeros(3)

   r1 = np.zeros(3)
   r2 = np.zeros(3)
   
   theCount1 = 0
   theCount2 = 0
   nopairCount = 0

   # calculate mean and cov matrix for alpha and beta density (of the pair electrons)
   
   n,counter,coords,spin = w.readNextWalkerWithSpin(nalpha)
   while (n != 0):
      if (n==nRef):
         r1 = np.array(coords[n1-1])
         r2 = np.array(coords[n2-1])
         if spin[n1-1] > 0 and spin[n2-1] < 0:
            mean1 += r1
            mean2 += r2
            for i in range(3):
               for j in range(i,3):
                  covmat1[i,j] += r1[i]*r1[j]
                  covmat2[i,j] += r2[i]*r2[j]
            theCount1 += 1
         elif spin[n1-1] < 0 and spin[n2-1] > 0:
            mean1 += r2
            mean2 += r1
            for i in range(3):
               for j in range(i,3):
                  covmat1[i,j] += r2[i]*r2[j]
                  covmat2[i,j] += r1[i]*r1[j]
            theCount2 += 1
         else:
            nopairCount += 1
      n,counter,coords,spin = w.readNextWalkerWithSpin(nalpha)

   theCount = theCount1 + theCount2
   mean1 = mean1 / theCount
   mean2 = mean2 / theCount
   # note the unbiased covariance estimator
   covmat1 = covmat1 / (theCount - 1)
   covmat2 = covmat2 / (theCount - 1)
   for i in range(3):
      for j in range(0,i):
         covmat1[i,j] = covmat1[j,i]
         covmat2[i,j] = covmat2[j,i]
   covmat1 = covmat1 - float(theCount)/float(theCount-1)*np.outer(mean1,mean1)
   covmat2 = covmat2 - float(theCount)/float(theCount-1)*np.outer(mean2,mean2)

   diffmax = np.amin(abs(covmat1-covmat2))
   print " pair, no pair count=",theCount1,theCount2,theCount,nopairCount
   print " mean1 = ",mean1
   print " mean2 = ",mean2
   print " covmat1:"
   for i in range(3):
      print covmat1[i,:]
   print " covmat2:"
   for i in range(3):
      print covmat2[i,:]

   print " min(abs(covmat1-covmat2))=",diffmax

   # diagonalize cov matrix to obtain rotation matrix

   mean = ( mean1 + mean2 ) / 2
   covmat = ( covmat1 + covmat2 ) / 2
   # for i in range(3):
   #    for j in range(0,i):
   #       covmat[i,j] = covmat[j,i]

   print " mean = ",mean
   print " covmat:"
   for i in range(3):
      print covmat[i,:]

   lmda,U = numpy.linalg.eig(covmat)

   print " lmda = ",lmda
   print " U:"
   for i in range(3):
      print U[i,:]

   # elegantly sort eigenvalues and eigenvectors with argsort!
   idx = lmda.argsort()
   lmda = lmda[idx]
   U = U[:,idx]

   print " after sorting: lmda = ",lmda
   print " U:"
   for i in range(3):
      print U[i,:]

   # start again reading xyz-file
   # this time put origin into mean and rotate into main "inertia" axes
   # calculate cov matrices between electrons, separately for x,y,z coordinates

   w = Walkers.Walkers(xyzfile)

   theCount = 0
   nopairCount = 0
   UT = np.transpose(U)
   mean1 = np.zeros(3)
   mean2 = np.zeros(3)
   covmatx = np.zeros((2,2))
   covmaty = np.zeros((2,2))
   covmatz = np.zeros((2,2))
   
   n,counter,coords,spin = w.readNextWalkerWithSpin(nalpha)
   while (n != 0):
      if (n==nRef):
         r1 = np.array(coords[n1-1])
         r1 = np.dot(UT,(r1 - mean))
         r2 = np.array(coords[n2-1])
         r2 = np.dot(UT,(r2 - mean))
         if spin[n1-1] > 0 and spin[n2-1] < 0:
            mean1 += r1
            mean2 += r2
            covmatx[0,0] += r1[0]*r1[0]
            covmatx[0,1] += r1[0]*r2[0]
            covmatx[1,1] += r2[0]*r2[0]
            covmaty[0,0] += r1[1]*r1[1]
            covmaty[0,1] += r1[1]*r2[1]
            covmaty[1,1] += r2[1]*r2[1]
            covmatz[0,0] += r1[2]*r1[2]
            covmatz[0,1] += r1[2]*r2[2]
            covmatz[1,1] += r2[2]*r2[2]
            theCount += 1
         elif spin[n1-1] < 0 and spin[n2-1] > 0:
            mean1 += r2
            mean2 += r1
            covmatx[0,0] += r2[0]*r2[0]
            covmatx[0,1] += r1[0]*r2[0]
            covmatx[1,1] += r1[0]*r1[0]
            covmaty[0,0] += r2[1]*r2[1]
            covmaty[0,1] += r1[1]*r2[1]
            covmaty[1,1] += r1[1]*r1[1]
            covmatz[0,0] += r2[2]*r2[2]
            covmatz[0,1] += r1[2]*r2[2]
            covmatz[1,1] += r1[2]*r1[2]
            theCount += 1
         else:
            nopairCount += 1
      n,counter,coords,spin = w.readNextWalkerWithSpin(nalpha)

   mean1 = mean1 / theCount
   mean2 = mean2 / theCount

   covmatx[1,0] = covmatx[0,1]
   covmaty[1,0] = covmaty[0,1]
   covmatz[1,0] = covmatz[0,1]

   # note the unbiased covariance estimator
   covmatx = covmatx / (theCount - 1)
   covmaty = covmaty / (theCount - 1)
   covmatz = covmatz / (theCount - 1)
   covmatx = covmatx - float(theCount)/float(theCount-1)*(mean1[0]*mean2[0])
   covmaty = covmaty - float(theCount)/float(theCount-1)*(mean1[1]*mean2[1])
   covmatz = covmatz - float(theCount)/float(theCount-1)*(mean1[2]*mean2[2])
   rhox = covmatx[0,1] / np.sqrt(covmatx[0,0]*covmatx[1,1])
   rhoy = covmaty[0,1] / np.sqrt(covmaty[0,0]*covmaty[1,1])
   rhoz = covmatz[0,1] / np.sqrt(covmatz[0,0]*covmatz[1,1])

   print " pair, no pair count = ",theCount,nopairCount
   print " mean1 = ",mean1
   print " mean2 = ",mean2
   print " covmatx:"
   for i in range(2):
      print covmatx[i,:]
   print " covmaty:"
   for i in range(2):
      print covmaty[i,:]
   print " covmatz:"
   for i in range(2):
      print covmatz[i,:]
   print " rho (x,y,z) = ",rhox,rhoy,rhoz 


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
   NOTE: electron and reference counting start with 1!
   '''
   sys.exit(0)


n1 = 0
n2 = 0
nRef = 1
args = sys.argv[1:]
while (len(args)>1):
   if args[0] == '-ep':
      n1 = int(args[1]) 
      n2 = int(args[2])
      args = args[3:]
   if args[0] == '-n':
      nRef = int(args[1])
      args = args[2:]

xyzfile = args[0]
charge = 0
mult = 1

correlationAnalysis(xyzfile,n1,n2,nRef,charge,mult)
