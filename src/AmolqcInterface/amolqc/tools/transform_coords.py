#!/usr/bin/env python
#
#
import numpy as np
import sys,getopt

def readGeomFromGaussian(comFileName,coordtype='cart'):
    '''read geometry from Gaussian input file. Return geometry
       and lines before and after coords.'''
    try:
        comFile = open(comFileName,'r')
    except IOError as e:
        print 'Could not open com file',comFileName
        print e
        raise

    iniLines = []
    endLines = []
    while True:
        line = comFile.readline()
        if not line: raise EOFError('no coordinates found')
        if '%' not in line:
            break
        iniLines.append(line)
    while True:
        if '#' not in line:
            break
        iniLines.append(line)
        line = comFile.readline()
    iniLines.append(line)
    for i in range(3):
       line = comFile.readline()
       iniLines.append(line)

    coords = []
    if coordtype=='cart':
        while True:
            line = comFile.readline()
            words = line.split()
            if len(words)==0:
                break
            coords.append([words[0],
                           np.array([float(words[1]),float(words[2]),float(words[3])])])
    else:
        print "unknown coordtype ",coordtype
        sys.exit(1)

    endLines.append(line)
    while True:
        line = comFile.readline()
        if not line:
            break
        endLines.append(line)

    comFile.close()

    return (coords,iniLines,endLines)

def writeGaussianComFile(iniLines,coords,endLines,comFileName):
    try:
        comFile = open(comFileName,'w')
    except IOError as e:
        print 'Could not open com file',comFileName
        print e
        raise
    for line in iniLines:
        comFile.write(line)
    for coord in coords:
        s = coord[0]+"  {0:11.6f} {1:11.6f} {2:11.6f}".format(coord[1][0],coord[1][1],coord[1][2])
        comFile.write(s+"\n")
    for line in endLines:
        comFile.write(line)
    comFile.close()


def translateCoords(coords,vec):
    for i in range(len(coords)):
        coords[i][1] = coords[i][1] + vec


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print '''    transform_coords.py 
    usage:
    for Gaussian com files:
    transform_coords.py infile outfile [-t x y z]
    '''
        sys.exit(1)

    infile = sys.argv[1]
    outfile = sys.argv[2]
    args = sys.argv[3:]
    typ = 'translate'
    while len(args)>0:
        if args[0]=='-t':
            typ = 'translate'
            x = float(args[1])
            y = float(args[2])
            z = float(args[3])
            vec = np.array([x,y,z])
            args = args[4:]
        else:
            print 'unrecognized option %s' % args[0]
            sys.exit(1)

    if typ=='translate':
        (coords,start,end) = readGeomFromGaussian(infile)
        translateCoords(coords,vec)
        writeGaussianComFile(start,coords,end,outfile)
    
