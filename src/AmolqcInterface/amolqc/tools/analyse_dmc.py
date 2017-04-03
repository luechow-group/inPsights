#!/usr/bin/env python
from numpy import *
import sys,getopt
#
if len(sys.argv) < 2:
    print 'arguments missing. use --help for help'
    sys.exit(0)

if '-h' in sys.argv or '--help' in sys.argv:
    print '''
    read and analyse the block output of a dmc run (saved as .dat)
    usage:
    analyse_dmc.py [-bn blockn] [-var var] [-k K] datfile
    -bn blockn: block length (number of steps in block)
    -var var: variance of E_loc as given in output
    -k K=5: calculate running averages using K subsequent data

    calculates the block sigma_i ( sigma_i^(bl) = < E_block**2 > - < E_block >^2 )
    and block mean ( < E_block >) as running average of K subsequent blocks.
    Each block energy is average of size * blockn data. current size is read from output.
    If block energies are independent, the number of independent data in the block energy
    is:
    N_i = var / sigma_i**2 + 1 
    global correlation length is: size*blockn / N_i
    '''
    sys.exit(0)


K = 5
blockn = 1
vari = 0
calcNcorr = False
args = sys.argv[1:]
while (len(args)>1):
    if args[0] == '-k':
        K = int(args[1])
        args = args[2:]
    if args[0] == '-bn':
        blockn = int(args[1])
        args = args[2:]
    if args[0] == '-var':
        vari = float(args[1])
        args = args[2:]
        calcNcorr = True
fname = args[0]

datFile = open(fname,"r")
blockE = []
size = []
while True:
    line = datFile.readline()
    if not line: break
    words = line.split()
    if len(words) < 3: break
    blockE.append(float(words[2]))
    size.append(float(words[1]))

print fname,": ",len(blockE)," block energies read"
if calcNcorr:
    print "using var=",vari," sigma_i=",sqrt(vari)," steps in block: ",blockn
print "running mean and sigma_i [and N_i, N_corr] for ",K," steps:"

for i in range(len(blockE)-K):
    sum1 = 0
    sum2 = 0
    ssum = 0
    for j in range(K):
        sum1 += blockE[i+j]
        sum2 += blockE[i+j]**2
        ssum += size[i+j]
    mean = sum1 / K
    sigma = sqrt(sum2/K - mean**2)
    means = ssum / K
    if calcNcorr:
        print '{0:5d} {1:12.5f} {2:10.5f} {3:15.1f} {4:10.1f}'.format(i,mean,sigma,vari / sigma**2 + 1, means*blockn / (vari / sigma**2 + 1))
    else:
        print '{0:5d} {1:12.5f} {2:10.5f}'.format(i,mean,sigma)

datFile.close()

tsum1 = 0
tsum2 = 0
tssum = 0
for i in range(len(blockE)):
    tsum1 += blockE[i]
    tsum2 += blockE[i]**2
    tssum += size[i]
mean = tsum1 / len(blockE)
sigma = sqrt(tsum2/len(blockE) - mean**2)
means = tssum / len(blockE)
print "mean and sigma_i for all blocks:"
if calcNcorr:
    print '{0:12.5f} {1:10.5f} {2:15.1f} {3:10.1f}'.format(mean,sigma,vari / sigma**2 + 1, means*blockn / (vari / sigma**2 + 1))
else:
    print '{0:12.5f} {1:10.5f}'.format(mean,sigma)
print "mean and sigma_m (if blocks independent):"
print '{0:12.5f} {1:10.5f}'.format(mean,sigma / sqrt(len(blockE)-1))



