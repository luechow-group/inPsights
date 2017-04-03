from pylab import *
import sys

if (len(sys.argv) < 2):
    print 'usage: python ',sys.argv[0],' histfilename'
    sys.exit(0)

#mean = float(sys.argv[1])
#sigma = float(sys.argv[2])
fname = sys.argv[1]

x,y1,y2 = loadtxt(fname,usecols=([0,1,2]),unpack=True)

width = x[1] - x[0]  # assuming equidistant points

p = bar(x,y1,width,color='g',align='center')
ylabel('N')
xlabel('E')
title('energy distribution')

plot(x,y2,'b-',label='gauss')

show()
