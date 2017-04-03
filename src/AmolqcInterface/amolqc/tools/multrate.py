#!/usr/bin/env python
import References
import sys,math

if  len(sys.argv) < 3:
    print '''
    reads a .ref file of a maximum calculation and calculates the multiplicity rates
    usage:    
    multrate.py [reffile.ref] [core_index] [number_of_electrons_to_be_observed]
    
    outputs:
    1. MULTIPLICITY RATE:
    calculates the multiplicity rate, observing the nearest electron maximum positions
    round a core (number given in input).
    2. MULTIPLICITY VALUE:
    calculates the maximum of PSI of the multiplicities by averaging over the maximum 
    values of each value of m_s. The values are divided by the maximum of all values 
    calculated this way.
    3. MULTIPLICITY AVERAGED CORRELATION:
    the averaged correlation of a multiplicity is independent from the given data. 
    Its just printed out as a reference to the overall averaged correlation.
    4. OVERALL AVERAGED CORRELATION:
    the overall averaged correlation is calculated only with regard to maxmultiplicity-1 
    electron maximum positions, disregarding those, that are always 'paired' (most common 
    inside the nucleus). It can be compared to the ep_analysis, as it is average of the 
    average of all correlation numbers between these observed positions.
    '''
    sys.exit(0)
#print sys.argv[1:]
r = References.References(sys.argv[1])
ne = r.getNElec()
nalpha = r.getNAlpha()
nelec = int(sys.argv[3])#default should be equal to the charge of the core
d1='n'+sys.argv[2]
mult_summary = [0]*(nelec+3)
multmax=0
maxtotal=0
mmaxvalue=0
maxvalue = [0]*(nelec+1)
#Rate Calculation
for n in range(1,r.getNMax()+1):
    for m in range(1,r.getMMax()+1):
        if r.DoesExist(n,m):
            dist=[-1.0]*ne
            for i in range(1,ne+1):
                d2='e'+str(i)
                dist[i-1] = r.getDistance(d1,d2,n,m)
            corr=0
            ind = sorted(range(len(dist)), key=lambda x:dist[x])
            for i in range(nelec):
                if ind[i]<nalpha:
                    corr = corr+1
                else:
                    corr = corr-1
            mult=abs(corr)+1
            if mult>multmax:
                multmax=mult
            #print n,m,mult
            mult_summary[mult]+= r.getValueAndFreq(n,m)[1]
            maxtotal+=r.getValueAndFreq(n,m)[1]
            if math.sqrt(math.exp(-r.getValueAndFreq(n,m)[0]))>maxvalue[mult]:
                maxvalue[mult]=math.sqrt(math.exp(-r.getValueAndFreq(n,m)[0]))
rate_summary=[0.0]*(multmax+1)
for i in range(multmax,1,-2):
    rate_summary[i] = float(i)*(mult_summary[i]-mult_summary[i+2])/(2*maxtotal)
if multmax%2!=0:
    rate_summary[1]=1.0*(mult_summary[1]-(mult_summary[3]/2))/maxtotal
#Correlation Calculation
corr_avg=[0.0]*(multmax+1)
overall_corr=0
if multmax>2:
    edges=0
    for i in range(multmax-1):
        edges+=i
    corrs=[0.0]*(multmax+1)
    for i in range(multmax,0,-2):
        alpha=(i-1+multmax-1)/2
        beta=(multmax-1)-alpha
        samea=0
        sameb=0
        for j in range(alpha):
            samea+=j
        for j in range(beta):
            sameb+=j
        corrs[i]=1.0*(edges-2*(samea+sameb))/edges
    for i in range(multmax,0,-2):
        if i==1:
            corr_avg[i]=corrs[i]
        else:
            for j in range(i,1,-2):
                corr_avg[i]+=2*corrs[j]
            corr_avg[i]+=corrs[1]
            corr_avg[i]=corr_avg[i]/i
    for i in range(multmax,0,-2):
        overall_corr+=mult_summary[i]*corrs[i]
    overall_corr=overall_corr/maxtotal
#Value Calculation
value_summary=[0.0]*(multmax+1)
for i in range(multmax,0,-2):
    if i==1:
        value_summary[i]=maxvalue[i]
    else:
        for j in range(i,1,-2):
           # print j            
            value_summary[i]+=2*maxvalue[j]
        value_summary[i]+=maxvalue[1]
        value_summary[i]=value_summary[i]/i
    if value_summary[i]>mmaxvalue:
        mmaxvalue=value_summary[i]    
for i in range(multmax,0,-2):
    #print value_summary[i]
    value_summary[i]=value_summary[i]/mmaxvalue    
#print rate_summary
print 'Pfad:  ',sys.path[0]
print '  Mult    Rate     Value  Correlation'
print '-------------------------------------'
for i in range(multmax,0,-2):
    if corr_avg[i]<0:
        print '{0:4}     {1:1.5f}  {2:1.5f}  {3:1.5f}'.format(i,rate_summary[i],value_summary[i],corr_avg[i])    
    else:
        print '{0:4}     {1:1.5f}  {2:1.5f}   {3:1.5f}'.format(i,rate_summary[i],value_summary[i],corr_avg[i])
print '-------------------------------------'
print 'Overall averaged correlation:',overall_corr
print ''
