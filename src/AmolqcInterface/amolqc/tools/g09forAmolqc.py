#!/usr/bin/env python
import sys,os
import subprocess
from gaussian2wf import gaussian2wf

pse = ['X','H','He',
       'Li','Be','B','C','N','O','F','Ne',
       'Na','Mg','Al','Si','P','S','Cl','Ar',
       'K','Ca',
           'Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn',
                'Ga','Ge','As','Se','Br','Kr']

pseDict = {'H':'Hydrogen','He':'Helium',
   'Li':'Lithium','Be':'Beryllium','B':'Boron','C':'Carbon','N':'Nitrogen','O':'Oxygen','F':'Fluorine','Ne':'Neon',
   'Na':'Sodium','Mg':'Magnesium','Al':'Aluminum','Si':'Silicon','P':'Phosphorous','S':'Sulfur','Cl':'Chlorine','Ar':'Argon',
   'K':'Potassium','Ca':'Calcium',
   'Sc':'Scandium','Ti':'Titanium','V':'Vanadium','Cr':'Chromium','Mn':'Manganese','Fe':'Iron','Co':'Cobalt','Ni':'Nickel',
   'Cu':'Copper','Zn':'Zinc',
   'Ga':'Gallium','Ge':'Germanium','As':'Arsenic','Se':'Selenium','Br':'Bromium','Kr':'Krypton'}

defaultDict = {}

def parseFile(fname):
    infile = open(fname,"r")
    inDict = defaultDict
    lines = infile.readlines()
    for line in lines:
        if line=='': continue
        words = line.split("=")
        key = words[0].strip()
        value = words[1].strip()
        inDict[key] = value
    return inDict

def readMolFile(fname):
    gfile = open(fname,"r")
    alllines = gfile.read().split("\n")
    geom = []
    for line in alllines:
        words = line.split()
        if (len(words)>=4):
            geom.append([words[0],words[1],words[2],words[3]])
    for e in geom:
        try:
            elemTyp = pse[int(e[0])]
            elemIdx = int(e[0])
            e.insert(0,elemTyp)
        except:
            elemTyp = e[0]
            elemIdx = pse.index(e[0])
            e.insert(1,elemIdx)
    return geom

def uniqueAtoms(geom):
    uAtoms = []
    for e in geom:
        if e[0] not in uAtoms:
            uAtoms.append(e[0])
    return uAtoms

def readOptLogFile(fname,mstdorient):
    try:
        logfile = open(fname,"r")
    except IOError as e:
        print 'Could not open log file ',fname
        print e
        raise
    if mstdorient:
        orient = 'Standard orientation'
    else:
        orient = 'Input orientation'
    while True:
        try:
            line = logfile.readline()
        except IOError as e:
            print 'IO-Error while reading ',fname
            print e
            raise
        if line == '':
            print 'Charge not found in log-file ',fname
            raise
        if 'Charge =' in line:
            words = line.split()
            charge = int(words[2])
            mult = int(words[5])
            break
    while True:
        try:
            line = logfile.readline()
        except IOError as e:
            print 'IO-Error while reading ',fname
            print e
            raise
        if line == '':
            print 'log file appears not to have converged'
            raise
        if '-- Stationary point found' in line:
            break
    while True:
        try:
            line = logfile.readline()
        except IOError as e:
            print 'IO-Error while reading ',fname
            print e
            raise
        if line == '':
            print orient,' not found in log file'
            raise
        if orient in line:
            for i in range(5):
                line = logfile.readline()
            geometry = []
            while '----' not in line:
                words = line.split()
                geometry.append([pse[int(words[1])],int(words[1]),words[3],words[4],words[5]])
                line = logfile.readline()
            break
    logfile.close()
    return charge,mult,geometry

def writeGeom(geom,handle):
    for e in geom:
       handle.write(e[0]+' '+e[2]+' '+e[3]+' '+e[4]+'\n')

def writeSortedGeom(geom,handle):
    geom = sorted(geom,key=lambda e: e[1],reverse=True)
    for e in geom:
       handle.write(e[0]+' '+e[2]+' '+e[3]+' '+e[4]+'\n')

def readBasisFromGbs(basis,uniqueAtoms):
    basispath = os.environ['AMOLQC']
    basisfile = open(basispath+"/bib/"+basis,"r")
    print "basisfile=",basispath+"/bib/"+basis
    gbasis = []
    for atom in uniqueAtoms:
        read = False
        start = '-'+atom
        print "start =",start
        basisfile.seek(0)
        for line in basisfile:
            if start in line:
                gbasis.append(atom+'  0'+"\n")
                read = True
            if '****' in line:
                if read:
                    gbasis.append(line)
                    read = False
                    break
            if read==True:
                if start not in line:
                    gbasis.append(line)
    return gbasis

def readBasisFromBaslib(basis,uniqueAtoms):
    basispath = os.environ['AMOLQC']
    basisfile = open(basispath+"/bib/"+basis,"r")
    print "basisfile=",basispath+"/bib/"+basis
    gbasis = []
    for atom in uniqueAtoms:
        read = False
        start = '####'+pseDict[atom]
        print "start =",start
        basisfile.seek(0)
        for line in basisfile:
            if start in line:
                gbasis.append(atom+' 0'+"\n")
                read = True
                basisfile.next()
            if '****' in line:
                if read:
                    gbasis.append(line)
                    read = False
                    break
            if read==True:
                if start not in line:
                    gbasis.append(line)
    return gbasis

def readEcpBasisFromGbs(basis,ecpAtoms,skip):
    basispath = os.environ['AMOLQC']
    basisfile = open(basispath+"/bib/"+basis,"r")
    print "basisfile=",basispath+"/bib/"+basis
    ebasis = []
    for atom in ecpAtoms:
        read = False
        start = atom+" 0"
        if not(skip): start='-'+atom+' 0'
        print "start =",start
        basisfile.seek(0)
        for line in basisfile:
            if start in line:
                ebasis.append(atom+' 0'+"\n")
                read = True
            if '****' in line:
                if read:
                    if not(skip): ebasis.append(line)
                    read = False
                    break
            if read==True:
                if start not in line:
                    ebasis.append(line)
    return ebasis

def getGBSFile(basisName):
    basispath = os.environ['AMOLQC']
    gbsEntry = basispath+"/bib/"+basisName
    basisfile = open(gbsEntry,"r")
    basisfile.close()
    return gbsEntry


ecptype = 'none'
ecpAtoms = []
withGaussian = True
if len(sys.argv) > 1:
    if sys.argv[1]=='-ng':
        withGaussian = False
        sys.argv.remove('-ng')

if len(sys.argv) > 1:
    fname = sys.argv[1]
    interactive = False
    inDict = parseFile(fname)
else:
    interactive = True

print
print "   * * *  g09forAmolqc  * * *"
print
print "composing and running a single point g09 calculation to generate MOs for amolqc"
print
print "version 0.91"
print
if interactive:
    print "Enter name for calculation (without extension) for .com/.log/.wf files"
    basename = str(raw_input("name for calculation: "))
else:
    basename = inDict['name']
print "base name = ",basename
if interactive:
    print "Which method to use for the orbitals? (Use legal g09 method with prefix 'R','RO' or 'U', e.g. RB3LYP)"
    method = raw_input("method for orbitals: ")
else:
    method = inDict['method']
if not (method[0]=='R' or method[0]=='U'):
    print "method name needs to start with 'R' or 'U'"
    sys.exit(1)
print "method = ",method
if interactive:
    print "Geometry is taken from gaussian log file (only converged opt calculations!) or a 'geom' file."
    geomFile = raw_input("log file name or 'geom' (in gaussian geometry format): ")
else:
    geomFile = inDict['geometry']

useLogFile = False
if '.log' in geomFile:
    useLogFile = True
if useLogFile:
    print "geometry from gaussian log file ",geomFile
else:
    print "geometry from file ",geomFile

stdorient=True
if useLogFile:
   if interactive:
       yn = raw_input("Use Input orientation (i.e. nosymm used for calculation)? (y/n) ")
   else:
      try:
         yn = inDict['nostdorient']
      except:
        yn='n'
   if yn=='n':
       stdorient=True
       print 'Searching for Standard Orientation in log-file.'
   else:
       stdorient=False
       print 'Searching for Input Orientation in log-file.'

if interactive:
    print "Enter Units, which are used within the geometry (Note: fort.7 Gaussian-Files uses Bohr)"
    yn = raw_input("Geometry given in Bohr (y: Geometry in Bohr; default:n): ")
else:
    try:
        yn = inDict['bohr']
    except:
        yn='n'
units=''
if yn=='y':
    units='units=au'
    print 'Atomic Units are used for input geometry'
else:
    units=''
if interactive:
    print "Now enter the name of the basis set or 'read' for reading the basis from the file 'basis'"
    print "Note that this is the basis set for gaussian calculation and the amolqc wavefunction."
    print "Only those basis sets that are available as '.gbs' files in the $AMOLQC/bib directory are allowed."
    print "Common examples are 'cc-pVTZ-f' and 'TZPAE'"
    basisName  = raw_input("Which basis in g09? (e.g.  'read' for reading 'basis'): ")
else:
    basisName = inDict['basis']
print "basis for orbitals: ",basisName

ecp = True
if interactive:
   if raw_input("Do you want to use Effective Core Potentials? (y/n)") == 'y':
      ecp = True
      ecpType = raw_input("Now enter the name of the ECP-type you want to use (for example: BFD-TZP).")
      atomType = raw_input("Now enter the atoms you want to use ECPs for (form: Si,O,C).")
   else:
      ecp = False
else:
   try:
      ecpType = inDict['ecptype']
      atomType = inDict['ecpatoms']
   except:
      ecp=False
   if ecp:
      if '-' in ecpType:
        ecpAtoms = str.split(atomType, ",")
if ecp:
    print "ecp type: ",ecpType
    print "ecp atoms: ",atomType

if useLogFile==False:
    if interactive:
        charge = int(raw_input("Charge? "))
        multiplicity = int(raw_input("Multiplicity? "))
    else:
        charge = int(inDict['charge'])
        multiplicity = int(inDict['multiplicity'])

if interactive:
    yn = raw_input("Input orientation (i.e. nosymm)? (y/n) ")
else:
    yn = inDict['nosymm']
if yn=='y':
    nosymm='NoSymm'
else:
    nosymm=''
print 'nosymm=',nosymm

if interactive:
    moType = raw_input("Do you want to localize MO (Boys localization)? (y/n) ")
else:
    moType = inDict['lmo']

if interactive:
    print "Give a title (in ' ') containing the system and the origin of the geometry (e.g. B3LYP/6-31G(d))"
    title = raw_input("title: ")
else:
    title = inDict['title']
print "title: ",title

nProcs=4
if interactive:
    nProcs = int(raw_input("How may CPUs to use for wf-generation?"))
else:
   try:
      nProcs = int(inDict['nProcs'])
   except:
      nProcs = 4
print 'Using nProcs=', str(nProcs)

compath = basename+".com"
basispath="basis"

if useLogFile:
    charge,multiplicity,geom = readOptLogFile(geomFile,stdorient)
else:
    geom = readMolFile(geomFile)

gbs = False
gbsEntry = ""
if basisName=="read":
    try:
        basisfile = open(basispath,"r")
        gbasis = basisfile.readl()
    except IOError as e:
        print basispath," could not be opened"
        print e
        raise
else:
   if ecp==False:
       try:
           gbs = True
           gbsEntry = getGBSFile(basisName+'.gbs')
       except IOError:
           print basisName+'.gbs', " not found in $AMOLQC/bib "
           gbs = False
           uAtoms = uniqueAtoms(geom)
           print 'unique atoms:',uAtoms
           try:
               gbasis = readBasisFromBaslib(basisName,uAtoms)
           except IOError as e:
               print basisName, "not found in $AMOLQC/bib "
               raise
   else:
      print "Using ECPs"
      gbs = False
      tbasisName = basisName+'.gbs'
      try:
         gbsEntry = getGBSFile(basisName+'.gbs')
      except IOError:
         print tbasisName+'.gbs', " not found in $AMOLQC/bib "
      try:
         ecpEntry = getGBSFile(ecpType+'.gbs')
         ecptype = ecpType
         ecpType = ecpType+'.gbs'
      except IOError:
         print ecpType+'.gbs', " not found in $AMOLQC/bib "
      uAtoms = uniqueAtoms(geom)
      nonecpAtoms = set(uAtoms)^set(ecpAtoms)
      nonecpAtoms = list(nonecpAtoms)
      print 'non-ECP atoms:',nonecpAtoms
      print 'ECP atoms',ecpAtoms
      try:
         gbasis = readBasisFromGbs(tbasisName,nonecpAtoms)
      except IOError as e:
         print tbasisName, "not found in $AMOLQC/bib "
         raise
      try:
         ebasis = readEcpBasisFromGbs(ecpType,ecpAtoms,False)
      except IOError as e:
         print tbasisName, "not found in $AMOLQC/bib "
         raise
      try:
         ecpType = str.split(ecpType,"-")
         ecpDef = readEcpBasisFromGbs(ecpType[0]+'.gbs',ecpAtoms,True)
      except IOError:
         print ecpType[0]+'.gbs', "not found in $AMOLQC/bib "

##################  WRITE GAUSSIAN INPUT  ############################################
gaussian = open(compath,"w");
gaussian.write("%chk="+basename+".chk"+"\n");
gaussian.write("%NProcShared="+str(nProcs)+"\n");
gaussian.write("%mem=5000MB"+"\n");
if ecp==False:
   gaussian.write("#p gfinput %s/gen 6D 10F int=nobasistransform punch=mo scf=tight %s %s\n" % (method, nosymm, units));
else:
   gaussian.write("#p gfinput %s/genecp 6D 10F int=nobasistransform punch=mo scf=tight %s %s\n" % (method, nosymm, units));
gaussian.write("\n");
gaussian.write(title+"\n");
gaussian.write("\n");
gaussian.write(str(charge)+" "+str(multiplicity)+"\n");
writeSortedGeom(geom,gaussian)
gaussian.write("\n");
if ecp == False:
   if gbs:
       gaussian.write("@"+gbsEntry+"\n")
   else:
       for item in gbasis:
           gaussian.write(item);
else:
   for item in gbasis:
      gaussian.write(item);
   for item in ebasis:
      gaussian.write(item);
   gaussian.write("\n");
   for item in ecpDef:
      gaussian.write(item)
gaussian.write("\n");
if moType=='y':
    gaussian.write('--Link1--'+"\n")
    gaussian.write("%chk="+basename+".chk"+"\n");
    gaussian.write("%NProcShared=4\n");
    gaussian.write("#p gen gfinput 6D 10F geom=check guess(save,read,local,only) int=nobasistransform pop=full punch=MO "+nosymm+" "+units+"\n");
    gaussian.write("\n");
    gaussian.write(title+' LOCAL-BOYS'+"\n");
    gaussian.write("\n");
    gaussian.write(str(charge)+" "+str(multiplicity)+"\n");
    gaussian.write("\n");
    if gbs:
        gaussian.write("@"+gbsEntry+"\n")
    else:
        for item in gbasis:
            gaussian.write(item);
    gaussian.write("\n");
gaussian.close();
print "gaussian input file has been written"
if withGaussian:
    print "now running g09 ..."
    cmd = ['g09',compath]
    subprocess.check_call(cmd)
    print "g09 has terminated"
    print "now creating wf file ..."
    if nosymm=='NoSymm':
        gaussian2wf(basename,basisName,method=method,title=title,stdorient=False,ecp=ecptype,ecpatoms=ecpAtoms)
    else:
        gaussian2wf(basename,basisName,method=method,title=title,stdorient=True,ecp=ecptype,ecpatoms=ecpAtoms)
print "done"
