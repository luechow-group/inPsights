#!/usr/bin/env python
#
#

import os
import shutil
import sys
import subprocess
#


class AmolqcMasterError(Exception):
   def __init__(self,value="amolqc master script error"):
      self.value = value
   def __str__(self):
      return repr(self.value)

class Amolqc:
   """ class for managing amolqc calculations """
   def __init__(self,mol,orbs,cfg='amolqc.cfg',verbose=True):
      self.mol = mol
      self.cwd = os.getcwd()
      self.orbs = orbs
      self.verbose = verbose
      self.name = self.mol+"-"+self.orbs
      self.jasdir = "Jastrow"
      self.gaudir = "Gaussian"
      self.maxdir = "Max"
      self.seddir = "SED/"+self.orbs
      self.job = False
      self.jobLines = []
      ifile = open(cfg,'r')
      self.runcmd = ifile.readline().rstrip('\n')
      self.jobTemplate = ifile.readlines()
      self.jobFilename = 'jobp'

   def __str__(self):
      s = "amolqc master object:\n"
      s += "name="+self.name+"\n"
      s += self.cwd+"\n"
      s += self.jasdir+"\n"
      s += self.gaudir+"\n"
      s += self.maxdir+"\n"
      s += self.seddir+"\n"
      s += str(self.job)+"\n"
      s += self.runcmd+"\n"
      s += self.jobFilename+"\n"
      return s

   def __repr__(self):
      s = "amolqc master object:\n"
      s += "name="+self.name+"\n"
      s += self.cwd+"\n"
      s += self.jasdir+"\n"
      s += self.gaudir+"\n"
      s += self.maxdir+"\n"
      s += self.seddir+"\n"
      s += str(self.job)+"\n"
      s += self.runcmd+"\n"
      s += self.jobFilename+"\n"
      return s

   def make_jastrow_dir(self):
      if self.verbose: print "creating dir ",self.jasdir
      try:
         os.makedirs(self.jasdir)
      except OSError as e:
         print e

   def make_max_dir(self):
      if self.verbose: print "creating dir ",self.maxdir
      try:
         os.makedirs(self.maxdir)
      except OSError as e:
         print e

   def make_sed_dir(self):
      if self.verbose: print "creating dir ",self.seddir
      try:
         os.makedirs(self.seddir)
      except OSError as e:
         print e

   def create_job(self):
      if self.verbose: print "collecting lines for job file"
      self.job = True
      self.jobLines = []

   def write_job(self):
      if self.verbose: print "writing job file"
      assert self.job and len(self.jobLines)>0, " no job data have been collected"
      of = open(self.jobFilename,"w")
      # replace jobname, cputime, queue ...
      for line in self.jobTemplate:
         of.write(line)
      of.write("cd "+self.cwd+"\n")
      for line in self.jobLines:
         of.write(line+"\n")
      of.write("date\n")
      of.write("echo '===   amolqc job done   ==='\n")
      of.close()
      self.job = False

   def copy_wf(self,suffix,source,dest):
      """ copy wf file to Jastrow """
      if suffix == "":
      	 wfname = self.name+".wf"
      else:
      	 wfname = self.name+"-"+suffix+".wf"
      if source=="Jastrow" or source=="j":
         src = self.jasdir+"/"+wfname
      elif source=="Gaussian" or source=="g":
         src = self.gaudir+"/"+wfname
      elif source=="Max" or source=="m":
         src = self.maxdir+"/"+wfname
      elif source=="SED" or source=="s":
         src = self.seddir+"/"+wfname
      else:
         raise AmolqcMasterError("unknown source")
      if dest=="Jastrow" or dest=="j":
         dst = self.jasdir
      elif dest=="Gaussian" or dest=="g":
         dst = self.gaudir
      elif dest=="Max" or dest=="m":
         dst = self.maxdir
      elif dest=="SED" or dest=="s":
         dst = self.seddir
      else:
         raise AmolqcMaster("unknown dest")
      if self.job:
         if self.verbose: print "adding cp ",src," ",dst," to job"
         self.jobLines.append("cp "+src+" "+dst)
      else:
      	 if self.verbose: print "copying ",src," to ",dst
         try:
            shutil.copy(src,dst)
         except IOError as e:
            print e

   def create_vm(self,suffix,jastrow,E_ref,seed=101,size=500,wf=""):
      """ create varmin in-file """
      assert jastrow[0:2]=='sm' or jastrow[0:2]=='ic',"illegal jastrow given"
      if wf=="":
         wf0name = self.name+".wf"        # start wf file
      else:
         wf0name = self.name+"-"+wf+".wf"
      wfname = self.name+"-"+suffix+".wf"     # optimized wf file
      inname = self.jasdir+"/"+self.name+"-"+suffix+".in"
      try:
         of = open(self.jasdir+"/"+wf0name,"r")
         of.close()
      except IOError as e:
         print "please note: start wf does not (yet) exist"
         ##raise
      try:
         inf = open(inname,"w")
      except IOError as e:
         print e
         raise
      inf.write("$gen(seed="+str(seed)+")\n")
      inf.write("$wf(read,file='"+wf0name+"')\n")
      inf.write("$change_jastrow(new_jastrow="+jastrow+")\n")
      inf.write("$generate_sample(size="+str(size)+")\n")
      inf.write("$jastrow_varmin1(E_ref="+str(E_ref)+")\n")
      inf.write("$wf(write,file='"+wfname+"')\n")
      inf.close()
      if self.verbose: print ".in file "+inname+" for Jastrow variance minimization written"
      if self.job:
         self.jobLines.append("cd "+self.jasdir)
	 self.jobLines.append(self.runcmd+" "+self.name+"-"+suffix)
	 self.jobLines.append("cd ..")
         
 
   def create_vm2(self,suffix,jastrow,E_ref,cnt=5,seed=101,size=500,blocks=20,wf="",stddev=0.001,testwf=True):
      """ create varmin in-file in two steps for sm3 or ic885 Jastrow"""
      assert jastrow=='sm0' or jastrow=='sm1' or jastrow=='ic',"illegal jastrow given"
      if wf=="":
         wf0name = self.name+".wf"        # start wf file
      else:
         wf0name = self.name+"-"+wf+".wf"
      wfname = self.name+"-"+suffix+".wf"     # optimized wf file
      inname = self.jasdir+"/"+self.name+"-"+suffix+".in"
      if testwf:
         try:
            of = open(self.jasdir+"/"+wf0name,"r")
            of.close()
         except IOError as e:
            print "please note: start wf does not (yet) exist"
            ##raise
      try:
         inf = open(inname,"w")
      except IOError as e:
         print e
         raise
      inf.write("$gen(seed="+str(seed)+")\n")
      inf.write("$wf(read,file='"+wf0name+"')\n")
      if jastrow=='sm0':
         inf.write("$change_jastrow(new_jastrow=sm1)\n")
      elif jastrow=='ic220':
         inf.write("$change_jastrow(new_jastrow=ic220)\n")
      inf.write("$sample(create,start=density,generate=random,size="+str(size)+")\n")
      inf.write("$sample(remove_outliers)\n")
      inf.write("$qmc(vmc,move=umr,steps=10,blocks=3,tau=0.01,accept_ratio=0.5,persist=9,discard_all)\n")
      inf.write("$qmc(vmc,move=umr,steps=10,blocks="+str(blocks)+",discard_all)\n")
      inf.write("!\n")
      inf.write("$sample(remove_outliers,tol=5.0,no_replace)\n")
      inf.write("$optimize_parameters(jastrow,method=varmin1,E_ref="+str(E_ref)+",max_iter=3)\n")
      inf.write("$sample(remove_outliers,no_replace)\n")
      inf.write("$sample(change_size,new_size=init_size)\n")
      inf.write("$qmc(vmc,move=umr,steps=10,blocks=3,persist=9,discard_all)\n")
      inf.write("$qmc(vmc,move=umr,steps=10,blocks="+str(blocks)+",discard_all)\n")
      inf.write("!\n")
      if jastrow[0:2]=='sm':
         jas2 = 'sm3'
      else:
         jas2 = 'ic885'
      inf.write("$change_jastrow(new_jastrow="+jas2+")\n")
      inf.write("$begin_loop(count="+str(cnt)+")\n")
      inf.write("$sample(remove_outliers,tol=5.0,no_replace)\n")
      inf.write("$optimize_parameters(jastrow,method=varmin1,E_ref="+str(E_ref)+",max_iter=3)\n")
      inf.write("$sample(remove_outliers,no_replace)\n")
      inf.write("$sample(change_size,new_size=init_size)\n")
      inf.write("$qmc(vmc,move=umr,steps=10,blocks=3,persist=9,discard_all)\n")
      inf.write("$qmc(vmc,move=umr,steps=10,blocks="+str(blocks)+",discard_all)\n")
      inf.write("$end_loop()\n")
      inf.write("$wf(write,file='"+wfname+"')\n")
      
      inf.write("!\n")
      inf.write("$sample(change_size,new_size=100)\n")
      inf.write("$qmc(vmc,steps=300,blocks=1000,discard=5,accept_ratio=0.5,stddev="+str(stddev)+")\n")
      inf.close()
      if self.verbose: print ".in file "+inname+" for two-step Jastrow variance minimization written"
      if self.job:
         self.jobLines.append("cd "+self.jasdir)
	 self.jobLines.append(self.runcmd+" "+self.name+"-"+suffix)
	 self.jobLines.append("cd ..")
 
   def create_em(self,suffix,eqb,seed=101,size=500,cnt=5,wf="",testwf=True):
      """ create em in-file starting from wf using eqb equilib blocks """
      assert eqb==20 or eqb==40 or eqb==80,"only 20, 40 or 80 allowed for eqb"
      if wf=="":
         wf0name = self.name+".wf"        # start wf file
      else:
         wf0name = self.name+"-"+wf+".wf"
      wfname = self.name+"-"+suffix+".wf"     # optimized wf file
      inname = self.jasdir+"/"+self.name+"-"+suffix+".in"
      if testwf:
         try:
            of = open(self.jasdir+"/"+wf0name,"r")
            of.close()
         except IOError as e:
            print "please note: start wf does not (yet) exist"
            ##raise
      try:
         inf = open(inname,"w")
      except IOError as e:
         print e
         raise
      inf.write("$gen(seed="+str(seed)+")\n")
      inf.write("$wf(read,file='"+wf0name+"')\n")
      inf.write("$generate_sample(size="+str(size)+")\n")
      if eqb==40:
         inf.write("$jastrow_emin40(count="+str(cnt)+")\n")
      elif eqb==80:
         inf.write("$jastrow_emin80(count="+str(cnt)+")\n")
      elif eqb==20:
         inf.write("$jastrow_emin20(count="+str(cnt)+")\n")
      inf.write("$wf(write,file='"+wfname+"')\n")
      inf.close()
      if self.verbose: print ".in file "+inname+" for Jastrow energy minimization written"
      if self.job:
         self.jobLines.append("cd "+self.jasdir)
	 self.jobLines.append(self.runcmd+" "+self.name+"-"+suffix)
	 self.jobLines.append("cd ..")
 
   def create_vmc(self,suffix,seed=101,size=100,stddev=0.001,wf=""):
      """ create vmc in-file starting from wf """
      if wf=="":
         wf0name = self.name+".wf"        # start wf file
      else:
         wf0name = self.name+"-"+wf+".wf"
      inname = self.jasdir+"/"+self.name+"-"+suffix+".in"
      try:
         of = open(self.jasdir+"/"+wf0name,"r")
         of.close()
      except IOError as e:
         print "please note: start wf does not (yet) exist"
         ##raise
      try:
         inf = open(inname,"w")
      except IOError as e:
         print e
         raise
      inf.write("$gen(seed="+str(seed)+")\n")
      inf.write("$wf(read,file='"+wf0name+"')\n")
      inf.write("$generate_sample(size="+str(size)+")\n")
      inf.write("$qmc(vmc,move=umr,steps=200,blocks=1000,discard=5,accept_ratio=0.5,stddev="+str(stddev)+")\n")
      inf.close()
      if self.verbose: print ".in file "+inname+" for VMC calc written"
      if self.job:
         self.jobLines.append("cd "+self.jasdir)
	 self.jobLines.append(self.runcmd+" "+self.name+"-"+suffix)
	 self.jobLines.append("cd ..")
    
   def create_maxima(self,suffix,seed=101,size=50,blocks=200,itmax=250,tol=0.001,wf=""):
      """ create vmc in-file for maxima search """
      if wf=="":
         wf0name = self.name+".wf"        # start wf file
      else:
         wf0name = self.name+"-"+wf+".wf"
      inname = self.maxdir+"/"+self.name+"-"+suffix+".in"
      try:
         of = open(self.maxdir+"/"+wf0name,"r")
         of.close()
      except IOError as e:
         print "please note: start wf does not exist"
         ##raise
      try:
         inf = open(inname,"w")
      except IOError as e:
         print e
         raise
      inf.write("$gen(seed="+str(seed)+")\n")
      inf.write("$wf(read,file='"+wf0name+"')\n")
      inf.write("$generate_sample(size="+str(size)+")\n")
      inf.write("$qmc(vmc,move=umr,steps=200,blocks="+str(blocks)+",discard=5,accept_ratio=0.5,\n")
      inf.write("     searchmax,mode=bfgs,itmax="+str(itmax)+",nmax=30,tol="+str(tol)+")\n")
      inf.close()
      if self.verbose: print ".in file "+inname+" for maxima calc written"
      if self.job:
         self.jobLines.append("cd "+self.maxdir)
	 self.jobLines.append(self.runcmd+" "+self.name+"-"+suffix)
	 self.jobLines.append("cd ..")
  
   def create_coc(self,suffix,max="",iters=3,seed=101,size=200,blocks=50,tol=0.001,wf=""):
      """ create in-file for coc search """
      if wf=="":
         wf0name = self.name+".wf"        # start wf file
      else:
         wf0name = self.name+"-"+wf+".wf"
      inname = self.maxdir+"/"+self.name+"-"+suffix+".in"
      try:
         of = open(self.maxdir+"/"+wf0name,"r")
         of.close()
      except IOError as e:
         print "please note: start wf does not exist"
         ##raise
      try:
         inf = open(inname,"w")
      except IOError as e:
         print e
         raise
      inf.write("$gen(seed="+str(seed)+")\n")
      inf.write("$wf(read,file='"+wf0name+"')\n")
      inf.write("$generate_sample(size="+str(size)+")\n")
      inf.write("$iterate_coc(vmc,move=umr,steps=200,blocks="+str(blocks)+",discard=3,accept_ratio=0.5,\n")
      if max=="":
         inf.write("     iters="+str(iters)+",readref,ref=1,tol="+str(tol)+")\n")
      else:
         maxfile = self.name+"-"+max+".ref"
         inf.write("     iters="+str(iters)+",ref_file="+maxfile+",ref=1,tol="+str(tol)+")\n")
      inf.close()
      if self.verbose: print ".in file "+inname+" for coc calc written"
      if self.job:
         self.jobLines.append("cd "+self.seddir)
	 self.jobLines.append(self.runcmd+" "+self.name+"-"+suffix)
	 self.jobLines.append("cd ..")
 
   def create_seds(self,suffix,max="",ref=1,seed=101,size=200,blocks=200,stddev=0.0002,grid=4.0,wf=""):
      """ create vmc in-file for maxima search """
      if wf=="":
         wf0name = self.name+".wf"        # start wf file
      else:
         wf0name = self.name+"-"+wf+".wf"
      inname = self.seddir+"/"+self.name+"-"+suffix+".in"
      try:
         of = open(self.seddir+"/"+wf0name,"r")
         of.close()
      except IOError as e:
         print "please note: start wf does not (yet) exist"
         ##raise
      try:
         inf = open(inname,"w")
      except IOError as e:
         print e
         raise
      inf.write("$gen(seed="+str(seed)+")\n")
      inf.write("$wf(read,file='"+wf0name+"',epart)\n")
      inf.write("$generate_sample(size="+str(size)+")\n")
      inf.write("$qmc(vmc,move=umr,steps=200,blocks="+str(blocks)+",discard=5,accept_ratio=0.5,stddev="+str(stddev)+",\n")
      if max=="":
         inf.write("     sed,epart,readref,ref="+str(ref)+",grid="+str(grid)+")\n")
      else:
         maxfile = self.name+"-"+max+".ref"
         inf.write("     sed,epart,ref_file="+maxfile+",ref="+str(ref)+",grid="+str(grid)+")\n")
      inf.close()
      if self.verbose: print ".in file "+inname+" for sed and epart calc written"
      if self.job:
         self.jobLines.append("cd "+self.seddir)
	 self.jobLines.append(self.runcmd+" "+self.name+"-"+suffix)
	 self.jobLines.append("cd ..")

   def isFinished(self,dir,suffix):
      """ return True if calculation with suffix has successfully finished """
      if dir=="Jastrow" or dir=="j":
         outname = self.jasdir
      elif dir=="Gaussian" or dir=="g":
         outname = self.gaudir
      elif dir=="SED" or dir=="s":
         outname = self.seddir
      if suffix=="":
         outname += "/"+self.name+".out"
      else:
         outname += "/"+self.name+"-"+suffix+".out"
      try:
         of = open(outname,"r")
      except IOError as e:
         print "out file does not exist",e
         raise
      lines = of.readlines()
      of.close()
      return lines[-3][0:4]=="Bye!"
 
   def result(self,dir,suffix):
      """ get results from calculations """
      if dir=="Jastrow" or dir=="j":
         outname = self.jasdir
      elif dir=="Gaussian" or dir=="g":
         outname = self.gaudir
      elif dir=="SED" or dir=="s":
         outname = self.seddir
      if suffix=="":
         outname += "/"+self.name+".out"
      else:
         outname += "/"+self.name+"-"+suffix+".out"
      try:
         of = open(outname,"r")
      except IOError as e:
         print "out file does not exist",e
         raise
      lines = of.readlines()
      normalTermination = self.isFinished(dir,suffix)
      word = lines[-5].split()
      date = word[5]+"-"+word[6]+"-"+word[7]+"-"+word[8]
      word = lines[-8].split()
      cputime = word[6]+":"+word[7]
      word = lines[-9].split()
      walltime = word[6]+":"+word[7]
      word = lines[7].split()
      version = word[1]
      word = lines[10].split()
      nproc = int(word[-2])
      energy = []
      stddev = []
      variance = []
      ncorr = []
      for i in range(len(lines)):
         if "FINAL RESULT" in lines[i]:
            word = lines[i+1].split()
            energy.append(float(word[3]))
            stddev.append(float(word[5]))
            word = lines[i+3].split()
            variance.append(float(word[5]))
            word = lines[i+4].split()
            ncorr.append(float(word[2]))
      return (energy,stddev,variance,ncorr,nproc,cputime,walltime,version,date,outname)
            
 
if __name__ == "__main__":
   print "testing Amolqc master module"
   a = Amolqc("n2h4","N2H4","rhf-tzp")
   print a


