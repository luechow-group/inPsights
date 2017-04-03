#!/usr/bin/python2.7
import glob
from optparse import OptionParser
import os
import re
import sys

modules = {}
uses = {}
includes = {}
ignore = ("omp_lib")

module_re = re.compile(r"\s*module\s+(\S+)\s*$", re.IGNORECASE)
use_re = re.compile(r"\s*use\s+(\w+)", re.IGNORECASE)

def findUse(name, line):
  m = use_re.match(line)
  if m:
    uselist = uses.setdefault(name, [])
    mod = m.group(1).lower()
    if mod not in uselist:
      log("Found use directive for %s in %s" % (mod, name))
      uselist.append(mod)
    uses[name] = uselist

def findModule(name, line, mods=modules):
  m = module_re.match(line)
  if m:
    mod = m.group(1).lower()
    log("Found module declaration for %s in %s" % (mod, name))
    if mod in mods:
      print >> sys.stderr, "Duplicate module name %s in %s and %s" % (mod, mods[mod], name)
      sys.exit(1)

    mods[m.group(1).lower()] = name

def findInclude(name, line):
  findModule(name, line, includes)

def findBoth(name, line):
  findModule(name, line)
  findUse(name, line)

def readFile(name, cb):
  with open(name, 'r') as f:
    for line in f:
      cb(name, line)

def globFiles(path, cb):
  for name in glob.glob(os.path.join(path, "*.f")):
    readFile(name, cb)

  for name in glob.glob(os.path.join(path, "*.f90")):
    readFile(name, cb)

def toObject(name):
  return re.sub(r"\.f(?:90)?$", ".o", name)

parser = OptionParser()
parser.add_option("-I" , "--include", action="store", type="string", dest="include")
parser.add_option("-b", "--basename", action="store_true", dest="basename")
parser.add_option("-d", "--dirprefix", action="store_true", dest="dirprefix")
parser.add_option("-v", "--verbose", action="store_true", dest="verbose")
parser.add_option("-f", "--file", action="store_true", dest="file")
(options, args) = parser.parse_args()

def log(msg):
  global options
  if not options.verbose: return
  print >> sys.stderr, msg

if options.include:
  if not os.path.isdir(options.include):
    print >> sys.stderr, "Include path has to be directory"
    sys.exit(1)

  path = os.path.normpath(options.include)
  globFiles(path, findInclude)

if options.basename:
  processName = lambda x: os.path.basename(x)
  processInclude = lambda x: x
elif options.dirprefix:
  processName = lambda x: "$(dir)/%s" % (os.path.basename(x))
  processInclude = lambda x: "$(utilsdir)/%s" % (os.path.basename(x))
else:
  processName = lambda x: x
  processInclude = lambda x: x

if not args:
  print "Usage: python %s [-I include-path] [-b|--basename] [-d|--dirprefix] [-v|--verbose] [-f|--file] path" % sys.argv[0]
  sys.exit(1)

name, = args
if options.file:
  if not os.path.isfile(name):
    print >> sys.stderr, "--file specified but target path is no file"
    sys.exit(1)

  readFile(name, findBoth)
else:
  if not os.path.isdir(name):
    print >> sys.stderr, "Target path has to be directory"
    sys.exit(1)

  globFiles(name, findBoth)

print "# dependency list (created with moduledepends.py)"

for path in sorted(uses):
  name = processName(path)
  deps = []
  missing = []
  for use in uses[path]:
    if use in modules:
      if modules[use] != path:
        deps.append(toObject(processName(modules[use])))
    elif use in includes:
      deps.append(toObject(processInclude(includes[use])))
    elif use in ignore:
      pass
    else:
      missing.append(use)

  deplist = list(set(deps))

  if deplist:
    print toObject(name) + ":", name, ' '.join(deplist)
  if missing:
    missinglist = ' '.join(list(set(missing)))
    print "# Missing in %s: %s" % (path, missinglist)
    print >> sys.stderr, "Missing modules in %s: %s" % (path, missinglist)
