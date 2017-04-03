#!/usr/bin/python
# -*- coding: utf-8 -*-

import cookielib, getopt, urllib2, sys
from xml.dom.ext.reader import HtmlLib

HOST = 'http://cccbdb.nist.gov'

opts, (_, datatype, formula) = getopt.getopt(sys.argv, 'b:')

jar = cookielib.CookieJar()
processor = urllib2.HTTPCookieProcessor()
opener = urllib2.build_opener(processor)

# Basissatz auswählen
# 6-31G*: 1, 3-21G: 2, 6-31G**: 3, 6-31+G**: 4, cc-pVDZ: 5, cc-pVTZ: 6, 6-311G*: 7, 3-21G*: 8, 6-31G: 9, 6-311+G(3df,2p): 10, CEP-31G: 11, CEP-31G*: 12, CEP-121G: 13, CEP-121G*: 14, LANL2DZ: 15, SDD: 16, aug-cc-pVDZ: 17, aug-cc-pVTZ: 18, 6-31G(2df,p): 19, STO-3G: 20, 6-311G**: 21, 6-311+G(3df,2pd): 22, cc-pCVDZ: 23, cc-pCVTZ: 24, cc-pVQZ: 25, aug-cc-pVQZ: 26, cc-pV(T+d)Z: 27, aug-cc-pV(T+d)Z: 28
for o, a in opts :
  if o == '-b':
    print "Selecting basis set"
    request = urllib2.Request(HOST + '/selectbasis2.asp',
                              'bchoice=%s&submit1=Submit' % a)
    jar.add_cookie_header(request)
    response = opener.open(request)
    jar.extract_cookies(response, request)

# Landingpage (bspw. geom1.asp, energy1.asp etc.) aufrufen,
# um das Cookie setzen zu lassen.
request = urllib2.Request(HOST + '/%s.asp' % datatype)
response = opener.open(request)
jar.extract_cookies(response, request)

# Cookie sowie Molekül übergeben
request = urllib2.Request(HOST + '/getform.asp',
                         'formula=%s&submit1=Submit' % formula)
jar.add_cookie_header(request)
response = opener.open(request)
jar.extract_cookies(response, request)

reader = HtmlLib.Reader()
doc = reader.fromString(response.read())

if datatype == "expgeom1":
  geomtable = doc.getElementsByTagName("table")[4]
  rows = geomtable.getElementsByTagName("tr")
  for row in rows:
    tds = row.getElementsByTagName("td")
    if not tds: continue

    print '\t'.join(map(lambda el: el.firstChild.nodeValue.strip(), tds))

elif datatype == "energy1":
  # TODO
  print doc.getElementsByTagName("a")
