#! /usr/bin/python

import sys,os
import glob
import re
from pylab import *
from pprint import *

dirList = glob.glob('./*.fit')
i=0
for fname in dirList:
  #print fname
  m = re.search('PbSe_(.)',fname)
  size = m.group(1)
  m = re.search('4f - (...)',fname)
  energy = m.group(1)

  #print size + "  " + energy

  f = open(fname, 'r')
  exec(f.read())
  if i==0:
    p = spectrum['peaks'].items()[0][1]
    sys.stdout.write("size, energy, ")
    for var in p['variables']:
      sys.stdout.write(var + ", ")
    sys.stdout.write('\n')
  i += 1

  max_amp = 0
  for key,val in spectrum['peaks'].items():
    if val['values'][0] > max_amp:
      max_amp = val['values'][0]

  for key,val in spectrum['peaks'].items():
    sys.stdout.write("%s, %s, "%(size,energy))
    val['values'][0] = val['values'][0]/max_amp
    for v in val['values']:
      sys.stdout.write(str(v) + ", ")
    sys.stdout.write('\n')
