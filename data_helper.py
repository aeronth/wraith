import os
import re
from pylab import *
from scipy.optimize import leastsq
from spectra_fitting import *
from pprint import *
from data_formats import *


def load_BL62_text_files(file_name_filter):
  """helper function to load a bunch of spectra from SSRL BL6-2 text files in the current dir"""
  files = os.listdir('.')

  spectra = {}
  
  for file in files:
    m = re.search('('+file_name_filter+'.*)\.dat',file)
    if (m):
      filename = file
      name = filename
      spectra[name] = Spectrum()
      spectra[name].name = name
      data = np.genfromtxt(filename,skip_header=37)
      #print data
      spectra[name].EE = data[:,1]
      spectra[name].data = data[:,1:]
    
  return spectra

def load_SUPER_text_files(file_name_filter):
  """helper function to load a bunch of spectra from AugerScan text files in the current dir"""
  files = os.listdir('.')

  spectra = {}
  
  for file in files:
    m = re.search('('+file_name_filter+'.*)\.\d\d\d',file)
    if (m):
      filename = file
      name = filename
      spectra[name] = Spectrum()
      spectra[name].name = name
      data = np.genfromtxt(filename,skip_header=12)
      #print data
      spectra[name].EE = data[:,0]
      spectra[name].data = data[:,1:]
    
  return spectra

def load_BL7_XES_files(file_name_filter):
  """helper function to load a bunch of spectra from ALS BL 7.0.1.1 XES text files in the current dir"""
  files = os.listdir('.')
  
  spectra = {}
  
  for file in files:
    m = re.search('('+file_name_filter+'.*)\_spec.txt',file)
    if (m):
      filename = file
      name = filename
      spectra[name] = Spectrum()
      spectra[name].name = name
      data = np.genfromtxt(filename,skip_header=0)
      #print data
      spectra[name].EE = data[:-1,0]
      spectra[name].data = data[:-1,1]
      
  return spectra

def load_BL7_XAS_files(file_name_filter):
  """helper function to load a bunch of spectra from ALS BL 7.0.1.1 XAS text files in the current dir"""
  files = os.listdir('.')
  
  spectra = {}
  
  for f in files:
    m = re.search('('+file_name_filter+'.*)\d\.txt$',f)
    if (m):
      filename = f
      print filename
      name = filename
      nametey = name+'TEY'
      nametfy = name+'TFY'
      spectra[nametey] = Spectrum()
      spectra[nametey].name = nametey

      spectra[nametfy] = Spectrum()
      spectra[nametfy].name = nametfy

      data = read_bl7011_xas(filename)

      spectra[nametey].EE = data['MonoEnergy']
      spectra[nametey].data = data['Counter1']/data['Izero']

      spectra[nametfy].EE = data['MonoEnergy']
      spectra[nametfy].data = data['Counter2']/data['Izero']

  return spectra

def load_AugerScan_text_files(file_name_filter):
  """helper function to load a bunch of spectra from AugerScan text files in the current dir"""
  files = os.listdir('.')

  special_chars = r'[ \t\(\"\(\)\*\&\^\%\$\#\@\!\_\+\-\=\[\]\{\}\|\;\'\"\:\/\.\,\<\>\\\?]'
  
  spectra = {}
  EE=[]
  data=[]
  i = 0;
  for file in files:
    m = re.search('('+file_name_filter+')\.txt',file)
    if (m):
      filename = file
      print filename
      f = open(filename, 'r')
      name = m.group(1)
      header = True
      for line in f:
        if (re.search('Element',line)):
          print '+e'
          header = False
          m = re.search('Element\s(.*);.*;\sDepth Cycle (\d+).*;\sTime Per Step (\d+);\sSweeps (\d+);.*',line)
          #print m.group(0)
          label = i#name+'-'+line
          i = i+1
          #print '+' + label
          spectrum = Spectrum()
          spectrum.element = m.group(1)
          spectrum.depth_cycle = double(m.group(2))
          spectrum.time_per_step = double(m.group(3))
          spectrum.sweeps = double(m.group(4))
          spectrum.name = name
          spectrum.EE=array([])
          spectrum.data=array([])
          spectra[label] = spectrum
        elif header:
          continue
        elif (re.search('^[\d\.]+\s\d+\s*$',line)):
          m = re.search('^([\d\.]+)\s(\d+)',line)
          #print m.group(0)
          spectra[label].EE = append(spectra[label].EE, double(m.group(1)))
          spectra[label].data = append(spectra[label].data, double(m.group(2)))

  return spectra
