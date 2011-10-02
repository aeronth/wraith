import os
import re
from pylab import *
from scipy.optimize import leastsq
from fitting_machinery import *
from pprint import *


class Spectrum:
  """A place to keep a spectrum"""
  
  def __init__(self):
    """Initialize stuff"""
    self.data = ones(1)
    self.EE = ones(1)
    self.offset = 0
    self.bg = lambda E, data: zeros(E.size)
    self.abg = lambda E: zeros(E.size)
    self.peaks = Peaks()
    self.name = "unnamed"

  def clear_peaks(self):
    self.peaks = Peaks()

  def clear_bg(self):
    self.bg = lambda E, data: zeros(E.size)

  def clear_abg(self):
    self.abg = lambda E: zeros(E.size)

  def get_spec(self):
    spec = {}
    spec['offset'] = self.offset
    spec['bg'] = {}
    spec['peaks'] = {}
    try:
      spec['bg'] = self.bg.get_spec()
    except:
      pass

    try:
      spec['peaks'] = self.peaks.get_spec()
    except:
      pass

    return spec

  def set_spec(self, spec):
    self.clear_abg()
    self.bg = Background()
    self.peaks = Peaks()
    self.offset = spec['offset']
    self.bg.set_spec(spec['bg'])
    self.peaks.set_spec(spec['peaks'])
 
  def E(self):
    return self.EE - self.offset

  def write_fits(self):
    file = open(self.name + '.py', 'w')
    pprint(self.get_spec(), file)
    file.close()

  def crop(self, range):
    """Crop to fit within the energy window given in range"""
    upper = self.E()>=range[0]
    self.EE = self.EE[upper]
    self.data = self.data[upper]
  
    lower = self.E()<=range[1]
    self.EE = self.EE[lower]
    self.data = self.data[lower]

  def sg(self,points):
    ker = {}
    ker[5] = r_[-3,12,17,12,-3]/35
    ker[7] = r_[-2,3,6,7,6,3,-2]/21
    ker[9] = r_[-21,14,39,54,59,54,39,14,-21]/231
    ker[11] = r_[-36,9,44,69,84,89,84,69,44,9,-36]/429
    ker[13] = r_[-11,0,9,16,21,24,25,24,21,16,9,0,-11]/143
    ker[15] = r_[-78,-13,42,87,122,147,162,167,162,147,122,87,42,-13,-78]/1105
    ker[17] = r_[-21,-6,7,18,27,34,39,42,43,42,39,34,27,18,7,-6,-21]/323
    ker[19] = r_[-136,-51,24,89,144,189,224,249,264,269,264,249,224,189,144,89,24,-51,-136]/2261
    ker[21] = r_[-171,-76,9,84,149,204,249,284,309,324,309,284,204,149,84,9,-76,-171]/3059
    ker[23] = r_[-42,-21,-2,15,30,43,54,63,70,75,78,79,78,75,70,63,54,43,30,15,-2,-21,-42]/8059
    ker[25] = r_[-253,-138,-33,62,147,222,287,322,387,422,447,462,467,462,447,422,387,322,287,222,147,62,-33,-138,-253]/5175

    sg = -ker[points]
    pad_start = array([])
    pad_end = array([])
    for i in range(0,points/2):
      pad_start = r_[self.data[0], pad_start]
      pad_end = r_[pad_end, self.data[-1]]
    sg1 = r_[pad_start, self.data, pad_end]
    sg1 = convolve(sg1, sg, 'valid')
    return sg1

  def sg1(self,points):
    ker = {}
    ker[5] = r_[-2:2.1]/10
    ker[7] = r_[-3:3.1]/28
    ker[9] = r_[-4:4.1]/60
    ker[11] = r_[-5:5.1]/110
    ker[13] = r_[-6:6.1]/182
    ker[15] = r_[-7:7.1]/280
    ker[17] = r_[-8:8.1]/408
    ker[19] = r_[-9:9.1]/570
    ker[21] = r_[-10:10.1]/770
    ker[23] = r_[-11:11.1]/1012
    ker[25] = r_[-12:12.1]/1300

    sg = ker[points]/(self.E()[0]-self.E()[1])
    pad_start = array([])
    pad_end = array([])
    for i in range(0,points/2):
      pad_start = r_[self.data[0], pad_start]
      pad_end = r_[pad_end, self.data[-1]]
    sg1 = r_[pad_start, self.data, pad_end]
    sg1 = convolve(sg1, sg, 'valid')
    return sg1

  def noabg(self):
    """return spectrum with bg subtracted"""
    return self.data - self.abg(self.E())

  def nobg(self):
    """return spectrum with bg subtracted"""
    return self.data - self.abg(self.E()) - self.bg(self.E(), self.data)

  def residuals(self):
    """return residuals from current fit"""
    spectrum = self.nobg() - self.peaks(self.E())
    return spectrum

  def full_fit(self):
    """return the current fit including peaks and bg"""
    return self.peaks(self.E())+self.bg(self.E(),self.data)

  def add_peak(self, name, variables, values, penalties, f):
    """Add a peak with fitting params, using function f"""
    self.peaks.add_peak(Peak(name, variables, values, penalties, f))
    self.peaks.optimize_fit(self.E(), self.nobg())

  def guess_abg_from_spec(self, spec):
    values = copy(spec['values'])
    ranges = copy(spec['ranges'])
    self.guess_abg(spec['name'], spec['variables'], values, ranges, eval(spec['function']), eval(spec['penalty_function']))

  def guess_abg(self, name, variables, values, ranges, function, penalty_function):
    penalties = []
    for range in ranges:
      penalties.append(Penalty(range,penalty_function))

    self.abg = AnalyticBackground(name, variables, values, penalties, function)
    self.abg.optimize_fit(self.E(), self.data)
 
  def guess_bg_from_spec(self, spec):
    values = copy(spec['values'])
    ranges = copy(spec['ranges'])
    self.guess_bg(spec['name'], spec['variables'], values, ranges, eval(spec['function']), eval(spec['penalty_function']))

  def guess_bg(self, name, variables, values, ranges, function, penalty_function):
    penalties = []
    for range in ranges:
      penalties.append(Penalty(range,penalty_function))

    self.bg = Background(name, variables, values, penalties, function, values[1]*2)
    self.bg.optimize_fit(self.E(), self.noabg())
  
  def guess_peak_from_spec(self, spec):
    values = copy(spec['values'])
    ranges = copy(spec['ranges'])
    
    self.guess_peak(spec['name'], spec['variables'], values, ranges, eval(spec['function']), eval(spec['penalty_function']))

  def guess_peak(self, name, variables, values, ranges, function, penalty_function):
    penalties = []
    for range in ranges:
      penalties.append(Penalty(range,penalty_function))

    self.add_peak(name, variables, values, penalties, function)
  
  def plot_f_of_data(self,f,scale=1.0,axes=None):
    if axes==None:
      axes=gca()
    lines = axes.plot(self.E(), f(self.data/scale),'--',label=f.func_name + "(" + self.name + ")")
    
  def plot(self,scale=1.0,axes=None):
    """Plot the spectra"""
    if axes==None:
      axes = gca()
    lines = axes.plot(self.E(), self.data/scale,':o', markersize=3, label=self.name,picker=2)

    for line in lines:
      line.spectrum = self

  def plot_sg1(self,points=5,scale=1.0,axes=None):
    """Plot the Savitski Golay derviative of the spectra"""
    if axes==None:
      axes = gca()
    lines = axes.plot(self.E(), self.sg1(points)/scale,':o', markersize=3, label=self.name+'-sg1(%d)'%(points,))

  def plot_nobg(self,scale=1.0, axes=None):
    """Plot sepctrum with bg subtracted"""
    if axes==None:
      axes = gca()
    lines = axes.plot(self.E(), self.nobg()/scale,':o', markersize=3, label=self.name+'-nobg',picker=2)
    for line in lines:
      line.spectrum = self

  def plot_peaks(self,scale=1.0, axes=None):
    """Plot the peaks summed together"""
    if axes==None:
      axes = gca()
    lines = axes.plot(self.E(), self.peaks(self.E())/scale,label=self.name+'-peaks')
    for line in lines:
      line.spectrum = self

  def plot_individual_peaks(self,scale=1.0, axes=None):
    """Plot the peaks as individual functions"""
    if axes==None:
      axes = gca()
    i = 0
    colors = ['b','g','r','c','m','y','k']
    for peak in self.peaks.peak_list:
      line = axes.fill_between(self.E(),0, peak(self.E())/scale,label=self.name+'-peak%d'%i,alpha=0.4,color=colors[i])
      line.spectrum = self
      peak.line = line
      i += 1

  def plot_individual_peaks_bg(self,scale=1.0, axes=None):
    """Plot the peaks as individual functions"""
    if axes==None:
      axes = gca()
    b = self.bg(self.E(), self.data)
    i = 0
    colors = ['b','g','r','c','m','y','k']
    for peak in self.peaks.peak_list:
      line = axes.fill_between(self.E(), b/scale, (peak(self.E())+b)/scale,label=self.name+'-peak%d'%i,alpha=0.4,color=colors[i])
      line.spectrum = self
      peak.line = line
      i += 1

  def plot_residuals(self,scale=1.0, axes=None):
    """plot the residuals"""
    if axes==None:
      axes = gca()
    lines = axes.plot(self.E(), self.residuals()/scale,label=self.name+'-residuals')
    for line in lines:
      line.spectrum = self

  def plot_bg(self,scale=1.0, axes=None):
    """Plot the bg"""
    if axes==None:
      axes = gca()
    lines = axes.plot(self.E(),(self.abg(self.E()) + self.bg(self.E(), self.data))/scale,label=self.name+'-bg')
    for line in lines:
      line.spectrum = self
    self.bg.line = lines[0]

  def plot_abg(self,scale=1.0, axes=None):
    """Plot the bg"""
    if axes==None:
      axes = gca()
    lines = axes.plot(self.E(),(self.abg(self.E()))/scale,label=self.name+'-abg')
    for line in lines:
      line.spectrum = self

  def plot_kbg(self,scale=1.0, axes=None):
    """Plot the bg"""
    if axes==None:
      axes = gca()
    lines = axes.plot(self.E(),(self.bg(self.E(), self.data))/scale,label=self.name+'-kbg')
    for line in lines:
      line.spectrum = self

  def plot_full_fit(self,scale=1.0, axes=None):
    """Plot the current fit including peaks and bg"""
    if axes==None:
      axes = gca()
    lines = axes.plot(self.E(),(self.peaks(self.E())+self.bg(self.E(),self.data))/scale,label=self.name+'-fullfit')
    for line in lines:
      line.spectrum = self

  def plot_full_summary(self,scale=1.0, axes=None,displayParams=True):
    """Plot the spectrum, the bg, and the full fit"""
    if axes==None:
      axes = gca()
    self.plot(scale, axes)
    self.plot_bg(scale, axes)
    self.plot_abg(scale, axes)
    self.plot_kbg(scale, axes)
    self.plot_full_fit(scale, axes)
    self.plot_individual_peaks_bg(scale, axes)

    if displayParams:
      if isinstance(self.bg, Background):
        spec = self.bg.get_spec()
        label = spec['name'] + "\n"
        label += spec['function']
        for var, val, range in zip(spec['variables'], spec['values'], spec['ranges']):
          label += "\n$%s = %0.2f\in\, (%0.2f,%0.2f)$"%(var,val,range[0],range[1])
        label += ""
        mu = self.E()[0]
        A = self.bg(self.E(),self.data)[0]
        EE = self.E()[::-1]
        col = self.bg.line.get_color()
        an = axes.annotate(label,
                      xy=(mu, A), xycoords='data',
                      xytext=(0,60), textcoords='offset points',
                      arrowprops=dict(arrowstyle='fancy',connectionstyle="arc3,rad=0.0",fc='%s'%col,alpha=0.6),
                      bbox=dict(boxstyle="roundtooth",fc='%s'%col,alpha=0.2),
                      fontsize='small',
                      family='sans-serif')
        an.draggable()
        an.fit_object = self.bg
        
      for peak in self.peaks:
        spec = peak.get_spec()
        label = spec['name'] + '\n'
        label += spec['function']
        for var, val, range in zip(spec['variables'], spec['values'], spec['ranges']):
          label += "\n$%s = %0.2f\in\, (%0.2f,%0.2f)$"%(var,val,range[0],range[1])
        label += ""
        mu = spec['values'][1]
        A = spec['values'][0]/scale
        EE = self.E()[::-1]
        back = self.bg(self.E(), self.data)[::-1]
        if mu > EE[-1]:
          b = back[-1]
        else:
          b = back[EE>=mu][0]
        col = peak.line.get_facecolor()
        an = axes.annotate(label,
                      xy=(mu, A + b), xycoords='data',
                      xytext=(0,60), textcoords='offset points',
                      arrowprops=dict(arrowstyle='fancy',connectionstyle="arc3,rad=0.0"),
                      bbox=dict(boxstyle="round"),
                      fontsize='small',
                      family='sans-serif')
        an.get_bbox_patch().set_facecolor(col[0])
        an.get_bbox_patch().set_alpha(col[0][3])
        an.draggable()
        an.fit_object = peak
  
  def plot_full_summary_nobg(self,scale=1.0, axes=None, displayParams=True):
    """Plot the spectrum, the bg, and the full fit"""
    if axes==None:
      axes = gca()
    self.plot_nobg(scale, axes)
    self.plot_peaks(scale, axes)
    self.plot_individual_peaks(scale, axes)

    if displayParams:
      for peak in self.peaks:
        spec = peak.get_spec()
        label = spec['name'] + '\n'
        label += spec['function']
        for var, val, range in zip(spec['variables'], spec['values'], spec['ranges']):
          label += "\n$%s = %0.2f\in\, (%0.2f,%0.2f)$"%(var,val,range[0],range[1])
        label += ""
        mu = spec['values'][1]
        A = spec['values'][0]/scale
        EE = self.E()[::-1]
        col = peak.line.get_facecolor()
        an = axes.annotate(label,
                      xy=(mu, A), xycoords='data',
                      xytext=(0,60), textcoords='offset points',
                      arrowprops=dict(arrowstyle='fancy',connectionstyle="arc3,rad=0.0"),
                      bbox=dict(boxstyle="round"),
                      fontsize='small',
                      family='sans-serif')
        an.get_bbox_patch().set_facecolor(col[0])
        an.get_bbox_patch().set_alpha(col[0][3])
        an.draggable()
        an.fit_object = peak
  
  def __call__(self):
    return self.data
  
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
      data = np.genfromtxt(filename,skip_header=12)
      #print data
      spectra[name].EE = data[:,0]
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
  
  for file in files:
    m = re.search('('+file_name_filter+'.*)\.xas$',file)
    if (m):
      filename = file
      name = filename
      spectra[name] = Spectrum()
      spectra[name].name = name
      data = np.genfromtxt(filename,skip_header=5)
      #print data
      spectra[name].EE = data[:,0]
      spectra[name].data = data[:,1:]
      
  return spectra

def load_AugerScan_text_files(file_name_filter):
  """helper function to load a bunch of spectra from AugerScan text files in the current dir"""
  files = os.listdir('.')

  special_chars = r'[ \t\(\"\(\)\*\&\^\%\$\#\@\!\_\+\-\=\[\]\{\}\|\;\'\"\:\/\.\,\<\>\\\?]'
  
  spectra = {}
  E=[]
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
          spectrum.E=array([])
          spectrum.data=array([])
          spectra[label] = spectrum
        elif header:
          continue
        elif (re.search('^[\d\.]+\s\d+\s*$',line)):
          m = re.search('^([\d\.]+)\s(\d+)',line)
          #print m.group(0)
          spectra[label].EE = append(spectra[label].E, double(m.group(1)))
          spectra[label].data = append(spectra[label].data, double(m.group(2)))

  return spectra
