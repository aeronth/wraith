import os
import re
from pylab import *
from scipy.optimize import leastsq
from fitting_machinery import *

PHI_sensitivity = {'Li 1s':    0.025, 'Be 1s':    0.074, 'B 1s':     0.156, 'C 1s':     0.296, 
                   'N 1s':     0.477, 'O 1s':     0.711, 'F 1s':     1.000, 'Ne 1s':    1.340, 
                   'Na 1s':    1.685, 'Mg 2p':    0.252, 'Al 2p':    0.193, 'Si 2p':    0.283, 
                   'P 2p':     0.412, 'S 2p':     0.570, 'Cl 2p':    0.770, 'Ar 2p':    1.011, 
                   'K 2p':     1.300, 'Ca 2p':    1.634, 'Sc 2p':    1.678, 'Ti 2p':    0.1789, 
                   'V 2p':     1.912, 'Cr 2p':    2.010, 'Mn 2p':    2.420, 'Fe 2p':    2.686, 
                   'Co 2p':    3.255, 'Ni 2p':    3.653, 'Cu 2p':    4.789, 'Zn 2p3/2': 3.3543, 
                   'Ga 2p3/2': 3.341, 'Ge 2p3/2': 3.100, 'As 3d':    0.570, 'Se 3d':    0.722, 
                   'Br 3d':    0.895, 'Kr 3d':    1.096, 'Rb 3d':    1.316, 'Sr 3d':    1.578, 
                   'Y 3d':     1.867, 'Zr 3d':    2.216, 'Nb 3d':    2.517, 'Mo 3d':    2.867, 
                   'Tc 3d':    3.266, 'Ru 3d':    3.696, 'Rh 3d':    4.179, 'Pd 3d':    5.198, 
                   'Ag 3d':    3.444, 'Cd 3d5/2': 3.777, 'In 3d5/2': 4.095, 'Sn 3d5/2': 4.473, 
                   'Sb 3d5/2': 4.925, 'Te 3d5/2': 4.925, 'I 3d5/2':  5.337, 'Xe 3d5/2': 5.702, 
                   'Cs 3d5/2': 6.032, 'Ba 3d5/2': 6.361, 'La 3d':    7.708, 'Ce 3d':    7.399, 
                   'Pr 3d':    6.356, 'Nd 3d':    4.597, 'Pm 3d':    3.754, 'Sm 3d5/2': 2.907, 
                   'Eu 4d':    2.210, 'Gd 4d':    2.207, 'Tb 4d':    2.201, 'Dy 4d':    2.198, 
                   'Ho 4d':    2.189, 'Er 4d':    2.184, 'Tm 4d':    2.172, 'Yb 4d':    2.169, 
                   'Lu 4f':    2.156, 'Hf 4f':    2.221, 'Ta 4f':    2.589, 'W 4f':     2.959, 
                   'Re 4f':    3.327, 'Os 4f':    3.747, 'Ir 4f':    4.217, 'Pt 4f':    4.674, 
                   'Au 4f':    5.240, 'Hg 4f':    5.797, 'Tl 4f':    6.447, 'Pb 4f':    6.968, 
                   'Bi 4f':    7.632, 'Th 4f7/2': 7.496, 'U 4f7/2':  8.476} 


PHI_sensitivity_a={'Li 1s':      0.025, 'Be 1s':      0.074, 'B 1s':       0.156, 'C 1s':       0.296, 
                   'N 1s':       0.477, 'O 1s':       0.711, 'F 1s':       1.000, 'Ne 1s':      1.340, 
                   'Na 1s':      1.685, 'Mg 2p':1/1.5*  0.252, 'Al 2p':1/1.5*  0.193, 'Si 2p':1/1.5*  0.283, 
                   'P 2p':1/1.5*   0.412, 'S 2p':1/1.5*   0.570, 'Cl 2p':1/1.5*  0.770, 'Ar 2p':1/1.5*  1.011, 
                   'K 2p':1/1.5*   1.300, 'Ca 2p':1/1.5*  1.634, 'Sc 2p':1/1.5*  1.678, 'Ti 2p':1/1.5*  0.1789, 
                   'V 2p':1/1.5*   1.912, 'Cr 2p':1/1.5*  2.010, 'Mn 2p':1/1.5*  2.420, 'Fe 2p':1/1.5*  2.686, 
                   'Co 2p':1/1.5*  3.255, 'Ni 2p':1/1.5*  3.653, 'Cu 2p':1/1.5*  4.789, 'Zn 2p':      3.3543, 
                   'Ga 2p':      3.341, 'Ge 2p':      3.100, 'As 3d':      0.570, 'Se 3d':      0.722, 
                   'Br 3d':1/1.666*0.895, 'Kr 3d':1/1.666*1.096, 'Rb 3d':1/1.666*1.316, 'Sr 3d':1/1.666*1/1.578, 
                   'Y 3d':1/1.666* 1.867, 'Zr 3d':1/1.666*2.216, 'Nb 3d':1/1.666*2.517, 'Mo 3d':1/1.666*2.867, 
                   'Tc 3d':1/1.666*3.266, 'Ru 3d':1/1.666*3.696, 'Rh 3d':1/1.666*4.179, 'Pd 3d':1/1.666*5.198, 
                   'Ag 3d':1/1.666*3.444, 'Cd 3d':      3.777, 'In 3d':      4.095, 'Sn 3d':      4.473, 
                   'Sb 3d':      4.925, 'Te 3d':      4.925, 'I 3d':       5.337, 'Xe 3d':      5.702, 
                   'Cs 3d':      6.032, 'Ba 3d':      6.361, 'La 3d':      7.708, 'Ce 3d':      7.399, 
                   'Pr 3d':1/1.666*6.356, 'Nd 3d':1/1.666*4.597, 'Pm 3d':1/1.666*3.754, 'Sm 3d':      2.907, 
                   'Eu 4d':1/1.666*2.210, 'Gd 4d':1/1.666*2.207, 'Tb 4d':1/1.666*2.201, 'Dy 4d':1/1.666*2.198, 
                   'Ho 4d':1/1.666*2.189, 'Er 4d':1/1.666*2.184, 'Tm 4d':1/1.666*2.172, 'Yb 4d':1/1.666*2.169, 
                   'Lu 4f':1/1.75* 2.156, 'Hf 4f':1/1.75* 2.221, 'Ta 4f':1/1.75* 2.589, 'W 4f':1/1.75*  2.959, 
                   'Re 4f':1/1.75* 3.327, 'Os 4f':1/1.75* 3.747, 'Ir 4f':1/1.75* 4.217, 'Pt 4f':1/1.75* 4.674, 
                   'Au 4f':1/1.75* 5.240, 'Hg 4f':1/1.75* 5.797, 'Tl 4f':1/1.75* 6.447, 'Pb 4f':1/1.75* 6.968, 
                   'Bi 4f':1/1.75* 7.632, 'Th 4f':      7.496, 'U 4f':       8.476} 


class Spectrum:
  """A place to keep a spectrum"""
  
  def __init__(self):
    """Initialize stuff"""
    self.data = ones(1)
    self.E = ones(1)
    self.bg = lambda E, data: zeros(E.size)
    self.peaks = Peaks()

  def crop(self, range):
    """Crop to fit within the energy window given in range"""
    upper = self.E>=range[0]
    self.E = self.E[upper]
    self.data = self.data[upper]
  
    lower = self.E<=range[1]
    self.E = self.E[lower]
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

    sg = ker[points]/(self.E[0]-self.E[1])
    pad_start = array([])
    pad_end = array([])
    for i in range(0,points/2):
      pad_start = r_[self.data[0], pad_start]
      pad_end = r_[pad_end, self.data[-1]]
    sg1 = r_[pad_start, self.data, pad_end]
    sg1 = convolve(sg1, sg, 'valid')
    return sg1

  def nobg(self):
    """return spectrum with bg subtracted"""
    return self.data - self.bg(self.E, self.data)

  def residuals(self):
    """return residuals from current fit"""
    spec = self.nobg() - self.peaks(self.E)
    return spec

  def full_fit(self):
    """return the current fit including peaks and bg"""
    return self.peaks(self.E)+self.bg(self.E,self.data)

  def add_peak(self, name, params, penalties, f):
    """Add a peak with fitting params, using function f"""
    self.peaks.add_peak(Peak(name, params, penalties, f))
    self.peaks.optimize_fit(self.E, self.nobg())

  def residual_m0(self):
    """return the zeroth moment of the residuals"""
    return sum(abs(self.residuals()))

  def residual_m1(self):
    """return the first moment of the residuals"""
    res = abs(self.residuals())
    return self.E[round(sum(r_[0:size(res)]*res)/sum(res))]

  def residual_m2(self):
    """return the second moment of the residuals"""
    res = abs(self.residuals())
    index_sigma = sqrt((sum((r_[0:size(res)]-self.residual_m1())**2 * res)/sum(res)))
    sigma = index_sigma * abs(mean(diff(self.E)))
    return sigma

  def simple_guess_s_peak(self):
    return

  def simple_guess_p_peak(self):
    return

  def simple_guess_d_peak(self):
    return

  def simple_guess_f_peak(self):
    return

  def simple_guess_peak(self):
    return

  def guess_bg_from_spec(self, bg_spec):
    initial_values = copy(bg_spec['initial_values'])
    ranges = copy(bg_spec['ranges'])
    self.guess_bg(bg_spec['name'], initial_values, ranges, eval(bg_spec['function']), eval(bg_spec['penalty_function']))

  def guess_bg(self, name, initial_values, ranges, function, penalty_function):
    penalties = []
    for range in ranges:
      penalties.append(Penalty(range,penalty_function))

    self.bg = Background(name, initial_values, penalties, function, initial_values[1]/2)
    self.bg.optimize_fit(self.E, self.data)
  
  def guess_peak_from_spec(self, peak_spec):
    initial_values = copy(peak_spec['initial_values'])
    initial_values[0] = max(self.nobg())*initial_values[0]

    ranges = copy(peak_spec['ranges'])
    ranges[0][0] = max(self.nobg())*ranges[0][0]
    ranges[0][1] = max(self.nobg())*ranges[0][1]
    
    self.guess_peak(peak_spec['name'], initial_values, ranges, eval(peak_spec['function']), eval(peak_spec['penalty_function']))

  def guess_peak(self, name, initial_values, ranges, function, penalty_function):
    penalties = []
    for range in ranges:
      penalties.append(Penalty(range,penalty_function))

    self.add_peak(name, initial_values, penalties, function)
  
  def plot(self,scale=1.0):
    """Plot the spectra"""
    plot(self.E, self.data/scale, 'k.')

  def plot_peaks(self,scale=1.0):
    """Plot the peaks summed together"""
    plot(self.E, self.peaks(self.E)/scale)

  def plot_individual_peaks(self,scale=1.0):
    """Plot the peaks as individual functions"""
    for peak in self.peaks.peak_list:
      plot(self.E, peak(self.E)/scale)

  def plot_individual_peaks_bg(self,scale=1.0):
    """Plot the peaks as individual functions"""
    b = self.bg(self.E, self.data)
    for peak in self.peaks.peak_list:
      plot(self.E, (peak(self.E)+b)/scale)

  def plot_residuals(self,scale=1.0):
    """plot the residuals"""
    plot(self.E, self.residuals()/scale)

  def plot_nobg(self,scale=1.0):
    """Plot sepctrum with bg subtracted"""
    plot(self.E,self.nobg()/scale,'k.')

  def plot_bg(self,scale=1.0):
    """Plot the bg"""
    plot(self.E,self.bg(self.E, self.data)/scale)

  def plot_full_fit(self,scale=1.0):
    """Plot the current fit including peaks and bg"""
    plot(self.E,(self.peaks(self.E)+self.bg(self.E,self.data))/scale, 'k:')

  def plot_full_summary(self,scale=1.0):
    """Plot the spectrum, the bg, and the full fit"""
    self.plot(scale)
    self.plot_bg(scale)
    self.plot_full_fit(scale)
    self.plot_individual_peaks_bg(scale)
    ylim(-0.1*max(self.data)/scale,1.2*max(self.data)/scale)
    i = 0
    for peak in self.peaks:
      p = peak.params
      x_factor = 0.8 - i/2 * 0.6
      y_factor = 0.75 - i%2 * 0.25
      x = xlim()[0] + x_factor * diff(xlim())
      y = ylim()[0] + y_factor * diff(ylim())
      #print "x_factor: %f\nx: %f\ny_factor: %f\ny: %f"%(x_factor,x,y_factor,y)
      label = "$Area=%0.0f$\n$\mu=%0.2f$\n$A=%0.2f$\n$\sigma=%0.2f$\n$m=%0.2f$"%(p[0]*p[2],p[1],p[0],p[2],p[3])
      if size(p)>4:
        label += "\n$\Delta=%0.2f$"%(p[6])
      text(x,
           y,
           label,
           bbox=dict(facecolor='red', alpha=0.2))
      plot(r_[x, p[1]], r_[y, (p[0]+min(self.data))/scale])
      i += 1

  def plot_full_summary_nobg(self,scale=1.0):
    """Plot the spectrum, the bg, and the full fit"""
    self.plot_nobg(scale)
    self.plot_individual_peaks(scale)
    ylim(-0.1*max(self.data)/scale,1.2*max(self.data)/scale)
    i = 0
    for peak in self.peaks:
      p = peak.params
      x_factor = 0.8 - i/2 * 0.6
      y_factor = 0.75 - i%2 * 0.25
      x = xlim()[0] + x_factor * diff(xlim())
      y = ylim()[0] + y_factor * diff(ylim())
      #print "x_factor: %f\nx: %f\ny_factor: %f\ny: %f"%(x_factor,x,y_factor,y)
      label = "$Area=%0.0f$\n$\mu=%0.2f$\n$A=%0.2f$\n$\sigma=%0.2f$\n$m=%0.2f$"%(p[0]*p[2],p[1],p[0],p[2],p[3])
      if size(p)>4:
        label += "\n$\Delta=%0.2f$"%(p[6])
      text(x,
           y,
           label,
           bbox=dict(facecolor='red', alpha=0.2))
      plot(r_[x, p[1]], r_[y, (p[0]+min(self.data))/scale])
      i += 1

  def summary_as_string(self):
    """Return a string that is representative of the full fit to the spectrum"""
    for peak in self.peaks:
      p = peak.params
      #print "x_factor: %f\nx: %f\ny_factor: %f\ny: %f"%(x_factor,x,y_factor,y)
      label = "Area=%0.0f\nPeak Position=%0.2f$\nA=%0.2f\nPeak Width=%0.2f\nGL factor=%0.2f"%(p[0]*p[2],p[1],p[0],p[2],p[3])
      if size(p)>4:
        label += "\n$\Delta=%0.2f$"%(p[6])

    return label
# 
  def __call__(self):
    return self.data

def crop_spectra(spectra,range):
  """helper function to crop a bunch of spectra"""
  for key in sorted(spectra.keys()):
    spectra[key].crop(range)

def tougaard_3_bg_subtract_spectrum(spectrum, B_range, C_range, D_range):
  """helper function to add a tougaard bg to a bunch of spectra"""
  params = [average(B_range), average(C_range), average(D_range)]
  penalties = [Penalty(B_range, exp_penalty), Penalty(C_range, exp_penalty), Penalty(D_range, exp_penalty)]
  spectrum.bg = Background("tougaard3", ['B','C','D'], params, penalties, K3, C_range[1]/2)
  spectrum.bg.optimize_fit(spectrum.E, spectrum.data)

def tougaard_3_bg_subtract_spectra(spectra, B_range, C_range, D_range):
  """helper function to add a tougaard bg to a bunch of spectra"""
  for key in sorted(spectra.keys()):
    spectrum = spectra[key]
    params = [average(B_range), average(C_range), average(D_range)]
    penalties = [Penalty(B_range, exp_penalty), Penalty(C_range, exp_penalty), Penalty(D_range, exp_penalty)]
    spectrum.bg = Background("tougaard3", ['B','C','D'], params, penalties, K3, C_range[1]/2)
    spectrum.bg.optimize_fit(spectrum.E, spectrum.data)

def tougaard_bg_subtract_spectrum(spectrum, B_range, C_range):
  """helper function to add a tougaard bg to a bunch of spectra"""
  params = [average(B_range), average(C_range)]
  penalties = [Penalty(B_range, exp_penalty), Penalty(C_range, exp_penalty)]
  spectrum.bg = Background("tougaard", ['B','C'], params, penalties, K, C_range[1]/2)
  spectrum.bg.optimize_fit(spectrum.E, spectrum.data)

def tougaard_bg_subtract_spectra(spectra, B_range, C_range):
  """helper function to add a tougaard bg to a bunch of spectra"""
  for key in sorted(spectra.keys()):
    spectrum = spectra[key]
    params = [average(B_range), average(C_range)]
    penalties = [Penalty(B_range, exp_penalty), Penalty(C_range, exp_penalty)]
    spectrum.bg = Background("tougaard", ['B','C'], params, penalties, K, C_range[1]/2)
    spectrum.bg.optimize_fit(spectrum.E, spectrum.data)

def plot_spectra_summary(spectra):
  """helper function to plot out fitting results for a bunch of spectra"""
  
  offset1 = 0

  keys = sorted(spectra.keys())

  for key in keys:
    spectrum = spectra[key]

    scale_factor = max(spectrum.nobg())
    normed = spectrum.nobg()/scale_factor
    plot(spectrum.E,normed+offset1,'b-')
    plot(spectrum.E,spectrum.peaks(spectrum.E)/scale_factor + offset1,'k-')
    for peak in spectrum.peaks.peak_list:
      plot(spectrum.E, peak(spectrum.E)/scale_factor + offset1)
    
    plot(spectrum.E, spectrum.residuals()/scale_factor + offset1 - 0.1,'r-')
    
    offset1 += 1.2
  
  grid('on')
  yticks(arange(0,size(keys)*1.2,1.2),keys)
  xlim([max(spectrum.E), min(spectrum.E)])

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
      spectra[name].E = data[:,0]
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
      spectra[name].E = data[:,0]
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
      spectra[name].E = data[:-1,0]
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
      spectra[name].E = data[:,0]
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
          spectra[label].E = append(spectra[label].E, double(m.group(1)))
          spectra[label].data = append(spectra[label].data, double(m.group(2)))

  return spectra
