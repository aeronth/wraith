import os
import re
from pylab import *
from scipy.optimize import leastsq

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
    self.simple_guess_peak()

  def simple_guess_p_peak(self):
    self.simple_guess_split_peak(r_[0.49,0.51], r_[0.1,15])
    
  def simple_guess_d_peak(self):
    self.simple_guess_split_peak(r_[0.66,0.68], r_[0.1,15])

  def simple_guess_f_peak(self):
    self.simple_guess_split_peak(r_[0.74,0.76], r_[0.1,15])

  def simple_guess_peak(self):
    a_guess = max(self.residuals())
    mu_guess = self.residual_m1()
    sigma_guess = self.residual_m2()
    m_guess = 0.5
    guess_parm = r_[a_guess,
                    mu_guess,
                    sigma_guess,
                    m_guess]
    print guess_parm
    a_penalty = Penalty(r_[0,20*a_guess],exp_bounded_penalty)
    mu_penalty = Penalty(r_[self.E[0],self.E[-1]],exp_bounded_penalty)
    sigma_penalty = Penalty(r_[0,2*sigma_guess],exp_bounded_penalty)
    m_penalty = Penalty(r_[0,1],exp_bounded_penalty)
    penalties = [a_penalty,
                 mu_penalty,
                 sigma_penalty,
                 m_penalty]

    self.add_peak(str(mu_guess), guess_parm, penalties, gl)

  def simple_guess_split_peak(self, ratio_window, split_window):
    a_guess = max(self.residuals())
    mu_guess = self.residual_m1()
    sigma_guess = self.residual_m2()
    m_guess = 0.5
    ratio_a_guess = average(ratio_window)
    ratio_area_guess = average(ratio_window)
    split_guess = average(split_window)
    guess_parm = r_[a_guess,
                    mu_guess,
                    sigma_guess,
                    m_guess,
                    ratio_a_guess,
                    ratio_area_guess,
                    split_guess]

    print guess_parm
    a_penalty = Penalty(r_[0,20*a_guess],exp_bounded_penalty)
    mu_penalty = Penalty(r_[self.E[0],self.E[-1]],exp_bounded_penalty)
    sigma_penalty = Penalty(r_[0,2*sigma_guess],exp_bounded_penalty)
    m_penalty = Penalty(r_[0,1],exp_bounded_penalty)
    ratio_a_penalty = Penalty(r_[0,1],exp_bounded_penalty)
    ratio_area_penalty = Penalty(ratio_window,exp_bounded_penalty)
    split_penalty = Penalty(split_window,exp_bounded_penalty)
    penalties = [a_penalty,
                 mu_penalty,
                 sigma_penalty,
                 m_penalty,
                 ratio_a_penalty,
                 ratio_area_penalty,
                 split_penalty]

    self.add_peak(str(mu_guess), guess_parm, penalties, spin_split_gl)

  def guess_peak_(self, name, ranges, function):
    guess_parm = array([])
    penalties = []
    for range in ranges:
      guess_parm = r_[guess_parm, average(range)]
      penalties.append(Penalty(range,exp_bounded_penalty))

    self.add_peak(name, guess_parm, penalties, function)
  
  def guess_peak(self, name, a_window, mu_window, sigma_window, m_window):
    a_guess = max(self.residuals())
    mu_guess = average(mu_window) 
    sigma_guess = average(sigma_window) 
    m_guess = average(m_window)
    guess_parm = r_[a_guess,
                    mu_guess,
                    sigma_guess,
                    m_guess]
                                   
    a_penalty = Penalty(a_window,exp_bounded_penalty)
    mu_penalty = Penalty(mu_window,exp_bounded_penalty)
    sigma_penalty = Penalty(sigma_window,exp_bounded_penalty)
    m_penalty = Penalty(m_window,exp_bounded_penalty)
    penalties = [a_penalty,
                 mu_penalty,
                 sigma_penalty,
                 m_penalty]
    
    self.add_peak(name, guess_parm, penalties, gl)

  def guess_split_peak(self, name, a_window, mu_window, sigma_window, m_window, ratio_a_window, ratio_area_window, split_window):
    a_guess = max(self.residuals())
    mu_guess = average(mu_window) 
    sigma_guess = average(sigma_window)
    m_guess = average(m_window)
    ratio_a_guess = average(ratio_a_window)
    ratio_area_guess = average(ratio_area_window)
    split_guess = average(split_window)
    guess_parm = r_[a_guess,
                    mu_guess,
                    sigma_guess,
                    m_guess,
                    ratio_a_guess,
                    ratio_area_guess,
                    split_guess]
                                   
    a_penalty = Penalty(a_window,exp_bounded_penalty)
    mu_penalty = Penalty(mu_window,exp_bounded_penalty)
    sigma_penalty = Penalty(sigma_window,exp_bounded_penalty)
    m_penalty = Penalty(m_window,exp_bounded_penalty)
    ratio_a_penalty = Penalty(ratio_a_window,exp_bounded_penalty)
    ratio_area_penalty = Penalty(ratio_area_window,exp_bounded_penalty)
    split_penalty = Penalty(split_window,exp_bounded_penalty)
    penalties = [a_penalty,
                 mu_penalty,
                 sigma_penalty,
                 m_penalty,
                 ratio_a_penalty,
                 ratio_area_penalty,
                 split_penalty]
    
    self.add_peak(name, guess_parm, penalties, spin_split_gl)

  def plot(self,scale=1.0):
    """Plot the spectra"""
    plot(self.E, self.data/scale, 'k.')
    #xlim([max(self.E), min(self.E)])

  def plot_peaks(self,scale=1.0):
    """Plot the peaks summed together"""
    plot(self.E, self.peaks(self.E)/scale)
    xlim([max(self.E), min(self.E)])

  def plot_individual_peaks(self,scale=1.0):
    """Plot the peaks as individual functions"""
    for peak in self.peaks.peak_list:
      plot(self.E, peak(self.E)/scale)
    xlim([max(self.E), min(self.E)])

  def plot_individual_peaks_bg(self,scale=1.0):
    """Plot the peaks as individual functions"""
    b = self.bg(self.E, self.data)
    for peak in self.peaks.peak_list:
      plot(self.E, (peak(self.E)+b)/scale)

    xlim([max(self.E), min(self.E)])

  def plot_residuals(self,scale=1.0):
    """plot the residuals"""
    plot(self.E, self.residuals()/scale)
    xlim([max(self.E), min(self.E)])

  def plot_nobg(self,scale=1.0):
    """Plot sepctrum with bg subtracted"""
    plot(self.E,self.nobg()/scale,'k.')
    xlim([max(self.E), min(self.E)])

  def plot_bg(self,scale=1.0):
    """Plot the bg"""
    plot(self.E,self.bg(self.E, self.data)/scale)
    xlim([max(self.E), min(self.E)])

  def plot_full_fit(self,scale=1.0):
    """Plot the current fit including peaks and bg"""
    plot(self.E,(self.peaks(self.E)+self.bg(self.E,self.data))/scale, 'k:')
    xlim([max(self.E), min(self.E)])

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
  penalties = [Penalty(B_range, exp_bounded_penalty), Penalty(C_range, exp_bounded_penalty), Penalty(D_range, exp_bounded_penalty)]
  spectrum.bg = Background("tougaard3", ['B','C','D'], params, penalties, K3, C_range[1]/2)
  spectrum.bg.optimize_fit(spectrum.E, spectrum.data)

def tougaard_3_bg_subtract_spectra(spectra, B_range, C_range, D_range):
  """helper function to add a tougaard bg to a bunch of spectra"""
  for key in sorted(spectra.keys()):
    spectrum = spectra[key]
    params = [average(B_range), average(C_range), average(D_range)]
    penalties = [Penalty(B_range, exp_bounded_penalty), Penalty(C_range, exp_bounded_penalty), Penalty(D_range, exp_bounded_penalty)]
    spectrum.bg = Background("tougaard3", ['B','C','D'], params, penalties, K3, C_range[1]/2)
    spectrum.bg.optimize_fit(spectrum.E, spectrum.data)

def tougaard_bg_subtract_spectrum(spectrum, B_range, C_range):
  """helper function to add a tougaard bg to a bunch of spectra"""
  params = [average(B_range), average(C_range)]
  penalties = [Penalty(B_range, exp_bounded_penalty), Penalty(C_range, exp_bounded_penalty)]
  spectrum.bg = Background("tougaard", ['B','C'], params, penalties, K, C_range[1]/2)
  spectrum.bg.optimize_fit(spectrum.E, spectrum.data)

def tougaard_bg_subtract_spectra(spectra, B_range, C_range):
  """helper function to add a tougaard bg to a bunch of spectra"""
  for key in sorted(spectra.keys()):
    spectrum = spectra[key]
    params = [average(B_range), average(C_range)]
    penalties = [Penalty(B_range, exp_bounded_penalty), Penalty(C_range, exp_bounded_penalty)]
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
    m = re.search('('+file_name_filter+'.*)\.xas',file)
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
 
###########################################
# Functions to feed into fitting routines 
###########################################

def K(p, E):
  """convolution kernel for Tougaard background"""
  B, C = p
  K_ = B * E / (C + E**2)**2
  K_ = K_*(K_>0)
  return K_
  
def K3(p, E):
  """convolution kernel for Tougaard background"""
  B, C, D = p
  K_ = B * E / ((C + E**2)**2 + D*E**2)
  K_ = K_*(K_>0)
  return K_
  
def voigt(p, E):
  """
    the voigt function = convolve(gaussian, lorentzian)
    p = A, mu, sigma
  """
  dE = E[1]-E[0]
  return dE * convolve( lorentzian(p,E), gaussian(p,E), 'same' )

def spin_split_gl(params, E):
  """
    spin split for the gl function
    a, mu, sigma are for peak1
    ratio_a * a == a for peak2
    ratio_sigma * sigma = sigma for peak2
    ratio_area == ratio_a * ratio_sigma can be fixed using a boundary penalty
      ratio_area -> 1/2 for p orbitals
      ratio_area -> 2/3 for d orbitals
      ratio_area -> 3/4 for f orbitals
    split is the spacing between peak1 and peak2
  """
  a, mu, sigma, m, ratio_a, ratio_area, split = params
  ratio_sigma = ratio_area / ratio_a
  return gl_(a, mu, sigma, m, E) + gl_(ratio_a * a, mu + split, ratio_sigma * sigma, m, E)

def gl_(a, mu, sigma, m, E):
  return a * exp(-2.772589 * (1 - m) * (E - mu)**2/sigma**2) / (1 + 4 * m * (E - mu)**2/sigma**2)

def gl(params, E):
  a, mu, sigma, m = params
  return gl_(a, mu, sigma, m, E)

def gl50(params, E):
  a, mu, sigma = params
  m = 0.5
  return gl_(a, mu, sigma, m, E)

def gls_(a, mu, sigma, m, E):
  return a * (1 - m) * exp(-2.772589 * (E - mu)**2/sigma**2) + m/(1 + 4 * (E - mu)**2/sigma**2) 

def gls(params, E):
  a, mu, sigma, m = params
  return gls_(a, mu, sigma, m, E)

def lorentzian_(a, mu, sigma, E):
  return a * 1/(1 + ((E - mu)/sigma)**2)

def lorentzian(params, E):
  a, mu, sigma = params
  return lorentzian(a, mu, sigma, E)

def gaussian_(a, mu, sigma, E):
  return a * exp(-((E-mu)/sigma)**2)

def gaussian(params, E):
  a, mu, sigma = params
  return gaussian_(a, mu, sigma, E)

###########################################
# Penalty functions
###########################################

def no_penalty(range, p):
  return 1.0

def exp_bounded_penalty(range, p):
  A = 60./(range[-1] - range[0])
  return 1+exp(A*-(p-range[0]))+exp(A*(p-range[-1]))

###########################################
# Fitting machinery                        
#   Penalty: a class to steer fitting
#   Background: generic convolution based background
#   Peak: a single fitting peak
#   Peaks: a collection of peaks
###########################################

class Penalty:
  """Encapsulates a penalty function to steer fitting for a parameter"""

  def __init__(self, params, f):
    """Initialize penalty function"""
    self.params = params
    self.f = f

  def __call__(self, p):
    """Penalty!"""
    return self.f(self.params, p)

class Background:
  """Implements a generic background based on a convolution kernel function"""

  def __init__(self, name, param_names, params, penalties, kernel, kernel_end):
    self.name = name
    self.param_names = param_names
    self.params = params
    self.penalties = penalties
    self.kernel = kernel
    self.kernel_end = kernel_end

  def f(self, params, E, spectrum):
    dE = E[1] - E[0]
    spectrum = spectrum - min(spectrum)
    bg = dE * convolve( spectrum, self.kernel( params, E)[::-1], 'full')
    return bg[bg.size-spectrum.size:]
  
  def residuals(self, params, E, spectrum):
    res = spectrum - self.f(params, E, spectrum)
    i = 0
    for p in params:
      res *= self.penalties[i](p)
      i += 1

    res[res<0] = res[res<0]*6
    return res
  
  def EE(self, dE):
    return arange(0, self.kernel_end, abs(dE))

  def optimize_fit(self, E, spectrum):
    offset = min(spectrum)
    spectrum = spectrum - offset
    self.dE = E[1]-E[0]
    plsq = leastsq(self.residuals, self.params, args=(self.EE(self.dE), spectrum))
    self.params = plsq[0]
    return self.f(self.params, self.EE(self.dE), spectrum) + offset
 
  def __call__(self, E, spectrum):
    offset = min(spectrum)
    spectrum = spectrum - offset
    self.dE = E[1]-E[0]
    return self.f(self.params, self.EE(self.dE), spectrum) + offset

class Peak:

  def __init__(self, name, params, penalties, f):
    self.name = name
    self.params = params
    self.penalties = penalties
    self.f = f

  def residuals(self, params, E, spectrum):
    res = spectrum - self.f(params, E)
    i = 0
    for p in params:
      res *= self.penalties[i](p)
      i += 1

    return res

  def __call__(self, E):
    return self.f(self.params, E)


class Peaks:

  def __init__(self):
    """something will go here"""
    self.peak_list = []

  def __getitem__(self, i):
    return self.peak_list[i]

  def residuals(self, params, E, spectrum):
    res = zeros(E.size)
    sum_up = zeros(E.size)
    param_eater = params
    for peak in self.peak_list:
      p = param_eater[:peak.params.size]
      param_eater = param_eater[peak.params.size:]
      sum_up += peak.f(p,E)

    res = spectrum-sum_up

    param_eater = params
    for peak in self.peak_list:
      i = 0
      p = param_eater[:peak.params.size]
      param_eater = param_eater[peak.params.size:]
      for param in p:
        res *= peak.penalties[i](param)
        i += 1

    return abs(res)

  def optimize_fit(self, E, spectrum):
    params = []
    for peak in self.peak_list:
      params = append(params, peak.params)

    plsq = leastsq(self.residuals, params, args=(E, spectrum))

    params = plsq[0]

    for peak in self.peak_list:
      peak.params = params[:peak.params.size]
      params = params[peak.params.size:]

    return self(E)

  def add_peak(self, peak):
    self.peak_list.append(peak)

  def __call__(self, E):
    sum_up = zeros(E.size)
    for peak in self.peak_list:
      sum_up += peak(E)
    return sum_up

