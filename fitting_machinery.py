import os
import re
from pylab import *
from scipy.optimize import leastsq
from fitting_functions import *

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

  def __init__(self, name, params, penalties, kernel, kernel_end):
    self.name = name
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
    self.res_call_count = 0

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

    if self.res_call_count % 200 == 0:
      print self.res_call_count
      figure()
      title(self.peak_list[0].name + ' - %f'%self.res_call_count )
      plot(E, res, 'k')
      plot(E, sum_up, 'b')
      plot(E, spectrum, 'r')

    param_eater = params
    for peak in self.peak_list:
      i = 0
      p = param_eater[:peak.params.size]
      param_eater = param_eater[peak.params.size:]
      for param in p:
        pen = peak.penalties[i](param)
        if self.res_call_count % 200 == 0:
          print 'param:   %f,   penalty: %f'% (param,pen)
        res *= pen
        i += 1

    self.res_call_count += 1

    return res

  def get_params(self):
    parm = []
    for peak in self.peak_list:
      parm = append(parm, peak.params)
    return parm

  def set_params(self, params):
    for peak in self.peak_list:
      peak.params = params[:peak.params.size]
      params = params[peak.params.size:]

  def optimize_fit(self, E, spectrum):
    params = self.get_params()

    p = leastsq(self.residuals, params, args=(E, spectrum))

    params = p[0]

    print 'Optimized parameters: %s' % (p,)

    self.set_params(params)

    return self(E)

  def add_peak(self, peak):
    self.peak_list.append(peak)

  def __call__(self, E):
    sum_up = zeros(E.size)
    for peak in self.peak_list:
      sum_up += peak(E)
    return sum_up

