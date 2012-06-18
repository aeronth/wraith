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

  def __init__(self, range, f):
    """Initialize penalty function"""
    self.range = range
    self.f = f

  def __call__(self, p):
    """Penalty!"""
    return self.f(self.range, p)

class AnalyticBackground:

  def __init__(self, spectrum, name, variables, values, penalties, f):
    self.spectrum = spectrum
    self.name = name
    self.variables = variables
    self.values = values
    self.penalties = penalties
    self.f = f

  def get_spec(self):
    ranges = []
    for pen in self.penalties:
      ranges.append(pen.range)

    spec = { 'name': self.name,
                'function': self.f.func_name,
                'penalty_function': self.penalties[0].f.func_name,
                'variables': self.variables,
                'values': self.values,
                'ranges': ranges
                }
    return spec

  def residuals(self, values, E, spectrum):
    res = spectrum - self.f(values, E)
    i = 0
    for p in values:
      res *= self.penalties[i](p)
      i += 1

    res[res<0] = res[res<0]*2
    return res

  def optimize_fit(self, E, spectrum):
    p = leastsq(self.residuals, self.values, args=(E, spectrum))

    self.values = p[0]
    print self.values

    return self(E)

  def __call__(self, E):
    return self.f(self.values, E)

   
class Background:
  """Implements a generic background based on a convolution kernel function"""

  def __init__(self, spectrum, name='tougaard', variables=['B','D'], values=r_[100,100], penalties=[Penalty(r_[0,100],no_penalty), Penalty(r_[0,100],no_penalty)], kernel=K, kernel_end=200):
    self.spectrum = spectrum
    self.name = name
    self.values = values
    self.variables = variables
    self.penalties = penalties
    self.kernel = kernel
    self.kernel_end = kernel_end
    self.optimization_history = ones((size(variables)+1,1))

  def set_spec(self, spec):
    self.penalties = []
    for range in spec['ranges']:
      self.penalties.append(Penalty(range, eval(spec['penalty_function'])))
    self.name = spec['name']
    self.variables = spec['variables']
    self.values = spec['values']
    self.kernel = eval(spec['function'])
    self.optimization_history = c_[self.optimization_history,r_[self.values,sum(abs(self.spectrum.residuals()))]]

  def get_spec(self):
    ranges = []
    for pen in self.penalties:
      ranges.append(pen.range)

    bg_spec = { 'name': self.name,
                'function': self.kernel.func_name,
                'penalty_function': self.penalties[0].f.func_name,
                'variables': self.variables,
                'values': self.values,
                'ranges': ranges
                }
    return bg_spec

  def f(self, values, E, spectrum):
    dE = E[1] - E[0]
    spectrum = spectrum - min(spectrum)
    bg = dE * convolve( spectrum, self.kernel( values, E)[::-1], 'full')
    return bg[bg.size-spectrum.size:]
  
  def residuals(self, values, E, spectrum):
    res = spectrum - self.f(values, E, spectrum)
    i = 0
    for p in values:
      res *= self.penalties[i](p)
      i += 1

    res[res<0] = res[res<0]*20
    self.optimization_history = c_[self.optimization_history,r_[values,sum(abs(res))]]
    return res
  
  def EE(self, dE):
    return arange(0, self.kernel_end, abs(dE))

  def optimize_fit(self, E, spectrum):
    offset = min(spectrum)
    spectrum = spectrum - offset
    self.dE = E[1]-E[0]
    plsq = leastsq(self.residuals, self.values, args=(self.EE(self.dE), spectrum))
    self.values = plsq[0]
    return self.f(self.values, self.EE(self.dE), spectrum) + offset
 
  def __call__(self, E, spectrum):
    offset = min(spectrum)
    spectrum = spectrum - offset
    self.dE = E[1]-E[0]
    return self.f(self.values, self.EE(self.dE), spectrum) + offset

class Peak:

  def __init__(self, spectrum, name='default', variables=['A','\mu','\sigma','m'], values=r_[100,100,1,0.5], penalties=[Penalty(r_[0,100],no_penalty)], f=gl):
    self.spectrum = spectrum
    self.name = name
    self.variables = variables
    self.values = values
    self.penalties = penalties
    self.f = f
    self.optimization_history = ones((size(variables)+1,1))
   
  def set_spec(self, spec):
    self.penalties = []
    for range in spec['ranges']:
      self.penalties.append(Penalty(range, eval(spec['penalty_function'])))
    self.name = spec['name']
    self.variables = spec['variables']
    self.values = spec['values']
    self.f = eval(spec['function'])
    self.optimization_history = c_[self.optimization_history,r_[self.values,sum(abs(self.spectrum.residuals()))]]

  def get_spec(self):
    ranges = []
    for pen in self.penalties:
      ranges.append(pen.range)

    spec = { 'name': self.name,
             'function': self.f.func_name,
             'penalty_function': self.penalties[0].f.func_name,
             'variables': self.variables,
             'values': self.values,
             'ranges': ranges
             }
    return spec

  def residuals(self, values, E, spectrum):
    res = spectrum - self.f(values, E)
    i = 0
    for p in values:
      res *= self.penalties[i](p)
      i += 1

    return res

  def __call__(self, E):
    return self.f(self.values, E)


class Peaks:

  def __init__(self, spectrum):
    """something will go here"""
    self.spectrum = spectrum
    self.peak_list = []

  def __getitem__(self, i):
    return self.peak_list[i]

  def set_spec(self,spec):
    for key,val in spec.items():
      self.peak_list.append(Peak())
      self.peak_list[-1].set_spec(val)

  def get_spec(self):
    spec = {}
    for peak in self.peak_list:
      spec['peak(%0.1f,%0.1f,...)'%(peak.values[0],peak.values[1])] = peak.get_spec()

    return spec

  def residuals(self, values, E, spectrum):

    res = zeros(E.size)
    sum_up = zeros(E.size)
    param_eater = values
    for peak in self.peak_list:
      p = param_eater[:peak.values.size]
      param_eater = param_eater[peak.values.size:]
      sum_up += peak.f(p,E)

    res = spectrum-sum_up

    param_eater = values
    for peak in self.peak_list:
      i = 0
      p = param_eater[:peak.values.size]
      param_eater = param_eater[peak.values.size:]
      peak.optimization_history = c_[peak.optimization_history,r_[p,sum(abs(res))]]
      for value in p:
        pen = peak.penalties[i](value)
        res *= pen
        i += 1

    return res

  def get_values(self):
    parm = []
    for peak in self.peak_list:
      parm = append(parm, peak.values)
    return parm

  def set_values(self, values):
    for peak in self.peak_list:
      peak.values = values[:peak.values.size]
      values = values[peak.values.size:]

  def optimize_fit(self, E, spectrum):
    values = self.get_values()

    p = leastsq(self.residuals, values, args=(E, spectrum))

    values = p[0]

    print 'Optimized parameters: %s' % (p,)

    self.set_values(values)

    return self(E)

  def add_peak(self, peak):
    self.peak_list.append(peak)

  def __call__(self, E):
    sum_up = zeros(E.size)
    for peak in self.peak_list:
      sum_up += peak(E)
    return sum_up

