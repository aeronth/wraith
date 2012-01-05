from pylab import *

###########################################
# Functions to feed into fitting routines 
###########################################

def analyzer_function(p, E):
    A, B = p
    BG = A/sqrt(abs(B-E))
    #print "A: %f, B: %f"%(A,B)
    return BG

def slope(p, E):
    m, b = p
    BG = b + m * E
    return BG

def fitting_arctan(p, E):
    A, B, u, b = p
    BG = b + A * arctan( B * (E - u) ) 
    return BG

def sloped_arctan(p, E):
    A, B, u, m, b = p
    BG = b + m * E + A * arctan( B * (E - u) ) 
    return BG
sloped_arctan.latex = r'$s_arctan(E) = b + mE + A {\rm arctan} ( B (E - u) )$'

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
    return gl_(a, mu, sigma, m, E) + gl_(ratio_a * a, mu - split, ratio_sigma * sigma, m, E)
spin_split_gl.latex = r'$ssgl(E) = gl(A,\mu,\sigma,m,E) + gl(R_A A, \mu + \Delta_E, \frac{R_{Area}}{R_A} \sigma, E)$'

def gl_(a, mu, sigma, m, E):
    return a * exp(-2.772589 * (1 - m) * (E - mu)**2/sigma**2) / (1 + 4 * m * (E - mu)**2/sigma**2)

def gl(params, E):
    a, mu, sigma, m = params
    return gl_(a, mu, sigma, m, E)
gl.latex = r'$gl(A,\mu,\sigma,m,E) = A e^{\frac{(-4ln(2) (1-m) (E - \mu)^2/\sigma^2)}{(1+4m(E-\mu)^2/\sigma^2)}}$'

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

#####################################################
# Penalty functions
#
# Returns a scale factor for residuals depending on
# the relationship between the parameter value and
# a declared range
#####################################################

#Residuals are unscaled regardless of range
def no_penalty(range, p):
    return 1.0

#Residuals are scaled by a constant factor of 100 if parameter value is outside of the range
def notch_penalty(range, p):
    value = 1.0
    if p < range[0] or p > range[1]:
        value = 100.0
    return value

#Residuals are scaled by an exponentially growing factor if parameter is outside of the range
def exp_penalty(range, p):
    A = 1.0/diff(range)
    value = 1.0
    if p < range[0]:
      value = 1.0 + exp(A*((range[0]+0.05*diff(range)) - p ))
    elif p > range[1]:
      value = 1.0 + exp(A*(p - (range[1]-0.05*diff(range)) ))
    return value

#Residuals are scaled by a quadratically growing factor if parameter value is outside of the range
def quad_penalty(range, p):
    A = 100.0/diff(range)
    value = 1.0
    lower_bound = (range[0]+0.02*diff(range))
    upper_bound = (range[1]-0.02*diff(range))
    if p < lower_bound:
        value = 1.0 + A*(lower_bound - p )**2
    elif p > upper_bound:
        value = 1.0 + A*(p - upper_bound )**2
    return value
    
