from pylab import *
import scipy.linalg
from scipy.misc import factorial
from scipy import signal

def savgol_smooth(A, points, degree):
  order = 4
  if points<5:
    order = 2
  kernel = savgol(points,points/2,points/2,degree,order)#-ker[points]
  pad_start = array([])
  pad_end = array([])
  for i in range(0,(points/2)-1):
    pad_start = r_[A[0], pad_start]
    pad_end = r_[pad_end, A[-1]]
  sg = r_[pad_start, A, pad_end]
  sg = signal.fftconvolve(sg, kernel, 'valid')
  return sg

def savgol(np, nl, nr, ld, m):
  a = zeros([m+1, m+1])
  b = zeros(m+1)
  b[ld] = 1
  sg = zeros(np)
  for ipj in range(0,2*m + 1):
    s = 0.
    if ipj==0:
      s = 1
    for k in range(1,nr + 1):
      s = s + (1.*k)**ipj

    for k in range(1,nl + 1):
      s = s + (-1.*k)**ipj

    mm = min(ipj, 2*m-ipj)

    for imj in range(-mm,mm+1,2):
      a[(ipj+imj)/2, (ipj-imj)/2] = s

  #print 'a = %s\nb = %s'%(a,b)

  c = scipy.linalg.solve(a,b)

  #print 'c = %s'%(c)

  for k in range(-nl,nr+1):
    s = c[0]
    fac = 1
    for mm in range(1,m+1):
      fac = fac*k
      s = s + c[mm]*fac
    kk = mod(k+np/2,np)
    sg[kk] = s

  #print 'sg = %s'%sg

  sg = sg * factorial(ld)

  return sg


