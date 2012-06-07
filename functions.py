from pylab import *
import scipy.linalg
from scipy.misc import factorial

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


