from math import sqrt
import warnings
c = 299792.458 # speed of light in km/s

def get_vel_shift(z1,z2,delz2=None,delz1=None):
  if z1==0. or z2==0.:
    warnings.warn('warning: redshift is 0.')
  vel = c*(z1-z2)/(1.+z1)
  if delz1==None and delz2==None:
    return [vel]
  else:
    vel_unc = (c/(1.+z1))*sqrt(delz2**2.+(delz1**2.)*((1.-z2)/(1.+z2))**2.)
    return [vel,vel_unc]
    
def vel_limit(vel,ident,delta=None):
  if delta is None:
    delta = 0.
  if fabs(vel)>1.0+delta:
    warnings.warn("vel("+ident+") = "+str(vel))
