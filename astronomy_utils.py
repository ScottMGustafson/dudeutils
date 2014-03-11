from math import sqrt
import warnings
import dude_xmlutils
c = 299792.458 # speed of light in km/s

class Absorber(dude_xmlutils.Data):
  def __init__(self, **kwargs):
    """
    accepted keywords:
    N, b, z, ion, iden
    """
    super(Absorber,self).Data(**kwargs)
    self.tag='Absorber'
    self.get_lines()

  def get_lines(self, fname='atom.dat'):
    """
    get all data for atom.dat and put into a list of dicts
    """
    linelst = []
    f=open(fname,'r')
    for line in f:
      line=line.split()
      ion, wave, f = line[0], float(line[1]), float(line[2])
      if ion==self.ion:
        linelst.append({'ion':ion,'wave':wave,'f':f})
    f.close()
    linelst.sort(key=lambda item:item['wave'] ,reverse=True)
    self.wave = [ item['wave'] for item in linelst ]
    self.f = [ item['f'] for item in linelst ]
    self.obs_wave = [ (1.+self.z)*item['wave'] for item in linelst ]
    
def get_vel_shift(z1,z2,delz2=None,delz1=None):
  """get vel shift between two redshifts"""
  if z1==0. or z2==0.:
    warnings.warn('warning: redshift is 0.')
  vel = c*(z1-z2)/(1.+z1)
  if delz1==None and delz2==None:
    return vel
  else:
    vel_unc = (c/(1.+z1))*sqrt(delz2**2.+(delz1**2.)*((1.-z2)/(1.+z2))**2.)
    return vel, vel_unc
    
def vel_limit(vel,ident,delta=None):
  if delta is None:
    delta = 0.
  if fabs(vel)>1.0+delta:
    warnings.warn("vel("+ident+") = "+str(vel))

def get_lines(fname='atom.dat'):
  """
  get all data for atom.dat and put into a list of dicts
  """
  linelst = []
  f=open(fname,'r')
  for line in f:
    line=line.split()
    ion, wave, f = line[0], float(line[1]), float(line[2])
    linelst.append({'ion':ion,'wave':wave,'f':f})
  return linelst

def get_series(lst=get_lines(),**kwargs):
  """
  parse out all data for a particular series of transitions for ion 'ion' 
  """
  ion = kwargs.get('ion','HI')
  retlst = []
  for item in lst:
    if item['ion']==ion:
      retlst.append(item)
  retlst = sorted(retlst, key=lambda item:item['wave'] ,reverse=True)
  if len(retlst)<chosen_transitions[-1]:
    raise Exception("chosen transitions not in data.")
  return retlst

def get_obs_wv(z,restwave):
  return restwave*(1.+z)

def z(obswv,restwv):
  return (obswv/restwv) -1.

def get_vel(wv1,wv2,restwv):
  """
  if a wv2 is redder than wv1, then vel should be positive

  gets velocity shift of wv2 to wv1

  input parameters:
  -----------------
  wv1, wv2 : the actual wavelength of the data
  restwv : rest wavelength of the desired center of the desired absorption line

  output:
  -------
  velocity shift to restwv 
  """

  assert(wv1>0. and restwv>0. and wv2>0.)
  z1=z(wv1,restwv)
  z2=z(wv2,restwv)
  return c*(z2-z1)/(1.+z1)
