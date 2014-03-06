
from math import sqrt
c = 299792.458 # speed of light in km/s

class Data(object):
  def __init__(self,**kwargs):
    for key, val in kwargs.iteritems():
      setattr(self,key,val)
    self.tag = 'ContinuumPoint' if 'x' in kwargs.keys() else 'Absorber'

  def __str__(self):
    if self.tag=='ContinuumPoint':
      return "%s %12.7lf %12.8E"%(self.id,self.x,self.y) 
    elif: self.tag=='Absorber'
      return "%s %9.6lf %9.6lf %10.8lf %8.4lf"%(self.id,self.N,self.b,self.z,self.vel)
    else:
      raise Exception('undefined tag: '+self.tag)

  def getData(self,root):
    attribute_list = ['N','b','z'] if self.tag=='Absorber' else ['x','y']
    flag = False
    for tag in root.findall('CompositeSpectrum'):
      for item in tag.findall(self.tag):
        if item.get('id')==self.id:
          flag = True
          for attr in attribute_list:
            item.set(attr,str(self.attr))
    if flag==False:
      raise Exception('id '+self.id+' not found')

  def writeData(self,root):
    attribute_list = ['N','b','z'] if self.tag=='Absorber' else ['x','y']
    flag = False
    for tag in root.findall('CompositeSpectrum'):
      for item in tag.findall(self.tag):
        if item.get('id')==self.id:
          flag = True
          for attr in attribute_list:
            self.attr=item.get(attr)
    if flag==False:
      raise Exception('id '+self.id+' not found')

def getVelocityShift(z1,z2,delz2=None,delz1=None):
  if z1==0. or z2==0.:
    print 'warning: redshift is 0.'
  vel = c*(z1-z2)/(1.+z1)
  if delz1==None and delz2==None:
    return [vel]
  else:
    delVel = (c/(1.+z1))*sqrt(delz2**2.+(delz1**2.)*((1.-z2)/(1.+z2))**2.)
    return [vel,delVel]
    
def vel_limit(vel,ident,delta=None):
  if delta is None:
    delta = 0.
  if fabs(vel)>1.0+delta:
    warnings.warn("vel("+ident+") = "+str(vel))






