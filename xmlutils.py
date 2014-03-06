import sys
import xml.etree.ElementTree as et
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

def parser(arg_list):
  datalst = []
  i=0
  action = 'write' if 'write' in arg_list else 'get'
  while i < len(arg_list):
    if "=" in arg_list[i]:
      key, val = arg_list[i].split('=')
      if key=='id':
        kwargs = {}
        kwargs['id']=val
        i+=1
        while arg_list[i][:2]!='id' and arg_list[i][:5]!='chisq':
          key, val = arg_list[i].split('=')
          kwargs[key] = val
          i+=1
        datalst.append(Data(**kwargs))
      elif key=='chisq':
        chisq=int(arg_list[i].split('=')[1])
    if '.xml' in arg_list[i]:
      xml_file = arg_list[i]
    i+=1

  return action, xml_file, datalst, chisq

action, xml_file, datalst, chisq = parser(sys.argv[1:])
tree = et.parse(xml_file) #parse xml file
root = tree.getroot()

for i in range(1,len(datalst)):
  setattr(datalst[i],'vel',getVelocityShift(datalst[0].z,datalst[i].z)[0])
  
if action=='get': 
  for item in datalst:
    item.getData(root)
  output_string = ''
  for item in datalst:
    output_string+=str(item)
  output_string+=' chisq='+str(chisq)
  f = open('absorberData.txt','a') # append read data to file
  f.write(output_string+'\n')
  f.close()

elif action=='write':
  for item in datalst:
    item.writeData(root)
  tree.write(xml_file)

else:
  raise Exception('unrecognized command: '+str(action))





