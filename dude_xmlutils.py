
class Data(object):
  def __init__(self,**kwargs):
    for key, val in kwargs.iteritems():
      setattr(self,key,val)
    #self.tag = 'ContinuumPoint' if 'x' in kwargs.keys() else 'Absorber'

  def __str__(self):
    if self.tag=='ContinuumPoint':
      return "%s %12.7lf %12.8E"%(self.iden,self.x,self.y) 
    elif: self.tag=='Absorber'
      return "%s %9.6lf %9.6lf %10.8lf %8.4lf"%(self.iden,self.N,self.b,self.z,self.vel)
    else:
      raise Exception('undefined tag: '+self.tag)

  def getData(self,root):
    attribute_list = ['N','b','z'] if self.tag=='Absorber' else ['x','y']
    flag = False
    for tag in root.findall('CompositeSpectrum'):
      for item in tag.findall(self.tag):
        if item.get('id')==self.iden:
          flag = True
          for attr in attribute_list:
            item.set(attr,str(self.attr))
    if flag==False:
      raise Exception('id '+self.iden+' not found')

  def writeData(self,root,attribute_list):
    #attribute_list = ['N','b','z'] if self.tag=='Absorber' else ['x','y']
    
    flag = False
    for tag in root.findall('CompositeSpectrum'):
      for item in tag.findall(self.tag):
        if item.get('id')==self.iden:
          flag = True
          for attr in attribute_list:
            self.attr=item.get(attr)
    if flag==False:
      raise Exception('id '+self.iden+' not found')



