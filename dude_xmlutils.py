import xml.etree.ElementTree as et

class _XMLFile(object):
  def __init__(self,name,assign_ids=False):
    self.name = name
    if name is None:
      raise Exception('no xml file specified')
    self.tree = et.parse(self.name)
    self.root = self.tree.getroot()
    if assign_ids:
      self.assign_ids()
  def assign_ids(self,tag='Absorber'):
    "assigns an id for each absorber where id not present"
    attribute_list = ['id','ionName','z']
    counter=0
    for tag in self.root.findall('CompositeSpectrum'):
      for item in tag.findall(tag):
        if item.get('id')=="" or item.get('id')=="null":
          ionName = item.get('ionName')
          for attr in attribute_list:
            item.set('id',ionName+'%d'%counter)
            counter+=1
    self.tree.write(self.xmlfile)

  def getData(self,iden,attribute_list,tag=None):
    for tag in self.root.findall('CompositeSpectrum'):
      for item in tag.findall(tag):
        if item.get('id')==iden:
          return [item.get(attr) for attr in attribute_list ]            
    raise Exception('id '+iden+' not found')

  def writeData(self,iden,tag=None,**kwargs):
    flag=False
    for tag in self.root.findall('CompositeSpectrum'):
      for item in tag.findall(tag):
        if item.get('id')==iden:
          flag = True
          for key, val in kwargs.iteritems():
            try:
              item.set(key,val)
            except:
              print("item "+key+" not found...skipping\n")
    if flag==False:
      raise Exception('id '+iden+' not found')
    self.tree.write(self.name)

class Data(object):
  def __init__(self,**kwargs):
    for key, val in kwargs.iteritems():
      setattr(self,key,val)
    assign_ids = kwargs.get('assign_ids',False)
    self.tag  = kwargs.get('tag','Absorber')
    self.xmlfile = _XMLFile(kwargs.get('xmlfile',None),assign_ids)     
    self.tree = et.parse(xmlfile)
    self.root = self.tree.getroot()
    self.iden = kwargs.get('iden',None)
    if not self.iden:
      raise Exception('no id specified')

  def __str__(self):
    if self.tag=='ContinuumPoint':
      return "%s %12.7lf %12.8E"%(self.iden,self.x,self.y) 
    elif: self.tag=='Absorber'
      return "%s %9.6lf %9.6lf %10.8lf %8.4lf"%(self.iden,self.N,self.b,self.z,self.vel)
    else:
      raise Exception('undefined tag: '+self.tag)

  def getData(self):
    lst = ['N','b','z'] if self.tag=='Absorber' else ['x','y']
    output = zip(lst, self.xmlfile.getData(self.iden, lst, self.tag))
    for item in output:
      setattr(self,item[0],item[1])
    
  def writeData(self,**kwargs):
    self.xmlfile.writeData(self.iden,self.tag,**kwargs)

