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
        attribute_list = ['id','ionName']
        counter=0
        for thetag in self.root.findall('CompositeSpectrum'):
            for item in thetag.findall(tag):
                if item.get('id')=="" or item.get('id')=="null":
                    item.set('id',item.get('ionName')+'%d'%counter)
                    counter+=1
        self.tree.write(self.name)

    def getData(self,iden,tag,attribute_list=None):
        for thetag in self.root.findall('CompositeSpectrum'):
            for item in thetag.findall(tag):
                if item.get('id')==iden:
                    if attribute_list is None:
                        attribute_list= [dict(it.attrib).keys()]
                    return list(zip(attribute_list, [item.get(attr) for attr in attribute_list ]))
        raise Exception('id '+iden+' not found in '+self.name)

    def getViewData(self,iden,tag,attribute_list=None):
        """needs separate get function"""
        if tag == 'Region':
            identityTag='start'
        else:
            identityTag='id'

        for thetag in self.root.findall(tag):
            if thetag.get(identityTag)==iden:
                if attribute_list is None:
                    attribute_list= list(dict(thetag.attrib).keys())
                return list(zip(attribute_list, [thetag.get(attr) for attr in attribute_list ]))
        raise Exception('id '+iden+' not found in '+self.name)

    def writeData(self,iden,tag,**kwargs):
        flag=False
        for thetag in self.root.findall('CompositeSpectrum'):
            for item in thetag.findall(tag):
                if item.get('id')==iden:
                    flag = True
                    for key, val in kwargs.items():
                        try:
                            item.set(key,val)
                        except:
                            print(("item "+key+" not found...skipping\n"))
        if flag==False:
            raise Exception('id '+iden+' not found')
        self.tree.write(self.name)

class Data(object):
    def __init__(self,**kwargs):
        for key, val in list(kwargs.items()):
            setattr(self,key,val)
        assign_ids  = kwargs.get('assign_ids',False)
        self.tag    = kwargs.get('tag','Absorber')
        self.xmlfile = _XMLFile(kwargs.get('xmlfile',None),assign_ids)         
        self.tree = et.parse(self.xmlfile.name)
        self.root = self.tree.getroot()
        self.iden = kwargs.get('iden',None).strip()
        self.vel=0.
        if not self.iden:
            raise Exception('no id specified')

    def getData(self,lst=None,function='getData',*args):
        func   = getattr(self.xmlfile,function)
        output = func(self.iden, self.tag, lst)
        for item in output:
            try:
                setattr(self,str(item[0]),float(item[1]))
            except:
                setattr(self,str(item[0]),str(item[1]))

    def getViewData(self,lst=None,*args):
        output = self.xmlfile.getViewData(self.iden, self.tag, lst)
        for item in output:
            try:
                setattr(self,str(item[0]),float(item[1]))
            except:
                setattr(self,str(item[0]),str(item[1]))
        
    def writeData(self,**kwargs):
        self.xmlfile.writeData(self.iden,self.tag,**kwargs)

class ContinuumPoint(Data):
    def __init__(self,**kwargs):
        super(ContinuumPoint, self).__init__(tag="ContinuumPoint",**kwargs)
        if kwargs.get('populate',True) is True:
            self.getData()
    def __str__(self):
        return "%s %12.7lf %12.8E"%(self.iden,self.x,self.y)
    def getData(self):
        super(ContinuumPoint, self).getData(['x','y']) 

class Absorber(Data):
    def __init__(self,**kwargs):
        super(Absorber, self).__init__(tag="Absorber",**kwargs)
        if kwargs.get('populate',True) is True:
            self.getData()
    def __str__(self):
        return "%s %9.6lf %9.6lf %10.8lf %8.4lf"%(self.iden,self.N,self.b,self.z,self.vel)
    def getData(self):
        super(Absorber, self).getData(['N','b','z'])

class VelocityView(Data):
    def __init__(self,**kwargs):
        super(VelocityView, self).__init__(tag="VelocityView",**kwargs)
        if kwargs.get('populate',True) is True:
            self.getData()
    def getData(self):
        super(VelocityView, self).getData(function='getViewData')

class SingleView(Data):
    def __init__(self,**kwargs):
        super(SingleView, self).__init__(tag="SingleView",**kwargs)
        if kwargs.get('populate',True) is True:
            self.getData()
    def getData(self):
        super(SingleView, self).getData(function='getViewData')

class Region(Data):
    def __init__(self,**kwargs):
        super(Region, self).__init__(tag="Region",**kwargs)
        if kwargs.get('populate',True) is True:
            self.getData()
    def getData(self):
        super(Region, self).getData(function='getViewData')

