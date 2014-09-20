"""
some tools to interface with the standard dude xml input/output
"""
import xml.etree.ElementTree as et
import datetime
from data_structures import ModelDB, read_in
from xml.dom import minidom

class XML_db(object):
    def __init__(self, db=None, filename=None):
        now = str(datetime.datetime.now())

        if filename is None:
            self.filename=now[0:9]+".xml"
        else:
            self.filename = filename
        self.db = db

        if db != None:
            self.root = self.load_from_db(self.db)
        else:
            self.root = None


    def load_from_db(self,db):
        filename = self.filename
        
        now = str(datetime.datetime.now())
        root = et.Element('modeldb')

        #set up header data
        head = et.SubElement(root, 'head')
        title = et.SubElement(head, 'title')
        title.text = 'Fitting Models'
        created = et.SubElement(head, 'dateCreated')
        created.text = now
        modified = et.SubElement(head, 'dateModified')
        modified.text = now

        models = et.SubElement(root, 'models')

        #load the model db
        for item in db.lst:
            current_group = None
            group_name = item.iden 
            if current_group is None or group_name != current_group.text:
                current_group = et.SubElement(root, 'model', {'id':group_name})

            children = []
            if item.absorbers!=None:
                children += get_children(item.absorbers) 
            if item.continuum_points!=None:
                children += get_children(item.continuum_points) 
            if item.regions!=None:
                children += get_children(item.regions)  

            if len(children)==0:
                raise Exception("no children are present")
            current_group.extend(children)
        return root

    def write(self):
        f = open(self.filename,'w')
        f.write(prettify(self.root))
        f.close()
        return

    def read(self):
        """read from xml, return ModelDB instance"""
        pass
        

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
        node = self.findNode(iden,tag)
        if attribute_list is None:
            attribute_list= list(dict(node.attrib).keys())
        return list(zip(attribute_list, [node.get(attr) for attr in attribute_list ]))

    def findNode(self,iden,tag):
        if iden is None:
            raise Exception('need to define an iden')
        for thetag in self.root.findall('CompositeSpectrum'):
            for item in thetag.findall(tag):
                if item.get('id')==iden:
                    return item
        raise Exception('id '+iden+' not found')

    def getNodeList(self,tag):
        if tag is None:
            raise Exception("tag must be defined")
        if tag=="Region":
            return list(self.root.findall(tag))
        else:
            thetag = self.root.findall('CompositeSpectrum')
            return [item for item in thetag.findall(tag)]
        
    def setVal(self,node,key,newval):
        node.set(key,newval)

    def writeOut(self):
        self.tree.write(self.name)

    def getDataList(self,tag):
        """return a list of all instances of tag"""
        temproot = self.root.findall('CompositeSpectrum')[0]
        lst = temproot.findall(tag)
        return [it for it in lst]

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
        node = self.findNode(iden,tag)
        for key, val in dict(kwargs).items():
            node.set(key,val)
        self.tree.write(self.name)

    def get_keys(self, tag):
        sp = self.root.findall('CompositeSpectrum')[0]
        children = sp.findall(tag)
        if children==[]:
            raise Exception(str(children)+"does not contain "+str(tag))
        return list(dict(children[0].attrib).keys())
        

class Data(object):
    def __init__(self,**kwargs):
        """
        mandatory kwargs:  
        ---------------
        xmlfile : (str) name of xmlfile
        tag:  (str) name of child class
        xmlnode: a node that xml.etree parsed from an xmlfile
        """
        for key, val in list(kwargs.items()):
            setattr(self,key,val)
        self.node = kwargs.get('xmlnode',None)  #to explicitly instantiate one node
        if not self.node is None:
            for key, val in dict(self.node.attrib).items():
                setattr(self,key,val)
        assign_ids  = kwargs.get('assign_ids',False)
        self.tag    = kwargs.get('tag')
        self.xmlfile = _XMLFile(kwargs.get('xmlfile'))       
        self.iden = kwargs.get('iden',None)

        if type(self.iden) is str:
            self.node = self.xmlfile.findNode(self.iden,self.tag)

    def getData(self,lst=None,function='getData',**kwargs):
        if function=='getDataList':
            return self.xmlfile.getDataList(self.tag,**kwargs)
        func   = getattr(self.xmlfile,function)
        output = func(self.iden, self.tag, lst)
        for item in output:
            try:
                setattr(self,str(item[0]),float(item[1]))
            except:
                setattr(self,str(item[0]),str(item[1]))
        
    def writeData(self,**kwargs):
        self.xmlfile.writeData(self.iden,self.tag,**kwargs)
    def writeOut(self):
        self.tree.write()
    def writeData(self, **kwargs):
        if not self.node is None:
            for key,val in dict(kwargs).items():
                setattr(self,key,val)
                self.node.set(key,val)
        else:
            #behavior to add new child node
            pass
        self.xmlfile.writeOut()

    def get_keys(self):
        return self.xmlfile.get_keys(self.tag)

class ContinuumPoint(Data):
    def __init__(self,**kwargs):
        super(ContinuumPoint, self).__init__(tag="ContinuumPoint",**kwargs)
        self.getData()
    def __str__(self):
        return "%s %12.7lf %12.8E"%(self.iden,self.x,self.y)
    def getData(self):
        super(ContinuumPoint, self).getData(['x','y']) 
    def write(self, **kwargs):
        super(ContinuumPoint, self).writeData(kwargs) 
        
class Absorber(Data):
    def __init__(self,**kwargs):
        super(Absorber, self).__init__(tag="Absorber",**kwargs)
        if not self.xmlfile:
            raise Exception("no xml fit file associated with this absorber: %s"%(self.iden))
        if kwargs.get('populate',True) is True:
            self.getData()
        self.get_lines()

    def __str__(self):
        return "iden=%6s N=%8.5lf b=%8.5lf z=%10.8lf"%(self.iden,self.N,self.b,self.z)
    
    def getData(self):
        super(Absorber, self).getData()
        self.ionName = self.ionName.replace(' ','')

    def locked(self,param):
        param_lock = {'N':'NLocked', 'b':'bLocked', 'z':'zLocked'}
        tf_lst = ['true', 'True', 'TRUE']
        try:
            ans = getattr(self,param_lock[param])
            if str(ans) not in tf_lst:
                return False
            else:
                return True
        except KeyError:
            raise Exception("no param named %s"%{param})
            
    def get_lines(self, filename='atom.dat'):
        """
        get all data for atom.dat and put into a list of dicts
        """
        linelst = []
        fname=open(filename,'r')
        for line in fname:
            line=line.split()
            ion, wave, f = line[0], float(line[1]), float(line[2])
            if ion==self.ionName:
                linelst.append({'ion':ion,'wave':wave,'f':f})
        fname.close()
        linelst.sort(key=lambda item:item['wave'] ,reverse=True)
        self.wave = [ item['wave'] for item in linelst ]
        self.f = [ item['f'] for item in linelst ]
        self.obs_wave = [ (1.+self.z)*item['wave'] for item in linelst ]
    def write(self, **kwargs):
        super(Absorber, self).writeData(kwargs) 
    def getShift(self, ref):
        """ 
        get velocity shift in km/s
         Inputs:
        --------
        ref: a reference Absorber instance

        returns:
        --------
        velocity (in km/s)
        """
        c=299792.458
        self.vel = (self.z-ref.z)*c/(1.+ref.z)
        return self.vel

    def __eq__(self,other):
        for item in ['N','b','z','iden','xmlfile']:
            if getattr(self,item)!=getattr(other,item):
                return False
        return True

    def __neq__(self,other):
        return not self.__eq__(other)
        
class VelocityView(Data):
    def __init__(self,**kwargs):
        super(VelocityView, self).__init__(tag="VelocityView",**kwargs)
        if kwargs.get('populate',True) is True:
            self.getData()
    def getData(self):
        super(VelocityView, self).getData(function='getViewData')
    def write(self, **kwargs):
        super(Absorber, self).writeData(kwargs) 

class SingleView(Data):
    def __init__(self,**kwargs):
        super(SingleView, self).__init__(tag="SingleView",**kwargs)
        if kwargs.get('populate',True) is True:
            self.getData()
    def getData(self):
        super(SingleView, self).getData(function='getViewData')
    def write(self, **kwargs):
        super(Absorber, self).writeData(kwargs) 

class Region(Data):
    def __init__(self,**kwargs):
        super(Region, self).__init__(tag="Region",**kwrgs)
        if kwargs.get('populate',True) is True:
            self.getData()
    def getData(self):
        """ gets and returns list of Region instances from xml src file"""
        super(Region, self).getData(function='getViewData')
    def write(self, **kwargs):
        super(Absorber, self).writeData(kwargs) 

def getContinuumPoints(xmlfile):
    xml = _XMLFile(xmlfile)
    conts = xml.getDataList('ContinuumPoint')
    return [ ContinuumPoint(xmlnode=item) for item in conts ]

def getList(xmlfile,classname):
    """ more general than the above.  """
    return _XMLFile(xmlfile).getDataList(name)
    

def consolidate_regions(lst):
    pass
    #TODO



def prettify(elem):
    reparsed = minidom.parseString(et.tostring(elem, 'utf-8'))
    return reparsed.toprettyxml(indent="  ")


def get_children(lst):
    elist = []
    for item in lst:
        attribs = item.get_keys()
        vals = [ str(getattr(item,it)) for it in attribs ]
        
        data= dict(zip(attribs, vals))

        elist.append( Element(item.tag, attrib=data) )
    return elist



    

