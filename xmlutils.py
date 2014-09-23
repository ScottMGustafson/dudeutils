"""
some tools to interface with the standard dude xml input/output
"""
#import abc
import xml.etree.ElementTree as et
import datetime
from xml.dom import minidom
import os.path

"""
the end goal here is that I want a tool that can with one command grab all 
current vals from a dude xml file.  


while running dude:
    fit and save
    grab and append to db

"""

class Basexml(object):
    #__metaclass__ = abc.ABCMeta
    @staticmethod
    def get_root(filename):
        if not os.path.exists(filename):
            return None
        tree = et.parse(filename)
        return tree.getroot()

    @staticmethod
    def read(filename):
        tree = et.parse(filename)
        return tree.getroot()
       
    #@abc.abstractmethod
    def write(self):
        """write the current data"""
        return

    @staticmethod
    def get_children(parent,child_tag):
        return parent.findall(child_tag)

    @staticmethod
    def get_node_data(node,attribute_list=None):
        if attribute_list is None:
            attribute_list= list(dict(node.attrib).keys())
        return dict(zip(attribute_list, [node.get(attr) for attr in attribute_list ]))
            
    @staticmethod
    def get_node(parent,tag,iden):
        for item in parent.findall(tag):
            if item.get('id')==iden:
                return item
        raise Exception('id '+iden+' not found')

    @staticmethod
    def get_node_list(parent,tag):
        return parent.findall(tag)

    @staticmethod
    def get_children(lst):
        elist = []
        for item in lst:
            attribs = item.get_keys()
            vals = [ str(getattr(item,it)) for it in attribs ]
            
            data= dict(zip(attribs, vals))

            elist.append( Element(item.tag, attrib=data) )
        return elist
        

class Model_xml(Basexml):
    """xml format for storing model data"""
    def __init__(self, db=None, **kwargs):
        """ 
        optional keywords:
        filename
        db
        """

        self.filename = kwargs.get("filename", Model_xml.default_filename() ) 
        self.db = kwargs.get("db", None)

    @classmethod
    def default_filename(cls):
        now = str(datetime.datetime.now())
        return now[0:10]+"model.xml"
        
    @staticmethod
    def create(filename,db):
        """create an xml file file structure.  returns root"""
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
            group_name = item.id 
            if current_group is None or group_name != current_group.text:
                current_group = et.SubElement(root, 'model', {'id':group_name})

            children = []
            if item.absorbers!=None:
                children += self.get_children(item.absorbers) 
            if item.continuum_points!=None:
                children += self.get_children(item.continuum_points) 
            if item.regions!=None:
                children += self.get_children(item.regions)  

            if len(children)==0:
                raise Exception("no children are present")
            current_group.extend(children)
        return root

    @staticmethod
    def write(filename, root):
        f = open(filename,'w')
        f.write(prettify(root))
        f.close()
        return



class Dudexml(Basexml):
    """
    the standard dude xml fit file
    """
    def __init__(self,name,assign_ids=False):
        self.name = name
        if name is None:
            raise Exception('no xml file specified')
        self.root = super(Dudexml,self).get_root(name)
        if assign_ids:
            self.assign_ids()

    def get_spectrum(self):
        return self.root.findall('CompositeSpectrum')[0]

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

    def get_node_data(self,**kwargs):
        """returns dict of data for node.  
        input:
        ------
        iden:  id to look for
        tag:   tag of xml node
        node:  if iden and tag not specified, may specify node instead.
        attribute_list:  list of attributes to return.  by default, returns all

        """
        iden=kwargs.get("iden",None)
        tag=kwargs.get("tag",None)
        node=kwargs.get("node",None)
        attribute_list=kwargs.get("attribute_list",None)

        if node == None:
            node = self.get_node(iden,tag)  

        if attribute_list!=None:
            args = [node, attribute_list] 
        else:
            args = [node]
        return super(Dudexml,self).get_node_data(*args)
  
    def get_node(self,iden,tag):
        if iden is None:
            raise Exception('need to define an iden')
        return super(Dudexml,self).get_node(self.root.findall('CompositeSpectrum'),tag,iden)

    def get_node_list(self,tag,parent=None):
        if tag is None:
            raise Exception("tag must be defined")
        if tag in ["ContinuumPoint","Absorber"]:
            return super(Dudexml,self).get_node_list(self.root.findall("CompositeSpectrum")[0],tag)
        elif parent:
            return super(Dudexml,self).get_node_list(parent,tag)
        else:
            return super(Dudexml,self).get_node_list(self.root,tag)

    def write(self,**kwargs):
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


        self.node = kwargs.get('xmlnode',None)  #to explicitly instantiate one node
        for key, val in list(kwargs.items()):
            setattr(self,key,val)
        
        if self.node != None:
            for key, val in dict(self.node.attrib).items():
                setattr(self,key,val)
        assign_ids  = kwargs.get('assign_ids',False)
        self.tag    = kwargs.get('tag')
        self.xmlfile = Dudexml(kwargs.get('xmlfile'))       
        self.id = kwargs.get('iden',None)

        if type(self.id) is str:
            self.node = self.xmlfile.get_node(self.id,self.tag)
        elif self.id is None and not self.node is None:
            self.parseNode(self.node)


#is there going to be an issue with id versus iden?
    def parseNode(self,node):
        data = node.attrib
        print(data)
        for key, val in data.items():
            try:
                setattr(self, str(key), float(val))
            except:
                setattr(self, str(key), str(val))
        
    def getData(self,**kwargs):
        if self.node==None:
            self.node=kwargs.get("node",None)
        if self.node!=None:
            dat=self.xmlfile.get_node_data(node=self.node)
        else:
            tag=kwargs.get("tag",self.tag)
            iden=kwargs.get("iden",self.id)
            dat = self.xmlfile.get_node_data(iden=iden, tag=tag)
        for key, val in dat.items():
            try:
                setattr(self,str(key),float(val))
            except:
                setattr(self,str(key),str(val))

    def write(self):
        self.writeNode()
        self.xmlfile.write()

    def writeNode(self):
        if self.node != None:
            for key,val in dict(kwargs).items():
                setattr(self,key,val)
                self.node.set(key,val)
        else:
            keys = self.xmlfile.get_keys(self.tag)  #get relevant keys to write
            vals = [getattr(self,key) for key in keys] #get this instances vals for them
            data = dict(zip(keys, vals))
            et.extend(self.xmlfile.get_spectrum(), self.tag, data)

    def get_keys(self):
        return self.xmlfile.get_keys(self.tag)

class ContinuumPoint(Data):
    def __init__(self,**kwargs):
        super(ContinuumPoint, self).__init__(tag="ContinuumPoint",**kwargs)
        self.getData()
    def __str__(self):
        return "%s %12.7lf %12.8E"%(self.id,self.x,self.y)
        
class Absorber(Data):
    def __init__(self,**kwargs):
        super(Absorber, self).__init__(tag="Absorber",**kwargs)
        if self.xmlfile is None and self.node is None:
            raise Exception("no xml fit file associated with this absorber: %s"%(self.id))
        if kwargs.get('populate',True) is True:
            self.getData()
        self.get_lines()

    def __str__(self):
        return "iden=%6s N=%8.5lf b=%8.5lf z=%10.8lf"%(self.id,self.N,self.b,self.z)

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
        for item in ['N','b','z','iden']:
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

class SingleView(Data):
    def __init__(self,**kwargs):
        super(SingleView, self).__init__(tag="SingleView",**kwargs)
        if kwargs.get('populate',True) is True:
            self.getData()

class Region(Data):
    def __init__(self,**kwargs):
        super(Region, self).__init__(tag="Region",**kwargs)
        if kwargs.get('populate',True) is True:
            self.getData()

def consolidate_regions(lst):
    pass
    #TODO

def prettify(elem):
    reparsed = minidom.parseString(et.tostring(elem, 'utf-8'))
    return reparsed.toprettyxml(indent="  ")


