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
    def __init__(self):
        pass

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
       
    def write(self,filename):
        self.tree.write(filename)

    @staticmethod
    def get_children(parent,child_tag):
        return parent.findall(child_tag)

    @staticmethod
    def get_node_data(node,attribute_list=None):
        if attribute_list is None:
            attribute_list= list(dict(node.attrib).keys())
        return dict(zip(attribute_list, [node.get(attr) for attr in attribute_list ]))
            
    @staticmethod
    def get_node(parent,tag,id):
        results=[item for item in parent.findall(tag) if item.get('id')==id]
        if len(results)>1:
            warnings.warn("more than one element contains the id \'%s\'\n returning the first one as default"%id)
            return results[0]
        elif len(results)==1:
            return results[0]
        else:
            raise Exception('id '+id+' not found')

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
        self.db = kwargs.get("db")

    @classmethod
    def default_filename(cls):
        now = str(datetime.datetime.now())
        return now[0:10]+"model.xml"

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
        self.tree = et.parse(name)
        self.root = self.tree.getroot()
        
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
        id:  id to look for
        tag:   tag of xml node
        node:  if id and tag not specified, may specify node instead.
        attribute_list:  list of attributes to return.  by default, returns all

        """
        id=kwargs.get("id",None)
        tag=kwargs.get("tag",None)
        node=kwargs.get("node",None)
        attribute_list=kwargs.get("attribute_list",None)

        if node == None:
            node = self.get_node(id,tag)  

        if attribute_list!=None:
            args = [node, attribute_list] 
        else:
            args = [node]
        return super(Dudexml,self).get_node_data(*args)
  
    def get_node(self,id,tag):
        if id is None:
            raise Exception('need to define an id')
        return super(Dudexml,self).get_node(self.root.findall('CompositeSpectrum')[0],tag,id)

    def get_node_list(self,tag,parent=None):
        if tag is None:
            raise Exception("tag must be defined")
        if tag in ["ContinuumPoint","Absorber"]:
            return super(Dudexml,self).get_node_list(self.root.findall("CompositeSpectrum")[0],tag)
        elif parent:
            return super(Dudexml,self).get_node_list(parent,tag)
        else:
            return super(Dudexml,self).get_node_list(self.root,tag)

    def set_node(self, id, tag,**kwargs):
        node=self.get_node(id,tag)
        for key, val in kwargs.items():
            node.set(key,str(val))

    def write(self):
        print("writing to "+str(self.name))
        self.tree.write(self.name)

    def get_keys(self, tag):
        sp = self.root.findall('CompositeSpectrum')[0]
        children = sp.findall(tag)
        if children==[]:
            raise Exception(str(children)+"does not contain "+str(tag))
        return list(dict(children[0].attrib).keys())

def prettify(elem):
    reparsed = minidom.parseString(et.tostring(elem, 'utf-8'))
    return reparsed.toprettyxml(indent="  ")


