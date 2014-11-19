import xmlutils
import warnings
import weakref
import uuid
import xml.etree.ElementTree as et
import collections
import builtins

tf = {"true":True, "false":False}

class ObjList(object):
    """this and associated subclasses are simply extended lists which implement a flyweight"""

    #store refs in a dict used to check is inst of data already exists.  
    #dict will be used as access point for the flyweight

    #_pool=weakref.WeakValueDictionary()  
    _pool = dict()

    def __new__(cls, objlist,id=None):
        """implements a flyweight pattern using ObjList._pool"""
        obj = ObjList._check_list_in_list(objlist, list(ObjList._pool.values()) )
        if obj is None:  #behold, we have a new element in our midst
            obj = object.__new__(cls)
            #obj.id = ObjList.generate_id() if id==None else id
            obj.id = str(builtins.id(obj)) if id==None else str(id)
            try:
                assert(obj.id not in ObjList._pool.keys())
            except AssertionError:
                obj.id = ObjList.generate_id()

            ObjList._pool[obj.id] = obj
        return obj
        
    @staticmethod
    def _check_two_lists(lst1,lst2):
        """items are not guaranteed to be in same order ... but this will do for now"""
        assert("List" not in lst1.__class__.__name__)
        if (len(lst1)!=len(lst2)): return False
    
        for i in range(len(lst1)):
            if lst1[i] != lst2[i]:
                return False
        return True

    @staticmethod
    def _check_list_in_list(lst, superlst):
        """superlst is a list of lists.  does is contain lst?"""

        for item in superlst:
            try:
                temp = list(item.objlist)
            except AttributeError:
                return None
            if ObjList._check_two_lists(temp,lst):
                return item
        return None

    def __init__(self,objlist,id=None):
        self.cls = objlist[0].__class__
        self.objlist = objlist
        self.nodelist = [item.node for item in self.objlist]

    def __iter__(self):
        for i in range(len(self.objlist)):
            yield self.objlist[i]               

    def __getitem__(self,i):
        return self.objlist[i]

    @staticmethod
    def generate_id():
        iden=str(uuid.uuid4())
        while iden in ObjList._pool.keys():
            iden=str(uuid.uuid4())  

        #ObjList.taken_names.append(iden)
        return iden

    def get_item(self,id):
        """get element from objlist"""
        for item in self.objlist:
            if id==item.id:
                return item

    @staticmethod
    def get(theID):
        """retrieve from _pool"""
        return ObjList._pool[theID]

    @staticmethod
    def set(value): 
        ObjList._pool[value.id] = value

    @staticmethod
    def factory(objlist,**kwargs):
        for cls in ObjList.__subclasses__():
            try:
                if cls.registrar_for(objlist[0].__class__.__name__):
                    return cls(objlist,**kwargs)
            except IndexError:
                return None
        raise TypeError("invalid type: "+objlist[0].__class__.__name__)

    def xml_rep(self,parent):
        """return the list of all relevant nodes in xml"""
        current = et.SubElement(parent,self.__class__.__name__,{"id":self.id})
        assert(len(self.nodelist)>0)
        current.extend(self.nodelist)
        return current

    @staticmethod
    def _xml_read(parent):
        """read from and ModelDb xml file"""
        for sublist in parent.findall(theTag):  #find all of one sub-type of objlist
            #make list of allelements (all individual absorbers for example)
            theID = sublist.get("id")
            objlist = [data_types.Data.factory(**{"node":item}) for item in sublist]
            yield cls(objlist,id=theID)

    @staticmethod
    def _read_node(node):
        """read from and ModelDb xml file"""
        for cls in ObjList.__subclasses__():
            theTag = cls.__name__
            if node.tag==theTag:  
                objlist = [Data.factory(**{"node":item}) for item in node]
                return cls(objlist,id=node.get("id"))

    @staticmethod
    def sublass_str():
        """return subclass names as list of strings"""
        return [item.__name__ for item in ObjList.__subclasses__() ]

    @staticmethod
    def get_from_xml(id,parent):
        for item in parent:
            if item.tag in [it+"s" for it in ObjList.subclass_str()]:
                if item.id==id:
                    return ObjList._read_node(item)
        return None

    @staticmethod
    def list_from_xml(parent):
        """returns a list of ObjList objects from a given parent"""
        return [ObjList._read_node(item) for item in parent]

    @staticmethod   
    def get_all_instances(subclass):
        lst = []
        if not type(subclass) is str:
            raise Exception("subclass input shoulf be a string")
        for val in dict(ObjList._pool).values():
            if val.__class__.__name__ == subclass:
                lst.append(val)
        return lst
        
class AbsorberList(ObjList):
    @classmethod
    def registrar_for(cls,tag):
        return tag=="Absorber"

class ContinuumPointList(ObjList):
    @classmethod
    def registrar_for(cls,tag):
        return tag=="ContinuumPoint"

class RegionList(ObjList):
    @classmethod
    def registrar_for(cls,tag):
        return tag=="Region"

class SingleViewList(ObjList):
    @classmethod
    def registrar_for(cls,tag):
        return tag=="SingleView"

class VelocityViewList(ObjList):
    @classmethod
    def registrar_for(cls,tag):
        return tag=="VelocityView"

class Data(object):
    node_attrib=[]
    def __init__(self,tag,**kwargs):
        self.tag=tag
        for key,val in kwargs.items():
            try:
                setattr(self,key,float(val))
            except:
                setattr(self,key,val)

    def __eq__(self,other):
        if type(self)!=type(other):
            return False
        return self.node.attrib==other.node.attrib

    def __neq__(self,other):
        return not self.__eq__(other)

    def __hash__(self):
        return builtins.id(self)

    @staticmethod
    def factory(**kwargs):
        tag=kwargs.get("tag")
        if tag==None:
            try:
                tag=kwargs.get("node").tag
                assert(tag!=None)
            except:
                raise Exception("need to specify either tag and id or node")

        for cls in Data.__subclasses__():
            if cls.registrar_for(tag):
                if "node" in kwargs.keys():
                    inst=cls.from_node(**kwargs)
                elif "xmlfile" in kwargs.keys():
                    inst=cls.from_file(**kwargs)
                else:
                    raise Exception("need to specify either an xml element node or an xml file")
                inst.keys=inst.node.attrib.keys()
                inst.parseNode()
                inst.set_data(**kwargs)
                try:  #if an additional constructor is specified, run it
                    inst.alt_init(**kwargs)
                except: 
                    pass
                return inst
        raise ValueError("%s not a valid data type"%(str(tag)))

    @classmethod
    def from_file(cls,**kwargs):
        """constructor from file"""
        tag=kwargs.pop("tag")
        id=kwargs.pop("id")
        #xmlfile = xmlutils.Dudexml(kwargs.pop("xmlfile"))
        #node = xmlfile.get_node(id=id, tag=tag)

        root = et.parse(kwargs.pop("xmlfile")).getroot()
        node = xmlutils.get_node(root.find('CompositeSpectrum'),tag,id)
            

        return cls(tag,id=id,xmlfile=xmlfile,node=node,**kwargs)


    @classmethod
    def from_node(cls,**kwargs):
        """constructor from node"""
        node=kwargs.get("node")
        try:
            tag = kwargs.pop("tag")
        except:
            tag = node.tag
        return cls(tag,**kwargs)

    @staticmethod
    def get_node_attrib(cls, kwargs):
        _kwargs = kwargs
        for key in kwargs.keys():
            if not key in cls.node_attrib:
                del(_kwargs[key]) 
        return _kwargs

    def locked(self,param):
        return getattr(self,param+"Locked")
    """
    def parse_kwargs(self,**kwargs): #TODO get rid of this since set_data does the same thing
        for key, val in list(kwargs.items()):  #this goes after to override any conflicts
            if "Locked" in key:
                val = True if val=='true' else False
            setattr(self,key,val)
        self.set_node(**kwargs)
    """
    def parseNode(self,node=None):
        """read from node, set attribs to self"""
        try:
            assert(type(self) in Data.__subclasses__())
        except AssertionError:
            raise Exception(str(type(self))+" is not recognized")
        if node==None:
            node=self.node
        
        data = node.attrib
        for key, val in data.items():
            if "Locked" in key:
                setattr(self,key,tf[val.lower()])
            else:
                try:
                    setattr(self, str(key), float(val))
                except:
                    setattr(self, str(key), str(val))

    def set_node(self,**kwargs):  
        """set values from self to node"""
        #kwargs=Data.get_node_attrib(kwargs)      
        for key, val in list(kwargs.items()):
            if key in self.node.attrib.keys():
                old=self.node.get(key)
                self.node.set(key, str(val))
                new=self.node.get(key)

    def set_data(self,**kwargs):
        """set kwargs to self and apply to node"""
        for key, val in list(kwargs.items()):  #this goes after to override any conflicts
            if "Locked" in key:
                val="true" if val else "false" 
            setattr(self,key,val)
        self.set_node(**kwargs)

class Absorber(Data):
    node_attrib=["id","N","NLocked","NError","b","bLocked","bError","z","zLocked","zError","ionName"]
    
    @classmethod
    def registrar_for(cls,tag):
        return tag=="Absorber"

    def __str__(self):
        return "%-5s id=%-6s N=%8.5lf b=%8.5lf z=%10.8lf"%(self.ionName,self.id,self.N,self.b,self.z)

    def alt_init(self,**kwargs):
        self.obs=[item.get_obs(self.z) for item in atomic_data[self.ionName.replace(" ","")]]

    def locked(self,param):
        param_lock = {'N':'NLocked', 'b':'bLocked', 'z':'zLocked'}
        try:
            ans = getattr(self,param_lock[param])
            if str(ans) not in ['true', 'True', 'TRUE']:
                return False
            else:
                return True
        except KeyError:
            raise Exception("no param named %s"%{param})

    def in_region(self,regions):
        """
        detects whether of not CENTER of line is in a optimization region.
        only needs to be in one region to pass
        """
        try:
            lines=self.obs
        except AttributeError:
            warnings.warn("\n\n  \'alt_init\' wasn\'t called for \'Absorber\'\n\n")
            self.alt_init()
            lines=self.obs
        for region in regions:
            for line in lines: 
                if line in region:
                    return True
        return False

class ContinuumPoint(Data):
    node_attrib=["id","x","xLocked","xError","y","yLocked","yError"]
    def __str__(self):
        return "id=%s x=%12.7lf y=%12.8E"%(self.id,self.x,self.y)
    @classmethod
    def registrar_for(cls,tag):
        return tag=="ContinuumPoint"
        
class Region(Data):
    node_attrib=["start","end"]

    def __contains__(self, key):
        return self.start<=float(key)<=self.end
        
    @classmethod
    def registrar_for(cls,tag):
        return tag=="Region"

    @staticmethod
    def consolidate_regions(lst):
        newlst=[]
        lst = sorted(lst,key=lambda item:item.start, reverse=False)
        while len(lst)>0:
            try:
                if lst[0].end>=lst[1].start:
                    lst[0].end=lst[1].end
                    del(lst[1])
                newlst.append(lst.pop(0))     
            except IndexError:
                newlst.append(lst.pop())
        for item in newlst:
            item.set_node(start=item.start, end=item.end)
        return newlst

class SingleView(Data):
    node_attrib=["id","centWave","waveRange","minFlux","maxFlux"]

    @classmethod
    def registrar_for(cls,tag):
        return tag=="SingleView"

class VelocityView(Data):
    node_attrib=["id","redshift","minWave","maxWave","minFlux","maxFlux","restWaves","labels"]

    @classmethod
    def registrar_for(cls,tag):
        return tag=="VelocityView"

class Singleton(type):
    _instances = {}
    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]

class AtomicData(object, metaclass=Singleton):
    def __init__(self):
        self.atomic_data = AtomicData.get_lines()

    @staticmethod
    def get_lines(fname='atom.dat'):
        all_lines={}
        f=open(fname,'r')
        for line in f:
            line=line.split()
            try:
                all_lines[line[0]].append(SpectralLine(**{'wave':float(line[1]),'f':float(line[2])}))
            except KeyError:
                all_lines[line[0]]=[SpectralLine(**{'wave':float(line[1]),'f':float(line[2])})]
        for k in all_lines.keys():
            all_lines[k] = sorted(all_lines[k], key=lambda item:item.wave, reverse=True)
        return all_lines


class SpectralLine(object):
    def __init__(self,**kwargs):
        for key, val in kwargs.items():
            setattr(self,key,val)

    def get_obs(self,z):
        return (1.+z)*self.wave


atomic_data=AtomicData().atomic_data  
#having this here and calling atomic_data from the module makes the singleton unnecessary

