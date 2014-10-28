import xml.etree.ElementTree as et
import data_types
import weakref
import uuid

class Flyweight(type):
    """a flyweight wrapper class"""
    def __init__(self,cls):
        self._instances = {}
        self._cls = cls
    def __call__(self,*args,**kwargs):
        self._instance.setdefault((*args,tuple(kwargs.items())), 
                                    self._cls(*args,**kwargs))
#class Singleton(type):
#    _instances = {}
#    def __call__(cls, *args, **kwargs):
#        if cls not in cls._instances:
#            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
#        return cls._instances[cls]


#@Flyweight
class ObjList(object):
    """this and associated subclasses are simply extended lists"""

    _pool=weakref.WeakValueDictionary()
    #key can be uuid, val will be object instance?

    def __new__(cls, **kwargs):
        obj=ObjList(**kwargs)
        if not obj in ObjList._pool.values(): #if already in pool
            obj = object.__new__(cls)
            ObjList._pool[obj.uuid] = obj
        else:
            obj=_pool[obj.uuid]  #return the old object, new one will be garbage collected
        return obj

    def __init__(self,objlst,taken_names=[],id=None):
        self.cls = objlst[0].__class__
        self.name = self.cls.__name__+"List"
        self.nodelist = [item.node for item in objlst]
        self.objlist = objlst
        if id==None:
            self.id=ObjList.generate_id(taken_names)    
        else:
            self.id=id  

    def __eq__(self,other):
        for i in range(len(self.objlst)):
            try:
                assert(self.objlst[i]==other.objlst[i])
            except IndexError:
                return False
            except AssertionError:
                if not self.objlst[i] in other.objlst:
                    return False
            except:
                raise
        return len(other.objlst)==len(self.objlst)  #last check     
        
    def __neq__(self,other):
        return not self.__eq__(other)

    def __iter__(self):
        for i in range(len(self.objlist)):
            yield self.objlist[i]               

    def __getitem__(self,i):
        return self.objlist[i]

    @staticmethod
    def factory(objlst,**kwargs):
        for cls in ObjList.__subclasses__():
            if cls.registrar_for(objlst[0].__class__.__name__):
                return cls(objlst,**kwargs)

    @staticmethod
    def generate_id(taken_names):
        iden=uuid.uuid4()
        while iden in taken_names:
            iden=uuid.uuid4()  
        return str(iden)

    def xml_rep(self,parent):
        """return the list of all relevant nodes in xml"""
        current = et.SubElement(parent,self.name,{"id":self.id})
        current.extend(self.nodelist)
        return parent

    @staticmethod
    def xml_read(parent):
        """read from and ModelDb xml file"""
        for cls in ObjList.__subclasses__():
            theTag = cls.__name__
            for sublist in parent.findall(theTag):  #find all of one sub-type of ObjLst
                #make list of allelements (all individual absorbers for example)
                theID = sublist.get("id")
                objlist = [data_types.Data.factory(**{"node":item}) for item in sublist]
                yield cls(objlist,id=theID)
            

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

    
