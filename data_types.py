import xmlutils
import warnings
import uuid
import xml.etree.ElementTree as et
import collections
import os
import sys
from scipy import constants

tf = {"true":True, "false":False}
c=constants.c/1000.  #speed of light in km/s

#TODO add functionality to directly change from absorber and have
#that propagate to all relevant entities...basically an observer pattern

class MissingPoolKey(Exception):
    """exception for when the pool is missing a key"""
    pass

class ObjList(object):

    _pool = dict() #stores objects already initialized

    def __init__(self, objlist, *args, **kwargs):
        self.objlist=objlist
        self._id=kwargs.pop("id",None)
        for key, val in kwargs.items():
            try:
                assert(type(val) in [str, float, int])
            except AssertionError:
                raise AssertionError("invalid type specified: %s:%s"%(
                    str(type(val)), str(val)))
            setattr(self,key,val)

    def __iter__(self):
        for i in range(len(self.objlist)):
            yield self.objlist[i]               

    def __getitem__(self,i):
        return self.objlist[i]

    def __len__(self):
        return len(self.objlist)

    def __eq__(self, other):
        """doensn't matter whether or not other is list or ObjList"""
        for item in other:
            if not item in self.objlist:
                return False
        return len(other) == len(self.objlist)

    def __neq__(self, other):
        return not self.__eq__(other)

    @staticmethod
    def clean_pool(lst):
        ObjList._pool = {k: v for (k, v) in ObjList._pool.items() if k in lst}

    @staticmethod
    def generate_id():
        iden=str(uuid.uuid4())
        while iden in list(ObjList._pool.keys()):
            iden=str(uuid.uuid4())  

        #ObjList.taken_names.append(iden)
        return iden

    @staticmethod
    def remove(iden):
        del(ObjList._pool[iden])

    def get_item(self,iden):
        """get element from objlist.  
        if AbsList, this will be the user assigned id, not the _pool key"""
        for item in self.objlist:
            if iden==item.id:
                return item

    @staticmethod
    def get(theID):
        """retrieve from _pool
        input:
        ------
        theID:  the objList's assigned ID

        output:
        -------
        Objlist or subclass of Objlist instance 

        raises:
        -------
        MissingPoolKey: when theID isn't recognized

        """
        try:
            return ObjList._pool[theID]
        except KeyError:
            msg = "\n  key not found: %s"%(str(theID))
            msg+="\n\n  available keys are:\n%s\n"%(str(sorted(ObjList._pool.keys())))
            raise MissingPoolKey(msg)
            #print("\n  key not found: %s"%(str(theID)))

    @staticmethod
    def set(value): 
        ObjList._pool[value.id] = value

    @property
    def id(self):
        return self._id

    @id.setter
    def id(self, value):
        """when setting a new id, changes should also be reflected in _pool"""
        if value in list(ObjList._pool.keys()):
            raise KeyError(str(value)+" already in pool")
        try:
            old = self._id
            assert(old in list(ObjList._pool.keys()))
        except AttributeError:
            raise Exception("trying to change id of data not yet registered in pool")
        except AssertionError:
            raise KeyError("attempting to access data not registered in pool")

        self._id = str(value)
        ObjList._pool[self._id] = ObjList._pool.pop(old)
       
    @staticmethod
    def factory(objlist,**kwargs):
        """
        factory method to produce instances of the chidren of ObjList.  If user 
        provides id=some value alread in _pool, then just return that data.
        new data will be registered in _pool

        Input:
        ------
        objlist : a list of data.  elements should be instances of some derived 
                  class of Data

        in **kwargs:
        ------------ 

        id : (default None) if specified, _pool will be searched for the 
             relevant data.  If not present, id will be set to specified value
             and a new instance will be created

        other entries will be suitable for any Data subclasses

        Output:
        -------
        object instance of desired subclass

        Raises:
        -------
        TypeError when cls isn't recognized

        """
        if objlist==None or objlist==[]:
            return None
        for key, value in ObjList._pool.items():
            if objlist==value:
                return value

        for cls in ObjList.__subclasses__():
            if cls.registrar_for(ObjList.classname(objlist[0].__class__)):
                iden=kwargs.pop('id',ObjList.generate_id()) 
                obj=cls(objlist,id=iden,**kwargs)

                #register data into the pool
                if not obj.id in list(ObjList._pool.keys()):
                    ObjList._pool[obj.id] = obj
                    return obj
                else:
                    return ObjList._pool[obj.id]
        raise TypeError("invalid type: "+ObjList.classname(objlist[0].__class__))

    @staticmethod
    def classname(cls):
        return cls.__name__.split('.')[-1]

    def xml_rep(self,parent):
        """return the list of all relevant nodes in xml"""
        try:
            assert(type(parent) in [str, unicode] )
            parent=str(parent)
        except:
            parent=ObjList.classname(parent.__class__)  #if a class is given instead of class name, just get the class name
            parent = parent.split('.')[-1]
        current = et.SubElement(parent,ObjList.classname(self.__class__),{"id":self.id})
        current.extend([item.node for item in self.objlist])
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
    def get_class(node):
        """
        get class of specified node.

        Input:
        ------
        node : xml node of data

        Output:
        -------
        the class that matches the specified node tag

        Raises:
        -------
        Exception if tag not recognized
        """
        for item in ObjList.__subclasses__():
            if node.tag==ObjList.classname(item):
                return item
        raise Exception("class %s not recognized"%(node.tag))

    @staticmethod
    def _read_node(node, verbose=False):
        """
        given a node that is an AbsorberList or similar type parse individual 
        datum

        input:
        ------
        node: node of one of AbsorberLists, etc...

        output:
        -------
        ObjList (or subclass) instance
        """ 
        if verbose:
            print(len(list(node)))
        objlist = [Data.factory(**{"node":item}) for item in list(node)]  #each item should be xml node
        return ObjList.factory(objlist,id=node.get("id"))


    @staticmethod
    def sublass_str():
        """return subclass names as list of strings"""
        return [item.__name__ for item in ObjList.__subclasses__() ]

    @staticmethod
    def get_from_xml(id,parent):
        """
        returns ObjList object instance from a given parent.  When _read_node 
        is called, data is automatically entered into ObjList._pool.  

        Input:
        ------
        parent : xml absorberLists node
        id : unique identifier of the data.

        Output:
        -------
        ObjList instance 

        """
        for item in list(parent):
            if item.tag in [it+"s" for it in ObjList.subclass_str()]:
                if item.id==id:
                    return ObjList._read_node(item)
        return None

    def refresh_list(xmlfile):
        ObjList._pool = {}
        tree = et.parse(xmlfile)
        root = tree.getroot()
        for tag in "AbsorberLists ContinuumPointLists RegionLists SingleViewLists VelocityViewLists".split():
            parent = root.find(tag)
            for objlist in parent:
                inst=ObjList._read_node(objlist)
                
    @staticmethod
    def list_from_xml(parent,verbose=False):
        """
        returns a list of ObjList objects from a given parent.  When _read_node 
        is called, data is automatically entered into ObjList._pool.  

        Input:
        ------
        parent : xml absorberLists node (list of absorberList instances)

        Output:
        -------
        list of class instances derived from ObjList 

        """
        if verbose:
            print(len(parent),' elements')
        return [ObjList._read_node(item, False) for item in parent]

    @staticmethod   
    def get_all_instances(subclass):
        lst = []
        if not type(subclass) is str:
            raise Exception("subclass input should be a string")
        for val in dict(ObjList._pool).values():
            if ObjList.classname(val.__class__) == subclass:
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

    def consolidate(self):
        """
        adapted from:
        http://codereview.stackexchange.com/questions/21307/consolidate-list-of-ranges-that-overlap
        """
        self.objlist = sorted(self.objlist, key=lambda x: x.start)
        ranges=[(item.start, item.end) for item in self.objlist]

        result = []
        current_start = -1
        current_end = -1 

        for start, end in sorted(ranges):
            if start > current_end:
                # this segment starts after the last segment stops
                # just add a new segment
                result.append( (start, end) )
                current_start, current_end = start, end
            else:
                # segments overlap, replace
                result[-1] = (current_start, end)
                # current_start already guaranteed to be lower
                current_end = max(current_end, end)

        for i in range(len(self.objlist)):
            if i<len(ranges):
                self.objlist[i].start=ranges[i][0]
                self.objlist[i].end=ranges[i][-1]

        while len(self.objlist)>len(ranges):
            del(self.objlist[-1])
            
     
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
        return id(self)

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
                inst.keys=list(inst.node.attrib.keys())
                inst.parse_node()
                #inst.set_data(**kwargs)
                try:  #if an additional constructor is specified, run it
                    inst.alt_init(**kwargs)
                except: 
                    pass
                return inst
        raise ValueError(
                "\n%s not a valid data type.  \nValid types are \n  %s"%(
                str(tag), str(Data.__subclasses__())
                ))

    @classmethod
    def from_file(cls,**kwargs):
        """constructor from file"""
        tag=kwargs.pop("tag")
        id=kwargs.pop("id")
        xmlfile = kwargs.pop("xmlfile")
        #node = xmlfile.get_node(id=id, tag=tag)

        root = et.parse(xmlfile).getroot()
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

    @staticmethod
    def read(filename,tag='Absorber', ids=None):
        
        """
        read the data.
        filename: input filename
        tag (default='Absorber'): the specified data type
        ids: which ids to read.  if None, then read all available.
        """

        if not type(filename) is str:
            filename.seek(0)
        etree=et.parse(filename)
        duderoot = etree.getroot()  ##should be SpecTool
 
        dudespec = duderoot.find("CompositeSpectrum")
        if dudespec is None:
            raise Exception("error reading fit file.  check %s to verify."%(filename))
        spectrum = dudespec.find("Spectrum")
        lst = []
        specdata=list(dudespec)+list(duderoot)

        if ids is None:
            for item in specdata:
                if item.tag == tag:
                    lst.append(Data.factory(node=item))
        else:
            try:
                if len(ids)==0: raise Exception('need at least one id specified')
            except:
                raise Exception("type of ids needs to be list of str")
            for item in specdata:
                if item.tag == tag and item.get('id') in ids:
                    lst.append(Data.factory(node=item))

        return lst

    def locked(self,param):
        return getattr(self,param+"Locked")

    def parse_node(self,node=None):
        """read from node, set attribs to self"""
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
        for key, val in dict(kwargs).items():
            self.node.set(key, str(val))

    def set_data(self,**kwargs):
        """set kwargs to self and apply to node"""
        for key, val in list(kwargs.items()): 
            if "Locked" in key:
                if not type(val) is bool:
                    assert(type(val) in [str, unicode])
                    assert(val in tf.keys())
                self.node.set(key, str(val).lower())
            setattr(self,key,val)
            self.node.set(key, str(val))


class Absorber(Data):
    node_attrib=["id","ionName",
                "N","NLocked","NError",
                "b","bLocked","bError",
                "z","zLocked","zError"]
    

    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)

    @classmethod
    def registrar_for(cls,tag):
        return tag=="Absorber"

    def __str__(self):
        return "%-5s id=%-6s N=%8.5lf b=%8.5lf z=%10.8lf"%(
                    self.ionName,self.id,self.N,self.b,self.z)



    def getShift(self, z):
        return (float(self.z) - z)*c/(1.+z)

    def get_wave(self, n=0):
        """get observed wave.  n=transition level such that 
            lya (n=2-->n=1)=0
            lyb (n=3-->n=1)=1
            and etc..
        """
        return (1.+self.z)*atomic_data[self.ionName][n].wave

    def get_obs_wave(self,n=0):
        return self.get_wave(n)

    def get_rest_wave(self, n=0):
        """get observed wave.  n=transition level such that 
            lya (n=2-->n=1)=0
            lyb (n=3-->n=1)=1
            and etc..
        """
        return atomic_data[self.ionName][n].wave

    def get_lines(self):
        """
        return list of SpectralLine instances containing relevant atomic data
        """
        lst=[]
        for item in atomic_data[self.ionName]:
            kwargs={}
            for key in ['f','gamma']:
                kwargs[key]=getattr(item,key)
            for key in Absorber.node_attrib:
                kwargs[key]=getattr(self,key)
            kwargs['wave']=float(getattr(item,'wave'))
            kwargs['obs_wave']=(1.+self.z)*float(getattr(item,'wave'))
            lst.append(SpectralLine(**kwargs))
        return lst

    @staticmethod
    def get_lyman_limit_tau(wave,z,N):
        """get tau for lyman limit systems"""
        lambda_l=911.8      
        if wave > lambda_l - 0.1: 
            wave = lambda_l - 0.1;

        acotz = np.atan(1./z)
        t2 = np.exp(-4.*z*acotz) / (1. - np.exp(-2.*np.pi*z))
        gaunt = 8. * np.pi * np.sqrt(3.) * (wave/lambda_l) * t2;
        return np.pow(10.,N) * ((wave/lambda_l)**3.) * 7.91e-18 * gaunt;
        
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
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
    def __str__(self):
        return "id=%s x=%12.7lf y=%12.8E"%(self.id,self.x,self.y)
    @classmethod
    def registrar_for(cls,tag):
        return tag=="ContinuumPoint"
        
class Region(Data):
    node_attrib=["start","end"]
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
    def __contains__(self, key):
        return self.start<=float(key)<=self.end
        
    @classmethod
    def registrar_for(cls,tag):
        return tag=="Region"

class SingleView(Data):
    node_attrib=["id","centWave","waveRange","minFlux","maxFlux"]
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
    @classmethod
    def registrar_for(cls,tag):
        return tag=="SingleView"

class VelocityView(Data):
    node_attrib=["id","labels",
                "minWave","maxWave","minFlux","maxFlux",
                "restWaves","redshift"]
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
    @classmethod
    def registrar_for(cls,tag):
        return tag=="VelocityView"

#class Singleton(type):
#    _instances = {}
#    def __call__(cls, *args, **kwargs):
#        if cls not in cls._instances:
#            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
#        return cls._instances[cls]

class AtomicData(object):
    def __init__(self):
        self.atomic_data = AtomicData.get_lines()

    @staticmethod
    def get_lines(fname='data/atom.dat'):
        all_lines={}
        f=open(fname,'r')
        for line in f:
            ion, line= line[0:5].strip(), line[5:]
            line=line.split()
            try:
                all_lines[ion].append(
                    SpectralLine(**{
                        'ionName':ion,
                        'wave': float(line[0]),
                        'f':    float(line[1]),
                        'gamma':float(line[2])
                    }))
            except KeyError:
                all_lines[ion]=[
                    SpectralLine(**{
                        'ionName':ion,
                        'wave': float(line[0]),
                        'f':    float(line[1]),
                        'gamma':float(line[2])
                    })]
        for k in all_lines.keys():
            all_lines[k] = sorted(all_lines[k], 
                                  key=lambda item:item.wave, reverse=True)
        return all_lines



class SpectralLine(object):
    mass_dict={
            'H':1.008,
            'D':2.014,
            'He':4.003,
            'Li':6.941,
            'Be':9.012,
            'B':10.081,
            'C':12.011,
            'N':14.007,
            'O':15.999,
            'F':18.998,
            'Ne':20.180,
            'Na':22.990,
            'Mg':24.305,
            'Al':26.982,
            'Si':28.086,
            'P':30.974,
            'S':32.066,
            'Cl':35.453,
            'Ar':39.948,
            'K':39.098,
            'Ca':40.078,
            'Sc':44.956,
            'Ti':47.867,
            'Cr':51.996,
            'Mn':54.938,
            'Fe':55.845,
            'Co':58.933,
            'Ni':58.693,
            'Cu':63.546,
            'Zn':65.38,
            'Ge':72.631
    }

    def __init__(self,**kwargs):
        for key, val in kwargs.items():
            setattr(self,key,val)

    def get_obs(self,z=None):
        try:
            if not z:
                z=self.z
            return (1.+z)*self.wave
        except:
            raise Exception(str(self.__dict__))

    @staticmethod
    def get_atom_name(ionName):
        ionName=ionName[:2]
        if ionName[-1] in 'IV ':
            return ionName[0]
        else:
            return ionName


#call on import
atomic_data=AtomicData().atomic_data  

