import xmlutils
import warnings

tf = {"true":True, "false":False}

atomic_data=SpectralLine.instance().atomic_data

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
        for item in self.__class__.node_attrib:
            if getattr(self,item)!=getattr(other,item):
                return False
        return True

    def __neq__(self,other):
        return not self.__eq__(other)

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
        xmlfile = xmlutils.Dudexml(kwargs.pop("xmlfile"))
        node = xmlfile.get_node(id=id, tag=tag)
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

    def parse_kwargs(self,**kwargs): #TODO get rid of this since set_data does the same thing
        """set kwargs to self and apply to node"""
        for key, val in list(kwargs.items()):  #this goes after to override any conflicts
            if "Locked" in key:
                val = True if val=='true' else False
            setattr(self,key,val)
        self.set_node(**kwargs)
           
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

class Singleton(object):
    """
    adapted from:
    http://stackoverflow.com/questions/42558/python-and-the-singleton-pattern
    """
    def __init__(self,decorated):
        self._decorated=decorated

    def singleton(self):
        try:
            return self._instance
        except AttributeError:
            self._instance=self._decorated()
            return self._instance

    def __call__(self):
        raise TypeError("singletons should be called as cls.singleton()")

    def __instancecheck__(self, inst):
        return isinstance(inst, self._decorated)


@Singleton
class SpectralLine(object):
    def __init__(self,**kwargs):
        self.atomic_data=SpectralLine.get_lines()
        for key, val in kwargs.items():
            setattr(self,key,val)
        self.obs = None

    def set_obs_wave(self,z):
        self.obs=self.wave*(1.+z)

    def get_obs(self,z):
        return (1.+z)*self.wave

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

