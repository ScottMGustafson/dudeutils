import xmlutils

tf = {"true":True, "false":False}

class Data(object):
    def __init__(self,tag):
        self.tag=tag

    def __eq__(self,other):
        for item in self.__dict__.keys():
            if getattr(self,item)!=getattr(other,item):
                return False
        return True

    def __neq__(self,other):
        return not self.__eq__(other)


    @staticmethod
    def factory(**kwargs):
        tag=kwargs.get("tag")
        if tag==None:
            if "node" in kwargs.keys():
                node=kwargs.get("node")
                tag=node.tag
                assert(tag!=None)
            else:
                raise Exception("need to specify either tag and id or node")
        
        for cls in Data.__subclasses__():
            if cls.registrar_for(tag):
                if "node" in kwargs.keys():
                    inst=cls(tag)
                    inst.from_node(**kwargs)
                    return inst
                elif "xmlfile" in kwargs.keys():
                    inst=cls(tag)
                    inst.from_file(**kwargs)
                    return inst
                else:
                    raise Exception("need to specify either an xml element node or an xml file")
        raise ValueError("%s not a valid data type"%(str(tag)))

    def from_node(self,**kwargs):
        """constructor from node"""
        self.node=kwargs.get("node")
        self.parseNode()
        self.set_data(**kwargs)
        self.keys=self.node.attrib.keys()

    def from_file(self,**kwargs):
        """constructor from file"""
        tag=kwargs.pop("tag")
        id=kwargs.pop("id")
        xmlfile = xmlutils.Dudexml(kwargs.pop("xmlfile"))
        self.node = xmlfile.get_node(id=id, tag=tag)
        self.keys=self.node.attrib.keys()
        self.parseNode()
        self.parse_kwargs(**kwargs)

    def parse_kwargs(self,**kwargs): #TODO get rid of this since set_data does the same thing
        """set kwargs to self and apply to node"""
        for key, val in list(kwargs.items()):  #this goes after to override any conflicts
            if "Locked" in key:
                val = True if val=='true' else False
            setattr(self,key,val)
        self.set_node(**kwargs)
           
    def set_node(self,**kwargs):  
        """set values from self to node"""      
        for key, val in list(kwargs.items()):
            if key in self.node.attrib.keys():
                old=self.node.get(key)
                self.node.set(key, str(val))
                new=self.node.get(key)

#is there going to be an issue with id versus id?
    def parseNode(self,node=None):
        """read from node, set attribs to self"""
        try:
            assert(type(self) in Data.__subclasses__())
        except:
            print("Data Subclasses:")
            print(str(Data.__subclasses__()))
            print("this instance:")
            print(str(type(self)))
            raise
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

    def set_data(self,**kwargs):
        """set kwargs to self and apply to node"""
        for key, val in list(kwargs.items()):  #this goes after to override any conflicts
            if "Locked" in key:
                val="true" if val else "false" 
            setattr(self,key,val)
        self.set_node(**kwargs)

    def get_keys(self):
        return self.keys

    def locked(self,param):
        return getattr(self,param+"Locked")

class ContinuumPoint(Data):
    def __str__(self):
        return "id=%s x=%12.7lf y=%12.8E"%(self.id,self.x,self.y)
    @classmethod
    def registrar_for(cls,tag):
        return tag=="ContinuumPoint"
        
class Absorber(Data):
    @classmethod
    def registrar_for(cls,tag):
        return tag=="Absorber"
    def __str__(self):
        return "%-5s id=%-6s N=%8.5lf b=%8.5lf z=%10.8lf"%(self.ionName,self.id,self.N,self.b,self.z)

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
        
class VelocityView(Data):
    @classmethod
    def registrar_for(cls,tag):
        return tag=="VelocityView"

class SingleView(Data):
    @classmethod
    def registrar_for(cls,tag):
        return tag=="SingleView"

class Region(Data):
    @classmethod
    def registrar_for(cls,tag):
        return tag=="Region"

def consolidate_regions(lst):
    pass
    #TODO
