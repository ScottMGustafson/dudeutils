"""
Implementation notes:

all xml fit file data (AbsorberList, ContinuumPointList, RegionList, 
SingleViewList, VelocityViewList) is written to the flyweight: 
data_types.ObjList._pool and will ONLY be accessed from there.

to access data here, we use Model.get(some id) where id is the key for our 
flyweight.

There shouldn't ever be the need to alter any model data, since it will be 
altered in dude, however, should the need arise, using Model.set_val will 
change the values and save as a new xml fit file (blah.xml_scratch.xml).  This 
is not written into data_types.ObjList._pool because it is expected that you 
will edit this in dude before using it as a viable model.

this will then be saved when you write it later.

"""



import xmlutils
import data_types
import warnings
import numpy as np
import matplotlib.pyplot as plt
import re
from constraints import Constraint
from numpy.random import random_sample
import data_types
import xml.etree.ElementTree as et



c = 299792.458

class Model(object): 
    #this dict maps ObjList subclasses to their associated data type
    model_classes = {"AbsorberList":"Absorber","ContinuumPointList":"ContinuumPoint",
                "RegionList":"Region","SingleViewList":"SingleView",
                "VelocityViewList":"VelocityView"}

    def __init__(self, **kwargs):
        """
        """ 
        self.id=kwargs.pop("id",data_types.ObjList.generate_id())
        #self.parse_kwargs(**kwargs)
        for key, val in dict(kwargs).items():
            setattr(self,key,val)  #data lists will be the id needed to fetch from the data pool
                 
        if self.xmlfile==None:
            try:
                self.xmlfile=Model.get(self.AbsorberList)[0].xmlfile.name
            except:
                raise Exception("must specify either xml filename or model data")

        #self.xml_fit = xmlutils.Dudexml(self.xmlfile)
        if not "AbsorberList" in kwargs.keys(): #no fitting data is specified in __init__, then read
            self.read()
        else:
            assert(type(self.AbsorberList) is str)
            assert(self.AbsorberList in data_types.ObjList._pool.keys())
        #self.test_chi2()

        for attr in ["chi2","pixels","params"]:
            try:
                assert(attr in kwargs.keys())
            except AssertionError:
                setattr(self,attr,0.)
                warnings.warn("%s wasn't defined.  setting to 0."%(attr))
        self.locked = {}
        if not self.params:
            self.count_params()

    def __eq__(self,other):
        for item in ['id','chi2','pixels','params','locked','AbsorberList']:
            if getattr(self,item)!=getattr(other,item):
                return False
        return True
    
    def __neq__(self,other):
        return not self.__eq__(other)


    def __str__(self):
        """
        string output for a model will be like:
    
        id=HI N=17.12345 b=12.345678 z=1.234567890
        id=SiI N=11.12345 b=12.345678 z=1.234567890
        id=OI N=11.12345 b=12.345678 z=1.234567890
        locked=OI:bLocked HI:zLocked
        chi2=123.4 pixels=154

        """
        string = "--------continuum points--------\n"
        for item in Model.get(self.ContinuumPointList):
            string+=str(item)+"\n"

        locked = {}
        string = "-----------AbsorberList------------\n"
        for item in Model.get(self.AbsorberList):
            string+=str(item)+"\n"
            for param in ['NLocked','bLocked','zLocked']:
                if getattr(item,param):
                    locked[item.absorber_id] = param
        string+="locked="
        for key, val in locked.items():
            if val:
                string+=str(key)+":"+str(val)+" "
        string+="\nchi2=%lf pixels=%lf params=%lf\n\n"%(float(self.chi2),float(self.pixels),float(self.params))
        return string

    """shouldn't need this
    def append(self,inst=None,tag=None,**kwargs):
        if not tag: raise Exception("must specify data type")
        if tag in Model.model_classes.values(): tag= inv_dict(tag)
        if inst:
            getattr(self,tag).append(inst)
            
        else:
            getattr(self,tag).append(data_types.Data.factory(**kwargs))
    """


    def build_xml(self,raw_data='',spname='',sptype=''):
        """build a dude-style xml for this model"""
        #define the boilerplate
        duderoot = et.Element("SpecTool",{"version":"1.0"})
        dudespec = et.SubElement(duderoot,"CompositeSpectrum", {"id":raw_data})
        et.SubElement(duderoot,"Spectrum",{"spec":spname,"spectype":sptype})

        #add spectrum stuff
        for attr in ["AbsorberList", "ContinuumPointList"]:
            try:
                dudespec.extend([item.node for item in Model.get(getattr(self,attr))])
            except:
                if getattr(self,attr) is None:
                    pass
                else: raise
                

        #the view stuff
        for attr in ["SingleViewList", "VelocityViewList", "RegionList"]:
            try:
                duderoot.extend([item.node for item in Model.get(getattr(self,attr))])
            except:
                if not getattr(self,attr):  #if the id is either '' or None
                    pass
                else: raise
        return duderoot

    """would be nice, but I wanna keep writing to _pool out of the picture except when instantiating new ObjList objects
    def consolidate_regions(self):
        self.RegionList = data_types.Region.consolidate_regions(self.RegionList)
        self.write()
    """
    def count_params(self):
        """get the number of params being optimized"""
        params=0
        for ab in Model.get(self.AbsorberList):
            if ab.in_region(Model.get(self.RegionList)):
                for lock in ["bLocked","zLocked","NLocked"]:
                    if getattr(ab,lock)==True:
                        params+=1
        self.params = params

    def get_datum(self,id,tag,param=None):
        """get an individual datum from model's data list.  (ex: an individual absorber)"""
        if tag in Model.model_classes.values(): tag=inv_dict(tag)
        for item in Model.get(getattr(self,tag)):
            if id==item.id:
                if param:
                    try:
                        return float(getattr(item,param))
                    except:
                        return getattr(item,param)
                else:
                    return item
        raise Exception("item not found: %s"%(id))
        
    @staticmethod
    def get(id):
        """just an alias for:"""
        return data_types.ObjList.get(id)

    def get_lst(self,attr):
        """just an alias for:"""
        return data_types.ObjList.get(getattr(self,attr))

    def check_vals(self):
        unphysical={"b":[1.0,50.],"N":[11.00,25.00]}
        abslist = self.get(self.AbsorberList)
        for item in abslist:
            assert(isinstance(item,data_types.Absorber))
            for key, val in unphysical.items():
                try:
                    assert(val[0]<=float(getattr(item,key))<=val[1])
                except: 
                    raise Exception("%s %s: %s has unphysical value of %lf"%(item.ionName,item.id,key,(getattr(item,key))))


    def get_model(self,**kwargs):  
        """get all model data from xml and set attribs for self"""
        for key, val in Model.model_classes.items():
            lst=[]
            for item in data_types.ObjList.get_all(val):
                assert(item.tag==val)
                inp = {"xmlfile":self.xmlfile,"node":item,"tag":item.tag}
                lst.append(data_types.Data.factory(**inp))  #creates list of absorbers, cont points, etc.
            if len(lst)>0:
                data=data_types.ObjList.factory(lst)  #store as objlist.  if already exists, does nothing by default
                setattr(self,key,data.id) #store key for flyweight
            else:
                setattr(self,key,None)  #else, store as none
        #self.parse_kwargs(**kwargs)
        self.test_chi2()

    def monte_carlo_set(self,id,tag,param,val_range):
        """set a param for Data `id` to a random value in val_range"""
        a=float(val_range[0])
        b=float(val_range[1])
        new = (b-a)*random_sample()+a
        self.set_val(id,tag,**{param:new})
    """
    def parse_kwargs(self,**kwargs):
        for key, val in kwargs.items():
            if val in data_types.ObjList.taken_names:
                setattr(self,key,data_types.ObjList.get(val))
            else:
                try:
                    setattr(self,key,float(val))
                except:
                    setattr(self,key,val)

    def parse_node(self,node):
        dat = self.xml.get_node_data(node=node)
        parse_kwargs(self,**dat)
    """

    """
    def pop(self,key,tag):
        if tag in Model.model_classes.values(): tag=inv_dict(tag)
        if type(key) == str:
            return self._pop_by_id(key,tag)
        elif type(key) == int:
            return self._pop_by_index(key,tag)
        else:
            lst = getattr(self,tag)
            for i in len(lst):
                if lst[i]==key:
                    ret = lst.pop(i)
                    setattr(self,tag,lst)
                    return ret
        raise Exception("item not found: %s"%(str(key)))

    def _pop_by_id(self,key,tag):
        lst = getattr(self,tag)
        for i in range(len(lst)):
            if lst[i].id==key:
                ret = lst.pop(i)
                setattr(self,tag,lst)
                return ret

    def _pop_by_index(self,i,tag):
        ret = getattr(self,tag).pop(i)
        setattr(self,tag,getattr(self,tag))
        return ret
    """

    def read(self,filename=None):
        """read from xml fit file, apply attribs to self"""
        if not filename:
            filename=self.xmlfile
        duderoot = et.parse(filename).getroot()  ##should be SpecTool
        dudespec = duderoot.find("CompositeSpectrum")
        if dudespec is None:
            raise Exception("error reading fit file.  check %s to verify."%(filename))
        spectrum = dudespec.find("Spectrum")
        spec, spectype = spectrum.get("spec"), spectrum.get("spectype")

        for key, val in Model.model_classes.items():
            lst = []
            for item in list(dudespec)+list(duderoot):
                #print(str(list(dudespec)+list(duderoot)))
                #raise AssertionError
                if item.tag == val:
                    lst.append(data_types.Data.factory(node=item))
            newobj = data_types.ObjList.factory(lst)
            assert(newobj.id in data_types.ObjList._pool.keys())  #test that data was added to pool
            setattr(self,key,newobj.id)

    def set_val(self,id,tag,**kwargs):
        """set the values of a given data element in xml file.  
           don't need to reset value in pool since data will be gotten
           later in get
        """
        if tag in Model.model_classes.values(): tag=inv_dict(tag)
        for item in getattr(self,tag):  #for item in abslist
            if item.id==id:  #if the id of the absorber matches
                item.set_data(**kwargs)
                print("setting to %s"%(str(kwargs)))
        new_lst= data_types.ObjList.factory(getattr(self,tag))  #dont need to return this, since being automatically written to _pool
        assert(new_lst.id in data_types.ObjList._pool.keys())  #test that data was added to pool
        self.write(kwargs.get("filename", self.xmlfile+"_scratch.xml"))

    def test_chi2(self):
        for item in ["chi2", "pixels", "params"]:
            try:
                float(getattr(self,item))
            except:
                if item=="params":
                    self.count_params()
                else:
                    warnings.warn(item+" not present...")
                    setattr(self,item,0.)
            
    def write(self,filename=None):
        #for item in attr:
            #get the node
        #    node=xmlutils.Dudexml.get_node(item.id,item.tag)
        #    for key in node.attrib.keys():
        #        node.set(key, str(item.key))
        if filename is None: filename=self.xmlfile
        root = self.build_xml()
        out = prettify(root)
        out.replace('\n\n','\n')
        f = open(filename,'w')
        f.write(out)
        
        #tree = et.ElementTree(root)
        #tree.write(self.xmlfile)


class ModelDB(object):
    def __init__(self, name=None, models=[], constraints=None,**kwargs): 
        """
        Model Database

        Inputs:
        -------
        models:  list of Model instances
        constraints: dict of dicts of tuples of floats.  (see ModelDB.constrain) 
        name: name of the xml models file.  (not the fit file)
        """

        for key, val in dict(kwargs).items():
            setattr(self,key,val)

        self.name=name
        #self.dbxml=xmlutils.Model_xml(filename=name)
        if len(models)>0:  #instantiate new db from models
            self.models = models
            for key in Model.model_classes.keys():
                setattr(self,key,[getattr(item,key) for item in self.models])
            
        elif name:   
            self.models = ModelDB.read(str(name), returndb=False)
        else:
            self.models = []

        if constraints:
            self.models = ModelDB.constrain(self.models,constraints)

    def __iter__(self):
        for i in range(len(self.models)):
            yield self.models[i]               

    def __getitem__(self,i):
        return self.models[i]

    def append(self, model):
        self.models.append(model)

    def get_all_abs(self,id,param,locked=False,constraints=None):
        """return a list of desired param values from all models"""
        x = []
        y = []
        if not constraints is None:
            lst=ModelDB.constrain(self.models,constraints)
            if len(lst)==0:
                raise Exception("no surviving models:\n%s"%(str(constraints)))
            if len(lst)==len(self.models):
                warnings.warn("everything passed")
        else:
            lst = self.models
        if locked:
            for item in lst:
                abslist = Model.get(item.AbsorberList)
                ab = abslist.get_item(id)
                if ab.locked(param):
                    x.append(float(getattr(ab,param)))
                    y.append(float(item.chi2))
        else:
            for item in lst:
                abslist = Model.get(item.AbsorberList)
                ab = abslist.get_item(id)
                x.append(float(getattr(ab,param)))
                y.append(float(item.chi2))


        if len(x)==0 or len(y)==0 or len(x)!=len(y):
            raise Exception("ill condittioned input: \n  x=%s\n  y=%s"%(str(x),str(y)))

        return x, y
            




    def best_fit(self,id,param,order,xmin=None,xmax=None, plot=True, constraints=None):
        """
        get a best fit of data with respect to `param'

        id: id of absorber
        param:  parameter name (N,b,z)
        order:  order of polynomial to fit
        xmax, xmin: range of values to consider
        locked:  get only locked parameters?
        plot:   plot the data?  otherwise return function, x, y
        """
        x, y = self.get_all_abs(self,id,param,locked=True,constraints=constraints)

        if xmin==None and xmax==None:
            xmax = max(x)
            xmin = min(x)

        x=np.array(x)
        y=np.array(y)

        coeffs=np.polyfit(x-x.mean(),y,int(order))
        f = np.poly1d(coeffs)

        if plot:
            xx = np.arange(xmin,xmax, np.abs(xmax-xmin)/100.)
            yy = f(xx-x.mean())
            plt.xlim(xmin,xmax)
            plt.plot(xx,yy,'b-')
            plt.plot(x,y,'ro')
            plt.show()
        return f, x, y

    @staticmethod
    def build_xml(models):
        """create an xmlfile file structure.  returns root"""
        import datetime
        now = str(datetime.datetime.now())
        #root = et.Element('modeldb')

        #set up header data
        root = et.Element('head')
        title = et.SubElement(root, 'title')
        title.text = 'Fitting Models'
        modified = et.SubElement(root, 'dateModified')
        modified.text = now
        modeldb = et.SubElement(root, 'ModelDB')

        #build individual models
        for item in models:
            current_group = None
            group_name = item.id 
            if current_group is None or group_name != current_group.text:
                data = {'id':group_name,'xmlfile':item.xmlfile,'chi2':str(item.chi2),
                    'pixels':str(item.pixels),'params':str(item.params)}
                for attr in Model.model_classes.keys():
                    inst = getattr(item,attr)
                    if not inst: continue
                    assert(type(inst) is str)
                    data[attr]=inst
                current_group = et.SubElement(modeldb, 'model', data)

        #build the actuall fitting data
        for attr in Model.model_classes.keys():  #write AbsorberList, cont points and views
            parent = et.SubElement(root,str(attr)+'s')
            instances = data_types.ObjList.get_all_instances(attr)
            #instances = [ item for item in data_types.ObjList._pool.values() if attr==str(type(item))]
            #print(attr, len(instances))
            if attr=="AbsorberList": assert(len(instances)>0)
            #parent.extend([ item.xml_rep(parent) for item in instances ])
            
            for item in instances:
                name=data_types.ObjList.classname(item.__class__)
                child=et.SubElement(parent,name,{"id":item.id})
                for it in item:
                    node=et.SubElement(child,it.node.tag,it.node.attrib)
                #child.extend(item.nodelist)
                #for it in item:
            
        return root

    def trim(self,constraints):
        self.models = ModelDB.constrain(self.models,constraints)

    @staticmethod
    def constrain(models,constraints):
        """
        example constraints:   
            constraints={"chi2":123,"params":3,"pixels":2345,"D":{"N":(12.3,14.3),"b":(15,16)}}
        """

        if type(constraints) is dict: 
            constraints=Constraint(**constraints)
        return [item for item in models if item in constraints]

    def get(self,xmlfile,chi2,pixels,params=None):
        """get from xml fit file"""
        mod = Model(xmlfile=xmlfile,chi2=chi2,pixels=pixels,params=params)
        #test for unphysical values 
        mod.check_vals()
            
        self.models.append(mod)

    def get_best_lst(self, id=None, param=None):
        """
        gets best fit.

        if `param' is specified, then gets best fit with a given parameter locked.
        param should be either bLocked, zLocked or NLocked
        """
        if not param is None:
            self.models=sorted(self.models, key=lambda x: x.chi2)
            return [(mod, mod.chi2) for mod in self.models]
        else:
            return self.get_locked(id, param)  #already sorted

    def get_locked(self, id, tag, param):
        tmp = []
       
        for mod in self.models:
            try:
                if to_bool(mod.get_datum(id,tag,param+"Locked")):
                    tmp.append(mod)
            except:
                pass
        tmp =  sorted(tmp, key=lambda x: x.chi2)
        return [(mod.get_datum(id,tag,param), mod.chi2) for mod in tmp]

    def get_model(self, id):
        for item in self.models:
            if item.id==id: 
                return item

    def get_vel_shift(self,id1,id2):
        return [item.get_vel(id1,id2) for item in self.models]

    def grab(self,xmlfile,chi2,pixels,params,**kwargs):
        """grab from xml file"""
        #need to reinstantiate xml file
        mod = Model(xmlfile=xmlfile,chi2=chi2,pixels=pixels,params=params,**kwargs)
        self.models.append(mod)
        return

    def pop(self,i):
        return self.models.pop(i)
    
    @staticmethod 
    def read(filename, returndb=True):
        """read from xml db, return inputs for Model"""
        
        tree = et.parse(filename)
        root = tree.getroot()
        #root=xmlutils.Model_xml.get_root(filename)
        data_types.ObjList._pool = {}  #clear out cache of data
        #build all data first before instantiating individual models
        for key in Model.model_classes.keys():
            parent = root.find(key+'s')
            objlist=data_types.ObjList.list_from_xml(parent) #instantiate all absorber/contpoint/view/etc data. data stored in data_types.ObjList._pool
        #print("\n\n\n"+str(data_types.ObjList._pool.keys())+"\n\n")
        model_list = []
        models = root.find('ModelDB').findall('model')

        for model in models:   
#get model data (includes an id mapping to something in objlisr)
            kwargs = {}
            for key, val in dict(model.attrib).items(): 
                kwargs[key] = val
            try:
                model_list.append(Model(**kwargs))
            except:
                raise
                raise Exception(str(kwargs))

        if len(models)==0:
            raise Exception("no models saved")

        if not returndb:
            return model_list
        else:
            return ModelDB(filename, models=model_list)
        

    def set_val(self,id,**kwargs):
        """set the values of a given data element"""
        new = self.models.pop(id)
        new.set_val(**kwargs)
        self.models.append(new)

    def write(self,filename=None):
        if filename==None:
            filename=self.name 
        root = ModelDB.build_xml(self.models)
        out = prettify(root)
        out.replace('\n\n','\n')
        f = open(filename,'w')
        f.write(out)
        #tree = et.ElementTree(root)
        #tree.write(filename)
        #self.dbxml.write(filename,root)

    @staticmethod
    def xml_in(**kwargs):
        """xml shorthand for an individual model"""
        pass
        

#    some helper functions:
#------------------------------
 
def inv_dict(tag, dic=Model.model_classes):
    if tag in dic.values():
        tmp = {v:k for k, v in dic.items()} 
        return tmp[tag]

def to_bool(string):
    string = string.lower()
    if string=="true":
        return True
    else:
        return False

def check_for_conflicts(root):
    """for some unknown reason, duplicates of a certain node will be printed.  
    this does not fix the underlying cause, but is a fix to prevent duplicate 
    printing.  If two nodes have the same id, but differing contents, an 
    exception will be raised"""
    ids = []
    for item in root:
        try:
            id = item.get("id")
            assert(id not in ids)
        except AssertionError:
            #code
            pass

    
def prettify(elem):
    """Return a pretty-printed XML string for the Element.
    copied and modified from:
    http://stackoverflow.com/questions/17402323/use-xml-etree-elementtree-to-write-out-nicely-formatted-xml-files
    """
    from xml.dom import minidom
    rough_string = et.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")




