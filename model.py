"""
Implementation notes:

all xml fit file data (AbsorberList, ContinuumPointList, RegionList, 
SingleViewList, VelocityViewList) is written to the flyweight: 
data_types.ObjList._pool and will ONLY be accessed from there.

to access data here, we use Model.get(some id) where id is the key for our 
dict of stored objects:  data_types.ObjList._pool
"""

import xmlutils
import data_types
import warnings
import os.path
import matplotlib.pyplot as plt
import re
from constraints import Constraint
from numpy.random import random_sample, randn
import data_types
import xml.etree.ElementTree as et
import copy
import observer
import io
import pickle

c = 299792.458 #speed of light in km/s
def tf(val):
    return "true" if val in [True,"true"] else "false"

class Model(object): 
    #this dict maps ObjList subclasses to their associated data type
    model_classes = {"AbsorberList":"Absorber",
                     "ContinuumPointList":"ContinuumPoint",
                     "RegionList":"Region",
                     "SingleViewList":"SingleView",
                     "VelocityViewList":"VelocityView"}

    def __init__(self, **kwargs):
        self.pixels=0
        self.chi2=0.
        self.params=0
        self.flux=None  #fits (or text) file with flux
        self.error=None #fits (or text) file with flux.  if test, then this is the same as flux and xml format will change a bit
        self.xmlfile=kwargs.pop("xmlfile",None)
        self.get_all=kwargs.pop("get_all",True)
        self.id=kwargs.pop("id",data_types.ObjList.generate_id())
        self.abs_ids=kwargs.pop("abs_ids",None)


        if "buff" in kwargs.keys():
            self.xmlfile="scratch.xml"
            buff=kwargs.pop("buff")
            import io
            assert(type(buff) is io.BytesIO)
            self.read(buff=buff)
            
        else:
            for key, val in dict(kwargs).items():
                setattr(self,key,val)  #data lists will be the id needed to fetch from the data pool
            if not self.xmlfile:
                try:
                    self.xmlfile=Model.get(self.AbsorberList)[0].xmlfile.name
                except:
                    raise Exception("must specify either xml filename or model data")
            if not "AbsorberList" in kwargs.keys(): #no fitting data is specified in __init__, then read
                self.read()

        for key, val in kwargs.items():
            setattr(self, key, val)

        self.chi2=float(self.chi2)

        self.pixels=int(float(self.pixels))
        self.params=int(float(self.params))

        

        self.locked = {} 
        self._dof=float(self.pixels)-float(self.params)

    def __eq__(self,other):
        attrs = list(['id','chi2','pixels','params'])+list(Model.model_classes.keys())
        for item in attrs:
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
        chi2=123.4 pixels=154 params=23

        """
        string = "-----------AbsorberList------------\n"
        for item in Model.get(self.AbsorberList):
            string+=str(item)+"\n"
        string+="\nchi2=%lf pixels=%lf params=%lf dh=%lf\n\n"%(
                float(self.chi2),float(self.pixels),float(self.params))
        return string

    @property
    def reduced_chi2(self):
        self._reduced_chi2 = float(self.chi2)/(self.dof)
        return self._reduced_chi2

    @reduced_chi2.setter
    def reduced_chi2(self,value):
        self._reduced_chi2 = float(value)

    @property
    def dof(self):
        self._dof=float(self.pixels-self.params)
        return self._dof

    @dof.setter
    def dof(self,value):
        warnings.warn("setting dof without adjusting pixels or params")
        self._dof=value

    @property
    def dh(self):
        d = self.get_datum('D','Absorber',param='N')
        h = self.get_datum('H','Absorber',param='N')
        self._dh=float(d)-float(h)
        return self._dh

    @dh.setter
    def dh(self,value):
        """shouldnt need this"""
        try:
            self._dh=value
        except AttributeError:
            raise Exception("no availble D/H")

    def get_vel(self,id1,id2):
        """
        get velocity shift between two absorbers

        input:
        ------
        id1, id2:  (string): ids of absorbers to use.  z2 will be our reference.
                    if vel>0, then id1 is red of id2
        output:
        -------
        velocity in km/s

        

        """
        z1 = float(self.get_datum(id1,'Absorber','z'))
        z2 = float(self.get_datum(id2,'Absorber','z'))
        return (z1-z2)*c/(1.+z2)

    def get_shift(self,id1,id2):
        """
        an alias for get_vel(id1, id2)

        """
        return self.get_vel(id1,id2)

    def copy(self):
        mod = copy.deepcopy(self)
        mod.id=data_types.ObjList.generate_id()
        return mod

    def build_xml(self,raw_data='',spname='',sptype=''):
        """build a dude-style xml for this model"""
        if raw_data=='':
            raw_data=self.flux
            if not raw_data:
                raise Exception("must specify location of fits source data")
        if spname=='':
            spname=raw_data
        sptype='fits' if spname.endswith('.fits') else 'ascii1'
        duderoot = et.Element("SpecTool",{"version":"1.0"})
        dudespec = et.SubElement(duderoot,"CompositeSpectrum", {"id":raw_data})
        spec= et.SubElement(dudespec,"Spectrum",{"spec":spname,"spectype":sptype})

        #add spectrum stuff
        for attr in ["ContinuumPointList","AbsorberList"]:
            try:
                lst = self.get_lst(attr)
            except:
                if getattr(self,attr) is None:
                    raise Exception("no valid continuum points or absorbers for %s"%(self.xmlfile))
                else: 
                    raise
            if lst is None:
                raise Exception(attr+" not found for "+self.xmlfile)
            dudespec.extend([item.node for item in lst])
                

        #the view stuff
        for attr in ["SingleViewList", "VelocityViewList", "RegionList"]:
            try:
                duderoot.extend([item.node for item in self.get_lst(attr)])
            except:
                pass
        return duderoot

    def get_datum(self,iden,tag,param=None):
        """get an individual datum from model's data list."""
        if tag in Model.model_classes.values(): 
            tag=inv_dict(tag)


        tst=self.get_lst(tag)
        assert(tst)

        for item in self.get_lst(tag):
            if iden==item.id:
                if param:
                    if "Locked" in param:
                        return bool(getattr(item,param))
                    try:
                        return float(getattr(item,param))
                    except:
                        return getattr(item,param)
                else:
                    return item
        raise Exception("item not found: %s"%(iden))
        
    @staticmethod
    def get(iden):
        """
        an alias for data_types.ObjList.get
        note that getting model_instance.*List only returns the id associated 
        with that particular object. to get the actual object, need to use
            Model.get(model.*List)

        further note that this is a staticmethod.


        input:
        ------
        iden: id (string) of object of interest

        output:
        -------
        an element from data_types.ObjList._pool, so a list of absorbers,
        continuumpoints, velocity views, region lists or single views

        raises:
        -------
        None

        """
        return data_types.ObjList.get(iden)

    def get_lst(self,attr):
        """just an alias for:"""
        return data_types.ObjList.get(getattr(self,attr))

    def check_vals(self):
        unphysical={"b":[1.0,50.],"N":[10.00,25.00]} 

        abslist = self.get(self.AbsorberList)
        for item in abslist:
            assert(isinstance(item,data_types.Absorber))
            for key, val in unphysical.items():
                try:
                    assert(val[0]<=float(getattr(item,key))<=val[1])
                except: 
                    raise Exception("%s %s: %s has unphysical value of %lf"%(item.ionName,item.id,key,(getattr(item,key))))

    def monte_carlo_set(self,iden, tag, val_range, param, gaussian=False):
        """set a param for Data `id` to a random value in val_range
        params:
        -------
        iden:  id or object instance of the item to be changed
        tag:   absorber, continuum point, etc..  see data_types.py
        val_range:  list or tuple of [min, max].  if gaussian=true, this will be 95% limits
        gaussian:  return random value on gaussian distribution?
        param: which param to alter?

        """
        a=float(val_range[0])
        b=float(val_range[1])
        if gaussian:
            twosigma=(val_range[1]-val_range[0])/2.
            sigma=twosigma/2.
            new = sigma*randn()+(a+b)/2.
        else:
            new = (b-a)*random_sample()+a
        self.set_val(iden,tag,**{param:new})

    def read(self,buff=None, filename=None):
        """read from xml fit file, apply attribs to self

        input:
        ------
        filename (default None)  filename to read.  if None, reads self.xmlfile

        output:
        -------
        None

        raises:
        -------
        Exception:  when xml parsing fails


        """
        if buff:
            filename=buff
        elif not filename:
            filename=self.xmlfile

            


        #read through the data:  
        #  if data is new, then create the new class instances and store data in
        #  ObjList._pool and point the model to the correct pool key. Otherwise,
        #  just point the model instance to the appropriate pool key

        for key, val in Model.model_classes.items():
            ids=self.abs_ids if val=="Absorber" else None #and not self.get_all else None  #if looking at absorbers, can specify which absorbers to include
            lst = data_types.Data.read(filename, tag=val, ids=ids)
            if key=='AbsorberList': assert(len(lst)>0)
            newobj = data_types.ObjList.factory(lst)
            if newobj:
                setattr(self,key,newobj.id)


        # read thorugh xml file for source data files
        try:
            if not type(filename) is str:
                filename.seek(0)
            duderoot = et.parse(filename).getroot()  ##should be SpecTool
        except:
            if type(filename) is str:
                raise Exception(filename+" failed to parse.")
            else:
                raise
        compositespec = duderoot.find("CompositeSpectrum")
        path=os.path.split(compositespec.get("id"))[0]

        #find the source data
        spectrum = compositespec.find("Spectrum")

        try:
            self.flux=spectrum.get('spec')
        except:
            self.flux=compositespec.get('id')

        try:
            self.error=spectrum.get('error')
            self.chi2=float(spectrum.get("chi2"))
            self.params=float(spectrum.get("params"))
            self.pixels=float(spectrum.get("pixels"))
        except:
            self.chi2=0.
            self.params=0.
            self.pixels=0.
            self.error=None
        try:
            self.flux  = os.path.join(path,spectrum.get("spec"))
            self.error = os.path.join(path,spectrum.get("error"))
        except:
            self.flux  = os.path.join(path,spectrum.get("spec"))
            self.error = self.flux
 
        #read through all data in xml.  create relevant classes
        # enable this if you want each model to store its own data

        #for key, val in Model.model_classes.items():
        #    lst = []
        #    for item in list(dudespec)+list(duderoot):
        #        if item.tag == val:
        #            lst.append(data_types.Data.factory(node=item))
        #    setattr(self,key,data_types.ObjList.factory(lst))


    def set_val(self,iden,tag="Absorber",**kwargs):
        """
        set value of item in model

        input:
        ------
        iden (unspecified type):  the item or the id of the item to set
        tag (string, optional):  the type of iden

        output:
        -------
        None

        raises:
        -------
        TypeError when type of iden or tag is not recognized

        """

        if type(iden)==str:
            if tag in Model.model_classes.values(): 
                tag=inv_dict(tag)
            else:
                raise TypeError("model.Model.set_val(): unrecognized type: %s"%(tag))

            #get the model's absorber list
            ab_lst=data_types.ObjList.get(getattr(self,tag))

            for item in ab_lst:  #for item in abslist
                if item.id==iden:  #if the id of the absorber matches
                    item.set_data(**kwargs)
                    #print("setting to %s %s"%(iden, str(kwargs)))
                    """for key, val in dict(kwargs).items():
                        try:
                            assert(key in item.node.attrib.keys())
                            assert(self.get_datum(iden,"Absorber",key)==val)
                        except AssertionError:
                            print("\n\n%s %s %s: %s != %s\nkeys=%s\n\n"%(
                                    self.xmlfile,iden,key,str(val),
                                    str(self.get_datum(iden,"Absorber",key)), 
                                    str(item.node.attrib.keys())))
                            raise"""
        else: #iden was the object, not just the id
            tag=iden.__class__.__name__.split('.')[-1]
            if tag in Model.model_classes.values():
                iden.set_data(**kwargs)
            elif tag in list(Model.model_classes.keys()):
                for key, val in dict(kwargs).items:
                    setattr(iden,key,val)
            else:
                raise TypeError("model.Model.set_val(): unrecognized type: %s"%(tag))


        #new_lst= data_types.ObjList.factory(ab_lst)  #dont need to return this, since being automatically written to _pool
        #assert(new_lst.id in data_types.ObjList._pool.keys())  #test that data was added to pool
        #self.write(kwargs.get("filename", self.xmlfile+"_scratch.xml"))
            
    def write(self,filename=None):
        if filename is None: 
            filename=self.xmlfile


        with open(filename,"w") as f:
            f.write("<?xml version=\"1.0\"?>\n")
            f.write("<SpecTool version=\"1.0\">\n")
            if self.flux.endswith('.fits'):
                f.write("<CompositeSpectrum id=\"%s\"><Spectrum spec=\"%s\" error=\"%s\" chi2=\"%lf\" pixels=\"%lf\" params=\"%d\"/>\n"%(
                      self.flux,os.path.split(self.flux)[-1],os.path.split(self.error)[-1],self.chi2,self.pixels,self.params
                   ))
            else:
                f.write("<CompositeSpectrum id=\"%s\"><Spectrum spec=\"%s\" spectype=\"ascii1\" chi2=\"%lf\" pixels=\"%lf\" params=\"%d\"/>\n"%(
                      self.flux,os.path.split(self.flux)[-1],self.chi2,self.pixels,self.params
                   ))
            for item in data_types.ObjList.get(getattr(self,"ContinuumPointList")):
                if not item.id: item.id="null"
                f.write("<ContinuumPoint x=\"%lf\" y=\"%E\" xError=\"0.0\" yError=\"0.0\" xLocked=\"%s\" yLocked=\"%s\" id=\"%s\"/>\n"%(
                    item.x,item.y,tf(item.xLocked),tf(item.yLocked),item.id
                ))

            for item in data_types.ObjList.get(getattr(self,"AbsorberList")):
                f.write("<Absorber ionName=\"%s\" N=\"%lf\" b=\"%lf\" z=\"%lf\" NError=\"0.0\" bError=\"0.0\" zError=\"0.0\" NLocked=\"%s\" bLocked=\"%s\" zLocked=\"%s\" id=\"%s\"/>\n"%(
                        item.ionName, item.N, item.b,item.z,tf(item.NLocked),tf(item.bLocked),tf(item.zLocked),item.id
                        ))
            f.write("</CompositeSpectrum>\n")

            for item in data_types.ObjList.get(getattr(self,"RegionList")):
                f.write("<Region start=\"%lf\" end=\"%lf\"/>\n"%(item.start, item.end))
            f.write("</SpecTool>")

    def get_spectral_line(self, iden, transition):
        if type(iden) is int:  #gave an index instead of id
            ab=Model.get(self.AbsorberList)[iden]
        else:
            ab=Model.get(self.AbsorberList).get_item(iden)
        try:
            return ab.get_lines()[transition]
        except AttributeError:
            raise Exception("absorber %s not found"%(str(iden)))
        except KeyError:
            raise KeyError("absorber %s has no transition %d"%(str(iden), transition))
        


class ModelDB(object):
    def __init__(self, name=None, models=[], constraints=None,**kwargs): 
        """
        Model Database

        Inputs:
        -------
        models:  list of Model instances
        constraints: see ModelDB.constrain
        name: name of the xml models file.  (not the fit file)
        """

        for key, val in dict(kwargs).items():
            setattr(self,key,val)

        self.name=name
        #self.dbxml=xmlutils.Model_xml(filename=name)
        if len(models)>0:  #instantiate new db from models
            self.models = models
            for key in list(Model.model_classes.keys()):
                try:
                    setattr(self,key,[getattr(item,key) for item in self.models])
                except:
                    setattr(self,key,[])
            
        elif name:   
            if name.endswith('.xml'):
                self.models = ModelDB.read(str(name), returndb=False)
                

            else:
                self.models=ModelDB.load_models(name).models
        else:
            self.models = []

        if constraints:
            self.constrain(constraints)
#this breaks encapsulation, but I need a reference to pool so the data will pickle ok
        self.pool=data_types.ObjList._pool

    def __iter__(self):
        for i in range(len(self.models)):
            yield self.models[i]  

    def __len__(self):
        return len(self.models)             

    def __getitem__(self,i):
        return self.models[i]

    def append(self, model):
        self.models.append(model)

    def append_lst(self, lst,constraints=None):
        for item in lst:
            self.models.append(item)
        if constraints:
            self.constrain(constraints)

    def get_all_abs(self,iden,param,locked=False,constraints=None):
        """return a list of desired param values from all models"""
        x = []
        y = []
        if not constraints is None:
            self.constrain(constraints)
        lst = self.models
        if locked:
            for item in lst:
                abslist = Model.get(item.AbsorberList)
                ab = abslist.get_item(iden)
                if ab.locked(param):
                    x.append(float(getattr(ab,param)))
                    y.append(float(item.chi2))
        else:
            for item in lst:
                abslist = Model.get(item.AbsorberList)
                ab = abslist.get_item(iden)
                if ab is None:
                    raise Exception(item.xmlfile+" has no absorber: "+iden)
                try:
                    x.append(float(getattr(ab,param)))
                    y.append(float(item.chi2))
                except:
                    pass


        if len(x)==0 or len(y)==0 or len(x)!=len(y):
            raise Exception("ill condittioned input: \n  x=%s\n  y=%s"%(str(x),str(y)))

        return x, y
            
    def remove_unused(self):
        lst=[]
        for mod in self.models:
            for key in list(Model.model_classes.keys()):
                try:
                    lst.append(getattr(mod,key))
                except AttributeError:
                    pass
        data_types.ObjList.clean_pool(list(set(lst)))

        for mod in list(self.models):
            for key in list(Model.model_classes.keys()):
                try: #delete the entire model if not in _pool.  This breaks encapsulation
                    if not getattr(mod,key) in list(data_types.ObjList._pool.keys()):
                        self.models.remove(mod)
                except (AttributeError, ValueError):
                    pass

        

    def remove(self, model):
        """
        remove a model from the database.
        also removes non-repeated entries in ObjList._pool

        input:
        ------
        model: which model to remove

        output:
        -------
        None

        raises:
        -------
        None

        """

        
        def check_for_conflicts():
            def get_keys(mod):
                out=[]
                for item in list(Model.model_classes.keys()):
                    try:
                        out.append(getattr(mod,item))
                    except AttributeError:
                        pass
                return out


            these_keys=get_keys(model)
            for mod in self.models:
                if mod is model:
                    continue
                for item in get_keys(mod):  #check if any of these_keys are repeated in mod's keys
                    if item in these_keys:
                        these_keys.remove(item)
                    if len(these_keys)==0:
                        return []
            return these_keys

        these_keys=check_for_conflicts()
        for key in these_keys: #remove items from ObjList._pool
            data_types.ObjList.remove(key)
     
        self.models.remove(model)

    def append_db(self,dbfile):
        """appends another xmldb from filename to the current db"""
        tree = et.parse(dbfile)
        root = tree.getroot()

        for key in Model.model_classes.keys():
            parent = root.find(key+'s')
            objlist=data_types.ObjList.list_from_xml(parent) #instantiate all absorber/contpoint/view/etc 
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
                raise Exception(str(kwargs))


    @staticmethod
    def filter(inp, filters):
        for key, val in filters.items():    
            if key=='chi2':
                inp=[item for item in inp if item.chi2<val]
            else:
                if type(val) is dict:
                    for param, rng in val.items():
                        if param=='shift':
                            inp=[item for item in inp if rng[0] <= item.get_shift(key,'H') <= rng[-1]] 
                        elif type(rng) in [str, bool]:
                            inp=[item for item in inp if item.get_datum(key,'Absorber',param)==rng] 
                        else:
                            inp=[item for item in inp if rng[0] <= float(item.get_datum(key,'Absorber',param)) <= rng[-1]] 
                else:
                    inp=[item for item in inp if getattr(item, key)==val]
        return inp
                




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
                for attr in list(Model.model_classes.keys()):
                    try:
                        inst = getattr(item,attr)
                        assert(type(inst) is str)
                        data[attr]=inst
                    except: 
                        if attr in ['AbsorberList', 'ContinuumPointList']:
                            warnings.warn(
                                "model "+group_name+"has no attribute "+attr
                            )
                        data[attr]="None"

                current_group = et.SubElement(modeldb, 'model', data)

        #build the actuall fitting data
        for attr in list(Model.model_classes.keys()):  #write AbsorberList, cont points and views
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

    def constrain(self,constraints):
        """
        example constraints:   
            constraints={"chi2":123,"params":3,"pixels":2345,"D":{"N":(12.3,14.3),"b":(15,16)}}
        """

        if type(constraints) is dict: 
            constraints=Constraint(**constraints)

        for item in [it for it in self.models if not it in constraints]:
            self.remove(item)

    def get(self,xmlfile,chi2,pixels,params=None):
        """get from xml fit file"""
        mod = Model(xmlfile=xmlfile,chi2=chi2,pixels=pixels,params=params)
        #test for unphysical values 
        mod.check_vals()
            
        self.models.append(mod)

    def get_lst_from_id(self,iden,attr):
        """get all models with a given continuum"""
        return [item for item in self.models if getattr(item,attr)==iden]

    def get_best_lst(self, iden=None, param=None):
        """
        gets best fit.

        if `param' is specified, then gets best fit with a given parameter locked.
        param should be either bLocked, zLocked or NLocked
        """
        if not param is None:
            self.models=sorted(self.models, key=lambda x: x.chi2)
            return [(mod, mod.chi2) for mod in self.models]
        else:
            return self.get_locked(iden, param)  #already sorted

    def get_locked(self, iden, tag, param):
        tmp = []
       
        for mod in self.models:
            try:
                if to_bool(mod.get_datum(iden,tag,param+"Locked")):
                    tmp.append(mod)
            except:
                pass
        tmp =  sorted(tmp, key=lambda x: x.chi2)
        return [(mod.get_datum(iden,tag,param), mod.chi2) for mod in tmp]

    def get_model(self, iden):
        for item in self.models:
            if item.id==iden: 
                return item

    def get_vel(self,id1,id2):
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
    def read(filename, returndb=True, verbose=False,time=False):
        """read from xml db, return inputs for Model"""

        if not filename.endswith('.xml'):
            try:
                return ModelDB.load_models(filename)
            except:
                msg= "ModelDB.read: input file %s is"%(filename)
                msg+=" of incorrect format. Please verify"
                raise Exception(msg)
                    

        import time 
        t=time.time()
        tree = et.parse(filename)
        t0=time.time()
        if time: print("time to parse file: %lf"%(t0-t))
        root = tree.getroot()
        #root=xmlutils.Model_xml.get_root(filename)
        data_types.ObjList._pool = {}  #clear out cache of data
        #build all data first before instantiating individual models
        for key in list(Model.model_classes.keys()):
            if verbose:
                print('reading '+key+'s instances...')
            parent = root.find(key+'s')
            objlist=data_types.ObjList.list_from_xml(parent,verbose) #instantiate all absorber/contpoint/view/etc data. data stored in data_types.ObjList._pool
        #print("\n\n\n"+str(data_types.ObjList._pool.keys())+"\n\n")

        t1=time.time()
        if time: print("time to get objlist: %lf"%(t1-t0))

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
                raise Exception(str(kwargs))

        t2=time.time()
        if time:print("instantiate models: %lf"%(t2-t1))

        if len(model_list)==0:
            raise Exception("no models saved")

        if not returndb:
            return model_list
        else:
            return ModelDB(models=model_list)
        

    def set_val(self,model,**kwargs):
        """set the values of a given data element"""
        if type(model) is Model:
            for item in self.models:
                if item is model:
                    item.set_val(**kwargs)
            raise Exception("model.ModelDB.set_val:  model not found")
        elif type(model) is str:
            model = self.get_model(iden)
            model.set_val(**kwargs)
    
    @staticmethod
    def dump_models(db,fname=None):
        if not fname:
            fname=db.name
        if not fname.endswith(".obj"):
            fname+=".obj"
        if len(db.pool.keys())==0:
            raise data_type.MissingPoolKey("no data to dump from pool")
        pickle.dump(db,open(fname, "wb"))

   
    @staticmethod
    def load_models(fname):
        db = pickle.load(open(fname, "rb"))
        data_types.ObjList._pool=db.pool


        """if len(data_types.ObjList._pool.keys())==0:
            print("need to refresh pooled data.")
            try:
                data_types.ObjList.refresh_list(fname.replace('.obj','.xml'))
            except: 
                data_types.ObjList.refresh_list(input('file to unpickle: '))
        db.pool=data_types.ObjList._pool"""
        return db

    def write(self,filename=None,verbose=False):

        if filename==None:
            if self.name == None:
                filename=input("name of database file to write: ")
            else:
                filename=self.name
                if not filename.endswith('.xml'):
                    filename=filename+'.xml'

        root = ModelDB.build_xml(self.models)
        out = prettify(root)
        out.replace('\n\n','\n')
                    
        if verbose:
            print("writing %s with %d models"%(filename,len(self.models)))
        f = open(filename,'w')
        f.write(out)
        #tree = et.ElementTree(root)
        #tree.write(filename)
        #self.dbxml.write(filename,root)

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
            iden = item.get("id")
            assert(iden not in ids)
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




