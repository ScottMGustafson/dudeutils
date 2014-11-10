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
        inputs:
        -------
        AbsorberList: data_types.ObjList() [effectively used as list(data_types.Absorber)]
        """ 
        self.id=kwargs.pop("id",data_types.ObjList.generate_id())
        self.parse_kwargs(**kwargs)

        if self.xmlfile==None:
            try:
                self.xmlfile=self.AbsorberList[0].xmlfile.name
            except:
                raise Exception("must specify either xml filename or model data")

        #self.xml_fit = xmlutils.Dudexml(self.xmlfile)
        self.read()
        self.test_chi2()

        """
        if not all(x in kwargs.keys() for x in Model.model_classes.keys()):
            if self.xmlfile:
                self.get_model()
            else:
                raise Exception("need to define at least one argument")
        """
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
        for item in self.ContinuumPointList:
            string+=str(item)+"\n"

        locked = {}
        string = "-----------AbsorberList------------\n"
        for item in self.AbsorberList:
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

    def append(self,inst=None,tag=None,**kwargs):
        if not tag: raise Exception("must specify data type")
        if tag in Model.model_classes.values(): tag= inv_dict(tag)
        if inst:
            getattr(self,tag).append(inst)
            
        else:
            getattr(self,tag).append(data_types.Data.factory(**kwargs))

    def build_xml(self,raw_data='',spname='',sptype=''):
        """build a dude-style xml for this model"""
        #define the boilerplate
        duderoot = et.Element("SpecTool",{"version":"1.0"})
        dudespec = et.SubElement(duderoot,"CompositeSpectrum", {"id":raw_data})
        et.SubElement(duderoot,"Spectrum",{"spec":spname,"spectype":sptype})

        #add spectrum stuff
        for attr in ["AbsorberList", "ContinuumPointList"]:
            try:
                dudespec.extend([item.node for item in getattr(self,attr)])
            except:
                if getattr(self,attr) is None:
                    pass
                else: raise
                

        #the view stuff
        for attr in ["SingleViewList", "VelocityViewList", "RegionList"]:
            try:
                duderoot.extend([item.node for item in getattr(self,attr)])
            except:
                if getattr(self,attr) is None:
                    pass
                else: raise
        return duderoot

    def consolidate_regions(self):
        self.RegionList = data_types.Region.consolidate_regions(self.RegionList)
        self.write()

    def count_params(self):
        """get the number of params being optimized"""
        params=0
        for ab in self.AbsorberList:
            
            if ab.in_region(self.RegionList):
                for lock in ["bLocked","zLocked","NLocked"]:
                    if getattr(ab,lock)==True:
                        params+=1
        self.params = params

    def get(self,id,tag,param=None):
        if tag in Model.model_classes.values(): tag=inv_dict(tag)
        for item in getattr(self,tag):
            if id==item.id:
                if param:
                    try:
                        return float(getattr(item,param))
                    except:
                        return getattr(item,param)
                else:
                    return item
        raise Exception("item not found: %s"%(id))
        
    def get_model(self,**kwargs):
        """get all model data from xml and set attribs for self"""
        for key, val in Model.model_classes.items():
            lst=[]
            for item in data_types.ObjList.get(val):
                assert(item.tag==val)
                inp = {"xmlfile":self.xmlfile,"node":item,"tag":item.tag}
                lst.append(data_types.Data.factory(**inp))
            if len(lst)>0:
                setattr(self,key,data_types.ObjList.factory(lst))
            else:
                setattr(self,key,None)
        self.parse_kwargs(**kwargs)
        self.test_chi2()

    def monte_carlo_set(self,id,tag,param,val_range):
        """set a param for Data `id` to a random value in val_range"""
        a=float(val_range[0])
        b=float(val_range[1])
        new = (b-a)*random_sample()+a
        self.set_val(id,tag,**{param:new})

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

    def read(self,filename=None):
        if not filename:
            filename=self.xmlfile
        duderoot = et.parse(filename).getroot()  ##should be SpecTool
        dudespec = duderoot.find("CompositeSpectrum")
        spectrum = dudespec.find("Spectrum")
        spec, spectype = spectrum.get("spec"), spectrum.get("spectype")

        for key, val in Model.model_classes.items():
            lst = []
            for item in list(dudespec)+list(duderoot):
                #print(str(list(dudespec)+list(duderoot)))
                #raise AssertionError
                if item.tag == val:
                    lst.append(data_types.Data.factory(node=item))
            setattr(self,key,data_types.ObjList.factory(lst))

    def set_val(self,id,tag,**kwargs):
        """set the values of a given data element"""
        if tag in Model.model_classes.values(): tag=inv_dict(tag)
        for item in getattr(self,tag):
            if item.id==id:
                item.set_data(**kwargs)
                print("setting to %s"%(str(kwargs)))
        new_lst= data_types.ObjList.factory(getattr(self,tag))
        self.tree.write(self.xmlfile)
        #self.xml_fit.write()

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
            
    def write(self):
        #for item in attr:
            #get the node
        #    node=xmlutils.Dudexml.get_node(item.id,item.tag)
        #    for key in node.attrib.keys():
        #        node.set(key, str(item.key))
        root = self.build_xml()
        self.tree._setroot(root)
        self.tree.write(self.xmlfile)
        #self.xml_fit.write()


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
        if len(models)>0:
            self.models = models
            for key in Model.model_classes.keys():
                setattr(self,key,[getattr(item,key) for item in self.models])
            self.root=self.build_xml()
        elif name:   
            self.models = ModelDB.read(str(name), return_db=False)
            self.tree = et.parse(self.filename)
            self.root = self.tree.getroot()
            #self.root=self.dbxml.read(name)
        else:
            self.models = []

        if constraints:
            self.models = ModelDB.constrain(self,constraints)

    def append(self, model):
        self.models.append(model)

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
        x = []
        y = []

        if constraints!=None:
            lst=ModelDB.constrain(self,constraints)
            if len(lst)==0:
                raise Exception("no surviving models:\n%s"%(str(constraints)))
            if len(lst)==len(self.models):
                warnings.warn("everything passed")
        else:
            lst = self.models
    
        for item in lst:
            ab = item.get(id,"Absorber")
            if ab.locked(param):
                x.append(float(getattr(ab,param)))
                y.append(float(item.chi2))

        if len(x)==0 or len(y)==0 or len(x)!=len(y):
            raise Exception("ill condittioned input: \n  x=%s\n  y=%s"%(str(x),str(y)))

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

    def build_xml(self):
        """create an xmlfile file structure.  returns root"""
        import datetime
        now = str(datetime.datetime.now())
        root = et.Element('modeldb')

        #set up header data
        head = et.SubElement(root, 'head')
        title = et.SubElement(head, 'title')
        title.text = 'Fitting Models'
        modified = et.SubElement(head, 'dateModified')
        modified.text = now
        models = et.SubElement(root, 'ModelDB')

        #build individual models
        for item in self.models:
            current_group = None
            group_name = item.id 
            if current_group is None or group_name != current_group.text:
                data = {'id':group_name,'xmlfile':item.xmlfile,'chi2':str(item.chi2),
                    'pixels':str(item.pixels),'params':str(item.params)}
                for attr in Model.model_classes.keys():
                    inst = getattr(item,attr)
                    if inst==None: continue
                    data[attr]=inst.id
                current_group = et.SubElement(models, 'model', data)

        #build the actuall fitting data
        for attr in Model.model_classes.keys():  #write AbsorberList, cont points and views
            parent = et.SubElement(root,str(attr)+'s')
            instances = data_types.ObjList.get_all_instances(attr)
            #instances = [ item for item in data_types.ObjList._pool.values() if attr==str(type(item))]
            print(attr, len(instances))
            if attr=="AbsorberList": assert(len(instances)>0)
            parent.extend([ item.xml_rep(parent) for item in instances ])
            
            """for item in instances:
                child=et.SubElement(parent,item.__class__.name,{"id":item.id})
                for it in item:
                    node=et.SubElement(child,item.__class__.name,item.node.attrib)
                #child.extend(item.nodelist)
                #for it in item:
            """  
        return root

    @staticmethod
    def constrain(obj,constraints):
        """
        example constraints:   
            constraints={"chi2":123,"params":3,"pixels":2345,"D":{"N":(12.3,14.3),"b":(15,16)}}
        """
        constraint=Constraint(**constraints)
        return [item for item in obj.models if constraint.compare(item)]

    def get(self,xmlfile,chi2,pixels):
        self.models.append(Model(xmlfile=xmlfile,chi2=chi2,pixels=pixels))

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

    def get_err(self, id, param_name):
        """get 1 sigma error from chi2 = chi2min + 1

        param_name is in [N,b,z]

        """
        tmp = self.get_best_lst(param=param_name+'Locked')
        lst = [item[0] for item in tmp]
        chi2 = [item[1] for item in tmp]
        chi2min = float(chi2[0])
        onesig=[]
        for item in lst:
            if item.chi2<chi2min+1.:
                onesig.append(getattr(item.getabs(id),param_name))
        return getattr(lst[0].getabs(id),param_name) ,max(onesig), min(onesig)

    def get_locked(self, id, tag, param):
        tmp = []
       
        for mod in self.models:
            try:
                if to_bool(mod.get(id,tag,param+"Locked")):
                    tmp.append(mod)
            except:
                pass
        tmp =  sorted(tmp, key=lambda x: x.chi2)
        return [(mod.get(id,tag,param), mod.chi2) for mod in tmp]

    def get_min_chi2(self):
        return np.amin(np.array([item.chi2 for item in self.models]))

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
    """
    def parse_kwargs(self,**kwargs):
        for key, val in kwargs.items():
            try:
                setattr(self,key,float(val))
            except:
                setattr(self,key,val)
    """
    def pop(self,i):
        return self.models.pop(i)
    
    @staticmethod 
    def read(filename):
        """read from xml db, return inputs for Model"""
        
        tree = et.parse(filename)
        root = tree.getroot()
        #root=xmlutils.Model_xml.get_root(filename)
        data_types.ObjList._pool = {}  #clear out cache of data
        #build all data first before instantiating individual models
        for key in Model.model_classes.keys():
            parent = root.find(key+'s')
            objlist=data_types.ObjList.list_from_xml(parent) #instantiate all absorber/contpoint/view/etc data. data stored in data_types.ObjList._pool
        print("\n\n\n"+str(data_types.ObjList._pool.keys())+"\n\n")
        model_list = []
        models = root.find('ModelDB').findall('model')

        for model in models:   
#get model data (includes an id mapping to something in objlisr)
            kwargs = {}
            for key, val in dict(model.attrib).items(): 
                if key in Model.model_classes.keys(): #key is classname,  val is an id
                    kwargs[key] = data_types.ObjList.get(val)
                else:
                    kwargs[key] = val
            try:
                model_list.append(Model(**kwargs))
            except:
                raise Exception(str(kwargs))

        if len(models)==0:
            raise Exception("no models saved")

        return ModelDB(filename, models=model_list)
        

    def set_val(self,id,**kwargs):
        """set the values of a given data element"""
        new = self.models.pop(id)
        new.set_val(**kwargs)
        self.models.append(new)

    def write(self,filename=None):
        if filename==None:
            filename=self.name 
        root = self.build_xml()
        self.tree._setroot(root)
        self.tree.write(filename)
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

    




