""" 
command line script to write new values to xml or to dump current vals into 
text file.
a typical command will be like:
    >>> python cmd_line_io.py get file.xml id=foobar id2=foo id3=bar
this will dump all the data from the listed absorbers to a file 
`absorberData.txt`.

    >>> python cmd_line_io.py write file.xml id=foobar N=12.5 b=8.3 z=4.0
will right the fitting parameters back into the xml
        

"""

import sys
import xmlutils
import data_types
import warnings
import numpy as np
import matplotlib as plt
import re
import astronomy_utils as astro


c = 299792.458
model_classes = {"absorbers":"Absorber","continuum_points":"ContinuumPoint","regions":"Region"}

class Model(object):
    
    def __init__(self, **kwargs):
        """
        inputs:
        -------
        absorbers: list(data_types.Absorber)  is a list of absorber references 
        """
        self.id=''
        self.parse_kwargs(**kwargs)
        if not all(x in kwargs.values() for x in ["absorbers", "continuum_points", "regions"]):
            if self.xmlfile:
                self.get_model(self.xmlfile)
            else:
                raise Exception("need to define at least one argument")
        if self.xmlfile==None:
            self.xmlfile=self.absorbers[0].xmlfile.name

        self.xml_fit = xmlutils.Dudexml(self.xmlfile)

        if "chi2" not in kwargs.keys() or "pixels" not in kwargs.keys():
            self.chi2=0.
            self.pixels=0.
            warnings.warn("chi2 and or pixels are 0.")
        self.locked = {}

    def __eq__(self,other):
        for item in ['id','chi2','pixels','locked','absorbers']:
            if getattr(self,item)!=getattr(other,item):
                return False
        return True
    
    def __neq__(self,other):
        return not self.__eq__(other)


    def __str__(self):
        """
        string output for a model will be like:
    
        iden=HI N=17.12345 b=12.345678 z=1.234567890
        iden=SiI N=11.12345 b=12.345678 z=1.234567890
        iden=OI N=11.12345 b=12.345678 z=1.234567890
        locked=OI:bLocked HI:zLocked
        chi2=123.4 pixels=154

        """
        string = ""
        locked = {}
        for item in self.absorbers:
            string+=str(item)+"\n"
            for param in ['NLocked','bLocked','zLocked']:
                if getattr(item,param):
                    locked[item.id] = param
        string+="locked="
        for key, val in locked.items():
            if val:
                string+=str(key)+":"+str(val)+" "
        string+="\nchi2=%lf pixels=%lf\n\n"%(float(self.chi2),float(self.pixels))
        return string

    def parse_kwargs(self,**kwargs):
        for key, val in kwargs.items():
            try:
                setattr(self,key,float(val))
            except:
                setattr(self,key,val)
        self.test_chi2

    def get_model(self,xmlfile=None,**kwargs):
        """get all model data from xml and set attribs for self"""
        xml = xmlutils.Dudexml(xmlfile) if xmlfile !=None else self.xml_fit
        for key, val in model_classes.items():
            lst=[]
            for item in xml.get_node_list(val):
                assert(item.tag==val)
                lst.append(data_types.Data.factory(xmlfile=xmlfile,node=item,tag=item.tag))
            if len(lst)>0:
                setattr(self,key,lst)
            else:
                setattr(self,key,None)
        self.parse_kwargs(**kwargs)
        self.test_chi2()

    def test_chi2(self):
        try:
            float(self.chi2)
            float(self.pixels)
        except:
            warnings.warn("chi2 and pixels not present...")
            self.pixels=0
            self.chi2=0.
            

    def constrain(self, constraints):
        """
        Inputs:
        -------
        constraints: a dict of constraints using Absorber attirbutes.  
        Each key-value pair is a string-tuple of floats.
        
        returns: boolean

        example:

        >>>ab=xmlutils.Absorber(N=12.2,b=10.,z=0.)
        >>>model=Model([ab])
        >>>model.constrain({ID:{N:(12,13), b:(9,11), z:(-1,1)}})
            True

        """
        for item in self.absorbers:
            try:
                for key, val in constraints[item.id].items():
                    if getattr(item,key)<val[0] or getattr(item,key)>val[1]:
                        return False
            except KeyError:
                pass
        return True

    def write(self):
        #ModelData to xml data:
        for key in model_classes.keys():
            for item in getattr(self,key):
                item.set_node(**item.__dict__) #set node values to current vals.
        self.xml_fit.write()

    def get(self,id,tag,param=None):
        if tag in model_classes.values(): tag= inv_dict(tag)
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
        
    def get_vel(self,iden1,iden2):
        z1 = self.get(iden1,'absorbers',"z")
        z2 = self.get(iden2,'absorbers',"z")
        return astro.get_vel_shift(z1,z2)
       
    def parse_node(self,node):
        dat = self.xml.get_node_data(node=node)
        parse_kwargs(self,**dat)

    def append(self,inst=None,tag=None,**kwargs):
        if not tag: tag=kwargs.get("tag")
        if not tag: raise Exception("must specify data type")
        if tag in model_classes.values(): tag= inv_dict(tag)
        if inst:
            getattr(self,tag).append(inst)
            
        else:
            getattr(self,tag).append(data_types.Data.factory(**kwargs))
            
    
    def pop(self,key,tag):
        if tag in model_classes.values(): tag=inv_dict(tag)
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

    def set(self,id,tag,**kwargs):
        """set the values of a given data element"""
        if tag in model_classes.values(): tag=inv_dict(tag)
        for item in getattr(self,tag):
            if item.id==id:
                item.set_node(**kwargs)
        

class ModelDB(object):
    def __init__(self, name=None, models=[], constraints=None,**kwargs): 
        """
        Model Database

        Inputs:
        -------
        models:  list of Model instances
        constraints: dict of dicts of tuples of floats.  (see Model.constrain) 
        name: name of the xml models file.  (not the fit file)
        """

        for key, val in dict(kwargs).items():
            setattr(self,key,val)

        self.name=name
        self.dbxml=xmlutils.Model_xml(filename=name)
        if len(models)>0:
            self.lst = models
            self.root=self.create(str(name))
        elif name:   
            self.lst = read_in(str(name), return_db=False)
            self.root=self.dbxml.read(name)
            
        else:
            self.lst = []

        if constraints:
            self.lst = [item for item in self.lst if item.constrain(constraints)]

    def pop(self,key):
        if type(key) == str:
            return self._pop_by_id(key)
        elif type(key) == int:
            return self._pop_by_index(key)
        else:
            for i in len(self.lst):
                if self.lst[i]==key:
                    return self.lst.pop(i)
        raise Exception("item not found: %s"%(str(key)))

    def _pop_by_id(self,key):
        for i in range(len(self.lst)):
            if self.lst[i].id==key:
                return self.lst.pop(i)

    def _pop_by_index(self,i,tag):
        return self.lst.pop(i)

    def set(self,id,**kwargs):
        """set the values of a given data element"""
        new = self.lst.pop(id)
        new.set(**kwargs)
        self.lst.append(new)

    def get_locked(self, iden, tag, param):
        tmp = []
       
        for mod in self.lst:
            try:
                if to_bool(mod.get(iden,tag,param+"Locked")):
                    tmp.append(mod)
            except:
                pass
        tmp =  sorted(tmp, key=lambda x: x.chi2)
        return [(mod.get(iden,tag,param), mod.chi2) for mod in tmp]

    def get_best_lst(self, iden=None, param=None):
        """
        gets best fit.

        if `param' is specified, then gets best fit with a given parameter locked.
        param should be either bLocked, zLocked or NLocked
        """
        if not param is None:
            self.lst=sorted(self.lst, key=lambda x: x.chi2)
            return [mod, mod.chi2 for mod in self.lst]
        else:
            return self.get_locked(iden, param)  #already sorted

    def get_min_chi2(self):
        return np.amin(np.array([item.chi2 for item in self.lst]))

    def best_fit(self,iden,param,order,xmin,xmax, locked=True, plot=True):
        """
        get a best fit of data with respect to `param'

        iden: id of absorber
        param:  parameter name (N,b,z)
        order:  order of polynomial to fit
        xmax, xmin: range of values to consider
        locked:  get only locked parameters?
        plot:   plot the data?  otherwise return function, x, y
        """
        x = []
        y = []
        for item in self.lst:
            ab = item.getabs(iden)
            if (locked and ab.locked(param)) or not locked:
                x.append(float(getattr(ab,param)))
                y.append(float(item.chi2))

        f = np.poly1d(np.polyfit(np.array(x),np.array(y),int(order)))
        xx = np.arange(xmin,xmax, (xmin-xmax)/(600.))
        if plot:
            plt.plot(x,y,'ro')
            plt.plot(xx,f(xx),'b-')
            plt.show()
        return f, x, y

    def get_err(self, iden, param_name):
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
                onesig.append(getattr(item.getabs(iden),param_name))
        return getattr(lst[0].getabs(iden),param_name) ,max(onesig), min(onesig)
        
    def append(self, model):
        self.lst.append(model)

    def write(self,filename=None):
        if filename==None:
            filename=self.name
        root = self.create(filename)
        self.dbxml.write(filename,root)

    def get(self,xmlfile,chi2,pixels):
        self.lst.append(Model(xmlfile=xmlfile,chi2=chi2,pixels=pixels))

    def get_model(self, iden):
        for item in self.lst:
            if item.id==iden: 
                return item

    def get_vel_shift(self,iden1,iden2):
        return [item.get_vel(iden1,iden2) for item in self.lst]

    def pop(self,i):
        return self.lst.pop(i)
    
    @staticmethod
    def read(filename):
        """read from xml, return inputs for Model"""
        root=xmlutils.Model_xml.get_root(filename)
        models = root.findall('model')
        if len(models)==0:
            raise Exception("no models saved")
        for model in models:
            lst = []
            kwargs={}
            for key, val in model_classes.items():
                tmplst=[]
                for item in model.findall(val):
                    assert(item.tag==val)
                    tmplst.append(data_types.Data.factory(node=item,tag=item.tag))
                if len(tmplst)>0:
                    kwargs[key]=tmplst
            for key, val in model.attrib.items():
                kwargs[key] = val
            lst.append(Model(**kwargs))
        return ModelDB(filename, models=lst)

    def parse_kwargs(self,**kwargs):
        for key, val in kwargs.items():
            try:
                setattr(self,key,float(val))
            except:
                setattr(self,key,val)

    def grab(self,xmlfile,chi2,pixels,**kwargs):
        """grab from xml file"""
        #need to reinstantiate xml file
        mod = Model(xmlfile=xmlfile,chi2=chi2,pixels=pixels,**kwargs)
        self.lst.append(mod)
        return

    def set(self,id,**kwargs):
        mod = self.get_model(id)

    def merge(self,other):
        self.lst += other.lst

    def create(self,filename):
        """create an xml file file structure.  returns root"""
        import xml.etree.ElementTree as et
        import datetime
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
        for item in self.lst:
            current_group = None
            group_name = item.id 
            if current_group is None or group_name != current_group.text:
                current_group = et.SubElement(root, 'model', {'id':group_name,'xmlfile':item.xmlfile,'chi2':str(item.chi2),'pixels':str(item.pixels)})

            children = []

            for att in model_classes.keys():
                if getattr(item,att) != None:
                    children+=[it.node for it in getattr(item,att)]
            if len(children)==0:
                raise Exception("no children are present")
            current_group.extend(children)
        return root


#some helper functions


def read_in(name,return_db=True):
    """
    read in a model database from text file
    Inputs:
    ------
    name: filename of model database (default = model_database.txt)

    returns:
    --------
    ModelDB instance
    """

    assert(type(name) is str)
    try:
        f = open(name,'r')
    except:
        raise Exception(str(name))
    models = []
    iden=0
    inp = f.readlines()
    assert(len(inp)>1)
    xmlfile = str(inp[0].strip())
    if ".xml" not in xmlfile:
        msg = "relevant dude .xml file must be notated as line 1 of %s\ngot %s instead"%(name,xmlfile)
        print(msg)
        xmlfile = input("enter relevant file:  ")
    i=1
    while i<len(inp):
        temp = []
        while i<len(inp):  #yes this looks weird, but is correct.
            if inp[i].strip()!='\n': 
                temp.append(inp[i])
                if 'chi2' in inp[i]: 
                    i+=1
                    break
                else:
                    i+=1
        if i>=len(inp):
            break
        newmod = parse_single_model(xmlfile, temp, iden=str(iden))
        models.append(newmod) 
        iden+=1
        i+=1
    if return_db:
        return ModelDB(name, models=models) 
    else:
        return models 



def parse_abs(xmlfile,data):
    """
    parse an individual absorber, either from a list of arguments or a raw string with single absorber
    """
    data = re.sub('=\s+','=', data.strip()).split()
    dct=dict([item.split('=') for item in data])
    dct["xmlfile"] = xmlfile
    return data_types.Absorber(**dct)

def parse_single_model(xmlfile, lines, iden=None):   
    """parse a string representation of a single model.

    an example model is:
    iden=HI    N=17.12345 b=12.345678 z=1.234567890
    iden=SiI   N=11.12345 b=12.345678 z=1.234567890
    iden=OI    N=11.12345 b=12.345678 z=1.234567890
    locked=OI:bLocked HI:zLocked
    chi2=123.4 pixels=154

    different models are differentiated by empty lines.

    Inputs:
    -------
    lines:  a list of a few lines of text

    returns:
    --------
    Model instance
    """

    try:
        str(xmlfile)
    except:
        xmlfile = str(xmlfile.name)

    absorbers = []
    for line in lines:
        if line[0:4]=='iden':
            ab = parse_abs(xmlfile, line)
            absorbers.append(ab)
        elif line[0:4]=='lock':
            lst=(line.split('=')[1]).split()
            dic = dict([item.split(':') for item in lst])
            for key, val in dic.items():
                for item in absorbers:
                    if item.id==key.strip():
                        setattr(item,val.strip(),True)
        elif 'chi2' in line:  #this should always be the last line
            lst=line.split()
            kw=dict([(item.strip()).split('=') for item in lst])
            return Model(absorbers=absorbers, iden=iden, **kw)
        else: 
            pass
    raise Exception("Input error for model database")

def to_bool(string):
    string = string.lower()
    if string=="true":
        return True
    else:
        return False
    

def inv_dict(tag, dic=model_classes):
    if tag in dic.values():
        tmp = {v:k for k, v in dic.items()} 
        return tmp[tag]

