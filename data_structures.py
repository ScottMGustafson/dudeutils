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
            string+="iden=%6s N=%8.5lf b=%8.5lf z=%lf\n"%(item.id,item.N,item.b,item.z)
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

    def get_model(self,xmlfile=None):
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

    def get(self,iden,param):
        try:
            return float(getattr(self.getabs(iden),param))
        except:
            return getattr(self.getabs(iden),param)

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

    def getabs(self, iden):
        for item in self.absorbers:
            if iden==item.id:
                return item

    def get_vel(self,iden1,iden2):
        z1 = float(self.getabs(iden=iden1).z)
        z2 = float(self.getabs(iden=iden2).z)
        return astro.get_vel_shift(z1,z2)
       
    def parse_node(self,node):
        dat = self.xml.get_node_data(node=node)
        parse_kwargs(self,**dat)

class ModelDB(object):
    def __init__(self, name=None, models=None, constraints=None,**kwargs): 
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

        if models:
            self.lst = models
            self.root=self.create(str(name))
        else:   
            self.lst = read_in(str(name), return_db=False)
            self.root=self.dbxml.read(filename)

        if constraints:
            self.lst = [item for item in self.lst if item.constrain(constraints)]

        self.dbxml=xmlutils.Model_xml()
        self.name = self.dbxml.filename if name==None else name

    def get_locked(self, iden, param):
        tmp = []
        for mod in self.lst:
            if getattr(mod.getabs(iden),param):
                tmp.append(mod)
        return sorted(tmp, key=lambda x: x.chi2)

    def get_best_lst(self, iden=None, param=None):
        """
        gets best fit.

        if `param' is specified, then gets best fit with a given parameter locked.
        param should be either bLocked, zLocked or NLocked
        """
        if not param is None:
            self.lst=sorted(self.lst, key=lambda x: x.chi2)
            return self.lst
        else:
            return self.get_locked(iden, param)

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
        lst = self.get_best_lst(param=param_name+'Locked')
        chi2min = float(lst[0].chi2)
        onesig=[]
        for item in lst:
            if item.chi2<chi2min+1.:
                onesig.append(getattr(item.getabs(iden),param_name))
        return getattr(lst[0].getabs(iden),param_name) ,max(onesig), min(onesig)
        
    def append(self, model):
        self.lst.append(self.model)

    def write(self,filename=None):
        if filename==None:
            filename=self.name
        root = self.create(filename)
        self.dbxml.write(filename,root)

    def write_to_db(self, clobber=False, name=None):
        import os.path
        if name is None:
            name=self.name
        if os.path.isfile(name) and not clobber:
            answer = input('ok to clobber? '+name+' y/n')
            if answer=='n':
                return
        f = open(name,'w')

        f.write(str(self.xmlfile)+'\n')
        for item in self.lst:
            f.write(str(item)+'\n')
        f.close()
        self.dbxml.write(self)  #should create xmldb in __init__

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
            kwargs["xmlfile"]=model.get("xmlfile")
            lst.append(Model(**kwargs))
        return ModelDB(filename, models=lst)

    def parse_kwargs(self,**kwargs):
        for key, val in kwargs.items():
            try:
                setattr(self,key,float(val))
            except:
                setattr(self,key,val)

    def set(self,id,tag,**kwargs):
        """set the values of a given data element"""
        for item in getattr(self,tag):
            if item.id==id:
                item.set_node(**kwargs)

    def grab(self,**kwargs):
        """grab from xml file"""
        #need to reinstantiate xml file
        for key, val in model_classes.items():
            lst=[]
            for item in self.xml.getDataList(val):
                lst.append(data_types.Data.factory(xmlfile=self.xml.name,node=item))
            kwargs[key] = lst

        kwargs["xmlfile"]=self.xml.name
             
        if kwargs.get("chi2",True):
            warnings.warn("no chi2 defined")
        return Model(**kwargs)

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
                current_group = et.SubElement(root, 'model', {'id':group_name,'xmlfile':item.xmlfile})

            children = []

            for att in model_classes.keys():
                if getattr(item,att) != None:
                    children+=[it.node for it in getattr(item,att)]
            if len(children)==0:
                raise Exception("no children are present")
            current_group.extend(children)
        return root


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
        print(str(newmod))
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

