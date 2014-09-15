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
import xml.etree.ElementTree as et
import dude_xmlutils
import astronomy_utils as astro
import warnings

class Model(object):
    def __init__(self,absorbers, **kwargs):
        """
        inputs:
        -------
        absorbers: list(dude_xmlutils.Absorber)  is a list of absorber references 
        """
        self.iden = kwargs.get('iden',None)
        self.ref = absorbers[0]
        self.absorbers = absorbers
        self.chisq = float(kwargs.get('chisq',0.))
        self.pixels= kwargs.get('pixels',0.)
        self.locked = {}

    def constrain(self, constraints):
        """
        Inputs:
        -------
        constraints: a dict of constraints using Absorber attirbutes.  Each key-value pair is a string-tuple of floats.
        

        returns: boolean

        example:

        >>>ab=dude_xmlutils.Absorber(N=12.2,b=10.,z=0.)
        >>>model=Model([ab])
        >>>model.constrain({ID:{N:(12,13), b:(9,11), z:(-1,1)}})
            True

        """
        for item in self.absorbers:
            try:
                for key, val in constraints[item.iden].items():
                    if getattr(item,key)<val[0] or getattr(item,key)>val[1]:
                        return False
            except KeyError:
                pass
        return True

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
            string+="iden=%6s N=%8.5lf b=%8.5lf z=%10.8lf\n"%(item.iden,item.N,item.b,item.z)
            for param in ['NLocked','bLocked','zLocked']:
                if getattr(item,param):
                    locked[item.iden] = param
        string+="locked="
        for key, val in locked.items():
            if val:
                string+=key+":"+val+" "
        string+="\nchi2=%lf pixels=%lf\n\n"%(float(self.chisq),float(self.pixels))
        return string
    def write(self):
        for item in self.absorbers:
            item.writeData()
    def getabs(self, iden):
        for item in self.absorbers:
            if iden==item.iden:
                return item


class ModelDB(object):
    def __init__(self,lst_of_Models=None, constraints=None, name='model_database.txt'): 
        """
        Model Database

        Inputs:
        -------
        lst_of_Models:  list of Model instances
        constraints: dict of dicts of tuples of floats.  (see Model.constrain) 


        """

        if lst_of_Models is None:
            lst_of_Models = read_in(name, return_db=False)
            
        if constraints:
            self.lst = [item for item in lst_of_Models if item.constrain(constraints)]
        else:
            self.lst = lst_of_Models
        self.filename=name

    def get_locked(self, iden, param):
        tmp = []
        for mod in self.lst:
            if getattr(mod.getabs(iden),param):
                tmp.append(mod)
        return sorted(tmp, key=lambda x: x.chisq)

    def best_fit(self, iden=None, param=None):
        """
        gets best fit.

        if `param' is specified, then gets best fit with a given parameter locked.
        param should be either bLocked, zLocked or NLocked
        """
        if not param is None:
            self.lst=sorted(self.lst, key=lambda x: x.chisq)
            return self.lst
        else:
            return self.get_locked(iden, param)
        
    def get_err(self, iden, param_name):
        """get 1 sigma error from chi2 = chi2min + 1

        param_name is in [N,b,z]

        """
        lst = self.best_fit(param=param_name+'Locked')
        chisqmin = float(lst[0].chisq)
        onesig=[]
        for item in lst:
            if item.chisq<chisqmin+1.:
                onesig.append(getattr(item.getabs(iden),param_name))
        return getattr(lst[0].getabs(iden),param_name) ,max(onesig), min(onesig)
        

    def append(self, model):
        self.lst.append(self.model)
        f = open(self.filename,'a')
        f.write(str(model))
        f.close()

    def write_to_db(self, clobber=False, name=None):
        import os.path
        if name is None:
            name=self.filename
        if os.path.isfile(name) and not clobber:
            answer = input('ok to clobber? '+name+' y/n')
            if answer=='n':
                return
        f = open(name,'w')
        for item in self.lst:
            f.write(str(item)+'\n')
        f.close()

    def get_model(self, iden):
        for item in self.lst:
            if item.iden==iden: 
                return item

def read_in(name='model_database.txt',return_db=False):
    """
    read in a model database from text file
    Inputs:
    ------
    name: filename of model database (default = model_database.txt)


    returns:
    --------
    ModelDB instance
    """

    f = open(name,'r')
    models = []
    iden=0
    inp = f.readlines()
    xmlfile = str(inp[0].strip())
    for i in range(1,len(inp)):
        temp = []
        while inp[i] != '\n':
            temp.append(inp[i])
            i+=1  
        models.append(parse_line(xmlfile,temp, iden=str(iden))) 
        iden+=1
    if return_db:
        return ModelDB(models) 
    else:
        return models 

def parse_abs(xml_file,data):
    """
    parse an individual absorber, either from a list of arguments or a raw string with single absorber

    """

    if type(data) is str:
        data = [item.strip() for item in data.split()]
    data=dict([item.split('=') for item in data])
    data["xmlfile"] = xml_file
    return dude_xmlutils.Absorber(**data)

def parse_line(xml_file,lines, iden=None):   
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

    absorbers = []
    for line in lines:
        if line[0:2]=='iden':
            absorbers.append(parse_abs(xml_file,line))
        if line[0:4]=='lock':
            line=(line.split('=')[1]).split()
            line = dict([item.split(':') for item in line])
            for key, val in line.items():
                for item in absorbers:
                    if item.iden==key.strip():
                        setattr(item,val.strip(),True)
        if line[0:4]=='chi2':  #this should always be the last line
            line=line.split()
            kw=dict([(item.strip()).split('=') for item in line])
            return Model(absorbers, iden=iden, **kw)
    raise Exception("Input error for model database")
           
def parse_old_format(input_file, xmlfile):
    """ used only for translating an old database.  should have 16 columns and really should only be run once ever by me"""
    f = open(input_file)
    count=0
    mods = []
    for line in f:
        try:
            assert(len(line.split())==16)
        except:
            if '.xml' in line:
                pass
            else:
                warnings.warn("this model is not the same as the rest: %d"%(count))
            continue
        id1,N1,b1,z1,id2,N2,b2,z2,v2,id3,N3,b3,z3,v3,_,chi2 = tuple(line.split())
        abs1 = dude_xmlutils.Absorber(iden=id1,N=N1,b=b1,z=z1,xmlfile=xmlfile)
        abs2 = dude_xmlutils.Absorber(iden=id2,N=N2,b=b2,z=z2,xmlfile=xmlfile)
        abs3 = dude_xmlutils.Absorber(iden=id3,N=N3,b=b3,z=z3,xmlfile=xmlfile)
        mods.append(Model([abs1,abs2,abs3], chisq=chi2,iden=str(count)))
        count+=1
    return ModelDB(mods,name=input_file+".new")
    
        

def parser(arg_list):
    """
    parses input arguments.
    """

    params = ['N','b','z']
    datalst = []
    ablst = []
    db = ModelDB()
    i=0
    
    while i<len(arg_list):
        item=arg_list[i]
        if item in ['write','get','get_best','get_err']:
            action = item
        elif '.xml' in item:
            datfile=item
        elif item in params:
            param=item
        elif 'id' in item or 'iden' in item:
            iden=item.split('=')[1].strip()
            if action=='get':  #write current model into database
                ablst.append(dude_xmlutils.Absorber(iden=ion, xmlfile=datfile))
            elif action=='write': #write commandline args into xml file
                absorber=dude_xmlutils.Absorber(iden=iden, xmlfile=datfile, populate=True)
                i+=1
                vals = {}
                while arg_list[i][0] in params:
                    item=arg_list[i].split('=')
                    vals[item[0]] = vals[item[-1]] 
                    i+=1
                absorber.writeData(**vals)
            abslst.append(absorber)
        elif 'chisq' in item:
            chisq=float(item.split('=')[-1])
        elif 'pixels' in item:
            pixels=float(item.split('=')[-1])
        else:
            raise Exception('unrecognized argument '+item)

    model = Model(abslst,chisq=chisq,pixels=pixels)

    if action=='write':
        model.write()    
    elif action=='get':
        db.append(model)
    elif action=='get_best':
        print(str(db.bestfit()))
    elif action=='get_err':
        best,mx,mn = db.get_err(iden,param)
        print("%s:  %s= %lf, %lf, %lf"%(iden,param,best,mx,mn))
        

if __name__ == '__main__':

    """
    input like:
    >>>python cmd_line_io.py write test.xml id=HI N=12 b=4 id=H2 N=17.8 z=3.2 
    """
    parser(sys.argv[1:])
