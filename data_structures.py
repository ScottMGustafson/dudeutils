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
import dude_xmlutils
import warnings
import numpy as np
import matplotlib as plt

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

    def get_best_lst(self, iden=None, param=None):
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

    def best_fit(self,iden,param,order,xmin,xmax, locked=True, plot=True):
        """
        get a best fit of data with respect to `param'

        iden: id of absorber
        param:  parameter name (N,b,z)
        order:  order of polynomial to fit
        xmax, xmin: range of values to consider
        locked:  get only locked parameters?
        plot:   plot the data?  otherwise return the minimum of the function of best fit
        """
        x = []
        y = []
        for item in self.lst:
            ab = item.getabs(iden)
            if (locked and ab.locked(param)) or not locked:
                x.append(float(getattr(ab,param))
                y.append(float(item.chisq))

        f = np.poly1d(np.polyfit(np.array(x),np.array(y),int(order)))
        xx = x.arange(xmin,xmax, (xmin-xmax)/(600.))
        if plot:
            plt.plot(x,y,'ro')
            plt.plot(xx,f(xx),'b-')
            plt.show()
        return np.amin( np.array(f(xx)) )

    def get_err(self, iden, param_name):
        """get 1 sigma error from chi2 = chi2min + 1

        param_name is in [N,b,z]

        """
        lst = self.get_best_lst(param=param_name+'Locked')
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

