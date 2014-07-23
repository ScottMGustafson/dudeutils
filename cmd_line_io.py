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

class Model(object):
    def __init__(self,absorbers, chisq=None, pixels=None):
        """
        inputs:
        -------
        absorbers: list(dude_xmlutils.Absorber)  is a list of absorber references 
        """
        self.ref = absorbers[0]
        self.absorbers = absorbers
        self.chisq = chisq  
        self.size = pixels
        self.locked = {}

    def constrain(self, constraints):
        for item in self.absorbers:
            for key, val in constraints[item.iden].items():
                if getattr(item,key)<val[0] or getattr(item,key)>val[1]
                    return False
        return True

    def __str__(self):
        pass

class ModelDB(object):
    def __init__(self,lst_of_Models, constraints=None): 
        if constraints:
            self.lst = [item for item in lst_of_Models if item.constrain(constraints)]
        else:
            self.lst = lst_of_Models

    def best_fit(self):
        return sorted(self.lst, key=lambda x: x.chisq)[0]
        
    def get_err(self, param_name):
        """get 1 sigma error from chi2 = chi2min + 1"""
        for item in self.lst:
            if item.chisq<=

    def print_out(self, name='model_database.txt', clobber=False):
        import os.path
        if os.path.isfile(name) and not clobber:
            answer = input('ok to clobber? '+name+' y/n')
            if answer=='n':
                return
        f = open(name,'w')
        for item in self.lst:
            f.write(str(item)+'\n')
        f.close()

def read_in(name='model_database.txt'):
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
    inp = f.readlines()
    for i in range(len(inp))
        temp = []
        while inp[i] != '\n':
            temp.append(inp[i])
            i+=1
        
        models.append(parse_line(temp)) 
    return ModelDB(models)  


def parse_line(lines):   
    """parse several lines of input
    Inputs:
    -------
    lines:  a list of a few lines of text

    returns:
    --------
    Model instance


    """
            
    pass
        

        


#  constraint is like {param1:[lower_lim,upper_lim], ... }
#  constraints is like {id1:constraint1, id2:constraint2, ... }




def parser(arg_list):
    datalst = []
    i=0
    action = 'write' if 'write' in arg_list else 'get'
    while i < len(arg_list):
        if '.xml' in arg_list[i]:
            xml_file = arg_list[i]
        if "=" in arg_list[i]:
            key, val = arg_list[i].split('=')
            if key=='id':
                kwargs = {}
                kwargs['iden']=val
                i+=1
                while arg_list[i][:2]!='id' and arg_list[i][:5]!='chisq':
                    key, val = arg_list[i].split('=')
                    kwargs[key] = val
                    i+=1
                datalst.append(dude_xmlutils.Absorber(**kwargs))
            elif key=='chisq':
                chisq=int(arg_list[i].split('=')[1])
        i+=1
    return action, xml_file, datalst, chisq

action, xml_file, datalst, chisq = parser(sys.argv[1:])

#for i in range(1,len(datalst)):
#    setattr(datalst[i],'vel',astro.get_vel_shift(datalst[0].z,datalst[i].z)[0])
    
if action=='get': 
    for item in datalst:
        item.getData()
        item.getShift[datalst[0]]
    output_string = ''
    for item in datalst:
        output_string+=str(item)
    output_string+=' chisq='+str(chisq)
    f = open('absorberData.txt','a') # append read data to file
    f.write(output_string+'\n')
    f.close()

elif action=='write':
    for item in datalst:
        item.writeData()
    tree.write(xml_file)

else:
    raise Exception('unrecognized command: '+str(action))
