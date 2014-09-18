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
from data_structures import Model, ModelDB
import warnings


           
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
