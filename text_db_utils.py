from model import Model, ModelDB

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
    id=0
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
        newmod = parse_single_model(xmlfile, temp, id=str(id))
        models.append(newmod) 
        id+=1
        i+=1
    if return_db:
        return ModelDB(name, models=models) 
    else:
        return models 



def parse_abs(xmlfile,data):
    """
    parse an individual absorber, either from a list of arguments or a raw 
    string with single absorber
    """
    data = re.sub('=\s+','=', data.strip()).split()
    dct=dict([item.split('=') for item in data])
    dct["xmlfile"] = xmlfile
    return data_types.Absorber(**dct)

def parse_single_model(xmlfile, lines, id=None):   
    """parse a string representation of a single model.

    an example model is:
    id=cont1 x=123   y=456
    id=cont2 x=2346   y=56896
    id=null  x=1234   y=4556
    id=cont2 x=2346   y=5966
    id=HI    N=17.12345 b=12.345678 z=1.234567890
    id=SiI   N=11.12345 b=12.345678 z=1.234567890
    id=OI    N=11.12345 b=12.345678 z=1.234567890
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
        if line[0:4]=='id':
            ab = parse_abs(xmlfile, line)
            absorbers.append(ab)
        elif line[0:4]=='lock':
            lst=(line.split('=')[1]).split()
            dic = dict([item.split(':') for item in lst])
            for key, val in dic.items(
                for item in absorbers:
                    if item.id==key.strip():
                        setattr(item,val.strip(),True)
        elif 'chi2' in line:  #this should always be the last line
            lst=line.split()
            kw=dict([(item.strip()).split('=') for item in lst])
            return Model(absorbers=absorbers, id=id, **kw)
        else: 
            pass
    raise Exception("Input error for model database")
