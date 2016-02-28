"""
Program to randomize continua.

"""

from model import * 
from dudeutils import run_optimize
import numpy as np


def get_lims(val,percent):
    """
    
    input:
    ------
    val: (float) initial value of continuum control point
    percent: percentage by which to vary continuum point

    output:
    -------
    list: [min, max] value

    raises:
    ------
    None

    """  
    return [(1.-percent)*val, (1.+percent)*val]

def construct_lims(model,fname,xpercent=None,ypercent=None):
    """
    
    input:
    ------
    model:  model.Model instance
    fname: (str) filename
    xpercent, ypercent: percentage by which to vary continuum control point

    output:
    -------
    dct (dict): dict of limits for each control point

    raises:
    -------
    None

    """  
    dct={}
    if fname=='':
        cont_lst=Model.get(mod.ContinuumPointList)
        pts=[item for item in cont_lst if item.id!='null']
        dct={}
        for item in pts:
            dct[item.id]={'ylim':ypercent,'xlim':xpercent}
        return dct
    
    with open(fname) as f:
        for line in f.readlines():
            try:
                key, xpercent, ypercent=line.split()
            except:
                key, ypercent=line.split()
                xpercent=0.
            dct[key]={'ylim':ypercent,'xlim':xpercent}
    return dct

def randomize_point(model,iden,ymax,ymin,xmax,xmin,gaussian=False):
    """
    
    input:
    ------
    model:  model.Model instance
    iden (str): id of continuum point
    ymax, ymin,xmax,xmin: limits for continuum point
    dct: dict of limits for continuum points
    gaussian (bool): if true, points are distributed along a gaussian distribution 

    output:
    -------
    None

    raises:
    ------
    None

    """  
    if ymax and not ymin==ymax:
        model.monte_carlo_set(iden, 'ContinuumPoint',[ymin,ymax],'y',gaussian)
    if xmax and not xmin==xmax:
        model.monte_carlo_set(iden, 'ContinuumPoint',[xmin,xmax],'x',gaussian)

def randomize_list_of_points(model,dct,gaussian=False):
    """
    
    input:
    ------
    model:  model.Model instance
    dct: dict of limits for continuum points
    gaussian (bool): if true, points are distributed along a gaussian distribution

    output: 
    -------
    None

    raises:
    -------
    None

    """  
    for key, val in dct.items():
        if type(val['xlim']) in [list,tuple]:
            xmin,xmax = val['xlim'][0],val['xlim'][-1]
        else:
            xmin,xmax = get_lims(model.get_datum(key,'ContinuumPoint','x'), val['xlim'])
        if type(val['ylim']) in [list,tuple]:
            ymin,ymax = val['ylim'][0],val['ylim'][-1]
        else:
            ymin,ymax = get_lims(model.get_datum(key,'ContinuumPoint','y'), val['ylim'])
        randomize_point(model,key,ymin,ymax,xmin,xmax,gaussian)

def get_continua(mod,n,limits_file='',gaussian=False,xpercent=0.,ypercent=0.02): 
    """
    
    input:
    ------
    mod:  (model.Model) model from which we get our continuum points
    n: (int) how many times to randomize continua
    limits_file: (string) name of file providing limits
    gaussian: (bool) if True, have point follow a gaussian distribution from the mean.
    xpercent: percentage by which to vary x componenents 
    ypercent: percentage by which to vary y componenets 

    output:
    -------
    list of models

    raises:
    ------
    None

    """     
    dct=construct_lims(mod,limits_file,xpercent,ypercent)
    models=[]
    path=os.path.split(mod.xmlfile)[0]
    for i in range(n):
        randomize_list_of_points(mod,dct,gaussian)
        fname=path+"/continuum_%d.xml"%(i)
        mod.write(fname)
        run_optimize(fname)
        models.append(Model(xmlfile=fname))
    return models
    
if __name__=='__main__':
    name="/home/scott/research/J0744+2059/2016-02-16.xml"
    n=5
    limits_file="data/cont_limits.dat"
    mod=Model(xmlfile=name)
    lst = mod.get_lst("ContinuumPointList")

    try:
        for item in lst:
            assert(item.xLocked and item.yLocked)
    except:
        raise Exception("need to set all continuumpoints to locked for %s\n"%(name))
    models = ModelDB(models=get_continua(mod,n,limits_file='',gaussian=True,xpercent=0.,ypercent=0.02))
    print(len(models))
    
