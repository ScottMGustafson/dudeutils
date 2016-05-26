from dudeutils.utilities import *
from dudeutils.model import *
import sys, os
from dudeutils.plot_distribution import plot_chi2
from configparser import ConfigParser
from scipy.constants import c
import numpy as np
from subprocess import TimeoutExpired
from dudeutils.get_data import get_data
c*=0.001   #convert to km/s

class ModelFailed(BaseException):
    def __init__(self,msg):
        self.msg=msg
    def __str__(self):
        return self.msg

def parse_config(config_file):
    """
    parse the config file to produce a dict of absorbers, parameters and 
    allowed ranges

    Input:
    ------
    config_file: config file to use.  data/random_sampling_config.cfg by default
    

    Output:
    -------
    dict of parameters and allowed ranges

    Raises:
    -------
    None

    """
    def set_config_defaults(dct):
        
        for key in "append step method gaussian vary_continuum to_buffer".split():
            if not key in list(dct['config'].keys()):
                dct['config'][key]=False

        if not "n" in list(dct['config'].keys()):
            dct['config']["n"]=2

        if not "method" in list(dct['config'].keys()):
            dct['config']["method"]="dude"
        

    def tf(val):
        if val.lower()=='true':
            return True
        elif val.lower()=='false':
            return False
        else:
            return val
    
    config = ConfigParser()
    config.optionxform=str
    config.read(config_file)
    dct={}
    for item in list(config.sections()):
        dct[item]=dict(config.items(item))

    for ab in dct.keys():
        if ab=='config':
            for key in dct[ab].keys():
                dct[ab][key]=tf(dct[ab][key])
                
        if ab=='continuum':
            for key, val in dct[ab].items():
                val=val.split(',')
                if len(val)==2:
                    dct[ab][key]={'xlim':float(val[0]),'ylim':float(val[1])}
                dct[ab][key] = {'ylim':float(val[0])}
        else:
            for k in ["N", "b", "z"]:
                try:
                    dct[ab][k]=dct[ab][k].replace(" ","")
                    dct[ab][k] = list(map(float,dct[ab][k].strip().split(',')))
                except KeyError:
                    pass

            cond=[]
            if "N" in dct[ab].keys():
                cond+=[ dct[ab]["N"][0]<10.,dct[ab]["N"][-1]>23.,
                        dct[ab]["N"][-1]<dct[ab]["N"][0]]

            if "b" in dct[ab].keys():
                cond+=[ dct[ab]["b"][0]<0.,dct[ab]["b"][-1]>50.,
                        dct[ab]["b"][-1]<dct[ab]["b"][0]]

            if "z" in dct[ab].keys():
                cond+=[ dct[ab]["z"][0]<0.,dct[ab]["z"][-1]<0.,
                        dct[ab]["z"][-1]<dct[ab]["z"][0]]

            
            if any(cond):
                raise ModelFailed("check your random sampling inputs")

    set_config_defaults(dct)
    return dct

#toggle lock controls  idens, params, locked,tag='Absorber'
def toggle_ab_locks(model, ab_cfg, locked,locked_dct=None):

    dct={}
    for key, val in ab_cfg.items():
        for attr in val.keys():
            if key in dct.keys():
                dct[key].append(attr)
            else:
                dct[key] =[attr]
    model.toggle_locks(dct,locked,'Absorber')
    if locked_dct:
        model.toggle_locks(locked_dct, True, 'Absorber')


def toggle_cont_ab_locks(model, ab_lock,ab_cfg,cont_cfg, 
                         cont_lock=None, locked_abs=None, locked_cont=None):

    if not cont_lock:
        cont_lock=not ab_lock
    toggle_ab_locks(model,ab_cfg,ab_lock,locked_dct=locked_abs)
    
    if locked_cont:
        cont_ids=[it for it in cont_cfg.keys() if not it in locked_cont.keys()]
    else:
        cont_ids=list(cont_cfg.keys())
    model.toggle_cont_lock(cont_ids,'y',locked=cont_lock)

def sort_by_z(model, idens):
    """
    Sometimes when optimizing a model, two absorbers close in velocity will 
    cross over each other in redshift.   This rearaanges the ids to the 
    appropriate absorbers.

    input:
    ------
    model: model.Model instance
    idens: list of absorber ids to consider, sorted by expected z.

    output:
    -------
    None

    """
    ab_lst = [Model.get_datum(iden,'Absorber') for iden in idens]
    aborbers=sorted([ab for item in ab_lst], key=lambda x: x.z ) 
    if idens!=[ab.id for item in absorbers]:
        for i in range(len(idens)):
            model.set_val(absorber[i].id, tag='Absorber', **{"id":idens[i]})



def perturb_absorbers(dct, model, **kwargs):
    """
    randomly perturb absorbers to explore the parameter-space

    Input:
    ------
    dct: (dict) a dict of parameters and ranges, output from parse_config
    model: a model.Model instance

    Output:
    -------
    a model.Model instance with perturbed absorbers

    Raises:
    -------
    None

    """

    for key, val in dct.items():
        if key in ['continuum', 'config']: continue
        if not type(val) is dict:
            raise Exception(str(val)+" should be a dict")
        for param, rng in val.items():
            model.monte_carlo_set(key,"Absorber",[ rng[0],rng[-1] ],param,kwargs.get('gaussian',True)) 


def perturb_values(model,ab_cfg,**kwargs):
    perturb_absorbers(ab_cfg, model, **kwargs)  
    cont_cfg=kwargs.get('cont_cfg',None)
    if cont_cfg:
        perturb_continua(model, cont_cfg, **kwargs)
    model.write(model.xmlfile)

def is_better(new_chi2,old_chi2,tol=4.): return new_chi2<old_chi2+tol

#the varying methods

def vary_ab(model,ab_cfg,locked_dct=None, **kwargs):
    """
    varies param values.

    params:
    -------
    model: Model instance
    ab_cfg: absorber config dict
    n: (int) number of iterations
    tol:  chi2 tolerance: model accepted if better than within prev_chi2+tol
    locked_keys,locked_params: list of keys, params to keep locked

    returns:
    --------
    list of new models

    raises:
    -------
    None
    """
    perturb_absorbers(ab_cfg, model, gaussian=True)
    toggle_ab_locks(model, ab_cfg, False,
                    locked_dct=locked_dct)
    model.write()
    if kwargs.get('step',False):
        buff=run_optimize(model.xmlfile,to_buffer=kwargs['to_buffer'],step=True,timeout=30)
        return populate_database(abs_ids=list(ab_cfg.keys()), return_list=True,
                            path=os.path.split(model.xmlfile)[0],
                            buff=buff)
    else:
        run_optimize(model.xmlfile,timeout=30,**kwargs)
        model.read()          
        return [model.copy()]

def vary_cont(model,ab_cfg,locked_abs=None, locked_cont=None, **kwargs):
    """
    varies param values AND continuum level at specified points.

    params:
    -------
    model: Model instance
    ab_cfg: absorber config dict
    cont_cfg: continuum config dict
    tol:  chi2 tolerance: model accepted if better than within prev_chi2+tol
    locked_keys,locked_params: list of keys, params to keep locked

    returns:
    --------
    list of new models

    raises:
    -------
    None
    """

    cont_cfg = kwargs.get('cont_cfg')
    def find_better_cont():
        """
        one iteration of our continuum varying procedure:
        toggle between optimizing continua and absorbers until find min chi2
        """
        toggle_cont_ab_locks(model,False,ab_cfg,cont_cfg,
                             locked_abs=locked_abs, locked_cont=locked_cont)
        model.write()
        run_optimize(model.xmlfile,timeout=30,**kwargs)
        model.read()
        toggle_cont_ab_locks(model,True,ab_cfg,cont_cfg,
                            locked_abs=locked_abs, locked_cont=locked_cont)
        model.write()

        return run_optimize(model.xmlfile,timeout=30,**kwargs)


    if kwargs.get('step',False):
        perturb_values(model,ab_cfg,**kwargs)
        buff=find_better_cont()
        return populate_database(abs_ids=list(ab_cfg.keys()),
                            path=os.path.split(model.xmlfile)[0],
                            buff=buff,return_list=True)
        
    else:
        tol=kwargs.get('tol',4.)
        old_model=model.copy()
        prev, new_chi2=old_model.chi2, 0.
        perturb_values(model,ab_cfg,**kwargs) #initial perturbation  
        count=0
        while np.fabs(new_chi2-prev)>tol and count<3:
            #while chi2 hasn't improved or optimization hasn't stalled
            #alternate between varying continua and absorbers
            prev=model.chi2
            find_better_cont()
            model.read()
            new_chi2=model.chi2
            if not is_better(new_chi2, prev,tol):
                count+=1
            else:
                count=0
        return model.copy()



def plt_sigmas(db,ax,x_attr,y_attr,xargs=[],yargs=[],**kwargs):
    min_chi2=min([item.chi2 for item in db])
    xfact=kwargs.get("xfact",1.)
    yfact=kwargs.get("yfact",1.)
    def one_sig(mod): return mod.chi2<=min_chi2+1.
    def two_sig(mod): return 1.<mod.chi2<=min_chi2+4.
    def morethan_2sig(mod): return 4.<mod.chi2

    ax.plot(  xfact*np.array(db.get_attr_lst(x_attr,morethan_2sig,*xargs)),
              yfact*np.array(db.get_attr_lst(y_attr,morethan_2sig,*yargs)),
              'ko'
           )

    ax.plot(  xfact*np.array(db.get_attr_lst(x_attr,two_sig,*xargs)),
              yfact*np.array(db.get_attr_lst(y_attr,two_sig,*yargs)),
              'bo'
           )

    ax.plot(  xfact*np.array(db.get_attr_lst(x_attr,one_sig,*xargs)),
              yfact*np.array(db.get_attr_lst(y_attr,one_sig,*yargs)),
              'co'
           )  

def perturb_continua(mod, dct, **kwargs): 
    """
    
    input:
    ------
    mod:  (model.Model) model from which we get our continuum points
    output:
    -------
    list of models

    raises:
    ------
    None

    """     
    for key, val in dct.items():
        y=mod.get_datum(key,"ContinuumPoint",'y')
        rng= (1.-float(val['ylim']))*y, (1.+float(val['ylim']))*y
        try:
            rng= list(map(float, rng))
        except:
            raise TypeError("error parsing continuum points: rng="+str(rng))
        mod.monte_carlo_set(key,"ContinuumPoint",rng,
                            'y',gaussian=kwargs.get('gaussian',True))
        if 'xlim' in val.keys():
            x=mod.get_datum(key,"ContinuumPoint",'x')
            rng=(1.+float(val['xlim']))*x, (1.-float(val['xlim']))*x
            mod.monte_carlo_set(key,"ContinuumPoint",rng,
                                'x',gaussian=kwargs.get('gaussian',True))


def filter_bad_models(models, dct, vel_pad=10.0,chi2pad=150.):
    if len(models)==0:
        return models
    min_chi2=min([float(item.chi2) for item in models])

    def _filter(model):

        for cnt in model.get(model.ContinuumPointList):
            for param in ["x","y"]:
                if str(getattr(cnt,param)).lower()=="nan":
                    return False

        if chi2pad:
            if float(model.chi2)>min_chi2+chi2pad or float(model.chi2)==0.:
                return False
        for iden, params in dct.items():
            if iden=='continuum':continue
            for param_name, param_range in params.items():
                val=model.get_datum(iden,"Absorber",param_name)
                if val==-1. or str(val).lower()=="nan":
                    return False 
                if param_name == "b":
                    if val<0.1:
                        return False
                    elif model.get_datum(iden,"Absorber","ionName")!="H I":
                        if val>60:   #a metal line with an extremely large 
                                     #b-value: not likely.
                            return False
                elif param_name == "N":
                    
                    if val<7.0 or val>23.:  #this means dude tried to throw it out.
                                  #if you are running this, then you've already decided that 
                                  #this absorber was needed, so will this is a bad model
                        return False
                    elif val>param_range[-1]+2 or val<param_range[0]-2:
                        return False 
                        #if model goes more than about 100 times bigger or smaller 
                        #outside your specified range
                elif param_name == "z":
                    vel_range=[c*(val-param_range[0])/(1.+param_range[-1]),
                               c*(val-param_range[-1])/(1.+param_range[-1])]
                    if (val>param_range[-1] and vel_range[-1]>vel_pad) or \
                       (val<param_range[0] and vel_range[0]<-1.*vel_pad):
                        return False
        return True

    if type(models) is ModelDB:
        for item in models.models:
            if not _filter(item):
                models.remove(item) 
        models.remove_unused()
        return models        
    else:
        return [item for item in models if _filter(item)]   

def plt_N_vs_cont(all_abs,ab_cfg, cont_cfg):    
    for cont_key in cont_cfg.keys():
        for ab in ab_cfg.keys():
            for attr in ['N']:
                x=[mod.get_datum(ab,tag="Absorber", param=attr) 
                    for mod in all_abs]
                y=[mod.get_datum(cont_key,tag="ContinuumPoint", param="y") 
                    for mod in all_abs]
                plt.plot(x,y,'ko')
                x_label=all_abs[0].get_datum(cont_key,tag="ContinuumPoint", param="x")
                plt.xlabel(r"log N(%s)"%(ab))
                plt.ylabel(r"cont level at $\lambda=$\ %5.1lf"%(x_label))
                plt.show()

def iterate_fn(fn,ab_cfg,n,kwargs):
    out_lst=[]
    for i in range(n):
        sys.stdout.write('iteration %d'%(i))
        sys.stdout.flush()        
        for key in ab_cfg.keys():
            for attr in ab_cfg[key].keys():
                if kwargs.get('_rawtxt',False):
                    model=kwargs["model"]
                    with open(model.xmlfile, 'w') as f:
                        f.writelines(kwargs.get('_rawtxt'))
                    model.read()
                    kwargs["model"]=model
                if kwargs['vary_continuum']:
                    kwargs['locked_keys']=list(ab_cfg[key].keys())
                    kwargs['locked_params']=['z']
                else:
                    kwargs['locked_keys']=[key]
                    kwargs['locked_params']=[attr]
                
                #a quick hack to ensure z cannot vary freely...ever.
                if not 'z' in kwargs['locked_params']:
                    kwargs['locked_params'].append('z')
                try:
                    out=fn(**kwargs)
                    out_lst+=out if type(out) is list else [out]
                    sys.stdout.write('.')
                    sys.stdout.flush()
                except TimeoutExpired:
                    sys.stdout.write('_')
                    sys.stdout.flush()
                    continue
                except KeyboardInterrupt:
                    return out_lst
                except:
                    raise
                
        sys.stdout.write('\n')
        sys.stdout.flush()
        out_lst=filter_bad_models(out_lst, ab_cfg)
    return out_lst

def random_sampling(db, model, ab_cfg, cont_cfg, **kwargs):
    """
    runs random sampling dude models to estimate errors.  This works best under 
    a small number of parameters and parameter-space to explore.  Do this when 
    you ALREADY have an idea of what the result should be.  This can then be
    used to get an estimate of the errors.

    input:
    ------
    model: (string or model.Model) string of xml file name
    iden:  (string) id of absorber intended to be varied
    param: (string) param to vary.  either N,b or z
    param_range: (list of tuple) range of value in the form [min,max].  
    n: (int) number of iterations
    abs_ids: (list of strings) which absorbers to keep in the model.
    constraints: (dict) model constraints.  see docstring for constraints.Constraint
    iden2: (string) if param==vel, then specify reference absorber id
    continuum: (bool) whether or not to vary continuum points specified in
               cont_pts.cfg as well

    output:
    -------
    db: model.ModelDB instance.  list of models

    Raises:
    -------
    AssertionError

    """

    def convert_to_vel():
        """
        if user instead specifies velocity from some other absorber, 
        convert to redshift, assuming the other is locked
        """
        zref=model.get_datum(iden2,"Absorber","z")
        try:
            assert(model.get_datum(iden2,"Absorber","zLocked"))
        except:
            model.set_val(iden2,"Absorber",zLocked='true')
            
        vel=param_range
        param_range=[(vel[0]/c)*(1.+zref)+zref, (vel[-1]/c)*(1.+zref)+zref]
        param='z'

        return param, param_range



    kwargs.update({'ab_cfg':ab_cfg,'cont_cfg':cont_cfg,'model':model})
    n=int(kwargs.pop('n',1))
    if kwargs['vary_continuum']:
        print('varying continua')
        model.lock_all_cont(tf=False)
        db.models+=iterate_fn(vary_cont,ab_cfg,n,kwargs)
    else:
        model.lock_all_cont()
        db.models+=iterate_fn(vary_ab,ab_cfg,n,kwargs)

    return db


def get_nsigma(db,n=1):
    """
    returns database of all models within n*sigma of best fitting model

    input:
    ------
    db: model.ModeDB database of models
    n: number of standard deviations.  n=1 (68% CI) by default

    output:
    -------
    db: model.ModelDB database of models containing only models within n*sigma
        of best fitting model

    Raises:
    -------
    None

    """
    def delta(): return float(n**2)
    lst=sorted(db.models, key=lambda x: x.chi2)
    return [item for item in lst if item.chi2<lst[0].chi2+delta()]

def main(config_file):
    #if not config_file:
    #    config_file=get_data('random_sampling_config.two_comp.cfg')
    ab_cfg=parse_config(config_file)
    cont_cfg=ab_cfg.pop('continuum',None)
    glob=ab_cfg.pop('config',None)
    if not glob:
        raise Exception
 
    if glob['append']:
        if type(glob['append']) is str:
            name=glob['append']
        else:
            name = input("path/name of db model to append: ")
        if name.endswith('.xml'):
            all_db=ModelDB(name=name)
        else:
            print('loading %s'%name)
            try:
                all_db=ModelDB.load_models(name)
            except:
                print('initializing empty DB')
                all_db=ModelDB(models=[])
    else:
        print('initializing new model db')
        all_db=ModelDB(models=[]) 

    

    model=Model(xmlfile=glob['source'])
    txt=open(glob['source'],'r').readlines()

    all_db=random_sampling( all_db, model, ab_cfg, cont_cfg, _rawtxt=txt, **glob)
    print("before filtering: %d"%len(all_db))
    all_db=filter_bad_models(all_db, ab_cfg, vel_pad=4.0,chi2pad=250.)

    open(glob['source'],'w').writelines(txt)   #in the case of keyboard interrupt, writes correct file back
    if len(all_db.models)==0:
        raise ModelFailed("no surviving models...\ncheck input for %s"%(config_file))

    print("after filtering: %d"%len(all_db))
    if type(glob['append']) is bool:
        name='latest.obj' 
    else:
        name=glob['append']
    name=os.path.splitext(name)[0]
    all_db.write(name+'.xml',True)
    print('writing %s'%(name+'.obj'))
    ModelDB.dump_models(all_db,name+'.obj')
    #models=get_nsigma(all_db,n=10)
    #for model in list(all_db.models):
    #    if not model in models:
    #        all_db.remove(model)



    if glob['plot']:    
        #plt_N_vs_cont(all_db,ab_cfg, cont_cfg)
        for key, val in ab_cfg.items():
            for attr in list(val.keys()):
                plot_chi2(all_db, iden=key, attr=attr,
                          xlabel=r"$%s(%s)$"%(attr,key),
                          constraints={})

if __name__=="__main__":
    main()



