from dudeutils import *
from model import *
from plot_distribution import plot_chi2
from configparser import ConfigParser
from scipy.constants import c
c*=0.001   #convert to km/s

def parse_config(config_file=os.path.join('data','random_sampling_config.cfg')):
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
    config = ConfigParser()
    config.optionxform=str
    config.read(config_file)
    dct={}
    for item in list(config.sections()):
        dct[item]=dict(config.items(item))

    for ab in dct.keys():
        if ab=='config':
            continue
        for k in ["N", "b", "z"]:
            dct[ab][k]=dct[ab][k].replace(" ","")
            dct[ab][k] = list(map(float,dct[ab][k].strip().split(',')))
        cond=[dct[ab]["N"][0]<10.,dct[ab]["N"][-1]>23.,
              dct[ab]["N"][-1]<dct[ab]["N"][0],
              dct[ab]["b"][0]<1.,dct[ab]["b"][-1]>50.,
              dct[ab]["b"][-1]<dct[ab]["b"][0],
              dct[ab]["z"][0]<0.,dct[ab]["z"][-1]<0.,
              dct[ab]["z"][-1]<dct[ab]["z"][0]]
        if any(cond):
            raise Exception("check your random sampling inputs")
    return dct

def perturb_absorbers(dct, model):
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
        for param, rng in val.items():
            
            model.monte_carlo_set(key,"Absorber",[ rng[0],rng[-1] ],param,False)
    model.write(model.xmlfile)



def filter_bad_models(models, dct, vel_pad=2.,chi2pad=50.):
    min_chi2=min([float(item.chi2) for item in models])
    def _filter(model):
        if float(model.chi2)>min_chi2+chi2pad:
            return False
        for iden, params in dct.items():
            for param_name, param_range in params.items():
                val=model.get_datum(iden,"Absorber",param_name)
                if val==-1.:
                    return False 

                if param_name == "b":
                    if val<0.1:
                        return False
                    elif model.get_datum(iden,"Absorber","ionName")!="H I":
                        if val>60:   #a metal line with an extremely large b-value: not 
                                     # likely for our line of work
                            return False
                elif param_name == "N":
                    if val<7.0 or val>23.:  #this means dude tried to throw it out.
                                  #if you are running this, then you've already decided that 
                                  #this absorber was needed, so will this is a bad model
                        return False

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
        return models        

    else:
        return [item for item in models if _filter(item)]
            
        
        
        

def random_sampling(model,iden,param,param_range,n,abs_ids,dct,constraints, 
                    method='dude',step=False, verbose=True, iden2=None):
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

    db=ModelDB(models=[]) 
    if type(model) is str:
        if step:
            model=Model(xmlfile=model, abs_ids=abs_ids)   
            #if stepping iteration, do not store all absorbers to conserve memory
        else:
            model=Model(xmlfile=model)
    
    if param=='vel':
        param, param_range=convert_to_vel()

    #a flag to determine whether or not to unlock our param after optimization
    already_locked = bool(model.get_datum(iden,'Absorber',param+"Locked"))

    #set parameter of interest to be locked during optimization
    model.set_val(iden,"Absorber",**{param+"Locked":True})
    assert(model.get_datum(iden,'Absorber',param+"Locked"))

    for i in range(n):
        #set value to random number within range
        perturb_absorbers(dct,model)
        #call dude to optimize
        try:
            run_optimize(model.xmlfile,step=step,verbose=verbose,method=method)
        except KeyboardInterrupt:
            continue
        #add to ModelDB database
        if step:
            newdb=populate_database(abs_ids, keep=False,
                                 path=os.path.split(model.xmlfile)[0],
                                 db=None, constraints=constraints)
            if len(newdb.models)==1:
                is_locked=[]
                for key, val in dct.items():
                    for param in val.keys():
                        if model.get_datum(key,"Absorber",param+"Locked"):
                            is_locked.append("%s: %s is locked"%(key,param))  

                if len(is_locked)>2:
                    msg="random_sampling.py: random_sampling():\n"
                    msg+="optimization procedure didn't run.\n"
                    msg+=str(is_locked)
                    raise Exception(msg)


    
            new_models=filter_bad_models(newdb.models, dct)
            if len(new_models)>0:
                print("appending new models")
                db.append_lst(new_models)
            else:
                print("skipping bad models")
        else:
            model.read(model.xmlfile) #model.xmlfile is rewritten by OptimizeXML
            db.append(model.copy())

    model.set_val(iden,"Absorber",**{param+"Locked":False}) 
    assert(not model.get_datum(iden,'Absorber',param+"Locked"))  


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

    return ModelDB(models=[item for item in lst if item.chi2<lst[0].chi2+delta()])

if __name__=="__main__":
    dct=parse_config()

    glob=dct.pop('config',None)
    if not glob: 
        raise Exception()



    source=glob['source']
    plot=True if glob['plot'].lower()=='true' else False
    append=True if glob['append'].lower()=='true' else False
    step=True if glob['step'].lower()=='true' else False
    step='dude' if glob['method'].lower()=='dude' else None
    n=int(glob['n'])


    
    if append:
        name=input("path/name of db model to append: ")
        if not name=="":
            if name.endswith('.xml'):
                all_db=ModelDB(name=name)
            else:
                all_db=ModelDB.load_models(name)
            
            print("initial filtering")
            all_db=filter_bad_models(all_db, dct)
        else:
            raise Exception("invalid name: %s"%(name)) 
    else:
        all_db=ModelDB(models=[]) 
        
    model=Model(xmlfile=source)
    constraints={}
    count=len(dct.keys())*3*n
    for key, val in dct.items():
        print("\n%d samplings remaining\n"%(count))
        for attr,rng in val.items():
            iden2=dct[iden]["iden2"] if attr=="vel" else None
            count-=1
            try:
                db=random_sampling( model, key, attr, rng, n,
                                    abs_ids=list(dct.keys()),dct=dct,
                                    constraints=constraints,
                                    iden2=iden2, step=step, verbose=False)

                min_chi2=min([item.chi2 for item in db.models])
                #added to remove irrelevant models more than like 10 sigma away                 
                all_db.append_db(get_nsigma(db, n=5.))
            except:
                print("failed either due to timeout, KeyboardInterrupt or other\n")
                pass

    if len(all_db.models)==0:
        raise Exception("no surviving models...")
     
    if plot:
        for key, val in dct.items():
            for attr,rng in val.items():
                plot_chi2(all_db, iden=key, attr=attr,
                          xlabel=r"$%s(%s)$"%(attr,key),
                          constraints={})
    if input("save db?").lower() in ['y', 'yes']:
        name=input("db path/name?") 
        if not name.endswith('.xml'):
            name+='.xml'
        all_db.write(name,True)
        name=name.strip('.xml')+'.obj'
        ModelDB.dump_models(all_db,name)

