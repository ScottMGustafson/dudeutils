from dudeutils import *
from model import *
from plot_distribution import plot_chi2
from configparser import ConfigParser

def parse_config(config_file=os.path.join('data','random_sampling_config.cfg')):
    config = ConfigParser()
    config.optionxform=str
    config.read(config_file)
    dct={}
    for item in list(config.sections()):
        dct[item]=dict(config.items(item))
    return dct

def perturb_absorbers(dct, model):
    for key, val in dct.items():
        for param, rng in val.items():
            model.monte_carlo_set(key,"Absorber",rng,param,False)
    model.write(model.xmlfile)
    return model

def random_sampling(model,iden,param,param_range,n,abs_ids,dct,constraints, step=False, iden2=None):
    """
    runs simmulated annealing algorithm on dude models

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

    """

    def convert_to_vel():
        """
        if if user instead specifies velocity from some other absorber, 
        convert to redshift, assuming the other is locked
        """
        zref=model.get_datum(iden2,"Absorber","z")
        try:
            assert(model.get_datum(iden2,"Absorber","zLocked"))
        except:
            model.set_val(iden2,"Absorber",zLocked='true')
            
        vel=param_range
        param_range=[(vel[0]/c)*(1.+zref)+zref, (vel[1]/c)*(1.+zref)+zref]
        param='z'

        return param, param_range

    db=ModelDB(models=[]) 
    if type(model) is str:
        if step:
            model=Model(xmlfile=model, abs_ids=abs_ids)   #if stepping iteration, do not store all absorbers
        else:
            model=Model(xmlfile=model)
    
    if param=='vel':
        param, param_range=convert_to_vel()

    #a flag to determine whether or not to unlock our param after optimization
    if model.get_datum(iden,'Absorber',param+"Locked")=='true':
        already_locked=True  
    else: 
        already_locked=False

    model.set_val(iden,"Absorber",**{param+"Locked":"true"})

    for i in range(n):
        #set value to random number within range
        perturb_absorbers(dct,model)
        #call dude to optimize
        run_optimize(model.xmlfile,step)
        #add to ModelDB database
        if step:
            db=populate_database(abs_ids, keep=False,
                                 path=os.path.split(model.xmlfile)[0],
                                 db=db, constraints=constraints)
        else:
            model.read(model.xmlfile) #model.xmlfile is rewritten by OptimizeXML
            db.append(model.copy())

    if not already_locked:
        model.set_val(iden,"Absorber",**{param+"Locked":"false"})
        model.write(model.xmlfile)

    return db


def get_nsigma(db,n=1):
    def delta(): return float(n**2)
    lst=sorted(db.models, key=lambda x: x.chi2)

    return ModelDB(models=[item for item in lst if item.chi2<lst[0].chi2+delta()])

if __name__=="__main__":
    dct=parse_config()

    glob=dct.pop('config',None)
    if not glob: 
        raise Exception()

    for key, val in dct.items():
        for k, rng in val.items():
            dct[key][k]=list(map(float,rng.strip().split(', ')))

    source=glob['source']
    plot=True if glob['plot'].lower()=='true' else False
    append=True if glob['append'].lower()=='true' else False
    step=True if glob['step'].lower()=='true' else False
    n=int(glob['n'])

    model=Model(xmlfile=source)
    constraints={}

    all_db=ModelDB(models=[])

    for key, val in dct.items():
        for attr,rng in val.items():
            iden2=dct[iden]["iden2"] if attr=="vel" else None
            db=random_sampling( model, key, attr, rng, n,
                                    abs_ids=list(dct.keys()),dct=dct,
                                    constraints=constraints,
                                    iden2=iden2, step=step)

            if append:
                all_db.append_lst(db.models)
    if plot:
        for key, val in dct.items():
            for attr,rng in val.items():
                plot_chi2(all_db, iden=key, attr=attr,
                          xlabel=r"$%s(%s)$"%(attr,key),
                          constraints={})

