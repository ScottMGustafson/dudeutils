from model import Model, ModelDB

"""a collection of helper functions"""

def newdb(xmlfile,chi2,pixels,dbfile=None,params=None,**kwargs):
    """get a model, append to new database"""
    model = Model(xmlfile=xmlfile,chi2=chi2,pixels=pixels)
    return ModelDB(dbfile,[model],**kwargs)

def load_from_db(dbxml):
    return ModelDB.read(dbxml)

def get_db(dbfile):
    return ModelDB.read(dbfile)

def random(xmlfile,itemid,tag,param,val_range,modelid=None):
    """get a new random setting a parameter to a random value in val_range (tuple or list)"""
    mod=Model(xmlfile=xmlfile,chi2=0,pixels=0,params=0,id=modelid)
    old = mod.get_datum(itemid,tag,"N")
    mod.monte_carlo_set(itemid,tag,param,val_range)
    new=mod.get_datum(itemid,tag,"N")
    assert(old!=new)
    print("model written to %s"%(xmlfile))

def all_conts(db=None):
    """get the set of all continua"""
    if type(db) is str:
        db = load_from_db(db)
    if db is None: 
        db = load_from_db('database.xml')
    conts=[]
    for item in db:
        conts.append((item.ContinuumPointList, item.xmlfile))
    return list(set(conts))

def cont_check_pipeline(reduced_chi2_limit=1.8,verbose=True, db=None):
    if db is None: 
        db = load_from_db('database.xml')
    conts = all_conts(db)

    ids = [item[0] for item in conts]
    for item in ids:
        assert(ids.count(item)==1)
    out=[]  #the 'good' continua
    for item in conts:
        models = db.get_lst_from_id(id=item[0],attr='ContinuumPointList')
        _min = min([float(mod.reduced_chi2) for mod in models])
        _min_chi2=min([float(mod.chi2) for mod in models])
        if _min <= reduced_chi2_limit:
            out.append((item[0],item[1],_min,_min_chi2))
    if verbose:
        print(len(out))
        for item in out:
            msg="%15s, id=%15s: \n    reduced chi2 = %6.4f, chi2=%6.4f"%(item[0],item[1],item[2],item[3])
            print(msg)
        print('\n')
    return out

def getDH(vr=(-0.5,0.5)):
    import matplotlib.pyplot as plt

    db = load_from_db('database.xml') 
    allmods = []
    _conts = cont_check_pipeline(db=db)
    conts = [item[0] for item in _conts]
    for item in db:
        if vr[0]<=item.get_shift('D','H')<=vr[1] and item.ContinuumPointList in conts:
            allmods.append(item)

    for cont in _conts:
        id = cont[0]
        name=cont[1]
        mods = []
        for mod in db.models:
            if vr[0]<=mod.get_shift('D','H')<=vr[1] and mod.ContinuumPointList==id:
                mods.append(mod)

        #dh = [{'dh':item.dh, 'chi2':item.chi2} for item in mods]
        plt.plot([item.dh for item in mods], [item.reduced_chi2 for item in mods],'ko')
        plt.title(name)
        plt.show()
        

    plt.plot([item.dh for item in allmods], [item.chi2 for item in allmods],'ko')
    plt.title("all continua")
    plt.show()

def check_for_cont_duplicates():
    import datetime
    import os
    def print_pretty(lst):
        for item in lst:
            print("  %15s  %15s"%(item[0],item[1]))
    def get_new_name(name):
        name, ext = os.path.splitext(name)
        num = 1
        new_name = name+str(num)+ext
        while os.path.isfile(new_name):
            num+=1
            new_name=name+str(num)+ext
        return new_name

    db = load_from_db('database.xml')
    conts = all_conts(db)
    names = [item[1] for item in conts]
    duplicates = list(set([(item[0],item[1]) for item in conts if names.count(item[1]) > 1 ]))
    print_pretty(duplicates)
    
    #rewrite duplicates to different xmlfiles
    for item in duplicates:
        models = db.get_lst_from_id(id=item[0],attr='ContinuumPointList')
        models = sorted( models, key=lambda x: x.chi2 )
        best_fit = models[0]
        new_name = get_new_name(best_fit.xmlfile)
        for i in range(0,len(db.models)):  #update db
            if db[i] in models:
                db[i].xmlfile=new_name  #update for model
            
        best_fit.write(filename=new_name)
                
    #write changes.  to avoid overwriting potentially desirable data, append date, time to fname
    db.write(str(datetime.datetime.now().isoformat()+"_db.xml"))

if __name__=='__main__':
    getDH()
        
    
    
    
