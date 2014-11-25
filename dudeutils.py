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

