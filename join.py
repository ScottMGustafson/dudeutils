from model import ModelDB, Model
import data_types
import os

def join(db1_name, db2_name, string='hi_', newdb="newdb.xml", outname='scratch.xml'):
    #Write something up for applying id changes to dict.

    #add some string identifier to every id, to avoid conflicts with other 
    #db.  the id doesn't actually matter, btw.
    f=open(db2_name,'r')
    _f = open(outname,'w')
    lines = f.readlines()
    for line in lines:
        keys=list(Model.model_classes.keys())+["id"]
        for _str in keys:
            _str = _str+"=\""
            if _str in line:
                line = line.replace(_str, _str+str(string))
        _f.write(line)
    f.close()
    _f.close()

    #now that there should be no conflicting ids, try combining dbs
    db1 = ModelDB.read(db1_name)
    pool1 = dict(data_types.ObjList._pool)
    data_types.ObjList._pool = {}  #this should already be done in ModelDB.read()


    db2 = ModelDB.read(outname)
    for key, val in data_types.ObjList._pool.items():
        if key in pool1.keys(): #correct key errors
            raise Exception("conflicting keys: "+str(key))

    #merge the data caches
    data_types.ObjList._pool=dict(list(data_types.ObjList._pool.items())+list(pool1.items()))
    db2.models += db1.models

    db2.write(newdb)
           
2014-11-24db.xml  2014-11-28db.xml  2014-11-28hidb.xml  
2014-11-28lodb.xml  2014-12-08db.xml  2014-12-14_hidb.xml  
2014-12-14_lodb.xml
