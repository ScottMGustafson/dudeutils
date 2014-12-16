from model import ModelDB

def join(db1_name, db2_name, db2string='hi', newdb="newdb.xml"):
    #Write something up for applying id changes to dict.

    #add some string identifier to every id, to avoid conflicts with other 
    #db.  the id doesn't actually matter, btw.
    f=open(db2string,'r')
    _f = open("scratch.xml",'w')
    for line in f.lines():
        line.replace("id=\"", "id=\""+str(db2string))
        _f.write(line)
    f.close()
    _f.close()


    db1 = ModelDB.read(db1_name)
    pool1 = data_types.ObjList._pool

    db2 = ModelDB.read("scratch.xml")
    for key, val in data_types.ObjList._pool.items():
        if key in pool1.keys(): #correct key errors
            raise Exception("conflicting keys")
    
    data_types.ObjList._pool+=pool1
    db2.models += db1.models

    db2.write(newdb)
            
    
            
    #there is a more elegant way to do this, but the lazyway is to just 
    #append all ids in an xmldb with some string, and then re-write the database
