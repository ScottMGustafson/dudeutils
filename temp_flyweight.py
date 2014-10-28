import weakref
import uuid


class ObjList(object):
    """this and associated subclasses are simply extended lists"""

    _pool=weakref.WeakValueDictionary()

    def __new__(cls, objlst,taken_names=[],id=None):
        lists = [item.objlst for item in ObjList._pool.values()]
        if not objlst in lists: #if not already in pool
            obj = object.__new__(cls)
            obj.id = ObjList.generate_id(taken_names) if id==None else id
            ObjList._pool[obj.id] = obj
        else:
            i=0
            while objlst!=lists[i].objlst:
                i+=1
            obj=_pool[lists[i].id]  #return the old object
        return obj

    def __init__(self,objlst,taken_names=[],id=None):
        self.cls = objlst[0].__class__
        self.objlist = objlst
        
    def __eq__(self,other):
        for i in range(len(self.objlst)):
            try:
                assert(self.objlst[i]==other.objlst[i])
            except IndexError:
                return False
            except AssertionError:
                if not self.objlst[i] in other.objlst:
                    return False
            except:
                raise
        return len(other.objlst)==len(self.objlst)  #last check     
        
    def __neq__(self,other):
        return not self.__eq__(other)

    def __iter__(self):
        for i in range(len(self.objlist)):
            yield self.objlist[i]               

    def __getitem__(self,i):
        return self.objlist[i]

    @staticmethod
    def generate_id(taken_names):
        iden=uuid.uuid4()
        while iden in taken_names:
            iden=uuid.uuid4()  

        taken_names.append(iden)
        return str(iden)

