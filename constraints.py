class Contraint(object):
    def __init__(self,**kwargs):
        for key, val in kwargs.items():
            try:
                setattr(self,key,tuple(map(float,val)))
            except:
                setattr(self,key,tuple(map(float,[0,val])))
        
    def compare(self):
        for key, val in self.__dict__.items():
            pass

class AbsConstraint(object):
    def __init__(self,**kwargs):
        self.id = kwargs.pop("id")
        for key, val in kwargs.items():
            try:
                setattr(self,key,tuple(map(float,val)))
            except:
                setattr(self,key,tuple(map(float,[0,val])))
            

