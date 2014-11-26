class Constraint(object):
    def __init__(self,**kwargs):
        self.abs=[]
        for key, val in kwargs.items():
            if type(val)==dict:  #only an absorber would have dict input
                self.abs.append(AbsConstraint(key,**val))
            elif key=='chi2':
                self.chi2 = ParamterConstraint("chi2",val)
            else:  #pixel, param   
                try: 
                    setattr(self,key,float(val))
                except:
                    if key in ["pixels","params"]:
                        raise
                    else:
                        setattr(self,key,val) 
                
    def __contains__(self,model):
        for key, val in self.__dict__.items():
            if key in ["pixels","params"]:
                if getattr(self,key)!=getattr(model,key):
                    return False
            elif key=='chi2':
                if not getattr(model,key) in getattr(self,key):
                    return False
            else:
                for ab_ in self.abs:
                    for ab in model.get_lst("AbsorberList"):
                        if ab.id==ab_.id:
                            if not ab in ab_:
                                return False
        return True


class AbsConstraint(object):
    def __init__(self,id,**kwargs):
        self.id = id
        for key, val in kwargs.items():
            setattr(self, str(key), ParameterConstraint(str(key),val))

    def __contains__(self,absorber):
        if absorber.id!=self.id:
            return False
        keys = [key for key in self.__dict__.keys() if key in ['N', 'b', 'z']]
        for key in keys:
            val = getattr(self,key)
            if not getattr(absorber,key) in getattr(self,key):
                return False
        return True


class ParameterConstraint(object):
    def __init__(self,id,rng):
        self.id=id
        self.rng = list(map(float,rng))
            
    def __contains__(self,val):
        try:
            return self.rng[0]<=val<=self.rng[1]
        except IndexError:
            raise IndexError("input range must be list or tuple of length 2.  Instead got:\n  %s"%(str(self.rng)))
   
