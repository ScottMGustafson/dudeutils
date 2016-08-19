"""
Some general constraint classes for constraining models.
"""

class Constraint(object):
    def __init__(self,**kwargs):
        """
        input:
        -----
        **kwargs:  dict of constraints.  
        format like :
            #for absorbrs
            {'id1': {param1: [min, max], param2:[min,max]}, 
             'id2': {param3: [min, max]}
            #for continuum points
            'continuum':{
                        'pt1': {x:[min,max], y:[min,max]}
                        }
           
            }
        """
        self.abs=[]
        self.cont=[]
        for key, val in kwargs.items():
            if type(val)==dict:  #only an absorber would have dict input
                if key=='continuum':
                    self.continuum.append(ContinuumConstraint(key,**val))
                else:
                    self.abs.append(AbsConstraint(key,**val))
            elif key=='chi2':
                if not type(val) in [list, tuple]:
                    val=[0,float(val)]
                self.chi2 = ParameterConstraint("chi2",val)
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
                if not float(getattr(model,key)) in getattr(self,key):
                    return False
            elif key=='continuum':
                for cont_ in self.cont:
                    for cont in model.get_lst('ContinuumPointList'):
                        if cont.id==cont_.id:
                            if not cont in cont_:
                                return False
            else:
                for ab_ in self.abs:
                    for ab in model.get_lst("AbsorberList"):
                        if ab.id==ab_.id:
                            if not ab in ab_:
                                return False
        return True

class AbstractConstraint(object):
    def __init__(self,id,**kwargs):
        self.id = id
        for key, val in kwargs.items():
            setattr(self, str(key), ParameterConstraint(str(key),val))

class ContinuumConstraint(AbstractConstraint):
    def __init__(self,id,**kwargs):
        super().__init__(id,**kwargs)
    def __contains__(self,absorber):
        if absorber.id!=self.id:
            return False
        keys = [key for key in self.__dict__.keys() if key in ['x','y']]
        for key in keys:
            val = getattr(self,key)
            if not getattr(absorber,key) in getattr(self,key):
                return False
        return True

class AbsConstraint(AbstractConstraint):
    def __init__(self,id,**kwargs):
        super().__init__(id,**kwargs)
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
        if type(rng) in [tuple, list]:
            self.rng = list(map(float,rng))
        else:
            self.rng=rng
            
    def __contains__(self,val):
        try:
            if type(val) in [float, int]:
                return self.rng[0]<=float(val)<=self.rng[1]
            elif type(val) is bool:
                print(val, self.rng)
                return val==self.rng
        except IndexError:
            raise IndexError("input range must be list or tuple of length 2.  Instead got:\n  %s"%(str(self.rng)))
   
