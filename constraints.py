class Constraint(object):
    def __init__(self,**kwargs):
        self.abs=[]
        for key, val in kwargs.items():
            if type(val)==dict:
                self.abs.append(AbsConstraint(key,**val))
            elif key=='chi2':
                Constraint.to_range(self,key,val) 
            else:  #pixel, param   
                try: 
                    setattr(self,key,float(val))
                except:
                    if key in ["pixels","params"]:
                        raise
                    else:
                        setattr(self,key,val)
                
    def compare(self,model):
        pix = self.__dict__.get("pixels",model.pixels)
        params=self.__dict__.get("params",model.params)
        xmlfile=self.__dict__.get("xmlfile")
        if model.pixels!=pix or params!=model.params or xmlfile!=model.xmlfile:
            return False
        if "chi2" in self.__dict__.keys():
            if not self.chi2[0]<=model.chi2<=self.chi2[1]:
                return False

        for ab in model.absorbers:
            for item in self.abs:
                if item.id==ab.id:
                    if not item.compare(ab):
                        return False
        return True

    @staticmethod
    def to_range(cls,key,val):
        """convert provided vales to range"""
        try:
            val=tuple(map(float,val))
        except:
            val=tuple(map(float,[0,val]))
        setattr(cls,key,val)

class AbsConstraint(object):
    def __init__(self,id,**kwargs):
        self.id = id
        for key, val in kwargs.items():
            Constraint.to_range(self,key,val)

    def compare(self,absorber):
        assert(absorber.id==self.id)
        for key in ["N","b","z"]:
            try:
                val = getattr(self,key)
                if not val[0]<=getattr(absorber,key)<=val[-1]:
                    return False
            except:
                pass
        return True

   
