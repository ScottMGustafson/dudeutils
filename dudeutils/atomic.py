from dudeutils.get_data import get_data

class AtomicData(object):
    def __init__(self):
        self.atomic_data = AtomicData.get_lines()

    @staticmethod
    def get_lines(fname=get_data('atom.dat')):
        all_lines={}
        f=open(fname,'r')
        for line in f:
            ion, line= line[0:5].strip(), line[5:]
            line=line.split()
            try:
                all_lines[ion].append(
                    SpectralLine(**{
                        'ionName':ion,
                        'wave': float(line[0]),
                        'f':    float(line[1]),
                        'gamma':float(line[2])
                    }))
            except KeyError:
                all_lines[ion]=[
                    SpectralLine(**{
                        'ionName':ion,
                        'wave': float(line[0]),
                        'f':    float(line[1]),
                        'gamma':float(line[2])
                    })]
        for k in all_lines.keys():
            all_lines[k] = sorted(all_lines[k], 
                                  key=lambda item:item.wave, reverse=True)
        return all_lines



class SpectralLine(object):
    mass_dict={
            'H':1.008,
            'D':2.014,
            'He':4.003,
            'Li':6.941,
            'Be':9.012,
            'B':10.081,
            'C':12.011,
            'N':14.007,
            'O':15.999,
            'F':18.998,
            'Ne':20.180,
            'Na':22.990,
            'Mg':24.305,
            'Al':26.982,
            'Si':28.086,
            'P':30.974,
            'S':32.066,
            'Cl':35.453,
            'Ar':39.948,
            'K':39.098,
            'Ca':40.078,
            'Sc':44.956,
            'Ti':47.867,
            'Cr':51.996,
            'Mn':54.938,
            'Fe':55.845,
            'Co':58.933,
            'Ni':58.693,
            'Cu':63.546,
            'Zn':65.38,
            'Ge':72.631
    }

    def __str__(self):
        return "%s: rest wave=%6.2lf A f=%lf"%(self.ionName, self.wave, self.f)

    def __init__(self,**kwargs):
        for key, val in kwargs.items():
            setattr(self,key,val)
        self.rest_wave=self.wave #set an alias for this, since I always forget the name

    def get_obs(self,z=None):
        try:
            if not z:
                z=self.z
            return (1.+z)*self.wave
        except:
            raise Exception(str(self.__dict__))

    @staticmethod
    def get_atom_name(ionName):
        ionName=ionName[:2]
        if ionName[-1] in 'IV ':
            return ionName[0]
        else:
            return ionName


#call on import
atomic_data=AtomicData().atomic_data  
