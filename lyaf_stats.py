import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import mlab
import spec_parser
from scipy.constants import c
c/=1000.

class DumpData(object):
    def __init__(self, **kwargs):
        self.metal=kwargs.get("metal")
        self.HI=kwargs.get("HI")
        self.metal.absorbers=[item for item in self.metal.absorbers if item.ionName!="HI"]
        self.HI.absorbers=[item for item in self.HI.absorbers if item.ionName=="HI"]
        self.name=os.path.split(self.HI.fname)[-1][0:5]
        self.HI_sp=kwargs.get("HI_sp",None)
        self.metal_sp=kwargs.get("metal_sp",None)
        for attr in 'metal HI name HI_sp metal_sp'.split():
            if not getattr(self,attr):
                print(attr,self.name)
                raise Exception()


    @staticmethod
    def factory(top_dir):
        def verify(dct,name):
            for key in 'HI metal metal_sp HI_sp'.split():
                try:
                    assert(key in dct.keys())
                except:
                    msg="name=%s\nkey=%s"%(str(name),str(key))
                    raise Exception(msg)
            assert(type(dct['HI'])==type(dct['metal'])==spec_parser.LineDump)
            assert(type(dct['HI_sp'])==type(dct['metal_sp'])==spec_parser.TextSpectrum)

        def update_dict(dct,name,key,val):
            if name in dct.keys():
                dct[name][key]=val
            else: 
                dct[name] = {key:val}


        dct={}
        for dirpath, dirname, filename in os.walk(top_dir):

            if not os.path.split(dirpath)[-1] in ["ascii (only HI)", 
                                                  "ascii (only Metal)"]:
                continue
            for f in filename:
                fname=os.path.join(dirpath,f)
                name=f[0:5]
                ab_type="HI" if os.path.split(dirpath)[-1]=="ascii\ \(only\ HI\)" else 'metal'
                if 'spec' in fname:
                    ab_type+='_sp'
                try:
                    data=spec_parser.Spectrum.sniffer(fname)
                except:
                    print('skipping ',fname)
                    continue
                update_dict(dct,name,ab_type,data)
            for key in dct.keys():
                verify(dct[key],fname)
        return [DumpData(**dct[name]) for name in dct.keys()]

    def set_SNR(self,vel=200.,ab_type='HI'):
        def _filter(line,vel=200.,ab_type='HI'):
            sp=getattr(self,ab_type+"_sp")
            if sp.waves[0]>item.obs_wave or line.obs_wave>sp.waves[-1]:
                return 0.
            ref=line.obs_wave

            cut=sp.get_cut(start=ref-vel*ref/c,end=ref+vel*ref/c)
            #if mostly absorption, i.e. within saturated HI line
            if np.mean(cut.abs)<0.3*np.mean(cut.cont):  
                return 0.
            
            return np.median(cut.cont)/np.median(cut.error)
 

        line_dump=getattr(self,ab_type)
        lines=[] 
        for ab in line_dump.absorbers:
            for item in ab.get_lines():
                lines.append(_filter(item))
            ab.SNR=max(lines)

        


def plot_hist(dat,attr, bin_rng=None, *args, **kwargs):
    num_bins=kwargs.get('num_bins',15)
    dat=[getattr(item,attr) for item in dat]
    n, bins, patches = plt.hist(dat, 
                                 num_bins, normed=1, facecolor='green', 
                                 alpha=0.75)
    y = mlab.normpdf( bins, np.mean(dat), np.std(dat))
    l = plt.plot(bins, y, 'r--', linewidth=1)
    plt.xlabel(attr)
    ab_type=kwargs.get("absorber_type","")
    plt.title(r'%s $%s=%3.2lf \pm %3.2lf  (n=%d)$'%(ab_type,attr,np.mean(dat), np.std(dat),len(dat)))
    if bin_rng:
        plt.ylabel("%s from z=%2.1lf to %2.1lf"%(attr, bin_rng[0],bin_rng[-1]))
    plt.show()

    return 

def bin_absorbers(dat, binsize, bin_start, bin_end, bin_attr='z'):
    output=[]
    b=bin_start
    for b in np.arange(bin_start, bin_end, binsize).tolist():
        output.append(
            [item for item in dat if b<=getattr(item,bin_attr)<b+binsize]
        )

    return output
        
def histogram_by_z(spec, binsize=0.3, bin_start=0.,bin_end=5.,*args,**kwargs):
    n_min=kwargs.get('n_min',30)
    all_abs=[]
    absorber_type=kwargs.get("absorber_type","HI")
    for sp in spec:
        all_abs+=getattr(sp,absorber_type).absorbers

    if kwargs.get("exclude",False):
        all_abs = [it for it in all_abs if it.ionName!=kwargs.get("exclude")]
    elif kwargs.get("ion",False):
        all_abs = [it for it in all_abs if it.ionName==kwargs.get("ion")] 

    data=bin_absorbers(all_abs,binsize,bin_start,bin_end)
    b=bin_start
    for sublst in data:
        if len(sublst)>30:
            for attr in ['N', 'b']:
                plot_hist(sublst,attr, bin_rng=[b,b+binsize],*args, **kwargs)
        b+=binsize

    z=[np.mean([item.z for item in it]) for it in data if len(it)>30]
    mu=[np.mean([item.N for item in it]) for it in data if len(it)>30]
    stdev=[np.std([item.N for item in it]) for it in data if len(it)>30]
    num_lines=[len(it) for it in data if len(it)>30]

    z,mu,stdev,num_lines=tuple(map(np.array, (z,mu,stdev,num_lines)))

    stderr=stdev/np.sqrt(num_lines)

    plt.clf()
    #plt.errorbar(z,mu,yerr=[stdev[i]/num_lines[i] for i in range(len(num_lines))], fmt='o')
    plt.errorbar(z,mu,yerr=stderr, fmt='o')    
    plt.xlabel('mean z')
    plt.ylabel('mean N (%s)'%(absorber_type))
    plt.show()

     
#TODO: scale bins for comoving volume  
#TODO: scatter binned N vs z and binned b vs z for met and HI 4 plots
#TODO: same as previous, but split by ionization:  mgII, CIV
  

if __name__=="__main__":
    spec=DumpData.factory('/home/scott/research/lyaf_stats')
    for item in 'HI metal HI_sp metal_sp name'.split():
        print("type=%s type=%s"%(item,type(getattr(spec[0],item))))

    """for item in spec:
        item.set_SNR(ab_type='HI')
        item.set_SNR(ab_type='metal')

    histogram_by_z(spec, absorber_type="HI")
    histogram_by_z(spec, absorber_type="metal")"""

    




 
        
