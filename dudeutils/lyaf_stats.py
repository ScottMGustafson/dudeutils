import pickle, sys, os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import mlab
import dudeutils.spec_parser as spec_parser
from scipy.constants import c
from dudeutils.lyaf_stats_dict import names
c/=1000.





class DumpData(object):
    def __init__(self, **kwargs):
        self.metal=kwargs.get("metal")
        self.HI=kwargs.get("HI")
        assert(type(self.metal)==type(self.HI)==spec_parser.LineDump)

        self.m1=[item for item in self.metal.absorbers 
                 if item.ionName.lower()=='m1'] #unidentified metals

        self.metal.absorbers=[item for item in self.metal.absorbers 
                              if item.ionName.lower() not in ["HI", "m1"] ]

        self.HI.absorbers=[item for item in self.HI.absorbers 
                            if item.ionName=="HI"]

        self.name=os.path.split(self.HI.fname)[-1][0:5]
        self.HI_sp=kwargs.get("HI_sp",None)
        self.metal_sp=kwargs.get("metal_sp",None)


    @staticmethod
    def factory(top_dir):
        def verify(dct,name):
            for key in 'HI metal metal_sp HI_sp'.split():
                try:
                    assert(key in dct.keys())
                except:
                    print(dct.keys())
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
                ab_type="HI" if os.path.split(dirpath)[-1]=="ascii (only HI)" else 'metal'
                if 'spec' in fname:
                    ab_type+='_sp'
                try:
                    data=spec_parser.Spectrum.sniffer(fname)
                except:
                    print('skipping ',fname)
                    continue
                update_dict(dct,name,ab_type,data)
        for key in dct.keys():
            verify(dct[key],key)
        return [DumpData(**dct[name]) for name in dct.keys()]

    def set_SNR(self,vel=200.,ab_type='HI',verbose=True):
        def _filter(line,vel=200.,ab_type='HI'):
            sp=getattr(self,ab_type+"_sp")
            if sp.waves[0]>item.obs_wave or line.obs_wave>sp.waves[-1]:
                return 0.
            ref=line.obs_wave

            cut=sp.get_cut(start=ref-vel*ref/c,end=ref+vel*ref/c)
            #if mostly absorption, i.e. within saturated HI line
            if not ab_type is 'HI':
                if np.mean(cut.abs)<0.3*np.mean(cut.cont):  
                    return 0.
            
            return np.median(cut.cont)/np.median(cut.error)
 

        line_dump=getattr(self,ab_type)
        lines=[] 
        if verbose: print("getting SNR for %s"%(self.name))
        for ab in line_dump.absorbers:
            for item in ab.get_lines():
                lines.append(_filter(item))

            ab.SNR=max(lines)
            if verbose:  
                sys.stdout.write('.')
                sys.stdout.flush()
        print('\n')
         

    @staticmethod
    def bin_absorbers(dat, bin_attr='z', binsize=None, bin_start=None, bin_end=None, bins=None):
        output=[]
        if not bins:   #if no custom bins specified, then take equallly spaced bins
            bins=np.arange(bin_start, bin_end, binsize)
            bins=list(zip(bins,bins+binsize))

        for b in bins:
            output.append(
                [item for item in dat if b[0]<=getattr(item,bin_attr)<b[1]]
            )
        return output


    @staticmethod
    def split_abs(obj1, absorbers):
        if type(absorbers) is str:
            return obj1.split(absorbers)
        elif type(absorbers) is list:
            return [obj1.split(item) for item in absorbers]
        else: 
            raise Exception("invalid input type: must be list or str")

def plot_hist(dat,attr, bin_rng=None, *args, **kwargs):
    if type(dat) is np.ndarray:
        dat=dat.tolist()
    assert(type(dat) is list)
    num_bins=kwargs.get('num_bins',15)
    dat=[getattr(item,attr) for item in dat]
    n, bins, patches = plt.hist(dat, 
                                 num_bins, normed=1, facecolor='green', 
                                 alpha=0.75)
    y = mlab.normpdf( bins, np.mean(dat), np.std(dat))
    l = plt.plot(bins, y, 'r--', linewidth=1)
    plt.xlabel(kwargs.get('xlabel', attr))

    ab_type=kwargs.get("absorber_type","")
    plt.title(r'%s $%s=%3.2lf \pm %3.2lf  (n=%d)$'%(ab_type,attr,np.mean(dat), np.std(dat),len(dat)))

    if kwargs.get("ylabel",None):
        plt.ylabel(kwargs.get("ylabel"))
    elif bin_rng:
        plt.ylabel("%s from z=%2.1lf to %2.1lf"%(attr, bin_rng[0],bin_rng[-1]))
    plt.show()

    return 



def flux_pdf(spec,rest_range=(1090,1170),**kwargs):
    if type(spec) is list:
        assert(type(spec[0]))
        for item in spec:
            print((1.+item.zem)*rest_range[0], (1.+item.zem)*rest_range[1])
            flux=item.HI_sp.get_cut(
                                    start=(1.+item.zem)*rest_range[0], 
                                    end=(1.+item.zem)*rest_range[1]
                                ).flux
            try:
                sp=np.concatenate(sp,flux)
            except:
                sp=flux

        print(sp.shape)
    else:
        assert(type(spec) is DumpData)
        zem, name = spec.zem, spec.name
        sp=spec.HI_sp.get_cut(start=(1.+zem)*rest_range[0], end=(1.+zem)*rest_range[1]).flux
        
    num_bins=kwargs.get('num_bins',40)
    n, bins, patches = plt.hist(sp, 
                                 num_bins, normed=1, facecolor='green', 
                                 alpha=0.75)
    #y = mlab.normpdf( bins, np.mean(sp.flux), np.std(sp.flux))
    #l = plt.plot(bins, y, 'r--', linewidth=1)
    plt.xlabel("flux")
    plt.title(kwargs.get("title",""))
    plt.show()
    
        
def histogram_by_z(spec, binsize=0.4, bin_start=0.,bin_end=5.,*args,**kwargs):
    n_min=kwargs.get('n_min',30)
    all_abs=[]
    absorber_type=kwargs.get("absorber_type","HI")
    for sp in spec:
        all_abs+=getattr(sp,absorber_type).absorbers

    if kwargs.get("exclude",False):
        all_abs = [it for it in all_abs if it.ionName!=kwargs.get("exclude")]
    elif kwargs.get("ion",False):
        all_abs = [it for it in all_abs if it.ionName==kwargs.get("ion")] 

    data=DumpData.bin_absorbers(all_abs,binsize=binsize,bin_start=bin_start,bin_end=bin_end)
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
    """spec=DumpData.factory('/home/scott/research/lyaf_stats')

    for key, val in names.items():
        for item in spec:
            if key in item.name:
                item.name=key
                item.zem=val
                continue
        
    for item in spec:
        item.set_SNR(ab_type='HI')
        item.set_SNR(ab_type='metal')"""

    #pickle.dump(spec,open("spec_dump.obj",'wb'))
    spec=pickle.load(open("spec_dump.obj",'rb'))
    for key, val in names.items():
        for item in spec:
            if key.strip() == item.name.strip():
                item.zem=val
                continue
    #flux_pdf([sp for sp in spec if hasattr(sp,'zem')],title="all spectra")

    met, HI,m1=[],[],[]
    for item in spec:
        if hasattr(item,'zem'):
            print(item,item.zem,item.name)
            #flux_pdf(item,title=item.name+" flux")
            met+=[it for it in item.metal.absorbers]
            HI+=[it for it in item.HI.absorbers]
            m1+=[it for it in item.m1]

    #bin by z
    zbins=np.arange(0,5,0.4)
    zbins=list(zip(zbins,zbins+0.4))
    met=DumpData.bin_absorbers(met,bins=zbins)
    HI=DumpData.bin_absorbers(HI,bins=zbins)

    #bin by SNR
    #bins=[(0,2.5),(2.5,5),(5,10),(10,20),(20,40),(40,80),(80,160)]
    bins=list(zip(np.arange(0,100,10),np.arange(0,100,10)+10.))
    met=[DumpData.bin_absorbers(it,bin_attr='SNR',bins=bins) for it in met]
    HI=[DumpData.bin_absorbers(it,bin_attr='SNR',bins=bins) for it in HI]
    #m1=DumpData.bin_absorbers(m1,bin_attr='SNR',bins=bins)

    i,j,k=0,0,0


    mean_metals=[]
    mean_HI=[]
    """
    for lst in [met, HI]:
        for i in range(len(zbins)):
            if i>len(lst): continue
            zbin=lst[i]
            for j in range(len(bins)):
                if j>len(zbin): continue
                SNRbin=zbin[j]
                if len(SNRbin)<30:                 
                    continue
                plot_hist(SNRbin,"N", 
                            xlabel="%s: z=%3.1lf--%3.1lf,   SNR=%4.1lf--%4.1lf"%(
                              "N",zbins[i][0],zbins[i][1],bins[j][0],bins[j][1]
                            ),
                            ylabel="metal" if lst is met else "HI",
                          )
                if lst is met:
                    mean_metals.append((zbins[i][0], bins[j][0], SNRbin))
                elif lst is HI:
                    mean_HI.append((zbins[i][0], bins[j][0], SNRbin))
    """
    #i= zbin, j=snrbin
    

    all_met, all_HI=[],[]
    for sp in spec:
        all_met+=sp.metal.absorbers
        all_HI+=sp.HI.absorbers

    x=[(item[0]+item[1])/2. for item in bins]
    met=DumpData.bin_absorbers(all_met,bin_attr='SNR',bins=bins)
    HI=DumpData.bin_absorbers(all_HI,bin_attr='SNR',bins=bins)
    

    
    plt.errorbar(x,
                 [np.mean([item.N for item in y]) for y in met],
                 yerr=[np.std([item.N for item in y]) for y in met], 
                 fmt='o') 
    plt.title("N vs SNR: metals")
    plt.show() 
    plt.errorbar(x,
                 [np.mean([item.N for item in y]) for y in HI],
                 yerr=[np.std([item.N for item in y]) for y in HI], 
                 fmt='o')  
    plt.title("N vs SNR: HI")
    plt.show() 

    #histogram_by_z(spec, absorber_type="HI")
    #histogram_by_z(spec, absorber_type="metal")


    




 
        
