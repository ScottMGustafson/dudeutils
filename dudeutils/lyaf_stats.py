import pickle, sys, os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import mlab
import dudeutils.spec_parser as spec_parser
from scipy.constants import c
from dudeutils.lyaf_stats_dict import names
from dudeutils.model import Model
import numpy.ma as ma
c/=1000.


plt.rc('text', usetex=True)

class AbSystem(object):
    def __init__(self, ab_lst, name=None):
        self.HI=ab_lst.pop(0)
        self.metals=ab_lst
        self.name=name
        self.z=self.HI.z
        self.data=[]
        self.split_by_ionization()
        
    def split_by_ionization(self):
        for i, ionization in enumerate("I II III IV V VI".split()):
            tmp=[]
            for ab in self.metals:
                if ab.ionName.endswith(ionization):
                    tmp.append(ab)
            self.data.append(tmp)

    @staticmethod
    def factory(ab_list, vel_tol=20.):
        """
        input:
        ------
        ab_list: list of Aborber Instances, without m1

        output:
        -------
        AbSystem instance
        """
        HIlst= [item for item in ab_list if item.ionName=="HI"]
        metlist=[item for item in ab_list if item.ionName not in "HI m1".split()]
        out=[]
        for HI in HIlst:
            newlst=[HI]
            for met in ab_list:
                if np.fabs(met.z-HI.z)*c/(1.+HI.z) < vel_tol:
                    newlst.append(met)
            out.append(AbSystem(newlst))
        return out
            
                     

class DumpData(object):
    def __init__(self, **kwargs):
        self.model=kwargs.get("model",None)
        self.metal=kwargs.get("metal",
                              spec_parser.LineDump(self.model))
        self.HI=kwargs.get("HI",
                           spec_parser.LineDump(self.model))
        self.HIonly=kwargs.get("HI",
                           spec_parser.LineDump(self.model))

        assert(type(self.metal)==type(self.HI)==spec_parser.LineDump)

        self.m1=[item for item in self.metal.absorbers 
                 if item.ionName.lower()=='m1'] #unidentified metals

        self.metal.absorbers=[item for item in self.metal.absorbers 
                              if item.ionName.lower() not in ["HI", "m1"] ]

        self.HI.absorbers=[item for item in self.HI.absorbers 
                            if item.ionName=="HI"]

        self.HIonly.absorbers=[item for item in self.HIonly.absorbers 
                            if item.ionName=="HI"]

        self.name=os.path.split(self.HI.fname)[-1][0:5]
        self.zem=names[self.name]
        self.HI_sp=kwargs.get("HI_sp",None)
        self.HIonly_sp=kwargs.get("HIonly_sp",None)
        self.metal_sp=kwargs.get("metal_sp",None)

    
    @staticmethod
    def factory(top_dir):
        dct={}
        models={}
        HI={}
        HIonly={}
        metals={}
        HIabs={}
        HIonly_abs={}
        metalabs={}
        for dirpath, dirname, filename in os.walk(top_dir):
            if str(os.path.split(dirpath)[-1]).strip() in ["ascii (only HI)", 
                                                  "ascii (only Metal)",
                                                  "ascii (HI)",
                                                  "xmls"]:


                for f in filename:
                    fname=os.path.join(dirpath,f)
                    name=f[0:5]

                    if f.endswith('.txt'):
                        if os.path.split(dirpath)[-1]=="ascii (HI)":
                            if "spec" in os.path.split(fname)[-1]:
                                HI[name]=fname
                            elif "abs" in os.path.split(fname)[-1]:
                                HIabs[name]=spec_parser.LineDump(fname)
                        if os.path.split(dirpath)[-1]=="ascii (only HI)":
                            if "spec" in os.path.split(fname)[-1]:
                                HIonly[name]=fname
                            elif "abs" in os.path.split(fname)[-1]:
                                HIonly_abs[name]=spec_parser.LineDump(fname)
                        elif os.path.split(dirpath)[-1]=="ascii (only Metal)":
                            if "spec" in os.path.split(fname)[-1]:
                                metals[name]=fname 
                            elif "abs" in os.path.split(fname)[-1]:
                                metalabs[name]=spec_parser.LineDump(fname)
                    elif f.endswith('.xml'):
                        if "HI" in fname and not "only" in fname:
                            print(fname)
                            models[name]=Model(xmlfile=fname)
                            models[name].read()
                    else:
                        pass
                

        for key in models.keys():
            try:
                dct[key]={"HI_sp":spec_parser.Spectrum.sniffer(HI[key]),
                          "HIonly_sp":spec_parser.Spectrum.sniffer(HIonly[key]),
                          "metal_sp":spec_parser.Spectrum.sniffer(metals[key]), 
                          "HI":HIabs[key],
                          "HIonly":HIabs[key],
                          "metal":metalabs[key], 
                          "model":models[key]}
            except KeyError:
                print("skipping %s"%key)
        return [DumpData(**dct[name]) for name in dct.keys()]

    def set_SNR(self,vel=200.,ab_type='HI',verbose=True):
        def _filter(line):
            sp=getattr(self,ab_type+"_sp")
            if sp.waves[0]>item.obs_wave or line.obs_wave>sp.waves[-1]:
                return 0.
            ref=line.obs_wave

            cut=sp.get_cut(start=ref-vel*ref/c,end=ref+vel*ref/c)
            #if mostly absorption, i.e. within saturated HI line
            if not ab_type is 'HI':
                if np.mean(cut.abs)<0.2*np.mean(cut.cont):  
                    return 0.
            
            return np.median(cut.cont)/np.median(cut.error)
 

        line_dump=getattr(self,ab_type)
        lines=[] 
        if verbose: print("getting SNR for %s"%(self.name))
        for ab in line_dump.absorbers:
            try:
                for item in ab.get_lines():
                    lines.append(_filter(item))
            except KeyError:
                print("offending sp: %s"%(line_dump.fname))
                print(str(ab))
                raise

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

    @staticmethod
    def get_DA(spec, zem,  exclude_ab=False, rest_range=None):
        """
        Get DA for a spec.  if exlcude, can exclude contribution from what you 
        specify by getting total DA and then adding suspected absorption from 
        whatever class of absorbers.

        input:
        ------
        spec:  Spectrum subclass object instance
        rest_range: rest wavelength region to consider.
        exclude: either `metal` or `HI` or None
                 if None, get simple DA from total spec.

        output:
        -------
        da:  mean(1-f[i]/c[i]) where f is flux at pixel i and c is the continuum 
             level
        """

        indices=None
        if rest_range:
            rest_range=((1.+zem)*rest_range[0], (1.+zem)*rest_range[1])
            indices=spec_parser.Spectrum.get_indices(spec.waves, rest_range)

        flux, ab, cont = tuple(map(ma.masked_invalid,(spec.flux[indices], spec.abs[indices], spec.cont[indices])))
       
        
        if exclude_ab:
            da=(1.0-(ab-flux)/cont)
        else:
            da=(1.0-flux/cont)     

        da=ma.masked_invalid(da)
        da=ma.masked_where(np.fabs(da)>10.,da)
        #raise Exception("break")
        return np.mean(da)

def plot_hist(dat,attr, bin_rng=None, *args, **kwargs):
    
    if not 'ax' in kwargs.keys():
        fig, ax = plt.subplots()
    else:
        ax=kwargs.get("ax")
    if type(dat) is np.ndarray:
        dat=dat.tolist()
    assert(type(dat) is list)
    num_bins=kwargs.get('num_bins',15)
    dat=[getattr(item,attr) for item in dat]
    n, bins, patches = ax.hist(dat, 
                                 num_bins, normed=1, facecolor='green', 
                                 alpha=0.75)
    y = mlab.normpdf( bins, np.mean(dat), np.std(dat))
    if attr=="N":
        attr=r"log N"
    ax.plot(bins, y, 'r--', linewidth=1)
    ax.set_xlabel(kwargs.get('xlabel', attr))

    ab_type=kwargs.get("absorber_type","")

    ax.set_title(r'%s $\langle $%s$ \rangle=%3.2lf \pm %3.2lf$ ($n=%d$)'%(ab_type,attr,np.mean(dat), np.std(dat),len(dat)))

    if kwargs.get("ylabel",None):
        ax.set_ylabel(kwargs.get("ylabel"))
    elif bin_rng:
        ax.set_ylabel(r"%s from z=%2.1lf to %2.1lf"%(attr, bin_rng[0],bin_rng[-1]))
    #plt.show()

    return ax

def get_DA(spec, rest_range=(1090,1170), exclude=None):
    if exclude=="HI":
        sp=spec.metal_sp
    elif exclude=="metals":
        sp=spec.HIonly_sp
    else:
        sp=spec.HI_sp #don't care which sp to use
    exclude_ab=True if exclude else False
    return DumpData.get_DA(sp, spec.zem, exclude_ab=exclude_ab, rest_range=rest_range)


def flux_pdf(spec,rest_range=(1090,1170),**kwargs):
    def prep_flux(spec):
        flux=spec.HI_sp.get_cut(
                                start=(1.+spec.zem)*rest_range[0], 
                                end=  (1.+spec.zem)*rest_range[1]
                                ).flux
        cont=spec.HI_sp.get_cut(
                                start=(1.+spec.zem)*rest_range[0], 
                                end=  (1.+spec.zem)*rest_range[1]).cont
        flux=ma.masked_invalid(flux/cont)
        flux=ma.masked_where(flux>1.4, flux)
        flux=ma.masked_where(flux<-0.2,flux)
        return flux

    if type(spec) is list:
        for item in spec:
            flux=prep_flux(item)
            try:
                sp=np.concatenate(sp,flux)
            except:
                sp=flux

    else:
        flux=prep_flux(spec)
    
    ax =kwargs.get('ax', None)
    fig=kwargs.get('fig', None)
    if not ax:
        fig, ax=plt.subplots()
        close=True
    num_bins=kwargs.get('num_bins',40)
    n, bins, patches = ax.hist(flux, 
                                 num_bins, normed=True, facecolor='green', 
                                 alpha=0.75)
    fig.canvas.draw()
    #y = mlab.normpdf( bins, np.mean(sp.flux), np.std(sp.flux))
    #l = plt.plot(bins, y, 'r--', linewidth=1)
    ax.set_xlabel(kwargs.get(xlabel,r"normalised flux"))
    ax.set_ylabel(kwargs.get(ylabel,r"prob. density"))
    if close:
        plt.show()
        plt.clf()

    

    #plt.title(kwargs.get("title",""))
    #plt.show()
    return ax
        
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
    fig, axes = plt.subplots(len(data),2)
    for i in range(len(data)):
        sublst=data[i]
        #if len(sublst)>30:
        #for attr in ['N', 'b']:
        if len(sublst)<5:
            continue
        axes[i][0]=plot_hist(sublst,'N', bin_rng=[b,b+binsize],*args, ax=axes[i][0], **kwargs)
        axes[i][1]=plot_hist(sublst,'b', bin_rng=[b,b+binsize],*args, ax=axes[i][1], **kwargs)
        b+=binsize

    plt.show()
    z=[np.mean([item.z for item in it]) for it in data if len(it)>30]
    mu=[np.mean([item.N for item in it]) for it in data if len(it)>30]
    stdev=[np.std([item.N for item in it]) for it in data if len(it)>30]
    num_lines=[len(it) for it in data if len(it)>30]

    z,mu,stdev,num_lines=tuple(map(np.array, (z,mu,stdev,num_lines)))

    stderr=stdev/np.sqrt(num_lines)

    plt.clf()
    #plt.errorbar(z,mu,yerr=[stdev[i]/num_lines[i] for i in range(len(num_lines))], fmt='o')
    absorber_type=r"\textrm{H}\,\textsc{i}" if absorber_type=="HI" else "metal"
    plt.errorbar(z,mu,yerr=stderr, fmt='o')    
    plt.xlabel('mean z')
    plt.ylabel('mean N (%s)'%(absorber_type))
    plt.show()

#TODO: scale bins for comoving volume  
#TODO: scatter binned N vs z and binned b vs z for met and HI 4 plots
#TODO: same as previous, but split by ionization:  mgII, CIV

if __name__=="__main__":
    """
    spec=DumpData.factory('/home/scott/research/lyaf_stats')

    for key, val in names.items():
        for item in spec:
            if key in item.name:
                item.name=key
                item.zem=val
                continue
        
    for item in spec:
        item.set_SNR(ab_type='HI')
        item.set_SNR(ab_type='metal')

    pickle.dump(spec,open("../saved_objects/spec_dump.obj",'wb'))

    
    #flux_pdf([sp for sp in spec if hasattr(sp,'zem')],title="all spectra")


    for item in spec:
        item.da=get_DA(item)
        item.metal_da=get_DA(item,exclude="HI")
        item.HI_da=get_DA(item,exclude="metals")
        print(item.name,item.da,item.metal_da, item.HI_da)
    """
    spec=pickle.load(open("../saved_objects/spec_dump.obj",'rb'))

    flux_pdf([sp for sp in spec if hasattr(sp,'zem')],title="")
    met, HI,m1=[],[],[]

    print(len([sp for sp in spec if hasattr(sp,'zem')]))
    
    raise Exception("breakoint")  
    for item in spec:
        if hasattr(item,'zem'):
            flux_pdf(item,ylabel=item.name)
            met+=[it for it in item.metal.absorbers if it.ionName!="HI"]
            HI+=[it for it in item.HI.absorbers if it.ionName=="HI"]
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
    for lst in [met, HI]:
        #fig, axes=plt.subplots(len(zbins),len(bins))
        for i in range(len(zbins)):
            if i>len(lst): continue
            zbin=lst[i]
            for j in range(len(bins)):
                if j>len(zbin): continue
                SNRbin=zbin[j]
                #if len(SNRbin)<30:                 
                #    continue
                try:
                    
                    """plot_hist(SNRbin,"N", 
                            xlabel=r"%s: z=$%3.1lf--%3.1lf$, SNR=$%4.1lf--%4.1lf$"%(
                              "N",zbins[i][0],zbins[i][1],bins[j][0],bins[j][1]
                            ),
                            ylabel="metal" if lst is met else r"\textrm{H}\,\textsc{i}",
                            ax=axes[i][j])"""
                except:
                    pass
                if lst is met:
                    mean_metals.append((zbins[i][0], bins[j][0], SNRbin))
                elif lst is HI:
                    mean_HI.append((zbins[i][0], bins[j][0], SNRbin))
        #fig.show()

    #i= zbin, j=snrbin
    

    all_met, all_HI=[],[]
    for sp in spec:
        all_met+=sp.metal.absorbers
        all_HI+=sp.HI.absorbers

    xx=[(item[0]+item[1])/2. for item in bins]
    met=DumpData.bin_absorbers(all_met,bin_attr='SNR',bins=bins)
    HI=DumpData.bin_absorbers(all_HI,bin_attr='SNR',bins=bins)

    print(len(met), len(HI))
    
    x,y=[],[]
    for i in range(len(met)):
        if len(met[i])>1:
            y.append(met[i])
            x.append(xx[i])

    plt.errorbar(x,
                 [np.mean([ab.N for ab in it]) for it in y],
                 yerr=[np.std([ab.N for ab in it]) for it in y], 
                 fmt='o') 
    plt.title(r"N vs SNR: metals")
    plt.xlabel(r"SNR")
    plt.ylabel(r"log N(metals)")
    plt.show() 


    x,y=[],[]
    for i in range(len(HI)):
        if len(met[i])>1:
            y.append(met[i])
            x.append(xx[i])

    plt.errorbar(x,
                 [np.mean([ab.N for ab in it]) for it in y],
                 yerr=[np.std([ab.N for ab in it]) for it in y], 
                 fmt='o')  
    plt.xlabel(r"SNR")
    plt.ylabel(r"log N(\textrm{H}\,\textsc{i})")
    plt.title(r"N vs SNR: \textrm{H}\,\textsc{i}")
    plt.show() 

    
    histogram_by_z(spec, absorber_type="HI")
    histogram_by_z(spec, absorber_type="metal")

"""
    plt.xcorr()
    plt.acorr()
"""



 
        
