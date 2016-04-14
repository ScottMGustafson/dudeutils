#import dudeutils
from spec_parser import *
from atomic import *
import matplotlib.pyplot as plt
from model import *
import os
import numpy as np
import numpy.ma as ma



def get_lines(model, dct):
    return [model.get_spectral_line(key, val) for key, val in dct.items()]


def plot_line(spec, model, fig, subplot_key, line,**kwargs):

    velocity=kwargs.pop("velocity",True)
    ax=kwargs.pop("ax",None)
    sharex=kwargs.get("sharex",None)
    sharey= kwargs.get("sharey",None)
    multiplier=kwargs.pop('multiplier',None)
    xlabel=kwargs.pop('xlabel',None)
    topmost=kwargs.pop('topmost',False)
    ref=kwargs.pop('ref',None)


    wave_r, ind, ref = prep_data(spec,model,ref,**kwargs)
    absorbers=[item for item in line] if type(line) is list else [line]
    waves, flux, error, ab0, cont = get_ab(spec,model,ref, ind, absorbers[0],return_ab=False)
    ab=[ab0]
    if len(absorbers)>1:
        ab+=list(get_ab(spec,model,ref, ind, absorbers[1:],return_ab=True))

    if multiplier:
        def f(x): return multiplier*x
        cont,flux,error= tuple( multiplier*x for x in (cont, flux, error) )
        ab=[f(it) for it in ab]
        #ab=[multiplier*np.array(it) for it in ab]

    tol=kwargs.pop('tol',1E12)

    cont=ma.masked_outside(cont,-1.*tol,tol)
    flux=ma.masked_outside(flux,-1.*tol,tol)
    error=ma.masked_outside(error,-1.*tol,tol)

    assert(flux.shape[0]==ab[0].shape[0]==cont.shape[0]==waves.shape[0])

    if not ax:
        ax=fig.add_subplot(subplot_key,sharex=sharex, sharey=sharey)
        ax.plot(waves, flux, 'k', linestyle='steps-mid')
        ax.plot(waves, np.zeros(waves.shape[0]), 'k--')
        ax.plot(waves, cont, 'r--', linestyle='steps-mid')
        try:
            ax.plot(waves, error, 'g', linestyle='steps-mid')
        except:
            pass

    def plot_ab(ab):
        """
        plot absorption separate from continuum so we can plot the indidual 
        absorption from diff absorbers
        """
        if type(ab) is list:
            assert(type(ab[0]) is np.ndarray)
            for item in ab:
                plot_ab(item)
        else:
            if not np.array_equal(ab,cont):
                ab=ma.masked_outside(ab,-1.*tol,tol)
                ax.plot(waves, ab, 'b--') 
    


    plot_ab(ab)


    x,labels, y=[],[],[]
    import lineid_plot
    for item in prep_absorbers(spec, model, absorbers):
        for line in item:
            obs_wave=Spectrum.convert_to_vel(line.obs_wave,ref)
            if np.amin(waves)<obs_wave<np.amax(waves):
                x.append(obs_wave)
                labels.append(line.ionName)
                y.append(np.mean(cont))

                
    box_loc=lineid_plot.get_box_loc(fig, ax, x, y, box_axes_space=0.06)
    if topmost:
        for i in range(len(x)):  
            ax.annotate(line.ionName, xy=(x[i], y[i]),
            xytext=(box_loc[i][0],
                    box_loc[i][1]),
            xycoords="data", textcoords="data",
            rotation=90, horizontalalignment="center",
            verticalalignment="center",
            fontsize=12,
            arrowprops=dict(arrowstyle="-",
                            relpos=(0.5, 0.0)),
            label=line.ionName)
    for i in range(len(x)): 
        ax.axvline(x=x[i], ymin=0, ymax=y[i], linewidth=1, color='k',linestyle='--')

    if xlabel:
        ax.set_xlabel(xlabel)
    else:
        ax.axes.get_xaxis().set_visible(False)
        ax.set_xlabel(None)
    ax.set_xlim(kwargs.get("xlims"),None)
    ax.set_ylim(kwargs.get("ylims"),None)
    return ax
        
def parse_spectrum(model):
    error=None
    if type(model) is str:
        if model.endswith('.xml'):
            model=Model(xmlfile=model)
    else:
        assert(type(model) is Model)

    try:
        spec, error =model.flux, model.error
    except:
        spec=model.flux
    if not os.path.sep in spec:
        pth, f_ = os.path.split(model.id)
        assert(f_==spec)
        spec=os.path.join(pth, spec)
        error=os.path.join(pth, error)


#TODO: error lies somewhere fit_absorption returns flux (and presumably also error)
#very subtle bug.
    #spec=Spectrum.sniffer(spec, error=error)
    sp_dump='test_SiII_out.dat'
    spec=Spectrum.sniffer(sp_dump)

    #if not hasattr(spec, 'abs')
    #    spec.cont, spec.abs, chi2 = Spectrum.fit_absorption(spec,model)
    return model, spec

def prep_absorbers(spec, model, absorbers):
    if type(absorbers) is list:
        return [prep_absorbers(spec, model, ab) for ab in absorbers]
    if not absorbers:
        return []
    
    if type(absorbers) is dict:
        absorbers=get_lines(model, absorbers)
        #print(str(absorbers[-1]),str(absorbers[-1].z))
    elif type(absorbers) is SpectralLine:
        pass
    else:
        raise TypeError("type of absorbers must either be list of "+\
                        "SpectralLine instances or dict")
    return absorbers


def get_wr_from_vr(ref,model,vel_r):
    return [(1+vel_r[0]/c)*ref, (1+vel_r[1]/c)*ref]

def prep_data(spec,model,ref=None,wave_r=None, vel_r=None,**kwargs):

    if ref:
        if type(ref) is dict:
            assert(len(ref.keys())==1)  #should only be one element long
            key=list(ref.keys())[0]
            ref=model.get_spectral_line(key, ref[key]).obs_wave
        elif type(ref) is SpectralLine:
            ref=model.get_spectral_line(key, ref[key]).obs_wave
        elif type(ref) is float:
            pass
        else:
            raise Exception("poorly defined input: ref")

    if vel_r:
        if wave_r or not ref:
            raise Exception("poorly defined input")
        else:
            wave_r=get_wr_from_vr(ref,model,vel_r)
 
    if wave_r:
        ind=Spectrum.get_indices(spec.waves, wave_r) 
    else:
        raise Exception("test fail")# np.arange(0,spec.waves.shape[0])


    if len(ind)==0:
        raise Exception("No data selected")

    return wave_r, ind, ref


def get_ab(spec,model,ref, ind, absorbers, vel_r=True, return_ab=True):

    if type(absorbers) is list:
        print('should be dict or spectral line',type(absorbers[0]), len(absorbers[0]))
        return [get_ab(spec,model,ref, ind, ab, vel_r, return_ab) for ab in absorbers]


    absorbers=prep_absorbers(spec, model, absorbers)
    waves, flux, error, ab, cont, chi2= Spectrum.fit_absorption(spec,
                                                        model,
                                                        ab_to_plot=absorbers,
                                                        indices=ind)  
    if type(spec) is TextSpectrum:  
        #needed because for some reason, fit_absorption is not giving right flux
        waves, flux, error = spec.waves[ind], spec.flux[ind], spec.error[ind] 

    if vel_r:
        waves=Spectrum.convert_to_vel(waves, ref)

    if return_ab:
        return ab
    else: 
        return waves, flux, error, ab, cont
    
def label_lines(ax, x, y, label):
    for i in range(len(label)):
        ax.annotate(label[i], xy=(x[i], y[i]), xytext=(x[i], 0.8),ha='center', va='center',
                    arrowprops=dict(arrowstyle="-",facecolor='black') )
        #plt.axvline(x=x[i], ymin=0., ymax=0.68, linewidth=1, color='k')
    return ax

    



if __name__=="__main__":
    #this is an example for plotting a few transitions of SiII

    xmlfile,  vel_r = "/home/scott/research/J0744+2059/SiII.xml", [-120., 120.]
    model, spec=parse_spectrum( Model(xmlfile=xmlfile) )
    fig = plt.figure()

    axes=[]

    axes.append(plot_line(  spec, model,
                            fig, 411,  
                            [{"SiII2":2},{"SiII1":2},{"SiII3":2}], 
                            vel_r=vel_r,velocity=True, ax=None,   
                            ref={"SiII2":2},
                            multiplier=10E14,topmost=True))
    
    axes.append(plot_line(  spec, model,
                            fig, 412,  
                            [{"SiII2":3},{"SiII1":3},{"SiII3":3}], 
                            vel_r=vel_r,velocity=True, ax=None,   
                            ref={"SiII2":3},
                            multiplier=10E14,
                            sharex=axes[0]))


    axes.append(plot_line(  spec, model,
                            fig, 413,  
                            [{"SiII1":5},{"SiII2":5},{"SiII3":5}], 
                            vel_r=vel_r,velocity=True, ax=None,   
                            ref={"SiII2":5},
                            multiplier=10E14,
                            sharex=axes[0]))

    axes.append(plot_line(  spec, model,
                            fig, 414,  
                            [{"SiII1":6},{"SiII2":6},{"SiII3":6}], 
                            vel_r=vel_r,velocity=True, ax=None,   
                            ref={"SiII2":6},
                            multiplier=10E14,
                            sharex=axes[0],
                            xlabel=r"velocity km s$^{-1}$"))



    fig.text(0.04, 0.5, 'Flux', va='center', rotation='vertical')

    plt.show()
    
