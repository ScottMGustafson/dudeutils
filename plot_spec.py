#import dudeutils
from spec_parser import *
from atomic import *
import matplotlib.pyplot as plt
from model import *
import os

def get_lines(model, dct):
    return [model.get_spectral_line(key, val) for key, val in dct.items()]


def plot_line(spec, model, fig, subplot_key, line, velocity=True, ax=None,**kwargs):
    sharex, sharey= kwargs.get("sharex",None), kwargs.get("sharey",None)
    waves, ab, cont, flux, error = get_data(spec,model,absorbers=line,**kwargs)
    assert(flux.shape[0]==ab.shape[0]==cont.shape[0]==waves.shape[0])
    if not ax:
        ax=fig.add_subplot(subplot_key,sharex=sharex, sharey=sharey)
        ax.plot(waves, flux, 'k', linestyle='steps-mid')
        if error: ax.plot(waves, error, 'g', linestyle='steps-mid')
        ax.plot(waves, cont, 'r', linestyle='steps-mid')
    ax.plot(waves, ab, 'b', linestyle='steps-mid') #this needs to be outside of if 
                                               #in case appending ab to 
                                               #pre-existing ax

    ax.set_xlim(kwargs.get("xlims"),None)
    ax.set_ylim(kwargs.get("ylims"),None)
    return ax

def get_data(spec,model,wave_r=None, vel_r=None, absorbers=None, ref_ab=None,**kwargs):
    if not hasattr(spec,'abs') or absorbers:       

        if type(absorbers) is dict:
            absorbers=get_lines(model, absorbers)
            print(str(absorbers[-1]))
        elif type(absorbers) is list:
            assert(type(absorbers[0]) is SpectralLine)
        else:
            raise TypeError("type of absorbers must either be list of SpectralLine instances or dict")
 
        ab, cont, chi2= Spectrum.fit_absorption(spec,
                                                model,
                                                ab_to_plot=absorbers)        
    else:
        ab, cont, error = spec.abs, spec.cont, spec.error

    if vel_r and ref_ab:
        if type(ref_ab) is dict:
            assert(len(ref_ab.keys())==1)  #should only be one element long
            key=list(ref_ab.keys())[0]
            ref_ab=model.get_spectral_line(key, ref_ab[key])
        else:
            assert(type(ref_ab) is SpectralLine)
            
        ind=Spectrum.get_indices(
                Spectrum.convert_to_vel(spec.waves,ref_ab.obs_wave), 
                vel_r) #get indices within velocity range
    elif wave_r:
        ind=Spectrum.get_indices(spec.waves, wave_r)
    else:
        ind=np.arange(0,spec.waves.shape[0])

    ab, cont, waves, flux =ab[ind], cont[ind], spec.waves[ind], spec.flux[ind]
    if spec.error:
        error=spec.error[ind]
    else:
        error=None
    if vel_r: 
        waves=Spectrum.convert_to_vel(waves, ref_ab.obs_wave)
       
    return waves, ab, cont, flux, error
        
def parse_spectrum(model):
    error=None
    if type(model) is str:
        if model.endswith('.xml'):
            model=dudeutils.get_model(model) 
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

    spec=Spectrum.sniffer(spec)
    if error:
        error=Spectrum.sniffer(error)
    #if not hasattr(spec, 'abs')
    #    spec.cont, spec.abs, chi2 = Spectrum.fit_absorption(spec,model)

    return model, spec

if __name__=="__main__":
    #this is an example for plotting a few transitions of SiII

    xmlfile,  vel_r = "/home/scott/research/J0744+2059/SiII.xml", [-120., 121.]
    model, spec=parse_spectrum( Model(xmlfile=xmlfile) )
    absorbers=[{"SiII3":2, "SiII1":2},{"SiII3":3, "SiII1":3},{"SiII3":4,"SiII1":4},{"SiII3":5},{"SiII3":6}]
    fig = plt.figure()

    axes=[]

    axes.append(plot_line(  spec, model,
                            fig, 411,  
                            absorbers[0], 
                            vel_r=vel_r,velocity=True, ax=None,   
                            ref_ab={"SiII3":absorbers[0]["SiII3"]}))


    sharex=axes[-1]
    
    axes.append(plot_line(  spec, model,
                            fig, 412,  
                            absorbers[1], 
                            vel_r=vel_r,velocity=True, ax=None,   
                            ref_ab={"SiII3":absorbers[1]["SiII3"]},
                            sharex=sharex))


    axes.append(plot_line(  spec, model,
                            fig, 413,  
                            absorbers[2], 
                            vel_r=vel_r,velocity=True, ax=None,   
                            ref_ab={"SiII3":absorbers[2]["SiII3"]},
                            sharex=sharex))

    axes.append(plot_line(  spec, model,
                            fig, 414,  
                            absorbers[3], 
                            vel_r=vel_r,velocity=True, ax=None,   
                            ref_ab={"SiII3":absorbers[3]["SiII3"]},
                            sharex=sharex))

    """axes.append(plot_line(  spec, model,
                            fig, 415,  
                            absorbers[4], 
                            vel_r=vel_r,velocity=True, ax=None,   
                            ref_ab={"SiII3":absorbers[4]["SiII3"]},
                            sharex=sharex))"""


    

    plt.show()
    
