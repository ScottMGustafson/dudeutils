import dudeutils
from spec_parser import *
import matplotlib.pyplot as plt
from model import *
import os

def get_lines(lst, dct):
    ablst=[]
    for key, transition in dct.items():
        ab_line=lst.get_item(key).get_lines()[transition]
        ablst.append(ab_line)
    return ablst

def plot_line(spec, model, fig, subplot_key, line, velocity=True, ax=None,**kwargs):
    sharex, sharey= kwargs.get("sharex",None), kwargs.get("sharey",None)
    if not ax:
        ax=fig.add_subplot(subplot_key,sharex=sharex, sharey=sharey)
        assert(ax)
    waves, ab, cont, flux, error = get_data(spec,model,absorbers=[line],**kwargs)

    ax.plot(waves, flux, 'k', linestyle='steps')
    ax.plot(waves, error, 'g', linestyle='steps')
    ax.plot(waves, cont, 'r', linestyle='steps')
    ax.plot(waves, ab, 'b', linestyle='steps')
    ax.set_xlim(kwargs.get("xlims"),None)
    ax.set_ylim(kwargs.get("ylims"),None)
    return ax

def get_data(spec,model,wave_r=None, vel_r=None, absorbers=None, ref_ab=None,**kwargs):
    if not hasattr(spec,'abs') or absorbers:        
        if not type(absorbers) is list:
            absorbers=[absorbers]
        if type(absorbers[0]) is dict:
            # given list of one-element dicts, get single spectral line
            # list of dicts:
            # [ {abs_id (str): transition number (int)} ]
            # 
            tmp=[]
            lst=Model.get(model.AbsorberList)
            for i in range(len(absorbers)): 
                tmp+=get_lines(lst, absorbers[i])
            absorbers=tmp
        ab, cont, chi2= Spectrum.fit_absorption(spec,
                                                model,
                                                ab_to_plot=absorbers)        
    else:
        ab, cont, error = spec.abs, spec.cont, spec.error



    if vel_r and ref_ab:
        if not type(ref_ab) is SpectralLine:
            ref_ab=get_lines(Model.get(model.AbsorberList), ref_ab)[0]
            
        ind=Spectrum.get_indices(
                Spectrum.convert_to_vel(spec.waves,ref_ab.obs_wave), 
                vel_r)
    elif wave_r:
        ind=Spectrum.get_indices(spec.waves, wave_r)
    else:
        ind=np.arange(0,spec.waves.shape[0])

    ab, cont, waves, flux =ab[ind], cont[ind], spec.waves[ind], spec.flux[ind]
    if spec.error:
        error=spec.error[ind]
    else:
        error=np.zeros(len(ind))
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
    #this is an example for plotting a few lyman alpha transitions plus 
    #an intervening line

    xmlfile,  vel_r = "test.xml", [-50., 50.]
    model, spec=parse_spectrum( Model(xmlfile=xmlfile) )
    absorbers=[{"test_ab",2},{"HI_1":0},{"HI_1":1},{"HI_1":5},{"HI_1":10}]
    fig = plt.figure()

    axes=[]

    axes.append(plot_line(spec, model,
                        fig, subplot_key, 
                        model.get_spectral_line("HI_1",0), 
                        velocity=True, sharex=None, 
                        ref_ab=model.get_spectral_line("HI_1",0)))


    sharex=axes[0]

    axes[0]=plot_line(spec, model,  #to append another line to this ax
                        fig, subplot_key, 
                        model.get_spectral_line("test_ab",2), 
                        velocity=True, ax=axes[0],   #to share the ax with previous
                        ref_ab=model.get_spectal_line("HI_1",0))

    subplot_key=100*len(absorbers)+12   #212, since loop starts with second plot
    for item in absorbers[2:]:
        axes.append(plot_line(spec, model,
                        fig, subplot_key, 
                        item, 
                        velocity=True, sharex=sharex, ref_ab=item))
        sharex=axes[0]
        subplot_key+=1  #aling them vertically
    plt.show()
    
    
