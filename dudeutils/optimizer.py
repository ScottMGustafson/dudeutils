import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import _spectrum
from dudeutils.spec_parser import Spectrum
from dudeutils.data_types import *
from dudeutils.model import *



def is_big_enough(line, lya_min=11.25, f_lya=0.416):
    """
    Tests whether a given line is optically thick enough to matter for the fit
    by comparing oscillator strength and column density to lyman-alpha . The
    limit is 11.25 for lyman-alpha at f=0.416 by default.

    cwave = (1.0+line.z)*line.rest
    vdopp = line.b/cwave*1.0e13
    alpha = line.gamma/(4.0*np.pi*vdopp)/(1.0+line.z)
    factor = (10.0**line.N)*2.647E-2*line.f/(sqrt(np.pi)*vdopp)*1.0/(1.0+line.z)

    Input:
    ------
    line: data_types.Absober instance.
    lya_min: minimum column for lya to matter.  (float)
    f_lya: oscillator strength of lya (float)
    
    Output:
    -------
    bool: True if given line is optically thick enough to be included.

    Raises:
    -------
    None
    """

    return line.f*10.0**line.N > (10.0**lya_min) * f_lya

def is_locked(param, ab, starts, ends):
    """
    Input:
    ------
    param : Param instance
    ab: an absorber.  data_types.Absorber type

    Ouput:
    ------
    True if locked, false if absorber should be optimized
    """
    if not param.locked:
        lst=ab.get_lines()
        for line in lst:
            for j in range(starts.shape[0]):
                #if there is at least one transition that is both in an 
                #optimization region and has an absorptive enough line
                if starts[j]<line.obs_wave<ends[j] and is_big_enough(line):
                    return False
    return True


def optimize(spec, model):

    """
    calculate best fit parameters using scipy.optimize.curve_fit

    Input:
    ------
    spec : observed_data.observedData instance

    Output:
    -------
    results from scipy.optimize.curve_fit.  uses levenberg marquardt optimization

    Raises:
    -------
    AssertionError
    """

    cont_points=model.get_lst("ContinuumPointList")
    cont_points = sorted(cont_points, key=lambda pt: pt.x)
    x=np.array([float(item.x) for item in cont_points])
    y=np.array([float(item.y) for item in cont_points])

    absorbers=[]
    for item in model.get_lst("AbsorberList"):
        absorbers+=item.get_lines()

    absorbers=list(filter(
                   lambda x: 
                        spec.waves[4]<x.get_obs(x.z)<spec.waves[-4], 
                        absorbers))
    regions=Model.get(model.RegionList)
    starts=np.array([item.start for item in regions],dtype=np.float)
    ends=np.array([item.end for item in regions],dtype=np.float)
    
    #convert into numpy array for c extension use


    N=np.array([i.N for i in absorbers], dtype=np.float)
    b=np.array([i.N for i in absorbers], dtype=np.float)
    z=np.array([i.N for i in absorbers], dtype=np.float)
    wave=np.array([i.wave for i in absorbers], dtype=np.float)
    gamma=np.array([i.gamma for i in absorbers], dtype=np.float)
    f=np.array([i.f for i in absorbers], dtype=np.float)




    #get indices relevant for optimization
    indices=[]
    for r in range(starts.shape[0]):
        indices+=Spectrum.get_indices(spec.waves,[starts[r],ends[r]])   

    assert(len(indices)<spec.waves.shape[0])

    abs_lst = list(model.get_lst('AbsorberList'))

    param_lst=[]
    for i in range(len(abs_lst)):
        a=abs_lst[i]
        #param_name,value, locked, error, parent_id,index=None)
        Nargs=["N",a.N,a.NLocked,a.NError,a,i]
        bargs=["b",a.b,a.bLocked,a.bError,a,i]
        zargs=["z",a.z,a.zLocked,a.zError,a,i]

        param_lst+=[Param(*Nargs), Param(*bargs), Param(*zargs)]
        

    locked, unlocked = [],[]
    for item in param_lst:
        if is_locked(item, abs_lst[i],starts,ends):
            locked.append(item)
        else:
            unlocked.append(item)

    all_params=list(unlocked)+list(locked)
    #values of the params to vary this gets passed as an argument to the 
    #function to optimize
    params = [item.value for item in unlocked] 
    #unlocked_indices = [item.index for item in unlocked] #index mapping to all_params

    absorbers = list(set([item.absorber for item in all_params]))

    assert(spec.waves[indices].shape[0]<spec.waves.shape[0])

    def absorption(waves, *params):

        """
        this is our function to pass into scipy.optimize

        calculate the best fit for given list of *params.  Needs to be in scope 
        of optimize so that it can inherit some of its objects, while keeping 
        inputs consistent with scipy.optimize.curve_fit

        Input:
        ------
        waves : 1d numpy array of wavelengths 
        *params : parameters to be optimized.  determined in the parent function

        Output:
        -------
        chi2 : chi-square calculated from regions specified in model.RegionList

        Raises:
        -------
        AssertionError
        """
     
        #update params to locked
        for i in range(len(unlocked)):
            unlocked[i].value=params[i]  #set param value
            all_params[unlocked[i].index] = unlocked[i]  #write back to all_params
        #now consolidate locked and unlocked parameters into lists for N,b,z

        _N,_b,_z=[],[],[]
        rest,gamma,f=[],[],[]
        for ab in absorbers:
            lines=ab.get_lines()
            for line in lines:
                if waves[4]<line.get_obs()<waves[-4]:
                    _N.append(line.N)
                    _b.append(line.b)
                    _z.append(line.z)
                    rest.append(float(line.wave))
                    gamma.append(float(line.gamma))
                    f.append(float(line.f))
        N,b,z,rest,gamma,f= tuple(map(np.array,(_N,_b,_z,rest,gamma,f)))
        #there is an entry of N, b, z for each individual absorption line
        #so if an absorber has multiple lines in atom.dat, the number will 
        #appear multiple times

        cont, absorption_sp, chi2= _spectrum.spectrum(spec.waves,
                                                     spec.flux,
                                                     spec.error,
                                                     x,y,
                                                     N,b,z,rest,gamma,f,
                                                     starts,ends)

        assert(absorption_sp[indices].shape[0] == waves.shape[0])
        return absorption_sp[indices]

    p0=[item.guess for item in unlocked]

    popt, pcov = curve_fit(absorption, 
                           spec.waves[indices], 
                           spec.flux[indices], 
                           p0=p0, 
                           sigma=spec.error[indices])
    assert(len(popt)==len(unlocked))

    
    j=0
    for i in [item.index for item in unlocked]:
        all_params[i].value=popt[j] 
        all_params[i].error=np.sqrt(pcov[j][j]) #one stdev
        j+=1

    #changes should also propagate to unlocked, bc all_params inlcudes unlocked

    for item in all_params:
        model.set_val(item.absorber,**{item.name:item.value, item.name+"Error":item.error})

    return unlocked, popt, pcov
        
    #rewrite new values to model

if __name__ == "__main__":
    import dudeutils
    test_xml = "test.xml"

    model=dudeutils.get_model(test_xml)
    sp = Spectrum.sniffer(model.flux, model.error)

    print("testing optimizer.fit_absorption")
    cont, ab, chi2 = Spectrum.fit_absorption(sp, model)

    continuum_points = sorted(model.get_lst('ContinuumPointList'),key=lambda u:u.x)
    
    #split x,y into lists
    x,y = map(list,zip(*[(item.x, item.y) for item in continuum_points]))

    plt.plot(src_data.waves,src_data.flux,"-k",linestyle='steps')
    plt.plot(x,y,"co")
    plt.plot(src_data.waves,cont,"-b",linestyle='steps')
    plt.plot(src_data.waves,ab,"-r",linestyle='steps')
    plt.ylim([-1E-14,1E-13])
    plt.show()

    print("optimizing spec")

    params, popt, pcov = optimizer.optimize(sp, model)

    print(str(popt))
    print("\n"+str(pcov))
