import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import _spectrum
from spec_parser import Spectrum

class Param(object):
    """
    class for the manipulation of individual absorption parameters, N,b and z"""
    def __init__(self,param_name,value, locked, error=0., parent=None, index=None, bounds=None):
        self.name=param_name
        self.absorber=parent
        self.locked=locked
        self.error=float(error)
        self.value=float(value)
        self.index=int(index)
        if bounds:
            self.guess=Param.random_initial_cond(bounds)
        else:
            self.guess=float(value)     

    @staticmethod
    def random_initial_cond(bounds):
        """
        provide randomized initial conditions given some set of bounds

        Input:
        ------
        bounds: length=2 list of floats

        Output:
        -------
        (float)

        Raises:
        -------
        None
        """
        return(bounds[1]-bounds[0]) * np.random.random_sample()+bounds[0]
            
    def __str__(self):
        return"%s=%lf, %sError=%lf, %sLocked=%s"%(
            self.name,self.value,
            self.name,self.error,
            self.name,str(self.locked).lower())
         

def consolidate(ranges):
    """
    consolidate list of ranges to combine overlapping ranges as one.

    Input:
    ------
    ranges:  list of ranges  (list of length=2 lists of floats)

    Output:
    -------
    result: list of consolidated ranges  (list of length=2 lists of floats)

    Raises:
    -------
    None

    """
    result = []
    current_start = -1
    current_stop = -1 

    for start, stop in sorted(ranges):
        if start > current_stop:
            # this segment starts after the last segment stops
            # just add a new segment
            result.append( (start, stop) )
            current_start, current_stop = start, stop
        else:
            # segments overlap, replace
            result[-1] = (current_start, stop)
            # current_start already guaranteed to be lower
            current_stop = max(current_stop, stop)

    return result
   

def fit_absorption(dat, model):
    """
    Input:
    ------
    dat : specParser.Spectrum instance
    model : model.Model instance

    Output:
    -------
    cont : the continuum of the spectrum
    absorption : absorption of the spectrum
    chi2 : chi-square calculated for spectral regions specified in model.RegionList 

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
                lambda x: dat.waves[4]<x.get_obs(x.z)<dat.waves[-4], 
                absorbers))

    regions=[(item.start, item.end) for item in list(model.get_lst("RegionList"))]
    regions=consolidate(regions) 
    starts=np.array([item[0] for item in regions],dtype=np.float)
    ends=np.array([item[1] for item in regions],dtype=np.float)

    assert(starts.shape[0]==ends.shape[0])
    assert(x.shape==y.shape)
    assert(dat.waves.shape==dat.flux.shape==dat.error.shape)
    
    N=np.array([item.N for item in absorbers], dtype=np.float)
    b=np.array([item.b for item in absorbers], dtype=np.float)
    z=np.array([item.z for item in absorbers], dtype=np.float)
    rest=np.array([item.wave for item in absorbers], dtype=np.float)
    gamma=np.array([item.gamma for item in absorbers], dtype=np.float)
    f=np.array([item.f for item in absorbers], dtype=np.float)


    cont, absorption, chi2= _spectrum.spectrum( dat.waves,dat.flux,dat.error,
                                                 x,y,
                                                 N,b,z,rest,gamma,f,
                                                 starts,ends )

    for i in range(len(absorbers)):
        model.set_val(absorbers[i],**{"N":N[i],"b":b[i],"z":z[i]}) #this propagates all the way down
    return cont, absorption, chi2



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

    regions=[(item.start, item.end) for item in list(model.get_lst("RegionList"))]
    regions=consolidate(regions) 
    starts=np.array([item[0] for item in regions],dtype=np.float)
    ends=np.array([item[1] for item in regions],dtype=np.float)

    #check that arrays are equaly sized
    assert(starts.shape==ends.shape)
    assert(x.shape==y.shape)
    assert(spec.waves.shape==spec.flux.shape==spec.error.shape)


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
    test_xml = "/home/scott/research/test.xml"

    model=dudeutils.get_model(test_xml)
    sp = spectrum.Spectrum.sniffer(model.flux, model.error)

    print("testing optimizer.fit_absorption")
    cont, ab, chi2 = optimizer.fit_absorption(sp, model)

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
