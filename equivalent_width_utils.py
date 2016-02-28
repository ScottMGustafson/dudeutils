from data_types import *
from model import *
import numpy as np
import spec_parser
from scipy.constants import pi,e,c,Planck,Boltzmann,epsilon_0,m_e

def measured_area(mx,mn,fwhm,fwhm_inst=6.5):
    #units = km/s   
    height=mx-mn
    sigma= fwhm/2.35
    return height*sigma*np.sqrt(2.*constants.pi)/mx


def expected_area(ab, trans, continuum_level):
    """
    given some absorber, calculated the expected absorption as 
    equivalent_width*continuum


    input:
    ------
    ab: data_types.Absorber instance
    trans: transition to consider.  0 is n=1
    continuum_level: mean continuum level about line center (float)

    output:
    -------
    absorption (float)

    raises:
    TypeError when trans is incorrectly specified

    """

    #get the appropriate line transition
    atom_name=SpectralLine.get_atom_name(ab.ionName)
    lines=ab.get_lines()
    if type(trans) is float:            
        for item in lines: 
            if trans==item.wave:
                line=item
                break
    elif type(trans) is int:
        line=lines[i]
    else:
        raise TypeError("ambiguous declaration of transition: %s"%(str(trans)))
        
    #given the transition, calculate the expected ew and the area of absorption
    
    return equivalent_width(  ab.N, line.f, line.obs_wave,
                                           units="pixels",
                                           )*continuum_level


def equivalent_width(logN,f,lmda,optically_thin=True,
                     units="Angstroms", vdisp=2.14):
    """
    calculate the equivalent width of an optically thin spectral line.

    input:
    ------
    logN (float):  log10 column density
    f (float): oscillator strength, i.e. a factor to account for quantum effects
    lmda (float): wavelength of line-center in Angstroms
    optically_thin (bool): is the line optically thin?  if false raises exception
    units (str):  units in which to  return EW.  accepts either "Angstroms" 
                  or "pixels".
    vdisp (float):  pixel width in km/s

    output:
    -------
    equivalent width (float)   

    """

    if not optically_thin:
        raise Exception("optically thick absorption not yet supported")
    lmda/=10.**10. #convert from angstroms to meters
    a=lmda**2. * e**2./(4.*epsilon_0*m_e*c*c)  #conglomerate of constants

    #when inputing N in cm-2, mass in kg, lambda in A, get out 10E-16 meters.
    #convert to Angstroms by mulitplying by 10.**6
    W=10.**6. * 10.**(logN)*f*a  
    if units=='Angstrom':
        return W
    if units=="pixels":
        #v/c = delL/L.  delL=W, L=central wave
        return (0.001*c/vdisp)*(W/lmda)  



def logN_from_ew(ew,f,lmda,optically_thin=True,
                     ew_units="Angstroms", vdisp=2.14):
    """
    calculate the logN from equivalent width of an optically thin spectral line.

    input:
    ------
    ew (float): equivalent width
    f (float): oscillator strength, i.e. a factor to account for quantum effects
    lmda (float): wavelength of line-center in Angstroms
    optically_thin (bool): is the line optically thin?  if false raises exception
    units (str):  units in which to  return EW.  accepts either "Angstroms" 
                  or "pixels".
    vdisp (float):  pixel width in km/s

    output:
    -------
    equivalent width (float)   


    raises:
    -------
    Exception: if optically_thick is True, since not yet supported
    """

    if not optically_thin:
        raise Exception("optically thick absorption not yet supported")
    lmda/=10.**10. 
    a=lmda**2. * e**2./(4.*epsilon_0*m_e*c*c)

    #when inputing N in cm-2, mass in kg, lambda in A, get out 10E-16 meters.
    #convert to Angstroms by mulitplying by 10.**6
    #W=10.**6. * 10.**(logN)*f*a 


    if ew_units=='Angstrom':
        return np.log10(W) - 6. -np.log10(f*a)

    elif ew_units=="pixels":
        def pixels_to_angstroms(x):
            return x*0.001*vdisp*lmda/c
        ew = pixels_to_angstroms(ew)
        return np.log10(ew)-np.log10(f*a)-6. 

def flux_to_counts(SNR, continuum_level, flux):
    """
    get an estimate of number of counts

    input:
    ------
    SNR: (float) signal-to-noise ratio
    continuum_level: (float) mean continuum level
    flux: (float, or list of floats or np.array of floats) observed flux

    output:
    -------
    photon counts (float)

    """
    return SNR**2. * flux/continuum_level

def abs_to_counts(SNR, continuum_level, flux):
    """
    get expected absorption
    
    input:
    ------
    SNR: (float) signal-to-noise ratio
    continuum_level: (float) mean continuum level
    flux: (float, or list of floats or np.array of floats) observed flux

    output:
    -------
    absortion in photon counts (float)

    raises:
    -------
    None
    """
    return SNR**2. * (1.-flux/continuum_level)

def mean_counts(SNR, continuum_level, flux):
    """ 
    input:
    ------
    SNR (float)
    continuum_level (float)  mean continuum level
    flux (list)  observed flux values

    output:
    -------
    float

    """
    return abs_to_counts(SNR, continuum_level, np.mean(flux))
           
def missing_counts(ab,trans,continuum_level, SNR, flux):
    """
    get missing counts assuming there is no absorption
    
    input:
    ------
    ab: data_types.Absorber instance
    trans: transition to consider.  0 is n=1
    SNR: (float) signal-to-noise ratio
    continuum_level: (float) mean continuum level
    flux: (float, or list of floats or np.array of floats) observed flux

    output:
    -------
    absortion in photon counts (float)

    raises:
    -------
    None
    """
    
    n_obs=0.
    for item in flux:
        n_obs+=abs_to_counts(SNR, continuum_level, item)

def counts(ab,trans,continuum_level, SNR, flux):
    n_obs=0.
    for item in flux:
        n_obs+=flux_to_counts(SNR, continuum_level, item)
    n_exp = expected_area(ab, trans, continuum_level)
    uncert = np.sqrt(mean_counts(SNR, continuum_level, flux))
    return n_obs-n_exp, uncert

def maxN(lmda,f,m,n_obs,mean_n,vdisp=2.14):

    lmda=lmda*10.**(-10.) #convert to meters  


    #W = (N/(0.01**2.))*f*e*e*lmda*lmda/(4*epsilon_0*m_e*c*c)   #ew in meters

    # W = (0.001*c/vdisp)*W/lmda   #  W in pixels
    # n_exp = W*SNR**2.  #expected pixels missing  SNR**2 = mean counts in unabsorbed region

    # H0:  W*(SNR**2) <2sigma   # assume nothing there

    #if W*SNR**2 > 2sigma(n_obs), somthing there.
    #  W(in pixels) = (0.001*c/vdisp)*(N/(0.01**2.))*f*e*e*lmda/(4*epsilon_0*m_e*c*c)


    #W<2sigma/SNR**2
    #Nmax=(2sigma/SNR**2)*(1000.*vdisp)*(0.001**2.)*(4*espilon_0*m_e*c)/(f*e*e*lmda)
    #    =(2./SNR)*(1000.*vdisp)*(0.001**2.)*(4*espilon_0*m_e*c)/(f*e*e*lmda) 
    #n_obs = flux_to_counts(continuum_level)

    #calculate total pixels observed:  n_exp/npix=SNR**2
    #                                  n_obs=SNR**2 * sum(Fi/cont)
    #                                  n_absorbed-n_exp-n_obs
    
    #
                                        
    #convert to an EW:                 EW =  n_absorbed/SNR**2
    #                                     =(nexp-n_obs)/SNR**2
    #                                     =(npix-n_obs/SNR**2)
    #
    #                                    EW=SNR**2-n_obs/npix
    #                                      =SNR**2(1-1/npix * sum(flux/cont))
    #  EW=total_counts_absorbed/mean_counts
    # in npix bins,  nexp expected counts, n_obs observed counts, nabs absorbed
    # assume our absorber is only one
    #  so in any num pix, get n_abs absorbed pix
    #  so npix doesn't matter for n_absorbed.  
    #   EW = total_counts_absorbed/mean_n
    #      = (nexp-n_obs)/mean_n

    #   N = 0.01**2.  * px * (4*epsilon_0*mc*c)/(f*e*e*lmda) * (npix-n_obs/SNR**2)
    #physical constants 

   # EW=mean_n - n_obs/len(flux)


    a=(4.*epsilon_0*m_e*c**2.)/(e*e*lmda*f)  
    #conversion to pixels:
    conversion_factors=0.01**2. * 1000.*vdisp/c

    maxN=   conversion_factors *a* (n_obs/mean_n) #=const * n_obs/mean_n * 1/W(N)
#convert to common log
    #return np.log10(maxN), np.sqrt(mean_n)/(n_obs*np.log(10.))    #logN, uncert in logN
    return np.log10(maxN), np.log10(conversion_factors*a/(np.sqrt(mean_n)))

def wave_to_pixels(wave,wave0,vdisp=2.14):
    #EW/L=1000.*v(km/s)/c(m/s)
    #-->  v = (EW/L) * 0.001*c 
    #v(pixels) = v*pixels/vdisp = (EW/L) * 0.001*c /vdisp
    return (wave/wave0) * 0.001*c/vdisp

def pixels_to_wave(npix,wave0,vdisp=2.14):
    #pixels = (EW/L) * 0.001*c /vdisp
    
    return  npix*wave0*1000.*vdisp/c

def get_maxN(ab,trans,continuum_level, SNR, flux, err, nsigma=2.,vdisp=2.14):


    n_obs=0.
    for item in flux:
        n_obs+=flux_to_counts(SNR, continuum_level, item)

    uncert=np.sqrt(flux_to_counts(SNR, continuum_level, np.mean(flux)))

    line=ab.get_lines()[trans]

    lmda=ab.get_wave(trans)*10.**(-10.) #convert from angstroms to meters
    rest_lmda=ab.get_rest_wave(trans)*10.**(-10.)
    f=line.f
    mean_n=SNR**2.
    print("SNR=%lf, mean=%lf, cont=%lf"%(SNR, SNR**2., cont*10.**(14.)))
    #physical constants and parameters
    a=(4.*epsilon_0*m_e*c*c)/(f*rest_lmda*rest_lmda*e*e)

    #conversion from pixels to wave, m-2 to cm-2
    const = rest_lmda*vdisp/(10.*c)

    #expected n if no absorption at all.  if less than n_obs, then there is emission 
    #or our continuum was poorly chosen
    n_exp=float(len(flux))*mean_n
    
    #equivalent width in pixels as function of n_exp and n_obs
    W=np.fabs(n_exp-n_obs)/mean_n
    #dW=W*nsigma/np.sqrt(mean_n)

    #dW1 = 2.*np.sqrt(sum( [(it/continuum_level)**2. for it in err] ))#/np.sqrt(len(err))

    dW = 2.*np.std(np.array(flux))/continuum_level
    print("ew, d(ew)=%lf+/-%lf pixels"%(W,dW))
    print("ew=%lf+/-%lf Angstroms"%(10.**10 * W*rest_lmda*1000.*vdisp/c,10.**10 * dW*rest_lmda*1000.*vdisp/c))


    N, dN = np.log10(a*const*W), dW/(W*np.log(10.))

    print("logN, dlogN",N,dN)
    
    #max N + nsigma standard deviations
    return N+2.*dN


if __name__ == "__main__":

    
    dct = {
            "CI":[ {'start':5299.5 , 'end':5300.5, 'trans':2},{'start':4761.6, 'end':4762.4,'trans':7},{'start':4545.5, 'end':4546.,'trans':13}],
            "CIV":[ {'start':6174.8 , 'end':6175.2, 'trans':1},{'start':6184.7, 'end':6185.1,'trans':0}],
            "SiIV":[{'start':5558.4 , 'end':5559.0, 'trans':1},{'start':5594.4,  'end':5595.,'trans':0}],
            "OVI":[{'start':4115.6 , 'end':4115.9, 'trans':1},{'start':4138.1,  'end':4138.4,'trans':0}],
            "CII":[{'start':5322.3 , 'end':5323.1, 'trans':0}],
            "SiII":[{'start':6088.7 , 'end':6089.2, 'trans':2},{'start':5026.5 , 'end':5027.2, 'trans':4},{'start':4759. , 'end':4759.3, 'trans':5}],
            "OI":[{'start':5193. , 'end':5193.6, 'trans':1}],
            "SiIII_":[{'start':4811.9 , 'end':4812.0, 'trans':0}],
            "CIII_":[{'start':3896.67 , 'end':3896.82, 'trans':0}],




            }



    model=Model(xmlfile="/home/scott/research/J1201+0116/2015-09-23Extraction/metal.xml")
    for key, val in dct.items():
        ab=None
        for item in Model.get(model.AbsorberList):
            if item.id==key:
                ab=item
                break
        
        if ab is None:
            print([item.id for item in Model.get(model.AbsorberList)])
            raise Exception()

        for item in val:

            spec=spec_parser.Spectrum.sniffer("/home/scott/research/J1201+0116/2015-09-23Extraction/combined.dat")
            #_s=spec.get_cut(start=item['start'], end=item['end'])  #should be separate, in case there is abs.
            spec=spec.get_cut(start=item['start'], end=item['end'])

            err=np.median(spec.error)
            cont=np.median(spec.cont)

            SNR=np.median(cont)/np.median(err)

            print("%s %s logNmax=%lf \n"%(key,ab.get_rest_wave(item['trans']),get_maxN(ab,item['trans'],cont,SNR,spec.flux.tolist(), spec.error.tolist() )))
    
       
