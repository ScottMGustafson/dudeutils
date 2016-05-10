import numpy as np
from scipy.stats import chisquare
import astropy.io.fits as fits
import wavelength
import searchArray

Ryd  = 13.60569253  #(eV)
hc   =  0.000123984193  # eV Angstroms
c = 299792.458 #in km/s
pi = np.pi

class FitData(object): 
    def __init__(self,dumpfile,*lst):
        if len(lst)!=5:
            try:
                lst = np.loadtxt(dumpfile,unpack=True)
            except:
                raise Exception(dumpfile)
            lst = np.delete(lst,0,0)  #first column is always a bunch of zeroes
        set_data(*lst)
        name=dumpfile
        assert(len(lst)==5)

    def set_data(self,*lst):  
        waves = lst[0]
        flux  = lst[1]
        error = lst[2]
        abs   = lst[3]
        cont  = lst[4]

    @staticmethod
    def get_ind(waves,start,end):
        lst1 = list(np.where(start<=waves)[0])
        lst2 = list(np.where(end>=waves)[0])
        return list(set(lst1+lst2))

    def get_cut(self,start=None,end=None,indices=None):
        if indices is None:
            indices=FitData.get_ind(start,end)

        return FitData('_.txt', waves[indices], flux[indices], 
            error[indices], cont[indices], abs[indices])

class Spectrum(object):
    def __init__(fits,dump=None):
        self.hdu=fits.open(fits)
        self._wave=wavelength.Wavelength(fits) #store all other wave data, such as wavescales here.
        self.wave=wavelength.get_waves(fits) #array of wavelengths as function of pixel
        self.header=hdu[0].header
        self.flux=hdu[0].data
        if dump:
            self.continuum=np.loadtxt(dump,unpack=True,usecols=(4,))  #for now get from dump.  eventually from ContinuumPointList
            assert(wave.shape[0]==continuum.shape[0])

    def chisquare(self, abs_list, wavelims):
        """
        abs_list:  in instance of data_types.absorberList
        wavelims:  either list of 2 ints, or list of lists of 2 ints
                   for the latter, we are looking at several absorption regions
                   simultaneously
        """

        cont, h, l = self.subtract_all_abs(abs_list)

        if type(wavelims[0]) is int:
            l=searchArray.find_nearest(wavelims[0])
            h=searchArray.find_nearest(wavelims[1])
            chi2 = chisquare(self.flux[l:h], self.continuum[l:h])[0]

        elif type(wavelims[0]) is list and type(wavelims[0][0]) is int:
            chi2=0.
            for item in wavelims:
                l=searchArray.find_nearest(item[0])
                h=searchArray.find_nearest(item[1])
                chi2+=chisquare(self.flux[l:h], self.continuum[l:h])[0]

        return chi2

    def subtract_all_abs(self,absorberlist):
        cont=self.continuum        
        for absorber in abs_list: 
            for ab in absorber.get_lines(): #break into individual absorption lines
                cont=subtract_abs_flux(self.wave,cont,ab)
        return cont

    @staticmethod
    def subtract_abs_flux(wave, cont, ab, subFlag=True):
        """
        subtract the flux from a single absorption line
        
        returns tau, low_pix, hi_pix
        """
        N=ab.N
        b=ab.b
        z=ab.z
        f=ab.f
        gamma=ab.gamma
        restWave=ab.rest
        tau_threshold = 0.001

        cwave = (1.0+z)*restWave
        vdopp = b/cwave*1.0e13
        alpha = gamma/(4.0*np.pi*vdopp)/(1.0+z)    
        fact = np.pow(10.0, N) * 2.647E-2 * \
          f/(np.sqrt(np.pi)*vdopp) * 1.0/(1.0+z)
        
        cpix = searchArray.find_nearest(wave,cwave)
        if (cpix >= wave.shape[0]-3) return cont, cpix, cpix
        if (cpix < 1)              return cont, cpix, cpix

        wave_high = (wave[cpix] + wave[cpix+1])/2
        wave_low  = (wave[cpix] + wave[cpix-1])/2

        nsamp = 10
        delta_wave = (wave_high - wave_low)/nsamp

        tau_sum = 0.0

        for i in range(nsamp):
            wave_i = wave_low + i*delta_wave
            vbar = c/b * (wave_i/cwave - 1.0)
            tau_sum += fact*voigt(vbar, alpha)

        tau = tau_sum / nsamp

        if subFlag: 
            continuum[cpix] -= tau
        else:         
            continuum[cpix] += tau

        lpix = cpix
        while (tau > tau_threshold and lpix > 1):
            lpix -= 1
            wave_high = (wave[lpix] + wave[lpix+1])/2
            wave_low  = (wave[lpix] + wave[lpix-1])/2

            vbar = c/b * (wave_low/cwave - 1.0)
            tau_low = fact*voigt(vbar, alpha)

            vbar = c/b * (wave_high/cwave - 1.0)
            tau_high = fact*voigt(vbar, alpha)

            tau = (tau_low + tau_high)/2

            if subFlag: continuum[lpix] -= tau
            else: continuum[lpix] += tau

        #Now walk the other way ...
        hpix = cpix

        vbar = c/b * (wave[hpix]/cwave - 1.0)
        tau = fact*voigt(vbar, alpha)

        while (tau > tau_threshold and hpix < size()-3):
            hpix += 1

            wave_high = (wave[hpix] + wave[hpix+1])/2
            wave_low  = (wave[hpix] + wave[hpix-1])/2

            vbar = c/b * (wave_low/cwave - 1.0)
            tau_low = fact*voigt(vbar, alpha)

            vbar = c/b * (wave_high/cwave - 1.0)
            tau_high = fact*voigt(vbar, alpha)

            tau = (tau_low + tau_high)/2

            if subFlag: continuum[hpix] -= tau
            else :      continuum[hpix] += tau

        return cont,lpix, hpix


