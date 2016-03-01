import numpy as np
from astropy.io import fits
import wavelength

class Spectrum(object):
    @staticmethod
    def sniffer(filename, *args, **kwargs):
        """detects raw text format, returns num columns"""
        if filename.endswith('.fits'):
            return FitsSpectrum(filename, error=kwargs.pop('error',None), *args, **kwargs)
        with open(filename) as f:
            if len(f.readline().split())==5:
                return TextSpectrum(filename, *args, **kwargs)
            elif len(f.readline().split())==6:
                return TextSpectrum(filename, *args, **kwargs)
            else:
                raise Exception("unrecognied filetype: %s"%(filename))

    @staticmethod
    def get_indices(lst,xr):
        """gets the indices of interest"""
        lst=np.array(lst)
        ind = []
        start, end = xr[0], xr[-1]
        lst1 = list(np.where(start<=lst)[0])
        lst2 = list(np.where(end>=lst)[0])
        ind = list(set(lst1) & set(lst2))
        return ind

    @staticmethod
    def convert_to_vel(waves,ref_wave):
        return c*(waves-ref_wave)/ref_wave

    @staticmethod        
    def truncate(arr, mod=20):
        ind = list(np.arange(0,arr.shape[0]-1,mod))
        return arr[ind]

    @staticmethod        
    def rebin(arr, binsize=20):
        ind = list(np.arange(0,arr.shape[0]-1,binsize))
        arr=[np.median(arr[i:i+binsize])  if i+binsize<arr.shape[0]-1 else np.median(arr[i:-1]) for i in ind]
        return np.array(arr)

class FitsSpectrum(Spectrum):
    def __init__(self, filename, error=None):
        self.hdu=fits.open(filename)
        if self.hdu[0].header['NAXIS']!=1:
            raise Exception("multidim spec not yet supported")
        self.waves, self.flux = wavelength.xy(filename)
        if error:
            _, self.error=wavelength.xy(error)
        else:
            self.error=None
            
class TextSpectrum(Spectrum): 
    def __init__(self,dumpfile,*lst):
        if len(lst)!=5:
            try:
                lst = np.loadtxt(dumpfile,unpack=True)
            except:
                raise Exception(dumpfile)
            if len(lst)==6:
                lst = np.delete(lst,0,0)  #first column is always a bunch of zeroes
            elif len(lst)==5:
                pass
            else:
                raise Exception("bad formatting for "+dumpfile)
        self.set_data(*lst)
        self.name=dumpfile


    def set_data(self,*lst):  
        assert(len(lst)==5)
        self.waves = lst[0]
        self.flux  = lst[1]
        self.error = lst[2]
        self.abs  = lst[3]
        self.cont  = lst[4]

    @staticmethod
    def get_ind(waves,beg,end):
        if beg>end:  #swap if entered in wrong order
            beg, end = end, beg 
        ind1=np.where(waves<end)[0]
        ind2=np.where(waves>beg)[0]
        return list(set(ind1).intersection(ind2))

    def get_cut(self,start=None,end=None,indices=None):
        if indices is None:
            indices=TextSpectrum.get_ind(self.waves,start,end)
        return TextSpectrum('_.txt', self.waves[indices], self.flux[indices], 
            self.error[indices], self.abs[indices], self.cont[indices])

    def get_ind_at(self,x,vel=False):
        if vel:
            return np.argmin(np.fabs(x-self.vel))
        else:
            return np.argmin(np.fabs(x-self.wave))

    def get_at(self,attr,x):
        assert(attr in ['flux','error','abs','cont'])
        return getattr(self,attr)[self.get_ind_at(x)]

    @staticmethod
    def shift_pixel(x, justification="center"):
        """shift a list by half a pixel"""

        if not justification=="center":
            raise Exception(justification+" type binning not yet supported")

        x=np.array(x).tolist()
        x.append(x[-1]+(x[-1]-x[-2]))  #append a pixel


        xx=[]
        for i in range(1,len(x)):
            try:
                xx.append((x[i]+x[i-1])/2.)
            except:
                print("failure in pixel %d"%(i))
                xx.append(x[i])

        return np.array(xx)

