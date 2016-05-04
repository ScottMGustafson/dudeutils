import numpy as np
from astropy.io import fits
import wavelength
from data_types import *
from atomic import *
import sys
from model import Model
import _spectrum

class Spectrum(object):
    @staticmethod
    def sniffer(filename, *args, **kwargs):
        """detects raw text format, returns correct spectrum class instance"""
        if type(filename) in [FitsSpectrum, TextSpectrum]:
            return filename
        elif type(filename) is Model:
            if filename.flux.endswith('.fits'):
                return FitsSpectrum(filename.flux, *args, **kwargs)   
            else: 
                return TextSpectrum(filename, *args, **kwargs)  
        
        elif filename.endswith('.fits'):
            return FitsSpectrum(filename, *args, **kwargs)
        elif filename.endswith('.xml'):
            return Model(xmlfile=filename)
        with open(filename) as f:
            try:
                cols=len(f.readline().split())
            except:
                raise Exception("cannot sniff ",filename)
            if cols==5 or cols==6:
                return TextSpectrum(filename, *args, **kwargs)
            elif cols in [7,8,10,11]:  
                #7, 10 columns for line dump.  cols+1 if ionName has a space in it
                return LineDump(filename, *args, **kwargs)
            
            else:
                raise Exception(
                    "unrecognized filetype: %s \nwith shape %d\nfirst line:\n%s"%(
                    filename, cols,f.readline()))

    @staticmethod
    def get_indices(lst,xr):
        """gets the indices of interest"""
        lst=np.array(lst)
        ind = []
        start, end = xr[0], xr[-1]
        lst1 = list(np.where(start<=lst)[0])
        lst2 = list(np.where(end>=lst)[0])
        ind = list(set(lst1) & set(lst2))
        return sorted(ind)

    @staticmethod
    def convert_to_vel(waves,ref_wave):
        """
        input:
        ------
        waves (float or array of floats): wavelengths
        ref_wave (float or atomic.SpectralLine): a reference absorption line
        """
        if type(ref_wave) is SpectralLine: #then is SpectralLine instance
            ref_wave=ref_wave.obs_wave
        elif type(ref_wave) is not float:
            raise Exception("convert_to_vel only accepts SpectralLine or float types") 
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




    @staticmethod
    def fit_absorption(spec, model, ab_to_plot=None, indices=None):
        """
        Input:
        ------
        spec : specParser.Spectrum instance
        model : model.Model instance

        Output:
        -------
        cont : the continuum of the spectrum
        absorption : absorption of the spectrum
        chi2 : chi-square calculated for spectral regions specified in
               model.RegionList 

        Raises:
        -------
        AssertionError

        """
        def truncate(indices=None):
            try:
                assert(spec.error.shape[0]==spec.waves.shape[0])
            except:
                spec.error=np.zeros(spec.waves.shape[0])
            if indices:
                return spec.waves[indices],spec.flux[indices],spec.error[indices]
            else:
                return spec.waves,spec.flux,spec.error

        wv,flux,e=truncate(indices)

        cont_points=model.get_lst("ContinuumPointList")
        cont_points = sorted(cont_points, key=lambda pt: pt.x)
        x=np.array([float(item.x) for item in cont_points])
        y=np.array([float(item.y) for item in cont_points])

        absorbers=[]

        abslst=ab_to_plot if ab_to_plot else Model.get(model.AbsorberList) 

        if type(abslst[0]) is  SpectralLine:
            absorbers=abslst
        elif type(abslst[0]) is Absorber:
            if not abslst[0]:
                raise Exception("no absorbers specified")
            for item in abslst:
                    absorbers+=item.get_lines()
        else: 
            raise Exception("type of ab_to_plot must either be list of SpectralLine instances or Absorber instances ")

        absorbers=[it for it in absorbers if wv[4]<it.get_obs(it.z)<wv[-4]]


        regions= RegionList.consolidate_list( 
                                  [ (item.start, item.end) 
                                  for item in list(model.get_lst("RegionList"))
                                  if item.start<wv[-1] and item.end>wv[0] ]
                 )

        starts=np.array([item[0] for item in regions],dtype=np.float)
        ends=np.array([item[1] for item in regions],dtype=np.float)

        #convert into numpy array for c extension use
        N=np.array([item.N for item in absorbers], dtype=np.float)
        b=np.array([item.b for item in absorbers], dtype=np.float)
        z=np.array([item.z for item in absorbers], dtype=np.float)
        rest=np.array([item.wave for item in absorbers], dtype=np.float)
        gamma=np.array([item.gamma for item in absorbers], dtype=np.float)
        f=np.array([item.f for item in absorbers], dtype=np.float)



        cont, absorption, chi2= _spectrum.spectrum(  wv,flux,e,
                                                     x,y,
                                                     N,b,z,rest,gamma,f,
                                                     starts,ends )

        return wv,flux,e, absorption, cont, chi2

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
        self.dump=filename.replace('.fits','.dat')
            
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
        self.dump=dumpfile  #this is an alias for self.name

    def set_data(self,*lst):  
        assert(len(lst)==5)
        
        i=0
        for it in "waves flux error abs cont".split():
            setattr(self,it,lst[i])
            i+=1

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
        assert(attr in 'flux error abs cont'.split())
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

class LineDump(object):
    def __init__(self,fname):
        self.fname=fname
        self.absorbers=[]
        for line in open(fname,'r').readlines():
            ionName=line[0:6].replace(' ','')    
            s=line[5:-1].strip().split()
            self.absorbers.append(
                Absorber(ionName=ionName, N=s[0], b=s[1], z=s[2])
            )

    def get_bin(self,rng):
        return [item for item in self.absorbers if rng[0]<=item.z<rng[-1] ]

    def get_ion(self, ionName):
        return [item for item in self.absorbers if item.ionName==ionName]


