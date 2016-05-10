import astropy.io.fits as fits
import numpy as np
from operator import attrgetter
import sys
import math
c = 299792.458     # speed of light in km/s

class WaveUtils1D(object):
    def __init__(self, crval, crpix, cdelt, loglin):
        #specify types as a check for type
        self.crval=float(crval)
        self.crpix=int(crpix)
        self.crdelt=float(cdelt)
        self.loglin=bool(loglin)

    def iterate(i,f):
        curr = self.get_wave(i)
        if i<f:
            while i<f:
                if self.loglin:
                    curr*=10.**self.cdelt
                else:
                    curr+=self.cdelt
                yield curr
                i+=1
        elif i>f:
            while i>f:
                if self.loglin:
                    curr/=10.**self.cdelt
                else:
                    curr-=self.cdelt
                yield curr
                i-=1
        else:
            raise StopIteration

    def get_wave(self, i):
        return WaveUtils1D.wave(self.crval, self.crpix, 
                                self.cdelt, i, self.loglin)

    @staticmethod
    def wave(crval1,crpix1,cdelt1,i,loglin=True, angstroms=True):
        """get wavelength of individual pixel"""
        crval1,crpix1,cdelt1=tuple(map(float, (crval1,crpix1,cdelt1)))
        if type(i) is np.ndarray:
            i=i.astype(np.float32)
        elif type(i) is int:
            i=float(i)
        else:
            pass
            
        if angstroms:
            crval1*=1.E3

        if not loglin:
            return crval1+cdelt1*(i-crpix1)
        else:
            return crval1*10.**(cdelt1*(i-crpix1))
    @staticmethod
    def is_loglin(head):
        if not head["CTYPE1"]=="LINEAR":
            raise Exception("only linear and loglinear wavescales accepted.")

        return True if head['DC-FLAG']==1 else False
        

class Wavelength(object):
    """wavelength utilities for 1D or 2D spectum"""
    def __init__(self,filename=None, shift=None):
        self.filename=filename
        self.hdu=fits.open(filename)
        self.head=self.hdu[0].header
        self.data=self.hdu[0].data
        self.size=self.hdu[0].data.shape[0]
        try:
            self.shift = float(shift)
        except:
            self.shift=0.0
        self.get_1D_params()
        self.loglin = WaveUtils1D.is_loglin(self.head)

    def shift_wave(self,make_permanent=False,shift=None,output=None):
        """shift=velocity shift in km/s
           make permanent:  write out to original file"""
        if not output:
            output=self.filename
        if shift:
            try:
                shift=float(shift)
            except:
                shift=float(self.shift)
            if self.hdu[0].header['CRVAL1']<1000.:
                shift/=1000.  #sometimes CRVAL not in angstroms
            if self.loglin:  
                if make_permanent:
                    self.hdu[0].header['CRVAL1']*=(1.+shift/c)
                    self.hdu.writeto(output,output_verify='fix',clobber=True)
                return WaveUtils1D.wave(float(self.crval1)*(1.+shift/c), self.crpix1, 
                                            self.cdelt1, np.arange(0,self.size), 
                                            self.loglin)
            else:
                if make_permanent:
                    self.hdu[0].header['CRVAL1']+=shift
                    self.hdu.writeto(output,output_verify='fix',clobber=True)
                return WaveUtils1D.wave(float(self.crval1)+shift, self.crpix1, 
                                            self.cdelt1, np.arange(0,self.size), 
                                            self.loglin)
        else:
            return WaveUtils1D.wave(float(self.crval1), self.crpix1, 
                                        self.cdelt1, np.arange(0,self.size), 
                                        self.loglin)
        
    def get_1D_params(self):
        for key in ['CTYPE1', 'CDELT1', 'CRPIX1', 'CRVAL1']:
            setattr(self,key.lower(),self.head[key])

    def xy(self,shift=None, make_permanent=False):
        x=self.shift_wave(make_permanent=make_permanent,shift=shift)
        return x, self.data

    def shiftCoeffs(self):
        try:
            self.shift = float(self.shift) 
        except TypeError:    #should not happen
            raise TypeError("TypeError: WaveData.__init__().    "+type(self.shift))
        except:
            msg = "unknown error: WaveData.__init__()\n"
            msg+= "shift type: %s \n" % str(type(self.shift))
            msg+= "shift = %lf\n" % str(self.shift)
            msg+= "coeff shape: %d    by    %d\n" % (self.coeffs.shape[0], self.coeffs.shape[1])
            msg+= "coeffs = %lf ...\n"%(self.coeffs[0])
            raise Exception(msg)
#is this the right thng to do?
        if self.loglin:
            return self.coeffs[0]*(1.+self.shift/c) #shift the wavescale by const. 
        else:
            return self.coeffs[0]+self.shift


    def closeInstance(self):
        self.hdu.close()

#aliases
def get_waves(filename):  
    return Wavelength.get_wave(filename)

def xy(filename, shift=None, make_permanent=False):
    inst=Wavelength(filename)
    x, y = inst.xy(shift, make_permanent=make_permanent)
    inst.hdu.close()
    return x, y

def shift(filename, shift, output=None):
    if not output:
        output=filename
    Wavelength(filename,shift).shift_wave(make_permanent=True,output=output)

