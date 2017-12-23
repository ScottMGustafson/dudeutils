import gc

import _spectrum
from scipy.signal import convolve

import dudeutils.wavelength as wavelength
from dudeutils.data_types import *
from dudeutils.model import Model


class Spectrum(object):
    @staticmethod
    def sniffer(filename, *args, **kwargs):
        """
        detects raw text format, returns correct spectrum class instance

        Parameters
        ----------
        filename
        args
        kwargs

        Returns
        -------
        Model or Spectrum subclass
        """
        if type(filename) in [FitsSpectrum, TextSpectrum]:
            return filename
        elif type(filename) is Model:
            if filename.flux.endswith('.fits'):
                return FitsSpectrum(filename.flux, error=filename.error)
            else:
                return TextSpectrum(filename.flux, *args, **kwargs)

        elif filename.endswith('.fits'):
            error = kwargs.pop('error', filename.replace('.fits', '_e.fits'))
            return FitsSpectrum(filename, error=error)
        elif filename.endswith('.xml'):
            return Model(xmlfile=filename)
        with open(filename) as f:
            try:
                cols = len(f.readline().split())
            except:
                raise Exception("cannot sniff ", filename)
            if cols == 5 or cols == 6:
                return TextSpectrum(filename, *args, **kwargs)
            elif cols in [7, 8, 10, 11]:
                # 7, 10 columns for line dump.  cols+1 if ionName has a space in it
                return LineDump(filename)

            else:
                raise Exception(
                    "unrecognized filetype: %s \nwith shape %d\nfirst line:\n%s" % (
                        filename, cols, f.readline()))

    @staticmethod
    def get_indices(lst, xr):
        """gets the indices of interest"""
        lst = np.array(lst)
        ind = []
        start, end = xr[0], xr[-1]
        lst1 = list(np.where(start <= lst)[0])
        lst2 = list(np.where(end >= lst)[0])
        ind = list(set(lst1) & set(lst2))
        return sorted(ind)

    @staticmethod
    def convert_to_vel(waves, ref_wave):
        """
        input:
        ------
        waves (float or array of floats): wavelengths
        ref_wave (float or atomic.SpectralLine): a reference absorption line
        """
        if type(ref_wave) is SpectralLine:  # then is SpectralLine instance
            ref_wave = ref_wave.obs_wave
        elif type(ref_wave) is not float:
            raise Exception("convert_to_vel only accepts SpectralLine or float types")
        return c * (waves - ref_wave) / ref_wave

    @staticmethod
    def truncate(arr, mod=20):
        ind = list(np.arange(0, arr.shape[0] - 1, mod))
        return arr[ind]

    @staticmethod
    def rebin(arr, binsize=20):
        ind = list(np.arange(0, arr.shape[0] - 1, binsize))
        arr = [np.median(arr[i:i + binsize]) if i + binsize < arr.shape[0] - 1 else np.median(arr[i:-1]) for i in ind]
        return np.array(arr)

    @staticmethod
    def get_regions(model, xr, wv, xregion):
        if xregion:
            if not xr:
                regions = RegionList.consolidate_list([(item.start, item.end)
                                                       for item in list(model.region_list)
                                                       if item.start < wv[-1] and item.end > wv[0]])
            else:
                regions = RegionList.consolidate_list([[xr[0], xr[-1]]])
        else:
            regions = RegionList.consolidate_list([[wv[0], wv[-1]]])
        return regions

    @staticmethod
    def get_chi2(flux, absorption, e, ind):
        tmp = (flux[ind] - absorption[ind]) / e[ind]
        tmp = tmp[~np.isnan(tmp)]
        tmp = tmp[~np.isinf(tmp)]
        return np.sum(tmp * tmp)

    @staticmethod
    def truncate(spec, indices=None):
        try:
            assert (spec.error.shape[0] == spec.waves.shape[0])
        except:
            spec.error = np.zeros(spec.waves.shape[0])
        if indices is not None:
            return spec.waves[indices], spec.flux[indices], spec.error[indices]
        else:
            return spec.waves, spec.flux, spec.error

    # from memory_profiler import profile
    # @staticmethod
    # @profile
    @staticmethod
    def fit_absorption(spec, model, vsig=3.4, vdisp=2.14,
                       ab_to_fit=None, cont_to_fit=None, get_all=False):
        """
        Input:
        ------
        spec : specParser.Spectrum instance
        model : model.Model instance
        ab_to_fit: included absorbers.  to include no absorbers, include as []

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
        try:
            vdisp, vsig = float(vdisp), float(vsig)
        except:
            print(vdisp, vsig)
            raise Exception()
        wv, flux, e = spec.waves, spec.flux, spec.error

        # wv, flux, e = Spectrum.truncate(spec, indices)
        abslst = ab_to_fit if ab_to_fit else model.absorber_list
        cont_points = cont_to_fit if cont_to_fit else model.cont_point_list
        cont_points = sorted(cont_points, key=lambda pt: pt.x)
        x = np.array([float(item.x) for item in cont_points])
        y = np.array([float(item.y) for item in cont_points])

        spec_lines = []

        if type(abslst[0]) is SpectralLine:
            spec_lines = abslst
        elif type(abslst[0]) is Absorber:
            if not abslst[0]:
                raise Exception("no absorbers specified")
            for item in abslst:
                spec_lines += item.get_lines()
        else:
            raise Exception("type of ab_to_fit must either be list of SpectralLine instances or Absorber instances ")

        spec_lines = [it for it in spec_lines if
                      wv[4] < it.get_obs(it.z) < wv[-4] and it.get_equiv_width(it.N, it.z, pixels=True) > 0.29]

        # convert into numpy array for c extension use
        N = np.array([item.N for item in spec_lines], dtype=np.float)
        b = np.array([item.b for item in spec_lines], dtype=np.float)
        z = np.array([item.z for item in spec_lines], dtype=np.float)
        rest = np.array([item.wave for item in spec_lines], dtype=np.float)
        gamma = np.array([item.gamma for item in spec_lines], dtype=np.float)
        f = np.array([item.f for item in spec_lines], dtype=np.float)

        cont, absorption = _spectrum.spectrum(wv,
                                              x, y,
                                              N, b, z, rest, gamma, f)

        ind = model.get_indices()
        width = int((6.0 * float(vsig) / float(vdisp))) + 1
        # kernel = 1.0 / (vsig * np.sqrt(2.0 * np.pi))
        kernel = np.exp(-0.5 * ((np.arange(width) * vdisp) / vsig) ** 2.)
        kernel /= np.sum(kernel)
        absorption = convolve(absorption, kernel, mode='same')
        absorption = np.concatenate((np.zeros(3), absorption[:-3]))  # looks like the convolve offsets abs by 5 or so

        model.update_dof()
        # ind = model.get_indices()

        chi2 = Spectrum.get_chi2(flux, absorption, e, ind)
        gc.collect()

        return absorption, cont, chi2


class FitsSpectrum(Spectrum):
    def __init__(self, filename, error=None):
        self.waves, self.flux = wavelength.xy(filename)
        if error:
            _, self.error = wavelength.xy(error)
        else:
            raise Exception('error is None.  flux file is %s' % (filename))
        self.dump = filename.replace('.fits', '.dat')


class TextSpectrum(Spectrum):
    attributes = ["waves", "flux", "error", "abs", "cont"]

    def __init__(self, dumpfile, *lst, **kwargs):
        if len(lst) == 0:
            try:
                lst = [kwargs[it] for it in TextSpectrum.attributes]
            except:
                Exception("bad formatting:\n" + str(kwargs))
        if len(lst) != 5:
            try:
                lst = np.loadtxt(dumpfile, unpack=True)
            except:
                raise
                # raise Exception(dumpfile)
            if len(lst) == 6:
                lst = np.delete(lst, 0, 0)  # first column is always a bunch of zeroes
            elif len(lst) == 5:
                pass
            else:
                raise Exception("bad formatting for " + dumpfile)
        self.set_data(*lst)
        self.name = dumpfile
        self.dump = dumpfile  # this is an alias for self.name

    def set_data(self, *lst):
        assert (len(lst) == 5)

        i = 0
        for it in "waves flux error abs cont".split():
            setattr(self, it, lst[i])
            i += 1

    @staticmethod
    def get_ind(waves, beg, end):
        if beg > end:
            beg, end = end, beg
        ind1 = np.where(waves < end)[0]
        ind2 = np.where(waves > beg)[0]
        return list(set(ind1).intersection(ind2))

    def get_cut(self, start=None, end=None, indices=None):
        if indices is None:
            indices = TextSpectrum.get_ind(self.waves, start, end)
        return TextSpectrum('_.txt', self.waves[indices], self.flux[indices],
                            self.error[indices], self.abs[indices], self.cont[indices])

    def get_ind_at(self, x, vel=False):
        if vel:
            return np.argmin(np.fabs(x - self.vel))
        else:
            return np.argmin(np.fabs(x - self.wave))

    def get_at(self, attr, x):
        assert (attr in 'flux error abs cont'.split())
        return getattr(self, attr)[self.get_ind_at(x)]

    @staticmethod
    def shift_pixel(x, justification="center"):
        """shift a list by half a pixel"""

        if not justification == "center":
            raise Exception(justification + " type binning not yet supported")

        x = np.array(x).tolist()
        x.append(x[-1] + (x[-1] - x[-2]))  # append a pixel

        xx = []
        for i in range(1, len(x)):
            try:
                xx.append((x[i] + x[i - 1]) / 2.)
            except:
                print("failure in pixel %d" % (i))
                xx.append(x[i])

        return np.array(xx)


class LineDump(object):
    def __init__(self, fname=None):
        self.absorbers = []
        if fname:
            self.parse(fname)

    def parse(self, fname, excluded_ions=[]):
        """
        parse a line dump and set absorbers attribute

        input: 
        ------
        fname:  string or model.Model instance

        output:
        -------
        None
        
        raises:
        -------
        TypeError when fname is of incorrect type

        """
        if type(fname) is Model:
            self.absorbers = Model.get(fname.AbsorberList)
            self.fname = fname.xmlfile
        elif type(fname) is str:
            self.fname = fname
            if fname.endswith('.xml'):
                mod = Model(xmlfile=fname)
                self.absorbers = [ab for ab in Model.get(mod.AbsorberList)
                                  if not ab in excluded_ions]
            else:
                for line in open(fname, 'r').readlines():
                    ionName = line[0:6].strip()  # .replace(' ','')
                    s = line[5:-1].strip().split()
                    self.absorbers.append(
                        Absorber(ionName=ionName, N=s[0], b=s[1], z=s[2])
                    )
        else:
            raise TypeError("expected str or model.Model instance.  Instead got type %s" % (str(type(fname))))

    def get_bin(self, rng):
        return [item for item in self.absorbers if rng[0] <= item.z < rng[-1]]

    def get_ion(self, ionName):
        return [item for item in self.absorbers if item.ionName == ionName]

    def split(self, ionName):
        from copy import deepcopy
        cpy = deepcopy(self)
        cpy.absorbers = [item for item in self.absorbers if item.ionName == ionName]
        self.absorbers = [item for item in self.absorbers if item.ionName != ionName]
        return cpy


if __name__ == '__main__':
    xmlfile = '/home/scott/J0744/J0744+2059.xml'
    config = '/home/scott/J0744/config.cfg'

    model = Model(xmlfile=xmlfile)
    model.update_dof()
    Config.configure('/home/scott/J0744/config.cfg')
    model = Model(xmlfile=Config.glob['source'])
    model.update_dof()
    spec = Spectrum.sniffer(model)
    vals = Spectrum.fit_absorption(spec, model)
