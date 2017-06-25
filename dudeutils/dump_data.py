import pickle, sys, os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import mlab
from dudeutils.spec_parser import Spectrum, LineDump, TextSpectrum, FitsSpectrum
from scipy.constants import c
from dudeutils.model import Model
import numpy.ma as ma

lya = 1216.668  # lyman alpha wavelength (A)
c /= 1000.  # speed of light in km/s

plt.rc('text', usetex=True)


class DumpData(object):
    def __init__(self, **kwargs):
        for key, val in kwargs.items():
            setattr(self, key, val)

        self.absorbers = [item for item in self.metal.absorbers + self.HI.absorbers]
        self.m1s = [item for item in self.metal.absorbers
                    if item.ionName.lower() == 'm1']  # unidentified metals

        self.metals = [item for item in self.metal.absorbers
                       if item.ionName.lower() not in ["HI", "m1"]]

        self.HIs = [item for item in self.HI.absorbers
                    if item.ionName == "HI"]

        self.name = os.path.split(self.HI.fname)[-1][0:5]
        names = kwargs.get("names")
        self.zem = names[self.name]
        self.HI_sp = kwargs.get("HI_sp", None)
        self.metal_sp = kwargs.get("metal_sp", None)

        self.set_SNR(ab_type='metal')
        self.set_SNR(ab_type='HI')

        self.da, self.da_error = self.get_DA(rest_range=(1090., 1170.))
        self.metal_da = self.get_DA(rest_range=(1090., 1170.), exclude_type="HI")[0]
        self.HI_da = self.get_DA(rest_range=(1090., 1170.), exclude_type="metal")[0]

    @staticmethod
    def format_correction(spec):
        """
        for some reason, in many of the lyaf TextSpectrum outputs of the 
        dataset that I'm looking at, there appears some other data at the end, 
        which I don't understand.  this truncates that.  If your dataset 
        doesn't have the same issue, then this isn't necessary.

        input:
        ------
        spec: dudeutils.spec_parser.Spectrum subclass instance

        output:
        ------
        spec:  truncated to not include data points with lambda<2900 \AA

        """
        indices = np.where(spec.waves > 2900.)[0]
        for attr in "waves flux abs cont error".split():
            setattr(spec, attr, getattr(spec, attr)[indices])
        return spec

    @staticmethod
    def factory(top_dir):
        """
        factory method for DumpData objects.

        input:
        ------
        top_dir: (str) path to top directory containing data.

        output:
        -------
        list of DumpData instances

        Raises:
        -------
        None

        """
        dct = {}
        models = {}
        HI = {}
        metals = {}
        HIabs = {}
        metalabs = {}
        for dirpath, dirname, filename in os.walk(top_dir):
            if str(os.path.split(dirpath)[-1]).strip() in ["ascii (only HI)",
                                                           "ascii (only Metal)",
                                                           "xmls"]:
                for f in filename:
                    fname = os.path.join(dirpath, f)
                    name = f[0:5]

                    if f.endswith('.txt'):
                        if os.path.split(dirpath)[-1] == "ascii (only HI)":
                            if "spec" in os.path.split(fname)[-1]:
                                HI[name] = fname
                            elif "abs" in os.path.split(fname)[-1]:
                                HIabs[name] = LineDump(fname)
                        elif os.path.split(dirpath)[-1] == "ascii (only Metal)":
                            if "spec" in os.path.split(fname)[-1]:
                                metals[name] = fname
                            elif "abs" in os.path.split(fname)[-1]:
                                metalabs[name] = LineDump(fname)
                    elif f.endswith('.xml'):
                        if "HI" in fname and not "only" in fname:
                            models[name] = Model(xmlfile=fname)
                            models[name].read()
                    else:
                        pass

        for key in models.keys():
            try:
                dct[key] = {"HI_sp": Spectrum.sniffer(HI[key]),
                            "metal_sp": Spectrum.sniffer(metals[key]),
                            "HI": HIabs[key],
                            "metal": metalabs[key],
                            "model": models[key]}
            except KeyError:
                print("skipping %s" % key)

            for attr in "HI_sp metal_sp".split():
                if type(dct[key][attr]) is TextSpectrum:
                    dct[key][attr] = DumpData.format_correction(dct[key][attr])

        return [DumpData(**dct[name]) for name in dct.keys()]

    @staticmethod
    def is_red(ab, sp):
        """
        tells whether a given absorber has at least one absorber red of the 
        lyman emission peak.
        
        input:
        ------
        ab: dudeutils.data_types.Absorber instance
        sp: DumpData instance

        output:
        -------
        True if has line red of lya, else False
        """
        lya_em = (1. + sp.zem) * 1216.
        z = ab.z
        for line in ab.get_lines():
            if line.get_obs(z) > lya_em and line.get_obs(z) < sp.HI_sp.waves[-1]:
                return True
        return False

    @staticmethod
    def split_red_blue(spec):
        """
        splits metals into those red of lya and those blue of lya peak

        input:
        ----------
        spec: list of DumpData instances

        output:
        -------
        bluemetals, redmetals: (tuple), lists of dudeutils.data_types.Absorber
                    instances of metals with at least one line blue or red of 
                    the lya emission peak, respectively
        """

        redmetals, bluemetals = [], []
        for item in spec:
            for ab in item.metal.absorbers:
                if ab.ionName != 'm1' and ab.b > 0.317:  # through out really small b values.  .317 chosen b/c 0.316 is a default small value in dude
                    if DumpData.is_red(ab, item):
                        redmetals.append(ab)
                    else:
                        bluemetals.append(ab)
        return bluemetals, redmetals

    def set_SNR(self, vel=200., ab_type='HI', verbose=True):
        """
        sets SNR for each absorber by estimating median of cont/error

        input:
        ------
        vel: velocity limits (+/-) to sample continuum and error
        ab_type: (str) one of 'HI', 'metal' or 'm1'
        verbose: (bool) prints details duriing progress

        output:
        -------
        None.  sets SNR for each absorber
        """

        def get_SNR(line):
            sp = getattr(self, ab_type + "_sp")
            if sp.waves[0] > item.obs_wave or line.obs_wave > sp.waves[-1]:
                return 0.
            ref = line.obs_wave

            cut = sp.get_cut(start=ref - vel * ref / c, end=ref + vel * ref / c)
            # if mostly absorption, i.e. within saturated HI line
            if not ab_type is 'HI':
                if np.median(cut.abs) < 0.6 * np.median(cut.cont):
                    # continuum is not a good indicator of SNR,
                    # bc blended with absorption
                    # take one third of the samples of highest value and use
                    # this to aprox the SNR
                    flux = np.sort(cut.flux)[int(-0.3 * cut.flux.shape[0]):-1]
                    return np.median(flux) / np.median(cut.error)
            return np.median(cut.cont) / np.median(cut.error)

        line_dump = getattr(self, ab_type)

        if verbose: print("getting SNR for %s" % (self.name))
        for ab in line_dump.absorbers:
            lines = []
            try:
                for item in ab.get_lines():
                    lines.append(get_SNR(item))
            except KeyError:
                print("offending sp: %s" % (line_dump.fname))
                print(str(ab))
                raise

            ab.SNR = max(lines)
            if ab.SNR == 0.:
                msg = "\n!!!all zero SNR!!!\n" + str(ab) + "\n  lambda=" + str(
                    [line.obs_wave for line in ab.get_lines()])
                print(msg)
            else:
                if verbose:
                    # print(str(ab)+" SNR=%3.2lf"%(ab.SNR))
                    sys.stdout.write('.')
                    sys.stdout.flush()
        if verbose: print('\n')

    @staticmethod
    def bin_absorbers(dat, bin_attr='z', binsize=None, bin_start=None,
                      bin_end=None, bins=None):
        """
        bins absorbers to specified bins

        input:
        ------
        dat: list of DumpData instances
        bin_attr: attribute over which to bins
        binsize: size of each bin in units of bin_attr
        bin_start, bin_end: start and end values 
        bins: custom bins to specify

        output:
        -------
        output: list of lists, one sublist corresponds to one bin

        """
        output = []
        if not bins:  # if no custom bins specified, then take equallly spaced bins
            bins = np.arange(bin_start, bin_end, binsize)
            bins = list(zip(bins, bins + binsize))

        for b in bins:
            output.append(
                [item for item in dat if b[0] <= getattr(item, bin_attr) < b[1]]
            )
        return output

    def get_DA(self, rest_range=None, exclude_type=None, mask_metals=False):
        """
        Get DA for a spec.  if exlcude, can exclude contribution from what you 
        specify by getting total DA and then adding suspected absorption from 
        whatever class of absorbers.

        input:
        ------
        spec:  Spectrum subclass object instance
        rest_range: rest wavelength region to consider.
        exclude: either `metal` or `HI` or None
                 if None, get simple DA from total spec.

        output:
        -------
        da:  mean(1-f[i]/c[i]) where f is flux at pixel i and c is the continuum 
             level
        """

        # ["waves", "flux","error", "abs","cont"]



        sp = self.metal_sp if exclude_type == "HI" else self.HI_sp
        waves, flux, err, ab, cont = DumpData.normalize(sp, mask_metals=mask_metals)

        indices = None
        if rest_range:
            rest_range = ((1. + self.zem) * rest_range[0], (1. + self.zem) * rest_range[1])
            indices = Spectrum.get_indices(spec.waves, rest_range)

            flux, ab, err = tuple(ma.masked_invalid(it[indices])
                                  for it in (flux, ab, err))
        else:
            flux, ab, err = tuple(ma.masked_invalid(it)
                                  for it in (flux, ab, err))

        if exclude_type:
            da = (1.0 - ab)

        else:
            da = (1.0 - flux)

        dda = err

        da = ma.masked_invalid(da)
        da = ma.masked_where(np.fabs(da) > 10., da)

        dda = ma.masked_invalid(dda)
        dda = ma.masked_where(np.fabs(dda) > 10., da)
        mse = np.sqrt(np.sum(dda ** 2.)) / float(dda.shape[0])
        return np.mean(da), mse

    @staticmethod
    def mask_metals(spec):
        def mask_range(b, z, lam_em, lam_cent):
            """
            this follows the technique from calura et al 2012
            """
            dz = 2. * b * np.sqrt(np.log(2.)) / c
            del_lam = lam_cent * dz
            return lam_cent - del_lam, lam_cent + del_lam

        waves, flux, err, ab, cont = tuple(getattr(spec.HI_sp, attr)
                                           for attr in TextSpectrum.attributes)

        for ab in spec.metals + spec.m1s:
            for line in ab.get_lines():
                beg, end = mask_range(ab.b, ab.z, line.wave, line.get_obs())
                waves = ma.masked_inside(waves, beg, end)

        return tuple(ma.masked_array(getattr(spec.HI_sp, attr), waves.mask)
                     for attr in TextSpectrum.attributes)

    @staticmethod
    def normalize(spec, mask_metals=False, rest_range=None, **kwargs):
        if mask_metals:
            waves, flux, err, ab, cont = DumpData.mask_metals(spec)
        else:
            _type = "metal_sp" if kwargs.get("exclude_type", "metal") == "HI" else "HI_sp"
            waves, flux, err, ab, cont = tuple(ma.masked_invalid(
                getattr(getattr(spec, _type), it))
                                               for it in TextSpectrum.attributes)

        if rest_range:
            waves = ma.masked_outside(waves,
                                      (1. + spec.zem) * rest_range[0],
                                      (1. + spec.zem) * rest_range[-1])

        zrng = kwargs.get('zrng', None)
        snr_rng = kwargs.get('snr_rng', None)

        if zrng:
            xmin, xmax = lya * (1. + zrng[0]), lya * (1. + zrng[1])
            waves = ma.masked_outside(waves, xmin, xmax)
        if snr_rng:
            snr = ma.masked_outside(ma.masked_invalid(cont / err), snr_rng[0], snr_rng[-1])
            waves.mask = np.logical_and(waves.mask, snr.mask)  # combine masks

        flux = ma.masked_array(flux, waves.mask)
        cont = ma.masked_array(cont, waves.mask)
        err = ma.masked_array(err, waves.mask)
        ab = ma.masked_array(ab, waves.mask)

        flux = ma.masked_outside(ma.masked_invalid(flux / cont), -0.2, 1.4)
        err = ma.masked_invalid(err / cont)

        all_masks = waves.mask
        for it in [err, ab, flux]:
            all_masks = np.logical_and(it.mask, all_masks)

        return tuple(ma.masked_array(item, all_masks)
                     for item in (waves, flux, err, ab))
