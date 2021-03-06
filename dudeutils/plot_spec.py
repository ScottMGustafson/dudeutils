# import dudeutils
from dudeutils.spec_parser import *
from dudeutils.atomic import *
import matplotlib.pyplot as plt
from matplotlib import rc
from dudeutils.model import *
import os
import numpy as np
import numpy.ma as ma
import lineid_plot

lineid_fontsize = 16
rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size': 16})
rc('text', usetex=True)


def get_lines(model, dct):
    return [model.get_spectral_line(key, val) for key, val in dct.items()]


def plot_line(spec, model, fig, line, **kwargs):
    velocity = kwargs.pop("velocity", True)
    ax = kwargs.pop("ax", None)
    multiplier = kwargs.pop('multiplier', None)
    xlabel = kwargs.pop('xlabel', None)
    topmost = kwargs.pop('topmost', False)
    ref = kwargs.pop('ref', None)
    ylabel = kwargs.pop('ylabel', None)
    debug = kwargs.pop('debug', False)
    if debug: print("prepping data")
    wave_r, ind, ref = prep_data(spec, model, ref, **kwargs)
    absorbers = [item for item in line] if type(line) is list else [line]
    if debug:
        print("getting ab")
        print(spec, model.flux)
    waves, flux, error, all_ab, cont = get_ab(spec, model, ref, ind, absorbers,
                                              return_ab=False, all_ab=True)

    ab = [all_ab]
    if len(absorbers) > 0:
        if debug:
            print(spec, model.flux)
        ab += list(get_ab(spec, model, ref, ind, absorbers, return_ab=True, all_ab=False))

    if multiplier:
        def f(x): return multiplier * x

        cont, flux, error = tuple(multiplier * x for x in (cont, flux, error))
        ab = [f(it) for it in ab]
        # ab=[multiplier*np.array(it) for it in ab]

    tol = kwargs.pop('tol', 1000.)

    cont = ma.masked_outside(cont, -1. * tol, tol)
    flux = ma.masked_outside(flux, -1. * tol, tol)
    error = ma.masked_outside(error, -1. * tol, tol)

    assert (flux.shape[0] == ab[0].shape[0] == cont.shape[0] == waves.shape[0])
    if debug: print("plotting")
    ax.plot(waves, flux, 'k', linestyle='steps-mid')
    ax.plot(waves, np.zeros(waves.shape[0]), 'k--')
    ax.plot(waves, cont, 'r--', linestyle='steps-mid')
    try:
        ax.plot(waves, error, 'g', linestyle='steps-mid')
    except:
        pass

    def plot_ab(ab):
        """
        plot absorption separate from continuum so we can plot the indidual 
        absorption from diff absorbers
        """
        if type(ab) is list:
            assert (type(ab[0]) is np.ndarray)
            for item in ab:
                plot_ab(item)
        else:
            if not np.array_equal(ab, cont):
                ab = ma.masked_outside(ab, -1. * tol, tol)
                ax.plot(waves, ab, 'b--')

    plot_ab(ab)

    x, labels, y = [], [], []

    for item in prep_absorbers(spec, model, absorbers):
        for line in item:
            obs_wave = Spectrum.convert_to_vel(line.obs_wave, ref)
            if np.amin(waves) < obs_wave < np.amax(waves):
                x.append(obs_wave)
                if hasattr(line, "custom_label"):
                    labels.append(line.custom_label)
                else:
                    labels.append(line.ionName)
                y.append(np.mean(cont))

    box_loc = lineid_plot.get_box_loc(fig, ax, x, y, box_axes_space=0.06)
    if topmost:
        for i in range(len(x)):
            ax.annotate(labels[i], xy=(x[i], y[i]),
                        xytext=(box_loc[i][0],
                                0.7 * box_loc[i][1]),
                        xycoords="data", textcoords="data",
                        rotation=90, horizontalalignment="center",
                        verticalalignment="center",
                        fontsize=lineid_fontsize,
                        arrowprops=dict(arrowstyle="-",
                                        relpos=(0.5, 0.0)),
                        label=labels[i])
            ax.axvline(x=x[i], ymin=0, ymax=0.05 * y[i], linewidth=1,
                       color='k', linestyle='--')

    else:
        for i in range(len(x)):
            ax.axvline(x=x[i], ymin=0, ymax=y[i], linewidth=1, color='k', linestyle='--')
    if xlabel:
        ax.set_xlabel(xlabel)
    else:
        ax.axes.get_xaxis().set_visible(False)
        ax.set_xlabel(None)
    # ax.set_xlim(kwargs.get("xlims"),None)
    ax.set_ylim(kwargs.get("ylims"), None)

    if ylabel:
        ax.yaxis.set_label_position("right")
        ax.set_ylabel(ylabel)

    return ax


def parse_spectrum(model):
    error = None
    if type(model) is str:
        if model.endswith('.xml'):
            model = Model(xmlfile=model)
    else:
        assert (type(model) is Model)

    spec = Spectrum.sniffer(model)

    #    spec.cont, spec.abs, chi2 = Spectrum.fit_absorption(spec,model)
    return model, spec


def prep_absorbers(spec, model, absorbers):
    """
    converts Absorber to list of SpectralLine instances
    """
    if type(absorbers) is list:
        return [prep_absorbers(spec, model, ab) for ab in absorbers]
    if not absorbers:
        return []

    if type(absorbers) is dict:
        absorbers = get_lines(model, absorbers)

    elif type(absorbers) is SpectralLine:
        pass
    else:
        raise TypeError("type of absorbers must either be list of " + \
                        "SpectralLine instances or dict")
    return absorbers


def get_wr_from_vr(ref, model, vel_r):
    return [(1 + vel_r[0] / c) * ref, (1 + vel_r[1] / c) * ref]


def prep_data(spec, model, ref=None, wave_r=None, vel_r=None, **kwargs):
    if ref:
        if type(ref) is dict:
            assert (len(ref.keys()) == 1)  # should only be one element long
            key = list(ref.keys())[0]
            ref = model.get_spectral_line(key, ref[key]).obs_wave
        elif type(ref) is SpectralLine:
            ref = model.get_spectral_line(key, ref[key]).obs_wave
        elif type(ref) is float:
            pass
        else:
            raise Exception("poorly defined input: ref")

    if vel_r:
        if wave_r or not ref:
            raise Exception("poorly defined input")
        else:
            wave_r = get_wr_from_vr(ref, model, vel_r)

    if wave_r:
        ind = Spectrum.get_indices(spec.waves, wave_r)
    else:
        raise Exception("test fail")  # np.arange(0,spec.waves.shape[0])

    if len(ind) == 0:
        raise Exception("No data selected")

    return wave_r, ind, ref


def get_ab(spec, model, ref, ind, absorbers, vel_r=True, return_ab=True, all_ab=False):
    if all_ab:
        absorbers = None
    else:
        if type(absorbers) is list:
            return [get_ab(spec, model, ref, ind, ab, vel_r, return_ab) for ab in absorbers]
        absorbers = prep_absorbers(spec, model, absorbers)
    if type(spec) is TextSpectrum:
        waves, flux, error, ab, cont = spec.waves[ind], spec.flux[ind], spec.error[ind], spec.abs[ind], spec.cont[ind]
        if not all_ab:
            waves, flux, error, ab, cont, chi2 = Spectrum.fit_absorption(spec,
                                                                         model,
                                                                         ab_to_plot=absorbers,
                                                                         indices=ind)

    else:
        waves, flux, error, ab, cont, chi2 = Spectrum.fit_absorption(spec,
                                                                     model,
                                                                     ab_to_plot=absorbers,
                                                                     indices=ind)

    if vel_r:
        waves = Spectrum.convert_to_vel(waves, ref)

    if return_ab:
        return ab
    else:
        return waves, flux, error, ab, cont


def label_lines(ax, x, y, label):
    for i in range(len(label)):
        ax.annotate(label[i], xy=(x[i], y[i]), xytext=(x[i], 0.8), ha='center', va='center',
                    arrowprops=dict(arrowstyle="-", facecolor='black'))
        # plt.axvline(x=x[i], ymin=0., ymax=0.68, linewidth=1, color='k')
    return ax


def plot_no_lines(center_wave, wave_range, exponent=14):
    model, spec = parse_spectrum(Model(xmlfile="/home/scott/research/J0744+2059/DH_1comp.copy.xml"))
    waves, flux, error, ab, cont, chi2 = Spectrum.fit_absorption(spec,
                                                                 model)
    flux *= 10 ** float(exponent)
    x, y = [], []
    for i in range(0, waves.shape[0]):
        if center_wave - np.fabs(wave_range[0]) < waves[i] < center_wave + np.fabs(wave_range[-1]):
            x.append(waves[i])
            y.append(flux[i])
    plt.plot(x, y, linestyle='steps', color='k')
    plt.xlabel(r'\r{A}ngstroms')
    plt.ylabel(r'$F_{\lambda}$ ($10^{%d}$ ergs s$^{-1}$ cm$^{-2}$ \r{A}$^{-1}$)' % (int(exponent)))
    plt.show()


def plot_continua(spec, db, fact=13, xlims=[4842, 4860.], wavelength=4847.25):
    indices = Spectrum.get_indices(spec.waves, xlims)

    count, ind = 0, 0
    for mod in db:
        waves, flux, error, ab, cont, chi2 = Spectrum.fit_absorption(spec, mod, indices=indices)
        if count == 0:
            plt.plot(waves, flux * 10. ** fact, color='k', linestyle='steps')
            ind = np.argmin(np.fabs(waves - wavelength))
        assert (ind != 0)
        mod.cont_level = cont[ind]
        mod.ab_level = ab[ind]
        if count == 0:
            plt.axvline(x=wavelength, ymin=0., ymax=mod.cont_level,
                        linewidth=1, color='k', linestyle='--')

        plt.plot(waves, ab * 10. ** fact, 'r-', alpha=0.1)
        plt.plot(waves, cont * 10. ** fact, 'b-', alpha=0.1)
        count += 1
    plt.axvline(x=wavelength, ymin=0., ymax=mod.cont_level, linewidth=1, color='k', linestyle='--')
    plt.xlabel(r"Wavelength (\r{A})")
    plt.ylabel(r"$F_{\lambda}$($10^{%d}$ erg s$^{-1}$ cm$^{-2}$ \r{A}$^{-1}$)" % (int(fact)))
    plt.xlim(xlims)
    plt.ylim([-0.05, 10.])
    plt.minorticks_on()
    plt.show()


if __name__ == "__main__":
    # this is an example for plotting a few transitions of SiII

    xmlfile, vel_r = "/home/scott/research/J0744+2059/SiII.xml", [-120., 120.]
    model, spec = parse_spectrum(Model(xmlfile=xmlfile))
    fig, axes = plt.subplots(3, sharex=True)

    axes[0] = plot_line(spec, model,
                        fig,
                        [{"SiII2": 2}, {"SiII1": 2}, {"SiII3": 2}],
                        vel_r=vel_r, velocity=True, ax=axes[0],
                        ref={"SiII2": 2},
                        multiplier=10E14,
                        topmost=True)

    axes[1] = plot_line(spec, model,
                        fig,
                        [{"SiII2": 3}, {"SiII1": 3}, {"SiII3": 3}],
                        vel_r=vel_r, velocity=True, ax=axes[1],
                        ref={"SiII2": 3},
                        multiplier=10E14)

    axes[2] = plot_line(spec, model,
                        fig,
                        [{"SiII1": 5}, {"SiII2": 5}, {"SiII3": 5}],
                        vel_r=vel_r, velocity=True, ax=axes[2],
                        ref={"SiII2": 5},
                        multiplier=10E14)

    fig.text(0.04, 0.5, 'Flux', va='center', rotation='vertical', size=22)

    plt.show()
