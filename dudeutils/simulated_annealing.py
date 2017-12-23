import copy

import numpy as np
from scipy.stats import chi2 as sp_chi2

from dudeutils.config import Config
from dudeutils.spec_parser import Spectrum


class OptConst(object):
    # default optimzation constraints
    N_epsilon, b_epsilon, z_epsilon = 0.00001, 0.0001, 10e-8
    x_epsilon, y_epsilon = 0.015, 10e-18  # abt half a pixel at 5000 Angstrom
    Nrng, brng, zrng = 0.7, 2.0, 10e-5
    alpha = 0.7
    max_steps = 50
    chi2_pad = 16.
    starting_step = 1.
    Tmin = 0.1
    min_step = 0.03

    @staticmethod
    def set_optimization_constants(kwargs):
        if not kwargs:
            return
        for key, val in kwargs.items():
            try:
                setattr(OptConst, key, float(val))
            except (TypeError, ValueError) as e:
                pass


def attr_vary(ab, attr, rng, step):
    if bool(getattr(ab, attr + 'Locked')):
        return
    _rng = step * float(rng[1] - rng[0]) / 2.
    current = getattr(ab, attr)
    setattr(ab, attr, np.fabs(np.random.normal(loc=current, scale=np.fabs(_rng))))
    if ab.id in dict(Config.ab_cfg, **Config.cont_cfg).keys():  # to prevent drift away from desired models
        if getattr(ab, attr) < rng[0]:
            setattr(ab, attr, current + np.fabs(np.random.normal(loc=0., scale=np.fabs(_rng))))
        elif getattr(ab, attr) > rng[-1]:
            setattr(ab, attr, current - np.fabs(np.random.normal(loc=0., scale=np.fabs(_rng))))


def ab_in_regions(ab, regions):
    for line in ab.get_lines():
        if regions.in_regions(line.get_obs(ab.z)):
            return True
    return False


def accept_prob(old_chi2, new_chi2, T):
    if new_chi2 < old_chi2:
        return 1.0
    else:
        delta = new_chi2 - old_chi2
        return min(np.exp(-delta / T), 1.0)


def set_ab_rng(ab, attr):
    try:
        rng = tuple(map(float, Config.ab_cfg[ab.id][attr]))
    except KeyError:
        rng = (getattr(ab, attr) - getattr(OptConst, attr + 'rng'),
               getattr(ab, attr) + getattr(OptConst, attr + 'rng'))
    return rng


def set_cnt_rng(cnt, attr):
    if cnt.id in Config.cont_cfg.keys():
        _rng = float(Config.cont_cfg[cnt.id]['ylim']) / 2.
        rng = ((1. - _rng) * getattr(cnt, attr), (1. + _rng) * getattr(cnt, attr))
    else:  # offset percentage
        rng = (1. - getattr(OptConst, attr + 'rng') * getattr(cnt, attr),
               1. + getattr(OptConst, attr + 'rng') * getattr(cnt, attr))
    return rng


def write_model(model):
    return [x for x in model.absorber_list if x.id not in ['null', '']], \
           [x for x in model.cont_point_list if x.id != 'null'], \
           model.chi2


class Anneal(object):
    def __init__(self, spec, model, writr=None, **kwargs):
        OptConst.set_optimization_constants(kwargs)
        if kwargs:
            for k, v in kwargs.items():
                setattr(self, k, v)
        self.spec, self.model = spec, model
        self.chi2crit = sp_chi2.ppf(0.99, int(model.dof))
        self.ab_lst, self.cnt_lst = self.model.absorber_list, self.model.cont_point_list

        self.writr = writr
        self.regions = self.model.region_list
        self.step = float(OptConst.starting_step)
        self.annealing_temp = sp_chi2.ppf(0.99, int(self.model.dof)) + OptConst.chi2_pad

        self.model.chi2 = self.get_chi2()
        self.best_model = self.model.copy()
        self.ctr = 0
        self.tie = kwargs.get('tie', None)

    def reset_model(self):
        self.model = None
        self.model = self.best_model.copy()

    def reset_best(self):
        self.best_model = None
        self.best_model = self.model.copy()

    def guess_ab(self, ab):
        if not ab_in_regions(ab, self.regions):
            return ab
        Nrng = set_ab_rng(ab, 'N')
        brng = set_ab_rng(ab, 'b')
        zrng = set_ab_rng(ab, 'z')

        attr_vary(ab, 'N', Nrng, self.step)  # start param range abt a factor of 2 (in log)
        attr_vary(ab, 'b', brng, self.step)
        attr_vary(ab, 'z', zrng, self.step)

        return ab

    def guess_cont(self, cont_point):
        if not bool(getattr(cont_point, 'yLocked')):
            if self.regions.in_regions(cont_point.x):
                y_rng = set_cnt_rng(cont_point, 'y')
                attr_vary(cont_point, 'y', y_rng, self.step)
        return cont_point

    def randomize_parameters(self):
        """
        
        Parameters
        ----------
        tie : (str, str), (optional)
            tuple of '_' separated ids to tie together.
            e.g.:  ('H_z', 'D_z') will tie together an absorber named 'D''s redshift to that of 'H' 
            param which is tied must be numeric type

        Returns
        -------

        """
        self.model.set_absorbers([self.guess_ab(copy.copy(ab)) for ab in self.model.absorber_list])
        self.model.set_cont_points([self.guess_cont(copy.copy(cnt)) for cnt in self.model.cont_point_list])

        if self.tie:
            ab1_id, ab2_id = self.tie[0].split('_')[0], self.tie[1].split('_')[0]
            param = self.tie[0].split('_')[1]
            i1 = self.model.get_ab_index(ab1_id)
            i2 = self.model.get_ab_index(ab2_id)
            if getattr(self.model.absorber_list[i1], param + 'Locked'):
                set_to = float(getattr(self.model.absorber_list[i1], param))
                setattr(self.model.absorber_list[i2], param, set_to)
            else:
                set_to = float(getattr(self.model.absorber_list[i2], param))
                setattr(self.model.absorber_list[i1], param, set_to)

    def get_chi2(self):
        return Spectrum.fit_absorption(self.spec, self.model,
                                       vdisp=Config.vdisp, vsig=Config.vsig,
                                       ab_to_fit=self.model.absorber_list, cont_to_fit=self.model.cont_point_list)[-1]

    def is_bad_model(self, chi2):
        return chi2 > max([self.best_model.chi2 + OptConst.chi2_pad, self.chi2crit])


def simulated_annealing(spec, model, writr=None, **kwargs):
    """
    
    Parameters
    ----------
    spec : spec_parser.Spectrum instance
        spectrum to fit
    model : model.Model instance
        the model to optimize
    writr : model_csv_io.ModelIO (optional)
        writer to which we may write
    kwargs : keyword args for Anneal()

    Returns
    -------
    model.Model instance
        The optimized model
    """
    annealr = Anneal(spec, model, writr, **kwargs)

    while annealr.annealing_temp > OptConst.Tmin and annealr.ctr < OptConst.max_steps:
        i = 0
        while i < OptConst.max_steps:
            old_chi2 = annealr.model.chi2
            annealr.randomize_parameters()
            annealr.model.chi2 = annealr.get_chi2()
            if accept_prob(old_chi2, annealr.model.chi2, annealr.annealing_temp) > np.random.uniform():
                if annealr.model.chi2 < annealr.best_model.chi2:
                    i -= OptConst.max_steps  # you get more iterations!!
                    annealr.reset_best()
                if annealr.step < OptConst.min_step:
                    annealr.annealing_temp *= (1. - float(annealr.ctr / OptConst.max_steps))
                    break
                else:
                    annealr.step *= OptConst.alpha
            else:
                annealr.reset_model()
                i += 1
                continue
        annealr.ctr += 1

    if annealr.writr:  # todo this is just temporary to see how it behaves.  eventually, just return the model
        annealr.writr.write([write_model(annealr.best_model)])
    else:
        return annealr.best_model
