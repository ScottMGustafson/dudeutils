import multiprocessing
import time
from copy import deepcopy
from random import shuffle

import matplotlib.pyplot as plt

from dudeutils.config import Config
from dudeutils.model import *
from dudeutils.model_csv_io import ModelIO
from dudeutils.simulated_annealing import simulated_annealing
from dudeutils.spec_parser import Spectrum
from dudeutils.utilities import *

c *= 0.001  # convert to km/s
default_iterations = 15


class ModelFailed(BaseException):
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return self.msg


def plt_sigmas(db, ax, x_attr, y_attr, xargs=[], yargs=[], **kwargs):
    min_chi2 = min([item.chi2 for item in db])
    xfact = kwargs.get("xfact", 1.)
    yfact = kwargs.get("yfact", 1.)

    def one_sig(mod): return mod.chi2 <= min_chi2 + 1.

    def two_sig(mod): return 1. < mod.chi2 <= min_chi2 + 4.

    def morethan_2sig(mod): return 4. < mod.chi2

    ax.plot(xfact * np.array(db.get_attr_lst(x_attr, morethan_2sig, *xargs)),
            yfact * np.array(db.get_attr_lst(y_attr, morethan_2sig, *yargs)),
            'ko'
            )

    ax.plot(xfact * np.array(db.get_attr_lst(x_attr, two_sig, *xargs)),
            yfact * np.array(db.get_attr_lst(y_attr, two_sig, *yargs)),
            'bo'
            )

    ax.plot(xfact * np.array(db.get_attr_lst(x_attr, one_sig, *xargs)),
            yfact * np.array(db.get_attr_lst(y_attr, one_sig, *yargs)),
            'co'
            )


def get_nsigma(db, n=1):
    """
    returns database of all models within n*sigma of best fitting model

    input:
    ------
    db: model.ModeDB database of models
    n: number of standard deviations.  n=1 (68% CI) by default

    output:
    -------
    db: model.ModelDB database of models containing only models within n*sigma
        of best fitting model

    Raises:
    -------
    None

    """

    lst = sorted(db.models, key=lambda x: x.chi2)
    return [item for item in lst if item.chi2 < lst[0].chi2 + n * n]


def _get_item(_id, lst):
    for i, it in enumerate(lst):
        if it.id == _id:
            return i
    raise Exception(_id, ' not in list ', [str(it) for it in lst])


class RandomSample(object):
    def __init__(self, sp, mod, writr, kwargs):
        self.spec = sp
        self.model = mod
        self.writr = writr
        self.kwargs = kwargs

    def __call__(self):
        simulated_annealing(self.spec, self.model, self.writr, **self.kwargs)


def random_sampling_multiproc(model, ab_cfg, cont_cfg, **kwargs):
    nproc = max([1, int(kwargs.get('num_processes', multiprocessing.cpu_count() - 2))])

    args = assemble_jobs(model, ab_cfg, cont_cfg, **kwargs)
    print('%d jobs with %d processes.' % (len(args), nproc))
    while True:
        if not multiprocessing.active_children():
            if len(args) == 0:
                break
            for _ in range(min([nproc, len(args)])):
                p = multiprocessing.Process(target=RandomSample(*args.pop()))
                print(p.name)
                p.start()
        else:
            time.sleep(1)

    print('done.  now joining data')
    model_outputs = kwargs.get('model_outputs')
    path, fname = os.path.split(model_outputs)

    for prefix in generate_prefixes(ab_cfg, cont_cfg):
        ModelIO(path=path, output_name=prefix + '.csv', overwrite=True).join_csv(prefix=prefix)


# def random_sampling_pool(model, ab_cfg, cont_cfg, **kwargs):
#     """
#     run random sampling in multiple processes, inspired by https://pymotw.com/3/multiprocessing/communication.html
#     Parameters
#     ----------
#     model
#     ab_cfg
#     cont_cfg
#     kwargs
#
#     Returns
#     -------
#
#     """
#     nproc = max([1, int(kwargs.get('num_processes', multiprocessing.cpu_count() - 2))])
#
#     args = assemble_jobs(model, ab_cfg, cont_cfg)
#     print('%d jobs with %d processes.' % (len(args), nproc))
#     p = multiprocessing.Pool(processes=nproc, maxtasksperchild=1)
#     p.starmap(simulated_annealing, [tuple(arg) for arg in args], chunksize=1)
#     p.close()
#     p.join()
#
#     print('done.  now joining data')
#     src = kwargs.get('source')
#     path, _ = os.path.split(src)
#     ModelIO(path=path).join_csv()
#     model = ModelIO.read_best(model, os.path.join(path, 'output.csv'))
#     return model


# def random_sampling(mod, ab_cfg, cont_cfg, **kwargs):
#     """
#     random sampling in one process
#
#     Parameters
#     ----------
#     spec
#     model
#     ab_cfg
#     cont_cfg
#     kwargs
#
#     Returns
#     -------
#
#     """
#
#     # go through each interesting param one-by-one.  note that unless otherwise specified,
#     # both cont and abs will be varied simultaneously
#     mod.update_dof()
#     sp = Spectrum.sniffer(mod)
#     path, _ = os.path.split(kwargs.get('source'))
#     writr = ModelIO(path, output_name="output.csv")
#     for _ in range(int(kwargs.get('num_iterations', default_iterations))):
#         ab_lst, cnt_lst = get_region_objects(mod)
#         for _id, dct in ab_cfg.items():
#             i = _get_item(_id, ab_lst)
#             for param, rng in dct.items():
#                 print(param, rng)
#                 old_val = getattr(ab_lst[i], param)
#                 _rng = float(rng[-1]) - float(rng[0])
#                 setattr(ab_lst[i], param + 'Locked', True)
#                 print('iteration ', str(int(_) + 1))
#                 setattr(ab_lst[i], param, np.fabs(np.random.normal(loc=old_val, scale=np.fabs(_rng))))
#                 simulated_annealing(sp, mod, ab_lst, cnt_lst, writr, kwargs)
#
#             break
#         break


def assemble_jobs(model, ab_cfg, cont_cfg, **kwargs):
    tasks = []
    # ab_lst, cnt_lst = get_region_objects(model)  # returns a copy of the lists
    sp = Spectrum.sniffer(model)
    path, _ = os.path.split(kwargs.get('source'))
    ctr = 0
    num_it = int(kwargs.get('num_iterations', default_iterations))
    iterations = list(range(num_it))
    shuffle(iterations)

    for _id, dct in ab_cfg.items():
        i = _get_item(_id, model.absorber_list)
        for param, rng in dct.items():
            all_writr_prefix = "output_%s_%s" % (_id.replace(' ', ''), str(param))
            _rng = float(rng[-1]) - float(rng[0])
            setattr(model.absorber_list[i], param + 'Locked', True)

            for k in Config.ab_cfg.keys():  # set the others to unlocked
                if k != _id:
                    _i = _get_item(k, model.absorber_list)
                    for _param in Config.ab_cfg[k].keys():
                        setattr(model.absorber_list[_i], _param + 'Locked', False)

            for _j in iterations:
                mod = model.copy()
                writr = ModelIO(path, output_name=all_writr_prefix + "_%d_%d.csv" % (ctr, _j))
                setattr(mod.absorber_list[i], param, rng[0] + _rng * (1. - float(_j / num_it)))

                msg = ""
                for x in mod.absorber_list:
                    msg += "%6.4lf %s %6.4lf %s %12.11lf %s " % (
                        x.N, "l" if x.NLocked else " ",
                        x.b, "l" if x.bLocked else " ",
                        x.z, "l" if x.zLocked else " ")
                print(msg)
                tasks.append(tuple([sp, mod, writr, kwargs]))  # model needs to be treated as if immutable
                ctr += 1
            setattr(model.absorber_list[i], param + 'Locked', False)

    for _id, dct in cont_cfg.items():
        i = _get_item(_id, model.cont_point_list)
        for param, rng in dct.items():
            param = param.replace('lim', '')  # xlim --> x, ylim --> y
            all_writr_prefix = "output_%s_%s" % (_id.replace(' ', ''), param)
            # all_writr = ModelIO(path, output_name=all_writr_prefix + '.csv')
            old_val = getattr(model.cont_point_list[i], param)
            oldlock = getattr(model.cont_point_list[i], param + 'Locked')
            setattr(model.cont_point_list[i], param + 'Locked', True)
            for _j in iterations:
                mod = model.copy()
                writr = ModelIO(path, output_name=all_writr_prefix + "_%d_%d.csv" % (ctr, _j))
                ctr += 1
                setattr(mod.cont_point_list[i], param, old_val * (1. - 2. * rng * float(_j / num_it - 0.5)))
                tasks.append(tuple([sp, mod, writr, kwargs]))
                setattr(model.cont_point_list[i], param + 'Locked', oldlock)

    return tasks


def get_region_objects(model):
    cont_pts = deepcopy(model.cont_point_list)  # need all so continuum is correct
    ab_lst = [ab for ab in deepcopy(model.absorber_list) if ab.in_region(model.region_list)]
    return ab_lst, cont_pts


def plot_spec(model):
    spec = Spectrum.sniffer(model)
    ab, cont, chi2 = Spectrum.fit_absorption(spec, model, vdisp=Config.vdisp, vsig=Config.vsig, get_all=True)
    waves, flux, error = spec.waves, spec.flux, spec.error
    plt.plot(waves, flux, linestyle='steps', color='k')
    plt.plot(waves, ab, linestyle='steps', color='r')
    plt.plot(waves, cont, linestyle='steps', color='b')
    plt.show()


def run(cfg_file, **kwargs):
    Config.configure(cfg_file)
    model = Model(xmlfile=Config.glob['source'])
    model.update_dof()
    kwargs = dict(kwargs, **Config.glob)
    random_sampling_multiproc(model, Config.ab_cfg, Config.cont_cfg, **kwargs)
    # plot_spec(model)
    # random_sampling(model, ab_cfg, cont_cfg, **glob)


def generate_prefixes(ab_cfg, cont_cfg):
    for _id, dct in dict(ab_cfg, **cont_cfg).items():
        for param, rng in dct.items():
            param = param.replace('lim', '')
            yield "output_%s_%s" % (_id.replace(' ', ''), param)


if __name__ == "__main__":
    Config.configure('/home/scott/J1201/config.cfg')
    print(generate_prefixes(Config.ab_cfg, Config.cont_cfg))
