import corner
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
from scipy import stats
from dudeutils.data_types import Absorber


def fit_power_law(x, y, ye, return_log=False):
    """
    if x is redshift, z, remember to enter arg as 1+z
    """
    x = np.log10(x)
    y, ye = lin2log(y, ye)
    m, b, r, p, stderr, db, dm = linear_regression(x, y, ye=ye)
    if return_log:
        def logfn(xx):
            return m * xx + b, np.sqrt((dm * xx) ** 2. + db ** 2.)

        return m, dm, b, db, logfn
    else:

        b, db = log2lin(b, db)

        def linfn(xx):
            return b * xx ** m, np.sqrt(xx ** (2. * m) * (db ** 2. + (dm * np.log(xx)) ** 2.))

        return m, dm, b, db, linfn


def linear_regression(x, y, verbose=True, clip=True, xe=None, ye=None):
    m, b, r, p, stderr = stats.linregress(x, y)
    assert (len(list(x)) == len(list(y)))
    if not ye is None:
        def f(x, m, b):
            return m * x + b

        assert (len(list(x)) == len(list(ye)))

        popt, pcov = curve_fit(f, x, y, sigma=ye, absolute_sigma=True)
        m = popt[0]
        b = popt[1]
        dm = np.sqrt(pcov[0][0])
        db = np.sqrt(pcov[1][1])
        print(pcov)
        return m, b, r, p, stderr, db, dm
    else:
        mx = np.mean(x)
        sx2 = np.sum(((x - mx) ** 2))
        db = np.sqrt(len(x)) * stderr * np.sqrt(1. / len(x) + mx * mx / sx2)
        dm = np.sqrt(len(x)) * stderr * np.sqrt(1. / sx2)
    return m, b, r, p, stderr, db, dm


def mad(data, axis=None):
    return np.mean(np.absolute(data - np.mean(data, axis)), axis)


def opt_bin_width(lst):
    """
    optimum bin width: freedman-diaconis rule
    h=2*IQR*n**0.3333

    input:
    ------
    lst:  list of floats
    output: 
    -------
    bin width, num bins
    """
    lst = np.sort(lst)
    q75, q25 = np.percentile(lst, [75, 25])
    iqr = q75 - q25
    h = 2. * iqr * np.array(lst).shape[0] ** -0.33333
    num_bins = (lst[-1] - lst[0]) / h
    return h, num_bins


def get_hist(values, nbins=False, attr='b'):
    """
    make a histogram.

    input:
    ------
    values: list of items.  
    nbins: number of bins.
    attr:  (optional.  default='b')  dudeutils.data_types.Absorber attribute

    output:
    -------
    x: bin centers
    hist: histogram value
    """
    if type(values[0]) is Absorber:
        values = np.array([getattr(item, attr) for item in values])

    if not nbins:
        _, nbins = opt_bin_width(values)
    hist, bin_edges = np.histogram(values, int(nbins))
    hist = hist / len(list(values))  # normalize hist by shape
    binsize = bin_edges[1] - bin_edges[0]
    x = bin_edges[:-1] + 0.5 * binsize
    return x, hist


def hui_rutledge(b, bsig):
    """
    returns number density per b: dN/db

    see hui and rutledge, 1999, APJ 517:541-548

    bsig is a parameter describing the data set's line widths.

    bsig^4 := 2/<\delta''^2> where \delta'' is a second derivative
    for the deviation of optical depth from the mean optical depth.
    in other words, it characterizes the curvature of a line at the line center.
    
    """
    b4 = (bsig / b) ** 4.
    return 4. * (b4 / b) * np.exp(-1. * b4)


def log2lin(logb, dlogb):
    b = 10. ** np.array(logb)
    db = b * np.log(10.) * np.array(dlogb)
    return b, db


def lin2log(b, db):
    return np.log10(b), np.array(db) / (np.array(b) * np.log(10.))


def get_sigma(cov, x, mu, alpha=0.05):
    """
    get an elipse centered in mu with axis over the eigenvectors of
    cov matrix and effective radius -2ln(alpha):
    
    (x-mu).T * np.linalg.inv(cov) * (x-mu) = -2*np.log(alpha)

    """
    pass


def get_cov(db, params):
    """get covariance matrix with weights from chi2 
    db:  ModelDB instance
    params: list of dicts, {absorber_id:param}
    data: data formatted such that list or np array of lists or np arrays.  
    each sub array is a set of measurements for one parameter.
    """

    # likelihood=np.array([np.exp(0.5*item.chi2) for item in db])
    # aweights=likelihood/np.sum(likelihood)


    data = []
    param_names = []
    for item in params:
        assert (len(item.keys()) == 1)
        for key, val in item.items():
            data.append([it.get_datum(key, 'Absorber', val) for it in db])
    return np.cov(data)


def spearman(db, id1, id2, attr1, attr2):
    x = np.array([item.get_datum(id1, "Absorber", attr1) for item in db.models])
    y = np.array([item.get_datum(id2, "Absorber", attr2) for item in db.models])
    return stats.spearmanr(x, y)


def cov(db, id1, id2, attr1, attr2, best1, best2):
    x = np.array([item.get_datum(id1, "Absorber", attr1) for item in db.models])
    y = np.array([item.get_datum(id2, "Absorber", attr2) for item in db.models])

    return np.sum((x - best1) * (y - best2) / x.shape[0])


def corr(db, id1, id2, attr1, attr2, best1, best2, err1, err2):
    return cov(db, id1, id2, attr1, attr2, best1, best2) / (err1 * err2)


def plot_xy(db, id1, id2, attr1, attr2):
    min_chi2 = min([item.chi2 for item in db.models])

    models_68, models_95, models = [], [], []
    for item in db.models:
        if item.chi2 <= min_chi2 + 1.0:
            models_68.append(item)
        elif item.chi2 <= min_chi2 + 4.:
            models_95.append(item)
        else:
            models.append(item)

    x_68 = np.array([item.get_datum(id1, "Absorber", attr1) for item in models_68])
    y_68 = np.array([item.get_datum(id2, "Absorber", attr2) for item in models_68])

    x_95 = np.array([item.get_datum(id1, "Absorber", attr1) for item in models_95])
    y_95 = np.array([item.get_datum(id2, "Absorber", attr2) for item in models_95])

    x = np.array([item.get_datum(id1, "Absorber", attr1) for item in models])
    y = np.array([item.get_datum(id2, "Absorber", attr2) for item in models])

    plt.plot(x, y, 'ko')
    plt.plot(x_95, y_95, 'bo')
    plt.plot(x_68, y_68, 'ro')

    plt.xlabel(attr1 + "(" + id1 + ")")
    plt.ylabel(attr2 + "(" + id2 + ")")

    plt.show()


def plt_corners(db, params, labels, title="", fname="corner.png"):
    """
    data should be arranges such that it is of shape (nsamples, ndim)

    """
    # Plot it.

    # likelihood=np.array([item.chi2 for item in db])
    # weights=likelihood/np.sum(likelihood)

    best = sorted(db.models, key=lambda x: x.chi2)[0]

    truths = []
    data = []
    for item in params:
        assert (len(item.keys()) == 1)
        for key, val in item.items():
            data.append([it.get_datum(key, 'Absorber', val) for it in db])
            truths.append(best.get_datum(key, 'Absorber', val))

    data = np.array(data).T

    nsamples, ndim = data.shape
    # Plot it.
    figure = corner.corner(data, labels=labels,
                           truths=truths,
                           quantiles=[0.16, 0.5, 0.84], \
                           show_titles=True, title_kwargs={"fontsize": 12},
                           verbose=True)
    figure.gca().annotate(title, xy=(0.5, 1.0), xycoords="figure fraction",
                          xytext=(0, -5), textcoords="offset points",
                          ha="center", va="top")
    figure.savefig(fname)


def grubbs(samp, alpha=0.05, attr=None, **kwargs):
    """
    from: 
        Grubbs, Frank E. (1950). "Sample criteria for testing outlying
        observations". Annals of Mathematical Statistics 21 (1): 27–58.
        doi:10.1214/aoms/1177729885.

    """

    def G(sample, val):
        return np.fabs((val - np.mean(sample))) / np.std(sample)

    if attr:
        samp = sorted(list(samp), key=lambda x: getattr(x, attr))
        _samp = [getattr(it, attr) for it in samp]
        Gmax = G(_samp, max(_samp))
        Gmin = G(_samp, min(_samp))
    else:
        samp = sorted(list(samp))
        Gmax = G(samp, max(samp))
        Gmin = G(samp, min(samp))

    interval = stats.t.interval(alpha / float(2 * len(samp)), len(samp) - 2)

    if Gmax > interval[-1]:
        print(samp[-1])
        del (samp[-1])
    if Gmax < interval[0]:
        print(samp[0])
        del (samp[0])

    return samp


def weights(lst):
    _tot = float(sum([1. / item ** 2. for item in lst]))
    return [(1 / item ** 2.) / _tot for item in lst]


def bootstrap(data):
    """
    returns a resampling of the dataset with length of data
    """
    data = np.array(data)
    indices = np.random.randint(len(data), size=len(data))
    return data[indices]


def bootstrap_resampling(data, stat=np.mean, samples=1000, attr=None, **kwargs):
    bootstrap_dist = []
    if attr:
        data = [getattr(it, attr) for it in list(data)]
    for i in range(samples):
        shuffled = bootstrap(data)
        bootstrap_dist.append(stat(shuffled))
    return bootstrap_dist
