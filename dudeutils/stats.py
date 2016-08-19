import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from dudeutils.utilities import *
import corner
from dudeutils.random_sampling import filter_bad_models, parse_config

#HIA  b<7(68%), N=20.583 (.62,.56) (68%)  (.52--.68 95%)
#DIa  b=8.7(8.5-9.1 68%), N=14.92 (14.68--15.05,  14.56--15.15 (95%))
#DIb  b=12.0 \pm0.7,  N=15.99\pm0.03 (68%),  \pm0.07 (95%)
#HIB  b=7 +10 -4, N=19.9 -0.4 + 0.1 (68%)  19.3--20.2 (95%)
#DIC  b=8.0 -2.2+0.9, N=15.02 -0.24+0.17
#HIC: b=18,6+0.3-0.2,  N=20.15\pm0.08  19.86--20.34 (95%)

def get_sigma(cov, x, mu,alpha=0.05):
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
    
    #likelihood=np.array([np.exp(0.5*item.chi2) for item in db])
    #aweights=likelihood/np.sum(likelihood)
    

    data=[]
    param_names=[]
    for item in params:
        assert(len(item.keys())==1)
        for key, val in item.items(): 
            data.append([it.get_datum(key,'Absorber',val) for it in db])     
    return np.cov(data)


def spearman(db,id1, id2, attr1, attr2):
    x=np.array([item.get_datum(id1,"Absorber",attr1) for item in db.models])
    y=np.array([item.get_datum(id2,"Absorber",attr2) for item in db.models])
    return stats.spearmanr(x,y)

def cov(db,id1, id2, attr1, attr2, best1, best2):
    x=np.array([item.get_datum(id1,"Absorber",attr1) for item in db.models])
    y=np.array([item.get_datum(id2,"Absorber",attr2) for item in db.models])
    
    return np.sum( (x-best1)*(y-best2)/x.shape[0] )

def corr(db,id1, id2, attr1, attr2, best1, best2, err1, err2):
    return cov(db,id1, id2, attr1, attr2, best1, best2) /(err1*err2)

def plot_xy(db,id1, id2, attr1, attr2):  
    min_chi2=min([item.chi2 for item in db.models])

    models_68, models_95, models=[],[],[]
    for item in db.models:
        if item.chi2<=min_chi2+1.0:
            models_68.append(item)
        elif item.chi2<=min_chi2+4.:
            models_95.append(item)
        else:
            models.append(item)

    x_68=np.array([item.get_datum(id1,"Absorber",attr1) for item in models_68])
    y_68=np.array([item.get_datum(id2,"Absorber",attr2) for item in models_68])

    x_95=np.array([item.get_datum(id1,"Absorber",attr1) for item in models_95])
    y_95=np.array([item.get_datum(id2,"Absorber",attr2) for item in models_95])

    x   =np.array([item.get_datum(id1,"Absorber",attr1) for item in models])
    y   =np.array([item.get_datum(id2,"Absorber",attr2) for item in models])

    plt.plot(x,y,'ko')
    plt.plot(x_95,y_95,'bo')
    plt.plot(x_68,y_68,'ro')

    plt.xlabel(attr1+"("+id1+")")
    plt.ylabel(attr2+"("+id2+")")

    plt.show()

def plt_corners(db,params,labels,title="",fname="corner.png"):
    """
    data should be arranges such that it is of shape (nsamples, ndim)

    """
# Plot it.

    #likelihood=np.array([item.chi2 for item in db])
    #weights=likelihood/np.sum(likelihood)

    best=sorted(db.models, key= lambda x:x.chi2)[0]

    truths=[]
    data=[]
    for item in params:
        assert(len(item.keys())==1)
        for key, val in item.items(): 
            data.append([it.get_datum(key,'Absorber',val) for it in db]) 
            truths.append(best.get_datum(key,'Absorber',val))

    data=np.array(data).T

    nsamples, ndim=data.shape
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

def mad(data, axis=None):
    return np.mean(np.absolute(data - np.mean(data, axis)), axis)

def grubbs(samp,alpha=0.05, attr=None,**kwargs):
    """
    from: 
        Grubbs, Frank E. (1950). "Sample criteria for testing outlying
        observations". Annals of Mathematical Statistics 21 (1): 27–58.
        doi:10.1214/aoms/1177729885.

    """
    def G(sample,val):
        return np.fabs((val-np.mean(sample)))/np.std(sample)
    if attr:
        samp=sorted(list(samp), key=lambda x: getattr(x,attr))
        _samp=[getattr(it,attr) for it in samp]
        Gmax=G(_samp,max(_samp))
        Gmin=G(_samp,min(_samp))
    else:
        samp=sorted(list(samp))
        Gmax=G(samp,max(samp))
        Gmin=G(samp,min(samp))
        
    interval=stats.t.interval(alpha/float(2*len(samp)),len(samp)-2)

    
    if Gmax>interval[-1]:
        print(samp[-1])
        del(samp[-1])
    if Gmax<interval[0]:
        print(samp[0])
        del(samp[0])

    return samp

def bootstrap(data):
    """
    returns a resampling of the dataset with length of data
    """
    data=np.array(data)
    indices = np.random.randint(len(data),size=len(data))
    return data[indices]

def bootstrap_resampling(data,stat=np.mean,samples=1000,attr=None,**kwargs):
    bootstrap_dist = []
    if attr:
        data=[getattr(it,attr) for it in list(data)]
    for i in range(samples):
        shuffled = bootstrap(data)
        bootstrap_dist.append(stat(shuffled))
    return bootstrap_dist

