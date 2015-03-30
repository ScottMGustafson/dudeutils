#determine if there is a significant difference in dh vals for diff continua choices.  if so, gives you median val for each range.  From this, you can estimate the continuum contribution to systematic error.


import model
import data_types
import scipy.stats as stats
import numpy as np
import dudeutils
import histogram
import math
import matplotlib
import matplotlib.pyplot as plt


dbname = 'database.xml'
c = 299792.458
velrange=(-0.5,0.5)
alpha = 0.05

plt.rcParams['xtick.major.pad']='10'
plt.rcParams['ytick.major.pad']='10'


def get_CI(data, CI=0.95):
    data = sorted(data)
    alpha=(1.-CI)/2.
    length=len(data)
    middle_pts = int(math.floor(0.95*length))
    pass
    

def get_best_fits(db,continuum,num_models=20):
    #constraints = Constraint(**{"D":{"NLocked":False},"H":{"NLocked":False},"chi2":(0.,1810.)})
    data = [item for item in db if item.ContinuumPointList==continuum] 
    #data = [item for item in db if item.get_datum(id,'Absorber','zLocked')==True]
    #data = [item for item in db if item in constraints] 
    data = [item for item in data if velrange[0]<=item.get_shift("D","H")<=velrange[1] and item.get_shift("H2","H")<-58.] 
    data=sorted(data, key=lambda x: x.chi2)
    
    return data[0:num_models]
    #return data

def get_best_cont(continua, db):
    for cont in continua:
        mods = [item for item in db if item.ContinuumPointList==cont[0] and velrange[0]<=item.get_shift('D','H')<=velrange[-1]]
        mods = sorted(mods, key=lambda x: x.chi2)
        _min = mods[0].chi2
        if _min<old_min:
            old_min=_min
            _min = min([float(item.chi2) for item in mods])
    return best


def compare_conts():
    db = model.ModelDB.read(dbname)
    continua = dudeutils.all_conts(db)
    #continua = dudeutils.cont_check_pipeline(reduced_chi2_limit=1.8,verbose=False,db=db)
    #continua = list(set([item.ContinuumPointList for item in db.models]))

    dh = []
    _all=[]
    _continua = []
    for cont in continua:
        data = get_best_fits(db, cont[0])
        #data = [item for item in db if item.ContinuumPointList==cont[0] and velrange[0]<=item.get_shift('D','H')<=velrange[-1] ]
        _continua.append(cont)
        dh.append([item.dh for item in data])
        _all.append(data)
        histogram.histogram(np.array([item.dh for item in data]), title=cont[-1])
        histogram.bootstrap_resampling(np.array([item.dh for item in data]), samples=5000,title=cont[-1]+" bootstrap")
        
    continua = _continua
    #each element of dh is a list of D/H vals for each continuum
        
    W, p = stats.levene(*dh)   #defaults to median as center.
    print(p)  #determines whether variance is significant.
    if p>alpha:
        print("no significant difference in DH values between chosen continua in chosen velocity range")
    for i in range(len(continua)):
        print('%20s med(dh)=%6.3f best dh=%6.3f red chi2=%6.3f datapoints=%4d name=%s'%(continua[i][-1],  np.median(dh[i]), dh[i][0], _all[i][0].reduced_chi2, len(dh[i]), _all[i][0].xmlfile)) 

#    _best_per_cont = [item[0] for item in dh]
#    histogram.histogram(array(_best_per_cont), title='all continua')
#    histogram.bootstrap_resampling(array(_best_per_cont), samples=5000,title='all continua bootstrap distribution (5000 samples)')

    _data = []
    for item in dh:
        _data+=list(item)
    histogram.histogram(np.array(_data), title='all continua')
    histogram.bootstrap_resampling(np.array(_data), samples=5000,title='all continua bootstrap distribution (5000 samples)')

    #print(get_CI(_data))

#    best = get_best_cont(continua, db)
#    dh = [item.dh for item in best]
#    histogram.histogram(array(dh),title=best[0].xmlfile)
#    histogram.bootstrap_resampling(array(dh), samples=5000,title=str(best[0].xmlfile+" bootstrap distribution (5000 samples)"))

def weight_models(lst):
    import numpy as np


    _tot = float(sum([1./item.reduced_chi2 for item in lst]))
    for item in lst:
        item.weight=(1./item.reduced_chi2)/_tot

    #now make a new list with a multiple instances of each element 
    #proportional'ish to their weight.  a bigger weight will increase chance of 
    #being chosen for a given bootstrap
    out = []
    for item in lst:
        mult = int(np.round(1000.*item.weight))
        out+=list( item.dh*np.ones(mult) )
    
    return out

def make_bootstrap(data, stat=np.mean, samples=1500):
    bootstrap_dist = []
    for i in range(samples):
        shuffled = histogram.bootstrap(data)
        bootstrap_dist.append(stat(shuffled))
        #total=np.concatenate( (total, shuffled, axis=0 )
    print('n=',len(data),'SE = ',np.std(bootstrap_dist),'mu = ',np.mean(bootstrap_dist),'med = ',np.median(bootstrap_dist),1.96*np.std(bootstrap_dist)*np.sqrt(len(data)))
    return bootstrap_dist

def plot_bootstrap(dist, samples=5000):
    mean_dist = make_bootstrap(dist, stat=np.mean, samples=samples)
    std_dist  = make_bootstrap(dist, stat=np.std, samples=samples)

    f, (ax0, ax1, ax2) = plt.subplots(1, 3, sharey=True,figsize=[8,3.])
    x0, hist = histogram.histogram(dist, density=True)
    x1, mean = histogram.histogram(mean_dist, density=True)
    x2, std = histogram.histogram(std_dist, density=True)

    ax0.locator_params(nbins=5)
    ax1.locator_params(nbins=5)
    #ax2.locator_params(nbins=5)

    ax0.plot(x0, hist, linestyle='steps')
    ax0.get_xaxis().get_major_formatter().set_useOffset(False)
    ax0.get_yaxis().get_major_formatter().set_useOffset(False)
    ax0.set_ylabel("frequency")
    ax0.set_xlabel("$D/H$")
    ax0.tick_params(axis='y', labelsize=13)
    dhticks=list(np.linspace(-4.63,-4.52,4))
    ax0.tick_params(axis='x', labelsize=13)
    #ax0.set_xticks(dhticks)
    #ax0.set_xticklabels(list(map(str, dhticks)))#, rotation='vertical')


    ax1.plot(x1, mean,linestyle='steps')
    ax1.get_xaxis().get_major_formatter().set_useOffset(False)
    ax1.get_yaxis().get_major_formatter().set_useOffset(False)
    #for tl in ax1.get_yticklabels(): tl.set_visible(False)
    ax1.set_xlabel("$D/H$")
    dhticks=list(np.linspace(-4.59,-4.555,4))
    ax1.set_xlim([-4.59, -4.55])
    ax1.tick_params(axis='x', labelsize=13)
    #ax1.set_xticks(dhticks)
    #ax1.set_xticklabels(list(map(str, dhticks)))#, rotation='vertical')
    ax1.yaxis.set_visible(False)

    #dhticks=list(np.linspace(-4.59,-4.552,5))

    ax2.plot(x2, std, linestyle='steps')
    ax2.tick_params(axis='x', labelsize=13)
    ax2.get_xaxis().get_major_formatter().set_useOffset(False)
    ax2.set_xlabel("$\sigma~$")
    ax2.yaxis.set_visible(False)
    ax2.set_xlim([0.01, 0.035])
    sigticks=list(np.linspace(0.005,0.035,4))
    #sigticks=list(np.linspace(0.0208,0.0250,5))
    ax2.set_xticks(sigticks)
    ax2.set_xticklabels(list(map(str,sigticks)))#, rotation='vertical')

    #plt.setp(plt.xticks()[0], rotation=90)
    #plt.setp(plt.xticks()[1], rotation=90)
    #plt.setp(plt.xticks()[2], rotation=90)

    #ax2.set_xticks([-200,-150,-100,-50,0,50,100,150,200])
    #ax2.set_xticklabels(['-200','','-100','','Zero',    '','100','','200'])

    #plt.subplots_adjust(top=0.87,bottom=0.13,left=0.1,right=0.95)
    #plt.savefig(output_fname)    
    #print('n=',len(data))
    f.subplots_adjust(wspace=0.05)
    plt.tight_layout()

    plt.show()

if __name__ == '__main__':
    db = model.ModelDB.read("all_contsdb.xml")
    continua = dudeutils.all_conts(db)
    continua = list(set([cont[-1] for cont in continua]))
    
    best_models=[]

    for cont in continua:
        models = [item for item in db.models if item.xmlfile==cont]
        models = sorted(models, key=lambda x: float(x.chi2))
        best_models.append(models[0])
        
    #weighted = weight_models(best_models)

    dh = [item.dh for item in best_models]
    ND = [float(item.get_datum("D","Absorber","N")) for item in best_models] 
    NH = [float(item.get_datum("H","Absorber","N")) for item in best_models] 
    print(min(dh), max(dh))
    

    plot_bootstrap(dh, samples=5000)
    #plot_bootstrap(ND, samples=5000)
    #plot_bootstrap(NH, samples=5000)
    #plot_bootstrap(weighted, samples=5000)




    #histogram.histogram(array(dh),xlabel='$log(D/H)$', ylabel='', density=True)
    #histogram.bootstrap_resampling(array(weighted), samples=5000, xlabel='$log(D/H)$ (5000 bootstrap resamples)',ylabel='',density=True)
    #histogram.bootstrap_resampling(array(dh), samples=5000, xlabel='$log(D/H)$ (5000 bootstrap resamples)',ylabel='',density=True)
    
    

