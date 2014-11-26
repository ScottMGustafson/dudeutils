import dudeutils
from model import c
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from constraints import Constraint
from matplotlib import rc

def xy(x_id, y_id, x,y, dbfile, constraints=None, chi2=False):
    if type(dbfile) is str:
        db = dudeutils.getdb(dbfile)
    else:
        db=dbfile
    _x, chi2 = db.get_all_abs(x_id, str(x), constraints=constraints)
    if y=='chi2' or y_id=="chi2":
        _y=chi2
    else:
        _y, _    = db.get_all_abs(y_id, str(y), constraints=constraints)
        assert(chi2==_ and len(_x)==len(_y))
    return _x, _y, chi2
    

def make_subplot(x,y, style='bo', xname=None,yname=None, ax=None,xlim=None,ylim=None):
    if ax is None:
        ax = plt.gca()
    line, = ax.plot(x, y, style)
    
    if yname: ax.set_ylabel(yname)
    if xname: ax.set_xlabel(xname)
    if not xlim is None: ax.set_xlim(xlim)
    if not ylim is None: ax.set_ylim(ylim)
    ax.ticklabel_format(useOffset=False)
    return line
       
def plotDH(db,constraints=None):
    ND, NH, chi2 = xy("D","H","N","N",db,constraints=constraints)
    zH, zD, _ =xy("D","H","z","z",db,constraints=constraints) 

    DH = [float(ND[i])-float(NH[i]) for i in range(len(NH))]
    vel = [(float(zD[i])-float(zH[i]))*c/(1.+float(zH[i])) for i in range(len(zH))]
    chi2 = [float(chi2[i]) for i in range(len(zH))]
    #chi2 = list(map(float, chi2))
    assert(len(DH)==len(vel)==len(chi2))
    fig1, (ax1, ax2) = plt.subplots(nrows=2)
    
    _1 = make_subplot(vel, DH, ax=ax1, yname='D/H',xlim=[-1.5,1.5])
    _2 = make_subplot(vel, chi2, ax=ax2, xname='D velocity (km/s)', yname='chi2',xlim=[-1.5,1.5])

    #fig2 = plt.figure()
    #plot(x, np.cos(x))
    plt.show()

def plotND_NH(db,constraints=None):
    ND, NH, chi2 = xy("D","H","N","N",db,constraints=constraints)
    zH, zD, _ =xy("D","H","z","z",db,constraints=constraints) 

    ND = [float(ND[i]) for i in range(len(ND))]
    NH = [float(NH[i]) for i in range(len(NH))]
    chi2 = [float(chi2[i]) for i in range(len(zH))]
    fig1, (ax1, ax2) = plt.subplots(nrows=2)
    
    _1 = make_subplot(ND, NH, ax=ax1, yname='NH',xlim=(12.7,13.1))
    _2 = make_subplot(ND, chi2, ax=ax2, xname='ND', yname='chi2',xlim=(12.7,13.1))

    #fig2 = plt.figure()
    #plot(x, np.cos(x))

    plt.show()

def plotchi2():
    hi = model.ModelDB(name="chi2_hi_database.txt.new", params=11, datapts=957, chi2min=1980., chi2lim=2000.)
    hi_onesigdb, hidb = get_onesig(hi,hi.chi2min)   
    data["hi"] = get_data(hidb)
    data["hionesig"] = get_data(hi_onesigdb)

    best = model.ModelDB(name="test_database.txt.new", params=8, datapts=957, chi2lim=1825., chi2min=1801., best=True)
    before_clip = get_data(best)
    best_onesigdb, bestdb = get_onesig(best,best.chi2min) 
    data["best"] = get_data(bestdb)
    data["bestonesig"] = get_data(best_onesigdb) 

    lo = model.ModelDB(name="chi2_lo_database.txt.new", params=11, datapts=957, chi2min=1980., chi2lim=2000.)
    lo_onesigdb, lodb = get_onesig(lo,lo.chi2min)  
    data["lo"] = get_data(lodb)
    data["loonesig"] = get_data(lo_onesigdb) 



#now set up the plots

    f = plt.figure(figsize=(4.1,9.3))

#hi
    gs1 = GridSpec(2,1)
    gs1.update(left=0.2, right=0.9, top=0.95, bottom=0.68, hspace=0.1)
    ax1 = plt.subplot(gs1[:-1, 0])
    ax2 = plt.subplot(gs1[-1, 0])
    #ax1.set_ylim([1980.,2000.])
    #ax2.set_ylim([-5.0,-4.4])
    #plt.setp( ax2.get_xticklabels(), visible=False)

    plot_one_pair(ax1,ax2,data["hi"], data["hionesig"])

#best
    gs2 = GridSpec(2,1)
    gs2.update(left=0.2, right=0.9, top=0.64, bottom=0.36, hspace=0.1)
    ax3 = plt.subplot(gs2[:-1, 0])
    ax4 = plt.subplot(gs2[-1, 0])
    #ax3.set_yticks(np.arange(1800.,1825.,5.))
    #ax3.set_ylim([1795.,1825.])
    #ax4.set_ylim([-5.0,-4.4])
    #plt.setp( ax4.get_xticklabels(), visible=False)
    plot_one_pair(ax3,ax4,data["best"], data["bestonesig"])
#lo
    gs3 = GridSpec(2,1)
    gs3.update(left=0.2, right=0.9, top=0.32, bottom=0.05, hspace=0.1)
    ax5 = plt.subplot(gs3[:-1, 0])
    ax6 = plt.subplot(gs3[-1, 0])
    #ax5.set_ylim([1980.,2000.])
    #ax6.set_ylim([-5.0,-4.4])
    plot_one_pair(ax5,ax6,data["lo"], data["loonesig"]) 

    #ax6.set_xlabel(r"velocity from {\textrm{H}\,\textsc{i}~} (km s$^{-1}$)")

    plt.show()
    #plt.savefig('/home/scott/Desktop/hi_best_lo.png',dpi=600)
    plt.clf()


if __name__ == '__main__':

    constraints = Constraint(**{"H2":{"z":(2.98755,2.9876)},"chi2":(2130.,2142.)})

    db=dudeutils.load_from_db('2014-11-24db.xml')
    db.trim(constraints)
    plotDH(db)#,constraints={"H2":{"z":(2.9875,2.98765)}})
    plotND_NH(db)
