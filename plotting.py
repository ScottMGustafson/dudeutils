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
    if y_id=="chi2":
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
       
def plotDH(db,constraints=None,the_title=None,dh_lim=None,chi2_lim=None):
    ND, NH, chi2 = xy("D","H","N","N",db,constraints=constraints)
    zH, zD, _ =xy("D","H","z","z",db,constraints=constraints) 

    DH = [float(ND[i])-float(NH[i]) for i in range(len(NH))]
    vel = [(float(zD[i])-float(zH[i]))*c/(1.+float(zH[i])) for i in range(len(zH))]
    chi2 = [float(chi2[i]) for i in range(len(zH))]
    #chi2 = list(map(float, chi2))
    assert(len(DH)==len(vel)==len(chi2))
    fig1, (ax1, ax2) = plt.subplots(nrows=2)
    
    _1 = make_subplot(vel, DH, ax=ax1, yname='D/H',ylim=dh_lim)
    _2 = make_subplot(vel, chi2, ax=ax2, xname='D velocity (km/s)', yname='chi2',ylim=chi2_lim)

    if the_title: plt.title(the_title)
    #fig2 = plt.figure()
    #plot(x, np.cos(x))
    plt.show()


def plot_cont_vel(db,constraints=None,the_title=None):
    zH2, zH, chi2 = xy("H2","H","z","z",db,constraints=constraints)
    zD, zH, _ = xy("D","H","z","z",db,constraints=constraints)

    vel2 = [(float(zH2[i])-float(zH[i]))*c/(1.+float(zH[i])) for i in range(len(zH))]
    velD = [(float(zD[i])-float(zH[i]))*c/(1.+float(zH[i])) for i in range(len(zH))]
    chi2 = [float(chi2[i]) for i in range(len(zH))]
    fig1, (ax1, ax2) = plt.subplots(nrows=2)
    
    _1 = make_subplot(vel2, velD, ax=ax1, yname='D velocity offset (km/s)')#,xlim=[-1.5,1.5])
    _2 = make_subplot(vel2, chi2, ax=ax2, xname='H2 velocity (km/s)', yname='chi2')#,xlim=[-1.5,1.5])
    if the_title: plt.title(the_title)
    #fig2 = plt.figure()
    #plot(x, np.cos(x))
    plt.show()

def plotND_NH(db,constraints=None,the_title=None):
    ND, NH, chi2 = xy("D","H","N","N",db,constraints=constraints)
    zH, zD, _ =xy("D","H","z","z",db,constraints=constraints) 

    ND = [float(ND[i]) for i in range(len(ND))]
    NH = [float(NH[i]) for i in range(len(NH))]
    chi2 = [float(chi2[i]) for i in range(len(zH))]
    fig1, (ax1, ax2) = plt.subplots(nrows=2)
    
    _1 = make_subplot(ND, NH, ax=ax1, yname='NH',xlim=(12.8,12.95))
    _2 = make_subplot(ND, chi2, ax=ax2, xname='ND', yname='chi2',xlim=(12.8,12.95))
    if the_title: plt.title(the_title)
    #fig2 = plt.figure()
    #plot(x, np.cos(x))

    plt.show()

def plotchi2():
    pass
    #plot 

if __name__ == '__main__':


    constraints = Constraint(**{"D":{"z":(2.988397,2.988426)}, "H2":{"z":(2.98759,2.987637), "b":(0.,40.)}})
    #constraints = Constraint(**{"H2":{"z":(2.98759,2.987637)}})

    #hidb=dudeutils.load_from_db('2014-11-28hidb.xml')
    #db.trim(constraints)
    #plotDH(hidb,constraints=constraints,the_title="hi")
    #plotND_NH(hidb,constraints=constraints,the_title="hi")
    #plot_cont_vel(hidb,constraints=constraints,the_title="hi")

    db=dudeutils.load_from_db('2014-12-08db.xml')
    plotDH(db,constraints=constraints,the_title="best",dh_lim=[-4.65,-4.4])
    plotND_NH(db,constraints=constraints,the_title="best")
    plot_cont_vel(db,constraints=constraints,the_title="best")

    #lodb=dudeutils.load_from_db('2014-11-28lodb.xml')
    #plotDH(lodb,constraints=constraints,the_title="lo")
    #plotND_NH(lodb,constraints=constraints,the_title="lo")
    #plot_cont_vel(lodb,constraints=constraints,the_title="lo")



