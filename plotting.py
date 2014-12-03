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
       
def plotDH(db,constraints=None):
    ND, NH, chi2 = xy("D","H","N","N",db,constraints=constraints)
    zH, zD, _ =xy("D","H","z","z",db,constraints=constraints) 

    DH = [float(ND[i])-float(NH[i]) for i in range(len(NH))]
    vel = [(float(zD[i])-float(zH[i]))*c/(1.+float(zH[i])) for i in range(len(zH))]
    chi2 = [float(chi2[i]) for i in range(len(zH))]
    #chi2 = list(map(float, chi2))
    assert(len(DH)==len(vel)==len(chi2))
    fig1, (ax1, ax2) = plt.subplots(nrows=2)
    
    _1 = make_subplot(vel, DH, ax=ax1, yname='D/H')#,xlim=[-1.5,1.5])
    _2 = make_subplot(vel, chi2, ax=ax2, xname='D velocity (km/s)', yname='chi2')#,xlim=[-1.5,1.5])

    #fig2 = plt.figure()
    #plot(x, np.cos(x))
    plt.show()


def plot_cont_vel(db,constraints=None):
    zH2, zH, chi2 = xy("H2","H","z","z",db,constraints=constraints)
    zD, zH, _ = xy("D","H","z","z",db,constraints=constraints)

    vel2 = [(float(zH2[i])-float(zH[i]))*c/(1.+float(zH[i])) for i in range(len(zH))]
    velD = [(float(zD[i])-float(zH[i]))*c/(1.+float(zH[i])) for i in range(len(zH))]
    chi2 = [float(chi2[i]) for i in range(len(zH))]
    fig1, (ax1, ax2) = plt.subplots(nrows=2)
    
    _1 = make_subplot(vel2, velD, ax=ax1, yname='D velocity offset (km/s)')#,xlim=[-1.5,1.5])
    _2 = make_subplot(vel2, chi2, ax=ax2, xname='H2 velocity (km/s)', yname='chi2')#,xlim=[-1.5,1.5])

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
    pass
    #plot 

if __name__ == '__main__':


    constraints = Constraint(**{"D":{"z":(2.988397,2.988426)}, "H2":{"z":(2.98759,2.987637)}})
    #constraints = Constraint(**{"H2":{"z":(2.98759,2.987637)}})

    hidb=dudeutils.load_from_db('2014-11-28hidb.xml')
    #db.trim(constraints)
    plotDH(hidb,constraints=constraints)
    #plotND_NH(hidb,constraints=constraints)
    #plot_cont_vel(hidb,constraints=constraints)

    db=dudeutils.load_from_db('2014-11-28db.xml')
    plotDH(db,constraints=constraints)
    #plotND_NH(db,constraints=constraints)
    #plot_cont_vel(db,constraints=constraints)

    lodb=dudeutils.load_from_db('2014-11-28lodb.xml')
    plotDH(lodb,constraints=constraints)
    #plotND_NH(lodb,constraints=constraints)
    #plot_cont_vel(lodb,constraints=constraints)



