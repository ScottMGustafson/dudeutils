import dudeutils
from model import c
import matplotlib.pyplot as plt

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
    
    return line
       
def plotDH(db):
    ND, NH, chi2 = xy("D","H","N","N",db)
    zH, zD, _ =xy("D","H","z","z",db) 

    DH = [float(ND[i])-float(NH[i]) for i in range(len(NH))]
    vel = [(float(zD[i])-float(zH[i]))*c/(1.+float(zH[i])) for i in range(len(zH))]
    chi2 = [float(chi2[i]) for i in range(len(zH))]
    #chi2 = list(map(float, chi2))
    assert(len(DH)==len(vel)==len(chi2))
    fig1, (ax1, ax2) = plt.subplots(nrows=2)
    
    _1 = make_subplot(vel, DH, ax=ax1, yname='D/H',xlim=[-1.2,1.2])
    _2 = make_subplot(vel, chi2, ax=ax2, xname='D velocity (km/s)', yname='chi2',xlim=[-1.2,1.2])

    #fig2 = plt.figure()
    #plot(x, np.cos(x))

    plt.show()

if __name__ == '__main__':

    db=dudeutils.load_from_db('2014-11-24db.xml')
    plotDH(db)
