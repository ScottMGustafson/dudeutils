import dudeutils
import matplotlib.pyplot as plt

def xy(x_id, y_id, x,y, dbfile, constraints=None, chi2=False):
    if type(dbfile) is str:
        db = dudeutils.getdb(dbfile)
    else:
        db=dbfile
    _x, chi2 = db.get_all_abs(x_id, str(x), constraints=constraints)
    if y=='chi2' or y_id=="chi2":
        y=chi2
    else:
        _y, _    = db.get_all_abs(y_id, str(y), constraints=constraints)
        assert(chi2==_ and len(_x)==len(_y))
    if chi2:
        return _x, _y, chi2
    else:
        return _x, _y
    

def make_subplot(x,y, style='bo', xname=None,yname=None, ax=None):
    if ax is None:
        ax = plt.gca()
    line, = ax.plot(x, y, style)
    if yname: ax.set_ylabel(yname)
    if xname: ax.set_xlabel(xname)
    
    return line
       

if __name__ == '__main__':

    db = dudeutils.getdb("2014-11-03db.xml")
    ND, NH = xy("D","H","N","N",db)
    velH2, chi2=xy("H2",None,"z","chi2",db)
    

    DH = [float(ND[i])-float(NH[i]) for i in range(len(NH))]
    velH2 = list(map(float, velH2))
    chi2 = list(map(float, chi2))
    assert(len(DH)==len(velH2)==len(chi2))
    fig1, (ax1, ax2) = plt.subplots(nrows=2)
    
    plot(velH2, DH, ax1)
    plot(velH2, chi2, ax2)

    #fig2 = plt.figure()
    #plot(x, np.cos(x))

    plt.show()
