from dudeutils import *
import numpy as np
import matplotlib.pyplot as plt
import warnings

def best_fit_DH(moddb,order,param="N",xmin=None,xmax=None, plot=True, constraints=None,thecolor='ro'):
    """
    get a best fit of data with respect to `param'

    id: id of absorber
    param:  parameter name (N,b,z)
    order:  order of polynomial to fit
    xmax, xmin: range of values to consider
    locked:  get only locked parameters?
    plot:   plot the data?  otherwise return function, x, y
    """
    x = []
    y = []

    if constraints!=None:
        lst=ModelDB.constrain(moddb,constraints)
        if len(lst)==0:
            raise Exception("no surviving models:\n%s"%(str(constraints)))
        if len(lst)==len(moddb.lst):
            warnings.warn("everything passed")
    else:
        lst = moddb.lst

    for item in lst:

        D = item.get("D","Absorber")
        H = item.get("H","Absorber")
        if D.locked(param) or H.locked(param):
            x.append(float(getattr(D,param))-float(getattr(H,param)))
            y.append(float(item.chi2))

    title=lst[-1].xmlfile

    if len(x)==0 or len(y)==0 or len(x)!=len(y):
        raise Exception("ill condittioned input: \n  x=%s\n  y=%s"%(str(x),str(y)))

    if xmin==None and xmax==None:
        xmax = max(x)
        xmin = min(x)

    x=np.array(x)
    y=np.array(y)

    #coeffs=np.polyfit(x-x.mean(),y,int(order))
    #f = np.poly1d(coeffs)
    if plot:
        #xx = np.arange(xmin,xmax, np.abs(xmax-xmin)/100.)
        #yy = f(xx-x.mean())
        plt.xlim(xmin,xmax)
        #plt.plot(xx,yy,'b-')
        plt.plot(x,y,thecolor)
        plt.title(title)
        plt.show()
    #return f, x, y

if __name__=="__main__":
    db = load_from_db("2014-10-15db.xml")
    z=2.98841195
    vel=2.0
    delz = (vel/299792.458)*(1.+z)
    best_constr = {"D":{"z":(z-delz,z+delz)},"H2":{"z":(2.9874,2.9877),"b":(5.0,30.)},"chi2":1820,"xmlfile":"best_2014-10-15.xml"}
    hi_constr = {"D":{"z":(z-delz,z+delz)},"H2":{"z":(2.9874,2.9877),"b":(5.0,30.)},"chi2":1900,"xmlfile":"hi_2014-10-15.xml"}
    lo_constr = {"D":{"z":(z-delz,z+delz)},"H2":{"z":(2.9874,2.9877),"b":(5.0,30.)},"chi2":1900,"xmlfile":"lo_2014-10-15.xml"}
    #db.best_fit("D","z",6,constraints=constraints)
    best_fit_DH(db,6,constraints=best_constr)
    #best_fit_DH(db,6,constraints=hi_constr,thecolor='bo')
    #best_fit_DH(db,6,constraints=lo_constr,thecolor='go')

