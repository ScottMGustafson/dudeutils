import numpy as np
import dudeutils.histogram as histogram
import matplotlib.pyplot as plt
from dudeutils.utilities import load_from_db
from dudeutils.model import ModelDB

#constraints_1sig = {'D':{'b':(13.0,14.6), 'N':(12.826,12.852), 'shift':(-0.3,0.3)}, 
 #              'H2':{'b':(8.9,11.), 'shift':(-61.4,-60.3), 'N':(12.62,12.68)}, 
 #              'xmlfile':"2015-01-22_best.xml"}

#constraints = {'D':{'b':(12.,16.), 'N':(12.8,12.88), 'shift':(-0.5,0.5)}, 
#               'H2':{'b':(0.,12.), 'shift':(-62.,-59.5), 'N':(12.55,12.75)}, 
#              'xmlfile':"2015-01-22_best.xml",'chi2':1790.}

constraints={}

def plot_distribution(db, dist_name="$D/H$",constraints=constraints):
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size=18)
    if type(db) is str:
        db = dudeutils.load_from_db(db)
    data = ModelDB.filter(db.models, constraints)
    f, ax1 = plt.subplots(1, figsize=[6,5])

    x1, mean = histogram.histogram(dist, density=True)

    ax1.plot(x1, mean,linestyle='steps')
    ax1.get_xaxis().get_major_formatter().set_useOffset(False)
    ax1.get_yaxis().get_major_formatter().set_useOffset(False)
    #for tl in ax1.get_yticklabels(): tl.set_visible(False)
    ax1.set_xlabel("$D/H$")
    ax1.set_ylabel("frequency")

    print(np.std(dist))
    plt.tight_layout()
    plt.show()

def plot_chi2(db, attr, xlabel="$N$", iden=None,  iden2=None, constraints=constraints):    
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size=18)
    if type(db) is str:
        db = dudeutils.load_from_db(db)
    data = ModelDB.filter(db.models, constraints)
    f, ax1 = plt.subplots(1, figsize=[6,5])

    if attr in ["N","b","z"]:
        x = [item.get_datum(iden=iden,tag="Absorber",param=attr) for item in data]
    elif attr in ["vel", "velocity"]:
        x = [item.get_vel(iden,iden2) for item in data]
    else:
        x = [getattr(item,attr) for item in data]
    y = [item.chi2 for item in data]

    ax1.plot(x,y,'ro')
    ax1.get_xaxis().get_major_formatter().set_useOffset(False)
    ax1.get_yaxis().get_major_formatter().set_useOffset(False)
    #for tl in ax1.get_yticklabels(): tl.set_visible(False)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel("$\chi^2$")
    ax1.minorticks_on()
    #plt.minorticks_on()
    plt.tight_layout()
    plt.show()

def chi2_xy(db, attr, xlabel="$N$", filename="out.txt",iden=None, iden2=None, constraints=constraints):    
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size=18)
    if type(db) is str:
        db = dudeutils.load_from_db(db)
    data = ModelDB.filter(db.models, constraints)
    f, ax1 = plt.subplots(1, figsize=[6,5])

    if attr in ["N","b","z"]:
        x = [item.get_datum(iden=iden,tag="Absorber",param=attr) for item in data]
    elif attr in ["vel", "velocity"]:
        x = [item.get_vel_shift(iden,iden2) for item in data]
    else:
        x = [getattr(item,attr) for item in data]
    y = [item.chi2 for item in data]

    ax1.plot(x,y,'ro')
    ax1.get_xaxis().get_major_formatter().set_useOffset(False)
    ax1.get_yaxis().get_major_formatter().set_useOffset(False)
    #for tl in ax1.get_yticklabels(): tl.set_visible(False)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel("$\chi^2$")

    plt.minorticks_on()
    plt.tight_layout()
    plt.show()
    with open(filename,"w") as f:
        for i in range(len(x)):
            f.write("%15E %15E \n"%(x[i], y[i]))
    return

    
  
  

