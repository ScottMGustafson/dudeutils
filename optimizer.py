import dude_xmlutils as dude
import numpy as np
from scipy.interpolate import splrep,splprep, splev
import dump_parser as dump
import matplotlib.pyplot as plt
xmlfile='test.xml'
dat_file='test.dat'


dat = dump.FitData(dat_file,xmlfile)

def init_absorber(iden,lockedParam,N_range,b_range,z_range,identity):
    """
    give randomized starting conditions for a given absorber.
    
    input params:
    -------------
    iden:  (str) iden for dude_xmlutils.Absorber class
    lockedParam:  (str) name of param to lock. (str)
    *_range: range of allowed values for given params

    """

    def randomize_param(val_range):
        return (val_range[1] - val_range[0]) * np.random.random_sample() + val_range[0]
    absorber = dude.Absorber(iden=identity,tag='Absorber',xmlfile=xmlfile)
    for param, rnge in [('N',N_range),('b',b_range),('z',z_range)]:
        if lockedParam != param:
            kw[param] = randomize_param(rnge) 
    absorber.write(**kw)
            

def fit_cont(dat):
    #x = np.array([float(item.x) for item in dat.contPoints])
    #y = np.array([float(item.y) for item in dat.contPoints])

    tck = splrep(dat.waves, dat.cont)
    continuum = splev( dat.waves, tck)
    return dat.waves, continuum

def abs_fit(abs_lst):
    """
    get a fit of the continuum with absorption.  voigt profile
    """
    pass

def absorber_fit(N,b,z):
    
    return #some fn to add/subtract from cont
    


#only need to get cont for highlighted regions xb, xe options


xcont = np.array([float(item.x) for item in dat.contPoints])
ycont = np.array([float(item.y) for item in dat.contPoints])
wave, continuum = fit_cont(dat)
fit=continuum + abs_fit(abs_lst)

plt.plot(dat.waves,fit,'g')
plt.plot(dat.waves,dat.cont,'g')
plt.plot(dat.waves,dat.flux,linestyle='steps',color='k')
plt.plot(xcont,ycont,'ro')
plt.xlim([4840.,4890.])
plt.ylim([0,10E-14])
plt.show()
    
    

