import numpy as np
from cloudy.cloudy_utils import calc_local_jnu

Ryd  = 13.60569253  #(eV)
hc   =  0.000123984193  # eV Angstroms
pi = 3.141592654

class FitData(object): 
    def __init__(self,dumpfile):
        out = np.loadtxt(dumpfile,unpack=True)
        self.waves = out[1]
        self.flux  = out[2]
        self.error = out[3]
        self.cont  = out[4]
        self.abs   = out[5]

    def get_cut(self,start,end):
        ind = []
        for i in range(0,self.waves.shape[0]):
            if start < self.waves[i] < end:
                ind.append(i)
        return self.waves[ind], self.flux[ind], self.error[ind], self.cont[ind], self.abs[ind]
        

