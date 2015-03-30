import numpy as np

Ryd  = 13.60569253  #(eV)
hc   =  0.000123984193  # eV Angstroms
pi = 3.141592654

class FitData(object): 
    def __init__(self,dumpfile,*lst):
        if len(lst)!=5:
            try:
                lst = np.loadtxt(dumpfile,unpack=True)
            except:
                raise Exception(dumpfile)
            lst = np.delete(lst,0,0)  #first column is always a bunch of zeroes
        self.set_data(*lst)
        self.name=dumpfile
        assert(len(lst)==5)

    def set_data(self,*lst):  
        self.waves = lst[0]
        self.flux  = lst[1]
        self.error = lst[2]
        self.abs   = lst[3]
        self.cont  = lst[4]

    @staticmethod
    def get_ind(waves,start,end):
        lst1 = list(np.where(start<=waves)[0])
        lst2 = list(np.where(end>=waves)[0])
        return list(set(lst1+lst2))

    def get_cut(self,start=None,end=None,indices=None):
        if indices is None:
            indices=FitData.get_ind(start,end)

        return FitData('_.txt', self.waves[indices], self.flux[indices], 
            self.error[indices], self.cont[indices], self.abs[indices])

