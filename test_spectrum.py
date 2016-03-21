import dudeutils
import _spectrum
#import matplotlib.pyplot as plt
import spec_parser
import sys
import unittest
import numpy as np
import optimizer
from data_types import *
from model import *

class TestSpectrum(unittest.TestCase):
    def setUp(self):
        self.model=dudeutils.get_model(
                        "/home/scott/research/J0744+2059/J0744+2059.xml")

        self.cont_points=self.model.get_lst("ContinuumPointList")
        self.cont_points = sorted(self.cont_points, key=lambda pt: pt.x)
        self.x=np.array([float(item.x) for item in self.cont_points])
        self.y=np.array([float(item.y) for item in self.cont_points])
        self.spec = spec_parser.Spectrum.sniffer(
                                    self.model.flux, error=self.model.error)

        absorbers=[]
        for item in self.model.get_lst("AbsorberList"):
            absorbers+=item.get_lines()

        absorbers=list(filter(
                       lambda x: 
                            self.spec.waves[4]<x.get_obs(x.z)<self.spec.waves[-4], 
                            absorbers))
        regions=Model.get(self.model.RegionList)
        self.starts=np.array([item.start for item in regions],dtype=np.float)
        self.ends=np.array([item.end for item in regions],dtype=np.float)
        
        #convert into numpy array for c extension use

        for item in 'N b z wave gamma f'.split():
            setattr(self,item,
                        np.array(
                                [getattr(i,item) for i in absorbers], 
                                dtype=np.float)
                    )


    def test_find_absorption(self):
        cont, ab, chi2= _spectrum.spectrum(self.spec.waves,
                                     self.spec.flux,
                                     self.spec.error,
                                     self.x,self.y,
                                     self.N,self.b,self.z,
                                     self.wave,self.gamma,self.f,
                                     self.starts,self.ends)
        self.assertTrue(cont.shape[0]==ab.shape[0])


    def test_optimizer(self):
        self.assertTrue(self.model.error)
        cont, ab, chi2= _spectrum.spectrum(self.spec.waves,
                                     self.spec.flux,
                                     self.spec.error,
                                     self.x,self.y,
                                     self.N,self.b,self.z,
                                     self.wave,self.gamma,self.f,
                                     self.starts,self.ends)

        self.assertTrue(cont.shape[0]==ab.shape[0])
        popt, pcov = optimizer.optimize(self.spec, self.model)
        #typically, optimizer doesn't converge on a solution.  likely something
        #up with how I pass input into scipy.optimize.
        #need a better way to optimize spectra.

if __name__ == '__main__':
    unittest.main()


