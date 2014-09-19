import numpy as np
import xmlutils as dude

Ryd  = 13.60569253  #(eV)
hc   =  0.000123984193  # eV Angstroms
pi = 3.141592654

class FitData(object): 
    def __init__(self,dumpfile=dumpfile,xmlfile=xmlfile):
        out = np.loadtxt(dumpfile,unpack=True)
        self.allwaves = out[1]
        self.flux  = out[2]
        self.error = out[3]
        self.cont  = out[4]
        self.abs   = out[5]
        self.contPoints = dude.getContinuumPoints(xmlfile)
        regions=self.get_highlighted_regions(xmlfile)
        temp = []
        for item in regions:
            temp.append(np.where(self.allwaves in [item.start, item.end])[0])
        #convert to single np array.
        
    def get_highlighted_regions(self,xmlfile)
        regions = [ dude.Region(xmlnode=it) for it in dude.getList(xmlfile,"Regions") ]
        regions.sort(key=lambda x: x.start)
        consolidated = []
        for i in range(len(regions)):
            j=i
            while regions[j].start<regions[i].end: j+=1
            consolidated.append(dude.Region(start=regions[i].start,end=regions[j].end))   
          
        return consolidated
        

    def get_cut(self,start,end):
        ind = []
        for i in range(0,self.waves.shape[0]):
            if start < self.waves[i] < end:
                ind.append(i)
        return self.waves[ind], self.flux[ind], self.error[ind], self.cont[ind], self.abs[ind]

        
        

