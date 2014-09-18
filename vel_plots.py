import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import data_structures

unc=3.0     #D uncertainty in km/s
dh_unc = 0.5
dh_accepted = -4.5    #just used to help decide on a scaling
vel_tol = 5.0
vel_guess=62.0
filenames = ["../data/chi2_lo_database.txt","../data/chi2_best_database.txt","../data/chi2_hi_database.txt"]
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.grid(False)
c = 299792.458

class Data(object):
    """
    Parameters:
    --------------------
    filename : name of file containing datapoints
    nh : HI column
    nd : D column
    vel: contaminnant velocity from HI
    chi2: chi squared
    d_vel: velocity shift of D
    params: #of params
    datapts: #of datapoints
    chi2lim: max chi2 to include in plot
    """
    def __init__(self, modeldb, **kwargs):
        self.best=False

        for key,val in kwargs.items():
            setattr(self,key,val)


        self._filter()

        self.dh = np.array([item.getabs(iden='D').N-item.getabs(iden='H').N for item in self.modeldb.lst ])
        self.d_vel = np.array([self.get_vel(item,'D','H') for item in self.modeldb.lst])
        self.vel = np.array([self.get_vel(item,'D','H') for item in self.modeldb.lst])
        
        
        
    @staticmethod
    def get_vel(mod,iden1,iden2):
        z1 = mod.getabs(iden=iden1).z
        z2 = mod.getabs(iden=iden2).z
        return (z1-z2)*c/(1.+z1)

    def _filter(self):
        ind = []    # what to delete from arrays
        for i in range(self.vel.shape[0]):
            if np.fabs(self.d_vel[i])>unc and \
                     np.fabs(self.nd[i]-self.nh[i]-dh_accepted)>dh_unc and \
                     self.chi2[i] > self.chi2lim:    
                ind.append(i)

        np.delete(self.nh,ind)
        np.delete(self.nd,ind)
        np.delete(self.vel,ind)
        np.delete(self.d_vel,ind)
        np.delete(self.chi2,ind)

        #get min chi square for a given set of models
        near_min = np.where(np.fabs(self.vel-vel_guess)<vel_tol)[0].tolist()
        self.chi2 = np.array([item.chi2 for item in modeldb.lst])
        self.chi2min = np.amin(self.chi2[near_min])
        #self.onesig, self.modeldb = get_onesig(modeldb,np.amin(chi2))

        """
        onesig = []
        for i in range(self.vel.shape[0]):
            if np.fabs(self.vel[i]-vel_guess)<vel_tol and self.chi2[i]<chi2min+1.01:
                onesig.append(i)

        assert(len(onesig)>0)

        ret = {'nh':self.nh[onesig],'nd':self.nd[onesig],'d_vel':self.d_vel[onesig], \
                            'vel':self.vel[onesig],'chi2':self.chi2[onesig], \
                            'datapts':self.datapts,'params':self.params}
        np.delete(self.nh,onesig)
        np.delete(self.nd,onesig)
        np.delete(self.vel,onesig)
        np.delete(self.d_vel,onesig)
        np.delete(self.chi2,onesig)
        return Data(nh=ret['nh'], nd=ret['nd'], d_vel=ret['d_vel'], vel=ret['vel'],\
                                chi2=ret['chi2'],datapts=ret['datapts'],params=ret['params'])

        """
def get_onesig(modeldb,chi2min):
    onesig = []
    for i in range(len(modeldb.lst)):
        if modeldb.lst[i].chi2 < chi2min+1.:
            onesig.append(modeldb.pop(i))

    onesigdb = data_structures.ModelDB(onesig)
    return onesigdb, modeldb
    

def plot_chi2(ax,x,y,onesigx=None,onesigy=None):
    ax.set_ylabel(r'$\chi^2$')
    plt.setp( ax.get_xticklabels(), visible=False)
    ax.plot(x,y,     'ko')
    if onesigx is not None:
        ax.plot(onesigx, onesigy,     'co')

    return

def plot_dh(ax,x,y,onesigx=None,onesigy=None):
    ax.set_ylabel(r'log\,N$_{\rm D I}$/N$_{\rm H I}$')
    ax.plot(x,y,     'ro')
    if onesigx is not None:
        ax.plot(onesigx, onesigy,     'co')
    return

def plot_one_pair_old(data,ax1,ax2):
    if data.best:
        plot_chi2(ax1, data.vel, data.chi2, data.onesig.vel, data.onesig.chi2)
        plot_dh(    ax2, data.vel, data.nd-data.nh, data.onesig.vel, data.onesig.nd-data.onesig.nh) 
    else:
        plot_chi2(ax1, data.vel, data.chi2)
        plot_dh(    ax2, data.vel, data.nd-data.nh)
    return

def plot_one_pair(data,ax1,ax2):
    plot_chi2(ax1, data.vel, data.chi2, data.onesig.vel, data.onesig.chi2)
    plot_dh(    ax2, data.vel, data.nd-data.nh, data.onesig.vel, data.onesig.nd-data.onesig.nh) 
    return

if __name__ == "__main__":
    # instantiate three instances of Data()


    #TODO this needs to be cleaned up
    hi = Data(data_structures.read_in("chi2_hi_database.txt.new"), params=11, datapts=957, chi2lim=2000.)
    hi_onesigdb, hidb = get_onesig(hi.modeldb,hi.chi2min)    
    hi = Data(hidb, params=11, datapts=957, chi2lim=2000.)
    hi_onesig = Data(hi_onesigdb, params=11, datapts=957, chi2lim=2000.)

    best = Data(data_structures.read_in("test_database.txt.new"), params=8, datapts=957, chi2lim=1825., best=True)
    best_onesigdb, bestdb = get_onesig(best.modeldb,best.chi2min) 
    best = Data(bestdb, params=8, datapts=957, chi2lim=1825., best=True)
    best_onesig = Data(best_onesigdb, params=8, datapts=957, chi2lim=1825., best=True)
 
    lo = Data(data_structures.read_in("chi2_lo_database.txt.new"), params=11, datapts=957, chi2lim=2000.)
    lo_onesigdb, lodb = get_onesig(lo.modeldb,lo.chi2min)  
    lo = Data(lodb, params=11, datapts=957, chi2lim=2000.)
    lo_onesig = Data(lo_onesigdb, params=11, datapts=957, chi2lim=2000.)

    f = plt.figure(figsize=(4.1,9.3))

    gs1 = GridSpec(2,1)
    gs1.update(left=0.2, right=0.9, top=0.95, bottom=0.68, hspace=0.1)
    ax1 = plt.subplot(gs1[:-1, 0])
    ax2 = plt.subplot(gs1[-1, 0])
    ax1.set_ylim([1980.,2000.])
    ax2.set_ylim([-5.0,-4.4])
    plt.setp( ax2.get_xticklabels(), visible=False)
    plot_one_pair(hi,ax1,ax2)

    gs2 = GridSpec(2,1)
    gs2.update(left=0.2, right=0.9, top=0.64, bottom=0.36, hspace=0.1)
    ax3 = plt.subplot(gs2[:-1, 0])
    ax4 = plt.subplot(gs2[-1, 0])
    ax3.set_yticks(np.arange(1800.,1825.,5.))
    ax3.set_ylim([1795.,1825.])
    ax4.set_ylim([-5.0,-4.4])
    plt.setp( ax4.get_xticklabels(), visible=False)
    plot_one_pair(best,ax3,ax4)

    gs3 = GridSpec(2,1)
    gs3.update(left=0.2, right=0.9, top=0.32, bottom=0.05, hspace=0.1)
    ax5 = plt.subplot(gs3[:-1, 0])
    ax6 = plt.subplot(gs3[-1, 0])
    ax5.set_ylim([1980.,2000.])
    ax6.set_ylim([-5.0,-4.4])
    plot_one_pair(lo,ax5,ax6)

    ax6.set_xlabel(r"velocity from {\textrm{H}\,\textsc{i}~} (km s$^{-1}$)")

    plt.savefig('~/Desktop/hi_best_lo.png',dpi=600)
    plt.clf()

    #best
    ######################################################################
    f, (chi2_ax, dh_ax) = plt.subplots(2, sharex=True, figsize=[3.5,4])
    chi2_ax.get_yaxis().get_major_formatter().set_useOffset(False)
    dh_ax.set_xticks(np.arange(61.80,62.15,0.1))
    #chi2_ax.set_xticks(np.arange(61.80,62.15,0.1))
    #plt.setp( chi2_ax.get_xticklabels(), visible=False)
    dh_ax.set_ylim( [-4.57,-4.535])
    #chi2_ax.set_xlim([61.80,62.15])
    dh_ax.set_xlim(    [61.80,62.10])
    chi2_ax.set_ylim([1802.75,1804.5])
    chi2_ax.set_yticks(np.arange(1803.,1804.5,0.25))

    plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
    dh_ax.set_ylabel(r'log\,N$_{\rm D I}$/N$_{\rm H I}$')
    chi2_ax.set_ylabel(r'$\chi^2$')
    dh_ax.set_xlabel(r"velocity from {\rm H I} ($km s^{-1}$)")
    #plt.setp( chi2_ax.get_xticklabels(), visible=False)
    #plt.setp( dh_ax.get_xticklabels(), visible=True)


    chi2_ax.plot(best.vel,best.chi2,     'ko')
    chi2_ax.plot(best.onesig.vel, best.onesig.chi2,     'co')

    dh_ax.plot(best.vel,best.nd-best.nh,     'ko')
    dh_ax.plot(best.onesig.vel, best.onesig.nd-best.onesig.nh,     'co')
    plt.gcf().tight_layout()
    plt.savefig("onesig_best.png")


 #lo
    ######################################################3

    plt.clf()
    f, (chi2_ax, dh_ax) = plt.subplots(2, sharex=True, figsize=[3.5,4])
    chi2_ax.get_yaxis().get_major_formatter().set_useOffset(False)
    #dh_ax.set_xticks(np.arange(61.80,62.15,0.1))
    #chi2_ax.set_xticks(np.arange(61.80,62.15,0.1))
    #plt.setp( chi2_ax.get_xticklabels(), visible=False)
    dh_ax.set_ylim( [-4.58,-4.55])
    #chi2_ax.set_xlim([61.80,62.15])
    dh_ax.set_xlim(    [61.80,62.40])
    chi2_ax.set_ylim([1980.,1982.])
    #chi2_ax.set_yticks(np.arange(1803.,1804.5,0.25))

    plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
    dh_ax.set_ylabel(r'$log(N_D/N_H)$')
    chi2_ax.set_ylabel(r'$\chi^2$')
    dh_ax.set_xlabel(r"velocity from HI ($km s^{-1}$)")
    #plt.setp( chi2_ax.get_xticklabels(), visible=False)
    #plt.setp( dh_ax.get_xticklabels(), visible=True)


    chi2_ax.plot(lo.vel,lo.chi2,     'ko')
    chi2_ax.plot(lo.onesig.vel, lo.onesig.chi2,     'co')

    dh_ax.plot(lo.vel,lo.nd-lo.nh,     'ko')
    dh_ax.plot(lo.onesig.vel, lo.onesig.nd-lo.onesig.nh,     'co')
    plt.gcf().tight_layout()
    plt.savefig("onesig_lo.png")


    #hi
    plt.clf()
    f, (chi2_ax, dh_ax) = plt.subplots(2, sharex=True, figsize=[3.5,4])
    chi2_ax.get_yaxis().get_major_formatter().set_useOffset(False)
    #dh_ax.set_xticks(np.arange(61.80,62.15,0.1))
    #chi2_ax.set_xticks(np.arange(61.80,62.15,0.1))
    #plt.setp( chi2_ax.get_xticklabels(), visible=False)
    dh_ax.set_ylim( [-4.58,-4.54])
    #chi2_ax.set_xlim([61.80,62.15])
    dh_ax.set_xlim(    [61.60,62.50])
    chi2_ax.set_ylim([1980.,1982.])
    #chi2_ax.set_yticks(np.arange(1803.,1804.5,0.25))

    plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
    dh_ax.set_ylabel(r'$log(N_D/N_H)$')
    chi2_ax.set_ylabel(r'$\chi^2$')
    dh_ax.set_xlabel(r"velocity from HI (km s$^{-1}$)")
    #plt.setp( chi2_ax.get_xticklabels(), visible=False)
    #plt.setp( dh_ax.get_xticklabels(), visible=True)


    chi2_ax.plot(hi.vel,hi.chi2,     'ko')
    chi2_ax.plot(hi.onesig.vel, hi.onesig.chi2,     'co')

    dh_ax.plot(hi.vel,hi.nd-hi.nh,     'ko')
    dh_ax.plot(hi.onesig.vel, hi.onesig.nd-hi.onesig.nh,     'co')
    plt.gcf().tight_layout()
    plt.savefig("onesig_hi.png")

