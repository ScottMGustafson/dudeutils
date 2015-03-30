import data_types
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy import stats



b_max = 18.  #max allowed line width in km/s
xmlfiles = ['228.xml','091.xml','093.xml' ,'1177.xml','125.xml','227.xml',
            '230.xml','257.xml','549.xml','092.xml','1062.xml',
            '124.xml','202.xml','231.xml','548.xml','414.xml','combined.xml' ]  
#the first element will be the reference spectrum
#absorber ids to be used for the red, blue and green CCDs respectively


blu_rng = [0.,4285.]
grn_rng = [4285., 5290.]
red_rng = [5290.,7000.]
wave = data_types.atomic_data['HI'][0].wave  #wavelength all diagnostic lines will use


#test for non-zero slope
#H0: slope =0
#H0: slope is nonzero
#shuffle x, y coords,
# p = fraction of points with more extreme slopes



class CCD(object):
    def __init__(self, xml, wave_range):
        self.name = xml
        self.absorbers = CCD.get_absorbers(self.name, wave_range)
        self.wave_range = wave_range
    def shift(self,ref, stat=np.mean):
        all_shifts = [self.absorbers[i].getShift(ref.absorbers[i].z) for i in range(len(self.absorbers))]
        #print("mean, std:",np.mean(all_shifts), np.std(all_shifts))
        x=[1215.67*(1.+float(ref.absorbers[i].z)) for i in range(len(self.absorbers))]
        #print(stats.spearmanr(x,all_shifts))
        m, b, r, p, se = stats.linregress(x,all_shifts)
        if not CCD.CL(m,se, cl=.95):
            pass
            #print("slope=%lf, p=%lf"%(m,p))
            #plt.plot(x, all_shifts, 'ko')
            #plt.title(self.name+" "+str(self.wave_range))
            #plt.show()
        dist, mu, SE = CCD.bootstrap_resampling(all_shifts,stat=stat)
        std = CCD.bootstrap_resampling(all_shifts)[1]
        return mu, SE, max(np.fabs(dist)), len(all_shifts), std

    @staticmethod
    def CL(mu,SE,mu0=0.,cl=.95):
        """assuming our distribution is gaussian, is the mean significant to 95% CL?"""
        if cl==0.95:
            return -1.96*SE<=mu-mu0<=1.96*SE
        else:
            return SE<=mu-mu0<=SE

    @staticmethod
    def get_absorbers(name, wr):
        lst = data_types.Data.read(name, tag='Absorber')
        return [ item for item in lst if wr[0]<=item.get_wave()<=wr[1] ]

    @staticmethod
    def bootstrap(data):
        data=np.array(data)
        indices = np.random.randint(len(data),size=len(data))
        return data[indices]

    @staticmethod
    def bootstrap_resampling(data,samples=5000, stat=np.mean):
        
        #total = []
        bootstrap_dist = []
        for i in range(samples):
            shuffled = CCD.bootstrap(data)
            #total=np.concatenate( (total, shuffled), axis=0 )
            bootstrap_dist.append(stat(shuffled))
        return bootstrap_dist, np.mean(bootstrap_dist), np.std(bootstrap_dist)


plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=18)
plt.rc('xtick.major', size=5, pad=14)

def histogram(data, num_bins=20, xlabel=None, ylabel=None, title=None, density=False):
    """plots a histogram.  If you need any filtering of data, do that before hand"""

    hist, _x = np.histogram(data, bins=num_bins, density=False)

    x=[]    
    for i in range(len(hist)):
        x.append((_x[i]+_x[i+1])/2.)

    if density:
        hist = normalize_hist(hist)

    return x, hist

def normalize_hist(hist):
    tot = 0
    for item in hist:
        tot+=float(item)
    return [item/tot for item in hist]





def plot_bootstrap(mean_dist, std_dist, color):

    f, (ax1, ax2) = plt.subplots(1, 2, sharey=True,figsize=[8,3.5])

    #ax1 = fig1.add_subplot(gs[:,1:50])
    #ax2 = fig1.add_subplot(gs[:,50:99])

    x1, mean = histogram(mean_dist, density=True)
    x2, std = histogram(std_dist, density=True)

    ax1.plot(x1, mean,linestyle='steps')
    ax1.get_xaxis().get_major_formatter().set_useOffset(False)
    #ax1.get_yaxis().get_major_formatter().set_useOffset(False)
    for tl in ax1.get_yticklabels(): tl.set_visible(False)
    ax1.set_xlabel("mean shift ($km~s^{-1}$)")
    #ax1.set_ylabel("frequency")
    ax1.set_ylabel(color+" CCD")

    ax2.plot(x2, std, linestyle='steps')
    ax2.get_xaxis().get_major_formatter().set_useOffset(False)
    ax2.set_xlabel("$\sigma~$ ($km~s^{-1}$)")
    ax2.yaxis.set_visible(False)

    #plt.subplots_adjust(top=0.87,bottom=0.13,left=0.1,right=0.95)
    #plt.savefig(output_fname)    
    #print('n=',len(data))
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':

    red_ref = CCD(xmlfiles[0],red_rng)
    blu_ref = CCD(xmlfiles[0],blu_rng)
    grn_ref = CCD(xmlfiles[0],grn_rng)

    #take out absorbers with an inappropriately large width
    exclude_red = [item.id for item in red_ref.absorbers if float(item.b) > b_max]
    exclude_blu = [item.id for item in blu_ref.absorbers if float(item.b) > b_max]
    exclude_grn = [item.id for item in grn_ref.absorbers if float(item.b) > b_max]
     

    #take out absorbers with an inappropriately large width


    for item in xmlfiles[1:]:
        blu = CCD(item, blu_rng).shift(blu_ref)
        grn = CCD(item, grn_rng).shift(grn_ref)
        red = CCD(item, red_rng).shift(red_ref)

        stdred = CCD(item, red_rng).shift(red_ref, stat=np.std)
        stdblu = CCD(item, blu_rng).shift(blu_ref, stat=np.std)
        stdgrn = CCD(item, grn_rng).shift(grn_ref, stat=np.std)
        
        print("%s:  b, g, r\n--------------------------------------------"%(item.strip('.xml')))
        #print("vccd  : %5.3lf, %5.3lf, %5.3lf"%(blu[0],grn[0],red[0]))
        #print("CL    : %5.3lf, %5.3lf, %5.3lf"%(1.96*blu[1],1.96*grn[1],1.96*red[1]))
        #print("max   : %5.3lf, %5.3lf, %5.3lf"%(blu[2],grn[2],red[2]))
        #print("n     : %d      %d      %d"%(blu[3],grn[3],red[3]))
        #print("shift?: %5s %5s %5s"%(not CCD.CL(blu[0],blu[1]), not CCD.CL(grn[0],grn[1]), not CCD.CL(red[0],red[1])))
        print("$%5.2lf \pm %5.2lf$ & $%5.2lf \pm %5.2lf$ & $%5.2lf \pm %5.2lf$ \\\\ \n" % (blu[0],stdblu[0], grn[0],stdgrn[0], red[0],stdred[0]))


    #print("dispersion for combined:  %lf, %lf, %lf"%(stdblu, stdgrn, stdred))

    #print("overall average error for each CCD:")

    #blu_mean, bmu, bstd=CCD.bootstrap_resampling([item[0] for item in std_errs], stat=np.std)
    #grn_mean, gmu, gstd=CCD.bootstrap_resampling([item[1] for item in std_errs], stat=np.std) 
    #red_mean, rmu, rstd=CCD.bootstrap_resampling([item[2] for item in std_errs], stat=np.std)

    #blu_std, _b=CCD.bootstrap_resampling([item[0] for item in std_errs],stat=np.std)[0:2]
    #grn_std, _g=CCD.bootstrap_resampling([item[1] for item in std_errs],stat=np.std)[0:2] 
    #red_std, _r=CCD.bootstrap_resampling([item[2] for item in std_errs],stat=np.std)[0:2]

    #plot_bootstrap(blu_mean, blu_std, 'blue')
    #plot_bootstrap(grn_mean, grn_std, 'green')
    #plot_bootstrap(red_mean, red_std,'red')


    
    #print("overall means: %lf %lf %lf"%(bmu, gmu, rmu))
    #print("stds: %lf %lf %lf"%(_b, _g, _r))
    #print("std error of mean: %lf %lf %lf"%(bstd, gstd, rstd))

    #print(np.mean([item[0] for item in std_errs]),np.mean([item[1] for item in std_errs]),np.mean([item[2] for item in std_errs]))





