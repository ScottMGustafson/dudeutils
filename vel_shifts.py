import data_types
import numpy as np



b_max = 18.  #max allowed line width in km/s
xmlfiles = ['228.xml','091.xml','093.xml' ,'1177.xml','125.xml','227.xml',
            '230.xml','257.xml','549.xml','092.xml','1062.xml',
            '124.xml','202.xml','231.xml','548.xml','414.xml' ]  
#the first element will be the reference spectrum
#absorber ids to be used for the red, blue and green CCDs respectively


blu_rng = [0.,4285.]
grn_rng = [4285., 5290.]
red_rng = [5290.,7000.]
wave = data_types.atomic_data['HI'][0].wave  #wavelength all diagnostic lines will use

class CCD(object):
    def __init__(self, xml, wave_range):
        self.name = xml
        self.absorbers = CCD.get_absorbers(self.name, wave_range)

    def shift(self,ref):
        all_shifts = [self.absorbers[i].getShift(ref.absorbers[i].z) for i in range(len(self.absorbers))]
        dist, mu, SE = CCD.bootstrap_resampling(all_shifts, plot=False)
        return mu, SE, max(np.fabs(dist)), len(all_shifts), np.std(all_shifts)

    @staticmethod
    def CL(mu,SE):
        """assuming our distribution is gaussian, is the mean significant to 95% CL?"""
        return -1.96*SE<=mu<=1.96*SE

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
    def bootstrap_resampling(data,samples=1500,plot=True):
        import histogram
        #total = []
        bootstrap_dist = []
        for i in range(samples):
            shuffled = CCD.bootstrap(data)
            #total=np.concatenate( (total, shuffled), axis=0 )
            bootstrap_dist.append(np.mean(shuffled))

        if plot: histogram.histogram(bootstrap_dist)

        return bootstrap_dist, np.mean(bootstrap_dist), np.std(bootstrap_dist)


if __name__ == '__main__':

    red_ref = CCD(xmlfiles[0],red_rng)
    blu_ref = CCD(xmlfiles[0],blu_rng)
    grn_ref = CCD(xmlfiles[0],grn_rng)

    #take out absorbers with an inappropriately large width
    exclude_red = [item.id for item in red_ref.absorbers if float(item.b) > b_max]
    exclude_blu = [item.id for item in blu_ref.absorbers if float(item.b) > b_max]
    exclude_grn = [item.id for item in grn_ref.absorbers if float(item.b) > b_max]
     

    #take out absorbers with an inappropriately large width

    std_errs = []

    for item in xmlfiles[1:]:
        red = CCD(item, red_rng).shift(red_ref)
        blu = CCD(item, blu_rng).shift(blu_ref)
        grn = CCD(item, grn_rng).shift(grn_ref)
        
        std_errs.append((blu[-1],grn[-1],red[-1]))
        print("%s:  b, g, r\n--------------------------------------------"%(item.strip('.xml')))
        print("vccd  : %5.3lf, %5.3lf, %5.3lf"%(blu[0],grn[0],red[0]))
        print("CL    : %5.3lf, %5.3lf, %5.3lf"%(1.96*blu[1],1.96*grn[1],1.96*red[1]))
        print("max   : %5.3lf, %5.3lf, %5.3lf"%(blu[2],grn[2],red[2]))
        print("n     : %d      %d      %d"%(blu[3],grn[3],red[3]))
        print("shift?: %5s %5s %5s"%(not CCD.CL(blu[0],blu[1]), not CCD.CL(grn[0],grn[1]), not CCD.CL(red[0],red[1])))
        print("%5.3lf \pm %5.3lf & %5.3lf \pm %5.3lf & %5.3lf \pm %5.3lf \\\\ \n" % (blu[0],1.96*blu[1], grn[0],1.96*grn[1], red[0],1.96*red[1]))


    print("overall average error for each CCD:")


    print(CCD.bootstrap_resampling([item[0] for item in std_errs])[1],CCD.bootstrap_resampling([item[1] for item in std_errs])[1],CCD.bootstrap_resampling([item[2] for item in std_errs])[1])


    print(np.mean([item[0] for item in std_errs]),np.mean([item[1] for item in std_errs]),np.mean([item[2] for item in std_errs]))





