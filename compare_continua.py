#determine if there is a significant difference in dh vals for diff continua choices.  if so, gives you median val for each range.  From this, you can estimate the continuum contribution to systematic error.


import model
import data_types
import scipy.stats as stats
from numpy import median
import dudeutils

dbname = 'database.xml'
c = 299792.458
velrange=(-.1,.1)
alpha = 0.05

if __name__ == '__main__':
    
    db = model.ModelDB.read(dbname)
    continua = dudeutils.cont_check_pipeline()
    #continua = list(set([item.ContinuumPointList for item in db.models]))

    data = []
    models = []
    for item in db:
        it = {}
        zD = float(item.get_datum('D','Absorber','z'))
        zH = float(item.get_datum('H','Absorber','z'))
        it['vel'] = (zD-zH)*c/(1.+zH)
        if velrange[0]<=it['vel']<=velrange[1]:
            it['dh'] = float(item.get_datum('D','Absorber',param='N')) - float(item.get_datum('H','Absorber',param='N'))
            it['continuum'] = item.ContinuumPointList
            it['reduced_chi2'] = item.reduced_chi2
            data.append(it)

    for item in data:
        assert(velrange[0]<=item['vel']<=velrange[1])

    dh = []
    chi2=[]
    for cont in continua:
        dh.append([item['dh'] for item in data if item['continuum']==cont[0]])
        chi2.append([item['reduced_chi2'] for item in data if item['continuum']==cont[0]])

    #each element of dh is a list of D/H vals for each continuum
        
    W, p = stats.levene(*dh)   #defaults to median as center.
    print(p)  #determines whether variance is significant.
    if p<=alpha:
        for i in range(len(continua)):
            ind = chi2[i].index(min(chi2[i])) #index of best fit
            print('%20s med(dh)=%6.3f best dh=%6.3f red chi2=%6.3f datapoints=%4d'%(continua[i][-1],median(dh[i]),dh[i][ind],chi2[i][ind],len(dh[i])))  #prints range of continua
    else:
        print("no significant difference in DH values between chosen continua in chosen velocity range")
        
    
        
        


    


