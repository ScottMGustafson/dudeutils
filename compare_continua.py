#determine if there is a significant difference in dh vals for diff continua choices.  if so, gives you median val for each range.  From this, you can estimate the continuum contribution to systematic error.


import model
import data_types
import scipy.stats as stats
from numpy import median, array
import dudeutils
import histogram

dbname = 'database.xml'
c = 299792.458
velrange=(-.5,.5)
alpha = 0.05

def get_best_cont(continua, db):
    old_min=10E10
    for cont in continua:
        mods = [item for item in db if item.ContinuumPointList==cont[0] and velrange[0]<=item.get_shift('D','H')<=velrange[-1]]
        _min = min([item.reduced_chi2 for item in mods])
        if _min<old_min:
            old_min=_min
            _min = min([float(item.chi2) for item in mods])
            best = [item for item in mods if float(item.chi2)<=float(_min)+5.]
            
    return best

if __name__ == '__main__':
    
    db = model.ModelDB.read(dbname)
    continua = dudeutils.cont_check_pipeline(reduced_chi2_limit=1.8,verbose=False)
    #continua = list(set([item.ContinuumPointList for item in db.models]))

    data = []
    models = []
    for item in db:
        it = {}
        if velrange[0]<=item.get_shift('D','H')<=velrange[1]:
            it['dh'] = item.dh
            it['continuum'] = item.ContinuumPointList
            it['reduced_chi2'] = item.reduced_chi2
            it['name'] = item.xmlfile 
            data.append(it)

    dh = []
    chi2=[]
    names = []
    for cont in continua:
        _dh=[item['dh'] for item in data if item['continuum']==cont[0]]
        dh.append(_dh)

        #histogram.histogram(array(_dh),title=cont[1])
        #histogram.bootstrap_resampling(array(_dh), samples=5000,title=cont[1])
        chi2.append([item['reduced_chi2'] for item in data if item['continuum']==cont[0]])
        names.append([item['name'] for item in data if item['continuum']==cont[0]])

    #each element of dh is a list of D/H vals for each continuum
        
    W, p = stats.levene(*dh)   #defaults to median as center.
    print(p)  #determines whether variance is significant.
    if p<=alpha:
        for i in range(len(continua)):
            ind = chi2[i].index(min(chi2[i])) #index of best fit
            print('%20s med(dh)=%6.3f best dh=%6.3f red chi2=%6.3f datapoints=%4d name=%s'%(continua[i][-1],median(dh[i]),dh[i][ind],chi2[i][ind],len(dh[i]),names[i][ind]))  #prints range of continua
    else:
        print("no significant difference in DH values between chosen continua in chosen velocity range")
        for i in range(len(continua)):
            ind = chi2[i].index(min(chi2[i])) #index of best fit
            print('%20s med(dh)=%6.3f best dh=%6.3f red chi2=%6.3f datapoints=%4d name=%s'%(continua[i][-1],median(dh[i]),dh[i][ind],chi2[i][ind],len(dh[i]),names[i][ind]))  #prints range of continua

    _data = []
    for item in dh:
        _data+=list(item)
    print(_data)
    histogram.bootstrap_resampling(array(_data), samples=5000,title='all data (5000 samples)')

    best = get_best_cont(continua, db)
    dh = [item.dh for item in best]
    histogram.histogram(array(dh),title=best[0].xmlfile)
    histogram.bootstrap_resampling(array(dh), samples=5000,title=best[0].xmlfile)
    
    
        

    
        
        


    



