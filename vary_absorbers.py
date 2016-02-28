import dudeutils
import spectrum

def get_ab(ab_lst,iden):
    for i in range(len(ab_lst)):
        if ab_lst[i].id==iden:
            return i
    return None

model_lst=dudeutils.newdb(xmlfile,0,0)
model=model_lst[0]
model.RegionList=data_types.RegionList.consolidate_regions(model.RegionList)

spectrum = spectrum.Spectrum(fits,dumpfile)
wavelims = [ [item.start, item.end] for item in model.RegionList]
chi2_original=spectrum.chi2(model.AbsorberList,wavelims)

model.chi2=chi2_original
model_lst[0]=model

#now vary the params:
param_space = {
                'absorber_id':{
                                'N':np.arange(12,15,0.2),
                                'b':np.arange(10,15,0.5),
                                'z':np.arange(2.7182,2.7184,0.000001)
                              },
                'absorber_id2':{
                                'N':np.arange(12,15,0.2),
                                'b':np.arange(10,15,0.5),
                                'z':np.arange(2.7182,2.7184,0.000001)
                              }
                }

## no step through the provided parameter space

#chi2, absorberList list:
varied_fits=[]

ab_lst=model.AbsorberList
for key, val in param_space.items():
    for k, v in val.items():
        i=get_ab(ab_lst,key)
        for i in range(len(v)):
            setattr(ab_lst,k,v)
            chi2=spectrum.chi2(ab_lst,wavelims)
            varied_fits.append([chi2, ab_lst])

