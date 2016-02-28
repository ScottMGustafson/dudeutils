from dudeutils import populate_database
from model import Model
import numpy as np
import matplotlib.pyplot as plt
ref_wave=4847.25

def filter_cont(wave,x,cont, range_allowed):
    ind = np.argmin(np.fabs(wave-x))
    if range_allowed[0]<cont[ind]<range_allowed[1]:
        return True
    else:
        return False

def label_lines(ax, x, y, label):
    for i in range(len(label)):
        ax.annotate(label[i], xy=(x[i], y[i]), xytext=(x[i], 0.8),ha='center', va='center',
                    arrowprops=dict(arrowstyle="-",facecolor='black') )
        #plt.axvline(x=x[i], ymin=0., ymax=0.68, linewidth=1, color='k')
    return ax

def get(wave,x,attr):
    ind = np.argmin(np.fabs(wave-x))
    return attr[ind]

def get_cont(pth,wv):
    x,ab,cont = np.loadtxt(pth,usecols=(1,4,5),unpack=True)
    return get(wv,x,cont)

db=None
plt.rc('text',usetex=True)
plt.rc('font',family='serif',size=16)
fig, ax = plt.subplots()
plt.minorticks_on()

pth='/home/scott/research/2015-09-23Extraction/'


line_wave1 = [4847.234, 4847.615, 4848.438, 4848.60]
line_label1 = [r'D\,\textsc{i}','X',r'H\,\textsc{i}',r'H\,\textsc{i}']
x,y,ab,cont = np.loadtxt(pth+"a.dat",usecols=(1,2,4,5),unpack=True)
ax.plot(x,(10.**13.)*y,'k-',linestyle='steps')
ax.ticklabel_format(useOffset=False)
files=[]

#added 2015-12-18 to remove some of the most extreme continua


#plot cont at 4845.9 versus N(D)
ND=[]
contND=[]
filelst=[str(l) for l in 'abcdefghijklmnopqrstuvwxyz123456789']
filelst+=['aa','bb','cc','dd','ee','gg','hh','ii','jj','kk','ll','mm','nn','oo','pp','qq','rr','ss','tt','uu','vv','ww','xx']
for item in filelst:
    dat = item+'.dat'

    try:
        x,ab,cont = np.loadtxt(pth+dat,usecols=(1,4,5),unpack=True)
    except:
        raise Exception(str(pth+dat))

    cond = [#filter_cont(4854.6,x,ab,[0.616E-13,0.68E-13]),
            #filter_cont(4854.75,x,ab,[0.616E-13,0.66E-13]),
            #filter_cont(4850.13,x,cont,[0.616E-13,0.679E-13]),
            #filter_cont(4839.95,x,cont,[0.68E-13,0.71E-13]),
            #filter_cont(4854.69,x,ab,[0.62E-13,0.646E-13]),
            #filter_cont(4847.3,x,ab,[0.402E-13,0.410E-13]),
            #get(4845.,x,cont)>get(4854.7,x,cont)
            ]


    if all(item for item in cond):
        try:
            ND.append(Model( 
                    xmlfile=pth+item+'.xml',abs_ids=['H','D'] 
                ).get_datum("D", tag="Absorber", param="N"))
        except:
            raise Exception(str(pth+item+'.xml'))
        contND.append(get(4847.25,x,cont))
        files.append(item+'.xml')
        ax.plot(x,(10.**13.)*ab,'b--',linestyle='default')
        ax.plot(x,(10.**13.)*cont,'b-',linestyle='default')
    else:
        print(item+'.xml')


#ax=label_lines(ax, line_wave1, [7E12, 7E12, 7E12, 7E12], line_label1)
    
print(str(len(files))+" continua survived.")



models=[Model( xmlfile=pth+f, abs_ids=['H','D'] ) for f in files]
tmp=sorted(models, key=lambda x: x.chi2, reverse=True)
print([ (item.xmlfile, item.chi2) for item in tmp[:5]],[ (item.xmlfile,item.chi2) for item in tmp[-5:]])

#print([item.xmlfile for item in models[0:5]])

#for mod in models:
#    print(mod.xmlfile, mod.chi2, mod.dh)
contND=np.array(contND)
ND=np.array(ND)
dh = np.array([item.dh for item in models])
lnhi=np.array([item.get_datum('H','Absorber',param='N') for item in models])
lndi=np.array([item.get_datum('D','Absorber',param='N') for item in models])
bhi=np.array([item.get_datum('H','Absorber',param='b') for item in models])
bdi=np.array([item.get_datum('D','Absorber',param='b') for item in models])
chi2=np.array([float(item.chi2) for item in models])

minchi2=np.amin(chi2)

lst1 = list(np.where(lndi<12.83)[0])
lst2 = list(np.where(lndi>12.77)[0])
#lst3 = list(np.where(chi2<(1.0+1.1*minchi2))[0])
ind_1sig = list(set(lst1) & set(lst2))# & set(lst3))
ind_not_1sig=list(set(np.arange(lndi.shape[0]))-set(ind_1sig))
#ind_not_1sig=list(set(ind_not_1sig) & set(lst3))

print(len(ind_1sig),len(ind_not_1sig))

best_ind = list(np.where((chi2-np.amin(chi2))==0)[0])

#print("max min:",np.amax(dh),np.amin(dh))
#print("dh",np.mean(dh), np.std(dh))
#print("lnhi",np.mean(lnhi), np.std(lnhi))
#print("lndi",np.mean(lndi), np.std(lndi))
#print("bhi",np.mean(bhi), np.std(bhi))
#print("bdi",np.mean(bdi), np.std(bdi))
#print('continuum:', np.mean(contND*10.**13.),np.std(contND*10.**13.))
#print('continuum:', np.amax(contND*10.**13.),np.amin(contND*10.**13.))


msg=""
for attr in [lndi, bdi, contND*10.**13.]:
    msg += " %s & %9.4lf & %9.4lf & %9.4lf & %9.4lf & %9.4lf \\\\ \n"%(
        "attr", attr[best_ind], np.mean(attr), 
         np.std(attr),np.amin(attr),np.amax(attr)
        )

print(msg)

plt.xlabel(r"Wavelength (\AA)")
plt.ylabel(r"$F_{\lambda}~(10^{13}~$erg~s$^{-1}~$cm$^{-2}~$\r{A}$^{-1})$")
plt.xlim(4846.,4852.)
plt.ylim(-0.05,1.0)
plt.minorticks_on()
plt.show()
plt.clf()

fig, ax = plt.subplots(2, sharex=True)

ax[0].scatter(contND[ind_1sig]*10.**13.,ND[ind_1sig],s=42, c='r', marker="o", label='a')
ax[0].scatter(contND[ind_not_1sig]*10.**13.,ND[ind_not_1sig],s=42, c='k', marker="^", label='b')
ax[1].scatter(contND[ind_1sig]*10.**13.,chi2[ind_1sig],s=42, c='r', marker="o", label='a')
ax[1].scatter(contND[ind_not_1sig]*10.**13.,chi2[ind_not_1sig],s=42, c='k', marker="^", label='b')

#ax[0].xlabel()
ax[0].set_ylabel(r"log N$_{DI}$~(cm$^{-2}$)")

ax[1].set_xlabel(r"Continuum Level at 4847.25 \AA~(flux units)")
ax[1].set_ylabel(r"$\chi^2$")
ax[0].ticklabel_format(useOffset=False)
ax[0].minorticks_on()
ax[1].ticklabel_format(useOffset=False)
ax[1].minorticks_on()
plt.show()
plt.clf()

plt.plot(np.array([item.get_datum('D','Absorber',param='b') for item in models]), chi2, 'ro')
plt.show()
plt.clf()

plt.plot(np.array([item.get_datum('D','Absorber',param='N') for item in models]), chi2, 'ro')
plt.show()
plt.clf()


