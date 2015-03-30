import dudeutils
from model import c, ModelDB
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from constraints import Constraint
from matplotlib import rc


#constraints_1sig = {'D':{'b':(13.8,15.2), 'N':(12.82,12.88), 'shift':(-0.5,0.5)}, 
#               'H2':{'b':(9.0,9.6), 'shift':(-61.3,-60.2), 'N':(12.62,12.68)}, 
#               'xmlfile':"2015-01-22_best.xml"}
constraints_1sig = {'D':{'b':(13.0,14.6), 'N':(12.826,12.852), 'shift':(-0.3,0.3)}, 
               'H2':{'b':(8.9,11.), 'shift':(-61.4,-60.3), 'N':(12.62,12.68)}, 
               'xmlfile':"2015-01-22_best.xml"}

constraints = {'D':{'b':(12.,16.), 'N':(12.8,12.88), 'shift':(-0.5,0.5)}, 
               'H2':{'b':(0.,12.), 'shift':(-62.,-59.5), 'N':(12.55,12.75)}, 
              'xmlfile':"2015-01-22_best.xml",'chi2':1790.}


def xy(x_id, y_id, x,y, dbfile, constraints=None, chi2=False):
    if type(dbfile) is str:
        db = dudeutils.getdb(dbfile)
    else:
        db=dbfile
    _x, chi2 = db.get_all_abs(x_id, str(x), constraints=constraints)
    if y_id=="chi2":
        _y=chi2
    else:
        _y, _    = db.get_all_abs(y_id, str(y), constraints=constraints)
        assert(chi2==_ and len(_x)==len(_y))
    return _x, _y, chi2
    
def make_subplot(x,y, style='bo', xname=None,yname=None, ax=None,xlim=None,ylim=None):
    if ax is None:
        ax = plt.gca()
    line, = ax.plot(x, y, style)
    
    if yname: ax.set_ylabel(yname)
    if xname: ax.set_xlabel(xname)
    if not xlim is None: ax.set_xlim(xlim)
    if not ylim is None: ax.set_ylim(ylim)
    ax.ticklabel_format(useOffset=False)
    return line
       
def plotDH(db,constraints=None,the_title=None,dh_lim=None,chi2_lim=None):
    ND, NH, chi2 = xy("D","H","N","N",db,constraints=constraints)
    zH, zD, _ =xy("D","H","z","z",db,constraints=constraints) 

    DH = [float(ND[i])-float(NH[i]) for i in range(len(NH))]
    vel = [(float(zD[i])-float(zH[i]))*c/(1.+float(zH[i])) for i in range(len(zH))]
    chi2 = [float(chi2[i]) for i in range(len(zH))]
    #chi2 = list(map(float, chi2))
    assert(len(DH)==len(vel)==len(chi2))
    fig1, (ax1, ax2) = plt.subplots(nrows=2)
    
    _1 = make_subplot(vel, DH, ax=ax1, yname='D/H',ylim=dh_lim)
    _2 = make_subplot(vel, chi2, ax=ax2, xname='D velocity (km/s)', yname='chi2',ylim=chi2_lim)

    if the_title: plt.title(the_title)
    #fig2 = plt.figure()
    #plot(x, np.cos(x))
    plt.show()


def plot_cont_vel(db,constraints=None,the_title=None):
    zH2, zH, chi2 = xy("H2","H","z","z",db,constraints=constraints)
    zD, zH, _ = xy("D","H","z","z",db,constraints=constraints)

    vel2 = [(float(zH2[i])-float(zH[i]))*c/(1.+float(zH[i])) for i in range(len(zH))]
    velD = [(float(zD[i])-float(zH[i]))*c/(1.+float(zH[i])) for i in range(len(zH))]
    chi2 = [float(chi2[i]) for i in range(len(zH))]
    fig1, (ax1, ax2) = plt.subplots(nrows=2)
    
    _1 = make_subplot(vel2, velD, ax=ax1, yname='D velocity offset (km/s)')#,xlim=[-1.5,1.5])
    _2 = make_subplot(vel2, chi2, ax=ax2, xname='H2 velocity (km/s)', yname='chi2')#,xlim=[-1.5,1.5])
    if the_title: plt.title(the_title)
    #fig2 = plt.figure()
    #plot(x, np.cos(x))
    plt.show()

def plotND_NH(db,constraints=None,the_title=None):
    ND, NH, chi2 = xy("D","H","N","N",db,constraints=constraints)
    zH, zD, _ =xy("D","H","z","z",db,constraints=constraints) 

    ND = [float(ND[i]) for i in range(len(ND))]
    NH = [float(NH[i]) for i in range(len(NH))]
    chi2 = [float(chi2[i]) for i in range(len(zH))]
    fig1, (ax1, ax2) = plt.subplots(nrows=2)
    
    _1 = make_subplot(ND, NH, ax=ax1, yname='NH',xlim=(12.8,12.95))
    _2 = make_subplot(ND, chi2, ax=ax2, xname='ND', yname='chi2',xlim=(12.8,12.95))
    if the_title: plt.title(the_title)
    #fig2 = plt.figure()
    #plot(x, np.cos(x))

    plt.show()

#vel = -61.3  --  -60.1
#N  12.61 -- 12.71
# b 8.7--11.4


# bD  = 13.02 -- 13.94
#ND   12.835 -- 12.840

def get_one_sig(data):
    #constraints = Constraint(**{"D":{"NLocked":False},"H":{"NLocked":False},"chi2":(0.,1810.)})

    if not type(data) is list:
        data = data.models#[item for item in db if int(item.params)==21] 
    data = ModelDB.filter(data, constraints)


    if len(data)<10:
        raise Exception('too few models survived: %d'%(len(data)))
    #data = [item for item in data if -0.5<=item.get_shift("D","H")<=0.5] 
    #data = [item for item in data if 13.017<=float(item.get_datum('D','Absorber','b'))<=13.95]
    #data = [item for item in data if 12.835<=float(item.get_datum('D','Absorber','b'))<=13.840]

    data = sorted(data, key=lambda x: float(x.reduced_chi2))
    #data = [item for item in data if item.xmlfile =="2015-01-22_best.xml"]

    #get a rough 1-sigma range before additional filtering
    get_range(data,'N','H2')
    get_range(data,'b','H2')
    get_range(data,'z','H2')
    get_range(data,'b','D')
    get_range(data,'N','D')

    try:
        min_ = float(data[0].reduced_chi2)
        minchi2=float(data[0].chi2)
    except IndexError:
        return data, data

    best = []
    out = []
    for item in data:
        df = float(item.pixels)-float(item.params)-1.
        
        if item.get_datum('D','Absorber','zLocked'):
            if float(item._chi2)<=minchi2+2.30:
                best.append(item)
            else:
                out.append(item)
        else:
            if float(item._chi2)<=minchi2+1.:
                best.append(item)
            else:
                out.append(item)
    
    print(str(best[0]))
    
    best = ModelDB.filter(best, constraints_1sig)

    return best, out




def get_range(lst,param,id):

    out=[]
    lst = [item for item in lst if item.get_datum(id,'Absorber',param+"Locked")]
    lst = sorted(lst, key=lambda x:float(x.chi2))
    minchi2=float(lst[0].chi2)
   
    for item in lst:
        if item.get_datum('D','Absorber','zLocked'):
            if float(item.chi2)<=minchi2+2.30:
                if param=='z':
                    out.append(float(item.get_shift(id,'H')))
                else:
                    out.append(float(item.get_datum(id,'Absorber',param)))
        else:
            if float(item.chi2)-1.<=minchi2+1.:
                if param=='z':
                    out.append(float(item.get_shift(id,'H')))
                else:
                    out.append(float(item.get_datum(id,'Absorber',param)))
        
    print(param,id,len(out),min(out),max(out))

    

def plot_abs_z(id,ref_id,vr=[-0.5,0.5]):
    """plots chi square versus velocity of absorber.  """
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size=18)
    from matplotlib.gridspec import GridSpec

    f = plt.figure(figsize=(6.,7.))
#hi
    gs1 = GridSpec(2,1)
    gs1.update(left=0.2, right=0.9, top=0.9, bottom=0.1, hspace=0.1)
    ax1 = plt.subplot(gs1[:-1, 0])
    ax2 = plt.subplot(gs1[-1, 0])

    db = dudeutils.load_from_db("database.xml")
#"z":(2.988399,2.98842),
    best, data = get_one_sig(db,'H2')
  

    ax1.plot([item.get_shift(id, ref_id) for item in data],[item.dh for item in data],'ko')
    ax1.plot([item.get_shift(id, ref_id) for item in best],[item.dh for item in best],'co')  
    ax2.plot([item.get_shift(id, ref_id) for item in data],[item.reduced_chi2 for item in data],'ko')
    ax2.plot([item.get_shift(id, ref_id) for item in best],[item.reduced_chi2 for item in best],'co')
    ax1.set_ylabel("$log(D/H)$")
    ax2.set_xlabel("shift from H I ($km~s^{-1}$)")
    ax2.set_ylabel("$\chi^2/df$")
    plt.setp( ax1.get_xticklabels(), visible=False)
    plt.show()

def plot_vel_v_dh(dbfile="database.xml"):
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size=18)
    db = dudeutils.load_from_db(dbfile)
    data = ModelDB.filter(db.models, constraints)
    for item in db.models:
        item._chi2=item.chi2

    data = sorted(data, key=lambda x: x.chi2)

    best, data = get_one_sig(data)

    assert(len(best)>1)

    #data = [item for item in data if 12.60<item.get_datum(tag="Absorber",id="H2",param="N")<12.69]
    #best = [item for item in best if 12.60<item.get_datum(tag="Absorber",id="H2",param="N")<12.69]
    #data = [item for item in data if item.get_datum(tag="Absorber",id="H2",param="b")<12.]
    #best = [item for item in best if item.get_datum(tag="Absorber",id="H2",param="b")<12.]


    """
    for item in data:
        if not item.get_datum('D','Absorber','zLocked'):
            if float(item.params)==26.:
                item._chi2=float(item.chi2)-1.
    for item in best:
        if not item.get_datum('D','Absorber','zLocked'):
            if float(item.params)==26.:
                item._chi2=float(item.chi2)-1.
    """

    dh_data=[item.dh for item in data]
    dh_best=[item.dh for item in best]
    chi2_data=[float(item._chi2) for item in data]
    chi2_best=[float(item._chi2) for item in best]
    plt.plot(dh_data, chi2_data, 'ko')
    plt.plot(dh_best, chi2_best, 'co')
    plt.xlabel('$D/H$')
    plt.ylabel('$\chi^2$')
    plt.gca().ticklabel_format(useOffset=False)
    plt.savefig('DH_chi2.png')
    plt.clf()


    for param in ['N', 'b', 'z']:
        for id in ['H', 'D', 'H2']:
            if param=='z':
                if id=='H':
                    continue
                else:
                    xdata= [float(item.get_shift(id,'H')) for item in data]
                    xbest=[float(item.get_shift(id,'H')) for item in best]
            else:
                xdata= [float(item.get_datum(tag="Absorber",id=id,param=param)) for item in data]
                xbest=[float(item.get_datum(tag="Absorber",id=id,param=param)) for item in best]

            plt.plot(xdata, chi2_data, 'ko')
            plt.plot(xbest, chi2_best, 'co')
            if id=='H2':
                id='contaminant'
            if param=='z':
                plt.xlabel('velocity shift of '+id+' ($km~s^{-1}$)')
            elif param=='N':
                plt.xlabel('$log('+param+'_{'+id+'})$')
            elif param=='b':
                plt.xlabel('$b_{'+id+'}~(km~s^{-1})$')
            else:
                plt.xlabel('$'+param+'_{'+id+'}$')
            plt.ylabel('$\chi^2$')
            plt.gca().ticklabel_format(useOffset=False)
            #plt.title('$'+param+'_{'+id+'}$ vs. $\chi^2$')
            plt.savefig(id+"_"+param+'_chi2.png')
            plt.clf()

            plt.plot(xdata, dh_data, 'ko')
            plt.plot(xbest, dh_best, 'co')
            #plt.title('$'+param+'_{'+id+'}$ vs. $D/H$')
            if param=='z':
                plt.xlabel('velocity shift of '+id+' ($km~s^{-1}$)')
            elif param=='N':
                plt.xlabel('$log('+param+'_{'+id+'})$')
            elif param=='b':
                plt.xlabel('$b_{'+id+'}~(km~s^{-1})$')
            plt.ylabel('$D/H$')
            plt.gca().ticklabel_format(useOffset=False)
            plt.savefig(id+"_"+param+'_dh.png')

            plt.clf()




if __name__ == '__main__':


    #constraints = Constraint(**{"D":{"z":(2.988399,2.988424),"b":{0.,25.}}, "H2":{"z":(2.98759,2.987637), "b":(0.,10.)},"chi2":(0.,1830.)})
    #constraints = Constraint(**{"H2":{"z":(2.98759,2.987637)}})

    #db=dudeutils.load_from_db('database.xml')
    #plotDH(db,constraints=constraints)
    #plotND_NH(db,constraints=constraints)
    #plot_cont_vel(db,constraints=constraints)
    #plot_vel_v_dh()
    #plot_abs_z("H2","H")


