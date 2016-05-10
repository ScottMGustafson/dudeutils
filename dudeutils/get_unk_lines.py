from atomic import SpectralLine
from dudeutils import newdb
from model import Model
from numpy import fabs, log10

c=299792.458
xmlfile="/home/scott/research/2015-09-23Extraction/best.xml"

vel_tol=0.5  #using diff in wv bc vels seem off

def get_vel(desired_obs, rest, z):
    z2=float(desired_obs)/float(rest)-1.
    return c*(z-z2)/(1.+z)
    
def vel_wv(wv1,wv2):
    return (wv2-wv1)*c/wv1

def get_rest(z, observed):
    return float(observed)/(1.+float(z))

def get_obs(z, rest):
    return (1.+float(z))*float(rest)

def get_z(rest, obs):
    return float(obs)/float(rest)-1.

def check_DH(line1, line2, vel_sep=81.55, vel_tol=5.0, dh=-4.6, dh_tol=0.5):
    if vel_sep-vel_tol<=c*(line2.z-line1.z)/(1.+line1.z)<=vel_tol+vel_sep:
        if fabs(line2.N-line1.N-dh)<=dh_tol:
            return True
    return False
    
def convert_to_column(abs1, abs2, N1):
    """
    convert the column of a line treated as HI to the column it would have at a 
    specified line.

    abs1,2 = SpectralLine instance
    typically, abs1=HI and abs2= our unidentified line

    ew=N*f*(e*lambda)**2 /(4*espilon_0*m*c**2)  
    so can roughly convert from one atom's column to another with
        
    (N1*f1*lambda1**2) / m1 = (N2*f2*lambda2**2) / m2

    so N2 should roughly be...
    N2=N1*(m2/m1)*(f1/f2)*(lambda1/lambda2)**2
    """
   

    if abs1.ionName ==abs2.ionName:
        return N1

    if abs1.ionName=="H I" and abs2.ionName in ["M1", "M I"] or \
        abs2.ionName=="H I" and abs1.ionName in ["M1", "M I"]:
        return N1

    m1=abs1.mass(units='amu')
    m2=abs2.mass(units='amu')

    return log10((10.**N1)*(m2/m1)*(abs1.f/abs2.f)*(abs1.wave/abs1.wave)**2.)


if __name__ == "__main__":
    model=Model(xmlfile=xmlfile)
    all_abs=Model.get(model.AbsorberList)
    
    all_HI = [item for item in all_abs if item.ionName=="H I"]
    all_M1 = [item for item in all_abs if item.ionName=="M1"]
    #rest wave if part of system at given z

    dh_candidates=[]
    #now check for possible D:
    for ref in all_HI:
        for ab in list(all_HI+all_M1):
            if check_DH(ref,ab):
                dh_candidates.append(ref,ab)

    if len(dh_candidates)>0:
        print('dh candidates:\n----------------------------------\n')
        for i in range(len(dh_candidates)):
            print(str(dh_candidates[i][0])+'\n',+str(dh_candidates[i][1])+'\n\n')

    candidates = []
    for m1 in all_M1:
        for hi in all_HI:
            for key, absorber in atomic_data.items():
                for line in absorber:
                    #get fitting lines assuming z of unk line is z_HI
                    wv=m1.get_wave(n=0)
                    
                    wv_=get_obs(hi.z,line.wave)
                    vel=c*(wv-wv_)/wv
                    column=convert_to_column(line,atomic_data['H I'][0],m1.N)
        
                    if fabs(wv-wv_)<=vel_tol and 10.<column<hi.N:
                        if line.ionName not in ["H I", "M1", "M 1", "D I"]:
                            candidates.append({
                                'name':line.ionName,
                                'rest':line.wave,
                                'z':hi.z,
                                'vel':wv-wv_,
                                'N':column,
                                'obs':get_obs(hi.z, line.wave),   #if vel is non zero, will be off expected line
                                'diff':hi.N-column
                                })

    candidates=sorted(candidates, key=lambda x:fabs(x["vel"]))
    for item in candidates:
        print("%6s:   rest=% 7.2lf     obs=% 7.2lf    z=% 10.8lf   vel=%  5.3lf   N=% 5.3lf  NHI-NX=% 5.3lf"%(item["name"], item["rest"], item['obs'], item["z"],item["vel"],item["N"],item["diff"]))


