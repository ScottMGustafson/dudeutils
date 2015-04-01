from data_types import atomic_data
from dudeutils import newdb
from model import Model
from numpy import fabs

c=299792.458
xmlfile="/home/scott/research/J0815+2710/J0815+2710.xml"

vel_tol=10.

def get_vel(desired_obs, rest, z):
    z2=float(desired_obs)/float(rest)-1.
    return c*(z-z2)/(1.+z)
    
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

if __name__ == "__main__":
    model=newdb(xmlfile=xmlfile,chi2=0.,pixels=0.,params=0.)[0]
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
                    if fabs(vel)<=vel_tol:
                        candidates.append({
                            'name':line.ionName,
                            'rest':line.wave,
                            'z':hi.z,
                            'vel':vel,
                            'N':m1.N,
                            'obs':get_obs(hi.z, line.wave)   #if vel is non zero, will be off expected line
                            })

    candidates=sorted(candidates, key=lambda x:fabs(x["vel"]))
    for item in candidates:
        print("%6s:   rest=% 7.2lf     obs=% 7.2lf    z=% 10.8lf   vel=%  5.3lf   N=%5.3lf"%(item["name"], item["rest"], item['obs'], item["z"],item["vel"],item["N"]))


