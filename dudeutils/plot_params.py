import random_sampling
import dudeutils
from plot_distribution import plot_chi2
import data_types


if __name__=="__main__":
    dct=random_sampling.parse_config()
    glob=dct.pop('config',None)

    db=random_sampling.get_nsigma(dudeutils.load_models("DHdb.obj"),n=4)


    for key, val in dct.items():
        for attr,rng in val.items():
            plot_chi2(db, iden=key, attr=attr,
                      xlabel=r"$%s(%s)$"%(attr,key))

            
