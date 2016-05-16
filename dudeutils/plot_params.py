import dudeutils.random_sampling as random_sampling
from dudeutils.utilities import load_models
from dudeutils.plot_distribution import plot_chi2
import dudeutils.data_types as data_types


if __name__=="__main__":
    dct=random_sampling.parse_config()
    glob=dct.pop('config',None)

    db=random_sampling.get_nsigma(load_models("DHdb.obj"),n=4)


    for key, val in dct.items():
        for attr,rng in val.items():
            plot_chi2(db, iden=key, attr=attr,
                      xlabel=r"$%s(%s)$"%(attr,key))

            
