from dudeutils import *
import numpy as np
import matplotlib.pyplot as plt
db = load_from_db("2014-10-06db.xml")

z=2.9884144
vel=4.0
delz = (vel/299792.458)*(1.+z)

constraints = {"D":{"z":(z-delz,z+delz)},"chi2":1890}
db.best_fit("H2","z",4,constraints=constraints)


