from dudeutils import *
import numpy as np
import matplotlib.pyplot as plt
db = load_from_db("2014-10-06db.xml")

z=2.9884144
vel=4.
delz = (vel/299792.458)*(1.+z)

constraints = {"D":{"z":(z-delz,z+delz)},"H2":{"z":(2.9865,2.9877)},"chi2":1876.,"params":8}
db.best_fit("H2","z",2,constraints=constraints)


