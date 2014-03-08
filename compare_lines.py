import dude_xmlutils
import astronomy_utils as astro
import xml.etree.ElementTree as et
import numpy.stats as stats
import os

"""
compare velocity shifts between absorption lines.
"""

reference='refspec.xml'

def get_files(path=os.getcwd(),ext='.xml')
  tmp = os.listdir(path)
  f = []
  for item in tmp:
    if item[-4:] == ext:
      f.append(item)

  return files

filelst = get_files()
filelst.pop(filelst.index(reference))

ids = ['M1'+i for i in range(49)]

datalst = [[dude_xmlutils.Data(iden=reference, xmlfile=reference, tag='Absorber') for ident in ids]]
for item in filelst:
  datalst.append([dude_xmlutils.Data(iden=ident, xmlfile=item, tag='Absorber') for ident in ids])

#get vel shifts wrt reference.xml

ref_dict = {}
for item in datalst[0]:
  ref_dict[item.iden]=float(item.z)

for i in range(1,len(datalst)):
  vel_dict = {}
  for ab in datalst[i]:
    vel_dict[ab.iden] = astro.get_vel_shift(ref_dict[ab.iden],ab.z)
  velshift.append(vel_dict)




