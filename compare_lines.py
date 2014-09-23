import xmlutils
import astronomy_utils as astro
import xml.etree.ElementTree as et
#import numpy.stats as stats
import os
from pprint import pprint

"""
compare velocity shifts between absorption lines.
"""

reference='1177.xml'

def get_files(path=os.getcwd(),ext='.xml'):
  tmp = os.listdir(path)
  f = []
  for item in tmp:
    if item[-4:] == ext:
      f.append(item)

  return f


filelst = get_files()
filelst.insert(0, filelst.pop(filelst.index(reference)))

ids = ['M1'+str(i) for i in range(49)]

datalst = []
for item in filelst:
  datalst.append([xmlutils.Data(iden=ident, xmlfile=item, tag='Absorber', assign_ids=True) for ident in ids])

#get vel shifts wrt reference.xml 

ref_dict = {}
for item in datalst[0]:
  ref_dict[item.id]=float(item.z)

velshift = []
for i in range(len(datalst)):
  vel_dict = {}
  vel_dict['aname'] = filelst[i]
  for ab in datalst[i]:
    vel_dict[ab.id] = astro.get_vel_shift(ref_dict[ab.id],float(ab.z))
  velshift.append(vel_dict)

count = 0.
avg=0.
for item in velshift:
  try:
    print('%s: %5.3lf %5.3lf\n'%(item['aname'],item['M110']-velshift[0]['M110'],\
                                          item['M128']-velshift[0]['M128']))
  except TypeError:
    print item['M117']
    print velshift[0]['M117']
  avg+=(item['M110']-velshift[0]['M110'])-(item['M128']-velshift[0]['M128'])
  count+=1.

avg/=count

print str(avg)
#for item in velshift:
#  pprint(item)



