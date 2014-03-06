import dude_xmlutils
import astronomy_utils as astro
import xml.etree.ElementTree as et

"""
compare velocity shifts between absorption lines.
"""

filelst = ['reference.xml','file1.xml','file2.xml','file3.xml']
ids = ['id1', 'id2', 'id3']

for item in filelst:
  tree = et.parse(item) #parse xml file
  root = tree.getroot()
  abslist = [dude_xmlutils.Data(iden=item, tag='Absorber').getData(root) for item in ids]
  datalst.append(abslist)

#get vel shifts wrt reference.xml

ref_dict = {}
for item in datalst[0]:
  ref_dict[item.iden]=float(item.z)

velshifts = [ref_dict]
for i in range(1,len(datalst):
  vel_dict = {}
  for ab in datalst[i]:
    vel_dict[ab.iden] = astro.get_vel_shift(ref_dict[ab.iden],ab.z)
  velshift.append(vel_dict)





