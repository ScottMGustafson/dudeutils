import dude_xmlutils as dude
import numpy as np

xmlfiles = ['228.xml','091.xml','093.xml' ,'1177.xml','125.xml','227.xml',\
    '230.xml','257.xml','549.xml','092.xml','1062.xml',\
    '124.xml','202.xml','231.xml','548.xml' ]  
#the first element will be the reference spectrum
#absorber ids to be used for the red, blue and green CCDs respectively
grn_abs  = ["M1%s"%(i) for i in range(7,22)]
red_abs  = ["M1%s"%(i) for i in range(6)] + ["M1%s"%(i) for i in range(36,40)]
blu_abs  = ["M1%s"%(i) for i in range(24,33)]


def str_rep(vb,vg,vr):
    return "vccd:  %5.3lf, %5.3lf, %5.3lf" % (vb,vg,vr)

class CCD(object):
    def __init__(self, xml, absorbers):
        self.name = xml.strip('.xml')
        self.absorbers = [ dude.Absorber(xmlfile=xml,iden=item) for item in absorbers ]
    def shift(self,ref):
        all_shifts = [self.absorbers[i].getShift(ref.absorbers[i]) for i in range(len(self.absorbers))]
        return np.mean(all_shifts), np.std(all_shifts)

red_ref = CCD(xmlfiles[0],red_abs)
blu_ref = CCD(xmlfiles[0],blu_abs)
grn_ref = CCD(xmlfiles[0],grn_abs)

for item in xmlfiles[1:]:
    red = CCD(item, red_abs).shift(red_ref)
    blu = CCD(item, blu_abs).shift(blu_ref)
    grn = CCD(item, grn_abs).shift(grn_ref)
    print("%s:  b, g, r\n--------------------------------------------"%(item.strip('.xml')))
    print(str_rep(blu[0], grn[0], red[0]))
    print("stdev: %5.3lf, %5.3lf, %5.3lf\n"%(blu[1],grn[1],red[1]))

