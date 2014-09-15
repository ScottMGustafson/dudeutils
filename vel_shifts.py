import dude_xmlutils as dude
import numpy as np

b_max = 18.  #max allowed line width in km/s
xmlfiles = ['228.xml','091.xml','093.xml' ,'1177.xml','125.xml','227.xml',\
    '230.xml','257.xml','549.xml','092.xml','1062.xml',\
    '124.xml','202.xml','231.xml','548.xml','548_post.xml','best_wave.xml','best_wave_2.xml']#,'414.xml' ]  
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
        mean = np.mean(all_shifts)
        std  =  np.std(all_shifts)
        return mean, std, max(np.fabs(np.array(all_shifts)))

red_ref = CCD(xmlfiles[0],red_abs)
blu_ref = CCD(xmlfiles[0],blu_abs)
grn_ref = CCD(xmlfiles[0],grn_abs)

#take out absorbers with an inappropriately large width
exclude_red = [item.iden for item in red_ref.absorbers if float(item.b) > b_max]
exclude_blu = [item.iden for item in blu_ref.absorbers if float(item.b) > b_max]
exclude_grn = [item.iden for item in grn_ref.absorbers if float(item.b) > b_max]

red_ref.absorbers = [ item for item in red_ref.absorbers if item.iden not in exclude_red ]
blu_ref.absorbers = [ item for item in blu_ref.absorbers if item.iden not in exclude_blu ]
grn_ref.absorbers = [ item for item in grn_ref.absorbers if item.iden not in exclude_grn ]

red_abs = [ item.iden for item in red_ref.absorbers]
blu_abs = [ item.iden for item in blu_ref.absorbers]
grn_abs = [ item.iden for item in grn_ref.absorbers]
 

#take out absorbers with an inappropriately large width

for item in xmlfiles[1:]:
    red = CCD(item, red_abs).shift(red_ref)
    blu = CCD(item, blu_abs).shift(blu_ref)
    grn = CCD(item, grn_abs).shift(grn_ref)
    
    print("%s:  b, g, r\n--------------------------------------------"%(item.strip('.xml')))
    print(str_rep(blu[0], grn[0], red[0]))
    print("stdev: %5.3lf, %5.3lf, %5.3lf"%(blu[1],grn[1],red[1]))
    print("max  : %5.3lf, %5.3lf, %5.3lf\n"%(blu[2],grn[2],red[2]))

