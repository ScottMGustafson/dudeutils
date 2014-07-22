import dude_xmlutils as dude

xmlfiles = []
red_abs  = []
blu_abs  = []
grn_abs  = []

def str_rep(vb,vg,vr):
    return "vccd:  %5.3lf, %5.3lf, %5.3lf" % (vb,vg,vr)

class CCD(object):
    def __init__(self, xml, absorbers):
        self.name = xml.strip('.xml')
        self.absorbers = [ dude.Absorber(iden=item) for item in absorbers ]
    def shift(self,ref)
        all_shifts = [self.absorbers[i].getShift(ref.absorbers[i]) for i in range(len(self.absorbers))]
        return np.mean(all_shifts), np.stddev(all_shifts)


red_ref = CCD(xmlfiles[0],red_abs)
blu_ref = CCD(xmlfiles[0],blu_abs)
grn_ref = CCD(xmlfiles[0],grn_abs)

for item in xmlfiles[1:]:
    red = CCD(item, red_abs).shift(red_ref)
    blu = CCD(item, blu_abs).shift(blu_ref)
    grn = CCD(item, grn_abs).shift(grn_ref)
    print("%s:  b, g, r\n--------------------------------------------"%(item.strip('.xml')))
    print(str_rep(blu[0], grn[0], red[0]))
    print("stdev: %5.3lf, %5.3lf, %5.3lf"%(blu[1],grn[1],red[1]))

