import unittest
from model import Model, ModelDB
import data_types


class TestWrite(unittest.TestCase):
    def test_write(self):
        self.model=Model(xmlfile='/home/scott/research/J0744+2059/SiII.xml')
        iden,tag,param,val="FeII3", "Absorber","NLocked",True
        old_val=self.model.get_datum(iden,tag,param)
        self.assertFalse(old_val==val)
        self.model.set_val(iden,tag,**{param:val})
        self.assertTrue(self.model.get_datum(iden,tag,param) is val)
        self.assertTrue(self.model.get_datum(iden,tag,param)==val)
        self.model.set_val(iden,tag,**{param:old_val})
        self.assertFalse(old_val is val)
        self.assertFalse(self.model.get_datum(iden,tag,param) is val)
        self.assertFalse(self.model.get_datum(iden,tag,param)==val)
        self.assertTrue(self.model.get_datum(iden,tag,param)==old_val)

    def test_absorber(self):        
        ab=data_types.Data.read("test.xml",tag='Absorber', ids="TheID")[0]
        #self.assertTrue(ab.NLocked)
        ab.set_data(**{"NLocked":False})
        self.assertFalse(ab.NLocked)
        ab.set_data(**{"NLocked":True})
        self.assertTrue(ab.NLocked)
        ab.set_data(**{"NLocked":False})
        self.assertFalse(ab.NLocked)

    def test_readDB(self):
        db=ModelDB(name="/home/scott/research/J0744+2059/SiIIdb.xml")
        
        
if __name__ == '__main__':
    unittest.main()


