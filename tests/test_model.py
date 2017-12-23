import unittest
from copy import deepcopy

from dudeutils.data_types import Data, Region, Absorber, ContinuumPoint, ObjList, AbsorberList, ContinuumPointList, \
    RegionList
from dudeutils.model import Model
from tests.mock_types import test_xml


class ModelTestCase(unittest.TestCase):
    def setUp(self):
        self.mod = Model(xmlfile=test_xml)

    def tearDown(self):
        self.mod = None

    # todo model copy isn't coping over exactly:
    # todo model tests for equality on
    # todo    list(['xmlfile', 'chi2', 'pixels', 'params']) + list(Model.model_classes.keys())
    # todo find where copy and mod are different.
    def test_copy(self):
        copy_mod = self.mod.copy()

        self.assertIsInstance(copy_mod, Model)
        self.assertEqual(self.mod.absorber_list.objlist[0], copy_mod.absorber_list.objlist[0],
                         'elements should be equal')
        self.assertFalse(self.mod.absorber_list.objlist[0] is copy_mod.absorber_list.objlist[0],
                         'elements should not be the same')
        self.assertFalse(self.mod.absorber_list.objlist is copy_mod.absorber_list.objlist,
                         'objlists should not be equal')
        self.assertNotEqual(id(self.mod), id(copy_mod), 'ids should not be equal')
        self.assertNotEqual(id(self.mod.absorber_list), id(copy_mod.absorber_list), 'ids should differ')
        self.assertEqual(len(self.mod.absorber_list), len(copy_mod.absorber_list), 'ids should differ')
        for i in range(len(copy_mod.absorber_list)):
            self.assertIsInstance(copy_mod.absorber_list[i], Absorber)
            self.assertEqual(type(copy_mod.absorber_list[i]), type(self.mod.absorber_list[i]))
            self.assertNotEqual(id(self.mod.absorber_list[i]), id(copy_mod.absorber_list[i]), 'ids should differ')
            self.assertEqual(self.mod.absorber_list[i], copy_mod.absorber_list[i], 'abs should be equal')
        self.assertEqual(copy_mod, self.mod, "%s \n!=\n %s" % (str(copy_mod), str(self.mod)))
        for i in range(len(copy_mod.absorber_list)):
            self.mod.absorber_list[i].N = -12.3
            self.assertNotEqual(self.mod.absorber_list[i], copy_mod.absorber_list[i], 'abs should be equal')
        self.assertNotEqual(copy_mod, self.mod, "%s \n==\n %s" % (str(copy_mod), str(self.mod)))

    def test_lst_decorator(self):
        self.assertEqual(len(self.mod.absorber_list), 1)
        self.assertEqual(len(self.mod.cont_point_list), 2)
        self.assertEqual(len(self.mod.region_list), 4)

        self.assertEqual(self.mod.absorber_list, self.mod.AbsorberList)
        self.assertEqual(self.mod.cont_point_list, self.mod.ContinuumPointList)
        self.assertEqual(self.mod.region_list, self.mod.RegionList)

    def test_type_in_model(self):
        lst = Data.read(test_xml)
        self.assertTrue(isinstance(lst[4], Region), str(type(lst[4])))

    def test__get_item(self):
        lst = [type('TempClass', (), {'id': val}) for val in ['1', '2', 'a', 'b']]
        for i, val in enumerate(lst):
            self.assertEqual(i, Model._get_item(getattr(val, 'id'), lst))

    def test_split_lst(self):
        self.assertIsInstance(test_xml, str)
        lst = Data.read(test_xml)
        lengths = [len([it for it in lst if isinstance(it, Absorber)]),
                   len([it for it in lst if isinstance(it, ContinuumPoint)]),
                   len([it for it in lst if isinstance(it, Region)])]
        self.assertIsInstance(lst, list)
        objlst_lst = ObjList.split_types(lst)
        self.assertTrue(len(objlst_lst[0]) == len(objlst_lst[0].objlist) > 0)
        # need to do this since no guaranteed ordering
        self.assertTrue(len(objlst_lst[0]) in lengths, str(objlst_lst[0].objlist))
        self.assertTrue(len(objlst_lst[1]) in lengths, str(objlst_lst[1].objlist))
        self.assertTrue(len(objlst_lst[2]) in lengths, str(objlst_lst[2].objlist))
        self.assertIsInstance(objlst_lst, list)
        self.assertTrue(len(objlst_lst) == 3)
        self.assertEqual(sum([len(x) for x in objlst_lst]), len(lst))
        self.assertTrue(type(objlst_lst[0]) in [ContinuumPointList, AbsorberList, RegionList])

    def test_set_absorbers(self):
        new_abs = deepcopy(self.mod.absorber_list)
        self.assertTrue(id(new_abs) != id(self.mod.absorber_list))
        self.assertTrue(id(new_abs[0]) != id(self.mod.absorber_list[0]))
        self.assertEqual(new_abs, self.mod.absorber_list)
        self.assertEqual(new_abs[0], self.mod.absorber_list[0])
        new_abs[0].z = -1.
        self.assertNotEqual(new_abs, self.mod.absorber_list)
        self.mod.set_absorbers(new_abs)
        self.assertEqual(new_abs, self.mod.absorber_list)
        self.assertEqual(self.mod.absorber_list[0].z, -1.)


if __name__ == '__main__':
    unittest.main()
