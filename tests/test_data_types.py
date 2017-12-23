import unittest
from copy import deepcopy

from dudeutils.data_types import Data, Region, Absorber, ContinuumPoint, ObjList
from tests.mock_types import test_xml


class DataTestCase(unittest.TestCase):
    def setUp(self):
        self.lst = Data.read(test_xml)

    def test_objlist_factory_raises_exception(self):
        self.assertRaises(Exception, ObjList.factory, [self.lst])

    def test_split_types(self):
        self.assertEqual(len(ObjList.split_types(self.lst)), 3)

    def test_type_in_model(self):
        self.assertEqual(len(self.lst), 7)
        self.assertIsInstance(self.lst[4], Region)
        self.assertIsInstance(self.lst[2], Absorber)
        self.assertIsInstance(self.lst[0], ContinuumPoint)

    def test_ab_constructor(self):
        ab = Absorber(id='id', ionName='C IV',
                      N=13.45, NLocked=True, NError=0,
                      b=13.45, bLocked=True, bError=0,
                      z=3.45, zLocked=True, zError=0)
        self.assertTrue(isinstance(ab, Absorber))
        self.assertTrue(len(ab.get_lines()) == 2)
        ab2 = Absorber(id='id', ionName='C IV',
                       N=13.45, NLocked=True, NError=0,
                       b=13.45, bLocked=True, bError=0,
                       z=3.45, zLocked=True, zError=0)
        self.assertEqual(ab, ab2)
        self.assertNotEqual(id(ab), id(ab2))

        # with incomplete signature...
        # ab = Absorber(ionName='C IV', N=13.45, b=13.45, z=3.45)
        # self.assertEqual(ab.NError, 0)
        # self.assertEqual(ab.zLocked, True)
        #
        # with self.assertRaises(Exception):
        #     Absorber(id='id',
        #              b=13.45, bLocked=True, bError=0,
        #              z=3.45, zLocked=True, zError=0)

    def test_constructor_region(self):
        reg1 = Region(start=123., end=456.)
        self.assertTrue(isinstance(reg1, Region))
        reg2 = Region(start=123., end=456.)
        self.assertEqual(reg1, reg2)
        self.assertNotEqual(id(reg1), id(reg2))

    def test_constructor_cont(self):
        cont1 = ContinuumPoint(x=123., y=456., xLocked=True, yLocked=True, xError=0, yError=0)
        self.assertTrue(isinstance(cont1, ContinuumPoint))
        cont2 = ContinuumPoint(x=123., y=456., xLocked=True, yLocked=True, xError=0, yError=0)
        self.assertEqual(cont1, cont2)
        self.assertNotEqual(id(cont1), id(cont2))

    def test_deep_copy(self):
        tmp = deepcopy(self.lst)
        for i in range(len(self.lst)):
            self.assertFalse(self.lst[i] is tmp[i])
            self.assertEqual(self.lst[i], tmp[i])

    def test_from_node(self):
        attrib_dct = {'ionName': 'C IV',
                      'N': float(11.2), 'b': float(2.02), 'z': float(1.92),
                      'NError': 0, 'bError': 0, 'zError': 0.0,
                      'NLocked': False, 'bLocked': False, 'zLocked': False, 'id': ""}
        test_node = type('testnode', (), {"tag": 'Absorber', 'attrib': attrib_dct})
        ab = Absorber.from_node(node=test_node)
        for attr in Absorber.node_attrib:
            self.assertTrue(hasattr(ab, attr))


class ObjListTestCase(unittest.TestCase):
    def setUp(self):
        self.obj = ObjList([1, 2, 3, 4], a=1, b=2)

    def test_instantiation(self):
        self.assertTrue(len(self.obj) == 4 == len(self.obj.objlist))
        self.assertTrue(self.obj.a == 1 and self.obj.b == 2)
        self.assertEqual(self.obj.objlist, [1, 2, 3, 4])
        self.assertTrue(2 in self.obj)
        self.assertEqual(self.obj.objlist, [x for x in self.obj])
        _new = ObjList([1, 2, 3, 4], a=1, b=2) + ObjList([5, 6, 7, 8], a=1, b=2)
        self.assertEqual(_new.objlist, [1, 2, 3, 4, 5, 6, 7, 8])

    def test_eq(self):
        other = ObjList([1, 2, 3, 4], a=1, b=2)
        self.assertEqual(self.obj, other)
        other_ = ObjList([1, 1, 3, 4], a=9, b=2)
        self.assertNotEqual(self.obj.objlist, other_.objlist)
        self.assertNotEqual(self.obj, other_)


if __name__ == '__main__':
    unittest.main()
