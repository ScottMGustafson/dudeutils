import unittest
import unittest.mock

from dudeutils.data_types import Absorber, Region, ContinuumPoint, ObjList
from dudeutils.model import Model
from dudeutils.simulated_annealing import *
from tests.mock_types import MockAb, MockSpec


class SimAnnealingTestCase(unittest.TestCase):
    def test_attr_vary(self):
        ab = MockAb(N=12, b=13, z=3, NLocked=False, bLocked=False, zLocked=True)
        old_z = ab.z
        attr_vary(ab, 'z', [5., 7.], 1.2)
        self.assertEqual(old_z, ab.z)

        old_N = ab.N
        attr_vary(ab, 'N', [5., 7.], 3.1)
        self.assertNotEqual(old_N, ab.N)

    def test_ab_in_regions(self):
        ab = Absorber(id='', N=15.1, b=45.3, z=0., ionName='H I')
        regions = ObjList.factory([Region(start=1000., end=1400.)])
        self.assertTrue(ab_in_regions(ab, regions))
        regions = ObjList.factory([Region(start=4000., end=5000.)])
        self.assertFalse(ab_in_regions(ab, regions))

    @unittest.mock.patch('dudeutils.spec_parser.Spectrum.fit_absorption')
    def test_reset_model(self, mock_ft_abs):
        mod = Model(spec=MockSpec(), pixels=75, params=8,
                    AbsorberList=[Absorber(id='id%d' % i, N=12.1, b=45.3, z=2.34, ionName='H I') for i in range(10)],
                    ContinuumPointList=[ContinuumPoint(x=1001, y=1003), ContinuumPoint(x=1010, y=1013),
                                        ContinuumPoint(x=1001, y=1103)],
                    RegionList=[Region(start=1200., end=1400.)])
        mock_ft_abs.return_value = np.array([]), np.array([]), 123.
        annealr = Anneal(mod.spec, mod)
        self.assertNotEqual(id(annealr.best_model), id(annealr.model))
        self.assertEqual(annealr.best_model.absorber_list, annealr.model.absorber_list)
        self.assertEqual(annealr.best_model, annealr.model, str(annealr.best_model) + "\n!=\n" + str(annealr.model))
        annealr.model.absorber_list = [Absorber(id='test_id', N=17.403, b=4.8, z=2.74, ionName='H I')]
        self.assertNotEqual(annealr.model.absorber_list, annealr.best_model.absorber_list)
        self.assertNotEqual(annealr.model, annealr.best_model,
                            str(annealr.best_model) + "\n==\n" + str(annealr.model))
        self.assertEqual('test_id', annealr.model.absorber_list[0].id)
        annealr.reset_model()
        self.assertEqual(annealr.model.absorber_list, annealr.best_model.absorber_list)
        self.assertNotEqual(id(mod.absorber_list[0]), id(annealr.best_model.absorber_list[0]))

    @unittest.mock.patch('dudeutils.spec_parser.Spectrum.fit_absorption')
    def test_reset_best(self, mock_ft_abs):
        mod = Model(spec=MockSpec(), pixels=75, params=8,
                    AbsorberList=[Absorber(id='id%d' % i, N=12.1, b=45.3, z=2.34, ionName='H I') for i in range(10)],
                    ContinuumPointList=[ContinuumPoint(x=1001, y=1003), ContinuumPoint(x=1010, y=1013),
                                        ContinuumPoint(x=1001, y=1103)],
                    RegionList=ObjList.factory([Region(start=1200., end=1400.)]))
        mock_ft_abs.return_value = np.array([]), np.array([]), 123.
        annealr = Anneal(mod.spec, mod)
        self.assertNotEqual(id(annealr.best_model), id(annealr.model))
        self.assertEqual(annealr.best_model, annealr.model)
        annealr.model.absorber_list[0].N = -12.3
        self.assertNotEqual(annealr.best_model, annealr.model)
        annealr.reset_best()
        self.assertEqual(annealr.best_model, annealr.model)
        self.assertEqual(annealr.model.absorber_list[0].N, -12.3)
        self.assertEqual(annealr.best_model.absorber_list[0].N, -12.3)

    @unittest.mock.patch('dudeutils.spec_parser.Spectrum.fit_absorption')
    def test_randomize_parameters(self, mock_ft_abs):
        mod = Model(spec=MockSpec(), pixels=75, params=8,
                    AbsorberList=[Absorber(id='id%d' % i, ionName='H I',
                                           N=12.1, b=45.3, z=2.34,
                                           NLocked=False, bLocked=False, zLocked=False) for i in range(10)],
                    ContinuumPointList=[ContinuumPoint(x=1001, y=1003), ContinuumPoint(x=1010, y=1013),
                                        ContinuumPoint(x=1001, y=1103)],
                    RegionList=ObjList.factory([Region(start=1000., end=6000.)]))
        mock_ft_abs.return_value = np.array([]), np.array([]), 123.
        annealr = Anneal(getattr(mod, 'spec'), mod)
        annealr.model.absorber_list[0].NLocked = True
        old_N = float(annealr.model.absorber_list[0].N)
        old_b = float(annealr.model.absorber_list[0].b)
        annealr.randomize_parameters()
        self.assertEqual(old_N, annealr.model.absorber_list[0].N)
        self.assertNotEqual(old_b, annealr.model.absorber_list[0].b)

    @unittest.mock.patch('dudeutils.spec_parser.Spectrum.fit_absorption')
    def test_randomize_parameters_dh_tied(self, mock_ft_abs):
        mod = Model(spec=MockSpec(), pixels=75, params=8,
                    AbsorberList=[Absorber(id='D', ionName='H I',
                                           N=12.1, b=45.3, z=2.34,
                                           NLocked=False, bLocked=False, zLocked=False),

                                  Absorber(id='H', ionName='H I',
                                           N=12.1, b=45.3, z=0.01,
                                           NLocked=False, bLocked=False, zLocked=True)
                                  ],
                    ContinuumPointList=[ContinuumPoint(x=1001, y=1003), ContinuumPoint(x=1010, y=1013),
                                        ContinuumPoint(x=1001, y=1103)],
                    RegionList=ObjList.factory([Region(start=1000., end=6000.)]))

        mock_ft_abs.return_value = np.array([]), np.array([]), 123.
        annealr = Anneal(getattr(mod, 'spec'), mod, tie=('H_z', 'D_z'))
        annealr.randomize_parameters()
        self.assertEqual(mod.absorber_list[0].z, mod.absorber_list[1].z)

        annealr.model.absorber_list[0].zLocked = True
        annealr.model.absorber_list[0].z = 1.23
        annealr.model.absorber_list[1].zLocked = False
        annealr.model.absorber_list[1].z = 2.63
        annealr.randomize_parameters()
        self.assertEqual(mod.absorber_list[0].z, mod.absorber_list[1].z)


    def test_deepcopy(self):
        a = [MockAb(id='test_id%d' % i) for i in range(7)]
        b = copy.deepcopy(a)
        self.assertNotEqual(id(a), id(b))
        self.assertEqual(a, b)


if __name__ == '__main__':
    unittest.main()
