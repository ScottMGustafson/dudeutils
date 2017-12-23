import unittest

from dudeutils.random_sampling import get_region_objects
from tests.mock_types import test_xml, MockCont, MockSpec, MockAb, MockRegion

ab_cfg = {"D": {"N": [12.82, 12.88], "b": [13, 14]}, "H": {"N": [1, 2]}}
cont_cfg = {"cont_1": [0.01]}
glob = {'source': test_xml, 'num_processes': 1, 'num_iterations': 3}


class FakeModel(object):
    def __init__(self):
        self.absorber_list = [MockAb(id='ab1'), MockAb(id='ab2'), MockAb(id='ab3'),
                              MockAb(id='ab4'), MockAb(id='ab5')]
        self.cont_point_list = [MockCont(id='cont1'), MockCont(id='cont2'),
                                MockCont(id='cont3'), MockCont(id='cont4'),
                                MockCont(id='cont5')]
        self.region_list = [MockRegion(0., 10000.)]
        self.spec = MockSpec()


class TestRandomSample(unittest.TestCase):
    def setUp(self):
        self.model = FakeModel()

    def test_get_region_objs(self):
        ablst, cntlst = get_region_objects(self.model)
        self.assertNotEqual(id(ablst), id(self.model.absorber_list))
        self.assertNotEqual(id(cntlst), id(self.model.cont_point_list))

        # def test_assemble(self):
        #     assemble_jobs(self.model, ab_cfg, cont_cfg, **glob)


if __name__ == '__main__':
    unittest.main()
