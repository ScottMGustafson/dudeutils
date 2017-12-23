import unittest

from dudeutils.model import Model
from tests.mock_types import test_xml


class RegionTestCase(unittest.TestCase):
    def setUp(self):
        self.model = Model(xmlfile=test_xml)

    def test_region(self):
        self.assertTrue(hasattr(self.model, 'RegionList'))
        self.assertTrue(len(self.model.region_list) == 4, str(self.model.region_list))

    def test_cmp(self):
        self.assertTrue(self.model.region_list[0] < self.model.region_list[1])


if __name__ == '__main__':
    unittest.main()
