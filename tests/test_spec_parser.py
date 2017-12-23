import unittest

from dudeutils.model import Model
from tests.mock_types import test_xml


class SpecParserTestCase(unittest.TestCase):
    def test_get_regions(self):
        model=Model(xmlfile=test_xml)
        pass

if __name__ == '__main__':
    unittest.main()
