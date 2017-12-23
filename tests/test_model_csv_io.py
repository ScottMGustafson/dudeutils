import os
import unittest

from dudeutils.model_csv_io import ModelIO
from tests import mock_types


class FakeModel(object):
    def __init__(self):
        self.absorber_list = [mock_types.MockAb(id='ab1'), mock_types.MockAb(id='ab2'), mock_types.MockAb(id='ab3'),
                              mock_types.MockAb(id='ab4'), mock_types.MockAb(id='ab5')]
        self.cont_point_list = [mock_types.MockCont(id='cont1'), mock_types.MockCont(id='cont2'),
                                mock_types.MockCont(id='cont3'), mock_types.MockCont(id='cont4'),
                                mock_types.MockCont(id='cont5')]
        self.region_list = [mock_types.MockRegion(0., 10000.)]
        self.spec = mock_types.MockSpec()


class MyTestCase(unittest.TestCase):
    def setUp(self):
        self.model = FakeModel()
        self.path = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'test_csv')

    def tearDown(self):
        for item in os.listdir(self.path):
            if os.path.isfile(os.path.join(self.path, item)) and item.endswith('.csv'):
                os.remove(os.path.join(self.path, item))

    def test_columns(self):
        len_ = 5
        lst = [type("AbsurdClassName", (), {"id": 'hello%d' % i}) for i in range(len_)]

        cols = ModelIO.ab_column_names(lst)
        self.assertEqual(4 * len_, len(cols))
        for i in range(4):
            self.assertTrue('hello0' in cols[i])

        cols = ModelIO.cnt_column_names(lst)
        self.assertEqual(2 * len_, len(cols))
        for i in range(2):
            self.assertTrue('hello0' in cols[i])

    def test_cnt_columns(self):
        lst = ModelIO.cnt_column_names(self.model.cont_point_list)
        self.assertTrue(len(lst) == 2 * len(self.model.cont_point_list))
        self.assertEqual(lst[0], self.model.cont_point_list[0].id + '_x')

    def test_ab_columns(self):
        lst = ModelIO.ab_column_names(self.model.absorber_list)
        self.assertTrue(len(lst) == 4 * len(self.model.absorber_list))
        self.assertEqual(lst[0], self.model.absorber_list[0].id + '_ionName')
        self.assertEqual(lst[1], self.model.absorber_list[0].id + '_N')
        self.assertEqual(lst[2], self.model.absorber_list[0].id + '_b')
        self.assertEqual(lst[3], self.model.absorber_list[0].id + '_z')

    def test_header(self):
        hdr = ModelIO.header_format(self.model.absorber_list, self.model.cont_point_list)
        self.assertEqual(hdr, 'chi2,' + ','.join(
            ModelIO.cnt_column_names(self.model.cont_point_list) + ModelIO.ab_column_names(
                self.model.absorber_list)) + '\n')

    def test_ab_data(self):
        dat = ModelIO.ab_data(self.model.absorber_list)


if __name__ == '__main__':
    unittest.main()
