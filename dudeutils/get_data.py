import os

_ROOT = os.path.abspath(os.path.dirname(__file__))


def get_data(data_file):
    return os.path.join(os.path.split(_ROOT)[0], 'data', data_file)

